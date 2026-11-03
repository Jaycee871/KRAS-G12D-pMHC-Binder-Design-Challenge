#!/usr/bin/env python3
"""Re-rank every scored design on the paired same-pose MM-GBSA ddG.

Run only after `scripts/ddg_control_7ow6.py` passes - a metric that cannot
recover the known answer on 7OW6 has no business ranking designs
(`results/validation/CONTROL_7OW6.md`).

Inputs, in the order a design's pose is taken from:

* `results/constrained-rescore/<backbone>_g12d.cif` - the pocket-constrained
  batch-03 poses, the only poses the 7OW6 control validated;
* `results/batch-0{1,2,3}/scored_structures/<sha1>_g12d.cif` - the
  unconstrained poses, taken for designs that have no constrained pose.

A design whose binder has no heavy atom within `LOCAL_RADIUS` of peptide
position 6 is reported as `no_site_contact`, not as ddg 0: the metric has
nothing to measure there, and that is itself the verdict (the binder is parked
off the peptide).

Usage:
    python3 scripts/ddg_rescore_designs.py --replicates 1 --workers 6
    python3 scripts/ddg_rescore_designs.py --only-top 12 --replicates 5   # finalists
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import sys
import time
from multiprocessing import get_context

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, "src"))

OUT_DIR = os.path.join(REPO, "results", "ddg-scoring")
BATCHES = ("batch-01", "batch-02", "batch-03")


def sha12(sequence: str) -> str:
    return hashlib.sha1(sequence.encode()).hexdigest()[:12]


def collect_designs():
    """Every scored design across batches 01-03, newest pose first."""
    constrained = {}
    cr = os.path.join(REPO, "results", "constrained-rescore")
    for line in open(os.path.join(cr, "checkpoint.jsonl")):
        row = json.loads(line)
        path = os.path.join(cr, f"{row['backbone']}_g12d.cif")
        if os.path.exists(path):
            constrained[row["sequence"]] = (row["backbone"], path)

    designs = {}
    for batch in BATCHES:
        ck = os.path.join(REPO, "results", batch, "checkpoint.jsonl")
        if not os.path.exists(ck):
            continue
        for line in open(ck):
            rec = json.loads(line)
            row = rec.get("row") or {}
            seq = row.get("sequence")
            if not seq:
                continue
            struct = rec.get("structure_file") or ""
            path = os.path.join(REPO, "results", batch, struct) if struct else ""
            name = row.get("name") or row.get("backbone") or sha12(seq)
            entry = designs.setdefault(seq, {
                "sequence": seq, "batches": [], "name": name,
                "backbone": row.get("backbone"), "strategy": row.get("strategy"),
                "pose": path, "pose_kind": "unconstrained",
                "iptm_g12d": row.get("iptm_g12d"),
                "specificity_gap": row.get("specificity_gap"),
                "peptide_buried_frac": row.get("peptide_buried_frac"),
                "pos6_contacts": row.get("pos6_contacts"),
                "asp6_hbond": row.get("asp6_hbond"),
                "asp6_salt_bridge": row.get("asp6_salt_bridge"),
                "net_interface_charge": row.get("net_interface_charge"),
            })
            if batch not in entry["batches"]:
                entry["batches"].append(batch)
            if not entry["pose"] and path:
                entry["pose"] = path
    for seq, (backbone, path) in constrained.items():
        entry = designs.setdefault(seq, {
            "sequence": seq, "batches": ["batch-03"], "name": backbone,
            "backbone": backbone, "strategy": backbone.rsplit("_", 1)[0],
            "iptm_g12d": None, "specificity_gap": None,
            "peptide_buried_frac": None, "pos6_contacts": None,
            "asp6_hbond": None, "asp6_salt_bridge": None,
            "net_interface_charge": None,
        })
        entry["pose"] = path
        entry["pose_kind"] = "constrained"
        entry["backbone"] = backbone
        entry["name"] = backbone
    return [d for d in designs.values() if d.get("pose")]


def _job(spec):
    os.environ.setdefault("OPENMM_CPU_THREADS", "2")
    import ddg_scoring as dd
    out = {k: spec[k] for k in (
        "name", "backbone", "strategy", "sequence", "pose", "pose_kind",
        "iptm_g12d", "specificity_gap", "peptide_buried_frac", "pos6_contacts",
        "asp6_hbond", "asp6_salt_bridge", "net_interface_charge", "batches")}
    t0 = time.time()
    try:
        model = dd.load_model(spec["pose"])
        groups = dd.groups_from_model(model, binder_sequence=spec["sequence"])
        result = dd.paired_ddg(model, groups, label=spec["name"],
                               replicates=spec["replicates"])
        out.update(ok=True, ddg_sel=result.ddg_sel, ddg_sd=result.ddg_sd,
                   ddg_values=result.ddg_values, dg_g12d=result.dg_g12d,
                   dg_wt=result.dg_wt, n_local_residues=result.n_residues,
                   seconds=result.seconds)
    except Exception as exc:                      # noqa: BLE001
        kind = "no_site_contact" if "does not touch the site" in str(exc) else "error"
        out.update(ok=False, status=kind, error=f"{type(exc).__name__}: {exc}",
                   seconds=time.time() - t0)
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--replicates", type=int, default=1)
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--only-top", type=int, default=0,
                    help="re-run only the N best designs of an existing pass")
    ap.add_argument("--resume", default="",
                    help="jsonl to resume from (skips labels already present)")
    ap.add_argument("--out", default=os.path.join(OUT_DIR, "designs_ddg.jsonl"))
    args = ap.parse_args()
    os.makedirs(OUT_DIR, exist_ok=True)

    designs = collect_designs()
    if args.only_top:
        prior = [json.loads(l) for l in open(args.resume or args.out)]
        scored = sorted((r for r in prior if r.get("ok")),
                        key=lambda r: r["ddg_sel"])[:args.only_top]
        keep = {r["sequence"] for r in scored}
        designs = [d for d in designs if d["sequence"] in keep]
    done = set()
    if args.resume and os.path.exists(args.resume) and not args.only_top:
        done = {json.loads(l)["sequence"] for l in open(args.resume)}
        designs = [d for d in designs if d["sequence"] not in done]
    for d in designs:
        d["replicates"] = args.replicates

    print(f"{len(designs)} designs to score "
          f"({len(done)} already done), {args.replicates} replicate(s) each",
          flush=True)
    mode = "a" if (args.resume == args.out and done) else "w"
    ctx = get_context("spawn")
    t0 = time.time()
    with open(args.out, mode) as fh, ctx.Pool(processes=args.workers) as pool:
        for i, row in enumerate(pool.imap_unordered(_job, designs), start=1):
            fh.write(json.dumps(row) + "\n")
            fh.flush()
            os.fsync(fh.fileno())
            if row.get("ok"):
                print(f"[{i}/{len(designs)}] {row['name']:<22} "
                      f"ddg {row['ddg_sel']:+7.2f} sd {row['ddg_sd']:.2f} "
                      f"({time.time()-t0:.0f}s)", flush=True)
            else:
                print(f"[{i}/{len(designs)}] {row['name']:<22} "
                      f"{row.get('status')}: {row.get('error','')[:70]}", flush=True)
    print(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
