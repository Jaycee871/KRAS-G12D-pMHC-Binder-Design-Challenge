#!/usr/bin/env python3
"""Does the energy function detect a carboxylate it *can* see? (7OW6 groove)

The 7OW6 TCR control failed, but a contact census showed why: the nearest TCR
heavy atom is 5.7 A from the Asp6 **side chain** and makes zero contacts within
5 A of it. Every "position-6 contact" the pipeline counts is to the peptide
*backbone*, which both states of a pair share atom-for-atom. So the TCR control
asks the metric to score an interaction that does not exist.

That failure is ambiguous on its own: it could mean the metric is ill-posed for
this epitope, or that the MM-GBSA energy function cannot see a carboxylate at
all. This script separates the two using a positive control that is *inside the
same crystal*: HLA-A*11:01 **Arg114 NH2 sits 2.95 A from Asp6 OD2**, with 15
heavy-atom contacts within 4 A. The groove grips the carboxylate the way a
selective binder is supposed to.

So: take the HLA heavy chain as the "binder" and the peptide as the target, and
run the identical paired calculation. If ddg_sel comes out clearly negative, the
energy function works and the TCR control is ill-posed. If it does not, the
energy function itself cannot score this chemistry and nothing built on it can
be trusted.

Usage:
    python3 scripts/ddg_groove_control.py [--replicates 3]
"""
from __future__ import annotations

import argparse, json, os, sys, time
from multiprocessing import get_context

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, "src"))
OUT_DIR = os.path.join(REPO, "results", "ddg-scoring")


def _job(spec):
    os.environ["OPENMM_CPU_THREADS"] = "2"
    os.environ["DDG_SOLUTE_DIELECTRIC"] = str(spec["eps_in"])
    import ddg_scoring as dd
    dd.SOLUTE_DIELECTRIC = float(spec["eps_in"])
    t0 = time.time()
    try:
        model = dd.load_model(os.path.join(REPO, "data/raw/7OW6.pdb"))
        # keep only the HLA heavy chain and the peptide
        for cid in [c.id for c in model]:
            if cid not in ("A", "C"):
                model.detach_child(cid)
        groups = dd.groups_from_model(model, binder_chains=["A"],
                                      peptide_position=spec["peptide_position"])
        res = dd.paired_ddg(model, groups, label=spec["label"],
                            replicates=spec["replicates"])
        row = res.as_dict(); row["ok"] = True
    except Exception as exc:                      # noqa: BLE001
        row = {"label": spec["label"], "ok": False,
               "error": f"{type(exc).__name__}: {exc}", "seconds": time.time()-t0}
    row.update({k: spec[k] for k in ("label", "eps_in", "peptide_position", "note")})
    return row


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--replicates", type=int, default=3)
    ap.add_argument("--out", default=os.path.join(OUT_DIR, "groove_control.jsonl"))
    args = ap.parse_args()
    os.makedirs(OUT_DIR, exist_ok=True)
    done = set()
    if os.path.exists(args.out):
        done = {json.loads(l)["label"] for l in open(args.out) if l.strip()}

    specs = [
        dict(label="groove_asp6_eps1", eps_in=1.0, peptide_position=6,
             note="HLA as binder, Arg114 2.95 A from Asp6 OD2 - the salt bridge is real"),
        dict(label="groove_asp6_eps4", eps_in=4.0, peptide_position=6,
             note="same, with the MM-GBSA buried-charge remedy"),
        dict(label="groove_offsite_p3_eps1", eps_in=1.0, peptide_position=3,
             note="negative: Val3, a buried apolar anchor, no salt bridge"),
    ]
    specs = [s for s in specs if s["label"] not in done]
    for s in specs:
        s["replicates"] = args.replicates
    print(f"{len(specs)} groove jobs, {args.replicates} replicates each", flush=True)

    ctx = get_context("spawn")
    with open(args.out, "a") as fh, ctx.Pool(processes=min(3, len(specs))) as pool:
        for row in pool.imap_unordered(_job, specs):
            fh.write(json.dumps(row) + "\n"); fh.flush(); os.fsync(fh.fileno())
            if row.get("ok"):
                print(f"{row['label']:<24} ddg {row['ddg_sel']:+6.2f} sd {row['ddg_sd']:.2f} "
                      f"{[round(v,2) for v in row['ddg_values']]}", flush=True)
            else:
                print(f"{row['label']:<24} FAILED {row['error']}", flush=True)
    print(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
