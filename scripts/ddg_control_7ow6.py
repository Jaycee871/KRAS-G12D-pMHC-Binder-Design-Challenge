#!/usr/bin/env python3
"""Positive and negative controls for the paired same-pose MM-GBSA ddG.

The gate on everything else, exactly as with the constrained Boltz-2 protocol
(`results/validation/CONTROL_7OW6.md`): the TCR of Poole et al. 2022 is
experimentally G12D-selective, so a metric that is going to rank designs has to
give it a clearly favourable (negative) `ddg_sel`, has to reproduce, and has to
be *silent* on controls where there is nothing to discriminate.

Jobs
----
positive
  `crystal_7ow6`        the deposited 7OW6 complex, TCR (chains D+E) as binder
  `predicted_7ow6`      the pocket-constrained Boltz-2 prediction of the same

negative / specificity
  `crystal_offsite_p3`  same crystal, but the paired mutation is peptide
                        position 3 (Val -> Gly) instead of position 6. The TCR
                        is selective for position 6, so a metric that reports
                        the same magnitude here is reporting "a side chain was
                        deleted", not selectivity.
  `decoy_7ow6`          the TCR on the unrelated A*11:01 decoy `AVFDRKSDAK`,
                        paired at position 6 (Lys -> Gly).
  `offpeptide_design`   a batch-01 design that Boltz-2 docked away from the
                        peptide. With no binder atom near position 6 the metric
                        must either refuse to score or return ~0.

The peptide-free (apo) state has no position 6 at all, so a *paired* ddG cannot
be defined for it; it is reported as not-applicable rather than faked.

Usage:
    python3 scripts/ddg_control_7ow6.py [--replicates 5] [--workers 5]
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
from multiprocessing import get_context

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, "src"))

OUT_DIR = os.path.join(REPO, "results", "ddg-scoring")


def _job(spec):
    os.environ.setdefault("OPENMM_CPU_THREADS", "2")
    import ddg_scoring as dd
    label = spec["label"]
    t0 = time.time()
    try:
        model = dd.load_model(os.path.join(REPO, spec["path"]))
        groups = dd.groups_from_model(
            model,
            binder_chains=spec.get("binder_chains"),
            peptide_position=spec.get("peptide_position", 6),
        )
        result = dd.paired_ddg(
            model, groups, label=label, replicates=spec["replicates"])
        row = result.as_dict()
        row["ok"] = True
        row["kind"] = spec["kind"]
        row["note"] = spec.get("note", "")
        return row
    except Exception as exc:                      # noqa: BLE001 - reported, not hidden
        return {
            "label": label, "ok": False, "kind": spec["kind"],
            "error": f"{type(exc).__name__}: {exc}",
            "note": spec.get("note", ""), "seconds": time.time() - t0,
        }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--replicates", type=int, default=5)
    ap.add_argument("--workers", type=int, default=5)
    ap.add_argument("--out", default=os.path.join(OUT_DIR, "control_7ow6_ddg.json"))
    args = ap.parse_args()
    os.makedirs(OUT_DIR, exist_ok=True)

    reps = args.replicates
    specs = [
        dict(label="crystal_7ow6", kind="positive",
             path="data/raw/7OW6.pdb", binder_chains=["D", "E"],
             peptide_position=6, replicates=reps,
             note="deposited crystal pose, the ground truth"),
        dict(label="predicted_7ow6", kind="positive",
             path="results/validation/7ow6_control_g12d_constrained.cif",
             binder_chains=["D", "E"], peptide_position=6, replicates=reps,
             note="pocket-constrained Boltz-2 prediction of the same complex"),
        dict(label="crystal_offsite_p3", kind="negative",
             path="data/raw/7OW6.pdb", binder_chains=["D", "E"],
             peptide_position=3, replicates=max(3, reps - 2),
             note="same crystal, paired mutation moved to peptide position 3"),
        dict(label="decoy_7ow6", kind="negative",
             path="results/validation/7ow6_control_decoy_constrained.cif",
             binder_chains=["D", "E"], peptide_position=6,
             replicates=max(3, reps - 2),
             note="TCR on the unrelated AVFDRKSDAK decoy, paired at position 6"),
        dict(label="offpeptide_design", kind="negative",
             path=os.environ.get(
                 "DDG_OFFPEPTIDE",
                 "results/batch-01/scored_structures/6cb2e86fe8c5_g12d.cif"),
             peptide_position=6, replicates=max(3, reps - 2),
             note="a batch-01 design docked away from the peptide"),
    ]

    ctx = get_context("spawn")
    with ctx.Pool(processes=min(args.workers, len(specs))) as pool:
        rows = pool.map(_job, specs)

    payload = {
        "generated": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        "method": "paired same-pose MM-GBSA (Amber ff14SB + GBn2), "
                  "src/ddg_scoring.py",
        "ddg_convention": "ddg_sel = dG_bind(G12D) - dG_bind(pos6->Gly); "
                          "negative favours G12D",
        "replicates": reps,
        "results": rows,
    }
    with open(args.out, "w") as fh:
        json.dump(payload, fh, indent=2)
    for row in rows:
        if row.get("ok"):
            print(f"{row['label']:>20}  ddg {row['ddg_sel']:+7.2f} "
                  f"sd {row['ddg_sd']:.2f}  reps {row['ddg_values']}  "
                  f"{row['seconds']:.0f}s")
        else:
            print(f"{row['label']:>20}  FAILED {row['error']}")
    print(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
