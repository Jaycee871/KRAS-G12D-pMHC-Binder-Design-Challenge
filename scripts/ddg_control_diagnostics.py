#!/usr/bin/env python3
"""Diagnose why the paired same-pose ddG fails its 7OW6 control.

The first control panel (`results/ddg-scoring/control_7ow6_ddg.json`) came out
inverted: the two positive controls scored +1.27 and +0.81 (wrong sign) and the
two negative controls scored -0.83 and -3.53, i.e. the metric "prefers" G12D
most strongly on an unrelated decoy peptide. Two mechanisms could do that, and
they are separable:

A. **Buried-charge over-penalty.** GB with a solute dielectric of 1 is known to
   over-pay the desolvation of a charge at a protein-protein interface. Here the
   charge is exactly the thing being measured (the Asp6 carboxylate), so the
   G12D state is penalised relative to a Gly that has no charge to desolvate.
   `src/ddg_scoring.py` already exposes `DDG_SOLUTE_DIELECTRIC` for this, with a
   note that it is only ever moved against this control.

B. **Unconverged minimisation.** dG_bind is ~-75 kcal/mol and the ddG is ~1, a
   1.5 % difference of two large numbers. At 500 L-BFGS iterations and a 1
   kJ/mol/nm tolerance a 2500-atom system is nowhere near a stationary point, so
   each replicate stops at a different place on the way down and the scatter
   (observed SD 0.8-1.9) is numerical, not physical.

This script runs the ground-truth crystal pose under the 2x2 of those two, plus
the decoy under the best corner, so the write-up can say which (if either) is
responsible. Every row is flushed to a jsonl checkpoint as it lands.

Usage:
    python3 scripts/ddg_control_diagnostics.py [--replicates 3] [--workers 4]
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
    dd.MIN_ITERATIONS = int(spec["min_iterations"])
    dd.MIN_TOLERANCE = float(spec["min_tolerance"])
    t0 = time.time()
    try:
        model = dd.load_model(os.path.join(REPO, spec["path"]))
        groups = dd.groups_from_model(
            model, binder_chains=spec.get("binder_chains"),
            peptide_position=spec.get("peptide_position", 6))
        res = dd.paired_ddg(model, groups, label=spec["label"],
                            replicates=spec["replicates"])
        row = res.as_dict()
        row.update(ok=True)
    except Exception as exc:                       # noqa: BLE001
        row = {"label": spec["label"], "ok": False,
               "error": f"{type(exc).__name__}: {exc}", "seconds": time.time()-t0}
    row.update({k: spec[k] for k in
                ("label", "kind", "eps_in", "min_iterations", "min_tolerance", "note")})
    return row


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--replicates", type=int, default=3)
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--out", default=os.path.join(OUT_DIR, "control_diagnostics.jsonl"))
    args = ap.parse_args()
    os.makedirs(OUT_DIR, exist_ok=True)

    done = set()
    if os.path.exists(args.out):
        done = {json.loads(l)["label"] for l in open(args.out) if l.strip()}

    CRYSTAL = dict(path="data/raw/7OW6.pdb", binder_chains=["D", "E"], peptide_position=6)
    DECOY = dict(path="results/validation/7ow6_control_decoy_constrained.cif",
                 binder_chains=["D", "E"], peptide_position=6)
    TIGHT = dict(min_iterations=5000, min_tolerance=0.05)
    LOOSE = dict(min_iterations=500, min_tolerance=1.0)

    specs = [
        dict(label="crystal_eps1_tight", kind="positive", eps_in=1.0, **CRYSTAL, **TIGHT,
             note="B alone: does tighter convergence remove the scatter?"),
        dict(label="crystal_eps4_loose", kind="positive", eps_in=4.0, **CRYSTAL, **LOOSE,
             note="A alone: does eps_in=4 restore the sign?"),
        dict(label="crystal_eps4_tight", kind="positive", eps_in=4.0, **CRYSTAL, **TIGHT,
             note="A+B together, the best case for the metric"),
        dict(label="decoy_eps4_tight", kind="negative", eps_in=4.0, **DECOY, **TIGHT,
             note="discrimination under the best corner"),
    ]
    specs = [s for s in specs if s["label"] not in done]
    for s in specs:
        s["replicates"] = args.replicates
    print(f"{len(specs)} diagnostic jobs ({len(done)} already done), "
          f"{args.replicates} replicates each", flush=True)

    ctx = get_context("spawn")
    t0 = time.time()
    with open(args.out, "a") as fh, ctx.Pool(processes=args.workers) as pool:
        for row in pool.imap_unordered(_job, specs):
            fh.write(json.dumps(row) + "\n"); fh.flush(); os.fsync(fh.fileno())
            if row.get("ok"):
                print(f"{row['label']:<20} eps{row['eps_in']:.0f} it{row['min_iterations']:<5} "
                      f"ddg {row['ddg_sel']:+6.2f} sd {row['ddg_sd']:.2f} "
                      f"{[round(v,2) for v in row['ddg_values']]} ({time.time()-t0:.0f}s)", flush=True)
            else:
                print(f"{row['label']:<20} FAILED {row['error']}", flush=True)
    print(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
