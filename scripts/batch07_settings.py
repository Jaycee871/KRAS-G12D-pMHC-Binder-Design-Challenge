#!/usr/bin/env python3
"""Batch-07 step 1: a systematic per-design Boltz-2 sampling sweep for bb54-g2.

`bb54-g2` (`B6-c6-46_w2_bb54_donor-anchored-g2`) is the only design in this
project that has ever shown wild-type discrimination on measured geometry
(-59.2 % peptide buried area, -73.1 % position-6 contacts; batch-06 Result 3),
and it has stalled at **2 of 3** on the Asp6 direct read-out gate. Batch-06
also left one specific, untried observation on the table: its *tightest pose*
came under `combo` (without_potentials + recycling_steps 8, 3.74/5.10 A) while
its *gate passes* came under `recyc8` (7.26/7.96 A). The best sampling setting
is design-specific and the two knobs are **not** additive, so this script walks
the grid around them for this one design.

Grid (each cell = ONE Boltz-2 call with `diffusion_samples: 3`, i.e. three
independent diffusion trajectories that cannot be byte-identical duplicates of
one another):

    recycling_steps in {1, 3, 5, 6, 8, 10}  x  without_potentials in {False, True}

plus two longer-trajectory cells at the two most promising recycling values and
one independent repeat of the best cell, to test whether a 3-of-3 is
reproducible or a lucky draw.

Nothing here is a new pipeline: the call body, the validated flat pocket
constraint, `step_scale` 1.638, the strictly-serial >= 8 s pacing and the whole
scoring path are `scripts/batch06_panel.run_draws` -> `run_batch04.score_candidate`,
re-used verbatim with `OUT` pointed at `results/batch-07`. No gate or threshold
is changed; `src/design_gates.py` is untouched. No ipTM gap and no MM-GBSA ddG
is computed anywhere - both are retired.
"""
from __future__ import annotations

import argparse, json, sys, time
from pathlib import Path
from typing import Any, Dict, List

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(REPO / "scripts"))

import batch06_panel as P  # noqa: E402
import run_batch04 as batch  # noqa: E402
import run_pipeline  # noqa: E402
import sampling_study as study  # noqa: E402
from nvidia_client import NvidiaClient  # noqa: E402

OUT = REPO / "results" / "batch-07"
P.OUT = OUT
batch.OUT = OUT

SUBJECT = "B6-c6-46_w2_bb54_donor-anchored-g2"
B06 = REPO / "results" / "batch-06"

# name: (recycling_steps, sampling_steps, without_potentials)
SWEEP: Dict[str, tuple] = {
    "r1":        (1, 50, False),
    "r1np":      (1, 50, True),
    "r3":        (3, 50, False),    # Boltz-2 default, batch-04/05 baseline
    "r3np":      (3, 50, True),     # = sampling-study `nopot`
    "r5":        (5, 50, False),
    "r5np":      (5, 50, True),
    "r6":        (6, 50, False),
    "r6np":      (6, 50, True),
    "r8":        (8, 50, False),    # = `recyc8`, batch-06's 2/3 cell
    "r8np":      (8, 50, True),     # = `combo`,  batch-06's 1/3 cell
    "r10":       (10, 50, False),
    "r10np":     (10, 50, True),
    "r8s100":    (8, 100, False),   # longer diffusion trajectory
    "r10s100np": (10, 100, True),
}
# cells re-run a second time, independently, to separate a reproducible
# setting from a lucky draw. Filled from --repeat on the command line.
REPEAT_SUFFIX = "_rep"

study.CONDITIONS.update(SWEEP)


def subject_design() -> Dict[str, Any]:
    rows = [json.loads(l) for l in
            (B06 / "panel_checkpoint.jsonl").read_text().splitlines() if l.strip()]
    src = next(r for r in rows if r["design"] == SUBJECT and r["arm"] == "g12d")
    row = src["row"]
    return {"name": SUBJECT, "sequence": src["sequence"],
            "backbone": src.get("backbone") or row.get("backbone", ""),
            "variant": src.get("variant") or row.get("variant", ""),
            "binder_length": row.get("binder_length", 46),
            "strategy": row.get("strategy", ""),
            "mpnn_score": row.get("mpnn_score", ""),
            "donor_bias_score": row.get("donor_bias_score", ""),
            "design_pos6_contacts": row.get("design_pos6_contacts", ""),
            "design_pos6_distance": row.get("design_pos6_distance", "")}


def main(argv=None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cells", default=",".join(SWEEP))
    parser.add_argument("--repeat", default="",
                        help="comma-separated cells to run a second, "
                             "independent time (suffix _rep)")
    parser.add_argument("--draws", type=int, default=3)
    parser.add_argument("--timeout", type=float, default=900.0)
    parser.add_argument("--max-boltz-calls", type=int, default=20)
    args = parser.parse_args(argv)

    (OUT / "structures").mkdir(parents=True, exist_ok=True)
    states = run_pipeline.load_target_states(REPO / "data" / "targets")
    client = NvidiaClient(timeout=args.timeout)
    design = subject_design()
    print(f"[step1] subject {design['name']} ({design['binder_length']} aa)")

    cells: List[str] = [c.strip() for c in args.cells.split(",") if c.strip()]
    for cell in [c.strip() for c in args.repeat.split(",") if c.strip()]:
        study.CONDITIONS[cell + REPEAT_SUFFIX] = SWEEP[cell]
        cells.append(cell + REPEAT_SUFFIX)

    calls, started = 0, time.time()
    for cell in cells:
        if calls >= args.max_boltz_calls:
            print("[stop] Boltz-2 call cap reached")
            break
        already = P.done_keys()
        recyc, samp, nopot = study.CONDITIONS[cell]
        print(f"\n[step1] === {cell}: recycling_steps={recyc} "
              f"sampling_steps={samp} without_potentials={nopot} ===", flush=True)
        calls += P.run_draws(client, design, "g12d", cell, states,
                             args.draws, already)
        rows = P.arm_rows(P.load_rows(), design["name"], "g12d", cell)
        passes = sum(1 for r in rows
                     if str(r.get("asp6_gate_pass")).lower() == "true")
        print(f"  -> {cell}: Asp6 gate {passes}/{len(rows)}", flush=True)

    print(f"\n[step1] {calls} Boltz-2 calls in {time.time() - started:.0f}s")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
