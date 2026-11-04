#!/usr/bin/env python3
"""Batch-08 stage 8: a diagnostic apo/decoy read on the most peptide-centred designs.

The normal gate order in `scripts/batch08_panel.py` runs apo and decoy only for
designs that clear the Asp6 3-of-3 gate *and* the WT criterion, and in batch-08
nothing did, so the batch would otherwise carry no apo arm at all.  The open
question of this project - every design buries MORE on apo peptide-free HLA
(1102-1239 A^2) than on the G12D complex (median 613 A^2) - is exactly the
question shape steering was meant to attack, so it is worth one read even
though no design is a candidate.

This changes **no gate and no threshold**.  It is a measurement on designs that
have already failed, chosen by a stated rule before the calls were made: the
two batch-08 designs with the highest median peptide share of buried area on
the G12D arm at r3np.  Nothing here can promote a design; `panel_results.json`
and the pass/fail verdicts are untouched.

apo is unconstrained (no chain C to point the pocket constraint at), exactly as
`results/validation/CONTROL_7OW6.md` did, so it is **not** like-for-like with
the constrained G12D/WT/decoy arms.  Same cell as the rest of batch-08: r3np.
"""
from __future__ import annotations

import argparse, json, sys, time
from pathlib import Path
from typing import Any, Dict, List

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(REPO / "scripts"))

import batch06_panel as p6  # noqa: E402
import run_batch04 as batch  # noqa: E402
import run_pipeline  # noqa: E402

OUT = REPO / "results" / "batch-08"
p6.OUT = OUT
batch.OUT = OUT
CONDITION = "nopot"   # r3np

# stated before the calls: the two most peptide-centred batch-08 designs
PROBE = ["B8-apex-46_bb46_donor-anchored-g2-dw",
         "B8-crest-46_bb46_donor-anchored-g2-dw"]


def sequences() -> Dict[str, str]:
    out: Dict[str, str] = {}
    for name in ("candidates.json", "candidates_w2.json"):
        for c in json.loads((OUT / name).read_text()):
            out[f"{c['backbone']}_{c['variant']}"] = c["sequence"]
    return out


def main(argv=None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--draws", type=int, default=3)
    parser.add_argument("--timeout", type=float, default=900.0)
    parser.add_argument("--max-boltz-calls", type=int, default=4)
    parser.add_argument("--max-retries", type=int, default=8)
    args = parser.parse_args(argv)

    (OUT / "structures").mkdir(parents=True, exist_ok=True)
    states = run_pipeline.load_target_states(REPO / "data" / "targets")
    client = p6.NvidiaClient(timeout=args.timeout, max_retries=args.max_retries)
    seqs = sequences()
    already, calls, started = p6.done_keys(), 0, time.time()
    summary: List[Dict[str, Any]] = []

    for name in PROBE:
        design = {"name": name, "sequence": seqs[name], "backbone": name.split("_donor")[0],
                  "variant": name.split("_", 2)[-1]}
        print(f"\n[probe] === {name} ===", flush=True)
        record: Dict[str, Any] = {"design": name, "sequence": design["sequence"]}
        for arm in ("apo", "decoy"):
            if calls >= args.max_boltz_calls:
                print("[stop] call cap reached"); break
            try:
                calls += p6.run_draws(client, design, arm, CONDITION, states,
                                      args.draws, already)
            except Exception as exc:  # noqa: BLE001
                print(f"  !! {arm} failed: {type(exc).__name__}: {exc}", flush=True)
                record[f"{arm}_error"] = f"{type(exc).__name__}: {exc}"
        rows = p6.load_rows()
        for arm in ("g12d", "apo", "decoy"):
            record[arm] = p6.arm_rows(rows, name, arm, CONDITION)
        summary.append(record)

    out = {"condition": CONDITION,
           "condition_alias": "r3np (recycling_steps 3, sampling_steps 50, "
                              "without_potentials True)",
           "note": "diagnostic only; apo is unconstrained and NOT like-for-like; "
                   "no gate, threshold or verdict is affected",
           "selection_rule": "the two batch-08 designs with the highest median "
                             "peptide share of buried area on the G12D arm",
           "boltz2_calls_this_segment": calls,
           "elapsed_seconds": round(time.time() - started, 1),
           "designs": summary}
    (OUT / "apo_probe.json").write_text(json.dumps(out, indent=2, default=str))
    print(f"\n[batch-08 apo probe] {calls} Boltz-2 calls -> {OUT / 'apo_probe.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
