#!/usr/bin/env python3
"""Batch-08 stage 7: de-wet the helix face of the *working* design, in place.

The shape-steered backbones of stages 1-3 never reached the Asp6 gate 3-of-3,
so the shape hypothesis is still untested on a design that works.  This is the
controlled version of the same test: take bb54-g2 - the one design that reaches
the gate reproducibly at r3np (8/9 pooled in batch-07) - and change nothing
except the chemistry it points at the alpha1/alpha2 helices.

Helix-facing positions are read off bb54-g2's own *predicted* r3np poses (a
binder residue whose closest heavy-atom approach to an alpha1/alpha2 helix is
under 6 A and closer than its approach to the peptide), the six positions
nearest Asp6 are protected, and the large hydrophobics among the rest become
serine.  Six positions change; the whole donor shell is untouched, the fold is
untouched, the length is untouched.  Parent and twin are run side by side in
the same cell in the same segment, so the comparison is paired.

Gate, criterion and arms are `scripts/batch06_panel.py` verbatim, as in
`scripts/batch08_panel.py`.
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

PARENT = "MMLPDDAIAELSGTPEARELARRLGGDERVAINYRAGRRGRVGVEA"
DEWET = "SSSPDDAIAELSGTPEARELARRLGGDERSAINYRAGRRGRSGSEA"
EDITS = [(1, "M"), (2, "M"), (3, "L"), (30, "V"), (42, "V"), (44, "V")]

DESIGNS: List[Dict[str, Any]] = [
    {"name": "B8-bb54-g2-parent", "sequence": PARENT, "variant": "parent-rerun",
     "backbone": "B6-c6-46_w2_bb54", "binder_length": 46},
    {"name": "B8-bb54-g2-dewet", "sequence": DEWET, "variant": "helix-face-dewet",
     "backbone": "B6-c6-46_w2_bb54", "binder_length": 46,
     "dewet_positions": [p for p, _ in EDITS]},
]


def main(argv=None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--draws", type=int, default=3)
    parser.add_argument("--timeout", type=float, default=900.0)
    parser.add_argument("--max-boltz-calls", type=int, default=10)
    parser.add_argument("--max-retries", type=int, default=8)
    args = parser.parse_args(argv)

    assert len(PARENT) == len(DEWET) == 46
    assert all(PARENT[p - 1] == aa for p, aa in EDITS)

    (OUT / "structures").mkdir(parents=True, exist_ok=True)
    states = run_pipeline.load_target_states(REPO / "data" / "targets")
    client = p6.NvidiaClient(timeout=args.timeout, max_retries=args.max_retries)
    already, calls, started = p6.done_keys(), 0, time.time()
    summary = []

    for design in DESIGNS:
        name = design["name"]
        print(f"\n[dewet] === {name} ===", flush=True)
        record: Dict[str, Any] = {k: design.get(k) for k in
                                  ("name", "sequence", "variant", "backbone",
                                   "binder_length", "dewet_positions")}
        if calls >= args.max_boltz_calls:
            break
        try:
            calls += p6.run_draws(client, design, "g12d", CONDITION, states,
                                  args.draws, already)
        except Exception as exc:  # noqa: BLE001
            record["error"] = f"{type(exc).__name__}: {exc}"
            summary.append(record); continue
        rows = p6.arm_rows(p6.load_rows(), name, "g12d", CONDITION)
        passes = sum(1 for r in rows if str(r.get("asp6_gate_pass")).lower() == "true")
        record["g12d"] = {"gate_passes": passes, "draws": len(rows), "rows": rows}
        print(f"  -> {CONDITION}: Asp6 gate {passes}/{len(rows)}", flush=True)
        record["gate_3of3"] = bool(passes >= args.draws)
        if not record["gate_3of3"]:
            summary.append(record); continue
        if calls < args.max_boltz_calls:
            calls += p6.run_draws(client, design, "wt", CONDITION, states,
                                  args.draws, already)
        all_rows = p6.load_rows()
        record["wt"] = p6.wt_verdict(p6.arm_rows(all_rows, name, "g12d", CONDITION),
                                     p6.arm_rows(all_rows, name, "wt", CONDITION))
        print(f"  -> WT arm: area loss {record['wt']['peptide_area_loss']}, "
              f"pos6 loss {record['wt']['pos6_contact_loss']}, "
              f"PASS={record['wt']['wt_gate_pass']}", flush=True)
        if record["wt"]["wt_gate_pass"]:
            for arm in ("apo", "decoy"):
                if calls >= args.max_boltz_calls:
                    break
                try:
                    calls += p6.run_draws(client, design, arm, CONDITION, states,
                                          args.draws, already)
                except Exception as exc:  # noqa: BLE001
                    print(f"  !! {arm} failed: {exc}", flush=True)
            all_rows = p6.load_rows()
            record["apo"] = p6.arm_rows(all_rows, name, "apo", CONDITION)
            record["decoy"] = p6.arm_rows(all_rows, name, "decoy", CONDITION)
        summary.append(record)

    out = {"condition": CONDITION, "wt_criterion": p6.WT_CRITERION,
           "edits": EDITS, "parent": PARENT, "dewet": DEWET,
           "boltz2_calls_this_segment": calls,
           "elapsed_seconds": round(time.time() - started, 1),
           "designs": summary}
    (OUT / "dewet_results.json").write_text(json.dumps(out, indent=2, default=str))
    print(f"\n[batch-08 dewet] {calls} Boltz-2 calls -> {OUT / 'dewet_results.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
