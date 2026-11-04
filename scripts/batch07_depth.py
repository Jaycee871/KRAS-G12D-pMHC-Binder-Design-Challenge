#!/usr/bin/env python3
"""Batch-07 stage 3c: how reproducible is the WT verdict itself, and apo/decoy.

Two independent 3-draw WT calls on bb54-g2 under the *same* cell (`r3np`) gave
opposite verdicts under the unretuned batch-06 criterion (-27.9 %/-25 % FAIL
vs -33.8 %/-50 % PASS). The criterion is not touched; instead more draws are
taken so the reported number has an error bar: two further G12D calls and two
further WT calls at `r3np`, giving n = 9 per state, reported both per-call (the
criterion as defined) and pooled (as context).

Because one WT call did clear, the apo and decoy arms are run for bb54-g2 under
`r3np` as the brief requires. The apo arm has no chain C, so it is run
UNCONSTRAINED, exactly as CONTROL_7OW6.md did, and is NOT like-for-like.
"""
from __future__ import annotations

import json, sys, time
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(REPO / "scripts"))

import batch06_panel as P  # noqa: E402
import batch07_settings as S  # noqa: E402
import batch07_panel as B7  # noqa: E402
import run_batch04 as batch  # noqa: E402
import run_pipeline  # noqa: E402
import sampling_study as study  # noqa: E402
from nvidia_client import NvidiaClient  # noqa: E402

OUT = REPO / "results" / "batch-07"
P.OUT = OUT
batch.OUT = OUT
for label in ("r3np_rep3", "r3np_rep4", "r3np_rep5"):
    study.CONDITIONS[label] = study.CONDITIONS["r3np"]


def main() -> int:
    states = run_pipeline.load_target_states(REPO / "data" / "targets")
    client = NvidiaClient(timeout=900.0)
    parent = B7.parent_design()
    jobs = [("g12d", "r3np_rep3"), ("wt", "r3np_rep3"),
            ("g12d", "r3np_rep4"), ("wt", "r3np_rep4"),
            ("apo", "r3np"), ("decoy", "r3np")]
    already, calls, started = P.done_keys(), 0, time.time()
    for arm, cond in jobs:
        print(f"\n[depth] {arm} @ {cond}", flush=True)
        calls += P.run_draws(client, parent, arm, cond, states, 3, already)

    rows = P.load_rows()
    name = parent["name"]
    def ar(arm, conds):
        return [r for c in conds for r in P.arm_rows(rows, name, arm, c)]
    g_cells = ["r3np", "r3np_rep3", "r3np_rep4"]
    w_cells = ["r3np", "r3np_rep2", "r3np_rep3", "r3np_rep4"]
    per_call = {}
    for c in w_cells:
        gc = c if c in g_cells else "r3np"
        per_call[c] = P.wt_verdict(P.arm_rows(rows, name, "g12d", gc),
                                   P.arm_rows(rows, name, "wt", c))
    pooled = P.wt_verdict(ar("g12d", g_cells), ar("wt", w_cells))
    out = {"design": name, "cell": "r3np", "boltz2_calls_this_segment": calls,
           "elapsed_seconds": round(time.time() - started, 1),
           "wt_criterion": P.WT_CRITERION,
           "per_call_wt_verdicts": per_call,
           "pooled_wt_verdict_context_only": pooled,
           "apo_rows": P.arm_rows(rows, name, "apo", "r3np"),
           "decoy_rows": P.arm_rows(rows, name, "decoy", "r3np")}
    (OUT / "depth_results.json").write_text(json.dumps(out, indent=2, default=str))
    for c, v in per_call.items():
        print(f"  WT {c}: area {v['peptide_area_loss']} pos6 {v['pos6_contact_loss']} "
              f"PASS={v['wt_gate_pass']}")
    print(f"  pooled (n=9/12): area {pooled['peptide_area_loss']} "
          f"pos6 {pooled['pos6_contact_loss']} PASS={pooled['wt_gate_pass']}")
    print(f"\n[depth] {calls} Boltz-2 calls")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
