#!/usr/bin/env python3
"""Batch-07 stage 3b: is the gate condition the same as the WT-discriminating one?

Stage 3 found both 3-of-3 designs fail the (verbatim) batch-06 WT criterion
under the very sampling cell that gave them 3-of-3. Batch-06, however, saw
bb54-g2 discriminate strongly under `recyc8` - where it was only 2-of-3. This
probe tests the obvious confound directly, changing no gate and no criterion:

  * run the G12D and WT arms of `bb54_rigid-gp-p2` under `recyc8`, the cell
    where bb54-g2 discriminated, and
  * repeat the WT arm of each 3-of-3 design once more, independently, under its
    own gate cell, to measure how noisy the WT verdict itself is.

A design still only passes if BOTH arms clear in the SAME cell: Asp6 gate 3/3
and WT loss by `batch06_panel.WT_CRITERION`, unretuned.
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

# independent repeats of an already-run arm get their own condition label so
# the checkpoint key stays unique; same Boltz-2 kwargs underneath.
study.CONDITIONS["r3_rep2"] = study.CONDITIONS["r3"]
study.CONDITIONS["r3np_rep2"] = study.CONDITIONS["r3np"]


def main() -> int:
    (OUT / "structures").mkdir(parents=True, exist_ok=True)
    states = run_pipeline.load_target_states(REPO / "data" / "targets")
    client = NvidiaClient(timeout=900.0)
    designs = {d["name"]: d for d in B7.screened_designs()}
    designs["B6-c6-46_w2_bb54_donor-anchored-g2"] = B7.parent_design()
    rigid = designs["bb54_rigid-gp-p2"]
    parent = designs["B6-c6-46_w2_bb54_donor-anchored-g2"]

    jobs = [
        (rigid, "g12d", "recyc8"), (rigid, "wt", "recyc8"),
        (rigid, "g12d", "r3_rep2"), (rigid, "wt", "r3_rep2"),
        (parent, "wt", "r3np_rep2"),
    ]
    already, calls, started = P.done_keys(), 0, time.time()
    for design, arm, cond in jobs:
        print(f"\n[cross] {design['name']} {arm} @ {cond}", flush=True)
        calls += P.run_draws(client, design, arm, cond, states, 3, already)

    rows = P.load_rows()
    out = {}
    for design in (rigid, parent):
        name = design["name"]
        for cond in ("recyc8", "r3_rep2", "r3np_rep2", design["condition"]):
            g = P.arm_rows(rows, name, "g12d", cond)
            w = P.arm_rows(rows, name, "wt", cond)
            if not g and not w:
                continue
            entry = {"g12d_draws": len(g),
                     "g12d_gate_passes": sum(1 for r in g if str(
                         r.get("asp6_gate_pass")).lower() == "true"),
                     "wt_draws": len(w)}
            if g and w:
                entry["wt"] = P.wt_verdict(g, w)
            out[f"{name}@{cond}"] = entry
            print(f"  {name}@{cond}: gate {entry['g12d_gate_passes']}/"
                  f"{entry['g12d_draws']}  wt={entry.get('wt', {}).get('wt_gate_pass')}"
                  f" area_loss={entry.get('wt', {}).get('peptide_area_loss')}"
                  f" pos6_loss={entry.get('wt', {}).get('pos6_contact_loss')}")
    (OUT / "cross_condition.json").write_text(json.dumps(
        {"boltz2_calls_this_segment": calls,
         "elapsed_seconds": round(time.time() - started, 1),
         "wt_criterion": P.WT_CRITERION, "cells": out}, indent=2, default=str))
    print(f"\n[cross] {calls} Boltz-2 calls")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
