#!/usr/bin/env python3
"""Batch-06 WT-arm DIAGNOSTIC for the two 2/3 near-misses.

This is **not** part of the gate and loosens nothing. No design cleared the
3-of-3 Asp6 read-out gate, so by `scripts/batch06_panel.py` no design was
entitled to a WT arm. This script runs the WT arm anyway, on the two designs
that reached 2 of 3, purely to answer the batch-06 scientific question --
did shrinking and re-centring the paratope buy ANY wild-type discrimination?
-- and the answer is reported as a diagnostic, never as a pass.

Same WT criterion as batch06_panel.WT_CRITERION, same sampling setting, same
pocket constraint, same 3-draws-in-one-call scheme.
"""
import json, sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts")); sys.path.insert(0, str(REPO / "src"))

import batch06_panel as P
import run_pipeline
from nvidia_client import NvidiaClient

NEAR_MISSES = [("B6-c6-46_w2_bb54_donor-anchored-g2", "recyc8"),
               ("B6-c6-46_w2_bb36_donor-anchored-g2", "recyc8")]


def main() -> int:
    states = run_pipeline.load_target_states(REPO / "data" / "targets")
    client = NvidiaClient(timeout=900.0)
    rows = P.load_rows()
    already = P.done_keys()
    out = []
    for name, cond in NEAR_MISSES:
        src = next(r for r in rows if r["design"] == name and r["arm"] == "g12d"
                   and r["condition"] == cond)
        design = {"name": name, "sequence": src["sequence"],
                  "backbone": src.get("backbone", ""), "variant": src.get("variant", "")}
        print(f"\n[wt-diagnostic] {name} @ {cond}", flush=True)
        P.run_draws(client, design, "wt", cond, states, 3, already)
        rows2 = P.load_rows()
        verdict = P.wt_verdict(P.arm_rows(rows2, name, "g12d", cond),
                               P.arm_rows(rows2, name, "wt", cond))
        verdict["design"], verdict["condition"] = name, cond
        verdict["status"] = ("DIAGNOSTIC ONLY - this design cleared the Asp6 "
                             "read-out gate in only 2 of 3 draws and is NOT a "
                             "batch-06 survivor regardless of this result")
        print(f"  -> peptide-area loss {verdict['peptide_area_loss']}, "
              f"pos6-contact loss {verdict['pos6_contact_loss']}, "
              f"would-pass={verdict['wt_gate_pass']}", flush=True)
        out.append(verdict)
    (P.OUT / "wt_diagnostic.json").write_text(json.dumps(
        {"criterion": P.WT_CRITERION, "designs": out}, indent=2, default=str))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
