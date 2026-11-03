#!/usr/bin/env python3
"""Batch-07 analysis: pose spreads, buried-area tables, distinct-draw counts.

Re-uses `scripts/batch06_analyze.pose_spread` (and through it
`scripts/replicate_rmsd_batch04.binder_rmsd`) verbatim, with `OUT` pointed at
`results/batch-07`: binder C-alpha RMSD in the **pMHC frame** with the binder
NOT re-fitted, and the frame's own RMSD reported alongside as the artefact
floor. In the apo arm Boltz-2 relabels the binder to chain C and there is no
peptide, so the frame is HLA chain A only - stated, not silently substituted.
No gate or threshold is changed; `src/design_gates.py` is untouched.
"""
from __future__ import annotations

import json, statistics, sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(REPO / "scripts"))

import batch06_analyze as A  # noqa: E402
import batch06_panel as P  # noqa: E402

OUT = REPO / "results" / "batch-07"
A.OUT = OUT
P.OUT = OUT


def num(v):
    return A.num(v)


def main() -> int:
    rows = A.load_rows(OUT / "panel_checkpoint.jsonl")
    screen = [json.loads(l) for l in
              (OUT / "screen_summary.jsonl").read_text().splitlines() if l.strip()]
    cells = sorted({(r["design"], r["arm"], r["condition"]) for r in rows})
    out = {"n_structures": len(rows), "n_cells": len(cells),
           "boltz2_calls_total": len(cells), "cells": {}}
    for design, arm, cond in cells:
        members = [r for r in rows if r["design"] == design and r["arm"] == arm
                   and r["condition"] == cond]
        vals = lambda k: [num(m["row"].get(k)) for m in members]  # noqa: E731
        def stat(k):
            xs = [v for v in vals(k) if v is not None]
            return {"values": xs,
                    "median": round(statistics.median(xs), 3) if xs else None}
        entry = {
            "draws": len(members),
            "distinct_structures": len({m["structure_sha1"] for m in members}),
            "asp6_gate_passes": sum(1 for m in members if str(
                m["row"].get("asp6_gate_pass")).lower() == "true"),
            "buried_area_peptide": stat("buried_area_peptide"),
            "buried_area_hla": stat("buried_area_hla"),
            "peptide_buried_frac": stat("peptide_buried_frac"),
            "pos6_contacts": stat("pos6_contacts"),
            "pose_spread": A.pose_spread(design, arm, cond, rows),
        }
        pbf = [v for v in vals("peptide_buried_frac") if v is not None]
        entry["peptide_buried_frac_ge_0.40"] = sum(1 for v in pbf if v >= 0.40)
        entry["peptide_buried_frac_ge_0.18"] = sum(1 for v in pbf if v >= 0.18)
        out["cells"][f"{design}|{arm}|{cond}"] = entry

    # the parent's independent WT repeat shares r3np's Boltz-2 kwargs, so it is
    # compared against the r3np G12D draws (same settings, different call).
    def arm_rows(design, arm, cond):
        return [r["row"] for r in rows if r["design"] == design
                and r["arm"] == arm and r["condition"] == cond]
    parent = "B6-c6-46_w2_bb54_donor-anchored-g2"
    out["parent_wt_repeat_vs_r3np_g12d"] = P.wt_verdict(
        arm_rows(parent, "g12d", "r3np"), arm_rows(parent, "wt", "r3np_rep2"))
    out["screen"] = {
        "designs": len(screen),
        "gate_3of3": [s["design"] for s in screen if s.get("gate_3of3")],
        "peptide_buried_frac_all_draws": None,
    }
    all_pbf = [num(r["row"].get("peptide_buried_frac")) for r in rows
               if r["arm"] == "g12d"]
    all_pbf = [v for v in all_pbf if v is not None]
    out["screen"]["peptide_buried_frac_all_draws"] = {
        "n": len(all_pbf), "median": round(statistics.median(all_pbf), 3),
        "max": round(max(all_pbf), 3),
        "ge_0.40": sum(1 for v in all_pbf if v >= 0.40),
        "ge_0.18": sum(1 for v in all_pbf if v >= 0.18)}
    (OUT / "analysis.json").write_text(json.dumps(out, indent=2, default=str))
    print(json.dumps({k: out[k] for k in
                      ("n_structures", "n_cells", "screen",
                       "parent_wt_repeat_vs_r3np_g12d")}, indent=2)[:3000])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
