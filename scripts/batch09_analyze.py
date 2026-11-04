#!/usr/bin/env python3
"""Analyse batch-09's matched regional-anchor counter-screen.

No new gate is introduced here. The analysis reports medians and simple,
control-free comparisons such as whether apo buries more total surface than
G12D under the SAME regional docking prior.
"""
from __future__ import annotations

import argparse
import json
import math
import statistics
import sys
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional

REPO = Path(__file__).resolve().parents[1]
OUT = REPO / "results" / "batch-09"
CHECKPOINT = OUT / "matched_anchor_checkpoint.jsonl"
sys.path.insert(0, str(REPO / "scripts"))
import batch09_matched_anchor as experiment  # noqa: E402


def median(values: Iterable[Any]) -> Optional[float]:
    clean: List[float] = []
    for value in values:
        if value is None or value == "":
            continue
        try:
            clean.append(float(value))
        except (TypeError, ValueError):
            continue
    return statistics.median(clean) if clean else None


def rows(historical=False) -> List[Dict[str, Any]]:
    if not CHECKPOINT.exists():
        raise FileNotFoundError(
            f"{CHECKPOINT} not found; run scripts/batch09_matched_anchor.py first"
        )
    data = [json.loads(line) for line in CHECKPOINT.read_text().splitlines()
            if line.strip()]
    validate_panel(data, historical=historical)
    return data


def validate_panel(data: List[Dict[str, Any]], historical=False) -> None:
    """Data-integrity checks, not biological pass/fail thresholds."""
    expected = {(key, arm, draw) for key in experiment.SUBJECT_KEYS
                for arm in experiment.ARMS for draw in range(3)}
    actual = [(r["design"], r["arm"], r["draw"]) for r in data]
    if len(set(actual)) != len(actual):
        raise ValueError("duplicate samples in matched-anchor panel")
    if set(actual) != expected:
        raise ValueError("matched-anchor panel is incomplete or contains unexpected samples; "
                         "expected 3 designs x 4 arms x 3 draws")
    for record in data:
        row = record.get("row") or {}
        if not historical:
            if row.get("readout_version") != experiment.panel.READOUT_VERSION:
                raise ValueError("legacy readout flags: run scripts/audit_carboxylate_readout.py; "
                                 "use --historical-recorded-metrics only to reproduce the old summary")
            if record["arm"] == "g12d" and row.get("asp6_buried_frac") is None:
                raise ValueError("corrected G12D readout requires Asp6 burial")
        if row.get("analysis_error"):
            raise ValueError("matched-anchor panel contains an analysis error")
        metrics = ["buried_area_total"]
        if record["arm"] != "apo":
            metrics += ["buried_area_peptide", "peptide_buried_frac", "pos6_contacts"]
        for metric in metrics:
            try:
                valid = math.isfinite(float(row[metric]))
            except (KeyError, TypeError, ValueError):
                valid = False
            if not valid:
                raise ValueError(f"missing/non-finite matched-anchor metric: {metric}")
    original_out = experiment.OUT
    try:
        experiment.OUT = OUT
        experiment.validate_records(data, experiment.load_subjects(), 3)
    finally:
        experiment.OUT = original_out


def arm_summary(records: List[Dict[str, Any]]) -> Dict[str, Any]:
    rs = [r.get("row") or {} for r in records]
    return {
        "n": len(rs),
        "median_total_buried_area": median(
            r.get("buried_area_total") for r in rs),
        "median_peptide_buried_area": median(
            r.get("buried_area_peptide") for r in rs),
        "median_peptide_buried_fraction": median(
            r.get("peptide_buried_frac") for r in rs),
        "median_pos6_contacts": median(r.get("pos6_contacts") for r in rs),
        "median_pos6_sidechain_contacts": median(
            r.get("pos6_sidechain_contacts") for r in rs),
        "asp6_gate_passes": sum(
            1 for r in rs if str(r.get("asp6_gate_pass")).lower() == "true"),
        "median_iptm": median(r.get("iptm") for r in rs),
        "analysis_errors": [
            r.get("analysis_error") for r in rs if r.get("analysis_error")
        ],
    }


def safe_ratio(numerator: Optional[float],
               denominator: Optional[float]) -> Optional[float]:
    if numerator is None or denominator in (None, 0):
        return None
    return numerator / denominator


def main(argv=()) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--historical-recorded-metrics', action='store_true',
                        help='rebuild legacy counts into a separate explicitly historical file')
    args = parser.parse_args(argv)
    data = rows(historical=args.historical_recorded_metrics)
    designs = sorted({r["design"] for r in data})
    result: Dict[str, Any] = {
        "purpose": "matched HLA-regional-anchor diagnostic only",
        "historical_recorded_metrics": args.historical_recorded_metrics,
        "readout_version": "legacy-recorded-flags" if args.historical_recorded_metrics else experiment.panel.READOUT_VERSION,
        "new_hard_gate": False,
        "panel_complete": True,
        "expected_draws_per_arm": 3,
        "designs": {},
    }

    for design in designs:
        per_arm: Dict[str, Any] = {}
        for arm in ("g12d", "wt", "apo", "decoy"):
            subset = [r for r in data
                      if r["design"] == design and r["arm"] == arm]
            if subset:
                per_arm[arm] = arm_summary(subset)

        g = per_arm.get("g12d", {})
        a = per_arm.get("apo", {})
        w = per_arm.get("wt", {})
        d = per_arm.get("decoy", {})

        g_total = g.get("median_total_buried_area")
        a_total = a.get("median_total_buried_area")
        g_pep = g.get("median_peptide_buried_area")
        w_pep = w.get("median_peptide_buried_area")
        d_pep = d.get("median_peptide_buried_area")

        result["designs"][design] = {
            "arms": per_arm,
            "matched_comparisons": {
                "apo_to_g12d_total_area_ratio": safe_ratio(a_total, g_total),
                "apo_buries_more_total_area_than_g12d": (
                    None if a_total is None or g_total is None
                    else a_total > g_total
                ),
                "wt_to_g12d_peptide_area_ratio_diagnostic_only":
                    safe_ratio(w_pep, g_pep),
                "decoy_to_g12d_peptide_area_ratio":
                    safe_ratio(d_pep, g_pep),
            },
        }

    output_name = ("historical_recorded_analysis_rebuilt.json" if args.historical_recorded_metrics
                   else "matched_anchor_analysis.json")
    (OUT / output_name).write_text(
        json.dumps(result, indent=2, default=str)
    )

    for design, info in result["designs"].items():
        comp = info["matched_comparisons"]
        print(
            f"{design}: apo/G12D total="
            f"{comp['apo_to_g12d_total_area_ratio']} "
            f"(apo>G12D={comp['apo_buries_more_total_area_than_g12d']}), "
            f"WT/G12D peptide="
            f"{comp['wt_to_g12d_peptide_area_ratio_diagnostic_only']}, "
            f"decoy/G12D peptide="
            f"{comp['decoy_to_g12d_peptide_area_ratio']}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
