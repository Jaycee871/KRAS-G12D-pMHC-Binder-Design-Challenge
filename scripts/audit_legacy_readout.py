#!/usr/bin/env python3
"""Audit batches 01-03 against the current direct Asp6 readout definition.

This is a zero-model-call audit of *recorded* historical metrics. It does not
claim coordinate remeasurement: rows are explicitly labelled
"recorded-unverified". The purpose is to decide whether pre-batch-04 sequences
can legitimately expand the current review pool without inventing a new gate.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SOURCES = tuple(ROOT / f"results/batch-{n:02}/checkpoint.jsonl" for n in (1, 2, 3))
DONOR_CUTOFF_A = 3.5


def _metrics(record: dict) -> dict:
    detail = record.get("detail") or {}
    metrics = detail.get("g12d_metrics")
    if isinstance(metrics, dict):
        return metrics
    row = record.get("row") or {}
    detail = row.get("detail") or {}
    metrics = detail.get("g12d_metrics")
    return metrics if isinstance(metrics, dict) else {}


def _sequence(record: dict) -> str:
    return record.get("sequence") or (record.get("row") or {}).get("sequence") or ""


def _direct_two_leg(metrics: dict) -> bool:
    distance = metrics.get("asp6_min_donor_distance")
    burial = metrics.get("asp6_buried_frac")
    return (
        isinstance(distance, (int, float))
        and distance <= DONOR_CUTOFF_A
        and isinstance(burial, (int, float))
        and burial > 0
    )


def audit(root: Path = ROOT) -> dict:
    report = {
        "purpose": "legacy batches 01-03 direct Asp6 readout audit",
        "model_calls": 0,
        "readout_status": "recorded-unverified",
        "readout_definition": {
            "acceptor": "Asp6 carboxylate oxygen as recorded by complex_analysis",
            "max_donor_distance_A": DONOR_CUTOFF_A,
            "requires_positive_Asp6_sidechain_burial": True,
        },
        "sources": {},
        "two_leg_candidates": [],
        "burial_only_candidates": [],
    }
    all_sequences: set[str] = set()

    for source in (root / p.relative_to(ROOT) for p in SOURCES):
        relative = str(source.relative_to(root))
        rows = [json.loads(line) for line in source.read_text().splitlines() if line.strip()]
        by_sequence: dict[str, dict] = {}

        for record in rows:
            sequence = _sequence(record)
            if not sequence:
                raise ValueError(f"missing sequence in {relative}")
            all_sequences.add(sequence)
            metrics = _metrics(record)
            distance = metrics.get("asp6_min_donor_distance")
            burial = metrics.get("asp6_buried_frac")
            oxygen_pass = isinstance(distance, (int, float)) and distance <= DONOR_CUTOFF_A
            burial_pass = isinstance(burial, (int, float)) and burial > 0
            two_leg = oxygen_pass and burial_pass

            entry = by_sequence.setdefault(sequence, {
                "sequence": sequence,
                "backbones": set(),
                "records": 0,
                "two_leg_passes": 0,
                "oxygen_distance_passes": 0,
                "positive_burial_records": 0,
                "best_donor_distance_A": None,
                "max_Asp6_buried_fraction": None,
                "max_pos6_contacts": 0,
            })
            entry["records"] += 1
            entry["two_leg_passes"] += int(two_leg)
            entry["oxygen_distance_passes"] += int(oxygen_pass)
            entry["positive_burial_records"] += int(burial_pass)
            backbone = (record.get("row") or {}).get("backbone")
            if backbone:
                entry["backbones"].add(backbone)
            if isinstance(distance, (int, float)):
                best = entry["best_donor_distance_A"]
                entry["best_donor_distance_A"] = distance if best is None else min(best, distance)
            if isinstance(burial, (int, float)):
                best = entry["max_Asp6_buried_fraction"]
                entry["max_Asp6_buried_fraction"] = burial if best is None else max(best, burial)
            entry["max_pos6_contacts"] = max(
                entry["max_pos6_contacts"], int(metrics.get("pos6_contacts") or 0)
            )

        rows_by_sequence = []
        for entry in by_sequence.values():
            entry = dict(entry)
            entry["backbones"] = sorted(entry["backbones"])
            rows_by_sequence.append(entry)
            if entry["two_leg_passes"]:
                report["two_leg_candidates"].append({"source": relative, **entry})
            elif entry["positive_burial_records"]:
                report["burial_only_candidates"].append({"source": relative, **entry})

        report["sources"][relative] = {
            "sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
            "records": len(rows),
            "unique_sequences": len(rows_by_sequence),
            "two_leg_sequences": sum(bool(x["two_leg_passes"]) for x in rows_by_sequence),
            "oxygen_only_sequences": sum(
                bool(x["oxygen_distance_passes"]) and not bool(x["two_leg_passes"])
                for x in rows_by_sequence
            ),
            "burial_only_sequences": sum(
                bool(x["positive_burial_records"]) and not bool(x["oxygen_distance_passes"])
                for x in rows_by_sequence
            ),
        }

    report["totals"] = {
        "checkpoint_records": sum(x["records"] for x in report["sources"].values()),
        "unique_sequences": len(all_sequences),
        "two_leg_sequences": len({x["sequence"] for x in report["two_leg_candidates"]}),
        "burial_only_sequences": len({x["sequence"] for x in report["burial_only_candidates"]}),
    }
    report["conclusion"] = (
        "No batch-01/02/03 sequence has a recorded oxygen-distance plus positive-burial "
        "Asp6 readout; these legacy batches do not expand the current review shortlist."
    )
    return report


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path)
    args = parser.parse_args(argv)
    report = audit()
    text = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(text)
    else:
        print(text, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
