#!/usr/bin/env python3
"""Apply the Asp6 read-out gate to batch-04 and rank the survivors.

Ranking uses measured geometry only, in the order agreed after both the ipTM
`specificity_gap` and the paired MM-GBSA `ddG_sel` failed the 7OW6 positive
control (see docs/PIPELINE.md section 6 and results/ddg-scoring/FINDINGS.md):

    1. peptide_buried_frac (floor 0.40)
    2. contacts to the Asp6 **side chain** (backbone contacts excluded)
    3. net interface charge out of strongly negative territory
    4. topological diversity (sequence-identity clustering at 70 %)

ipTM is reported as a pose-quality read-out only: never a gap, never a
tie-breaker.
"""
from __future__ import annotations

import csv, json, statistics, sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))
import run_batch04  # noqa: E402  (column list and gate constants)

OUT = REPO / "results" / "batch-04"
PEPTIDE_BURIED_FLOOR = 0.40
CHARGE_FLOOR = -2.0


def identity(a: str, b: str) -> float:
    n = min(len(a), len(b))
    if not n:
        return 0.0
    return sum(1 for x, y in zip(a[:n], b[:n]) if x == y) / n


def main() -> int:
    rows = [json.loads(line)["row"]
            for line in (OUT / "checkpoint.jsonl").read_text().splitlines()
            if line.strip()]
    # last write wins per sequence
    by_sequence = {r["sequence"]: r for r in rows}
    rows = list(by_sequence.values())

    survivors = [r for r in rows if r.get("gate_pass")]
    survivors.sort(key=lambda r: (
        -(float(r.get("peptide_buried_frac") or 0)),
        -(int(r.get("pos6_sidechain_contacts") or 0)),
        -(float(r.get("net_interface_charge") or 0)),
    ))

    # topological diversity: cluster survivors at 70 % identity, keep the
    # representative order but flag redundancy.
    clusters: list[list[dict]] = []
    for row in survivors:
        for cluster in clusters:
            if identity(cluster[0]["sequence"], row["sequence"]) >= 0.70:
                cluster.append(row)
                break
        else:
            clusters.append([row])
    for index, cluster in enumerate(clusters, start=1):
        for member in cluster:
            member["cluster"] = index

    ranked = survivors
    for rank, row in enumerate(ranked, start=1):
        row["rank"] = rank
        row["name"] = f"KRASG12D_B04_{rank:03d}_{row['backbone']}"

    charges = [float(r["net_interface_charge"]) for r in rows
               if r.get("net_interface_charge") not in ("", None)]
    donor_distances = [float(r["asp6_min_donor_distance"]) for r in rows
                       if r.get("asp6_min_donor_distance") not in ("", None)]
    summary = {
        "scored": len(rows),
        "gate_survivors": len(survivors),
        "clusters": len(clusters),
        "net_interface_charge": {
            "median": round(statistics.median(charges), 2) if charges else None,
            "mean": round(statistics.fmean(charges), 2) if charges else None,
            "min": min(charges) if charges else None,
            "max": max(charges) if charges else None,
            "n_negative": sum(1 for c in charges if c < 0),
            "n_zero_or_positive": sum(1 for c in charges if c >= 0),
            "batch02_median": -1.0, "batch03_median": -1.0,
        },
        "asp6_donor_distances": sorted(donor_distances),
        "designs_with_any_donor_in_range": len(donor_distances),
        "peptide_buried_frac": {
            "max": max((float(r["peptide_buried_frac"]) for r in rows
                        if r.get("peptide_buried_frac") not in ("", None)),
                       default=None),
            "n_above_floor": sum(1 for r in rows
                                 if float(r.get("peptide_buried_frac") or 0)
                                 >= PEPTIDE_BURIED_FLOOR),
            "floor": PEPTIDE_BURIED_FLOOR,
        },
        "survivors_meeting_peptide_floor": sum(
            1 for r in survivors
            if float(r.get("peptide_buried_frac") or 0) >= PEPTIDE_BURIED_FLOOR),
        "survivors_above_charge_floor": sum(
            1 for r in survivors
            if float(r.get("net_interface_charge") or 0) >= CHARGE_FLOOR),
    }
    (OUT / "gate_results.json").write_text(json.dumps(
        {"summary": summary, "survivors": survivors}, indent=2, default=str))

    columns = ["cluster"] + run_batch04.CSV_COLUMNS
    all_rows = sorted(rows, key=lambda r: (
        not r.get("gate_pass"),
        -(float(r.get("peptide_buried_frac") or 0)),
        -(int(r.get("pos6_sidechain_contacts") or 0))))
    with (OUT / "designs.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        for row in all_rows:
            writer.writerow({c: row.get(c, "") for c in columns})

    print(json.dumps(summary, indent=2))
    for row in ranked:
        print(f"#{row['rank']} {row['backbone']} ({row['variant']}) "
              f"cluster={row['cluster']} "
              f"pep_buried={row['peptide_buried_frac']} "
              f"pos6_sc={row['pos6_sidechain_contacts']} "
              f"donor={row['asp6_min_donor_distance']} A "
              f"({row['asp6_donor_partner']}) "
              f"asp6_buried={row['asp6_buried_frac']} "
              f"charge={row['net_interface_charge']} "
              f"area_pep={row['buried_area_peptide']} iptm={row['iptm_g12d']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
