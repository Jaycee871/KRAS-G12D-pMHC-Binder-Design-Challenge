#!/usr/bin/env python3
"""Build a diversified Track-1 submission portfolio from the frozen evidence pool.

This is a submission-allocation layer, not a new biological gate. It never
changes historical batch outcomes. The goal is to avoid spending several of
Track 1's top-20 wet-lab slots on highly correlated or already-known failure
modes.

Policy:
- ranks 1-3: the three corrected matched-panel leads (exploit)
- ranks 4-20: one sequence per new family, no >2-Cys burden, no repo sequence
  liability, no exact/family reuse from the batch-10 failed set; choose the
  most sequence-diverse candidate greedily (explore)
- ranks 21-30: reserve designs, allowing known failed families but still
  enforcing clean sequence liabilities and at most two sequences per family
  across the full file
"""
from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path
from typing import Dict, List, Mapping, Sequence

ROOT = Path(__file__).resolve().parents[1]
POOL = ROOT / "results" / "readout-correction" / "submission-review" / "evidence_pool.json"
OUTDIR = ROOT / "results" / "final-submission"

import sys
sys.path.insert(0, str(ROOT / "src"))
from run_pipeline import sequence_filter_failures  # noqa: E402
from validate_submission import validate_file  # noqa: E402

LEAD_ALIASES = (
    "B5-c6-62_bb50_donor-anchored-g1",
    "B4-c6-62_bb25_donor-anchored",
    "B6-c6-46_w2_bb54_donor-anchored-g2",
)

BATCH10_FAILED_FAMILIES = {
    "B6-c6-46_w2_bb36",
    "B5-c6-62_bb173",
    "bb54_parent_and_derivatives",
    "B5-c6-62_bb90",
    "B4-c67-65_bb26",
    "B6-c6-46_w2_bb38",
}

NO_DIRECT_REASON = "no recorded successful direct Asp6 readout"
TOP_N = 20
TOTAL_N = 30
MAX_CYS = 2


def _alias(candidate: Mapping) -> str:
    aliases = candidate.get("aliases") or []
    return str(aliases[0]) if aliases else ""


def _protocol_stats(candidate: Mapping) -> Dict[str, int]:
    g_records = g_success = coord_records = coord_success = 0
    for protocol in candidate.get("protocols") or []:
        if protocol.get("arm") != "g12d":
            continue
        g_records += int(protocol.get("records") or 0)
        g_success += int(protocol.get("success_records") or 0)
        if protocol.get("readout_status") == "coordinate-remeasured":
            coord_records += int(protocol.get("records") or 0)
            coord_success += int(protocol.get("success_records") or 0)
    return {
        "g12d_records": g_records,
        "g12d_success_records": g_success,
        "coord_g12d_records": coord_records,
        "coord_g12d_success_records": coord_success,
    }


def _feature(candidate: Mapping) -> Dict:
    seq = str(candidate["sequence"])
    return {
        "candidate": dict(candidate),
        "sequence": seq,
        "alias": _alias(candidate),
        "family": str(candidate.get("family") or ""),
        "length": len(seq),
        "cys_count": seq.count("C"),
        "sequence_filter_failures": sequence_filter_failures(seq),
        **_protocol_stats(candidate),
    }


def _levenshtein(a: str, b: str) -> int:
    previous = list(range(len(b) + 1))
    for i, ca in enumerate(a, start=1):
        current = [i]
        for j, cb in enumerate(b, start=1):
            current.append(
                min(
                    previous[j] + 1,
                    current[j - 1] + 1,
                    previous[j - 1] + (ca != cb),
                )
            )
        previous = current
    return previous[-1]


def normalized_edit_similarity(a: str, b: str) -> float:
    return 1.0 - _levenshtein(a, b) / max(len(a), len(b))


def min_distance_to_selected(item: Mapping, selected: Sequence[Mapping]) -> float:
    return min(
        1.0 - normalized_edit_similarity(item["sequence"], other["sequence"])
        for other in selected
    )


def _clean_for_portfolio(item: Mapping) -> bool:
    return (
        40 <= int(item["length"]) <= 100
        and int(item["cys_count"]) <= MAX_CYS
        and not item["sequence_filter_failures"]
    )


def select_portfolio(pool: Mapping) -> List[Dict]:
    features = [_feature(c) for c in pool["candidates"]]

    fixed: List[Dict] = []
    for alias in LEAD_ALIASES:
        matches = [
            item for item in features
            if alias in (item["candidate"].get("aliases") or [])
        ]
        if len(matches) != 1:
            raise RuntimeError(f"expected one lead for {alias}, got {len(matches)}")
        fixed.append(matches[0])

    selected = list(fixed)
    used_families = {item["family"] for item in selected}

    exploration = [
        item
        for item in features
        if item["candidate"].get("exclusion_reason") == NO_DIRECT_REASON
        and _clean_for_portfolio(item)
        and item["family"] not in used_families
        and item["family"] not in BATCH10_FAILED_FAMILIES
    ]

    while len(selected) < TOP_N:
        available = [
            item for item in exploration if item["family"] not in used_families
        ]
        if not available:
            raise RuntimeError("not enough family-diverse exploration candidates")

        # Submission-level portfolio heuristic:
        # 1) maximize distance from everything already selected
        # 2) prefer an explicitly donor-anchored design
        # 3) prefer zero Cys, then fewer Cys
        # 4) deterministic alias tie-break
        ranked = sorted(
            available,
            key=lambda item: (
                -min_distance_to_selected(item, selected),
                -int("donor-anchored" in item["alias"]),
                -int(item["cys_count"] == 0),
                int(item["cys_count"]),
                item["alias"],
                item["sequence"],
            ),
        )
        chosen = ranked[0]
        chosen["portfolio_min_distance"] = min_distance_to_selected(chosen, selected)
        selected.append(chosen)
        used_families.add(chosen["family"])

    family_counts: Dict[str, int] = {}
    for item in selected:
        family_counts[item["family"]] = family_counts.get(item["family"], 0) + 1

    selected_sequences = {item["sequence"] for item in selected}
    reserve_pool = [
        item
        for item in features
        if item["sequence"] not in selected_sequences
        and _clean_for_portfolio(item)
    ]

    # Reserve order intentionally puts direct-readout evidence first, but it is
    # outside the guaranteed top-20 screening window. This is where known
    # batch-10 failures and retired direct-readout families may live.
    reserve_pool.sort(
        key=lambda item: (
            -int(bool(item["candidate"].get("any_direct_readout"))),
            -int(item["coord_g12d_success_records"]),
            -int(item["g12d_success_records"]),
            -int(not bool(item["candidate"].get("exclusion_reason"))),
            -int(item["cys_count"] == 0),
            item["alias"],
            item["sequence"],
        )
    )

    for item in reserve_pool:
        if len(selected) >= TOTAL_N:
            break
        if family_counts.get(item["family"], 0) >= 2:
            continue
        selected.append(item)
        family_counts[item["family"]] = family_counts.get(item["family"], 0) + 1

    if len(selected) != TOTAL_N:
        raise RuntimeError(f"expected {TOTAL_N} portfolio designs, got {len(selected)}")

    for rank, item in enumerate(selected, start=1):
        if rank <= 3:
            bucket = "exploit-lead"
            rationale = "corrected matched-panel lead retained as exploitation slot"
        elif rank <= TOP_N:
            bucket = "exploration-top20"
            rationale = (
                "family-diverse exploration: no successful direct Asp6 readout "
                "recorded, but avoids known batch-10 family failure and sequence liabilities"
            )
        else:
            bucket = "reserve"
            rationale = (
                "reserve outside top-20 window; may include known failed/retired "
                "families if higher-ranked designs are filtered"
            )
        item["rank"] = rank
        item["bucket"] = bucket
        item["rationale"] = rationale
        item["name"] = (
            f"KRAS_PORT_{rank:02d}_"
            + hashlib.sha256(item["sequence"].encode()).hexdigest()[:8]
        )
    return selected


def _write_csv(path: Path, rows: Sequence[Mapping], columns: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(columns))
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row.get(column, "") for column in columns})


def build() -> List[Dict]:
    pool = json.loads(POOL.read_text())
    rows = select_portfolio(pool)

    _write_csv(
        OUTDIR / "track1_portfolio30.csv",
        [
            {
                "name": row["name"],
                "sequence": row["sequence"],
                "molecule_class": "protein",
            }
            for row in rows
        ],
        ("name", "sequence", "molecule_class"),
    )

    audit_rows = []
    for row in rows:
        candidate = row["candidate"]
        audit_rows.append(
            {
                "rank": row["rank"],
                "bucket": row["bucket"],
                "name": row["name"],
                "alias": row["alias"],
                "family": row["family"],
                "sequence": row["sequence"],
                "length": row["length"],
                "cys_count": row["cys_count"],
                "sequence_filter_failures": ";".join(row["sequence_filter_failures"]),
                "direct_readout": bool(candidate.get("any_direct_readout")),
                "coord_g12d_success_records": row["coord_g12d_success_records"],
                "coord_g12d_records": row["coord_g12d_records"],
                "all_g12d_success_records": row["g12d_success_records"],
                "all_g12d_records": row["g12d_records"],
                "evidence_tier": candidate.get("tier") or "",
                "exclusion_reason": candidate.get("exclusion_reason") or "",
                "review_warnings": ";".join(candidate.get("review_warnings") or []),
                "portfolio_min_distance": (
                    round(float(row["portfolio_min_distance"]), 6)
                    if row.get("portfolio_min_distance") is not None
                    else ""
                ),
                "rationale": row["rationale"],
            }
        )

    _write_csv(
        OUTDIR / "portfolio30_audit.csv",
        audit_rows,
        (
            "rank",
            "bucket",
            "name",
            "alias",
            "family",
            "sequence",
            "length",
            "cys_count",
            "sequence_filter_failures",
            "direct_readout",
            "coord_g12d_success_records",
            "coord_g12d_records",
            "all_g12d_success_records",
            "all_g12d_records",
            "evidence_tier",
            "exclusion_reason",
            "review_warnings",
            "portfolio_min_distance",
            "rationale",
        ),
    )

    summary = {
        "designs": len(rows),
        "top20_families": len({row["family"] for row in rows[:20]}),
        "max_family_count_top20": max(
            sum(other["family"] == row["family"] for other in rows[:20])
            for row in rows[:20]
        ),
        "max_family_count_full30": max(
            sum(other["family"] == row["family"] for other in rows)
            for row in rows
        ),
        "top20_gt2_cys": sum(row["cys_count"] > 2 for row in rows[:20]),
        "full30_gt2_cys": sum(row["cys_count"] > 2 for row in rows),
        "top20_sequence_filter_failures": sum(
            bool(row["sequence_filter_failures"]) for row in rows[:20]
        ),
        "top20_batch10_failed_family_exploration_slots": sum(
            row["rank"] > 3 and row["rank"] <= 20
            and row["family"] in BATCH10_FAILED_FAMILIES
            for row in rows
        ),
        "policy": {
            "top20_family_cap": 1,
            "full30_family_cap": 2,
            "max_cys": MAX_CYS,
            "ranks_1_3": "fixed corrected matched-panel leads",
            "ranks_4_20": "family-diverse exploration selected by max-min sequence distance",
            "ranks_21_30": "reserve, direct-readout evidence prioritized",
        },
    }
    (OUTDIR / "portfolio30_summary.json").write_text(json.dumps(summary, indent=2) + "\n")

    errors = validate_file(OUTDIR / "track1_portfolio30.csv")
    if errors:
        raise RuntimeError("generated portfolio failed validation: " + "; ".join(errors))
    return rows


def main() -> int:
    rows = build()
    print(f"WROTE {len(rows)} designs to {OUTDIR / 'track1_portfolio30.csv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
