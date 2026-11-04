#!/usr/bin/env python3
"""Build portfolio v3: evidence-risk hedges + family-diverse exploration.

This is a submission allocation layer only. It does not change any biological
gate or historical batch conclusion.

Top-20 allocation:
  1-3   corrected matched-panel leads
  4-8   evidence-bearing, explicitly high-risk hedges
  9-20  family-diverse exploration

The bb135 4-Cys sequence is the sole >2-Cys exception because it has the
strongest oxygen-specific Asp6 readout in the inventory. Its apo/decoy and
developability risks remain explicit in the audit.
"""
from __future__ import annotations

import csv
import hashlib
import importlib.util
import json
from pathlib import Path
from typing import Dict, List, Mapping, Sequence

ROOT = Path(__file__).resolve().parents[1]
OUTDIR = ROOT / "results" / "final-submission"
POOL = ROOT / "results" / "readout-correction" / "submission-review" / "evidence_pool.json"

spec = importlib.util.spec_from_file_location(
    "portfolio_v2", ROOT / "scripts" / "build_portfolio_submission.py"
)
v2 = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(v2)

LEAD_ALIASES = v2.LEAD_ALIASES
HEDGE_ALIASES = (
    "B5-c6-62_bb135_donor-anchored-g1",
    "B6-c6-46_w2_bb36_donor-anchored-g2",
    "B8-apex-54_bb06_donor-anchored-g2",
    "B8-apex-46_bb08_donor-anchored-g1-dw",
    "B8-span-50_bb28_donor-anchored-g1-dw",
)
BB135_ALIAS = HEDGE_ALIASES[0]
NO_DIRECT_REASON = v2.NO_DIRECT_REASON
BATCH10_FAILED_FAMILIES = v2.BATCH10_FAILED_FAMILIES
TOP_N = 20
TOTAL_N = 30

RISK_NOTES = {
    "B5-c6-62_bb135_donor-anchored-g1": (
        "high-evidence/high-risk: oxygen-specific G12D Asp6 7/7; WT direct "
        "readout 0/4; unrelated-decoy Asp8 donor contact 3/4; apo total buried "
        "area 1101 A^2 versus ~630 A^2 G12D in the earlier asymmetric panel; "
        "4 Cys. Two G12D predictions repeatedly place C20-C53 and C27-C42 as "
        "pair-like SG-SG contacts, but no disulfide was explicitly constrained."
    ),
    "B6-c6-46_w2_bb36_donor-anchored-g2": (
        "evidence-bearing/known-risk: oxygen-specific G12D 2/6 in batch-06; "
        "batch-10 matched G12D 1/3 and apo/G12D median total-area ratio 1.258, "
        "so it failed the frozen batch-10 advance rule."
    ),
    "B8-apex-54_bb06_donor-anchored-g2": (
        "limited positive hedge: oxygen-specific G12D 1/3; batch-08 family was "
        "retired after replication failure; matched apo counter-screen not run."
    ),
    "B8-apex-46_bb08_donor-anchored-g1-dw": (
        "limited positive hedge: oxygen-specific G12D 1/3; batch-08 family was "
        "retired after replication failure; matched apo counter-screen not run."
    ),
    "B8-span-50_bb28_donor-anchored-g1-dw": (
        "limited positive hedge: oxygen-specific G12D 1/3; batch-08 family was "
        "retired after replication failure; matched apo counter-screen not run."
    ),
}


def _find(features: Sequence[Mapping], alias: str) -> Dict:
    matches = [
        dict(item) for item in features
        if alias in (item["candidate"].get("aliases") or [])
    ]
    if len(matches) != 1:
        raise RuntimeError(f"expected exactly one candidate for {alias}, got {len(matches)}")
    return matches[0]


def _eligible_exploration(item: Mapping, used_families: set[str]) -> bool:
    return (
        item["candidate"].get("exclusion_reason") == NO_DIRECT_REASON
        and 40 <= int(item["length"]) <= 100
        and int(item["cys_count"]) <= 2
        and not item["sequence_filter_failures"]
        and item["family"] not in used_families
        and item["family"] not in BATCH10_FAILED_FAMILIES
    )


def _reserve_clean(item: Mapping) -> bool:
    if not (40 <= int(item["length"]) <= 100):
        return False
    if item["sequence_filter_failures"]:
        return False
    if item["alias"] == BB135_ALIAS:
        return True
    return int(item["cys_count"]) <= 2


def select_portfolio(pool: Mapping) -> List[Dict]:
    features = [v2._feature(c) for c in pool["candidates"]]

    selected: List[Dict] = []
    for alias in LEAD_ALIASES:
        item = _find(features, alias)
        item["bucket"] = "exploit-lead"
        item["risk_note"] = (
            "strongest corrected matched-panel evidence; unresolved matched apo "
            "over-burial remains a negative signal"
        )
        selected.append(item)

    for alias in HEDGE_ALIASES:
        item = _find(features, alias)
        item["bucket"] = "evidence-risk-hedge"
        item["risk_note"] = RISK_NOTES[alias]
        selected.append(item)

    if len({x["family"] for x in selected}) != len(selected):
        raise RuntimeError("lead/hedge block must use distinct families")

    used_families = {x["family"] for x in selected}
    exploration = [
        item for item in features if _eligible_exploration(item, used_families)
    ]

    while len(selected) < TOP_N:
        available = [
            item for item in exploration if item["family"] not in used_families
        ]
        if not available:
            raise RuntimeError("not enough family-diverse exploration candidates")
        ranked = sorted(
            available,
            key=lambda item: (
                -v2.min_distance_to_selected(item, selected),
                -int("donor-anchored" in item["alias"]),
                -int(item["cys_count"] == 0),
                int(item["cys_count"]),
                item["alias"],
                item["sequence"],
            ),
        )
        chosen = dict(ranked[0])
        chosen["portfolio_min_distance"] = v2.min_distance_to_selected(chosen, selected)
        chosen["bucket"] = "exploration-top20"
        chosen["risk_note"] = (
            "weak negative / under-tested: no successful direct Asp6 readout "
            "recorded; generally one initial G12D screen and no oxygen-specific "
            "coordinate remeasurement"
        )
        selected.append(chosen)
        used_families.add(chosen["family"])

    family_counts: Dict[str, int] = {}
    for item in selected:
        family_counts[item["family"]] = family_counts.get(item["family"], 0) + 1

    selected_sequences = {item["sequence"] for item in selected}
    reserve = [
        item for item in features
        if item["sequence"] not in selected_sequences and _reserve_clean(item)
    ]
    reserve.sort(
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
    for raw in reserve:
        if len(selected) >= TOTAL_N:
            break
        if family_counts.get(raw["family"], 0) >= 2:
            continue
        item = dict(raw)
        item["bucket"] = "reserve"
        item["risk_note"] = (
            "reserve outside top-20 allocation; historical negative evidence and "
            "exclusion reasons remain authoritative"
        )
        selected.append(item)
        family_counts[item["family"]] = family_counts.get(item["family"], 0) + 1

    if len(selected) != TOTAL_N:
        raise RuntimeError(f"expected {TOTAL_N} designs, got {len(selected)}")

    for rank, item in enumerate(selected, start=1):
        item["rank"] = rank
        item["name"] = (
            f"KRAS_V3_{rank:02d}_"
            + hashlib.sha256(item["sequence"].encode()).hexdigest()[:8]
        )
    return selected


def _write_csv(path: Path, rows: Sequence[Mapping], columns: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(columns))
        writer.writeheader()
        for row in rows:
            writer.writerow({c: row.get(c, "") for c in columns})


def build() -> List[Dict]:
    pool = json.loads(POOL.read_text())
    rows = select_portfolio(pool)

    upload_rows = [
        {"name": r["name"], "sequence": r["sequence"], "molecule_class": "protein"}
        for r in rows
    ]
    upload = OUTDIR / "track1_portfolio30_v3.csv"
    _write_csv(upload, upload_rows, ("name", "sequence", "molecule_class"))

    audit = []
    for r in rows:
        c = r["candidate"]
        audit.append(
            {
                "rank": r["rank"],
                "bucket": r["bucket"],
                "name": r["name"],
                "alias": r["alias"],
                "family": r["family"],
                "sequence": r["sequence"],
                "length": r["length"],
                "cys_count": r["cys_count"],
                "direct_readout": bool(c.get("any_direct_readout")),
                "coord_g12d_success_records": r["coord_g12d_success_records"],
                "coord_g12d_records": r["coord_g12d_records"],
                "all_g12d_success_records": r["g12d_success_records"],
                "all_g12d_records": r["g12d_records"],
                "evidence_tier": c.get("tier") or "",
                "exclusion_reason": c.get("exclusion_reason") or "",
                "review_warnings": ";".join(c.get("review_warnings") or []),
                "portfolio_min_distance": (
                    round(float(r["portfolio_min_distance"]), 6)
                    if r.get("portfolio_min_distance") is not None else ""
                ),
                "risk_note": r["risk_note"],
            }
        )
    _write_csv(
        OUTDIR / "portfolio30_v3_audit.csv",
        audit,
        (
            "rank","bucket","name","alias","family","sequence","length","cys_count",
            "direct_readout","coord_g12d_success_records","coord_g12d_records",
            "all_g12d_success_records","all_g12d_records","evidence_tier",
            "exclusion_reason","review_warnings","portfolio_min_distance","risk_note",
        ),
    )

    top20 = rows[:20]
    summary = {
        "designs": len(rows),
        "top20_families": len({r["family"] for r in top20}),
        "top20_direct_readout_designs": sum(bool(r["candidate"].get("any_direct_readout")) for r in top20),
        "top20_buckets": {
            "exploit_lead": sum(r["bucket"] == "exploit-lead" for r in top20),
            "evidence_risk_hedge": sum(r["bucket"] == "evidence-risk-hedge" for r in top20),
            "exploration": sum(r["bucket"] == "exploration-top20" for r in top20),
        },
        "top20_gt2_cys": [
            {"rank": r["rank"], "alias": r["alias"], "cys_count": r["cys_count"]}
            for r in top20 if r["cys_count"] > 2
        ],
        "max_family_count_top20": max(
            sum(other["family"] == r["family"] for other in top20) for r in top20
        ),
        "max_family_count_full30": max(
            sum(other["family"] == r["family"] for other in rows) for r in rows
        ),
        "bb135_policy": (
            "explicit high-evidence/high-risk exception; 4 Cys is a risk flag, "
            "not a hard exclusion"
        ),
    }
    (OUTDIR / "portfolio30_v3_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n"
    )

    errors = v2.validate_file(upload)
    if errors:
        raise RuntimeError("v3 validation failed: " + "; ".join(errors))
    return rows


if __name__ == "__main__":
    rows = build()
    print(f"WROTE {len(rows)} designs")
