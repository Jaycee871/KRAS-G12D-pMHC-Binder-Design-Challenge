import importlib.util
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "portfolio_v3", ROOT / "scripts" / "build_portfolio_v3.py"
)
v3 = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(v3)


def _rows():
    pool = json.loads(v3.POOL.read_text())
    return v3.select_portfolio(pool)


def test_v3_top20_is_balanced_and_family_diverse():
    rows = _rows()
    top20 = rows[:20]
    assert len(rows) == 30
    assert len({r["sequence"] for r in rows}) == 30
    assert len({r["family"] for r in top20}) == 20
    assert sum(bool(r["candidate"].get("any_direct_readout")) for r in top20) == 8
    assert [r["bucket"] for r in top20].count("exploit-lead") == 3
    assert [r["bucket"] for r in top20].count("evidence-risk-hedge") == 5
    assert [r["bucket"] for r in top20].count("exploration-top20") == 12


def test_v3_keeps_bb135_as_explicit_high_risk_exception():
    rows = _rows()
    top20 = rows[:20]
    bb135 = next(r for r in top20 if r["alias"] == v3.BB135_ALIAS)
    assert bb135["rank"] == 4
    assert bb135["cys_count"] == 4
    assert bb135["coord_g12d_success_records"] == 7
    assert bb135["coord_g12d_records"] == 7
    assert "high-evidence/high-risk" in bb135["risk_note"]
    assert [r["alias"] for r in top20 if r["cys_count"] > 2] == [v3.BB135_ALIAS]


def test_v3_hedges_are_exactly_predeclared():
    rows = _rows()
    assert [r["alias"] for r in rows[:3]] == list(v3.LEAD_ALIASES)
    assert [r["alias"] for r in rows[3:8]] == list(v3.HEDGE_ALIASES)


def test_v3_does_not_repopulate_top20_with_batch10_failures():
    rows = _rows()
    top_aliases = {r["alias"] for r in rows[:20]}
    exact_failed = {
        "B5-c6-62_bb173_donor-anchored-g1",
        "bb54_graft-p3",
        "B5-c6-62_bb90_donor-anchored-g2",
        "B4-c67-65_bb26_donor-anchored",
        "B6-c6-46_w2_bb38_donor-anchored-g1",
        "bb54_graft-p1-v3",
        "bb54_graft-p2",
    }
    # bb36 is the single deliberate batch-10 evidence-risk hedge.
    assert top_aliases.isdisjoint(exact_failed)
    assert "B6-c6-46_w2_bb36_donor-anchored-g2" in top_aliases


def test_v3_non_bb135_sequences_keep_two_cys_cap_and_sequence_filters():
    rows = _rows()
    for row in rows:
        if row["alias"] == v3.BB135_ALIAS:
            continue
        assert row["cys_count"] <= 2
        assert not row["sequence_filter_failures"]
