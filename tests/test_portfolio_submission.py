import importlib.util
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "build_portfolio_submission",
    ROOT / "scripts" / "build_portfolio_submission.py",
)
portfolio = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(portfolio)


def _rows():
    pool = json.loads(portfolio.POOL.read_text())
    return portfolio.select_portfolio(pool)


def test_portfolio_has_thirty_unique_sequences_and_twenty_unique_top_families():
    rows = _rows()
    assert len(rows) == 30
    assert len({row["sequence"] for row in rows}) == 30
    assert len({row["name"] for row in rows}) == 30
    assert len({row["family"] for row in rows[:20]}) == 20


def test_top20_family_cap_and_full30_family_cap_are_locked():
    rows = _rows()
    top20 = rows[:20]
    full = rows
    assert max(
        sum(other["family"] == row["family"] for other in top20)
        for row in top20
    ) == 1
    assert max(
        sum(other["family"] == row["family"] for other in full)
        for row in full
    ) <= 2


def test_top20_keeps_three_leads_then_uses_exploration_not_batch10_failures():
    rows = _rows()
    assert [row["alias"] for row in rows[:3]] == list(portfolio.LEAD_ALIASES)
    for row in rows[3:20]:
        assert row["candidate"]["exclusion_reason"] == portfolio.NO_DIRECT_REASON
        assert row["family"] not in portfolio.BATCH10_FAILED_FAMILIES


def test_portfolio_removes_high_cys_burden_and_repo_sequence_liabilities():
    rows = _rows()
    assert all(row["cys_count"] <= 2 for row in rows)
    assert all(not row["sequence_filter_failures"] for row in rows)
    aliases = {row["alias"] for row in rows}
    assert "B5-c6-62_bb135_donor-anchored-g1" not in aliases
    assert "B8-span-50_bb26_donor-anchored-g2" not in aliases


def test_known_batch10_exact_failures_are_outside_top20():
    rows = _rows()
    top_aliases = {row["alias"] for row in rows[:20]}
    failed_aliases = {
        "B6-c6-46_w2_bb36_donor-anchored-g2",
        "B5-c6-62_bb173_donor-anchored-g1",
        "bb54_graft-p3",
        "B5-c6-62_bb90_donor-anchored-g2",
        "B4-c67-65_bb26_donor-anchored",
        "B6-c6-46_w2_bb38_donor-anchored-g1",
        "bb54_graft-p1-v3",
        "bb54_graft-p2",
    }
    assert top_aliases.isdisjoint(failed_aliases)
