import csv
from pathlib import Path

from src.validate_submission import validate_file

ROOT = Path(__file__).resolve().parents[1]
SUBMISSION = ROOT / "results" / "final-submission" / "track1_fallback20.csv"
AUDIT = ROOT / "results" / "final-submission" / "fallback20_audit.csv"
SHORTLIST = (
    ROOT
    / "results"
    / "readout-correction"
    / "submission-review"
    / "review_shortlist.csv"
)


def _rows(path):
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def test_fallback20_is_schema_valid_and_exactly_twenty_unique_designs():
    assert validate_file(SUBMISSION) == []
    rows = _rows(SUBMISSION)
    assert len(rows) == 20
    assert len({row["name"] for row in rows}) == 20
    assert len({row["sequence"] for row in rows}) == 20
    assert all(row["molecule_class"] == "protein" for row in rows)


def test_first_eleven_are_exactly_the_curated_review_shortlist():
    submitted = _rows(SUBMISSION)
    shortlist = _rows(SHORTLIST)
    assert len(shortlist) == 11
    assert [row["sequence"] for row in submitted[:11]] == [
        row["sequence"] for row in shortlist
    ]


def test_tail_is_explicitly_fallback_only_in_audit():
    submitted = _rows(SUBMISSION)
    audit = _rows(AUDIT)
    assert len(audit) == 20
    assert [row["sequence"] for row in audit] == [
        row["sequence"] for row in submitted
    ]
    assert all(
        row["disposition"] == "retained-review" for row in audit[:11]
    )
    assert all(row["disposition"] == "fallback-only" for row in audit[11:])
    assert all(row["exclusion"].strip() for row in audit[11:])
