#!/usr/bin/env python3
"""Build a compliant Proteinbase submission CSV from a scored designs table.

Input: a CSV of scored designs with at least `name` and `sequence` columns.
`molecule_class` is preserved when supplied; otherwise a single-chain design
is emitted as the official Proteinbase class `protein`.
Scoring columns are optional; by default designs are ranked by the
`specificity_gap` column (descending) if present, else by `score`, else the
input order is preserved.

The output is validated with validate_submission.py before being written, so a
non-compliant file is never produced.

Usage:
    python3 src/make_submission.py scored.csv -o submission.csv -n 20
"""
from __future__ import annotations

import argparse
import csv
import sys
import tempfile
from pathlib import Path
from typing import List

try:  # allow both `python src/make_submission.py` and `from src...` imports
    from validate_submission import MAX_DESIGNS, validate_file, ValidationError
except ImportError:  # pragma: no cover - import shim
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from validate_submission import MAX_DESIGNS, validate_file, ValidationError

DEFAULT_RANK_COLUMNS = ("specificity_gap", "score")
OUTPUT_COLUMNS = ("name", "sequence", "molecule_class")
DEFAULT_TOP_N = 20  # Track 1 screens the top 20; files may contain up to 40.


def read_scored(path: Path) -> List[dict]:
    with path.open(newline="", encoding="utf-8-sig") as fh:
        reader = csv.DictReader(fh)
        if reader.fieldnames is None:
            raise ValidationError(f"{path} is empty: no header row")
        cols = [c.strip().lower() for c in reader.fieldnames if c]
        for required in ("name", "sequence"):
            if required not in cols:
                raise ValidationError(
                    f"{path} is missing required column '{required}' "
                    f"(found: {', '.join(cols) or 'none'})"
                )
        rows = []
        for raw in reader:
            row = {
                (k.strip().lower() if k else ""):
                    (v.strip() if isinstance(v, str) else "")
                for k, v in raw.items()
            }
            if row.get("name") or row.get("sequence"):
                rows.append(row)
    return rows


def pick_rank_column(rows: List[dict], explicit: str | None) -> str | None:
    if explicit:
        if not rows or explicit not in rows[0]:
            raise ValidationError(f"rank column '{explicit}' not found in input")
        return explicit
    if not rows:
        return None
    for col in DEFAULT_RANK_COLUMNS:
        if col in rows[0]:
            return col
    return None


def rank(rows: List[dict], column: str | None, ascending: bool) -> List[dict]:
    if column is None:
        return list(rows)

    def key(row):
        raw = row.get(column, "")
        try:
            return (0, float(raw))
        except (TypeError, ValueError):
            # Unscored rows sort last regardless of direction.
            return (1, 0.0)

    scored = [r for r in rows if key(r)[0] == 0]
    unscored = [r for r in rows if key(r)[0] == 1]
    scored.sort(key=lambda r: float(r[column]), reverse=not ascending)
    return scored + unscored


def write_submission(rows: List[dict], path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    cols = list(OUTPUT_COLUMNS)
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=cols)
        writer.writeheader()
        for row in rows:
            sequence = (row.get("sequence") or "").strip()
            molecule_class = (row.get("molecule_class") or "").strip().lower()
            if not molecule_class and ":" not in sequence:
                molecule_class = "protein"
            writer.writerow({
                "name": row.get("name", ""),
                "sequence": sequence,
                "molecule_class": molecule_class,
            })
    return path


def build(
    source: Path,
    out: Path,
    top_n: int = DEFAULT_TOP_N,
    rank_column: str | None = None,
    ascending: bool = False,
) -> Path:
    """Read, rank, truncate, validate and write. Returns the output path."""
    if top_n < 1:
        raise ValidationError(f"--top-n must be >= 1 (got {top_n})")
    if top_n > MAX_DESIGNS:
        raise ValidationError(
            f"--top-n {top_n} exceeds the Track 1 submission maximum of {MAX_DESIGNS}"
        )

    rows = read_scored(source)
    if not rows:
        raise ValidationError(f"{source} contains no designs")

    column = pick_rank_column(rows, rank_column)
    selected = rank(rows, column, ascending)[:top_n]

    # Validate a temp copy first so we never leave a bad file on disk.
    with tempfile.TemporaryDirectory() as tmp:
        probe = write_submission(selected, Path(tmp) / "probe.csv")
        errors = validate_file(probe)
    if errors:
        raise ValidationError(
            "generated submission would be invalid:\n  - " + "\n  - ".join(errors)
        )

    write_submission(selected, out)
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("source", type=Path, help="scored designs CSV")
    ap.add_argument("-o", "--out", type=Path, default=Path("submission.csv"))
    ap.add_argument("-n", "--top-n", type=int, default=DEFAULT_TOP_N)
    ap.add_argument("--rank-column", default=None,
                    help="column to rank by (default: specificity_gap, then score)")
    ap.add_argument("--ascending", action="store_true",
                    help="rank low-to-high (e.g. when ranking by KD)")
    args = ap.parse_args(argv)

    try:
        out = build(args.source, args.out, args.top_n,
                    args.rank_column, args.ascending)
    except ValidationError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    except OSError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2

    n = sum(1 for _ in out.open()) - 1
    print(f"WROTE: {out} ({n} design(s)) - validated against all submission rules")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
