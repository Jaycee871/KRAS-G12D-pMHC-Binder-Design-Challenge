#!/usr/bin/env python3
"""Validate a Proteinbase Challenge 5 submission file (CSV or FASTA).

Rules enforced (from the Anthropic x Adaptyv competition submission page):
  * file size <= 5 MB
  * CSV must have `name` and `sequence` columns; `molecule_class` optional
  * at most 20 designs per submission
  * each chain between 10 and 250 amino acids (inclusive)
  * only the 20 canonical amino acids
  * scFv / Fab written as {VH}:{VL} - exactly one colon, both parts valid
  * names must be present and unique

Exit code 0 = valid, 1 = invalid, 2 = usage/IO error.

Usage:
    python3 src/validate_submission.py submission.csv
"""
from __future__ import annotations

import argparse
import csv
import os
import sys
from pathlib import Path
from typing import List, Tuple

MAX_DESIGNS = 20
MIN_CHAIN_LEN = 10
MAX_CHAIN_LEN = 250
MAX_FILE_BYTES = 5 * 1024 * 1024  # 5 MB
VALID_AA = set("ACDEFGHIKLMNPQRSTVWY")
REQUIRED_COLUMNS = ("name", "sequence")
OPTIONAL_COLUMNS = ("molecule_class",)
MULTI_CHAIN_CLASSES = {"scfv", "fab"}


class ValidationError(Exception):
    pass


def parse_csv(path: Path) -> List[dict]:
    with path.open(newline="", encoding="utf-8-sig") as fh:
        reader = csv.DictReader(fh)
        if reader.fieldnames is None:
            raise ValidationError("CSV is empty: no header row found")
        cols = [c.strip().lower() for c in reader.fieldnames if c is not None]
        missing = [c for c in REQUIRED_COLUMNS if c not in cols]
        if missing:
            raise ValidationError(
                f"missing required column(s): {', '.join(missing)} "
                f"(found: {', '.join(cols) or 'none'})"
            )
        rows = []
        for raw in reader:
            row = {
                (k.strip().lower() if k else ""): (v.strip() if isinstance(v, str) else "")
                for k, v in raw.items()
            }
            if not any(row.get(c) for c in REQUIRED_COLUMNS):
                continue  # skip fully blank lines
            rows.append(row)
    return rows


def parse_fasta(path: Path) -> List[dict]:
    rows: List[dict] = []
    name = None
    chunks: List[str] = []
    with path.open(encoding="utf-8-sig") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    rows.append({"name": name, "sequence": "".join(chunks)})
                name = line[1:].strip()
                chunks = []
            else:
                if name is None:
                    raise ValidationError(
                        "FASTA starts with sequence data before any '>' header"
                    )
                chunks.append(line)
    if name is not None:
        rows.append({"name": name, "sequence": "".join(chunks)})
    if not rows:
        raise ValidationError("FASTA contains no records")
    return rows


def load_rows(path: Path) -> List[dict]:
    suffix = path.suffix.lower()
    if suffix == ".csv":
        return parse_csv(path)
    if suffix in {".fasta", ".fa", ".faa", ".fas"}:
        return parse_fasta(path)
    raise ValidationError(
        f"unsupported file type '{suffix or path.name}': expected .csv or .fasta"
    )


def validate_sequence(seq: str, molecule_class: str = "") -> List[str]:
    """Return a list of problems with one design's sequence (empty = ok)."""
    problems: List[str] = []
    seq = seq.strip().upper().replace(" ", "")
    mc = (molecule_class or "").strip().lower()

    if not seq:
        return ["sequence is empty"]

    n_colons = seq.count(":")
    if n_colons > 1:
        return [
            f"malformed multi-chain sequence: expected exactly one ':' "
            f"separating {{VH}}:{{VL}}, found {n_colons}"
        ]
    if n_colons == 0 and mc in MULTI_CHAIN_CLASSES:
        problems.append(
            f"molecule_class '{mc}' requires a two-chain sequence written as "
            "{VH}:{VL}, but no ':' was found"
        )

    chains = seq.split(":")
    labels = ["VH", "VL"] if n_colons == 1 else ["chain"]
    for label, chain in zip(labels, chains):
        if not chain:
            problems.append(f"{label} is empty (stray ':' in sequence)")
            continue
        bad = sorted(set(chain) - VALID_AA)
        if bad:
            problems.append(
                f"{label} contains invalid amino acid letter(s): {', '.join(bad)}"
            )
        if len(chain) < MIN_CHAIN_LEN:
            problems.append(
                f"{label} is {len(chain)} aa, below the minimum of {MIN_CHAIN_LEN}"
            )
        elif len(chain) > MAX_CHAIN_LEN:
            problems.append(
                f"{label} is {len(chain)} aa, above the maximum of {MAX_CHAIN_LEN}"
            )
    return problems


def validate_file(path: Path) -> List[str]:
    """Return a list of human-readable errors; empty list means valid."""
    errors: List[str] = []

    if not path.exists():
        raise ValidationError(f"file not found: {path}")

    size = path.stat().st_size
    if size > MAX_FILE_BYTES:
        errors.append(
            f"file is {size / 1024 / 1024:.2f} MB, above the 5 MB limit"
        )
    if size == 0:
        raise ValidationError(f"file is empty: {path}")

    rows = load_rows(path)

    if len(rows) == 0:
        errors.append("submission contains no designs")
    if len(rows) > MAX_DESIGNS:
        errors.append(
            f"submission contains {len(rows)} designs, above the limit of "
            f"{MAX_DESIGNS} per submission"
        )

    seen = {}
    for idx, row in enumerate(rows, start=1):
        name = (row.get("name") or "").strip()
        label = f"row {idx}" + (f" ('{name}')" if name else "")
        if not name:
            errors.append(f"{label}: missing name")
        elif name in seen:
            errors.append(f"{label}: duplicate name (also row {seen[name]})")
        else:
            seen[name] = idx
        for problem in validate_sequence(
            row.get("sequence", ""), row.get("molecule_class", "")
        ):
            errors.append(f"{label}: {problem}")

    return errors


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("path", type=Path, help="submission .csv or .fasta")
    ap.add_argument("--quiet", action="store_true", help="only print on failure")
    args = ap.parse_args(argv)

    try:
        errors = validate_file(args.path)
    except ValidationError as exc:
        print(f"INVALID: {exc}", file=sys.stderr)
        return 1
    except OSError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2

    if errors:
        print(f"INVALID: {args.path} failed {len(errors)} check(s):", file=sys.stderr)
        for err in errors:
            print(f"  - {err}", file=sys.stderr)
        return 1

    if not args.quiet:
        rows = load_rows(args.path)
        print(f"VALID: {args.path} ({len(rows)} design(s), "
              f"{os.path.getsize(args.path)} bytes)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
