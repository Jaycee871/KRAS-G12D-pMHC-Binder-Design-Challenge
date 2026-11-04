#!/usr/bin/env python3
"""Batch-05 replication: 3 independent constrained draws per screening passer.

Replication is part of the **gate** in batch-05, not a post-hoc check: a design
passes only if it clears the Asp6 direct read-out gate in >= 3 of 3 independent
constrained Boltz-2 draws (the screening call that selected it is reported
separately as draw 0 and does not count towards the 3).

This is batch-04's `scripts/replicate_batch04.py` and
`scripts/replicate_rmsd_batch04.py` re-used verbatim with their output
directory pointed at `results/batch-05` - including the binder C-alpha RMSD
spread measured in the pMHC frame with the binder *not* re-fitted.
"""
from __future__ import annotations

import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(REPO / "scripts"))

import run_batch04, replicate_batch04, replicate_rmsd_batch04  # noqa: E402

OUT = REPO / "results" / "batch-05"


def main() -> int:
    replicates = sys.argv[1] if len(sys.argv) > 1 else "3"
    for module in (run_batch04, replicate_batch04, replicate_rmsd_batch04):
        module.OUT = OUT
    (OUT / "replicate_structures").mkdir(parents=True, exist_ok=True)
    argv = sys.argv
    sys.argv = ["replicate_batch04", "--replicates", replicates]
    rc = replicate_batch04.main()
    sys.argv = ["replicate_rmsd_batch04"]
    rc = replicate_rmsd_batch04.main() or rc
    sys.argv = argv
    return rc


if __name__ == "__main__":
    raise SystemExit(main())
