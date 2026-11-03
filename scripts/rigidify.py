#!/usr/bin/env python3
"""Rigidification helpers for batch-07: make the donor contact constitutive.

Batches 04, 05 and 06 all stalled at 2-of-3 on the Asp6 direct read-out, and
batch-06's diagnosis was explicit: the contact is **intermittent, not weak**.
`bb54-g2` reads the carboxylate at 2.5-3.5 A when it reads it at all, so the
problem is not chemistry but the mobility of the element that carries the
donor. These helpers attack exactly that, in two ways, both deterministic and
both measured on the *design* backbone before any Boltz-2 call is spent:

* `secondary_structure` / `loop_burden` - a CA-geometry secondary-structure
  assignment (no DSSP dependency, no hydrogens on an RFdiffusion backbone) used
  to ask whether the donor-carrying residues sit on a regular helix/strand or
  dangle on a long loop;
* `rigidifying_edits` - deterministic sequence edits that stiffen the scaffold
  *away* from the paratope: Gly -> Ala at loop positions with a negative phi
  (left alone where phi is positive, since only Gly is allowed there),
  Pro at loop positions whose phi sits in the proline window, and at most one
  Cys-Cys staple across a pair of CB atoms 3.3-5.2 A apart and >= 5 residues
  apart in sequence.

The Asp6 donor shell is **protected**: no position within `shell_protect` of
position 6 is ever edited, so the donor geometry that works is kept exactly as
`run_batch04.donor_anchored_variant` placed it. Nothing here shrinks the
paratope - batch-06 measured that shrinking does not re-centre the interface.
"""
from __future__ import annotations

import sys
from pathlib import Path
from typing import Any, Dict, List, Sequence, Tuple

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

import complex_analysis  # noqa: E402

HELIX_13 = (4.8, 6.3)       # CA(i)-CA(i+3) in an alpha helix
HELIX_14 = (5.6, 7.2)       # CA(i)-CA(i+4) in an alpha helix
STRAND_12 = 6.4             # CA(i)-CA(i+2) above this is extended
PRO_PHI = (-95.0, -45.0)    # proline-tolerated phi window
SS_PAIR = (3.3, 6.5)        # CB-CB of a disulfide (3.5-4.5 A ideal; up to 6.5 A
                            # is taken as "formable after local relaxation" and is
                            # labelled as such - it is a hypothesis, not a bond)


def binder_residues(pdb_text: str, binder_chain: str = "A") -> List[Any]:
    structure = complex_analysis.parse_structure(pdb_text, fmt="pdb")
    model = complex_analysis._clean_model(next(iter(structure)))
    return [r for r in model[binder_chain]
            if complex_analysis.is_polymer_residue(r)]


def _coords(residues, name: str):
    out = []
    for r in residues:
        out.append(np.array(r[name].coord, dtype=float) if name in r else None)
    return out


def dihedral(p0, p1, p2, p3) -> float:
    b0, b1, b2 = p0 - p1, p2 - p1, p3 - p2
    b1n = b1 / np.linalg.norm(b1)
    v = b0 - np.dot(b0, b1n) * b1n
    w = b2 - np.dot(b2, b1n) * b1n
    x, y = np.dot(v, w), np.dot(np.cross(b1n, v), w)
    return float(np.degrees(np.arctan2(y, x)))


def phi_angles(residues) -> List[float]:
    """Backbone phi per residue (None-safe: 999.0 where it cannot be formed)."""
    N, CA, C = _coords(residues, "N"), _coords(residues, "CA"), _coords(residues, "C")
    out = [999.0] * len(residues)
    for i in range(1, len(residues)):
        if any(x is None for x in (C[i - 1], N[i], CA[i], C[i])):
            continue
        out[i] = dihedral(C[i - 1], N[i], CA[i], C[i])
    return out


def secondary_structure(residues) -> str:
    """Per-residue H / E / L string from CA geometry alone."""
    ca = _coords(residues, "CA")
    n = len(ca)
    ss = ["L"] * n
    def d(i, j):
        if 0 <= i < n and 0 <= j < n and ca[i] is not None and ca[j] is not None:
            return float(np.linalg.norm(ca[i] - ca[j]))
        return None
    for i in range(n):
        d13, d14 = d(i, i + 3), d(i, i + 4)
        if (d13 is not None and d14 is not None
                and HELIX_13[0] <= d13 <= HELIX_13[1]
                and HELIX_14[0] <= d14 <= HELIX_14[1]):
            for k in range(i, min(i + 5, n)):
                ss[k] = "H"
    for i in range(n):
        d12 = d(i, i + 2)
        if ss[i] == "L" and d12 is not None and d12 >= STRAND_12:
            d12b = d(i - 2, i)
            if d12b is not None and d12b >= STRAND_12:
                ss[i] = "E"
    return "".join(ss)


def loop_runs(ss: str) -> List[Tuple[int, int]]:
    """1-based inclusive (start, end) of every loop run."""
    runs, start = [], None
    for index, char in enumerate(ss, start=1):
        if char == "L" and start is None:
            start = index
        elif char != "L" and start is not None:
            runs.append((start, index - 1)); start = None
    if start is not None:
        runs.append((start, len(ss)))
    return runs


def donor_loop_metrics(pdb_text: str, shell: Sequence[Dict[str, Any]],
                       binder_chain: str = "A",
                       shell_depth: int = 4) -> Dict[str, Any]:
    """How rigid is the element carrying the Asp6 donor?

    `donor_loop_length` is the length of the loop run containing the nearest
    shell position (0 when that position is already on a helix or strand), and
    `shell_ss` the secondary structure under the `shell_depth` nearest shell
    positions. A *constitutive* contact wants the donor on, or one residue
    from, a regular element - not in the middle of a long flexible loop.
    """
    residues = binder_residues(pdb_text, binder_chain)
    ss = secondary_structure(residues)
    runs = loop_runs(ss)
    near = [entry["position"] for entry in shell[:shell_depth]]
    anchor = near[0] if near else None
    loop_len, loop_span = 0, None
    if anchor is not None and anchor <= len(ss) and ss[anchor - 1] == "L":
        for start, end in runs:
            if start <= anchor <= end:
                loop_len, loop_span = end - start + 1, (start, end)
                break
    helix_frac = ss.count("H") / len(ss) if ss else 0.0
    return {"ss": ss, "ss_helix_frac": round(helix_frac, 3),
            "ss_loop_frac": round(ss.count("L") / len(ss), 3) if ss else None,
            "anchor_position": anchor,
            "anchor_ss": ss[anchor - 1] if anchor and anchor <= len(ss) else None,
            "shell_ss": "".join(ss[p - 1] for p in near if p <= len(ss)),
            "donor_loop_length": loop_len, "donor_loop_span": loop_span,
            "n_loop_runs": len(runs),
            "longest_loop": max((e - s + 1 for s, e in runs), default=0)}


def rigidifying_edits(sequence: str, pdb_text: str,
                      shell: Sequence[Dict[str, Any]],
                      binder_chain: str = "A", shell_protect: int = 6,
                      do_proline: bool = True, do_staple: bool = False,
                      ) -> Tuple[str, List[str]]:
    """Deterministic stiffening edits away from the protected donor shell."""
    residues = binder_residues(pdb_text, binder_chain)
    ss = secondary_structure(residues)
    phi = phi_angles(residues)
    protected = {entry["position"] for entry in shell[:shell_protect]}
    out, notes = list(sequence), []
    n = min(len(out), len(ss))

    # 1. Gly -> Ala on loops, except where a positive phi needs a glycine.
    for i in range(1, n + 1):
        if i in protected or ss[i - 1] != "L" or out[i - 1] != "G":
            continue
        if phi[i - 1] > 0:
            notes.append(f"G{i} kept (phi {phi[i - 1]:.0f} positive)")
            continue
        out[i - 1] = "A"
        notes.append(f"G{i}A")

    # 2. Proline at loop positions inside the proline phi window, never at the
    #    first residue of a helix-following turn and never two in a row.
    if do_proline:
        last = -9
        for i in range(2, n):
            if i in protected or ss[i - 1] != "L":
                continue
            if out[i - 1] in "PGC" or i - last < 3:
                continue
            if PRO_PHI[0] <= phi[i - 1] <= PRO_PHI[1] and ss[i - 2] == "L":
                out[i - 1] = "P"
                notes.append(f"{sequence[i - 1]}{i}P (phi {phi[i - 1]:.0f})")
                last = i

    # 3. At most one disulfide staple between two non-shell positions.
    if do_staple:
        cb = _coords(residues, "CB")
        best = None
        for i in range(1, n + 1):
            for j in range(i + 5, n + 1):
                if i in protected or j in protected:
                    continue
                if cb[i - 1] is None or cb[j - 1] is None:
                    continue
                if out[i - 1] == "P" or out[j - 1] == "P":
                    continue
                distance = float(np.linalg.norm(cb[i - 1] - cb[j - 1]))
                if SS_PAIR[0] <= distance <= SS_PAIR[1]:
                    if best is None or distance < best[2]:
                        best = (i, j, distance)
        if best:
            i, j, distance = best
            out[i - 1] = out[j - 1] = "C"
            notes.append(f"staple C{i}-C{j} (CB-CB {distance:.2f} A)")
    return "".join(out), notes
