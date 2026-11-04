#!/usr/bin/env python3
"""Where on the groove a paratope sits: alpha1/alpha2 helix vs peptide.

Why this module exists
----------------------
Batch-06 measured that shrinking a binder scales the whole interface down and
moves the peptide share of buried area by 0.004 (0.209 -> 0.213); batch-07
measured that rigidifying it moves nothing at all, and that the HLA groove is
still ~75 % of bb54-g2's interface (peptide 208 A^2 vs HLA 613 A^2).  Both
batches partitioned buried area only as *peptide vs HLA chain*, which cannot
tell "lying along the alpha1/alpha2 helices" apart from "resting on the groove
rim beside the peptide".  Batch-08 attacks interface **shape**, so it needs the
finer partition and an explicit read-out of how the binder *approaches* the
groove:

* :func:`buried_area_partition` - buried area split peptide / alpha1 helix /
  alpha2 helix / rest of the HLA heavy chain / beta-2m, from the same
  Shrake-Rupley SASA difference `complex_analysis.interface_metrics` uses, so
  the totals are directly comparable with every earlier batch.
* :func:`contact_partition` - the same split on heavy-atom contact counts, and
  the version of it that works on an RFdiffusion output (binder chain A,
  concatenated target chain B) *before* any Boltz-2 call is paid for.
* :func:`groove_frame` / :func:`approach_metrics` - a local frame on the
  groove (peptide long axis, groove normal "up", helix-to-helix "lateral") and
  where the binder's centre of mass sits in it.  A paratope that comes down on
  the peptide from above has a large positive height and a small lateral
  offset; one that drapes along a helix has a large lateral offset.

Nothing here is a gate.  Thresholds live in the batch-08 driver and are stated
before the batch is run; this module only measures.

Helix definitions
-----------------
`ALPHA1_HELIX` = HLA heavy-chain residues 57-84, `ALPHA2_HELIX` = 138-180, the
canonical class-I groove helices.  Checked against 7OW6 itself with a
CA(i)-CA(i+4) < 6.5 A helicity test on the trimmed target, which returns
helical segments 57-65, 76-82 (alpha1) and 140-147, 156-165, 167-172
(alpha2) - i.e. the canonical ranges bracket every helical segment that lines
the groove, and the 66-75 kink and 148-155 turn are included deliberately
because they are part of the same two walls the paratope has to avoid.
"""
from __future__ import annotations

from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np

import complex_analysis as ca

ALPHA1_HELIX = (57, 84)
ALPHA2_HELIX = (138, 180)
# The beta-sheet floor of the groove, used only to orient the "up" vector.
GROOVE_FLOOR_RANGES = ((1, 56), (85, 137))
PEPTIDE_LENGTH = 10


def _in(resseq: int, span: Tuple[int, int]) -> bool:
    return span[0] <= resseq <= span[1]


def segment_of(resseq: int) -> str:
    """Which part of the HLA heavy chain residue `resseq` belongs to."""
    if _in(resseq, ALPHA1_HELIX):
        return "hla_alpha1"
    if _in(resseq, ALPHA2_HELIX):
        return "hla_alpha2"
    return "hla_other"


# --- buried-area partition on a predicted complex ---------------------------

def buried_area_partition(text: str, binder_sequence: str,
                          arm: str = "g12d", fmt: str = "mmcif") -> Dict[str, Any]:
    """Buried area split peptide / alpha1 / alpha2 / rest-of-HLA / b2m.

    Same Shrake-Rupley difference as `complex_analysis.interface_metrics`
    (1.4 A probe, binder present vs binder chain detached, same coordinates),
    so `buried_area_total` here equals the value the other scripts report.
    Works on the apo arm too, where there is no peptide chain.
    """
    model = ca._clean_model(next(iter(ca.parse_structure(text, fmt=fmt))))
    chain_map = ca.map_complex_chains(
        model, binder_sequence=binder_sequence, binder_chain_hint="D",
        require_peptide=(arm != "apo"))
    binder, peptide, hla, b2m = (chain_map.binder, chain_map.peptide,
                                 chain_map.hla, chain_map.b2m)
    bound = ca._atom_sasa(model)
    free_model = ca._clean_model(next(iter(ca.parse_structure(text, fmt=fmt))))
    free_model.detach_child(binder)
    free = ca._atom_sasa(free_model)

    parts: Dict[str, float] = {"peptide": 0.0, "hla_alpha1": 0.0,
                               "hla_alpha2": 0.0, "hla_other": 0.0, "b2m": 0.0,
                               "unassigned": 0.0}
    total = 0.0
    for (chain_id, resseq, _atom), free_value in free.items():
        delta = free_value - bound.get((chain_id, resseq, _atom), 0.0)
        if delta <= 0:
            continue
        total += delta
        if peptide is not None and chain_id == peptide:
            parts["peptide"] += delta
        elif hla is not None and chain_id == hla:
            parts[segment_of(resseq)] += delta
        elif b2m is not None and chain_id == b2m:
            parts["b2m"] += delta
        else:
            parts["unassigned"] += delta

    helix = parts["hla_alpha1"] + parts["hla_alpha2"]
    out: Dict[str, Any] = {
        "arm": arm,
        "chains": {"binder": binder, "peptide": peptide, "hla": hla, "b2m": b2m},
        "buried_area_total": round(total, 1),
        "buried_area_peptide": round(parts["peptide"], 1),
        "buried_area_hla_alpha1": round(parts["hla_alpha1"], 1),
        "buried_area_hla_alpha2": round(parts["hla_alpha2"], 1),
        "buried_area_hla_helices": round(helix, 1),
        "buried_area_hla_other": round(parts["hla_other"], 1),
        "buried_area_b2m": round(parts["b2m"], 1),
        "helix_area_share": round(helix / total, 3) if total > 0 else None,
        "peptide_area_share": (round(parts["peptide"] / total, 3)
                               if total > 0 else None),
    }
    return out


# --- contact partition ------------------------------------------------------

def _contacts(search, residues) -> int:
    hits = set()
    for residue in residues:
        for atom in ca.heavy_atoms(residue):
            for neighbour in search.search(atom.coord, ca.CONTACT_CUTOFF):
                hits.add(id(neighbour))
    return len(hits)


def design_shape_metrics(pdb_text: str, binder_chain: str = "A",
                         target_chain: str = "B",
                         peptide_length: int = PEPTIDE_LENGTH) -> Dict[str, Any]:
    """Shape read-out of an **RFdiffusion output** backbone, before any scoring.

    The NIM returns the designed binder as chain A and the kept target segments
    concatenated into chain B in contig order, so for the contig
    `A1-180/0 C1-10/0 L-L` chain B is HLA residues 1-180 followed by the
    10-mer.  The HLA segment is contiguous and gap-free in 7OW6 (checked), so
    the i-th residue of that segment is HLA residue i+1 and the helix spans
    above apply by index.

    Returns the alpha1/alpha2 contact split, the peptide contact fraction that
    batch-06 filtered on, and the approach geometry of :func:`approach_metrics`.
    """
    from Bio.PDB import NeighborSearch

    model = ca._clean_model(next(iter(ca.parse_structure(pdb_text, fmt="pdb"))))
    chains = {c.id: c for c in model}
    for needed in (binder_chain, target_chain):
        if needed not in chains:
            raise ca.StructureAnalysisError(
                f"chain {needed!r} missing from the designed complex "
                f"(chains {sorted(chains)})")
    target = [r for r in chains[target_chain] if ca.is_polymer_residue(r)]
    if len(target) <= peptide_length:
        raise ca.StructureAnalysisError(
            f"target chain {target_chain} has {len(target)} residues, which "
            f"cannot be the groove plus a {peptide_length}-mer")
    peptide = target[-peptide_length:]
    hla = target[:-peptide_length]
    binder_atoms = ca.heavy_atoms(chains[binder_chain])
    if not binder_atoms:
        raise ca.StructureAnalysisError("binder chain has no heavy atoms")
    search = NeighborSearch(binder_atoms)

    alpha1 = [r for i, r in enumerate(hla, start=1) if _in(i, ALPHA1_HELIX)]
    alpha2 = [r for i, r in enumerate(hla, start=1) if _in(i, ALPHA2_HELIX)]
    other = [r for i, r in enumerate(hla, start=1)
             if not (_in(i, ALPHA1_HELIX) or _in(i, ALPHA2_HELIX))]

    n_pep = _contacts(search, peptide)
    n_a1 = _contacts(search, alpha1)
    n_a2 = _contacts(search, alpha2)
    n_other = _contacts(search, other)
    n_target = n_pep + n_a1 + n_a2 + n_other

    frame = groove_frame([r["CA"].coord for r in hla if "CA" in r],
                         [r["CA"].coord for r in peptide if "CA" in r],
                         hla_resseqs=[i for i, _ in enumerate(hla, start=1)])
    approach = approach_metrics(
        [a.coord for a in binder_atoms if a.get_id() == "CA"], frame)

    return {
        "peptide_contacts": n_pep,
        "alpha1_contacts": n_a1,
        "alpha2_contacts": n_a2,
        "helix_contacts": n_a1 + n_a2,
        "hla_other_contacts": n_other,
        "target_contacts": n_target,
        "peptide_contact_frac": round(n_pep / n_target, 3) if n_target else 0.0,
        "helix_contact_frac": round((n_a1 + n_a2) / n_target, 3) if n_target else 0.0,
        **approach,
    }


# --- the groove frame -------------------------------------------------------

def groove_frame(hla_ca: Sequence, peptide_ca: Sequence,
                 hla_resseqs: Optional[Sequence[int]] = None) -> Dict[str, Any]:
    """Local orthonormal frame on the groove.

    * `axis` - the peptide long axis (first principal component of its CAs),
      N-to-C oriented.
    * `up` - the groove normal, from the beta-sheet floor centroid towards the
      peptide centroid, orthogonalised against `axis`.  Positive `up` is out of
      the groove, i.e. where a TCR (or a binder) comes from.
    * `lateral` - `up` x `axis`, which points from one helix wall to the other.
    """
    hla_ca = np.asarray(hla_ca, dtype=float)
    peptide_ca = np.asarray(peptide_ca, dtype=float)
    if len(peptide_ca) < 3 or len(hla_ca) < 10:
        raise ca.StructureAnalysisError("not enough CA atoms to build a frame")
    if hla_resseqs is None:
        hla_resseqs = list(range(1, len(hla_ca) + 1))
    floor = np.array([c for c, n in zip(hla_ca, hla_resseqs)
                      if any(_in(n, span) for span in GROOVE_FLOOR_RANGES)])
    if len(floor) < 10:
        raise ca.StructureAnalysisError("groove floor not found in this chain")

    centre = peptide_ca.mean(axis=0)
    centred = peptide_ca - centre
    axis = np.linalg.svd(centred, full_matrices=False)[2][0]
    if np.dot(peptide_ca[-1] - peptide_ca[0], axis) < 0:
        axis = -axis
    up = centre - floor.mean(axis=0)
    up = up - np.dot(up, axis) * axis
    up /= np.linalg.norm(up)
    lateral = np.cross(up, axis)
    lateral /= np.linalg.norm(lateral)
    return {"centre": centre, "axis": axis, "up": up, "lateral": lateral}


def approach_metrics(binder_ca: Sequence, frame: Dict[str, Any],
                     over_peptide_lateral: float = 7.0) -> Dict[str, Any]:
    """Where the binder sits in the groove frame.

    `approach_height` is the mean height of the binder CAs above the peptide
    centroid along the groove normal - how far "above the groove axis" the
    paratope sits.  `lateral_offset` is the mean *signed* displacement towards
    one helix wall and `abs_lateral_offset` its magnitude; a binder laid along
    alpha1 or alpha2 has a large magnitude.  `over_peptide_frac` is the share
    of binder CAs within `over_peptide_lateral` A of the groove's long axis,
    i.e. genuinely over the peptide rather than over a wall.
    """
    binder_ca = np.asarray(list(binder_ca), dtype=float)
    if len(binder_ca) == 0:
        raise ca.StructureAnalysisError("binder has no CA atoms")
    delta = binder_ca - frame["centre"]
    height = delta @ frame["up"]
    lateral = delta @ frame["lateral"]
    along = delta @ frame["axis"]
    return {
        "approach_height": round(float(height.mean()), 2),
        "min_height": round(float(height.min()), 2),
        "lateral_offset": round(float(lateral.mean()), 2),
        "abs_lateral_offset": round(float(abs(lateral.mean())), 2),
        "lateral_spread": round(float(lateral.std()), 2),
        "along_axis_offset": round(float(along.mean()), 2),
        "over_peptide_frac": round(
            float(np.mean(np.abs(lateral) <= over_peptide_lateral)), 3),
    }


def predicted_approach(text: str, binder_sequence: str, arm: str = "g12d",
                       fmt: str = "mmcif") -> Dict[str, Any]:
    """:func:`approach_metrics` on a predicted (Boltz-2) complex."""
    model = ca._clean_model(next(iter(ca.parse_structure(text, fmt=fmt))))
    chain_map = ca.map_complex_chains(
        model, binder_sequence=binder_sequence, binder_chain_hint="D",
        require_peptide=(arm != "apo"))
    if chain_map.peptide is None or chain_map.hla is None:
        return {}
    hla = [r for r in model[chain_map.hla] if ca.is_polymer_residue(r)]
    peptide = [r for r in model[chain_map.peptide] if ca.is_polymer_residue(r)]
    frame = groove_frame([r["CA"].coord for r in hla if "CA" in r],
                         [r["CA"].coord for r in peptide if "CA" in r],
                         hla_resseqs=[r.id[1] for r in hla if "CA" in r])
    binder = [r for r in model[chain_map.binder] if ca.is_polymer_residue(r)]
    return approach_metrics([r["CA"].coord for r in binder if "CA" in r], frame)


def helix_facing_positions(pdb_text: str, binder_chain: str = "A",
                           target_chain: str = "B",
                           peptide_length: int = PEPTIDE_LENGTH,
                           radius: float = 6.0) -> List[Dict[str, Any]]:
    """Binder positions (1-based) whose heavy atoms sit against alpha1/alpha2.

    Measured on the RFdiffusion backbone, which is backbone-only on the target
    side, so this is deliberately coarse - it is used to *select against* helix
    contact in sequence design, not to score an interface.
    """
    model = ca._clean_model(next(iter(ca.parse_structure(pdb_text, fmt="pdb"))))
    chains = {c.id: c for c in model}
    target = [r for r in chains[target_chain] if ca.is_polymer_residue(r)]
    hla = target[:-peptide_length]
    peptide = target[-peptide_length:]
    helix_atoms = np.array(
        [a.coord for i, r in enumerate(hla, start=1)
         if _in(i, ALPHA1_HELIX) or _in(i, ALPHA2_HELIX)
         for a in ca.heavy_atoms(r)])
    peptide_atoms = np.array([a.coord for r in peptide for a in ca.heavy_atoms(r)])
    out: List[Dict[str, Any]] = []
    binder = [r for r in chains[binder_chain] if ca.is_polymer_residue(r)]
    for index, residue in enumerate(binder, start=1):
        atoms = np.array([a.coord for a in ca.heavy_atoms(residue)])
        if not len(atoms):
            continue
        d_helix = float(np.min(np.linalg.norm(
            atoms[:, None, :] - helix_atoms[None, :, :], axis=-1)))
        d_peptide = float(np.min(np.linalg.norm(
            atoms[:, None, :] - peptide_atoms[None, :, :], axis=-1)))
        if d_helix <= radius and d_helix < d_peptide:
            out.append({"position": index, "helix_distance": round(d_helix, 2),
                        "peptide_distance": round(d_peptide, 2)})
    out.sort(key=lambda r: r["helix_distance"])
    return out
