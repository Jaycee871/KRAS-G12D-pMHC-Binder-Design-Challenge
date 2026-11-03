#!/usr/bin/env python3
"""Pose filters: is the binder actually sitting on the peptide?

Why this module exists
----------------------
Batch-01 spent four Boltz-2 calls on every ProteinMPNN sequence and 45 of 56
designs turned out to bury no peptide surface at all. Two different things
were going wrong, and they need two different filters:

1. **The designed backbone never pointed at the peptide.** RFdiffusion is told
   to use peptide hotspots, but nothing checks the backbone it returns. A
   backbone whose binder chain sits on the alpha1/alpha2 helices cannot produce
   a specific sequence, and it is free to detect: the designed complex is right
   there in the RFdiffusion output PDB. :func:`design_pose_metrics` /
   :func:`design_pose_ok` do that check *before any Boltz-2 call is made*.

2. **Boltz-2 re-docked the binder somewhere else.** Boltz-2 predicts from
   sequence alone and routinely parks the binder on the HLA heavy chain or the
   beta-2m face. Once the G12D prediction is in hand that is visible from the
   interface metrics, so the remaining three states (WT, apo, decoy) can be
   skipped: :func:`predicted_pose_ok`. That is three of four calls saved on
   every off-peptide prediction, which is where most of the batch-02 budget
   saving comes from.

Geometry notes
--------------
In an RFdiffusion NIM output the designed binder is chain ``A`` and every kept
target segment is concatenated into chain ``B`` in contig order, so for the
contig ``A1-180/0 C1-10/0 L-L`` the last ten residues of chain B are the
peptide and its position 6 is the G12D aspartate. The target is returned
**backbone only** (N, CA, C, O - no side chains), so this filter is
deliberately coarse: heavy-atom contact counts against the peptide backbone
plus the closest approach to position 6. It is a *cheap triage*, not the
interface analysis in :mod:`complex_analysis`.
"""
from __future__ import annotations

from typing import Any, Dict, List, Mapping, Optional, Tuple

import complex_analysis
from complex_analysis import StructureAnalysisError

# --- thresholds -------------------------------------------------------------
# Designed (RFdiffusion) complex. Calibrated on the 32 batch-01 backbones:
# see results/batch-02/SUMMARY.md.
CONTACT_CUTOFF = complex_analysis.CONTACT_CUTOFF        # 5.0 A, heavy atom
DESIGN_MIN_PEPTIDE_CONTACTS = 10       # binder heavy atoms within 5 A of the peptide
DESIGN_MIN_PEPTIDE_RESIDUES = 3        # distinct peptide residues contacted
DESIGN_MIN_POS6_CONTACTS = 1           # binder heavy atoms within 5 A of position 6
DESIGN_MAX_POS6_DISTANCE = 6.0         # A, closest binder heavy atom to position 6

# Predicted (Boltz-2) G12D complex, from complex_analysis.interface_metrics.
PREDICTED_MIN_PEPTIDE_AREA = 25.0      # A^2 buried on the peptide chain
PREDICTED_MIN_POS6_CONTACTS = 1        # at least one heavy atom near position 6

PEPTIDE_LENGTH = 10
MUTATED_PEPTIDE_POSITION = complex_analysis.MUTATED_PEPTIDE_POSITION
MUTATED_RESIDUE_NAME = complex_analysis.MUTATED_RESIDUE_NAME


def _heavy(atoms) -> List:
    return [a for a in atoms if a.element != "H"]


def design_pose_metrics(
    pdb_text: str,
    binder_chain: str = "A",
    target_chain: str = "B",
    peptide_length: int = PEPTIDE_LENGTH,
    peptide_position: int = MUTATED_PEPTIDE_POSITION,
    expected_residue: Optional[str] = MUTATED_RESIDUE_NAME,
) -> Dict[str, Any]:
    """Where the *designed* binder sits on the RFdiffusion output complex.

    Returns heavy-atom contact counts of the binder against the peptide
    segment and against the rest of the target, the number of distinct peptide
    residues contacted, and the closest approach to the mutated position.

    Raises `StructureAnalysisError` when the structure is not the expected
    binder + concatenated-target pair, or when the residue found at the
    mutated position is not the expected one - never returns a quiet zero.
    """
    from Bio.PDB import NeighborSearch

    structure = complex_analysis.parse_structure(pdb_text, fmt="pdb")
    model = complex_analysis._clean_model(next(iter(structure)))
    chains = {c.id: c for c in model}
    if binder_chain not in chains:
        raise StructureAnalysisError(
            f"binder chain {binder_chain!r} is not in the designed complex "
            f"(chains: {sorted(chains)})"
        )
    if target_chain not in chains:
        raise StructureAnalysisError(
            f"target chain {target_chain!r} is not in the designed complex "
            f"(chains: {sorted(chains)})"
        )

    target_residues = [r for r in chains[target_chain]
                       if complex_analysis.is_polymer_residue(r)]
    if len(target_residues) <= peptide_length:
        raise StructureAnalysisError(
            f"target chain {target_chain} has only {len(target_residues)} "
            f"residues, which cannot be the groove plus a {peptide_length}-mer "
            f"peptide"
        )
    peptide_residues = target_residues[-peptide_length:]
    other_residues = target_residues[:-peptide_length]

    pos6 = peptide_residues[peptide_position - 1]
    resname = pos6.get_resname().strip().upper()
    if expected_residue and resname != expected_residue.upper():
        raise StructureAnalysisError(
            f"position {peptide_position} of the peptide segment is {resname}, "
            f"expected {expected_residue.upper()}; the contig order or the "
            f"peptide length assumption is wrong (segment: "
            f"{''.join(complex_analysis.one_letter(r) for r in peptide_residues)})"
        )

    binder_atoms = _heavy(chains[binder_chain].get_atoms())
    if not binder_atoms:
        raise StructureAnalysisError(
            f"binder chain {binder_chain} has no heavy atoms")
    search = NeighborSearch(binder_atoms)

    def contacts(residues) -> Tuple[int, int]:
        """(binder heavy atoms in contact, residues contacted)."""
        hits = set()
        touched = 0
        for residue in residues:
            before = len(hits)
            for atom in _heavy(residue):
                for neighbour in search.search(atom.coord, CONTACT_CUTOFF):
                    hits.add(id(neighbour))
            if len(hits) > before:
                touched += 1
        return len(hits), touched

    peptide_contacts, peptide_residues_contacted = contacts(peptide_residues)
    hla_contacts, _ = contacts(other_residues)
    pos6_contacts, _ = contacts([pos6])

    pos6_distance = min(
        float(((atom.coord - other.coord) ** 2).sum() ** 0.5)
        for atom in _heavy(pos6) for other in binder_atoms
    )

    total = peptide_contacts + hla_contacts
    return {
        "binder_length": len([r for r in chains[binder_chain]
                              if complex_analysis.is_polymer_residue(r)]),
        "peptide_sequence": "".join(
            complex_analysis.one_letter(r) for r in peptide_residues),
        "peptide_contacts": peptide_contacts,
        "peptide_residues_contacted": peptide_residues_contacted,
        "hla_contacts": hla_contacts,
        "pos6_contacts": pos6_contacts,
        "pos6_distance": round(pos6_distance, 2),
        "peptide_contact_frac": (peptide_contacts / total) if total else 0.0,
    }


def design_pose_ok(
    metrics: Mapping[str, Any],
    min_peptide_contacts: int = DESIGN_MIN_PEPTIDE_CONTACTS,
    min_peptide_residues: int = DESIGN_MIN_PEPTIDE_RESIDUES,
    min_pos6_contacts: int = DESIGN_MIN_POS6_CONTACTS,
    max_pos6_distance: float = DESIGN_MAX_POS6_DISTANCE,
) -> Tuple[bool, str]:
    """Does this designed backbone sit on the peptide? `(ok, reason)`."""
    reasons: List[str] = []
    if metrics["peptide_contacts"] < min_peptide_contacts:
        reasons.append(
            f"only {metrics['peptide_contacts']} binder heavy atoms within "
            f"{CONTACT_CUTOFF:.0f} A of the peptide "
            f"(need >= {min_peptide_contacts})"
        )
    if metrics["peptide_residues_contacted"] < min_peptide_residues:
        reasons.append(
            f"only {metrics['peptide_residues_contacted']} peptide residues "
            f"contacted (need >= {min_peptide_residues})"
        )
    if metrics["pos6_contacts"] < min_pos6_contacts:
        reasons.append(
            f"only {metrics['pos6_contacts']} binder heavy atoms within "
            f"{CONTACT_CUTOFF:.0f} A of position 6 "
            f"(need >= {min_pos6_contacts})"
        )
    if metrics["pos6_distance"] > max_pos6_distance:
        reasons.append(
            f"closest approach to position 6 is "
            f"{metrics['pos6_distance']:.1f} A (need <= "
            f"{max_pos6_distance:.1f} A)"
        )
    return (not reasons), "; ".join(reasons)


def predicted_pose_ok(
    metrics: Mapping[str, Any],
    min_peptide_area: float = PREDICTED_MIN_PEPTIDE_AREA,
    min_pos6_contacts: int = PREDICTED_MIN_POS6_CONTACTS,
) -> Tuple[bool, str]:
    """Did Boltz-2 keep the binder on the peptide? `(ok, reason)`.

    Takes the output of :func:`complex_analysis.interface_metrics` on the G12D
    prediction. A design that fails here is not worth the WT, apo and decoy
    predictions: whatever its ipTM says, it is not reading the mutation.
    """
    if not metrics:
        return False, "no interface metrics (analysis failed)"
    reasons: List[str] = []
    area = float(metrics.get("buried_area_peptide") or 0.0)
    contacts = float(metrics.get("pos6_contacts") or 0.0)
    if area < min_peptide_area:
        reasons.append(
            f"only {area:.0f} A^2 buried on the peptide "
            f"(need >= {min_peptide_area:.0f} A^2)"
        )
    if contacts < min_pos6_contacts:
        reasons.append(
            f"only {contacts:.0f} binder heavy atoms within "
            f"{CONTACT_CUTOFF:.0f} A of peptide position 6 "
            f"(need >= {min_pos6_contacts})"
        )
    return (not reasons), "; ".join(reasons)
