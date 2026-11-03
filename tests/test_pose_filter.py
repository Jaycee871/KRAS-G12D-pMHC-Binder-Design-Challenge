"""Unit tests for the batch-02 pose filters.

Batch-01's core finding was that the binders were not sitting on the peptide:
45 of 56 scored designs buried no peptide surface at all, after four Boltz-2
calls each. These tests pin the two filters that stop that happening again -
one on the designed RFdiffusion complex (free), one on the Boltz-2 G12D
prediction (saves the other three states).

Nothing here makes an API call.
"""
from __future__ import annotations

import pytest

import pose_filter
from complex_analysis import StructureAnalysisError

PEPTIDE = ["VAL", "VAL", "VAL", "GLY", "ALA", "ASP", "GLY", "VAL", "GLY", "LYS"]


def backbone_pdb(binder_offset: float, peptide=PEPTIDE, binder_residues=20,
                 groove_residues=15) -> str:
    """An RFdiffusion-shaped output PDB.

    Chain A is the designed binder; chain B is the concatenated target, the
    groove first and the peptide last, exactly as the NIM returns it. The
    binder is a straight chain `binder_offset` Angstrom above the peptide, so
    the test can place it on the peptide or far away from it.
    """
    lines = []
    serial = 1

    def atoms(resname, chain, resseq, x, y):
        nonlocal serial
        for name in ("N", "CA", "C", "O"):
            lines.append(
                f"ATOM  {serial:>5} {name:<4} {resname:>3} {chain}{resseq:>4}"
                f"    {x:>8.3f}{y:>8.3f}{0.0:>8.3f}  1.00  0.00"
                f"           {name[0]}"
            )
            serial += 1

    for index in range(binder_residues):
        atoms("GLY", "A", index + 1, index * 3.8, binder_offset)
    resseq = binder_residues + 1
    for index in range(groove_residues):
        atoms("ALA", "B", resseq, index * 3.8, 20.0)
        resseq += 1
    for index, resname in enumerate(peptide):
        atoms(resname, "B", resseq, index * 3.8, 0.0)
        resseq += 1
    return "\n".join(lines) + "\nTER\nEND\n"


# --- filter 1: the designed complex -----------------------------------------

def test_a_binder_on_the_peptide_passes():
    metrics = pose_filter.design_pose_metrics(backbone_pdb(3.5))
    assert metrics["peptide_sequence"] == "VVVGADGVGK"
    assert metrics["peptide_contacts"] >= pose_filter.DESIGN_MIN_PEPTIDE_CONTACTS
    assert metrics["pos6_contacts"] >= 1
    ok, reason = pose_filter.design_pose_ok(metrics)
    assert ok and reason == ""


def test_a_binder_parked_away_from_the_peptide_is_rejected():
    metrics = pose_filter.design_pose_metrics(backbone_pdb(30.0))
    assert metrics["peptide_contacts"] == 0
    ok, reason = pose_filter.design_pose_ok(metrics)
    assert not ok
    assert "peptide" in reason


def test_a_binder_near_the_peptide_but_off_position_6_is_rejected():
    """Contacts the peptide terminus only: no read-out of the mutation."""
    metrics = pose_filter.design_pose_metrics(
        backbone_pdb(3.5, binder_residues=3))
    assert metrics["pos6_distance"] > pose_filter.DESIGN_MAX_POS6_DISTANCE
    ok, reason = pose_filter.design_pose_ok(metrics)
    assert not ok
    assert "position 6" in reason


def test_the_mutated_position_is_checked_not_assumed():
    """Wild-type peptide in the target segment: the filter must refuse it."""
    wt = list(PEPTIDE)
    wt[5] = "GLY"
    with pytest.raises(StructureAnalysisError, match="expected ASP"):
        pose_filter.design_pose_metrics(backbone_pdb(3.5, peptide=wt))


def test_a_missing_binder_chain_raises():
    pdb = "\n".join(line for line in backbone_pdb(3.5).splitlines()
                    if line[21:22] != "A")
    with pytest.raises(StructureAnalysisError, match="binder chain"):
        pose_filter.design_pose_metrics(pdb)


def test_thresholds_are_explicit():
    assert pose_filter.DESIGN_MIN_PEPTIDE_CONTACTS == 10
    assert pose_filter.DESIGN_MIN_PEPTIDE_RESIDUES == 3
    assert pose_filter.DESIGN_MIN_POS6_CONTACTS == 1
    assert pose_filter.DESIGN_MAX_POS6_DISTANCE == 6.0


# --- filter 2: the Boltz-2 prediction ---------------------------------------

def test_a_prediction_on_the_peptide_survives():
    ok, reason = pose_filter.predicted_pose_ok(
        {"buried_area_peptide": 215.9, "pos6_contacts": 9})
    assert ok and reason == ""


def test_a_prediction_parked_on_the_hla_is_rejected():
    ok, reason = pose_filter.predicted_pose_ok(
        {"buried_area_peptide": 0.0, "pos6_contacts": 0,
         "buried_area_total": 780.0})
    assert not ok
    assert "buried on the peptide" in reason
    assert "position 6" in reason


def test_a_prediction_with_a_sliver_of_peptide_area_is_rejected():
    """Batch-01's mean was 21 A^2 on the peptide: not a paratope."""
    ok, _ = pose_filter.predicted_pose_ok(
        {"buried_area_peptide": 21.0, "pos6_contacts": 2})
    assert not ok


def test_failed_analysis_is_not_a_pass():
    ok, reason = pose_filter.predicted_pose_ok({})
    assert not ok
    assert "analysis failed" in reason
