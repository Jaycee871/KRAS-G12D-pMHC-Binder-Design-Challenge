"""Unit tests for src/ddg_scoring.py - the paired same-pose ddG.

These cover the bookkeeping the metric stands on, not the force field:

* the two states of a pair must contain **exactly the same residues**, because
  the whole point of a paired calculation is that everything except peptide
  position 6 cancels;
* the WT state must be a real glycine - resname GLY, side chain gone, backbone
  untouched - and the G12D state must still carry the carboxylate;
* a binder that does not reach position 6 must *raise*, not quietly return
  ddg 0, since "no signal" and "nothing measured" are the exact pair of
  outcomes batch-01 could not tell apart;
* the two-chain binder form (the 7OW6 TCR) must be expressible.

The energy terms themselves are validated against the 7OW6 control in
`results/ddg-scoring/FINDINGS.md`, not here: a force-field evaluation is far too
slow for a unit test.
"""
from __future__ import annotations

import pytest

import complex_analysis as ca
import ddg_scoring as dd

from test_complex_analysis import BINDER_SEQUENCE, bound_complex, make_cif, \
    HLA_SEQUENCE, B2M_SEQUENCE, model_of


def groups_of(text, **kwargs):
    model = model_of(text)
    return model, dd.groups_from_model(model, **kwargs)


# --- group identification ---------------------------------------------------

def test_groups_from_model_maps_a_single_chain_binder_by_sequence():
    _model, groups = groups_of(bound_complex(), binder_sequence=BINDER_SEQUENCE)
    assert groups.binder_chains == ("D",)
    assert groups.peptide_chain == "C"
    assert set(groups.target_chains) == {"A", "B", "C"}


def test_groups_from_model_accepts_a_two_chain_binder():
    """The 7OW6 control's binder is the TCR alpha/beta pair, not one chain."""
    text = make_cif([
        ("A", HLA_SEQUENCE, -14.0),
        ("B", B2M_SEQUENCE, -24.0),
        ("C", ca.G12D_PEPTIDE, 0.0),
        ("D", BINDER_SEQUENCE, 7.0),
        ("E", BINDER_SEQUENCE, 12.0),
    ])
    _model, groups = groups_of(text, binder_chains=["D", "E"])
    assert groups.binder_chains == ("D", "E")
    assert groups.peptide_chain == "C"
    assert set(groups.target_chains) == {"A", "B", "C"}


def test_groups_from_model_rejects_an_absent_binder_chain():
    with pytest.raises(dd.DdgError):
        groups_of(bound_complex(), binder_chains=["Z"])


# --- truncation and the in-place mutation -----------------------------------

def test_both_states_contain_the_same_residues():
    model, groups = groups_of(bound_complex(), binder_sequence=BINDER_SEQUENCE)
    g12d = dd.build_local_system(model, groups)
    wt = dd.build_local_system(model, groups, mutate_site_to_gly=True)
    assert g12d.n_residues == wt.n_residues
    assert g12d.n_fragments == wt.n_fragments
    assert g12d.binder_new_chains == wt.binder_new_chains
    assert (g12d.site_new_chain, g12d.site_resseq) == \
           (wt.site_new_chain, wt.site_resseq)


def _site_residue(local):
    chain = local.structure[0][local.site_new_chain]
    return chain[(" ", local.site_resseq, " ")]


def test_g12d_state_keeps_the_carboxylate_and_wt_state_is_a_glycine():
    model, groups = groups_of(bound_complex(), binder_sequence=BINDER_SEQUENCE)
    g12d = _site_residue(dd.build_local_system(model, groups))
    wt = _site_residue(
        dd.build_local_system(model, groups, mutate_site_to_gly=True))

    assert g12d.get_resname() == "ASP"
    assert {"OD1", "OD2"} <= {a.get_name() for a in g12d}

    assert wt.get_resname() == "GLY"
    assert {a.get_name() for a in wt} == {"N", "CA", "C", "O"}
    # the backbone is untouched: a paired calculation may not move the peptide
    for name in ("N", "CA", "C", "O"):
        assert g12d[name].get_coord() == pytest.approx(wt[name].get_coord())


def test_truncation_drops_distant_residues_but_keeps_the_site():
    model, groups = groups_of(bound_complex(), binder_sequence=BINDER_SEQUENCE)
    full = sum(len(dd.polymer_residues(c)) for c in model)
    local = dd.build_local_system(model, groups, radius=8.0)
    assert local.n_residues < full
    assert _site_residue(local).get_resname() == "ASP"


def test_a_binder_that_does_not_reach_position_six_raises():
    """`no signal` and `nothing measured` must never look the same."""
    model, groups = groups_of(bound_complex(binder_y=60.0),
                              binder_sequence=BINDER_SEQUENCE)
    with pytest.raises(dd.DdgError) as excinfo:
        dd.build_local_system(model, groups)
    assert "does not touch the site" in str(excinfo.value)


def test_an_out_of_range_peptide_position_raises():
    model, groups = groups_of(bound_complex(), binder_sequence=BINDER_SEQUENCE)
    groups.peptide_position = 99
    with pytest.raises(dd.DdgError):
        dd.build_local_system(model, groups)


# --- conventions the write-up depends on ------------------------------------

def test_sign_convention_is_documented_as_negative_favours_g12d():
    assert "negative" in dd.paired_ddg.__doc__.lower()
    assert "ddG_sel" in dd.__doc__ or "ddg_sel" in dd.__doc__


def test_restraint_shell_is_inside_the_truncation_sphere():
    assert dd.FREE_RADIUS < dd.LOCAL_RADIUS
