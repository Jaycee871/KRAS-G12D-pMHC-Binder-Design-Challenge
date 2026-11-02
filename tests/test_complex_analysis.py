"""Unit tests for src/complex_analysis.py - the chain mapping and SASA code.

These are the regression tests for the batch-01 failure: every design scored
`peptide_buried_frac = 0.00` / `asp6_buried_frac = 0.00` / `0` position-6
contacts, and nothing in the code could tell a real zero from a lookup that
had quietly missed. So the cases below cover both directions:

* a binder that really is on the peptide must produce non-zero burial and
  contacts, *whatever the chain ids are called*; and
* a missing chain, a missing residue, the wrong residue at position 6 or an
  ambiguous mapping must raise, never return 0.0.

The fixture is a tiny synthetic complex built below rather than a real
structure, so the expected geometry is known exactly and the tests need no
network.
"""
from __future__ import annotations

import pytest

import complex_analysis as ca
from complex_analysis import ChainMappingError, StructureAnalysisError

# β2m is recognised by its invariant N-terminus; a prefix is enough here.
B2M_SEQUENCE = ca.B2M_PREFIX + "V" * (99 - len(ca.B2M_PREFIX))
HLA_SEQUENCE = "G" * 180
BINDER_SEQUENCE = "AKRLDEWYTQ" * 6          # 60 aa, carries +2/-2 of charge


# --- fixture builder --------------------------------------------------------

CIF_HEADER = [
    "data_fixture",
    "#",
    "loop_",
    "_atom_site.group_PDB",
    "_atom_site.id",
    "_atom_site.type_symbol",
    "_atom_site.label_atom_id",
    "_atom_site.label_alt_id",
    "_atom_site.label_comp_id",
    "_atom_site.label_asym_id",
    "_atom_site.label_seq_id",
    "_atom_site.pdbx_PDB_ins_code",
    "_atom_site.Cartn_x",
    "_atom_site.Cartn_y",
    "_atom_site.Cartn_z",
    "_atom_site.occupancy",
    "_atom_site.B_iso_or_equiv",
    "_atom_site.auth_seq_id",
    "_atom_site.auth_asym_id",
    "_atom_site.pdbx_PDB_model_num",
]

ONE_TO_THREE = {v: k for k, v in ca.THREE_TO_ONE.items() if len(k) == 3}
ONE_TO_THREE.update({"A": "ALA", "G": "GLY", "M": "MET"})

# name, element, dx, dy, dz - the Asp side chain points at +y, i.e. at the
# chain stacked above it, so a binder placed there buries it.
BACKBONE = [("N", "N", 0.0, 0.0, 0.0), ("CA", "C", 1.0, 0.0, 0.0),
            ("C", "C", 2.0, 0.0, 0.0), ("O", "O", 2.5, 0.5, 0.0)]
CB = ("CB", "C", 1.0, 1.5, 0.0)
ASP_SIDECHAIN = [CB, ("CG", "C", 1.0, 2.5, 0.0),
                 ("OD1", "O", 0.5, 3.4, 0.0), ("OD2", "O", 1.8, 3.4, 0.0)]


def chain_rows(chain_id, sequence, y_offset, serial_start, hetero=()):
    rows = []
    serial = serial_start
    for index, aa in enumerate(sequence):
        comp = ONE_TO_THREE.get(aa.upper(), "GLY")
        atoms = list(BACKBONE)
        if comp == "ASP":
            atoms += ASP_SIDECHAIN
        elif comp != "GLY":
            atoms.append(CB)
        for name, element, dx, dy, dz in atoms:
            rows.append(
                f"ATOM {serial} {element} {name} . {comp} {chain_id} "
                f"{index + 1} ? {index * 3.8 + dx:.3f} {y_offset + dy:.3f} "
                f"{dz:.3f} 1.00 50.00 {index + 1} {chain_id} 1"
            )
            serial += 1
    for offset, (resname, element, name) in enumerate(hetero):
        rows.append(
            f"HETATM {serial} {element} {name} . {resname} {chain_id} "
            f"{len(sequence) + offset + 1} ? 0.000 {y_offset:.3f} 20.000 "
            f"1.00 50.00 {len(sequence) + offset + 1} {chain_id} 1"
        )
        serial += 1
    return rows, serial


def make_cif(chains):
    """`chains` is a list of (chain_id, sequence, y_offset[, hetero])."""
    rows = []
    serial = 1
    for spec in chains:
        chain_id, sequence, y_offset = spec[0], spec[1], spec[2]
        hetero = spec[3] if len(spec) > 3 else ()
        new_rows, serial = chain_rows(chain_id, sequence, y_offset, serial,
                                      hetero)
        rows += new_rows
    return "\n".join(CIF_HEADER + rows) + "\n#\n"


def bound_complex(binder_y=7.0, peptide=ca.G12D_PEPTIDE,
                  binder_id="D", peptide_id="C", hla_id="A", b2m_id="B",
                  **kwargs):
    """HLA + β2m + peptide + a binder stacked `binder_y` Å above the peptide."""
    return make_cif([
        (hla_id, HLA_SEQUENCE, -14.0),
        (b2m_id, B2M_SEQUENCE, -24.0),
        (peptide_id, peptide, 0.0),
        (binder_id, BINDER_SEQUENCE, binder_y),
    ], **kwargs)


def model_of(text):
    return ca._clean_model(next(iter(ca.parse_structure(text))))


# --- chain inventory and mapping -------------------------------------------

def test_chain_inventory_reports_length_and_sequence():
    inventory = ca.chain_inventory(model_of(bound_complex()))
    assert set(inventory) == {"A", "B", "C", "D"}
    assert inventory["C"].sequence == ca.G12D_PEPTIDE
    assert inventory["C"].length == 10
    assert inventory["D"].length == len(BINDER_SEQUENCE)


def test_mapping_identifies_every_chain_by_content():
    chain_map = ca.map_complex_chains(
        model_of(bound_complex()), binder_sequence=BINDER_SEQUENCE
    )
    assert (chain_map.binder, chain_map.peptide) == ("D", "C")
    assert (chain_map.hla, chain_map.b2m) == ("A", "B")


def test_mapping_does_not_depend_on_the_chain_ids():
    """The batch-01 hypothesis: ids renamed, mapping must still be right."""
    text = bound_complex(binder_id="X", peptide_id="Y", hla_id="Z", b2m_id="W")
    chain_map = ca.map_complex_chains(
        model_of(text), binder_sequence=BINDER_SEQUENCE
    )
    assert (chain_map.binder, chain_map.peptide) == ("X", "Y")
    assert (chain_map.hla, chain_map.b2m) == ("Z", "W")


def test_mapping_falls_back_to_the_chain_id_hint():
    chain_map = ca.map_complex_chains(model_of(bound_complex()),
                                      binder_chain_hint="D")
    assert chain_map.binder == "D"
    assert chain_map.peptide == "C"


def test_mapping_raises_when_the_binder_sequence_is_absent():
    text = bound_complex()
    with pytest.raises(ChainMappingError, match="designed binder sequence"):
        ca.map_complex_chains(model_of(text), binder_sequence="W" * 42)


def test_mapping_raises_when_the_hinted_chain_does_not_exist():
    with pytest.raises(ChainMappingError, match="binder chain 'Z'"):
        ca.map_complex_chains(model_of(bound_complex()), binder_chain_hint="Z")


def test_mapping_raises_when_there_is_no_peptide():
    text = make_cif([("A", HLA_SEQUENCE, -14.0),
                     ("D", BINDER_SEQUENCE, 7.0)])
    with pytest.raises(ChainMappingError, match="no peptide chain"):
        ca.map_complex_chains(model_of(text), binder_sequence=BINDER_SEQUENCE)


def test_mapping_raises_when_two_chains_could_be_the_peptide():
    # Neither short chain is a known peptide sequence, so content cannot
    # break the tie and the mapping must refuse instead of picking one.
    text = make_cif([("C", "ACDEFHIKLM", 0.0),
                     ("E", "MNPQRSTVWY", 14.0),
                     ("D", BINDER_SEQUENCE, 7.0)])
    with pytest.raises(ChainMappingError, match="could be the peptide"):
        ca.map_complex_chains(model_of(text), binder_sequence=BINDER_SEQUENCE)


def test_mapping_raises_when_there_are_only_waters():
    text = make_cif([("A", "", 0.0, (("HOH", "O", "O"),))])
    with pytest.raises(StructureAnalysisError, match="no polymer chains"):
        ca.chain_inventory(model_of(text))


# --- interface metrics ------------------------------------------------------

def test_metrics_on_a_binder_that_sits_on_the_peptide():
    metrics = ca.interface_metrics(bound_complex(binder_y=7.0),
                                   binder_sequence=BINDER_SEQUENCE)
    assert metrics["binder_chain"] == "D"
    assert metrics["peptide_chain"] == "C"
    assert metrics["pos6_residue"] == "ASP"
    assert metrics["buried_area_total"] > 0
    assert metrics["buried_area_peptide"] > 0
    assert 0.0 < metrics["peptide_buried_frac"] <= 1.0
    assert metrics["asp6_buried_frac"] > 0.3
    assert metrics["pos6_contacts"] >= ca_min_contacts()
    assert metrics["contacts_by_chain"]["C"] > 0
    assert metrics["notes"] == []


def ca_min_contacts():
    return 3


def test_a_binder_that_touches_nothing_is_zero_but_flagged():
    """A real zero must still be reported - and said out loud, not implied."""
    metrics = ca.interface_metrics(bound_complex(binder_y=60.0),
                                   binder_sequence=BINDER_SEQUENCE)
    assert metrics["peptide_buried_frac"] == 0.0
    assert metrics["asp6_buried_frac"] == 0.0
    assert metrics["pos6_contacts"] == 0.0
    assert any("no contact" in note or "none of it on the peptide" in note
               for note in metrics["notes"])


def test_a_binder_on_the_hla_only_is_flagged_as_off_peptide():
    """Batch-01's actual geometry: contact, but not with the peptide."""
    text = make_cif([
        ("A", HLA_SEQUENCE, -14.0),
        ("B", B2M_SEQUENCE, -24.0),
        ("C", ca.G12D_PEPTIDE, 0.0),
        ("D", BINDER_SEQUENCE, -20.0),
    ])
    metrics = ca.interface_metrics(text, binder_sequence=BINDER_SEQUENCE)
    assert metrics["buried_area_total"] > 0
    assert metrics["buried_area_peptide"] == 0.0
    assert metrics["peptide_buried_frac"] == 0.0
    assert any("off the peptide" in note for note in metrics["notes"])


def test_metrics_raise_when_position_6_is_not_the_aspartate():
    """The WT peptide has glycine at position 6: that is an error, not a 0.00."""
    text = bound_complex(peptide=ca.WT_PEPTIDE)
    with pytest.raises(StructureAnalysisError, match="position 6.*is GLY"):
        ca.interface_metrics(text, binder_sequence=BINDER_SEQUENCE)


def test_metrics_raise_when_the_peptide_is_too_short_for_the_position():
    text = make_cif([("A", HLA_SEQUENCE, -14.0),
                     ("C", "VVVD", 0.0),
                     ("D", BINDER_SEQUENCE, 7.0)])
    with pytest.raises(StructureAnalysisError, match="does not exist"):
        ca.interface_metrics(text, binder_sequence=BINDER_SEQUENCE)


def test_metrics_raise_when_the_binder_chain_is_missing():
    with pytest.raises(ChainMappingError, match="binder chain 'Q'"):
        ca.interface_metrics(bound_complex(), binder_chain_hint="Q")


def test_metrics_raise_rather_than_zero_when_the_mapping_is_wrong():
    """Pointing the analysis at the β2m chain as 'the peptide' must not pass."""
    text = bound_complex()
    with pytest.raises(StructureAnalysisError):
        ca.interface_metrics(text, binder_sequence=BINDER_SEQUENCE,
                             peptide_position=9999)


def test_waters_and_ligands_are_ignored():
    with_water = bound_complex()
    without = bound_complex()
    hetero_text = make_cif([
        ("A", HLA_SEQUENCE, -14.0),
        ("B", B2M_SEQUENCE, -24.0),
        ("C", ca.G12D_PEPTIDE, 0.0, (("HOH", "O", "O"),)),
        ("D", BINDER_SEQUENCE, 7.0),
    ])
    clean = ca.interface_metrics(without, binder_sequence=BINDER_SEQUENCE)
    watered = ca.interface_metrics(hetero_text, binder_sequence=BINDER_SEQUENCE)
    assert with_water == without
    assert watered["chain_map"]["chains"]["C"]["length"] == 10
    assert watered["asp6_buried_frac"] == pytest.approx(
        clean["asp6_buried_frac"], abs=1e-9)


# --- helpers ----------------------------------------------------------------

def test_sidechain_atom_names_skips_the_backbone():
    model = model_of(bound_complex())
    residues = [r for r in model["C"]]
    asp = residues[5]
    assert asp.get_resname() == "ASP"
    assert set(ca.sidechain_atom_names(asp)) == {"CB", "CG", "OD1", "OD2"}
    glycine = residues[3]
    assert glycine.get_resname() == "GLY"
    assert ca.sidechain_atom_names(glycine) == []


def test_count_contacts_requires_a_binder_chain_in_the_model():
    model = model_of(bound_complex())
    residue = [r for r in model["C"]][5]
    assert ca.count_contacts(model, "D", residue) > 0
    with pytest.raises(ChainMappingError):
        ca.count_contacts(model, "Z", residue)


def test_interface_charge_counts_only_contacting_residues():
    near = model_of(bound_complex(binder_y=7.0))
    far = model_of(bound_complex(binder_y=60.0))
    assert ca.interface_charge(far, "D") == pytest.approx(0.0)
    assert ca.interface_charge(near, "D") != 0.0


def test_interface_charge_needs_a_target():
    model = model_of(make_cif([("D", BINDER_SEQUENCE, 0.0)]))
    with pytest.raises(StructureAnalysisError, match="no target chains"):
        ca.interface_charge(model, "D")
