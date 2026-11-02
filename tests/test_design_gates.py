import pytest

import design_gates as dg


def good_record(**overrides):
    """A design that clears every gate; override one field per test."""
    record = {
        "name": "d1",
        "iptm_g12d": 0.86,
        "iptm_wt": 0.31,
        "iptm_apo": 0.18,
        "iptm_decoy": 0.27,
        "peptide_buried_frac": 0.58,
        "asp6_buried_frac": 0.81,
        "pos6_contacts": 7,
        "pae_interface_g12d": 6.2,
        "net_interface_charge": 1.0,
    }
    record.update(overrides)
    return record


# --- thresholds match docs/PIPELINE.md --------------------------------------

def test_documented_thresholds():
    assert dg.ON_TARGET_IPTM_MIN == 0.75
    assert dg.OFF_TARGET_IPTM_MAX == 0.50
    assert dg.SPECIFICITY_GAP_MIN == 0.25
    assert dg.PEPTIDE_BURIED_FRACTION_MIN == 0.40
    assert dg.MIN_POS6_HEAVY_ATOM_CONTACTS == 3


def test_decoy_is_in_the_off_target_panel():
    assert dg.OFF_TARGET_STATES == ("wt", "apo", "decoy")
    assert dg.DECOY_PEPTIDE == "AVFDRKSDAK"


# --- specificity gap --------------------------------------------------------

def test_specificity_gap_uses_worst_off_target():
    rec = good_record(iptm_g12d=0.90, iptm_wt=0.20, iptm_apo=0.10,
                      iptm_decoy=0.45)
    assert dg.worst_off_target(rec) == 0.45
    assert dg.specificity_gap(rec) == pytest.approx(0.45)


def test_specificity_gap_none_when_on_target_missing():
    rec = good_record()
    del rec["iptm_g12d"]
    assert dg.specificity_gap(rec) is None


def test_specificity_gap_none_without_any_off_target():
    rec = {"iptm_g12d": 0.9}
    assert dg.worst_off_target(rec) is None
    assert dg.specificity_gap(rec) is None


def test_specificity_gap_tolerates_partial_panel():
    rec = {"iptm_g12d": 0.9, "iptm_wt": 0.4}
    assert dg.specificity_gap(rec) == pytest.approx(0.5)


def test_string_values_are_parsed():
    rec = {k: str(v) for k, v in good_record().items()}
    assert dg.evaluate(rec) == []


def test_non_numeric_value_raises():
    with pytest.raises(dg.GateError, match="not numeric"):
        dg.evaluate(good_record(iptm_decoy="n/a"))


# --- gates ------------------------------------------------------------------

def test_good_record_passes():
    assert dg.evaluate(good_record()) == []
    assert dg.passes(good_record())


def test_low_on_target_iptm_fails():
    failures = dg.evaluate(good_record(iptm_g12d=0.70, iptm_decoy=0.20))
    assert any("ipTM(G12D)" in f for f in failures)


def test_high_decoy_iptm_fails_the_counter_screen():
    failures = dg.evaluate(good_record(iptm_decoy=0.62))
    assert any("ipTM(decoy)" in f for f in failures)
    assert any("max(WT, apo, decoy)" in f for f in failures)


def test_high_wt_iptm_still_fails():
    assert any("ipTM(wt)" in f for f in dg.evaluate(good_record(iptm_wt=0.55)))


def test_high_apo_iptm_still_fails():
    assert any("ipTM(apo)" in f for f in dg.evaluate(good_record(iptm_apo=0.70)))


def test_decoy_binder_fails_even_with_clean_wt_and_apo():
    """A groove-rim binder: great on G12D, clean on WT/apo, binds the decoy."""
    failures = dg.evaluate(good_record(iptm_g12d=0.88, iptm_wt=0.10,
                                       iptm_apo=0.12, iptm_decoy=0.84))
    assert failures
    assert any("specificity_gap" in f for f in failures)


def test_narrow_gap_fails():
    failures = dg.evaluate(good_record(iptm_g12d=0.76, iptm_decoy=0.60))
    assert any("specificity_gap" in f for f in failures)


def test_low_peptide_buried_fraction_fails():
    failures = dg.evaluate(good_record(peptide_buried_frac=0.22))
    assert any("peptide-buried fraction" in f for f in failures)


def test_low_asp6_buried_fraction_fails():
    failures = dg.evaluate(good_record(asp6_buried_frac=0.35))
    assert any("Asp6 buried fraction" in f for f in failures)


def test_asp6_threshold_boundary_passes():
    assert dg.evaluate(
        good_record(asp6_buried_frac=dg.ASP6_BURIED_FRACTION_MIN)
    ) == []


def test_too_few_pos6_contacts_fails():
    failures = dg.evaluate(good_record(pos6_contacts=2))
    assert any("position 6" in f for f in failures)


def test_high_pae_interface_fails():
    failures = dg.evaluate(good_record(pae_interface_g12d=12.4))
    assert any("PAE-interface" in f for f in failures)


# --- interface charge -------------------------------------------------------

def test_excess_negative_interface_charge_fails():
    failures = dg.evaluate(good_record(net_interface_charge=-4))
    assert any("excess negative charge" in f for f in failures)


def test_excess_positive_interface_charge_fails():
    failures = dg.evaluate(good_record(net_interface_charge=6))
    assert any("non-specific electrostatic" in f for f in failures)


def test_charge_window_boundaries_pass():
    assert dg.evaluate(
        good_record(net_interface_charge=dg.NET_INTERFACE_CHARGE_MIN)) == []
    assert dg.evaluate(
        good_record(net_interface_charge=dg.NET_INTERFACE_CHARGE_MAX)) == []


# --- Asp6 buried fraction from SASA ----------------------------------------

def test_asp6_buried_fraction_from_sasa():
    assert dg.asp6_buried_fraction(100.0, 25.0) == pytest.approx(0.75)
    assert dg.asp6_buried_fraction(80.0, 0.0) == pytest.approx(1.0)
    assert dg.asp6_buried_fraction(80.0, 80.0) == pytest.approx(0.0)


def test_asp6_buried_fraction_rejects_zero_free_sasa():
    with pytest.raises(dg.GateError, match="zero SASA"):
        dg.asp6_buried_fraction(0.0, 0.0)


def test_asp6_buried_fraction_rejects_growth_on_binding():
    with pytest.raises(dg.GateError, match="exceeds the free pMHC"):
        dg.asp6_buried_fraction(50.0, 60.0)


def test_asp6_buried_fraction_rejects_negative_sasa():
    with pytest.raises(dg.GateError, match="non-negative"):
        dg.asp6_buried_fraction(-1.0, 0.0)


# --- completeness -----------------------------------------------------------

def test_missing_columns_reported():
    rec = good_record()
    del rec["iptm_decoy"]
    del rec["asp6_buried_frac"]
    assert set(dg.missing_columns(rec)) == {"iptm_decoy", "asp6_buried_frac"}
    assert not dg.is_complete(rec)
    assert not dg.passes(rec)
    # ...but the gates that *can* be checked still pass.
    assert dg.evaluate(rec) == []


def test_blank_counts_as_missing():
    assert dg.missing_columns(good_record(iptm_decoy="")) == ["iptm_decoy"]


def test_complete_record_is_complete():
    assert dg.is_complete(good_record())
