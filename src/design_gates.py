#!/usr/bin/env python3
"""Hard gates and ranking terms for Challenge 5 designs (see `docs/PIPELINE.md`).

This module is the single machine-readable copy of the thresholds documented in
`docs/PIPELINE.md` section 6, so the document and the tooling cannot drift
apart. It scores a *design record* - one row of the scored-designs table - and
reports which gates it fails.

The counter-screen panel is four-way: G12D (on-target), WT, apo HLA+B2M and an
unrelated same-allele decoy peptide (`AVFDRKSDAK`, see
`data/targets/decoy_a1101_complex.fasta`).

Usage:
    from design_gates import evaluate, specificity_gap
"""
from __future__ import annotations

from typing import Dict, List, Mapping, Optional

# --- ipTM gates -------------------------------------------------------------
ON_TARGET_IPTM_MIN = 0.75          # ipTM(G12D) >= 0.75
OFF_TARGET_IPTM_MAX = 0.50         # max(WT, apo, decoy) <= 0.50
SPECIFICITY_GAP_MIN = 0.25         # score(G12D) - max(off-targets) >= 0.25

# --- interface geometry gates -----------------------------------------------
PEPTIDE_BURIED_FRACTION_MIN = 0.40  # >= 40 % of buried area on the peptide chain
ASP6_BURIED_FRACTION_MIN = 0.60     # >= 60 % of the Asp6 side chain buried
ASP6_BURIED_FRACTION_PREFERRED = 0.75
MIN_POS6_HEAVY_ATOM_CONTACTS = 3    # binder heavy atoms within 5 A of position 6
PAE_INTERFACE_MAX = 10.0            # Angstrom, on the G12D run

# --- interface electrostatics (Xie et al. 2025: excess negative charge hurts) -
NET_INTERFACE_CHARGE_MIN = -2.0
NET_INTERFACE_CHARGE_MAX = 4.0

# Off-target states scored against the same binder, in panel order.
OFF_TARGET_STATES = ("wt", "apo", "decoy")

# Column names in the scored-designs table.
IPTM_COLUMNS = {
    "g12d": "iptm_g12d",
    "wt": "iptm_wt",
    "apo": "iptm_apo",
    "decoy": "iptm_decoy",
}
ASP6_BURIED_COLUMN = "asp6_buried_frac"
PEPTIDE_BURIED_COLUMN = "peptide_buried_frac"
NET_INTERFACE_CHARGE_COLUMN = "net_interface_charge"
POS6_CONTACTS_COLUMN = "pos6_contacts"
PAE_INTERFACE_COLUMN = "pae_interface_g12d"

DECOY_PEPTIDE = "AVFDRKSDAK"


class GateError(ValueError):
    """Raised when a design record cannot be interpreted at all."""


def _as_float(record: Mapping, key: str) -> Optional[float]:
    """Return record[key] as a float, or None when absent/blank/unparseable."""
    if key not in record:
        return None
    raw = record[key]
    if raw is None or (isinstance(raw, str) and not raw.strip()):
        return None
    try:
        return float(raw)
    except (TypeError, ValueError):
        raise GateError(f"column '{key}' is not numeric: {raw!r}")


def off_target_scores(record: Mapping) -> Dict[str, float]:
    """The off-target ipTMs present in `record`, keyed by state."""
    out: Dict[str, float] = {}
    for state in OFF_TARGET_STATES:
        value = _as_float(record, IPTM_COLUMNS[state])
        if value is not None:
            out[state] = value
    return out


def worst_off_target(record: Mapping) -> Optional[float]:
    """max(ipTM(WT), ipTM(apo), ipTM(decoy)) over whichever are present."""
    scores = off_target_scores(record)
    return max(scores.values()) if scores else None


def specificity_gap(record: Mapping) -> Optional[float]:
    """score(G12D) - max(score(WT), score(apo), score(decoy)).

    Returns None when the on-target score or every off-target score is missing.
    """
    on = _as_float(record, IPTM_COLUMNS["g12d"])
    off = worst_off_target(record)
    if on is None or off is None:
        return None
    return on - off


def asp6_buried_fraction(sasa_free: float, sasa_complex: float) -> float:
    """Buried fraction of the Asp6 side chain.

    `sasa_free` is the solvent-accessible surface area of the peptide position 6
    aspartate side chain in the *free* pMHC (binder absent); `sasa_complex` is
    the same quantity in the binder-bound complex. Both in square Angstrom.
    """
    if sasa_free < 0 or sasa_complex < 0:
        raise GateError("SASA values must be non-negative")
    if sasa_free == 0:
        raise GateError(
            "Asp6 side chain has zero SASA in the free pMHC: nothing to bury "
            "(check the chain/residue selection)"
        )
    if sasa_complex > sasa_free:
        raise GateError(
            f"Asp6 SASA in the complex ({sasa_complex}) exceeds the free pMHC "
            f"value ({sasa_free}): the binder cannot expose a buried residue"
        )
    return (sasa_free - sasa_complex) / sasa_free


def evaluate(record: Mapping) -> List[str]:
    """Return the gate failures for one design record; empty list means pass.

    Missing columns are not failures - a partially scored design simply has
    fewer gates applied - so this can be run mid-pipeline. Use `is_complete` to
    check that every gate could actually be evaluated.
    """
    failures: List[str] = []

    on = _as_float(record, IPTM_COLUMNS["g12d"])
    if on is not None and on < ON_TARGET_IPTM_MIN:
        failures.append(
            f"ipTM(G12D) {on:.3f} < {ON_TARGET_IPTM_MIN:.2f}"
        )

    for state, value in off_target_scores(record).items():
        if value > OFF_TARGET_IPTM_MAX:
            failures.append(
                f"ipTM({state}) {value:.3f} > {OFF_TARGET_IPTM_MAX:.2f} "
                f"(counter-screen gate: max(WT, apo, decoy) <= "
                f"{OFF_TARGET_IPTM_MAX:.2f})"
            )

    gap = specificity_gap(record)
    if gap is not None and gap < SPECIFICITY_GAP_MIN:
        failures.append(
            f"specificity_gap {gap:.3f} < {SPECIFICITY_GAP_MIN:.2f}"
        )

    pep = _as_float(record, PEPTIDE_BURIED_COLUMN)
    if pep is not None and pep < PEPTIDE_BURIED_FRACTION_MIN:
        failures.append(
            f"peptide-buried fraction {pep:.2f} < "
            f"{PEPTIDE_BURIED_FRACTION_MIN:.2f}"
        )

    asp6 = _as_float(record, ASP6_BURIED_COLUMN)
    if asp6 is not None and asp6 < ASP6_BURIED_FRACTION_MIN:
        failures.append(
            f"Asp6 buried fraction {asp6:.2f} < {ASP6_BURIED_FRACTION_MIN:.2f}"
        )

    contacts = _as_float(record, POS6_CONTACTS_COLUMN)
    if contacts is not None and contacts < MIN_POS6_HEAVY_ATOM_CONTACTS:
        failures.append(
            f"only {contacts:.0f} binder heavy atoms within 5 A of peptide "
            f"position 6 (need >= {MIN_POS6_HEAVY_ATOM_CONTACTS})"
        )

    pae = _as_float(record, PAE_INTERFACE_COLUMN)
    if pae is not None and pae >= PAE_INTERFACE_MAX:
        failures.append(
            f"PAE-interface(G12D) {pae:.1f} A >= {PAE_INTERFACE_MAX:.0f} A"
        )

    charge = _as_float(record, NET_INTERFACE_CHARGE_COLUMN)
    if charge is not None:
        if charge < NET_INTERFACE_CHARGE_MIN:
            failures.append(
                f"net interface charge {charge:+.1f} below "
                f"{NET_INTERFACE_CHARGE_MIN:+.1f}: excess negative charge at "
                f"the interface (Xie et al. 2025)"
            )
        elif charge > NET_INTERFACE_CHARGE_MAX:
            failures.append(
                f"net interface charge {charge:+.1f} above "
                f"{NET_INTERFACE_CHARGE_MAX:+.1f}: non-specific electrostatic "
                f"attraction to the groove"
            )

    return failures


REQUIRED_COLUMNS = (
    IPTM_COLUMNS["g12d"],
    IPTM_COLUMNS["wt"],
    IPTM_COLUMNS["apo"],
    IPTM_COLUMNS["decoy"],
    PEPTIDE_BURIED_COLUMN,
    ASP6_BURIED_COLUMN,
    NET_INTERFACE_CHARGE_COLUMN,
)


def missing_columns(record: Mapping) -> List[str]:
    """Gate inputs that `record` does not supply (blank counts as missing)."""
    out = []
    for col in REQUIRED_COLUMNS:
        raw = record.get(col)
        if raw is None or (isinstance(raw, str) and not raw.strip()):
            out.append(col)
    return out


def is_complete(record: Mapping) -> bool:
    """True when every gate in `REQUIRED_COLUMNS` can be evaluated."""
    return not missing_columns(record)


def passes(record: Mapping) -> bool:
    """True when the record is fully scored and clears every hard gate."""
    return is_complete(record) and not evaluate(record)
