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

# --- ipTM: ranking terms, not hard gates (batch-02) -------------------------
# Batch-01 produced a best ipTM(G12D) of 0.870 and a best specificity gap of
# +0.254 while burying a mean of 21 A^2 on the peptide: the reported ipTM is
# the strongest binder-to-*any*-chain interface, so a binder parked on the HLA
# alpha-3 face still scores 0.4-0.5. Until the pose filter is proven these
# thresholds are reported as *shortfalls* (`ranking_shortfalls`) and used for
# ranking; they no longer veto a design in `evaluate`.
ON_TARGET_IPTM_MIN = 0.75          # ipTM(G12D) >= 0.75 (ranking target)
OFF_TARGET_IPTM_MAX = 0.50         # max(WT, apo, decoy) <= 0.50 (ranking target)
SPECIFICITY_GAP_MIN = 0.25         # score(G12D) - max(off-targets) (ranking target)

# --- interface geometry gates -----------------------------------------------
PEPTIDE_BURIED_FRACTION_MIN = 0.40  # >= 40 % of buried area on the peptide chain
# `asp6_buried_frac` is now a *ranking* term only. Measured on 7OW6 itself the
# groove alone buries 0.90 of the Asp6 side chain, leaving 8.8 A^2 exposed, and
# the selective TCR of Poole et al. 2022 buries none of that remainder - so a
# 0.60 hard gate was unreachable even by a known selective binder.
ASP6_BURIED_FRACTION_PREFERRED = 0.75
# What replaces it: the binder has to *read the carboxylate*. Either a donor
# nitrogen or hydroxyl oxygen within 3.5 A of OD1/OD2 (an Asp6 hydrogen bond /
# salt bridge), or - the Poole escape clause - a position-6 heavy-atom contact
# footprint at least as dense as that selective TCR's, which reads the
# aspartate indirectly (electrostatically) and makes no hydrogen bond to it
# because HLA Arg114 already caps the carboxylate.
MIN_POS6_HEAVY_ATOM_CONTACTS = 3    # binder heavy atoms within 5 A of position 6
POS6_CONTACTS_WITHOUT_HBOND = 7     # contacts that substitute for the H-bond
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
ASP6_HBOND_COLUMN = "asp6_hbond"
ASP6_SALT_BRIDGE_COLUMN = "asp6_salt_bridge"
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


_TRUE = {"true", "1", "yes", "y", "t"}
_FALSE = {"false", "0", "no", "n", "f", ""}


def _as_bool(record: Mapping, key: str) -> Optional[bool]:
    """Return record[key] as a bool, or None when absent/blank."""
    if key not in record:
        return None
    raw = record[key]
    if raw is None:
        return None
    if isinstance(raw, bool):
        return raw
    text = str(raw).strip().lower()
    if text in _TRUE:
        return True
    if text in _FALSE:
        return None if text == "" else False
    raise GateError(f"column '{key}' is not a boolean: {raw!r}")


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


def asp6_readout_failures(record: Mapping) -> List[str]:
    """The Asp6 read-out criterion that replaced `asp6_buried_frac >= 0.60`.

    Passes when the binder hydrogen bonds or salt bridges to the Asp6
    carboxylate (a donor N or hydroxyl O within 3.5 A of OD1/OD2), or makes at
    least `POS6_CONTACTS_WITHOUT_HBOND` heavy-atom contacts with position 6 -
    the footprint of the selective TCR of Poole et al. 2022, which discriminates
    G12D from WT without an Asp6 hydrogen bond.
    """
    hbond = _as_bool(record, ASP6_HBOND_COLUMN)
    salt = _as_bool(record, ASP6_SALT_BRIDGE_COLUMN)
    contacts = _as_float(record, POS6_CONTACTS_COLUMN)
    if hbond is None and salt is None:
        return []          # not measured: no opinion (mid-pipeline records)
    if hbond or salt:
        return []
    if contacts is not None and contacts >= POS6_CONTACTS_WITHOUT_HBOND:
        return []
    found = f"{contacts:.0f}" if contacts is not None else "an unknown number of"
    return [
        f"no hydrogen bond or salt bridge to the Asp6 carboxylate (donor "
        f"within 3.5 A of OD1/OD2) and only {found} heavy-atom contacts with "
        f"position 6 (need a bond, or >= {POS6_CONTACTS_WITHOUT_HBOND} contacts)"
    ]


def ranking_shortfalls(record: Mapping) -> List[str]:
    """ipTM targets the record misses. Ranking information, not hard gates.

    Demoted from `evaluate` in batch-02: the reported ipTM is the strongest
    binder-to-any-target-chain interface, so it says little about *where* the
    binder sits, which is the failure mode batch-01 actually had.
    """
    shortfalls: List[str] = []
    on = _as_float(record, IPTM_COLUMNS["g12d"])
    if on is not None and on < ON_TARGET_IPTM_MIN:
        shortfalls.append(f"ipTM(G12D) {on:.3f} < {ON_TARGET_IPTM_MIN:.2f}")
    for state, value in off_target_scores(record).items():
        if value > OFF_TARGET_IPTM_MAX:
            shortfalls.append(
                f"ipTM({state}) {value:.3f} > {OFF_TARGET_IPTM_MAX:.2f} "
                f"(counter-screen target: max(WT, apo, decoy) <= "
                f"{OFF_TARGET_IPTM_MAX:.2f})"
            )
    gap = specificity_gap(record)
    if gap is not None and gap < SPECIFICITY_GAP_MIN:
        shortfalls.append(
            f"specificity_gap {gap:.3f} < {SPECIFICITY_GAP_MIN:.2f}")
    return shortfalls


def evaluate(record: Mapping, include_ranking: bool = False) -> List[str]:
    """Return the hard-gate failures for one design record; empty means pass.

    Hard gates are the *geometry and chemistry* of the interface: peptide share
    of the buried area, position-6 contacts, the Asp6 read-out criterion,
    interface charge and interface PAE/PDE. The ipTM thresholds are ranking
    terms since batch-02 and are only included when `include_ranking` is true
    (see `ranking_shortfalls`).

    Missing columns are not failures - a partially scored design simply has
    fewer gates applied - so this can be run mid-pipeline. Use `is_complete` to
    check that every gate could actually be evaluated.
    """
    failures: List[str] = []

    pep = _as_float(record, PEPTIDE_BURIED_COLUMN)
    if pep is not None and pep < PEPTIDE_BURIED_FRACTION_MIN:
        failures.append(
            f"peptide-buried fraction {pep:.2f} < "
            f"{PEPTIDE_BURIED_FRACTION_MIN:.2f}"
        )

    contacts = _as_float(record, POS6_CONTACTS_COLUMN)
    if contacts is not None and contacts < MIN_POS6_HEAVY_ATOM_CONTACTS:
        failures.append(
            f"only {contacts:.0f} binder heavy atoms within 5 A of peptide "
            f"position 6 (need >= {MIN_POS6_HEAVY_ATOM_CONTACTS})"
        )

    failures.extend(asp6_readout_failures(record))

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

    if include_ranking:
        failures.extend(ranking_shortfalls(record))
    return failures


REQUIRED_COLUMNS = (
    IPTM_COLUMNS["g12d"],
    IPTM_COLUMNS["wt"],
    IPTM_COLUMNS["apo"],
    IPTM_COLUMNS["decoy"],
    PEPTIDE_BURIED_COLUMN,
    POS6_CONTACTS_COLUMN,
    ASP6_HBOND_COLUMN,
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
