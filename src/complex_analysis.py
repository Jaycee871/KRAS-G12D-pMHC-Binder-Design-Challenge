#!/usr/bin/env python3
"""Chain mapping and interface geometry for a predicted binder/pMHC complex.

Why this module exists
----------------------
Batch-01 scored every design with ``peptide_buried_frac = 0.00``,
``asp6_buried_frac = 0.00`` and ``0`` position-6 contacts, which is exactly
what a *broken chain mapping* looks like - and also exactly what a *real*
docking failure looks like. The old code could not tell the two apart: it
located the peptide by the hard-coded chain id ``"C"``, located position 6 by
list index, never checked the residue it found was the aspartate, and returned
``0.0`` from every path where a chain, a residue or an atom was missing.

So this module does two things the old code did not:

1. **Maps chains by content, not by convention.** Chains are identified from
   their one-letter sequences and lengths (peptide by sequence match, β2m by
   its invariant sequence, HLA heavy chain as the long remaining chain, binder
   by the designed sequence), and the mapping is returned alongside the
   metrics so every number can be traced to a chain.
2. **Fails loudly.** A missing chain, an ambiguous mapping, a peptide too short
   for the mutated position, the wrong residue at that position, or a
   side chain with zero reference SASA all raise `ChainMappingError` /
   `StructureAnalysisError`. The only zero this module ever returns is a
   measured zero - a binder that really does bury no target surface - and that
   case is flagged explicitly in ``notes``.

Everything is computed on the *same* predicted coordinates: Shrake-Rupley SASA
with a 1.4 Å probe, once with the binder present and once with the binder chain
detached, so the two differ only by the binder.
"""
from __future__ import annotations

import io
from dataclasses import dataclass, field, asdict
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

# --- constants --------------------------------------------------------------
SASA_PROBE_RADIUS = 1.4      # Å, PIPELINE.md section 6
CONTACT_CUTOFF = 5.0         # Å, heavy-atom contact definition
INTERFACE_CUTOFF = 5.0       # Å, "interface residue" definition

# --- Asp6 hydrogen bond / salt bridge (PIPELINE.md section 6, batch-02) -----
# Replaces the unreachable `asp6_buried_frac >= 0.60` gate: in 7OW6 the groove
# alone buries ~90 % of the Asp6 side chain, and the selective TCR of Poole et
# al. 2022 scores 0.00 on buried fraction while making 7 contacts with
# position 6. What a *selective* paratope has to do is read the carboxylate,
# i.e. donate a hydrogen bond to it.
HBOND_DISTANCE = 3.5         # Å, donor heavy atom to OD1/OD2
SALT_BRIDGE_DISTANCE = 4.0   # Å, charged donor group to OD1/OD2
ASP_CARBOXYLATE_ATOMS = ("OD1", "OD2")
# Side-chain hydroxyl / sulfhydryl donors, by residue.
HYDROXYL_DONORS = {
    "SER": {"OG"}, "THR": {"OG1"}, "TYR": {"OH"}, "CYS": {"SG"},
}
# Side-chain nitrogen donors that also carry (part of) a positive charge.
CHARGED_DONORS = {
    "ARG": {"NE", "NH1", "NH2"},
    "LYS": {"NZ"},
    "HIS": {"ND1", "NE2"},
}

BACKBONE_ATOMS = frozenset({"N", "CA", "C", "O", "OXT"})

# Peptide chemistry of the four target states (PIPELINE.md section 0).
G12D_PEPTIDE = "VVVGADGVGK"
WT_PEPTIDE = "VVVGAGGVGK"
DECOY_PEPTIDE = "AVFDRKSDAK"
MUTATED_PEPTIDE_POSITION = 6     # 1-based; the G12D aspartate
MUTATED_RESIDUE_NAME = "ASP"

# β2-microglobulin is invariant, so it can be recognised outright.
B2M_PREFIX = "IQRTPKIQVYSRHPAENGKSNFLNCYVSGFHPSDIEVDLLKNGERIEK"
B2M_LENGTH_RANGE = (95, 102)
HLA_MIN_LENGTH = 150             # α1/α2(/α3) heavy chain
PEPTIDE_MAX_LENGTH = 15

# Formal charge per residue at pH 7.4; His counted at ~10 % protonation.
RESIDUE_CHARGE = {"ARG": 1.0, "LYS": 1.0, "ASP": -1.0, "GLU": -1.0, "HIS": 0.1}

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V", "SEC": "U", "PYL": "O", "MSE": "M",
}


class StructureAnalysisError(ValueError):
    """The structure cannot be analysed as a binder/pMHC complex."""


class ChainMappingError(StructureAnalysisError):
    """A chain could not be identified, or the mapping is ambiguous."""


# --- chain inventory --------------------------------------------------------

@dataclass
class ChainInfo:
    """One polymer chain of a parsed model."""
    chain_id: str
    length: int
    sequence: str
    first_resseq: int
    last_resseq: int

    def as_dict(self) -> Dict[str, Any]:
        return asdict(self)


@dataclass
class ChainMap:
    """Which chain of a complex is which molecule."""
    binder: str
    peptide: Optional[str] = None
    hla: Optional[str] = None
    b2m: Optional[str] = None
    chains: Dict[str, ChainInfo] = field(default_factory=dict)

    @property
    def target_chains(self) -> List[str]:
        return [cid for cid in self.chains if cid != self.binder]

    def as_dict(self) -> Dict[str, Any]:
        return {
            "binder": self.binder,
            "peptide": self.peptide,
            "hla": self.hla,
            "b2m": self.b2m,
            "chains": {cid: info.as_dict() for cid, info in self.chains.items()},
        }

    def describe(self) -> str:
        return ", ".join(
            f"{cid}:{info.length}aa" for cid, info in self.chains.items()
        )


def one_letter(residue) -> str:
    return THREE_TO_ONE.get(residue.get_resname().strip().upper(), "X")


def is_polymer_residue(residue) -> bool:
    """True for a standard amino-acid residue (not water, ion, ligand)."""
    hetflag = residue.id[0]
    name = residue.get_resname().strip().upper()
    if hetflag == "W":
        return False
    if hetflag != " " and name != "MSE":
        return False
    return name in THREE_TO_ONE


def chain_inventory(model) -> Dict[str, ChainInfo]:
    """Per-chain length, sequence and residue numbering of a parsed model."""
    inventory: Dict[str, ChainInfo] = {}
    for chain in model:
        residues = [r for r in chain if is_polymer_residue(r)]
        if not residues:
            continue
        inventory[chain.id] = ChainInfo(
            chain_id=chain.id,
            length=len(residues),
            sequence="".join(one_letter(r) for r in residues),
            first_resseq=residues[0].id[1],
            last_resseq=residues[-1].id[1],
        )
    if not inventory:
        raise StructureAnalysisError("structure contains no polymer chains")
    return inventory


def map_complex_chains(
    model,
    binder_sequence: Optional[str] = None,
    binder_chain_hint: Optional[str] = None,
    peptide_sequences: Sequence[str] = (G12D_PEPTIDE, WT_PEPTIDE, DECOY_PEPTIDE),
    require_peptide: bool = True,
) -> ChainMap:
    """Work out which chain is the binder, the peptide, the HLA and β2m.

    Identification is by *content*:

    * **binder** - the chain whose sequence equals `binder_sequence` when one is
      given, else the chain named by `binder_chain_hint`, else (when exactly one
      chain is left over after the target chains are identified) that chain.
    * **peptide** - a chain whose sequence is one of `peptide_sequences`, else
      the single short (≤ 15 aa) non-binder chain.
    * **β2m** - a ~99-residue chain starting with the invariant β2m sequence.
    * **HLA heavy chain** - the remaining long (≥ 150 aa) chain.

    Raises `ChainMappingError` - never guesses - when a required chain is
    absent or two chains match equally well.
    """
    inventory = chain_inventory(model)
    description = ", ".join(
        f"{cid} ({info.length} aa, {info.sequence[:12]}…)"
        for cid, info in inventory.items()
    )

    binder: Optional[str] = None
    if binder_sequence:
        wanted = binder_sequence.strip().upper()
        matches = [cid for cid, info in inventory.items() if info.sequence == wanted]
        if len(matches) > 1:
            raise ChainMappingError(
                f"binder sequence of {len(wanted)} aa matches {len(matches)} "
                f"chains {matches}; cannot tell which is the binder. "
                f"Chains present: {description}"
            )
        if matches:
            binder = matches[0]
        elif binder_chain_hint is None:
            raise ChainMappingError(
                f"no chain in the structure has the designed binder sequence "
                f"({len(wanted)} aa). Chains present: {description}"
            )
    if binder is None and binder_chain_hint is not None:
        if binder_chain_hint not in inventory:
            raise ChainMappingError(
                f"binder chain {binder_chain_hint!r} is not in the structure. "
                f"Chains present: {description}"
            )
        binder = binder_chain_hint

    remaining = {cid: info for cid, info in inventory.items() if cid != binder}

    # β2-microglobulin: invariant sequence, ~99 residues.
    b2m = None
    b2m_hits = [
        cid for cid, info in remaining.items()
        if info.sequence.startswith(B2M_PREFIX[:20])
        or (B2M_LENGTH_RANGE[0] <= info.length <= B2M_LENGTH_RANGE[1]
            and info.sequence[:4] == B2M_PREFIX[:4])
    ]
    if len(b2m_hits) == 1:
        b2m = b2m_hits[0]
    elif len(b2m_hits) > 1:
        raise ChainMappingError(
            f"{len(b2m_hits)} chains look like β2-microglobulin ({b2m_hits}); "
            f"chains present: {description}"
        )

    # Peptide: exact sequence match first, then "the one short chain".
    wanted_peptides = {p.strip().upper() for p in peptide_sequences if p}
    peptide_hits = [
        cid for cid, info in remaining.items()
        if info.sequence in wanted_peptides and cid != b2m
    ]
    if not peptide_hits:
        peptide_hits = [
            cid for cid, info in remaining.items()
            if cid != b2m and info.length <= PEPTIDE_MAX_LENGTH
        ]
    peptide: Optional[str] = None
    if len(peptide_hits) == 1:
        peptide = peptide_hits[0]
    elif len(peptide_hits) > 1:
        raise ChainMappingError(
            f"{len(peptide_hits)} chains could be the peptide ({peptide_hits}); "
            f"chains present: {description}"
        )
    elif require_peptide:
        raise ChainMappingError(
            f"no peptide chain in the structure: expected one of "
            f"{sorted(wanted_peptides)} or a chain of at most "
            f"{PEPTIDE_MAX_LENGTH} residues. Chains present: {description}"
        )

    # HLA heavy chain: the long chain that is left.
    hla_hits = [
        cid for cid, info in remaining.items()
        if cid not in (b2m, peptide) and info.length >= HLA_MIN_LENGTH
    ]
    hla = hla_hits[0] if len(hla_hits) == 1 else None
    if len(hla_hits) > 1:
        raise ChainMappingError(
            f"{len(hla_hits)} chains could be the HLA heavy chain ({hla_hits}); "
            f"chains present: {description}"
        )

    if binder is None:
        leftovers = [cid for cid in remaining if cid not in (b2m, peptide, hla)]
        if len(leftovers) != 1:
            raise ChainMappingError(
                f"cannot identify the binder chain: {len(leftovers)} unassigned "
                f"chains {leftovers}. Chains present: {description}"
            )
        binder = leftovers[0]

    return ChainMap(binder=binder, peptide=peptide, hla=hla, b2m=b2m,
                    chains=inventory)


# --- parsing ----------------------------------------------------------------

def parse_structure(text: str, fmt: str = "mmcif"):
    """Parse mmCIF (or PDB) text into a Bio.PDB structure."""
    if fmt == "mmcif":
        from Bio.PDB import MMCIFParser

        parser = MMCIFParser(QUIET=True)
    else:
        from Bio.PDB import PDBParser

        parser = PDBParser(QUIET=True)
    structure = parser.get_structure("complex", io.StringIO(text))
    try:
        next(iter(structure))
    except StopIteration:
        raise StructureAnalysisError("structure contains no models")
    return structure


def _clean_model(model):
    """Drop waters, ligands, hydrogens and alternate locations, in place.

    Experimental structures (the 7OW6 validation case) carry waters, glycans
    and altlocs that an SASA computed on a *predicted* model never has; leaving
    them in silently changes every buried-area number.
    """
    for chain in list(model):
        for residue in list(chain):
            if not is_polymer_residue(residue):
                chain.detach_child(residue.id)
                continue
            for atom in list(residue):
                if atom.element == "H":
                    residue.detach_child(atom.get_id())
                elif atom.is_disordered():
                    try:
                        residue.detach_child(atom.get_id())
                        residue.add(atom.disordered_get())
                    except Exception:  # noqa: BLE001
                        pass
        if not list(chain):
            model.detach_child(chain.id)
    return model


def _atom_sasa(model) -> Dict[Tuple[str, int, str], float]:
    from Bio.PDB.SASA import ShrakeRupley

    ShrakeRupley(probe_radius=SASA_PROBE_RADIUS).compute(model, level="A")
    return {
        (atom.get_parent().get_parent().id, atom.get_parent().id[1], atom.get_id()):
        atom.sasa
        for atom in model.get_atoms()
    }


def sidechain_atom_names(residue) -> List[str]:
    """Heavy side-chain atom names of `residue` (empty for glycine)."""
    return [
        atom.get_id() for atom in residue
        if atom.get_id() not in BACKBONE_ATOMS and atom.element != "H"
    ]


# --- metrics ----------------------------------------------------------------

def interface_metrics(
    text: str,
    binder_sequence: Optional[str] = None,
    binder_chain_hint: Optional[str] = None,
    peptide_position: int = MUTATED_PEPTIDE_POSITION,
    expected_residue: Optional[str] = MUTATED_RESIDUE_NAME,
    fmt: str = "mmcif",
    peptide_sequences: Sequence[str] = (G12D_PEPTIDE, WT_PEPTIDE, DECOY_PEPTIDE),
) -> Dict[str, Any]:
    """Interface geometry of a binder/pMHC complex, with the chain mapping.

    Returns `peptide_buried_frac`, `asp6_buried_frac`, `pos6_contacts` and
    `net_interface_charge` (PIPELINE.md section 6) plus the per-chain buried
    areas and contact counts that say *where* the binder actually sits, and the
    `chain_map` the numbers were computed from.

    Raises `ChainMappingError` / `StructureAnalysisError` rather than returning
    zero when something is missing. A genuine zero - a binder that buries no
    target surface at all - comes back as 0.0 with a note in ``notes``.
    """
    structure = parse_structure(text, fmt=fmt)
    model = _clean_model(next(iter(structure)))
    chain_map = map_complex_chains(
        model,
        binder_sequence=binder_sequence,
        binder_chain_hint=binder_chain_hint,
        peptide_sequences=peptide_sequences,
    )
    binder_chain, peptide_chain = chain_map.binder, chain_map.peptide
    if peptide_chain is None:
        raise ChainMappingError(
            f"no peptide chain identified in {chain_map.describe()}"
        )

    notes: List[str] = []
    bound = _atom_sasa(model)

    free_structure = parse_structure(text, fmt=fmt)
    free_model = _clean_model(next(iter(free_structure)))
    if binder_chain not in [c.id for c in free_model]:
        raise ChainMappingError(
            f"binder chain {binder_chain!r} vanished from the reference model"
        )
    free_model.detach_child(binder_chain)
    free = _atom_sasa(free_model)

    buried_total = 0.0
    buried_by_chain: Dict[str, float] = {}
    for key, free_value in free.items():
        delta = free_value - bound.get(key, 0.0)
        if delta <= 0:
            continue
        buried_total += delta
        buried_by_chain[key[0]] = buried_by_chain.get(key[0], 0.0) + delta

    buried_peptide = buried_by_chain.get(peptide_chain, 0.0)
    if buried_total <= 0.0:
        notes.append(
            f"binder chain {binder_chain} buries no target surface at all: it "
            f"makes no contact with the pMHC in this prediction"
        )
        peptide_buried_frac = 0.0
    else:
        peptide_buried_frac = buried_peptide / buried_total
        if buried_peptide <= 0.0:
            notes.append(
                f"binder buries {buried_total:.0f} Å² of the target but none of "
                f"it on the peptide chain {peptide_chain}: the paratope is off "
                f"the peptide"
            )

    metrics: Dict[str, Any] = {
        "chain_map": chain_map.as_dict(),
        "binder_chain": binder_chain,
        "peptide_chain": peptide_chain,
        "peptide_buried_frac": peptide_buried_frac,
        "buried_area_total": buried_total,
        "buried_area_peptide": buried_peptide,
        "buried_area_by_chain": {k: round(v, 1) for k, v in buried_by_chain.items()},
    }

    # --- the mutated position --------------------------------------------
    residues = [r for r in model[peptide_chain] if is_polymer_residue(r)]
    if len(residues) < peptide_position:
        raise StructureAnalysisError(
            f"peptide chain {peptide_chain} has {len(residues)} residues, so "
            f"position {peptide_position} does not exist "
            f"(sequence {chain_map.chains[peptide_chain].sequence})"
        )
    target_residue = residues[peptide_position - 1]
    resname = target_residue.get_resname().strip().upper()
    if expected_residue and resname != expected_residue.upper():
        raise StructureAnalysisError(
            f"peptide position {peptide_position} of chain {peptide_chain} is "
            f"{resname}, expected {expected_residue.upper()} "
            f"(sequence {chain_map.chains[peptide_chain].sequence}); the "
            f"chain or the numbering is wrong"
        )
    metrics["pos6_residue"] = resname
    metrics["pos6_resseq"] = target_residue.id[1]

    names = sidechain_atom_names(target_residue)
    if not names:
        raise StructureAnalysisError(
            f"peptide position {peptide_position} ({resname}) has no heavy "
            f"side-chain atoms, so its buried fraction is undefined"
        )
    resseq = target_residue.id[1]
    missing = [n for n in names if (peptide_chain, resseq, n) not in free]
    if missing:
        raise StructureAnalysisError(
            f"side-chain atoms {missing} of {resname}{resseq} are missing from "
            f"the binder-free reference model"
        )
    free_sc = sum(free[(peptide_chain, resseq, n)] for n in names)
    bound_sc = sum(bound.get((peptide_chain, resseq, n), 0.0) for n in names)
    if free_sc <= 0.0:
        raise StructureAnalysisError(
            f"{resname}{resseq} side chain has zero SASA in the binder-free "
            f"pMHC, so there is nothing for the binder to bury; the chain "
            f"selection or the structure is wrong "
            f"(chains: {chain_map.describe()})"
        )
    metrics["asp6_sasa_free"] = free_sc
    metrics["asp6_sasa_bound"] = bound_sc
    metrics["asp6_buried_frac"] = (free_sc - min(bound_sc, free_sc)) / free_sc

    partners = asp6_hbond_partners(model, binder_chain, target_residue)
    metrics["asp6_hbond_partners"] = partners
    metrics["asp6_hbond"] = any(p["hbond"] for p in partners)
    metrics["asp6_salt_bridge"] = any(p["salt_bridge"] for p in partners)
    metrics["asp6_min_donor_distance"] = (
        partners[0]["distance"] if partners else None)
    if not partners:
        notes.append(
            f"no binder hydrogen-bond donor within {HBOND_DISTANCE} Å of the "
            f"{resname}{resseq} carboxylate: the binder does not read the "
            f"mutated side chain"
        )

    contacts = count_contacts(model, binder_chain, target_residue)
    metrics["pos6_contacts"] = float(contacts)
    metrics["contacts_by_chain"] = contacts_by_chain(model, binder_chain)
    metrics["net_interface_charge"] = interface_charge(model, binder_chain)
    metrics["notes"] = notes
    return metrics


def is_donor_atom(residue, atom) -> bool:
    """True when `atom` of `residue` can donate a hydrogen bond.

    Any nitrogen (backbone amide or side chain) counts, as do the side-chain
    hydroxyls of Ser/Thr/Tyr and the Cys sulfhydryl. Hydrogens are not present
    in predicted structures, so donor identity is decided by chemistry rather
    than by an explicit H.
    """
    name = atom.get_id().strip().upper()
    resname = residue.get_resname().strip().upper()
    if atom.element == "N":
        return True
    return name in HYDROXYL_DONORS.get(resname, ())


def asp6_hbond_partners(
    model,
    binder_chain: str,
    residue,
    hbond_distance: float = HBOND_DISTANCE,
    salt_bridge_distance: float = SALT_BRIDGE_DISTANCE,
) -> List[Dict[str, Any]]:
    """Binder donors hydrogen bonding to the Asp6 carboxylate.

    One record per donor atom within `hbond_distance` of OD1/OD2, flagged as a
    salt bridge when the donor is a charged Arg/Lys/His group within
    `salt_bridge_distance`. Sorted by distance, closest first.
    """
    from Bio.PDB import NeighborSearch

    if binder_chain not in [c.id for c in model]:
        raise ChainMappingError(f"binder chain {binder_chain!r} is not in the model")
    carboxylate = [a for a in residue
                   if a.get_id().strip().upper() in ASP_CARBOXYLATE_ATOMS]
    if not carboxylate:
        raise StructureAnalysisError(
            f"{residue.get_resname()}{residue.id[1]} has no "
            f"{'/'.join(ASP_CARBOXYLATE_ATOMS)} carboxylate oxygens, so no "
            f"hydrogen bond to it can be measured"
        )
    binder_atoms = heavy_atoms(model[binder_chain])
    if not binder_atoms:
        raise StructureAnalysisError(f"binder chain {binder_chain} has no heavy atoms")
    search = NeighborSearch(binder_atoms)

    best: Dict[Tuple[str, int, str], Dict[str, Any]] = {}
    cutoff = max(hbond_distance, salt_bridge_distance)
    for oxygen in carboxylate:
        for atom in search.search(oxygen.coord, cutoff):
            donor_residue = atom.get_parent()
            if not is_donor_atom(donor_residue, atom):
                continue
            resname = donor_residue.get_resname().strip().upper()
            name = atom.get_id().strip().upper()
            distance = float(((atom.coord - oxygen.coord) ** 2).sum() ** 0.5)
            charged = name in CHARGED_DONORS.get(resname, ())
            is_hbond = distance <= hbond_distance
            is_salt_bridge = charged and distance <= salt_bridge_distance
            if not (is_hbond or is_salt_bridge):
                continue
            key = (resname, donor_residue.id[1], name)
            record = {
                "donor_residue": resname,
                "donor_resseq": donor_residue.id[1],
                "donor_atom": name,
                "acceptor_atom": oxygen.get_id().strip().upper(),
                "distance": round(distance, 2),
                "hbond": is_hbond,
                "salt_bridge": is_salt_bridge,
            }
            if key not in best or distance < best[key]["distance"]:
                best[key] = record
    return sorted(best.values(), key=lambda r: r["distance"])


def heavy_atoms(entity) -> List:
    return [a for a in entity.get_atoms() if a.element != "H"]


def count_contacts(model, binder_chain: str, residue) -> int:
    """Binder heavy atoms within `CONTACT_CUTOFF` of `residue`'s heavy atoms."""
    from Bio.PDB import NeighborSearch

    if binder_chain not in [c.id for c in model]:
        raise ChainMappingError(
            f"binder chain {binder_chain!r} is not in the model"
        )
    binder_atoms = heavy_atoms(model[binder_chain])
    if not binder_atoms:
        raise StructureAnalysisError(
            f"binder chain {binder_chain} has no heavy atoms"
        )
    search = NeighborSearch(binder_atoms)
    hits = set()
    for atom in residue:
        if atom.element == "H":
            continue
        for neighbour in search.search(atom.coord, CONTACT_CUTOFF):
            hits.add(id(neighbour))
    return len(hits)


def contacts_by_chain(model, binder_chain: str) -> Dict[str, int]:
    """Binder heavy atoms within the cutoff of each other chain."""
    from Bio.PDB import NeighborSearch

    binder_atoms = heavy_atoms(model[binder_chain])
    search = NeighborSearch(binder_atoms)
    out: Dict[str, int] = {}
    for chain in model:
        if chain.id == binder_chain:
            continue
        hits = set()
        for atom in heavy_atoms(chain):
            for neighbour in search.search(atom.coord, CONTACT_CUTOFF):
                hits.add(id(neighbour))
        out[chain.id] = len(hits)
    return out


def interface_charge(model, binder_chain: str) -> float:
    """Net formal charge of the binder residues forming the interface."""
    from Bio.PDB import NeighborSearch

    target_atoms = [
        atom for chain in model if chain.id != binder_chain
        for atom in heavy_atoms(chain)
    ]
    if not target_atoms:
        raise StructureAnalysisError(
            f"no target chains besides the binder {binder_chain!r}"
        )
    search = NeighborSearch(target_atoms)
    charge = 0.0
    for residue in model[binder_chain]:
        if any(search.search(atom.coord, INTERFACE_CUTOFF)
               for atom in residue if atom.element != "H"):
            charge += RESIDUE_CHARGE.get(residue.get_resname().strip().upper(), 0.0)
    return charge
