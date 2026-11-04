#!/usr/bin/env python3
"""Paired, same-pose MM-GBSA ddG for the G12D vs WT peptide substitution.

Why this module exists
----------------------
Until 2026-11-03 designs were ranked on

    specificity_gap = ipTM(G12D) - max(ipTM(WT), ipTM(apo), ipTM(decoy))

with each ipTM coming from an *independent* Boltz-2 co-folding run. The 7OW6
positive control killed that quantity (`results/validation/CONTROL_7OW6.md`):
on the experimentally G12D-selective TCR of Poole et al. 2022 the gap is
-0.006 - the wrong sign - and in the design regime byte-identical re-runs moved
ipTM(G12D) by 0.06-0.23, i.e. more than any gap the project has ever produced.
Two independent dockings differ by far more than one side chain is worth.

What replaces it is a *paired* calculation on **one** pose:

    ddG_sel = dG_bind(G12D pose) - dG_bind(same pose, peptide position 6 -> Gly)

Both terms come from the same predicted complex; the only difference between
them is the three side-chain atoms CB/CG/OD1/OD2 of peptide residue 6 and the
local relaxation around them. Everything the two states share - the docking
geometry, the HLA groove, beta2m, the rest of the peptide - cancels to machine
precision in the subtraction instead of being re-sampled.

**A negative `ddG_sel` means the binder prefers G12D**, by that many kcal/mol.

Method
------
* Energy model: Amber ff14SB (`amber14-all.xml`) with GBn2 implicit solvent
  (`implicit/gbn2.xml`, ACE non-polar surface term included), NoCutoff.
  Single-trajectory MM-GBSA: the complex, the isolated binder and the isolated
  target are evaluated at the *same* coordinates, so
  `dG_bind = E_complex - E_binder - E_target` keeps only interaction and
  desolvation terms and has no conformational-entropy term. For a *difference*
  of two dG_bind values that share a pose, the usual single-trajectory
  criticisms (no receptor relaxation, no entropy) largely cancel.
* Local truncation: only residues with a heavy atom within `LOCAL_RADIUS` of
  peptide position 6 are retained, each contiguous run of residues becoming its
  own chain with proper termini. The retained set is computed **once**, on the
  G12D coordinates, and reused for every state, so the two states of a pair are
  built from an identical residue list. Everything dropped is >= 16 A from the
  mutation and contributes equally to both states.
* Relaxation: restrained L-BFGS minimisation. Heavy atoms further than
  `FREE_RADIUS` from position 6 are harmonically restrained to their input
  coordinates (`RESTRAINT_K`); atoms inside that shell, and all hydrogens, move
  freely. This is the "repack or minimise locally" step - it lets the binder
  side chains and the peptide backbone respond to the mutation without letting
  the pose drift.
* Replicates: `replicate` seeds a Gaussian jitter (`JITTER_SIGMA`) applied to
  the free atoms before minimisation, identically in both states of the pair.
  Replicate 0 is unjittered. The spread over replicates is the reproducibility
  figure that the old ipTM gap failed.

Public API
----------
    groups_from_model(model, ...)      -> ComplexGroups
    paired_ddg(path, ..., replicates)  -> PairedResult
    binding_energy(...)                -> one dG_bind

Units are kcal/mol throughout the public API.
"""
from __future__ import annotations

import io
import math
import os
import random
import sys
from dataclasses import dataclass, field, asdict
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
from Bio.PDB import PDBIO, Chain as BioChain, Model as BioModel, Structure as BioStructure
from Bio.PDB.Residue import Residue

import openmm as mm
import openmm.app as app
from openmm import unit
from pdbfixer import PDBFixer

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import complex_analysis as ca  # noqa: E402

# --- knobs ------------------------------------------------------------------
LOCAL_RADIUS = 16.0        # A, residues kept around peptide position 6
FREE_RADIUS = 8.0          # A, heavy atoms minimised without restraint
RESTRAINT_K = 10.0         # kcal/mol/A^2 on the restrained heavy atoms
MIN_ITERATIONS = 500       # L-BFGS iterations
MIN_TOLERANCE = 1.0        # kJ/mol/nm
JITTER_SIGMA = 0.01        # nm, replicate perturbation of the free atoms
MIN_FRAGMENT_LENGTH = 2    # residues; shorter truncation fragments are dropped
FORCEFIELD = ("amber14-all.xml", "implicit/gbn2.xml")
KJ_PER_KCAL = 4.184
SOLUTE_DIELECTRIC = float(os.environ.get("DDG_SOLUTE_DIELECTRIC", "1.0"))

BACKBONE_KEEP = ("N", "CA", "C", "O", "OXT")


class DdgError(RuntimeError):
    """The ddG calculation could not be carried out."""


# --- group bookkeeping ------------------------------------------------------

@dataclass
class ComplexGroups:
    """Which chains of a parsed model are the binder, and which the target."""
    binder_chains: Tuple[str, ...]
    target_chains: Tuple[str, ...]
    peptide_chain: Optional[str]
    peptide_position: int = 6

    def as_dict(self) -> Dict[str, object]:
        return asdict(self)


def groups_from_model(
    model,
    binder_sequence: Optional[str] = None,
    binder_chains: Optional[Sequence[str]] = None,
    peptide_position: int = 6,
) -> ComplexGroups:
    """Identify binder/target chains, by content where possible.

    `binder_chains` overrides the content mapping, which is needed for the 7OW6
    control: its binder is the two-chain TCR (alpha + beta), which
    `complex_analysis.map_complex_chains` - written for one-chain designs -
    cannot express.
    """
    inventory = ca.chain_inventory(model)
    if binder_chains:
        binder = tuple(binder_chains)
        missing = [c for c in binder if c not in inventory]
        if missing:
            raise DdgError(f"binder chains absent from the model: {missing}")
        peptide = None
        for cid, info in inventory.items():
            if cid in binder:
                continue
            if info.length <= ca.PEPTIDE_MAX_LENGTH:
                peptide = cid
        target = tuple(c for c in inventory if c not in binder)
    else:
        cmap = ca.map_complex_chains(model, binder_sequence=binder_sequence)
        binder = (cmap.binder,)
        peptide = cmap.peptide
        target = tuple(c for c in (cmap.hla, cmap.b2m, cmap.peptide) if c)
    if peptide is None:
        raise DdgError("no peptide chain found; a peptide-free state has no ddG")
    return ComplexGroups(binder, target, peptide, peptide_position)


def polymer_residues(chain) -> List[Residue]:
    return [r for r in chain if ca.is_polymer_residue(r)]


def _heavy(residue) -> List:
    return [a for a in residue if a.element != "H" and not a.get_name().startswith("H")]


def peptide_site_residue(model, groups: ComplexGroups) -> Residue:
    """The peptide residue being mutated (1-based position within the peptide)."""
    res = polymer_residues(model[groups.peptide_chain])
    idx = groups.peptide_position - 1
    if idx >= len(res):
        raise DdgError(
            f"peptide chain {groups.peptide_chain} has {len(res)} residues, "
            f"no position {groups.peptide_position}"
        )
    return res[idx]


# --- truncation -------------------------------------------------------------

_CHAIN_ALPHABET = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"


@dataclass
class LocalSystem:
    """A truncated local complex, ready to be written out as PDB."""
    structure: BioStructure.Structure
    binder_new_chains: Tuple[str, ...]
    target_new_chains: Tuple[str, ...]
    site_new_chain: str
    site_resseq: int
    n_residues: int
    n_fragments: int
    center: np.ndarray

    def pdb_text(self) -> str:
        writer = PDBIO()
        writer.set_structure(self.structure)
        buf = io.StringIO()
        writer.save(buf)
        return buf.getvalue()


def build_local_system(
    model,
    groups: ComplexGroups,
    radius: float = LOCAL_RADIUS,
    mutate_site_to_gly: bool = False,
) -> LocalSystem:
    """Cut a sphere of radius `radius` around peptide position 6.

    The kept-residue set is decided on the *unmutated* coordinates so that both
    states of a pair contain exactly the same residues. When
    `mutate_site_to_gly` is true the site residue is rebuilt as glycine in
    place: resname GLY, side chain deleted, backbone untouched. The peptide
    position-6 substitution this project cares about (Asp -> Gly) is purely a
    deletion, so no side-chain rebuilding - and therefore no rotamer guess - is
    involved.
    """
    site = peptide_site_residue(model, groups)
    center_atoms = np.array([a.get_coord() for a in _heavy(site)], dtype=float)
    if len(center_atoms) == 0:
        raise DdgError("the mutation site has no heavy atoms")

    out = BioStructure.Structure("local")
    out_model = BioModel.Model(0)
    out.add(out_model)

    binder_new: List[str] = []
    target_new: List[str] = []
    site_new_chain = None
    site_resseq = None
    n_res = 0
    next_chain = 0

    for cid in list(groups.binder_chains) + list(groups.target_chains):
        if cid not in {c.id for c in model}:
            continue
        residues = polymer_residues(model[cid])
        kept_idx = []
        for i, res in enumerate(residues):
            coords = np.array([a.get_coord() for a in _heavy(res)], dtype=float)
            if len(coords) == 0:
                continue
            d = np.linalg.norm(
                coords[:, None, :] - center_atoms[None, :, :], axis=-1).min()
            if d <= radius:
                kept_idx.append(i)
        # contiguous fragments
        fragments: List[List[int]] = []
        for i in kept_idx:
            if fragments and i == fragments[-1][-1] + 1:
                fragments[-1].append(i)
            else:
                fragments.append([i])
        for frag in fragments:
            has_site = any(residues[i] is site for i in frag)
            if len(frag) < MIN_FRAGMENT_LENGTH and not has_site:
                continue
            if next_chain >= len(_CHAIN_ALPHABET):
                raise DdgError("too many truncation fragments for one PDB file")
            new_id = _CHAIN_ALPHABET[next_chain]
            next_chain += 1
            new_chain = BioChain.Chain(new_id)
            out_model.add(new_chain)
            for n, i in enumerate(frag, start=1):
                res = residues[i]
                is_site = res is site
                name = "GLY" if (is_site and mutate_site_to_gly) else res.get_resname()
                new_res = Residue((" ", n, " "), name, res.segid)
                for atom in res:
                    if atom.element == "H" or atom.get_name().startswith("H"):
                        continue
                    if is_site and mutate_site_to_gly and atom.get_name() not in BACKBONE_KEEP:
                        continue
                    new_atom = atom.copy()
                    new_atom.set_altloc(" ")
                    new_res.add(new_atom)
                new_chain.add(new_res)
                n_res += 1
                if is_site:
                    site_new_chain, site_resseq = new_id, n
            if cid in groups.binder_chains:
                binder_new.append(new_id)
            else:
                target_new.append(new_id)

    if site_new_chain is None:
        raise DdgError("the mutation site was dropped by truncation (bug)")
    if not binder_new:
        raise DdgError(
            f"no binder residue within {radius:.0f} A of peptide position "
            f"{groups.peptide_position}: this binder does not touch the site"
        )
    return LocalSystem(
        structure=out,
        binder_new_chains=tuple(binder_new),
        target_new_chains=tuple(target_new),
        site_new_chain=site_new_chain,
        site_resseq=site_resseq,
        n_residues=n_res,
        n_fragments=next_chain,
        center=center_atoms,
    )


# --- OpenMM plumbing --------------------------------------------------------

def _prepare(pdb_text: str) -> Tuple[app.Topology, unit.Quantity]:
    fixer = PDBFixer(pdbfile=io.StringIO(pdb_text))
    fixer.findMissingResidues()
    fixer.missingResidues = {}          # never build loops across a truncation
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    fixer.removeHeterogens(keepWater=False)
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.0)
    return fixer.topology, fixer.positions


def _chain_atom_indices(topology: app.Topology) -> Dict[str, List[int]]:
    out: Dict[str, List[int]] = {}
    for chain in topology.chains():
        out.setdefault(chain.id, []).extend(
            a.index for r in chain.residues() for a in r.atoms())
    return out


def _system(topology: app.Topology) -> mm.System:
    """The MM-GBSA system: Amber ff14SB + GBn2, no cutoff.

    `DDG_SOLUTE_DIELECTRIC` raises the solute dielectric above 1. GB with
    eps_in = 1 is known to over-penalise a charge buried at a protein-protein
    interface, which is exactly the situation here (the Asp6 carboxylate), and
    eps_in = 2-4 is the usual MM-GBSA remedy. It is a knob, so it is only ever
    moved against the 7OW6 control, never against a design.
    """
    ff = app.ForceField(*FORCEFIELD)
    kwargs = dict(nonbondedMethod=app.NoCutoff, constraints=None,
                  rigidWater=False)
    if SOLUTE_DIELECTRIC != 1.0:
        kwargs["soluteDielectric"] = SOLUTE_DIELECTRIC
    return ff.createSystem(topology, **kwargs)


def _energy(topology: app.Topology, positions) -> float:
    """Potential energy, kcal/mol, of `topology` at `positions`."""
    system = _system(topology)
    integrator = mm.VerletIntegrator(1.0 * unit.femtosecond)
    context = mm.Context(
        system, integrator, mm.Platform.getPlatformByName("CPU"))
    context.setPositions(positions)
    energy = context.getState(getEnergy=True).getPotentialEnergy()
    value = energy.value_in_unit(unit.kilojoule_per_mole) / KJ_PER_KCAL
    del context, integrator
    return value


def _subset(topology: app.Topology, positions, keep_chains: Sequence[str]):
    """Topology/positions restricted to `keep_chains` (coordinates unchanged)."""
    modeller = app.Modeller(topology, positions)
    drop = [r for c in modeller.topology.chains() if c.id not in keep_chains
            for r in c.residues()]
    modeller.delete(drop)
    return modeller.topology, modeller.positions


def _minimise(
    topology: app.Topology,
    positions,
    center: np.ndarray,
    free_radius: float,
    jitter_sigma: float,
    seed: int,
) -> Tuple[object, float, int]:
    """Restrained minimisation. Returns (positions, energy kcal/mol, n_free)."""
    system = _system(topology)
    pos = np.array(positions.value_in_unit(unit.nanometer), dtype=float)
    centre_nm = center / 10.0
    d = np.linalg.norm(pos[:, None, :] - centre_nm[None, :, :], axis=-1).min(axis=1)
    free_mask = d <= (free_radius / 10.0)
    atoms = list(topology.atoms())
    for atom in atoms:
        if atom.element is not None and atom.element.symbol == "H":
            free_mask[atom.index] = True

    restraint = mm.CustomExternalForce(
        "k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
    restraint.addGlobalParameter(
        "k", RESTRAINT_K * KJ_PER_KCAL * 100.0 * unit.kilojoule_per_mole / unit.nanometer**2)
    for name in ("x0", "y0", "z0"):
        restraint.addPerParticleParameter(name)
    n_restrained = 0
    for atom in atoms:
        if free_mask[atom.index]:
            continue
        restraint.addParticle(atom.index, pos[atom.index].tolist())
        n_restrained += 1
    system.addForce(restraint)

    if jitter_sigma > 0:
        rng = np.random.default_rng(seed)
        noise = rng.normal(0.0, jitter_sigma, size=pos.shape)
        noise[~free_mask] = 0.0
        pos = pos + noise

    integrator = mm.VerletIntegrator(1.0 * unit.femtosecond)
    context = mm.Context(system, integrator,
                         mm.Platform.getPlatformByName("CPU"))
    context.setPositions(pos * unit.nanometer)
    mm.LocalEnergyMinimizer.minimize(
        context, MIN_TOLERANCE, MIN_ITERATIONS)
    state = context.getState(getPositions=True, getEnergy=True)
    out_positions = state.getPositions(asNumpy=True)
    energy = state.getPotentialEnergy().value_in_unit(
        unit.kilojoule_per_mole) / KJ_PER_KCAL
    del context, integrator
    return out_positions, energy, int(free_mask.sum())


# --- the metric -------------------------------------------------------------

@dataclass
class StateResult:
    """One state (G12D or WT) of a pair, at one replicate."""
    state: str
    replicate: int
    dg_bind: float
    e_complex: float
    e_binder: float
    e_target: float
    n_atoms: int
    n_free_atoms: int
    n_residues: int

    def as_dict(self) -> Dict[str, object]:
        return asdict(self)


@dataclass
class PairedResult:
    label: str
    ddg_sel: float
    ddg_sd: float
    ddg_values: List[float]
    dg_g12d: List[float]
    dg_wt: List[float]
    groups: Dict[str, object]
    n_residues: int
    n_fragments: int
    seconds: float
    notes: List[str] = field(default_factory=list)

    def as_dict(self) -> Dict[str, object]:
        return asdict(self)


def binding_energy(
    local: LocalSystem,
    replicate: int = 0,
    free_radius: float = FREE_RADIUS,
    state: str = "state",
) -> StateResult:
    """Minimise one truncated state and return its single-trajectory dG_bind."""
    topology, positions = _prepare(local.pdb_text())
    sigma = 0.0 if replicate == 0 else JITTER_SIGMA
    minimised, _e, n_free = _minimise(
        topology, positions, local.center, free_radius, sigma, seed=1000 + replicate)
    e_complex = _energy(topology, minimised)
    b_top, b_pos = _subset(topology, minimised, local.binder_new_chains)
    t_top, t_pos = _subset(topology, minimised, local.target_new_chains)
    e_binder = _energy(b_top, b_pos)
    e_target = _energy(t_top, t_pos)
    return StateResult(
        state=state,
        replicate=replicate,
        dg_bind=e_complex - e_binder - e_target,
        e_complex=e_complex,
        e_binder=e_binder,
        e_target=e_target,
        n_atoms=topology.getNumAtoms(),
        n_free_atoms=n_free,
        n_residues=local.n_residues,
    )


def paired_ddg(
    model,
    groups: ComplexGroups,
    label: str = "",
    replicates: int = 1,
    radius: float = LOCAL_RADIUS,
    free_radius: float = FREE_RADIUS,
    state_sink=None,
) -> PairedResult:
    """ddG_sel = dG_bind(G12D) - dG_bind(position 6 -> Gly), same pose.

    Negative means the binder prefers the G12D peptide.
    """
    import time
    t0 = time.time()
    g12d = build_local_system(model, groups, radius=radius, mutate_site_to_gly=False)
    wt = build_local_system(model, groups, radius=radius, mutate_site_to_gly=True)
    if g12d.n_residues != wt.n_residues:
        raise DdgError("the two states do not contain the same residues (bug)")

    ddgs: List[float] = []
    dg_g12d: List[float] = []
    dg_wt: List[float] = []
    for rep in range(replicates):
        a = binding_energy(g12d, replicate=rep, free_radius=free_radius, state="g12d")
        b = binding_energy(wt, replicate=rep, free_radius=free_radius, state="wt")
        if state_sink is not None:
            state_sink(label, a)
            state_sink(label, b)
        dg_g12d.append(a.dg_bind)
        dg_wt.append(b.dg_bind)
        ddgs.append(a.dg_bind - b.dg_bind)

    mean = sum(ddgs) / len(ddgs)
    sd = (math.sqrt(sum((v - mean) ** 2 for v in ddgs) / (len(ddgs) - 1))
          if len(ddgs) > 1 else 0.0)
    return PairedResult(
        label=label,
        ddg_sel=mean,
        ddg_sd=sd,
        ddg_values=ddgs,
        dg_g12d=dg_g12d,
        dg_wt=dg_wt,
        groups=groups.as_dict(),
        n_residues=g12d.n_residues,
        n_fragments=g12d.n_fragments,
        seconds=time.time() - t0,
    )


def load_model(path: str):
    """Parse a predicted complex (mmCIF or PDB) into a cleaned Bio.PDB model."""
    text = open(path).read()
    fmt = "mmcif" if path.lower().endswith((".cif", ".mmcif")) else "pdb"
    structure = ca.parse_structure(text, fmt=fmt)
    model = list(structure.get_models())[0]
    return ca._clean_model(model)
