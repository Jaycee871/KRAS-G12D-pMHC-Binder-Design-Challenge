#!/usr/bin/env python3
"""Batch-11: peptide-only generation, then unbiased G12D/apo screening.

This batch changes the generative mechanism rather than retuning the old
HLA-plus-peptide shape filters. RFdiffusion and ProteinMPNN see only the bound
KRAS G12D 10-mer. The HLA groove is restored only for a steric-fit check and
for complete-complex Boltz-2 screening.

Stages:
  plan     Print the frozen call budget; no API calls.
  design   RFdiffusion + ProteinMPNN on peptide only.
  screen   Unconstrained G12D + apo, 3 diffusion samples per call.
  confirm  Unconstrained WT + decoy for screen advancers only.

All model-call stages are checkpointed. Re-running does not intentionally
re-spend a completed Boltz-2 cell.
"""
from __future__ import annotations

import argparse
import hashlib
import io
import json
import os
import statistics
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

import numpy as np
from Bio.PDB import NeighborSearch, PDBIO, PDBParser, Select
from Bio.SVDSuperimposer import SVDSuperimposer

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(REPO / "scripts"))

import complex_analysis as ca  # noqa: E402
import run_pipeline  # noqa: E402
import selectivity_panel as panel  # noqa: E402
from nvidia_client import NvidiaApiError, NvidiaClient  # noqa: E402

OUT = REPO / "results" / "batch-11"
TARGET_PDB = REPO / "data" / "targets" / "g12d_complex.pdb"

PEPTIDE_LENGTH = 10
ASP6 = 6
VAL8 = 8
BINDER_CHAIN = "A"
RFD_TARGET_CHAIN = "B"
SEVERE_CLASH_DISTANCE = 2.2
MAX_SEVERE_CLASH_ATOMS = 2
DESIGN_CONTACT_DISTANCE = 5.0
DESIGN_SIDECHAIN_DISTANCE = 4.5

STRATEGIES: Tuple[Dict[str, Any], ...] = (
    {"name": "B11-peptide40", "binder_length": 40,
     "hotspots": ["C5", "C6", "C8"], "sampling_temp": 0.20},
    {"name": "B11-peptide46", "binder_length": 46,
     "hotspots": ["C5", "C6", "C8"], "sampling_temp": 0.20},
    {"name": "B11-peptide52", "binder_length": 52,
     "hotspots": ["C5", "C6", "C8"], "sampling_temp": 0.25},
)

SCREEN_SETTINGS = {
    "recycling_steps": 8,
    "sampling_steps": 50,
    "diffusion_samples": 3,
    "without_potentials": False,
}

DONOR_AAS = set("RKHSTY")
STRONG_DONORS = set("RKH")
HYDROPHOBIC_AROMATIC = set("FWYILVM")
ACIDIC = set("DE")


class PeptideSelect(Select):
    def __init__(self, chain_id: str):
        self.chain_id = chain_id

    def accept_chain(self, chain):
        return int(chain.id == self.chain_id)

    def accept_residue(self, residue):
        return int(ca.is_polymer_residue(residue))


def _polymer_residues(chain) -> List:
    return [r for r in chain if ca.is_polymer_residue(r)]


def _chain_by_length(model, length: int) -> str:
    matches = [c.id for c in model if len(_polymer_residues(c)) == length]
    if len(matches) != 1:
        raise ca.StructureAnalysisError(
            f"expected one {length}-residue chain, found {matches}"
        )
    return matches[0]


def peptide_only_target(target_pdb: Path = TARGET_PDB) -> Tuple[str, str]:
    """Return the crystallographic 10-mer PDB text and its original chain id."""
    structure = PDBParser(QUIET=True).get_structure("target", str(target_pdb))
    model = next(iter(structure))
    peptide_chain = _chain_by_length(model, PEPTIDE_LENGTH)
    writer = PDBIO()
    writer.set_structure(structure)
    handle = io.StringIO()
    writer.save(handle, PeptideSelect(peptide_chain))
    text = handle.getvalue()
    if "ATOM" not in text:
        raise ca.StructureAnalysisError("peptide-only target contains no atoms")
    return text, peptide_chain


def peptide_only_contig(binder_length: int, peptide_chain: str = "C") -> str:
    return f"{peptide_chain}1-{PEPTIDE_LENGTH}/0 {binder_length}-{binder_length}"


def _design_model(pdb_text: str):
    return ca._clean_model(next(iter(ca.parse_structure(pdb_text, fmt="pdb"))))


def _design_chains(model) -> Tuple[Any, Any]:
    chains = {c.id: c for c in model}
    if BINDER_CHAIN not in chains or RFD_TARGET_CHAIN not in chains:
        raise ca.StructureAnalysisError(
            f"RFdiffusion output needs chains A/B, found {sorted(chains)}"
        )
    target = _polymer_residues(chains[RFD_TARGET_CHAIN])
    if len(target) != PEPTIDE_LENGTH:
        raise ca.StructureAnalysisError(
            f"peptide-only RF target has {len(target)} residues, expected 10"
        )
    return chains[BINDER_CHAIN], chains[RFD_TARGET_CHAIN]


def residue_contact_metrics(
    pdb_text: str,
    peptide_position: int,
    cutoff: float = DESIGN_CONTACT_DISTANCE,
) -> Dict[str, Any]:
    """Binder contact to one peptide residue in peptide-only RF output."""
    model = _design_model(pdb_text)
    binder, target_chain = _design_chains(model)
    target = _polymer_residues(target_chain)
    residue = target[peptide_position - 1]
    binder_atoms = ca.heavy_atoms(binder)
    residue_atoms = ca.heavy_atoms(residue)
    sc_names = set(ca.sidechain_atom_names(residue))
    sidechain_atoms = [
        a for a in residue if a.element != "H" and a.get_id() in sc_names
    ]
    if not binder_atoms or not residue_atoms or not sidechain_atoms:
        return {
            "position": peptide_position,
            "residue": residue.get_resname().strip().upper(),
            "min_distance": None,
            "sidechain_min_distance": None,
            "contact_atoms": 0,
            "sidechain_contact_atoms": 0,
            "binder_residues": [],
        }

    search = NeighborSearch(binder_atoms)
    atom_hits = set()
    side_hits = set()
    residue_hits = set()
    for atom in residue_atoms:
        for neighbour in search.search(atom.coord, cutoff):
            atom_hits.add(id(neighbour))
            parent = neighbour.get_parent()
            residue_hits.add(parent.id[1])
    for atom in sidechain_atoms:
        for neighbour in search.search(atom.coord, cutoff):
            side_hits.add(id(neighbour))
            residue_hits.add(neighbour.get_parent().id[1])

    def min_distance(atoms):
        return min(
            float(np.linalg.norm(a.coord - b.coord))
            for a in binder_atoms for b in atoms
        )

    return {
        "position": peptide_position,
        "residue": residue.get_resname().strip().upper(),
        "min_distance": round(min_distance(residue_atoms), 3),
        "sidechain_min_distance": round(min_distance(sidechain_atoms), 3),
        "contact_atoms": len(atom_hits),
        "sidechain_contact_atoms": len(side_hits),
        "binder_residues": sorted(residue_hits),
    }


def residue_shell_positions(
    pdb_text: str,
    peptide_position: int,
    radius: float = 10.0,
) -> List[Dict[str, Any]]:
    """Binder residue positions closest to one peptide side chain."""
    model = _design_model(pdb_text)
    binder, target_chain = _design_chains(model)
    residue = _polymer_residues(target_chain)[peptide_position - 1]
    sc_names = set(ca.sidechain_atom_names(residue))
    target_atoms = [
        a for a in residue if a.element != "H" and a.get_id() in sc_names
    ]
    if not target_atoms:
        return []
    out = []
    for index, binder_residue in enumerate(_polymer_residues(binder), start=1):
        atoms = ca.heavy_atoms(binder_residue)
        if not atoms:
            continue
        distance = min(
            float(np.linalg.norm(a.coord - b.coord))
            for a in atoms for b in target_atoms
        )
        if distance <= radius:
            out.append({"position": index, "distance": round(distance, 3)})
    out.sort(key=lambda x: (x["distance"], x["position"]))
    return out


def _reference_chains(target_pdb: Path = TARGET_PDB):
    structure = PDBParser(QUIET=True).get_structure("reference", str(target_pdb))
    model = next(iter(structure))
    peptide_id = _chain_by_length(model, PEPTIDE_LENGTH)
    hla_candidates = [
        c for c in model
        if 250 <= len(_polymer_residues(c)) <= 320
    ]
    if len(hla_candidates) != 1:
        raise ca.StructureAnalysisError(
            f"expected one HLA heavy chain, found {[c.id for c in hla_candidates]}"
        )
    return model, hla_candidates[0].id, peptide_id


def hla_reinsertion_metrics(
    pdb_text: str,
    target_pdb: Path = TARGET_PDB,
) -> Dict[str, Any]:
    """Align RF peptide to 7OW6, then measure binder/HLA steric fit."""
    design_model = _design_model(pdb_text)
    binder, target_chain = _design_chains(design_model)
    design_peptide = _polymer_residues(target_chain)

    ref_model, hla_id, peptide_id = _reference_chains(target_pdb)
    ref_peptide = _polymer_residues(ref_model[peptide_id])

    mobile = np.array([r["CA"].coord for r in design_peptide], dtype=float)
    reference = np.array([r["CA"].coord for r in ref_peptide], dtype=float)
    if mobile.shape != reference.shape or len(mobile) != PEPTIDE_LENGTH:
        raise ca.StructureAnalysisError("peptide CA alignment shape mismatch")

    sup = SVDSuperimposer()
    sup.set(reference, mobile)
    sup.run()
    rot, tran = sup.get_rotran()

    binder_coords = np.array(
        [a.coord for a in ca.heavy_atoms(binder)], dtype=float
    )
    binder_aligned = np.dot(binder_coords, rot) + tran

    hla_residues = [
        r for r in _polymer_residues(ref_model[hla_id]) if r.id[1] <= 180
    ]
    hla_coords = np.array(
        [a.coord for r in hla_residues for a in ca.heavy_atoms(r)], dtype=float
    )
    if not len(binder_aligned) or not len(hla_coords):
        raise ca.StructureAnalysisError("empty binder or HLA atom set")

    distances = np.linalg.norm(
        binder_aligned[:, None, :] - hla_coords[None, :, :], axis=-1
    )
    minima = distances.min(axis=1)
    return {
        "peptide_alignment_rmsd": round(float(sup.get_rms()), 4),
        "severe_clash_atoms": int(np.sum(minima < SEVERE_CLASH_DISTANCE)),
        "hla_contact_atoms_5A": int(np.sum(minima < 5.0)),
        "min_hla_distance": round(float(minima.min()), 3),
    }


def design_backbone_assessment(
    pos6: Mapping[str, Any],
    pos8: Mapping[str, Any],
    reinsertion: Mapping[str, Any],
) -> Tuple[bool, List[str]]:
    """Frozen batch-11 residue-level and steric-fit gate."""
    reasons: List[str] = []
    for label, metric in (("Asp6", pos6), ("Val8", pos8)):
        if int(metric.get("sidechain_contact_atoms") or 0) < 1:
            reasons.append(f"{label} has no side-chain contact")
        dist = metric.get("sidechain_min_distance")
        if dist is None or float(dist) > DESIGN_SIDECHAIN_DISTANCE:
            reasons.append(
                f"{label} side-chain distance {dist} > {DESIGN_SIDECHAIN_DISTANCE}"
            )
    union = set(pos6.get("binder_residues") or []) | set(
        pos8.get("binder_residues") or []
    )
    if len(union) < 2:
        reasons.append("Asp6/Val8 contacts do not span two binder residues")
    clashes = int(reinsertion.get("severe_clash_atoms") or 0)
    if clashes > MAX_SEVERE_CLASH_ATOMS:
        reasons.append(
            f"{clashes} severe HLA clash atoms > {MAX_SEVERE_CLASH_ATOMS}"
        )
    return not reasons, reasons


def _distance_weight(distance: float) -> float:
    return max(0.0, 1.0 - max(0.0, distance - 3.0) / 7.0)


def paired_sequence_score(
    sequence: str,
    shell6: Sequence[Mapping[str, Any]],
    shell8: Sequence[Mapping[str, Any]],
) -> float:
    score = 0.0
    for entry in shell6:
        pos = int(entry["position"])
        if not 1 <= pos <= len(sequence):
            continue
        aa = sequence[pos - 1]
        w = _distance_weight(float(entry["distance"]))
        if aa in STRONG_DONORS:
            score += 1.5 * w
        elif aa in DONOR_AAS:
            score += 1.0 * w
        elif aa in ACIDIC:
            score -= 1.5 * w
    for entry in shell8:
        pos = int(entry["position"])
        if not 1 <= pos <= len(sequence):
            continue
        aa = sequence[pos - 1]
        w = _distance_weight(float(entry["distance"]))
        if aa in HYDROPHOBIC_AROMATIC:
            score += 1.0 * w
        elif aa in "KRDE":
            score -= 0.5 * w
    return round(score, 4)


def paired_anchor_variant(
    sequence: str,
    shell6: Sequence[Mapping[str, Any]],
    shell8: Sequence[Mapping[str, Any]],
) -> str:
    """One deterministic bidentate variant: Asp6-facing Arg + Val8-facing Tyr."""
    if not shell6 or not shell8:
        return sequence
    p6 = int(shell6[0]["position"])
    p8 = next(
        (int(e["position"]) for e in shell8 if int(e["position"]) != p6),
        None,
    )
    if p8 is None or p6 > len(sequence) or p8 > len(sequence):
        return sequence
    residues = list(sequence)
    residues[p6 - 1] = "R"
    residues[p8 - 1] = "Y"
    return "".join(residues)


def design_one(
    strategy: Mapping[str, Any],
    index: int,
    peptide_pdb: str,
    peptide_chain: str,
    args,
) -> Dict[str, Any]:
    name = f"{strategy['name']}_bb{index:02d}"
    client = NvidiaClient(timeout=args.timeout)
    try:
        pdb = client.rfdiffusion_generate(
            input_pdb=peptide_pdb,
            contigs=peptide_only_contig(
                int(strategy["binder_length"]), peptide_chain
            ),
            hotspot_res=list(strategy["hotspots"]),
            diffusion_steps=args.diffusion_steps,
            random_seed=args.seed + index * 31 + int(strategy["binder_length"]),
        )
    except Exception as exc:  # noqa: BLE001
        return {
            "backbone": name,
            "strategy": strategy["name"],
            "ok": False,
            "stage": "rfdiffusion",
            "error": f"{type(exc).__name__}: {exc}",
            "candidates": [],
        }

    (OUT / "all_structures").mkdir(parents=True, exist_ok=True)
    (OUT / "all_structures" / f"{name}.pdb").write_text(pdb)

    try:
        pos6 = residue_contact_metrics(pdb, ASP6)
        pos8 = residue_contact_metrics(pdb, VAL8)
        fit = hla_reinsertion_metrics(pdb)
        ok, reasons = design_backbone_assessment(pos6, pos8, fit)
    except Exception as exc:  # noqa: BLE001
        return {
            "backbone": name,
            "strategy": strategy["name"],
            "ok": False,
            "stage": "geometry",
            "error": f"{type(exc).__name__}: {exc}",
            "candidates": [],
        }

    record: Dict[str, Any] = {
        "backbone": name,
        "strategy": strategy["name"],
        "binder_length": strategy["binder_length"],
        "hotspots": strategy["hotspots"],
        "ok": ok,
        "reasons": reasons,
        "pos6": pos6,
        "pos8": pos8,
        "hla_reinsertion": fit,
        "candidates": [],
    }
    if not ok:
        return record

    shell6 = residue_shell_positions(pdb, ASP6)
    shell8 = residue_shell_positions(pdb, VAL8)
    record["shell6"] = shell6[:10]
    record["shell8"] = shell8[:10]

    try:
        designed = client.proteinmpnn_design(
            input_pdb=pdb,
            input_pdb_chains=[BINDER_CHAIN],
            num_seq_per_target=args.seqs_per_backbone,
            sampling_temp=[float(strategy["sampling_temp"])],
            use_soluble_model=True,
            random_seed=args.seed + index + int(strategy["binder_length"]),
        )
    except Exception as exc:  # noqa: BLE001
        record["mpnn_error"] = f"{type(exc).__name__}: {exc}"
        return record

    samples: List[Dict[str, Any]] = []
    for entry in designed.get("sequences") or []:
        sequence = str(entry["sequence"]).split("/")[0]
        if run_pipeline.sequence_filter_failures(sequence):
            continue
        score = paired_sequence_score(sequence, shell6, shell8)
        samples.append(
            {
                "sequence": sequence,
                "mpnn_score": entry.get("score"),
                "paired_score": score,
                "variant": "mpnn",
            }
        )

    samples.sort(
        key=lambda x: (
            -float(x["paired_score"]),
            float(x["mpnn_score"]) if x["mpnn_score"] is not None else 9.0,
            x["sequence"],
        )
    )
    picked = samples[: args.seqs_kept_per_backbone]
    if picked:
        graft = paired_anchor_variant(picked[0]["sequence"], shell6, shell8)
        if (
            graft != picked[0]["sequence"]
            and not run_pipeline.sequence_filter_failures(graft)
        ):
            picked.append(
                {
                    "sequence": graft,
                    "mpnn_score": picked[0]["mpnn_score"],
                    "paired_score": paired_sequence_score(graft, shell6, shell8),
                    "variant": "paired-anchor",
                    "parent_sequence": picked[0]["sequence"],
                }
            )

    seen = set()
    for candidate in picked:
        seq = candidate["sequence"]
        if seq in seen:
            continue
        seen.add(seq)
        candidate.update(
            {
                "backbone": name,
                "strategy": strategy["name"],
                "binder_length": strategy["binder_length"],
                "design_pos6": pos6,
                "design_pos8": pos8,
                "hla_reinsertion": fit,
            }
        )
        record["candidates"].append(candidate)
    return record


def design_stage(args) -> int:
    if not TARGET_PDB.exists():
        raise FileNotFoundError(
            f"{TARGET_PDB} missing; run python src/fetch_target.py first"
        )
    OUT.mkdir(parents=True, exist_ok=True)
    peptide_pdb, peptide_chain = peptide_only_target()
    jobs = [
        (strategy, index)
        for strategy in STRATEGIES
        for index in range(args.backbones_per_strategy)
    ]
    records: List[Dict[str, Any]] = []
    started = time.time()
    with ThreadPoolExecutor(max_workers=args.design_workers) as pool:
        futures = [
            pool.submit(
                design_one, strategy, index, peptide_pdb, peptide_chain, args
            )
            for strategy, index in jobs
        ]
        for future in as_completed(futures):
            record = future.result()
            records.append(record)
            print(
                f"[design] {len(records)}/{len(jobs)} "
                f"bidentate={sum(1 for r in records if r.get('ok'))}",
                flush=True,
            )

    records.sort(key=lambda r: r["backbone"])
    candidates: List[Dict[str, Any]] = []
    seen = set()
    for record in records:
        for candidate in record.get("candidates") or []:
            if candidate["sequence"] in seen:
                continue
            seen.add(candidate["sequence"])
            candidates.append(candidate)
    candidates.sort(
        key=lambda c: (
            -float(c.get("paired_score") or 0.0),
            c["backbone"],
            c["sequence"],
        )
    )

    (OUT / "backbone_poses.json").write_text(
        json.dumps(
            [{k: v for k, v in r.items() if k != "candidates"} for r in records],
            indent=2,
            default=str,
        )
    )
    (OUT / "candidates.json").write_text(
        json.dumps(candidates, indent=2, default=str)
    )
    run = {
        "status": "completed",
        "stage": "design",
        "mechanism": "peptide-only RFdiffusion + ProteinMPNN; HLA reinserted after generation",
        "strategies": STRATEGIES,
        "backbones_requested": len(jobs),
        "backbones_bidentate_ok": sum(1 for r in records if r.get("ok")),
        "candidate_sequences": len(candidates),
        "design_gate": {
            "asp6_sidechain_contact_atoms": ">=1",
            "asp6_sidechain_distance": f"<={DESIGN_SIDECHAIN_DISTANCE}",
            "val8_sidechain_contact_atoms": ">=1",
            "val8_sidechain_distance": f"<={DESIGN_SIDECHAIN_DISTANCE}",
            "distinct_binder_residues": ">=2 across Asp6/Val8 contacts",
            "severe_hla_clash_atoms": f"<={MAX_SEVERE_CLASH_ATOMS} at {SEVERE_CLASH_DISTANCE} A",
        },
        "elapsed_seconds": round(time.time() - started, 1),
        "track1_slots_spent": 0,
    }
    (OUT / "design_run.json").write_text(json.dumps(run, indent=2))
    print(json.dumps(run, indent=2))
    return 0


def select_screen_candidates(
    candidates: Sequence[Mapping[str, Any]], limit: int
) -> List[Dict[str, Any]]:
    """Prefer backbone diversity, then fill by the frozen paired score."""
    ranked = sorted(
        (dict(c) for c in candidates),
        key=lambda c: (
            -float(c.get("paired_score") or 0.0),
            c["backbone"],
            c["sequence"],
        ),
    )
    out: List[Dict[str, Any]] = []
    used_backbones = set()
    used_sequences = set()
    for candidate in ranked:
        if candidate["backbone"] in used_backbones:
            continue
        out.append(candidate)
        used_backbones.add(candidate["backbone"])
        used_sequences.add(candidate["sequence"])
        if len(out) >= limit:
            return out
    for candidate in ranked:
        if candidate["sequence"] in used_sequences:
            continue
        out.append(candidate)
        used_sequences.add(candidate["sequence"])
        if len(out) >= limit:
            break
    return out


def _checkpoint_path() -> Path:
    return OUT / "panel_checkpoint.jsonl"


def _read_checkpoint() -> List[Dict[str, Any]]:
    path = _checkpoint_path()
    if not path.exists():
        return []
    return [
        json.loads(line)
        for line in path.read_text().splitlines()
        if line.strip()
    ]


def _append_checkpoint(record: Mapping[str, Any]) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    with _checkpoint_path().open("a") as handle:
        handle.write(json.dumps(dict(record), default=str) + "\n")
        handle.flush()
        os.fsync(handle.fileno())


def _cell_complete(
    rows: Sequence[Mapping[str, Any]], design: str, arm: str
) -> bool:
    draws = {
        int(r["draw"])
        for r in rows
        if r.get("design") == design and r.get("arm") == arm
    }
    return draws >= {0, 1, 2}


def _design_key(sequence: str) -> str:
    return "B11_" + hashlib.sha256(sequence.encode()).hexdigest()[:10]


def predicted_position_sidechain_contacts(
    mmcif_text: str,
    binder_sequence: str,
    arm: str,
    position: int,
) -> Dict[str, Any]:
    if arm == "apo":
        return {"position": position, "sidechain_contacts": None,
                "sidechain_min_distance": None}
    model = ca._clean_model(
        next(iter(ca.parse_structure(mmcif_text, fmt="mmcif")))
    )
    chain_map = ca.map_complex_chains(
        model,
        binder_sequence=binder_sequence,
        binder_chain_hint=run_pipeline.BINDER_CHAIN,
        require_peptide=True,
    )
    residues = _polymer_residues(model[chain_map.peptide])
    residue = residues[position - 1]
    sc_names = set(ca.sidechain_atom_names(residue))
    target_atoms = [
        a for a in residue if a.element != "H" and a.get_id() in sc_names
    ]
    binder_atoms = ca.heavy_atoms(model[chain_map.binder])
    if not target_atoms or not binder_atoms:
        return {
            "position": position,
            "residue": residue.get_resname().strip().upper(),
            "sidechain_contacts": 0,
            "sidechain_min_distance": None,
        }
    search = NeighborSearch(binder_atoms)
    hits = set()
    for atom in target_atoms:
        for neighbour in search.search(atom.coord, ca.CONTACT_CUTOFF):
            hits.add(id(neighbour))
    distance = min(
        float(np.linalg.norm(a.coord - b.coord))
        for a in binder_atoms for b in target_atoms
    )
    return {
        "position": position,
        "residue": residue.get_resname().strip().upper(),
        "sidechain_contacts": len(hits),
        "sidechain_min_distance": round(distance, 3),
    }


def _call_and_record(
    client: NvidiaClient,
    candidate: Mapping[str, Any],
    arm: str,
    states: Mapping[str, Sequence[Mapping[str, Any]]],
) -> int:
    """One three-draw Boltz-2 cell. Returns 1 only for a new API call."""
    sequence = str(candidate["sequence"])
    design = _design_key(sequence)
    responses = OUT / "responses"
    structures_dir = OUT / "structures"
    responses.mkdir(parents=True, exist_ok=True)
    structures_dir.mkdir(parents=True, exist_ok=True)
    response_path = responses / f"{design}_{arm}.json"

    if response_path.exists():
        payload = json.loads(response_path.read_text())
        api_calls = 0
    else:
        polymers = list(states[arm]) + [
            {
                "id": run_pipeline.BINDER_CHAIN,
                "molecule_type": "protein",
                "sequence": sequence,
            }
        ]
        payload = client.boltz2_predict(polymers, **SCREEN_SETTINGS)
        tmp = response_path.with_suffix(".tmp")
        tmp.write_text(json.dumps(payload))
        os.replace(tmp, response_path)
        api_calls = 1

    structures = payload.get("structures") or []
    if len(structures) != 3:
        raise RuntimeError(
            f"{design}/{arm}: expected 3 structures, got {len(structures)}"
        )
    hashes = [
        hashlib.sha1(str(item["structure"]).encode()).hexdigest()
        for item in structures
    ]
    if len(set(hashes)) != 3:
        raise RuntimeError(f"{design}/{arm}: duplicate diffusion structures")

    existing = _read_checkpoint()
    existing_keys = {
        (r["design"], r["arm"], int(r["draw"])) for r in existing
    }
    binder_index = len(states[arm])
    for draw, item in enumerate(structures):
        key = (design, arm, draw)
        if key in existing_keys:
            continue
        text = str(item["structure"])
        metrics = panel.panel_metrics(text, sequence, arm)
        p8 = predicted_position_sidechain_contacts(text, sequence, arm, VAL8)
        metrics["pos8_sidechain_contacts"] = p8["sidechain_contacts"]
        metrics["pos8_sidechain_min_distance"] = p8["sidechain_min_distance"]
        if p8.get("residue"):
            metrics["pos8_residue"] = p8["residue"]
        metrics["iptm"] = run_pipeline.binder_interface_iptm(
            payload, binder_index, sample=draw
        )
        pde = payload.get("complex_ipde_scores") or []
        metrics["pde_interface"] = pde[draw] if draw < len(pde) else None
        stem = f"{design}_{arm}_r{draw}"
        structure_file = structures_dir / f"{stem}.cif"
        structure_file.write_text(text)
        _append_checkpoint(
            {
                "design": design,
                "sequence": sequence,
                "backbone": candidate.get("backbone"),
                "variant": candidate.get("variant"),
                "paired_score": candidate.get("paired_score"),
                "arm": arm,
                "draw": draw,
                "sampling": SCREEN_SETTINGS,
                "constraint": None,
                "structure_sha1": hashes[draw],
                "row": metrics,
                "structure_file": str(structure_file.relative_to(OUT)),
                "scored_at": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
            }
        )
    return api_calls


def summarize_screen(
    rows: Sequence[Mapping[str, Any]],
    selected: Sequence[Mapping[str, Any]],
) -> Dict[str, Any]:
    summary: Dict[str, Any] = {}
    for candidate in selected:
        sequence = str(candidate["sequence"])
        design = _design_key(sequence)
        g = [
            r for r in rows
            if r.get("design") == design and r.get("arm") == "g12d"
        ]
        a = [
            r for r in rows
            if r.get("design") == design and r.get("arm") == "apo"
        ]
        complete = len(g) == 3 and len(a) == 3
        if not complete:
            summary[design] = {
                "backbone": candidate.get("backbone"),
                "complete": False,
                "advance_to_wt_decoy": False,
            }
            continue
        g_gate = sum(bool(r["row"].get("asp6_gate_pass")) for r in g)
        g_p8 = sum(
            int(r["row"].get("pos8_sidechain_contacts") or 0) >= 1
            for r in g
        )
        g_area = statistics.median(
            float(r["row"]["buried_area_total"]) for r in g
        )
        a_area = statistics.median(
            float(r["row"]["buried_area_total"]) for r in a
        )
        ratio = a_area / g_area if g_area > 0 else None
        advance = bool(
            g_gate >= 2
            and g_p8 >= 2
            and ratio is not None
            and ratio < 1.0
        )
        summary[design] = {
            "backbone": candidate.get("backbone"),
            "variant": candidate.get("variant"),
            "paired_score": candidate.get("paired_score"),
            "complete": True,
            "g12d_Asp6_gate_passes": g_gate,
            "g12d_Val8_contact_draws": g_p8,
            "g12d_median_total_buried_area": round(g_area, 3),
            "apo_median_total_buried_area": round(a_area, 3),
            "apo_to_g12d_total_area_ratio": round(ratio, 6),
            "advance_to_wt_decoy": advance,
        }
    return summary


def screen_stage(args) -> int:
    candidates_path = OUT / "candidates.json"
    if not candidates_path.exists():
        raise FileNotFoundError("batch-11 candidates.json is missing")
    candidates = json.loads(candidates_path.read_text())
    selected = select_screen_candidates(candidates, args.max_screen_candidates)
    states = run_pipeline.load_target_states(REPO / "data" / "targets")
    rows = _read_checkpoint()
    client = NvidiaClient(timeout=args.timeout)
    calls = 0
    started = time.time()

    for candidate in selected:
        design = _design_key(candidate["sequence"])
        for arm in ("g12d", "apo"):
            rows = _read_checkpoint()
            if _cell_complete(rows, design, arm):
                continue
            if calls >= args.max_boltz_calls:
                raise RuntimeError(
                    f"screen call cap reached before completion: {calls}"
                )
            calls += _call_and_record(client, candidate, arm, states)

    rows = _read_checkpoint()
    summary = summarize_screen(rows, selected)
    run = {
        "status": "completed",
        "stage": "screen",
        "selected_candidates": [
            {
                "design": _design_key(c["sequence"]),
                "sequence": c["sequence"],
                "backbone": c.get("backbone"),
                "variant": c.get("variant"),
                "paired_score": c.get("paired_score"),
            }
            for c in selected
        ],
        "arms": ["g12d", "apo"],
        "sampling": SCREEN_SETTINGS,
        "constraint": None,
        "advance_rule": (
            "G12D corrected Asp6 >=2/3 AND G12D Val8 side-chain contact >=2/3 "
            "AND apo/G12D median total buried area <1.0"
        ),
        "screen_summary": summary,
        "advancers": [
            key for key, value in summary.items()
            if value.get("advance_to_wt_decoy")
        ],
        "boltz2_calls_this_run": calls,
        "elapsed_seconds": round(time.time() - started, 1),
        "historical_gates_changed": False,
        "track1_slots_spent": 0,
    }
    (OUT / "screen_run.json").write_text(json.dumps(run, indent=2))
    print(json.dumps(run, indent=2))
    return 0


def _candidate_map() -> Dict[str, Dict[str, Any]]:
    candidates = json.loads((OUT / "candidates.json").read_text())
    return {_design_key(c["sequence"]): c for c in candidates}


def confirm_stage(args) -> int:
    screen_path = OUT / "screen_run.json"
    if not screen_path.exists():
        raise FileNotFoundError("screen_run.json missing")
    screen = json.loads(screen_path.read_text())
    advancers = list(screen.get("advancers") or [])
    candidate_map = _candidate_map()
    states = run_pipeline.load_target_states(REPO / "data" / "targets")
    client = NvidiaClient(timeout=args.timeout)
    calls = 0
    started = time.time()

    for design in advancers:
        candidate = candidate_map[design]
        for arm in ("wt", "decoy"):
            rows = _read_checkpoint()
            if _cell_complete(rows, design, arm):
                continue
            if calls >= args.max_boltz_calls:
                raise RuntimeError(
                    f"confirm call cap reached before completion: {calls}"
                )
            calls += _call_and_record(client, candidate, arm, states)

    rows = _read_checkpoint()
    confirm_summary: Dict[str, Any] = {}
    for design in advancers:
        item: Dict[str, Any] = {}
        for arm in ("wt", "decoy"):
            cell = [
                r for r in rows
                if r.get("design") == design and r.get("arm") == arm
            ]
            if len(cell) != 3:
                item[arm] = {"complete": False}
                continue
            item[arm] = {
                "complete": True,
                "median_total_buried_area": round(
                    statistics.median(
                        float(r["row"]["buried_area_total"]) for r in cell
                    ),
                    3,
                ),
                "median_peptide_buried_area": round(
                    statistics.median(
                        float(r["row"].get("buried_area_peptide") or 0.0)
                        for r in cell
                    ),
                    3,
                ),
                "acidic_donor_contact_draws": sum(
                    any(
                        c.get("peptide_residue") in {"ASP", "GLU"}
                        for c in (r["row"].get("donor_contacts") or [])
                    )
                    for r in cell
                ),
            }
        confirm_summary[design] = item

    run = {
        "status": "completed",
        "stage": "confirm",
        "advancers": advancers,
        "arms": ["wt", "decoy"],
        "sampling": SCREEN_SETTINGS,
        "constraint": None,
        "boltz2_calls_this_run": calls,
        "confirm_summary": confirm_summary,
        "interpretation": (
            "diagnostic only; no post-hoc WT/decoy hard threshold added"
        ),
        "elapsed_seconds": round(time.time() - started, 1),
        "historical_gates_changed": False,
        "track1_slots_spent": 0,
    }
    (OUT / "confirm_run.json").write_text(json.dumps(run, indent=2))
    print(json.dumps(run, indent=2))
    return 0


def plan_stage(args) -> int:
    plan = {
        "backbones": len(STRATEGIES) * args.backbones_per_strategy,
        "rfdiffusion_calls_max": len(STRATEGIES) * args.backbones_per_strategy,
        "proteinmpnn_calls_max": len(STRATEGIES) * args.backbones_per_strategy,
        "phase1_candidates_max": args.max_screen_candidates,
        "phase1_boltz2_calls_max": args.max_screen_candidates * 2,
        "phase1_structures_max": args.max_screen_candidates * 2 * 3,
        "phase2_boltz2_calls_max": args.max_screen_candidates * 2,
        "phase2_only_for_advancers": True,
    }
    print(json.dumps(plan, indent=2))
    return 0


def main(argv=None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--stage", choices=("plan", "design", "screen", "confirm"), required=True
    )
    parser.add_argument("--backbones-per-strategy", type=int, default=8)
    parser.add_argument("--seqs-per-backbone", type=int, default=12)
    parser.add_argument("--seqs-kept-per-backbone", type=int, default=2)
    parser.add_argument("--diffusion-steps", type=int, default=15)
    parser.add_argument("--design-workers", type=int, default=4)
    parser.add_argument("--seed", type=int, default=11105)
    parser.add_argument("--timeout", type=float, default=600.0)
    parser.add_argument("--max-screen-candidates", type=int, default=12)
    parser.add_argument("--max-boltz-calls", type=int, default=24)
    args = parser.parse_args(argv)

    if args.backbones_per_strategy < 1:
        parser.error("--backbones-per-strategy must be >=1")
    if not 1 <= args.max_screen_candidates <= 12:
        parser.error("--max-screen-candidates must be 1..12")

    if args.stage == "plan":
        return plan_stage(args)
    if args.stage == "design":
        return design_stage(args)
    if args.stage == "screen":
        return screen_stage(args)
    return confirm_stage(args)


if __name__ == "__main__":
    raise SystemExit(main())
