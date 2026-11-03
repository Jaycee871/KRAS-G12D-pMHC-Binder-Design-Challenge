#!/usr/bin/env python3
"""Multi-strategy design batch driver for Challenge 5.

Wraps `src/run_pipeline.py` so one batch can explore several *design
hypotheses* at once instead of a single contig/hotspot/temperature setting.
Every design is tagged with the strategy it came from, so the batch report can
say which hypothesis actually produced specific binders.

What it does, per strategy:

    1. RFdiffusion backbones against the trimmed a1/a2 groove + peptide, with
       that strategy's contig (binder length) and hotspot set.
    2. ProteinMPNN sequences for the binder chain at that strategy's sampling
       temperature, filtered with `run_pipeline.sequence_filter_failures`.
    3. The best `--seqs-scored-per-backbone` sequences per backbone (lowest
       ProteinMPNN score) go to Boltz-2 against all four target states.
    4. `design_gates.evaluate` on the result; CSV + details JSON + structures.

Rate limiting and durability (batch-01 post-mortem)
---------------------------------------------------
Batch-01 ran Boltz-2 from a 6-worker pool, drew 235 HTTP 429s and lost every
one of the 50 finished candidates when the run was interrupted, because the
CSV was only written at the end. Both are fixed here:

* Boltz-2 calls go through the process-wide throttle in `nvidia_client`
  (strictly serial, minimum interval between requests, global token bucket),
  so `--workers` above 1 no longer applies to scoring.
* **Every scored candidate is appended to `checkpoint.jsonl` and its structure
  written to disk the moment it is finished**, with an fsync, so an interrupt
  never costs more than the call in flight. `--resume` picks up from that file
  and `--rebuild` regenerates the CSV from it without any API call.

`--rescore-from candidates.json` skips backbone generation and sequence design
entirely and re-scores an existing candidate list - the way to re-run analysis
over backbones and sequences that have already been paid for.

The NVIDIA API key is read from `NVIDIA_API_KEY` by `nvidia_client`; it is
never taken as an argument, logged or written to disk.
"""
from __future__ import annotations

import argparse
import collections
import json
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional

sys.path.insert(0, str(Path(__file__).resolve().parent))

import complex_analysis
import design_gates
import nvidia_client as nvidia
import pose_filter
import run_pipeline
from nvidia_client import NvidiaApiError, NvidiaClient

REPO_ROOT = Path(__file__).resolve().parents[1]

# --- design hypotheses ------------------------------------------------------
# Hotspots are chain+residue labels on the *trimmed* design target, where the
# peptide keeps chain id C (residues 1-10, position 6 = the G12D aspartate) and
# the a1/a2 groove keeps chain id A.
PEPTIDE_CORE = ["C4", "C5", "C6", "C7", "C8"]
POS6_TIGHT = ["C5", "C6", "C7"]
# a1/a2 rim residues lining the groove over the peptide bulge (HLA-A*11:01
# numbering): the "full groove" control strategy of PIPELINE.md section 0,
# included precisely to show that rim-centred paratopes lose specificity.
GROOVE_RIM = ["A62", "A65", "A72", "A146", "A150", "A155", "A158"]

POS6_CENTRED = POS6_TIGHT  # C5-C7, the mutated Asp6 and its neighbours

# Batch-01 strategies, kept so a batch-01 run can be reproduced
# (`--strategy-set batch01`).
BATCH01_STRATEGIES: List[Dict[str, Any]] = [
    {
        "name": "S1-core-short",
        "rationale": "short minibinder arcing over the peptide; hotspots on the "
                     "solvent-exposed peptide core C4-C8, low MPNN temperature",
        "binder_length": 65,
        "hotspots": PEPTIDE_CORE,
        "sampling_temp": 0.1,
    },
    {
        "name": "S2-core-long",
        "rationale": "longer scaffold on the same peptide-core hotspots, higher "
                     "MPNN temperature for more sequence diversity",
        "binder_length": 95,
        "hotspots": PEPTIDE_CORE,
        "sampling_temp": 0.2,
    },
    {
        "name": "S3-pos6-tight",
        "rationale": "hotspots collapsed onto C5-C7, i.e. the mutated Asp6 and "
                     "its immediate neighbours, to force a position-6-centred "
                     "paratope (Xie et al. 2025 mutation-site hotspot)",
        "binder_length": 80,
        "hotspots": POS6_TIGHT,
        "sampling_temp": 0.15,
    },
    {
        "name": "S4-full-groove",
        "rationale": "control: whole peptide plus a1/a2 rim, the docking mode "
                     "PIPELINE.md section 0 predicts will bind all four states",
        "binder_length": 85,
        "hotspots": ["C1", "C3", "C5", "C6", "C8", "C10"] + GROOVE_RIM,
        "sampling_temp": 0.15,
    },
]

# Batch-02: S1-core-short was the only batch-01 strategy that produced
# peptide-contacting designs, so the whole budget goes to that family - small
# binders around 65 residues, hotspots on the peptide core or collapsed onto
# C5-C7 - spread over *many backbones* rather than many sequences per
# backbone (same-backbone sequences give correlated poses, and the pose is
# what fails).
BATCH02_STRATEGIES: List[Dict[str, Any]] = [
    {
        "name": "S1-core-65",
        "rationale": "the batch-01 S1-core-short setting exactly (65 aa, C4-C8, "
                     "T=0.1): the only strategy that produced a peptide-centred "
                     "design, re-run over many more backbones",
        "binder_length": 65,
        "hotspots": PEPTIDE_CORE,
        "sampling_temp": 0.1,
    },
    {
        "name": "S1-core-58",
        "rationale": "shorter S1 variant (58 aa): less scaffold to drape onto "
                     "the a1/a2 helices, so a larger share of the paratope has "
                     "to land on the peptide",
        "binder_length": 58,
        "hotspots": PEPTIDE_CORE,
        "sampling_temp": 0.1,
    },
    {
        "name": "S1-core-72",
        "rationale": "longer S1 variant (72 aa): enough scaffold for a stable "
                     "fold while staying in the batch-01 length band that "
                     "worked",
        "binder_length": 72,
        "hotspots": PEPTIDE_CORE,
        "sampling_temp": 0.1,
    },
    {
        "name": "S1-pos6-65",
        "rationale": "S1 length with hotspots centred on C5-C7, i.e. on the "
                     "mutated Asp6 and its two neighbours, to pull the paratope "
                     "onto the mutation itself",
        "binder_length": 65,
        "hotspots": POS6_CENTRED,
        "sampling_temp": 0.1,
    },
    {
        "name": "S1-pos6-68",
        "rationale": "C5-C7 hotspots at 68 aa and a slightly higher MPNN "
                     "temperature, for sequence diversity over the same "
                     "position-6-centred pose",
        "binder_length": 68,
        "hotspots": POS6_CENTRED,
        "sampling_temp": 0.15,
    },
]

STRATEGY_SETS = {
    "batch01": BATCH01_STRATEGIES,
    "batch02": BATCH02_STRATEGIES,
    "all": BATCH02_STRATEGIES + BATCH01_STRATEGIES,
}

STRATEGIES: List[Dict[str, Any]] = BATCH02_STRATEGIES

EXTRA_COLUMNS = [
    "strategy",
    "backbone",
    "binder_length",
    "hotspots",
    "mpnn_temp",
    "mpnn_score",
    # Where the binder actually sits: batch-01 could not tell "no peptide
    # contact" from "the analysis is broken", so the breakdown is now a column.
    "buried_area_total",
    "buried_area_peptide",
    "buried_area_hla",
    "contacts_peptide",
    "contacts_hla",
    # Does the binder read the Asp6 carboxylate? (batch-02 gate, replaces
    # the unreachable asp6_buried_frac >= 0.60.)
    "asp6_hbond",
    "asp6_salt_bridge",
    "asp6_min_donor_distance",
    "asp6_hbond_partner",
    # Pose filters: the designed RFdiffusion pose, and whether Boltz-2 kept it.
    "design_peptide_contacts",
    "design_peptide_residues",
    "design_pos6_contacts",
    "design_pos6_distance",
    "pose_ok",
    "pose_reject_reason",
    "boltz_calls",
    "analysis_error",
    "n_gate_failures",
    "gate_failures",
    "ranking_shortfalls",
]
BATCH_CSV_COLUMNS = run_pipeline.CSV_COLUMNS + EXTRA_COLUMNS


# --- transport bookkeeping --------------------------------------------------

class CountingClient(NvidiaClient):
    """`NvidiaClient` that records every status it had to back off on."""

    retry_statuses: "collections.Counter" = collections.Counter()
    _lock = threading.Lock()

    def _backoff(self, attempt: int, response=None) -> float:
        if response is not None:
            with CountingClient._lock:
                CountingClient.retry_statuses[response.status_code] += 1
            print(f"[retry] HTTP {response.status_code} on attempt {attempt}; "
                  f"backing off", file=sys.stderr)
        return super()._backoff(attempt, response)


_local = threading.local()


MAX_RETRIES = 4


def client(timeout: float) -> CountingClient:
    """One client per worker thread (requests.Session is not shared)."""
    existing = getattr(_local, "client", None)
    if existing is None:
        existing = CountingClient(timeout=timeout, max_retries=MAX_RETRIES)
        _local.client = existing
    return existing


BACKBONE_STATS: "collections.Counter" = collections.Counter()
CALL_TIMES: List[float] = []
CALL_COUNTS: "collections.Counter" = collections.Counter()
_stats_lock = threading.Lock()


def _record(kind: str, seconds: float) -> None:
    with _stats_lock:
        CALL_COUNTS[kind] += 1
        CALL_TIMES.append(seconds)


# --- stages -----------------------------------------------------------------

def make_backbone(strategy: Mapping[str, Any], index: int, target_pdb: str,
                  args) -> Optional[Dict[str, Any]]:
    contigs = run_pipeline.binder_contig(int(strategy["binder_length"]))
    started = time.time()
    try:
        pdb = client(args.timeout).rfdiffusion_generate(
            input_pdb=target_pdb,
            contigs=contigs,
            hotspot_res=list(strategy["hotspots"]),
            diffusion_steps=args.diffusion_steps,
            random_seed=args.seed + index,
        )
    except NvidiaApiError as exc:
        print(f"[rfdiffusion] {strategy['name']} backbone {index} FAILED: {exc}",
              file=sys.stderr)
        return None
    _record("rfdiffusion", time.time() - started)
    record = {
        "strategy": strategy["name"],
        "backbone": f"{strategy['name']}_bb{index:02d}",
        "contigs": contigs,
        "pdb": pdb,
    }
    record.update(backbone_pose(record, skip=args.no_pose_filter))
    return record


# --- pose filter 1: the designed RFdiffusion complex ------------------------

def backbone_pose(backbone: Mapping[str, Any], skip: bool = False) -> Dict[str, Any]:
    """Measure, in place, whether the designed binder sits on the peptide.

    This is the cheap half of the batch-02 budget saving: it costs no API call
    at all and it rejects whole backbones - and therefore every sequence they
    would have produced - before ProteinMPNN or Boltz-2 are asked for anything.
    """
    if skip:
        return {"pose_ok": True, "pose_reject_reason": "", "pose_metrics": {}}
    try:
        metrics = pose_filter.design_pose_metrics(backbone["pdb"])
    except complex_analysis.StructureAnalysisError as exc:
        print(f"[pose] {backbone['backbone']} unreadable: {exc}", file=sys.stderr)
        return {"pose_ok": False,
                "pose_reject_reason": f"design pose unreadable: {exc}",
                "pose_metrics": {}}
    ok, reason = pose_filter.design_pose_ok(metrics)
    if not ok:
        print(f"[pose] drop backbone {backbone['backbone']}: {reason}")
    return {"pose_ok": ok,
            "pose_reject_reason": "" if ok else f"design pose: {reason}",
            "pose_metrics": metrics}


def design_sequences(backbone: Mapping[str, Any], strategy: Mapping[str, Any],
                     args) -> List[Dict[str, Any]]:
    started = time.time()
    try:
        designed = client(args.timeout).proteinmpnn_design(
            input_pdb=backbone["pdb"],
            input_pdb_chains=[run_pipeline.RFD_BINDER_CHAIN],
            num_seq_per_target=args.seqs_per_backbone,
            sampling_temp=[float(strategy["sampling_temp"])],
            use_soluble_model=True,
            random_seed=args.seed,
        )
    except NvidiaApiError as exc:
        print(f"[proteinmpnn] {backbone['backbone']} FAILED: {exc}", file=sys.stderr)
        return []
    _record("proteinmpnn", time.time() - started)

    out: List[Dict[str, Any]] = []
    for record in designed["sequences"]:
        sequence = record["sequence"].split("/")[0]
        failures = run_pipeline.sequence_filter_failures(sequence)
        if failures:
            print(f"[filter] drop {backbone['backbone']} {sequence[:16]}...: "
                  f"{'; '.join(failures)}")
            continue
        pose = backbone.get("pose_metrics") or {}
        out.append({
            "sequence": sequence,
            "mpnn_score": record.get("score"),
            "strategy": strategy["name"],
            "backbone": backbone["backbone"],
            "binder_length": strategy["binder_length"],
            "hotspots": " ".join(strategy["hotspots"]),
            "mpnn_temp": strategy["sampling_temp"],
            "design_peptide_contacts": pose.get("peptide_contacts", ""),
            "design_peptide_residues": pose.get("peptide_residues_contacted", ""),
            "design_pos6_contacts": pose.get("pos6_contacts", ""),
            "design_pos6_distance": pose.get("pos6_distance", ""),
        })
    # Lowest ProteinMPNN score first (negative log-likelihood: lower is better).
    out.sort(key=lambda r: (r["mpnn_score"] is None, r["mpnn_score"] or 0.0))
    return out[: args.seqs_scored_per_backbone]


def _analysis_columns(row: Dict[str, Any], metrics: Mapping[str, Any],
                      chain_map: Mapping[str, Any]) -> None:
    """Copy the interface breakdown into the CSV row."""
    hla = chain_map.get("hla")
    peptide = chain_map.get("peptide")
    by_area = metrics.get("buried_area_by_chain") or {}
    by_contact = metrics.get("contacts_by_chain") or {}
    row["buried_area_total"] = run_pipeline._round(metrics.get("buried_area_total"), 1)
    row["buried_area_peptide"] = run_pipeline._round(
        metrics.get("buried_area_peptide"), 1)
    row["buried_area_hla"] = run_pipeline._round(by_area.get(hla), 1) if hla else ""
    row["contacts_peptide"] = by_contact.get(peptide, "") if peptide else ""
    row["contacts_hla"] = by_contact.get(hla, "") if hla else ""


def score_candidate(candidate: Mapping[str, Any], target_states, boltz_kwargs,
                    args) -> Optional[Dict[str, Any]]:
    """Boltz-2 + section 6 metrics for one sequence.

    The G12D state is predicted first and the binder's *pose* is checked
    (`pose_filter.predicted_pose_ok`) before the WT, apo and decoy states are
    paid for. A prediction that puts the binder on the HLA heavy chain or on
    beta-2m instead of the peptide cannot be mutation-specific whatever its
    ipTM says, so the remaining three calls are skipped - three of four calls
    saved on every off-peptide prediction.
    """
    sequence = candidate["sequence"]
    row: Dict[str, Any] = dict(candidate)
    row["molecule_class"] = "single_chain"
    row["analysis_error"] = ""
    row["boltz_calls"] = 0
    row["pose_ok"] = ""
    row["pose_reject_reason"] = ""
    detail: Dict[str, Any] = {"sequence": sequence, "states": {},
                              "strategy": candidate["strategy"],
                              "backbone": candidate["backbone"]}
    g12d_structure = None

    for state in run_pipeline.STATES:
        polymers = list(target_states[state]) + [{
            "id": run_pipeline.BINDER_CHAIN,
            "molecule_type": "protein",
            "sequence": sequence,
        }]
        started = time.time()
        try:
            payload = client(args.timeout).boltz2_predict(polymers, **boltz_kwargs)
        except NvidiaApiError as exc:
            print(f"[boltz2] {candidate['backbone']} {state} FAILED: {exc}",
                  file=sys.stderr)
            return None
        _record("boltz2", time.time() - started)
        row["boltz_calls"] = int(row["boltz_calls"]) + 1

        binder_index = len(polymers) - 1
        iptm = run_pipeline.binder_interface_iptm(payload, binder_index)
        row[design_gates.IPTM_COLUMNS[state]] = (
            round(iptm, 4) if iptm is not None else ""
        )
        detail["states"][state] = {
            "binder_interface_iptm": iptm,
            "iptm_scores": payload.get("iptm_scores"),
            "complex_plddt": run_pipeline._first(payload, "complex_plddt_scores"),
            "complex_iplddt": run_pipeline._first(payload, "complex_iplddt_scores"),
            "complex_ipde": run_pipeline._first(payload, "complex_ipde_scores"),
            "pair_chains_iptm_scores": payload.get("pair_chains_iptm_scores"),
        }
        if state == "g12d":
            g12d_structure = payload["structures"][0]["structure"]
            # The binder is located by its *sequence*, not by a chain id, and
            # anything missing raises rather than quietly scoring 0.00 - the
            # batch-01 failure mode. A raised error is recorded on the row so
            # it shows up in the CSV instead of looking like a real zero.
            try:
                metrics = complex_analysis.interface_metrics(
                    g12d_structure, binder_sequence=sequence,
                    binder_chain_hint=run_pipeline.BINDER_CHAIN,
                )
            except complex_analysis.StructureAnalysisError as exc:
                detail["metrics_error"] = f"{type(exc).__name__}: {exc}"
                row["analysis_error"] = f"{type(exc).__name__}: {exc}"
                print(f"[analysis] {candidate['backbone']} {type(exc).__name__}: "
                      f"{exc}", file=sys.stderr)
                metrics = {}
            detail["g12d_metrics"] = metrics
            if metrics:
                row["peptide_buried_frac"] = run_pipeline._round(
                    metrics.get("peptide_buried_frac"))
                row["asp6_buried_frac"] = run_pipeline._round(
                    metrics.get("asp6_buried_frac"))
                row["net_interface_charge"] = run_pipeline._round(
                    metrics.get("net_interface_charge"), 1)
                contacts = metrics.get("pos6_contacts")
                row["pos6_contacts"] = int(contacts) if contacts is not None else ""
                _analysis_columns(row, metrics, metrics.get("chain_map", {}))
                detail["notes"] = metrics.get("notes", [])
            if metrics:
                row["asp6_hbond"] = bool(metrics.get("asp6_hbond"))
                row["asp6_salt_bridge"] = bool(metrics.get("asp6_salt_bridge"))
                row["asp6_min_donor_distance"] = (
                    metrics.get("asp6_min_donor_distance") or "")
                partners = metrics.get("asp6_hbond_partners") or []
                row["asp6_hbond_partner"] = (
                    f"{partners[0]['donor_residue']}{partners[0]['donor_resseq']}"
                    f":{partners[0]['donor_atom']}@{partners[0]['distance']}A"
                    if partners else "")
            row["pae_interface_g12d"] = run_pipeline._round(
                detail["states"]["g12d"]["complex_ipde"], 2)

            # --- pose filter 2: did Boltz-2 keep the binder on the peptide? --
            if not args.no_pose_filter:
                ok, reason = pose_filter.predicted_pose_ok(metrics)
                row["pose_ok"] = ok
                detail["pose_ok"] = ok
                if not ok:
                    row["pose_reject_reason"] = f"predicted pose: {reason}"
                    detail["pose_reject_reason"] = row["pose_reject_reason"]
                    print(f"[pose] {candidate['backbone']} off-peptide after "
                          f"Boltz-2 ({reason}); skipping WT/apo/decoy",
                          flush=True)
                    break

    row["specificity_gap"] = run_pipeline._round(design_gates.specificity_gap(row))
    failures = design_gates.evaluate(row)
    if row["pose_reject_reason"]:
        failures.append(row["pose_reject_reason"])
    if row["analysis_error"]:
        failures.append(f"interface analysis failed: {row['analysis_error']}")
    row["gate_failures"] = "; ".join(failures)
    row["n_gate_failures"] = len(failures)
    shortfalls = design_gates.ranking_shortfalls(row)
    row["ranking_shortfalls"] = "; ".join(shortfalls)
    detail["gate_failures"] = failures
    detail["ranking_shortfalls"] = shortfalls
    return {"row": row, "detail": detail, "structure": g12d_structure}


# --- checkpointing ----------------------------------------------------------
# One JSON object per scored candidate, appended and fsynced the moment the
# candidate finishes. Nothing in a batch is ever only in memory.

CHECKPOINT_NAME = "checkpoint.jsonl"
_checkpoint_lock = threading.Lock()


def checkpoint_path(out_dir: Path) -> Path:
    return Path(out_dir) / CHECKPOINT_NAME


def structure_name(sequence: str) -> str:
    import hashlib

    return hashlib.sha1(sequence.encode()).hexdigest()[:12]


def append_checkpoint(out_dir: Path, result: Mapping[str, Any]) -> Path:
    """Write one finished candidate to disk, structure included, durably."""
    out_dir = Path(out_dir)
    structures = out_dir / "scored_structures"
    structures.mkdir(parents=True, exist_ok=True)
    sequence = result["row"]["sequence"]
    structure_file = ""
    if result.get("structure"):
        target = structures / f"{structure_name(sequence)}_g12d.cif"
        target.write_text(result["structure"])
        structure_file = str(target.relative_to(out_dir))
    record = {
        "sequence": sequence,
        "row": result["row"],
        "detail": result["detail"],
        "structure_file": structure_file,
        "scored_at": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
    }
    path = checkpoint_path(out_dir)
    with _checkpoint_lock:
        with path.open("a") as handle:
            handle.write(json.dumps(record, default=str) + "\n")
            handle.flush()
            import os as _os

            _os.fsync(handle.fileno())
    return path


def load_checkpoint(out_dir: Path) -> Dict[str, Dict[str, Any]]:
    """Scored candidates already on disk, keyed by sequence (last wins)."""
    path = checkpoint_path(out_dir)
    done: Dict[str, Dict[str, Any]] = {}
    if not path.exists():
        return done
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line:
            continue
        try:
            record = json.loads(line)
        except json.JSONDecodeError:
            print(f"[checkpoint] skipping a truncated line in {path}",
                  file=sys.stderr)
            continue
        done[record["sequence"]] = record
    return done


def results_from_checkpoint(out_dir: Path) -> List[Dict[str, Any]]:
    out_dir = Path(out_dir)
    results = []
    for record in load_checkpoint(out_dir).values():
        structure = None
        if record.get("structure_file"):
            candidate = out_dir / record["structure_file"]
            if candidate.exists():
                structure = candidate.read_text()
        results.append({"row": record["row"], "detail": record["detail"],
                        "structure": structure})
    return results


def finalize(out_dir: Path, results: List[Dict[str, Any]], args,
             started: float, backbones: int, candidates: int,
             strategies) -> Dict[str, Any]:
    """Write designs.csv / details / structures / run_stats from `results`."""
    import csv as _csv

    out_dir = Path(out_dir)
    results.sort(key=lambda r: (r["row"].get("specificity_gap") or -99,
                                r["row"].get("asp6_buried_frac") or 0),
                 reverse=True)
    for index, result in enumerate(results, start=1):
        result["row"]["name"] = run_pipeline.design_name(
            result["row"]["sequence"], index)
        result["detail"]["name"] = result["row"]["name"]

    csv_path = out_dir / "designs.csv"
    with csv_path.open("w", newline="") as handle:
        writer = _csv.DictWriter(handle, fieldnames=BATCH_CSV_COLUMNS,
                                 extrasaction="ignore")
        writer.writeheader()
        for result in results:
            writer.writerow({c: result["row"].get(c, "")
                             for c in BATCH_CSV_COLUMNS})

    (out_dir / "designs.details.json").write_text(
        json.dumps([r["detail"] for r in results], indent=2, default=str))

    (out_dir / "structures").mkdir(parents=True, exist_ok=True)
    for rank, result in enumerate(results, start=1):
        if result["structure"] is None or rank > args.top_structures:
            continue
        name = result["row"]["name"]
        (out_dir / "structures" / f"{rank:02d}_{name}_g12d.cif").write_text(
            result["structure"])

    elapsed = time.time() - started
    scored_four_state = [r for r in results
                         if int(r["row"].get("boltz_calls") or 0) >= 4]
    pose_survivors = [r for r in results
                      if str(r["row"].get("pose_ok")).lower() == "true"]
    stats = {
        "elapsed_seconds": round(elapsed, 1),
        "backbones": backbones,
        "backbones_generated": BACKBONE_STATS.get("generated", backbones),
        "backbones_pose_ok": BACKBONE_STATS.get("pose_ok", backbones),
        "backbones_pose_rejected": BACKBONE_STATS.get("pose_rejected", 0),
        "raw_sequences": candidates,
        "scored": len(results),
        "scored_four_state": len(scored_four_state),
        "pose_filter_survivors": len(pose_survivors),
        "pose_filter_rejected": len(results) - len(pose_survivors),
        "passing": sum(1 for r in results if not r["row"]["gate_failures"]),
        "analysis_errors": sum(1 for r in results
                               if r["row"].get("analysis_error")),
        "call_counts": dict(CALL_COUNTS),
        "total_calls": sum(CALL_COUNTS.values()),
        "retry_statuses": dict(CountingClient.retry_statuses),
        "throttle": (nvidia.default_throttle().stats()
                     if nvidia.default_throttle() else None),
        "strategies": [{k: v for k, v in s.items()} for s in strategies],
        "boltz_kwargs": {
            "recycling_steps": args.recycling_steps,
            "sampling_steps": args.sampling_steps,
            "diffusion_samples": 1,
            "without_potentials": True,
        },
        "diffusion_steps": args.diffusion_steps,
        "seqs_per_backbone": args.seqs_per_backbone,
        "seqs_scored_per_backbone": args.seqs_scored_per_backbone,
        "pose_filter": {
            "enabled": not args.no_pose_filter,
            "design_min_peptide_contacts": pose_filter.DESIGN_MIN_PEPTIDE_CONTACTS,
            "design_min_peptide_residues": pose_filter.DESIGN_MIN_PEPTIDE_RESIDUES,
            "design_min_pos6_contacts": pose_filter.DESIGN_MIN_POS6_CONTACTS,
            "design_max_pos6_distance": pose_filter.DESIGN_MAX_POS6_DISTANCE,
            "predicted_min_peptide_area": pose_filter.PREDICTED_MIN_PEPTIDE_AREA,
            "predicted_min_pos6_contacts": pose_filter.PREDICTED_MIN_POS6_CONTACTS,
        },
    }
    (out_dir / "run_stats.json").write_text(json.dumps(stats, indent=2,
                                                       default=str))
    return stats


# --- driver -----------------------------------------------------------------

def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--backbones-per-strategy", type=int, default=8)
    parser.add_argument("--seqs-per-backbone", type=int, default=8)
    parser.add_argument("--seqs-scored-per-backbone", type=int, default=4)
    parser.add_argument("--diffusion-steps", type=int, default=15)
    parser.add_argument("--recycling-steps", type=int, default=3)
    parser.add_argument("--sampling-steps", type=int, default=50)
    parser.add_argument("--workers", type=int, default=1,
                        help="Boltz-2 scoring workers. Kept at 1: the hosted "
                             "endpoint 429s under concurrency (see "
                             "scripts/calibrate_boltz2.py) and nvidia_client "
                             "serialises Boltz-2 anyway")
    parser.add_argument("--design-workers", type=int, default=4)
    parser.add_argument("--timeout", type=float, default=600.0)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--out-dir", default=str(REPO_ROOT / "results" / "batch-01"))
    parser.add_argument("--target-dir", default=str(run_pipeline.DEFAULT_TARGET_DIR))
    parser.add_argument("--top-structures", type=int, default=5)
    parser.add_argument("--strategies", default="",
                        help="comma-separated strategy names (default: all)")
    parser.add_argument("--max-retries", type=int, default=4,
                        help="retries per HTTP request after the first attempt")
    parser.add_argument("--rescore-from", default="",
                        help="skip backbone generation and sequence design and "
                             "re-score the candidates in this candidates.json")
    parser.add_argument("--resume", action="store_true",
                        help="skip candidates already in checkpoint.jsonl")
    parser.add_argument("--rebuild", action="store_true",
                        help="rebuild designs.csv and the summary inputs from "
                             "checkpoint.jsonl without making any API call")
    parser.add_argument("--strategy-set", default="batch02",
                        choices=sorted(STRATEGY_SETS),
                        help="which strategy table to run (default: batch02)")
    parser.add_argument("--no-pose-filter", action="store_true",
                        help="disable both pose filters (batch-01 behaviour): "
                             "score every sequence against all four states")
    parser.add_argument("--max-boltz-calls", type=int, default=0,
                        help="stop starting new candidates once this many "
                             "Boltz-2 calls have been made (0 = no cap)")
    parser.add_argument("--time-budget", type=float, default=0.0,
                        help="stop starting new candidates after this many "
                             "seconds of wall clock (0 = no limit)")
    parser.add_argument("--boltz-min-interval", type=float, default=None,
                        help="seconds between Boltz-2 requests (default: the "
                             "calibrated value in nvidia_client)")
    parser.add_argument("--from-backbones", default="",
                        help="reuse the RFdiffusion backbones already written to "
                             "this directory instead of generating new ones; the "
                             "file name encodes the strategy")
    args = parser.parse_args(argv)

    global MAX_RETRIES
    MAX_RETRIES = args.max_retries

    out_dir = Path(args.out_dir)
    (out_dir / "structures").mkdir(parents=True, exist_ok=True)
    (out_dir / "all_structures").mkdir(parents=True, exist_ok=True)

    if args.boltz_min_interval is not None:
        import os as _os

        _os.environ["NVIDIA_BOLTZ2_MIN_INTERVAL"] = str(args.boltz_min_interval)
        nvidia.configure_throttle(nvidia.build_throttle())

    started = time.time()

    if args.rebuild:
        results = results_from_checkpoint(out_dir)
        stats = finalize(out_dir, results, args, started, 0, len(results),
                         STRATEGIES)
        print(f"[batch] rebuilt {stats['scored']} designs from "
              f"{checkpoint_path(out_dir)} without any API call "
              f"({stats['passing']} clear every gate)")
        return 0

    target_dir = Path(args.target_dir)
    target_states = run_pipeline.load_target_states(target_dir)
    target_pdb = run_pipeline.trimmed_design_target(target_dir / "g12d_complex.pdb")

    strategies = STRATEGY_SETS[args.strategy_set]
    if args.strategies:
        wanted = {s.strip() for s in args.strategies.split(",")}
        strategies = [s for s in strategies if s["name"] in wanted]

    boltz_kwargs = {
        "recycling_steps": args.recycling_steps,
        "sampling_steps": args.sampling_steps,
        "diffusion_samples": 1,
        "without_potentials": True,
    }

    print(f"[batch] {len(strategies)} strategies x "
          f"{args.backbones_per_strategy} backbones, "
          f"{args.seqs_per_backbone} sequences each, "
          f"{args.seqs_scored_per_backbone} scored per backbone")

    # stage 1+2: backbones and sequences, in parallel over backbones
    jobs = [(s, i) for s in strategies for i in range(args.backbones_per_strategy)]
    backbones: List[Dict[str, Any]] = []
    candidates: List[Dict[str, Any]] = []

    if args.rescore_from:
        candidates = json.loads(Path(args.rescore_from).read_text())
        if args.strategies:
            wanted = {s["name"] for s in strategies}
            candidates = [c for c in candidates if c["strategy"] in wanted]
        print(f"[batch] re-scoring {len(candidates)} existing candidates from "
              f"{args.rescore_from}; no RFdiffusion or ProteinMPNN calls")
        return score_all(
            candidates, out_dir, target_states, boltz_kwargs, args, started,
            backbones=len({c.get("backbone") for c in candidates}),
            strategies=strategies)

    reused: Dict[str, Dict[str, Any]] = {}
    if args.from_backbones:
        by_name = {s["name"]: s for s in strategies}
        for path in sorted(Path(args.from_backbones).glob("*.pdb")):
            stem = path.stem
            strategy_name = stem.split("_bb")[0]
            if strategy_name not in by_name:
                continue
            record = {
                "strategy": strategy_name,
                "backbone": stem,
                "contigs": run_pipeline.binder_contig(
                    int(by_name[strategy_name]["binder_length"])),
                "pdb": path.read_text(),
            }
            record.update(backbone_pose(record, skip=args.no_pose_filter))
            reused[stem] = record
        jobs = [(by_name[b["strategy"]], b) for b in reused.values()]
        print(f"[batch] reusing {len(jobs)} backbones from {args.from_backbones}")

    def backbone_then_sequences(job):
        strategy, index = job
        if args.from_backbones:
            bb = index  # already a backbone record
        else:
            bb = make_backbone(strategy, index, target_pdb, args)
        if bb is None:
            return None, []
        if not bb.get("pose_ok", True):
            # No ProteinMPNN call, no Boltz-2 calls: the designed binder is not
            # on the peptide, so nothing it encodes can be mutation-specific.
            return bb, []
        return bb, design_sequences(bb, strategy, args)

    with ThreadPoolExecutor(max_workers=args.design_workers) as pool:
        for future in as_completed([pool.submit(backbone_then_sequences, j)
                                    for j in jobs]):
            bb, seqs = future.result()
            if bb is not None:
                BACKBONE_STATS["generated"] += 1
                if bb.get("pose_ok", True):
                    BACKBONE_STATS["pose_ok"] += 1
                    backbones.append(bb)
                else:
                    BACKBONE_STATS["pose_rejected"] += 1
                target_file = (out_dir / "all_structures" /
                               f"{bb['backbone']}.pdb")
                if not target_file.exists():
                    target_file.write_text(bb["pdb"])
            candidates.extend(seqs)
            print(f"[batch] backbones={BACKBONE_STATS['generated']} "
                  f"on-peptide={BACKBONE_STATS['pose_ok']} "
                  f"candidates={len(candidates)}")

    # de-duplicate identical sequences across backbones
    seen = set()
    unique: List[Dict[str, Any]] = []
    for candidate in candidates:
        if candidate["sequence"] in seen:
            continue
        seen.add(candidate["sequence"])
        unique.append(candidate)
    candidates = unique

    (out_dir / "candidates.json").write_text(json.dumps(candidates, indent=2))
    (out_dir / "backbone_poses.json").write_text(json.dumps(
        [{k: v for k, v in bb.items() if k != "pdb"} for bb in backbones],
        indent=2, default=str))
    print(f"[batch] {BACKBONE_STATS['generated']} backbones generated, "
          f"{len(backbones)} passed the design pose filter, "
          f"{len(candidates)} candidates to score "
          f"(<= {len(candidates) * len(run_pipeline.STATES)} Boltz-2 calls, "
          f"1 per off-peptide candidate)")

    return score_all(candidates, out_dir, target_states, boltz_kwargs, args,
                     started, backbones=len(backbones), strategies=strategies)


def score_all(candidates: List[Dict[str, Any]], out_dir: Path, target_states,
              boltz_kwargs, args, started: float, backbones: int,
              strategies) -> int:
    """Score every candidate, checkpointing each one as soon as it is done."""
    out_dir = Path(out_dir)
    already = load_checkpoint(out_dir) if args.resume else {}
    todo = [c for c in candidates if c["sequence"] not in already]
    if already:
        print(f"[batch] resuming: {len(already)} candidates already in "
              f"{checkpoint_path(out_dir).name}, {len(todo)} left to score")
    print(f"[batch] {len(todo)} candidates to score "
          f"({len(todo) * len(run_pipeline.STATES)} Boltz-2 calls, serial)")

    done = 0
    stopped = ""
    if args.workers > 1:
        print(f"[batch] WARNING: --workers {args.workers} requested, but "
              f"nvidia_client serialises Boltz-2; concurrency only adds 429s",
              file=sys.stderr)

    def budget_exhausted() -> str:
        """Why no further candidate should be started, or ''."""
        if args.max_boltz_calls and CALL_COUNTS["boltz2"] >= args.max_boltz_calls:
            return (f"Boltz-2 call budget reached "
                    f"({CALL_COUNTS['boltz2']}/{args.max_boltz_calls})")
        if args.time_budget and (time.time() - started) >= args.time_budget:
            return (f"time budget reached "
                    f"({time.time() - started:.0f}s/{args.time_budget:.0f}s)")
        return ""

    def record_result(result, index: int) -> None:
        # Durability first: on disk before anything else happens.
        append_checkpoint(out_dir, result)
        row = result["row"]
        print(f"[boltz2] {index}/{len(todo)} scored "
              f"{row.get('backbone')} calls={row.get('boltz_calls')} "
              f"pose_ok={row.get('pose_ok')} "
              f"gap={row.get('specificity_gap')} "
              f"pep_buried={row.get('peptide_buried_frac')} "
              f"pos6={row.get('pos6_contacts')} "
              f"asp6_hbond={row.get('asp6_hbond')} "
              f"({time.time() - started:.0f}s elapsed)", flush=True)

    if args.workers <= 1:
        # Serial is the only sane mode (the throttle serialises Boltz-2
        # anyway), and it is the only mode in which a call budget can actually
        # be respected - so the budget checks live here.
        for candidate in todo:
            stopped = budget_exhausted()
            if stopped:
                print(f"[batch] stopping early: {stopped}; "
                      f"{len(todo) - done} candidates left unscored",
                      file=sys.stderr)
                break
            result = score_candidate(candidate, target_states, boltz_kwargs, args)
            done += 1
            if result is None:
                print(f"[boltz2] {done}/{len(todo)} - candidate dropped "
                      f"(API error)", file=sys.stderr)
                continue
            record_result(result, done)
    else:
        with ThreadPoolExecutor(max_workers=args.workers) as pool:
            futures = [pool.submit(score_candidate, c, target_states,
                                   boltz_kwargs, args) for c in todo]
            for future in as_completed(futures):
                result = future.result()
                done += 1
                if result is None:
                    print(f"[boltz2] {done}/{len(todo)} - candidate dropped "
                          f"(API error)", file=sys.stderr)
                    continue
                record_result(result, done)

    results = results_from_checkpoint(out_dir)
    stats = finalize(out_dir, results, args, started, backbones,
                     len(candidates), strategies)
    stats["stopped_early"] = stopped
    stats["unscored_candidates"] = max(0, len(todo) - done)
    (out_dir / "run_stats.json").write_text(json.dumps(stats, indent=2,
                                                       default=str))
    print(f"[batch] done in {stats['elapsed_seconds']:.0f}s: "
          f"{stats['scored']} scored "
          f"({stats['pose_filter_survivors']} survived the pose filter, "
          f"{stats['scored_four_state']} four-state), "
          f"{stats['passing']} clear every gate, "
          f"{stats['analysis_errors']} analysis errors, "
          f"{stats['total_calls']} API calls, "
          f"retries {stats['retry_statuses'] or 'none'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
