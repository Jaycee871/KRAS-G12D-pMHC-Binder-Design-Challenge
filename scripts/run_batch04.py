#!/usr/bin/env python3
"""Batch-04: donor-biased, pos6-centred design under the Asp6 read-out gate.

What is different from batches 02/03, and why
---------------------------------------------
Batches 02 and 03 produced 88 scored designs, **zero** of which put any
hydrogen-bond donor within 3.5 A of the Asp6 carboxylate, with a median net
interface charge of -1.0 (see `results/gate-audit/ASP6_GATE_AUDIT.md`). Both
ranking quantities this project has tried - the ipTM `specificity_gap` and the
paired same-pose MM-GBSA `ddG_sel` - failed the 7OW6 positive control and are
retired. Nothing here ranks on either.

This batch therefore changes the *designing*, not the scoring:

1. RFdiffusion hotspots collapsed onto peptide position 6 (C6, C6-C7, C5-C7)
   and the design-stage pose filter tightened on position 6, so only backbones
   that actually reach the mutated side chain are paid for.
2. **Donor bias in sequence design.** The hosted ProteinMPNN NIM ignores
   `bias_AA_jsonl` (accepted, no effect on identical seeds) and answers HTTP 500
   for `pssm_jsonl`/`bias_by_res_jsonl`, so API-side biasing is not available.
   The bias is applied locally instead, which is measurable and reproducible:
   the binder positions pointing at Asp6 are identified on the designed
   backbone, ProteinMPNN is over-sampled, samples are *selected* on donor
   content at those positions, and a deterministic donor-anchored variant is
   derived (closest position -> Arg, second -> Tyr/Ser, acidic residues in the
   pos6 shell -> their neutral amides) to push the interface charge off
   negative.
3. Boltz-2 is called **constrained** with the validated flat pocket term
   `{"binder":"D","contacts":[{"id":"C","residue_index":6}],"max_distance":6.0}`,
   strictly serially at the calibrated >= 8 s interval, **once per candidate**
   (the G12D state). The other three states existed only to feed
   `specificity_gap`, which is retired; ipTM is used here as a pose-quality
   read-out, never as a gap and never as a tie-breaker.
4. Gate: a binder donor (Arg/Lys/His/Ser/Thr/Tyr) heavy atom within 3.5 A of
   Asp6 OD1/OD2 **and** `asp6_buried_frac > 0`. Survivors are ranked on
   measured geometry only.
"""
from __future__ import annotations

import argparse, json, os, sys, time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any, Dict, List, Optional

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

import complex_analysis, design_gates, pose_filter, run_pipeline  # noqa: E402
from nvidia_client import BOLTZ2_PATH, NvidiaApiError, NvidiaClient  # noqa: E402

OUT = REPO / "results" / "batch-04"

# --- the gate ---------------------------------------------------------------
GATE_DONOR_RESIDUES = {"ARG", "LYS", "HIS", "SER", "THR", "TYR"}
GATE_DONOR_ATOMS = {"NE", "NH1", "NH2", "NZ", "ND1", "NE2", "OG", "OG1", "OH"}
GATE_DONOR_DISTANCE = 3.5

# --- design hypotheses ------------------------------------------------------
STRATEGIES: List[Dict[str, Any]] = [
    {"name": "B4-c6-62", "binder_length": 62, "hotspots": ["C6"],
     "sampling_temp": 0.15,
     "rationale": "hotspot collapsed onto position 6 alone - the strongest "
                  "available pull of the paratope onto the mutated aspartate"},
    {"name": "B4-c67-65", "binder_length": 65, "hotspots": ["C6", "C7"],
     "sampling_temp": 0.2,
     "rationale": "position 6 plus its C-terminal neighbour at the batch-01/02 "
                  "length that produced peptide-centred poses"},
    {"name": "B4-pos6-68", "binder_length": 68, "hotspots": ["C5", "C6", "C7"],
     "sampling_temp": 0.2,
     "rationale": "the S1-pos6-68 setting that produced bb13, the only design "
                  "in three batches to read the carboxylate, re-run with the "
                  "donor bias"},
    {"name": "B4-pos6-72", "binder_length": 72, "hotspots": ["C5", "C6", "C7"],
     "sampling_temp": 0.25,
     "rationale": "longer scaffold on the same pos6 hotspots at a higher MPNN "
                  "temperature, for topological diversity among survivors"},
]

DONOR_AAS = set("RKHSTY")
ACID_TO_AMIDE = {"D": "N", "E": "Q"}

CSV_COLUMNS = [
    "rank", "name", "strategy", "backbone", "variant", "binder_length",
    "sequence", "mpnn_score", "donor_bias_score",
    "iptm_g12d", "pae_interface_g12d",
    "peptide_buried_frac", "asp6_buried_frac",
    "asp6_min_donor_distance", "asp6_donor_partner",
    "asp6_min_any_donor_distance",
    "pos6_contacts", "pos6_sidechain_contacts",
    "buried_area_total", "buried_area_peptide", "buried_area_hla",
    "net_interface_charge", "asp6_hbond", "asp6_salt_bridge",
    "pose_ok", "pose_reject_reason", "gate_pass", "gate_reason",
    "design_pos6_contacts", "design_pos6_distance", "analysis_error",
]


# --- local donor bias -------------------------------------------------------

def pos6_shell_positions(pdb_text: str, radius: float = 10.0,
                         binder_chain: str = "A", target_chain: str = "B",
                         peptide_length: int = 10) -> List[Dict[str, Any]]:
    """Binder positions (1-based) pointing at Asp6, closest first.

    Measured on the RFdiffusion output complex, which is backbone-only on the
    target side; distance is the closest heavy-atom approach of the binder
    residue to the position-6 residue.
    """
    structure = complex_analysis.parse_structure(pdb_text, fmt="pdb")
    model = complex_analysis._clean_model(next(iter(structure)))
    chains = {c.id: c for c in model}
    target = [r for r in chains[target_chain]
              if complex_analysis.is_polymer_residue(r)]
    pos6 = target[-peptide_length:][complex_analysis.MUTATED_PEPTIDE_POSITION - 1]
    pos6_atoms = complex_analysis.heavy_atoms(pos6)
    out = []
    binder = [r for r in chains[binder_chain]
              if complex_analysis.is_polymer_residue(r)]
    for index, residue in enumerate(binder, start=1):
        atoms = complex_analysis.heavy_atoms(residue)
        if not atoms:
            continue
        distance = min(float(((a.coord - b.coord) ** 2).sum() ** 0.5)
                       for a in atoms for b in pos6_atoms)
        if distance <= radius:
            out.append({"position": index, "distance": round(distance, 2)})
    out.sort(key=lambda r: r["distance"])
    return out


def donor_bias_score(sequence: str, shell: List[Dict[str, Any]]) -> float:
    """Donor content of the Asp6 shell, distance-weighted, acid-penalised."""
    score = 0.0
    for entry in shell:
        position = entry["position"]
        if position > len(sequence):
            continue
        weight = max(0.0, 1.0 - (entry["distance"] - 3.0) / 9.0)
        aa = sequence[position - 1]
        if aa in DONOR_AAS:
            score += weight * (1.5 if aa in "RKH" else 1.0)
        elif aa in ACID_TO_AMIDE:
            score -= weight * 1.5
    return round(score, 3)


def donor_anchored_variant(sequence: str, shell: List[Dict[str, Any]]) -> str:
    """Deterministic donor graft: Arg anchor at position 6's nearest residue.

    The nearest shell position becomes Arg, the second nearest Tyr (a neutral
    donor, so the anchor is not two full charges), and every Asp/Glu inside the
    shell becomes its amide - which is the specific correction for the -1.0
    median interface charge of batches 02 and 03.
    """
    residues = list(sequence)
    for rank, entry in enumerate(shell):
        position = entry["position"]
        if position > len(residues):
            continue
        aa = residues[position - 1]
        if rank == 0 and aa not in "RK":
            residues[position - 1] = "R"
        elif rank == 1 and aa not in DONOR_AAS:
            residues[position - 1] = "Y"
        elif aa in ACID_TO_AMIDE:
            residues[position - 1] = ACID_TO_AMIDE[aa]
    return "".join(residues)


# --- stage 1 + 2 ------------------------------------------------------------

def backbone_and_sequences(client: NvidiaClient, strategy, index, target_pdb,
                           args) -> Dict[str, Any]:
    name = f"{strategy['name']}{args.tag}_bb{index:02d}"
    try:
        pdb = client.rfdiffusion_generate(
            input_pdb=target_pdb,
            contigs=run_pipeline.binder_contig(int(strategy["binder_length"])),
            hotspot_res=list(strategy["hotspots"]),
            diffusion_steps=args.diffusion_steps,
            random_seed=args.seed + index * 17,
        )
    except NvidiaApiError as exc:
        print(f"[rfdiffusion] {name} FAILED: {exc}", file=sys.stderr)
        return {"backbone": name, "ok": False, "candidates": []}

    (OUT / "all_structures" / f"{name}.pdb").write_text(pdb)
    try:
        metrics = pose_filter.design_pose_metrics(pdb)
    except complex_analysis.StructureAnalysisError as exc:
        print(f"[pose] {name} unreadable: {exc}", file=sys.stderr)
        return {"backbone": name, "ok": False, "candidates": []}

    ok, reason = pose_filter.design_pose_ok(
        metrics, min_pos6_contacts=args.min_pos6_contacts,
        max_pos6_distance=args.max_pos6_distance)
    record = {"backbone": name, "strategy": strategy["name"], "ok": ok,
              "reason": reason, "pose": metrics, "candidates": []}
    if not ok:
        print(f"[pose] drop {name}: {reason}")
        return record

    shell = pos6_shell_positions(pdb)
    record["pos6_shell"] = shell[:8]
    try:
        designed = client.proteinmpnn_design(
            input_pdb=pdb, input_pdb_chains=[run_pipeline.RFD_BINDER_CHAIN],
            num_seq_per_target=args.seqs_per_backbone,
            sampling_temp=[float(strategy["sampling_temp"])],
            use_soluble_model=True, random_seed=args.seed + index)
    except NvidiaApiError as exc:
        print(f"[proteinmpnn] {name} FAILED: {exc}", file=sys.stderr)
        return record

    samples = []
    for entry in designed["sequences"]:
        sequence = entry["sequence"].split("/")[0]
        if run_pipeline.sequence_filter_failures(sequence):
            continue
        samples.append({"sequence": sequence, "mpnn_score": entry.get("score"),
                        "donor_bias_score": donor_bias_score(sequence, shell)})
    if not samples:
        return record
    # Donor-biased *selection*: best donor content first, MPNN score breaks ties.
    samples.sort(key=lambda s: (-s["donor_bias_score"],
                                s["mpnn_score"] if s["mpnn_score"] is not None else 9))
    picked = samples[: args.seqs_scored_per_backbone]
    for sample in picked:
        sample.update({"strategy": strategy["name"], "backbone": name,
                       "binder_length": strategy["binder_length"],
                       "variant": "mpnn-selected",
                       "design_pos6_contacts": metrics["pos6_contacts"],
                       "design_pos6_distance": metrics["pos6_distance"]})
    # ... plus one deterministic donor-anchored graft off the best sample.
    grafted = donor_anchored_variant(picked[0]["sequence"], shell)
    if (grafted != picked[0]["sequence"]
            and not run_pipeline.sequence_filter_failures(grafted)):
        picked.append({"sequence": grafted, "mpnn_score": picked[0]["mpnn_score"],
                       "donor_bias_score": donor_bias_score(grafted, shell),
                       "strategy": strategy["name"], "backbone": name,
                       "binder_length": strategy["binder_length"],
                       "variant": "donor-anchored",
                       "design_pos6_contacts": metrics["pos6_contacts"],
                       "design_pos6_distance": metrics["pos6_distance"]})
    record["candidates"] = picked
    return record


# --- stage 3: constrained Boltz-2 -------------------------------------------

POCKET = [{"binder": run_pipeline.BINDER_CHAIN,
           "contacts": [{"id": run_pipeline.PEPTIDE_CHAIN, "residue_index": 6}],
           "max_distance": 6.0}]


def boltz2_constrained(client: NvidiaClient, polymers, args) -> Dict[str, Any]:
    body = {"polymers": [dict(p) for p in polymers],
            "recycling_steps": args.recycling_steps,
            "sampling_steps": args.sampling_steps,
            "diffusion_samples": 1, "step_scale": 1.638,
            "without_potentials": False, "constraints": POCKET}
    return client.post(BOLTZ2_PATH, body)


def gate_donor_distance(metrics) -> Dict[str, Any]:
    """Closest gate-qualifying donor (Arg/Lys/His/Ser/Thr/Tyr side chain)."""
    best = None
    for partner in metrics.get("asp6_hbond_partners") or []:
        if (partner["donor_residue"] in GATE_DONOR_RESIDUES
                and partner["donor_atom"] in GATE_DONOR_ATOMS
                and partner["distance"] <= GATE_DONOR_DISTANCE):
            if best is None or partner["distance"] < best["distance"]:
                best = partner
    any_donor = metrics.get("asp6_min_donor_distance")
    return {"donor": best, "any_donor_distance": any_donor}


def pos6_sidechain_contacts(structure_text: str, metrics) -> Optional[int]:
    """Binder heavy atoms within 5 A of the Asp6 **side chain** only.

    Backbone contacts to position 6 do not discriminate G12D from WT - the
    backbone is identical in both - so this is the contact count the ranking
    uses.
    """
    from Bio.PDB import NeighborSearch

    structure = complex_analysis.parse_structure(structure_text, fmt="mmcif")
    model = complex_analysis._clean_model(next(iter(structure)))
    binder_chain = metrics["binder_chain"]
    peptide_chain = metrics["peptide_chain"]
    residues = [r for r in model[peptide_chain]
                if complex_analysis.is_polymer_residue(r)]
    pos6 = residues[5]
    names = set(complex_analysis.sidechain_atom_names(pos6))
    atoms = [a for a in pos6 if a.get_id() in names and a.element != "H"]
    if not atoms:
        return None
    search = NeighborSearch(complex_analysis.heavy_atoms(model[binder_chain]))
    hits = set()
    for atom in atoms:
        for neighbour in search.search(atom.coord, complex_analysis.CONTACT_CUTOFF):
            hits.add(id(neighbour))
    return len(hits)


def score_candidate(client, candidate, target_states, args) -> Dict[str, Any]:
    sequence = candidate["sequence"]
    row = {k: candidate.get(k, "") for k in
           ("strategy", "backbone", "variant", "binder_length", "mpnn_score",
            "donor_bias_score", "design_pos6_contacts", "design_pos6_distance")}
    row["sequence"] = sequence
    row["analysis_error"] = ""
    polymers = list(target_states["g12d"]) + [
        {"id": run_pipeline.BINDER_CHAIN, "molecule_type": "protein",
         "sequence": sequence}]
    payload = boltz2_constrained(client, polymers, args)
    structure = payload["structures"][0]["structure"]
    row["iptm_g12d"] = run_pipeline._round(
        run_pipeline.binder_interface_iptm(payload, len(polymers) - 1), 4)
    row["pae_interface_g12d"] = run_pipeline._round(
        run_pipeline._first(payload, "complex_ipde_scores"), 2)
    try:
        metrics = complex_analysis.interface_metrics(
            structure, binder_sequence=sequence,
            binder_chain_hint=run_pipeline.BINDER_CHAIN)
    except complex_analysis.StructureAnalysisError as exc:
        row["analysis_error"] = f"{type(exc).__name__}: {exc}"
        return {"row": row, "structure": structure, "metrics": {}}

    by_area = metrics.get("buried_area_by_chain") or {}
    hla = (metrics.get("chain_map") or {}).get("hla")
    row.update({
        "peptide_buried_frac": run_pipeline._round(metrics["peptide_buried_frac"]),
        "asp6_buried_frac": run_pipeline._round(metrics["asp6_buried_frac"]),
        "buried_area_total": run_pipeline._round(metrics["buried_area_total"], 1),
        "buried_area_peptide": run_pipeline._round(metrics["buried_area_peptide"], 1),
        "buried_area_hla": run_pipeline._round(by_area.get(hla), 1) if hla else "",
        "pos6_contacts": int(metrics["pos6_contacts"]),
        "net_interface_charge": run_pipeline._round(metrics["net_interface_charge"], 1),
        "asp6_hbond": bool(metrics["asp6_hbond"]),
        "asp6_salt_bridge": bool(metrics["asp6_salt_bridge"]),
    })
    row["pos6_sidechain_contacts"] = pos6_sidechain_contacts(structure, metrics)
    gate = gate_donor_distance(metrics)
    row["asp6_min_any_donor_distance"] = gate["any_donor_distance"] or ""
    row["asp6_min_donor_distance"] = gate["donor"]["distance"] if gate["donor"] else ""
    row["asp6_donor_partner"] = (
        f"{gate['donor']['donor_residue']}{gate['donor']['donor_resseq']}"
        f":{gate['donor']['donor_atom']}@{gate['donor']['distance']}A"
        if gate["donor"] else "")
    ok, reason = pose_filter.predicted_pose_ok(metrics)
    row["pose_ok"] = ok
    row["pose_reject_reason"] = "" if ok else reason

    reasons = []
    if not gate["donor"]:
        reasons.append(
            f"no Arg/Lys/His/Ser/Thr/Tyr donor within {GATE_DONOR_DISTANCE} A of "
            f"Asp6 OD1/OD2 (closest donor of any kind: "
            f"{gate['any_donor_distance'] or 'none within 4.0 A'})")
    if not (metrics["asp6_buried_frac"] > 0):
        reasons.append("asp6_buried_frac = 0")
    row["gate_pass"] = not reasons
    row["gate_reason"] = "; ".join(reasons)
    return {"row": row, "structure": structure, "metrics": metrics}


# --- driver -----------------------------------------------------------------

def append_checkpoint(record: Dict[str, Any]) -> None:
    path = OUT / "checkpoint.jsonl"
    with path.open("a") as handle:
        handle.write(json.dumps(record, default=str) + "\n")
        handle.flush()
        os.fsync(handle.fileno())


def main(argv=None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--backbones-per-strategy", type=int, default=24)
    parser.add_argument("--seqs-per-backbone", type=int, default=16)
    parser.add_argument("--seqs-scored-per-backbone", type=int, default=2)
    parser.add_argument("--diffusion-steps", type=int, default=15)
    parser.add_argument("--recycling-steps", type=int, default=3)
    parser.add_argument("--sampling-steps", type=int, default=50)
    parser.add_argument("--design-workers", type=int, default=4)
    parser.add_argument("--seed", type=int, default=4042)
    parser.add_argument("--timeout", type=float, default=600.0)
    parser.add_argument("--min-pos6-contacts", type=int, default=2)
    parser.add_argument("--max-pos6-distance", type=float, default=5.0)
    parser.add_argument("--max-boltz-calls", type=int, default=160)
    parser.add_argument("--score-time-budget", type=float, default=2700.0)
    parser.add_argument("--replicates", type=int, default=2)
    parser.add_argument("--max-replicated", type=int, default=8)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--rescore-from", default="",
                        help="skip stages 1-2 and score the candidates in this "
                             "candidates.json (no RFdiffusion, no ProteinMPNN)")
    parser.add_argument("--tag", default="",
                        help="suffix for backbone names and artefact "
                             "files, so a second wave does not "
                             "overwrite the first")
    args = parser.parse_args(argv)

    started = time.time()
    for sub in ("all_structures", "scored_structures"):
        (OUT / sub).mkdir(parents=True, exist_ok=True)

    target_dir = REPO / "data" / "targets"
    target_states = run_pipeline.load_target_states(target_dir)
    target_pdb = run_pipeline.trimmed_design_target(target_dir / "g12d_complex.pdb")

    client = NvidiaClient(timeout=args.timeout)

    # stage 1 + 2
    backbone_records: List[Dict[str, Any]] = []
    if args.rescore_from:
        candidates = json.loads(Path(args.rescore_from).read_text())
        print(f"[design] skipping stages 1-2: scoring {len(candidates)} "
              f"candidates from {args.rescore_from}")
        return score_stage(client, candidates, backbone_records, target_states,
                           args, started)
    candidates: List[Dict[str, Any]] = []
    jobs = [(s, i) for s in STRATEGIES for i in range(args.backbones_per_strategy)]
    with ThreadPoolExecutor(max_workers=args.design_workers) as pool:
        futures = [pool.submit(backbone_and_sequences, NvidiaClient(timeout=args.timeout),
                               s, i, target_pdb, args) for s, i in jobs]
        for future in as_completed(futures):
            record = future.result()
            backbone_records.append(record)
            candidates.extend(record.get("candidates") or [])
            print(f"[design] backbones={len(backbone_records)}/{len(jobs)} "
                  f"on-pos6={sum(1 for r in backbone_records if r['ok'])} "
                  f"candidates={len(candidates)} "
                  f"({time.time() - started:.0f}s)", flush=True)

    seen, unique = set(), []
    for candidate in candidates:
        if candidate["sequence"] in seen:
            continue
        seen.add(candidate["sequence"])
        unique.append(candidate)
    candidates = unique
    (OUT / f"candidates{args.tag}.json").write_text(json.dumps(candidates, indent=2))
    (OUT / f"backbone_poses{args.tag}.json").write_text(json.dumps(
        [{k: v for k, v in r.items() if k != "candidates"}
         for r in backbone_records], indent=2, default=str))
    print(f"[design] stage 1-2 done in {time.time() - started:.0f}s: "
          f"{len(backbone_records)} backbones, "
          f"{sum(1 for r in backbone_records if r['ok'])} reach position 6, "
          f"{len(candidates)} candidates")

    return score_stage(client, candidates, backbone_records, target_states,
                       args, started)


def score_stage(client, candidates, backbone_records, target_states, args,
                started) -> int:
    """Stage 3: constrained Boltz-2, one call per candidate, strictly serial."""
    done_sequences = set()
    if args.resume and (OUT / "checkpoint.jsonl").exists():
        done_sequences = {json.loads(l)["row"]["sequence"]
                          for l in (OUT / "checkpoint.jsonl").read_text().splitlines() if l.strip()}
    # donor-biased ordering: score the most donor-rich candidates first, so a
    # budget stop never costs the candidates most likely to pass the gate.
    candidates.sort(key=lambda c: -float(c.get("donor_bias_score") or 0))
    todo = [c for c in candidates if c["sequence"] not in done_sequences]
    calls = 0
    scored = 0
    stopped = ""
    score_started = time.time()
    for candidate in todo:
        if calls >= args.max_boltz_calls:
            stopped = f"Boltz-2 call cap reached ({calls}/{args.max_boltz_calls})"
            break
        if time.time() - score_started >= args.score_time_budget:
            stopped = f"scoring time budget reached ({time.time() - score_started:.0f}s)"
            break
        try:
            result = score_candidate(client, candidate, target_states, args)
        except (NvidiaApiError, Exception) as exc:  # noqa: BLE001
            print(f"[boltz2] {candidate['backbone']} dropped: "
                  f"{type(exc).__name__}: {exc}", file=sys.stderr)
            calls += 1
            continue
        calls += 1
        scored += 1
        name = f"{candidate['backbone']}_{candidate['variant']}"
        (OUT / "scored_structures" / f"{name}_g12d.cif").write_text(result["structure"])
        append_checkpoint({"row": result["row"], "structure_file":
                           f"scored_structures/{name}_g12d.cif",
                           "scored_at": time.strftime("%Y-%m-%dT%H:%M:%S%z")})
        row = result["row"]
        print(f"[boltz2] {scored}/{len(todo)} {name} iptm={row.get('iptm_g12d')} "
              f"pep_buried={row.get('peptide_buried_frac')} "
              f"pos6_sc={row.get('pos6_sidechain_contacts')} "
              f"donor={row.get('asp6_min_donor_distance') or '-'} "
              f"charge={row.get('net_interface_charge')} "
              f"GATE={row.get('gate_pass')} ({time.time() - started:.0f}s)",
              flush=True)

    stats = {"elapsed_seconds": round(time.time() - started, 1),
             "backbones_generated": len(backbone_records),
             "backbones_pose_ok": sum(1 for r in backbone_records if r["ok"]),
             "candidates": len(candidates), "scored": scored,
             "boltz2_calls": calls, "stopped_early": stopped,
             "strategies": STRATEGIES,
             "pose_filter": {"min_pos6_contacts": args.min_pos6_contacts,
                             "max_pos6_distance": args.max_pos6_distance},
             "constraint": POCKET,
             "boltz_kwargs": {"recycling_steps": args.recycling_steps,
                              "sampling_steps": args.sampling_steps,
                              "diffusion_samples": 1,
                              "without_potentials": False},
             "gate": {"donor_residues": sorted(GATE_DONOR_RESIDUES),
                      "donor_distance": GATE_DONOR_DISTANCE,
                      "asp6_buried_frac": "> 0"}}
    (OUT / f"run_stats{args.tag}.json").write_text(json.dumps(stats, indent=2, default=str))
    print(f"[batch-04] {scored} scored, {calls} Boltz-2 calls, "
          f"{stats['elapsed_seconds']:.0f}s {stopped}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
