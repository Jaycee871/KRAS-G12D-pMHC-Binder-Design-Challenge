#!/usr/bin/env python3
"""Batch-05: scale-up of strategy B4-c6-62 with a donor graft for EVERY backbone.

This is deliberately **not** a new pipeline. Every stage is
`scripts/run_batch04.py`'s code, imported and re-used; this driver changes only
the three things batch-04's SUMMARY recommended:

1. **One strategy, at volume.** `B4-c6-62` (62 aa binder, hotspot `C6` alone,
   ProteinMPNN `sampling_temp` 0.15) is the only strategy in batch-04 that
   produced a gate survivor whose pose replicated, so the whole budget goes to
   it. Here it is re-named `B5-c6-62`; the parameters are identical.
2. **A deterministic donor-anchored graft for every backbone** (and for every
   MPNN parent selected on that backbone), rather than batch-04's single graft
   off the single best sample, which some backbones never got at all. MPNN is
   kept only where batch-04 already used it - as the source of the parent
   sequence the graft is applied to. MPNN-selected sequences are *not* scored
   on their own: 0 of 38 passed the gate in batch-04, and the deterministic
   graft is what produced both passes.
3. **Unique artefact names.** Batch-04 named scored structures
   `{backbone}_{variant}`, so two candidates from one backbone overwrote each
   other. Variant labels here are unique per candidate
   (`donor-anchored-g1`, `-g2`, ...), which keeps the batch-04 replicate and
   RMSD scripts working unchanged while removing the collision.

Scoring is unchanged: one constrained Boltz-2 call per candidate, strictly
serial at the calibrated >= 8 s interval, with the validated flat pocket term
`{"binder":"D","contacts":[{"id":"C","residue_index":6}],"max_distance":6.0}`.
Replication is folded into the gate downstream (scripts/replicate_batch05.py):
a design passes only if it clears the Asp6 donor gate in 3 of 3 independent
draws. Nothing is ranked on ipTM, `specificity_gap` or `ddG_sel`.
"""
from __future__ import annotations

import argparse, json, sys, time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any, Dict, List

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(REPO / "scripts"))

import run_batch04 as batch  # noqa: E402
import run_pipeline  # noqa: E402
from nvidia_client import NvidiaApiError, NvidiaClient  # noqa: E402

OUT = REPO / "results" / "batch-05"
batch.OUT = OUT  # every batch-04 helper writes into results/batch-05 instead

STRATEGY: Dict[str, Any] = {
    "name": "B5-c6-62", "binder_length": 62, "hotspots": ["C6"],
    "sampling_temp": 0.15,
    "rationale": "batch-04's B4-c6-62 verbatim - the only strategy that "
                 "produced a gate survivor whose pose replicated - run at "
                 "volume with a deterministic donor graft per backbone",
}


def grafts_for_backbone(record: Dict[str, Any]) -> List[Dict[str, Any]]:
    """A deterministic donor-anchored graft for *every* selected MPNN parent.

    The pos6 shell is recomputed from the saved RFdiffusion backbone (the
    record only keeps its first 8 entries) so the graft is byte-identical to
    what batch-04's `donor_anchored_variant` would produce.
    """
    name = record["backbone"]
    pdb_path = OUT / "all_structures" / f"{name}.pdb"
    if not pdb_path.exists():
        return []
    shell = batch.pos6_shell_positions(pdb_path.read_text())
    out: List[Dict[str, Any]] = []
    seen = set()
    parents = [c for c in record.get("candidates") or []
               if c.get("variant") == "mpnn-selected"]
    for index, parent in enumerate(parents, start=1):
        grafted = batch.donor_anchored_variant(parent["sequence"], shell)
        if grafted in seen or run_pipeline.sequence_filter_failures(grafted):
            continue
        seen.add(grafted)
        out.append({**parent,
                    "sequence": grafted,
                    "parent_sequence": parent["sequence"],
                    "donor_bias_score": batch.donor_bias_score(grafted, shell),
                    "variant": f"donor-anchored-g{index}"})
    return out


def main(argv=None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--backbones", type=int, default=64)
    parser.add_argument("--seqs-per-backbone", type=int, default=16)
    parser.add_argument("--seqs-scored-per-backbone", type=int, default=2)
    parser.add_argument("--diffusion-steps", type=int, default=15)
    parser.add_argument("--recycling-steps", type=int, default=3)
    parser.add_argument("--sampling-steps", type=int, default=50)
    parser.add_argument("--design-workers", type=int, default=4)
    parser.add_argument("--seed", type=int, default=5050)
    parser.add_argument("--timeout", type=float, default=600.0)
    parser.add_argument("--min-pos6-contacts", type=int, default=2)
    parser.add_argument("--max-pos6-distance", type=float, default=5.0)
    parser.add_argument("--max-boltz-calls", type=int, default=220)
    parser.add_argument("--score-time-budget", type=float, default=7200.0)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--tag", default="")
    args = parser.parse_args(argv)

    started = time.time()
    for sub in ("all_structures", "scored_structures", "replicate_structures"):
        (OUT / sub).mkdir(parents=True, exist_ok=True)

    target_dir = REPO / "data" / "targets"
    target_states = run_pipeline.load_target_states(target_dir)
    target_pdb = run_pipeline.trimmed_design_target(target_dir / "g12d_complex.pdb")

    candidates_path = OUT / f"candidates{args.tag}.json"
    backbone_records: List[Dict[str, Any]] = []

    if args.resume and candidates_path.exists():
        candidates = json.loads(candidates_path.read_text())
        poses = OUT / f"backbone_poses{args.tag}.json"
        backbone_records = json.loads(poses.read_text()) if poses.exists() else []
        print(f"[design] resuming: {len(candidates)} candidates from "
              f"{candidates_path}")
    else:
        candidates = []
        with ThreadPoolExecutor(max_workers=args.design_workers) as pool:
            futures = [pool.submit(batch.backbone_and_sequences,
                                   NvidiaClient(timeout=args.timeout),
                                   STRATEGY, i, target_pdb, args)
                       for i in range(args.backbones)]
            for future in as_completed(futures):
                record = future.result()
                backbone_records.append(record)
                record["grafts"] = grafts_for_backbone(record)
                candidates.extend(record["grafts"])
                print(f"[design] backbones={len(backbone_records)}/{args.backbones} "
                      f"on-pos6={sum(1 for r in backbone_records if r['ok'])} "
                      f"grafts={len(candidates)} "
                      f"({time.time() - started:.0f}s)", flush=True)

        seen, unique = set(), []
        for candidate in candidates:
            if candidate["sequence"] in seen:
                continue
            seen.add(candidate["sequence"])
            unique.append(candidate)
        candidates = unique
        candidates_path.write_text(json.dumps(candidates, indent=2))
        (OUT / f"backbone_poses{args.tag}.json").write_text(json.dumps(
            backbone_records, indent=2, default=str))

    pose_ok = sum(1 for r in backbone_records if r.get("ok"))
    with_graft = len({c["backbone"] for c in candidates})
    print(f"[design] stage 1-2 done in {time.time() - started:.0f}s: "
          f"{len(backbone_records)} backbones, {pose_ok} reach position 6, "
          f"{with_graft} backbones with >= 1 graft, {len(candidates)} grafts")

    return score_stage(candidates, backbone_records, target_states, args, started)


def score_stage(candidates, backbone_records, target_states, args, started) -> int:
    client = NvidiaClient(timeout=args.timeout)
    checkpoint = OUT / "checkpoint.jsonl"
    done = set()
    if checkpoint.exists():
        done = {json.loads(l)["row"]["sequence"]
                for l in checkpoint.read_text().splitlines() if l.strip()}
        if done:
            print(f"[score] checkpoint holds {len(done)} scored sequences")
    candidates = sorted(candidates,
                        key=lambda c: -float(c.get("donor_bias_score") or 0))
    todo = [c for c in candidates if c["sequence"] not in done]
    calls = scored = 0
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
            result = batch.score_candidate(client, candidate, target_states, args)
        except (NvidiaApiError, Exception) as exc:  # noqa: BLE001
            print(f"[boltz2] {candidate['backbone']} dropped: "
                  f"{type(exc).__name__}: {exc}", file=sys.stderr)
            calls += 1
            continue
        calls += 1
        scored += 1
        row = result["row"]
        row["parent_sequence"] = candidate.get("parent_sequence", "")
        name = f"{candidate['backbone']}_{candidate['variant']}"
        (OUT / "scored_structures" / f"{name}_g12d.cif").write_text(result["structure"])
        batch.append_checkpoint({"row": row,
                                 "structure_file": f"scored_structures/{name}_g12d.cif",
                                 "scored_at": time.strftime("%Y-%m-%dT%H:%M:%S%z")})
        print(f"[boltz2] {scored}/{len(todo)} {name} iptm={row.get('iptm_g12d')} "
              f"pep_buried={row.get('peptide_buried_frac')} "
              f"pos6_sc={row.get('pos6_sidechain_contacts')} "
              f"donor={row.get('asp6_min_donor_distance') or '-'} "
              f"charge={row.get('net_interface_charge')} "
              f"GATE={row.get('gate_pass')} ({time.time() - started:.0f}s)",
              flush=True)

    stats = {"elapsed_seconds": round(time.time() - started, 1),
             "backbones_generated": len(backbone_records),
             "backbones_pose_ok": sum(1 for r in backbone_records if r.get("ok")),
             "grafts": len(candidates), "scored_this_segment": scored,
             "boltz2_calls_this_segment": calls, "stopped_early": stopped,
             "strategy": STRATEGY,
             "pose_filter": {"min_pos6_contacts": args.min_pos6_contacts,
                             "max_pos6_distance": args.max_pos6_distance},
             "constraint": batch.POCKET,
             "boltz_kwargs": {"recycling_steps": args.recycling_steps,
                              "sampling_steps": args.sampling_steps,
                              "diffusion_samples": 1,
                              "without_potentials": False},
             "gate": {"donor_residues": sorted(batch.GATE_DONOR_RESIDUES),
                      "donor_distance": batch.GATE_DONOR_DISTANCE,
                      "asp6_buried_frac": "> 0",
                      "replication": "3 of 3 independent constrained draws"}}
    suffix = args.tag or ("_resumed" if (OUT / "run_stats.json").exists() else "")
    (OUT / f"run_stats{suffix}.json").write_text(json.dumps(stats, indent=2, default=str))
    print(f"[batch-05] {scored} scored this segment, {calls} Boltz-2 calls, "
          f"{stats['elapsed_seconds']:.0f}s {stopped}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
