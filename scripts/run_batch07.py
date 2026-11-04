#!/usr/bin/env python3
"""Batch-07 step 2: a rigidifying redesign around bb54-g2's donor graft.

Not a new pipeline. Every API call, constraint, pacing rule and scoring path is
batch-04/05/06 code, imported and re-used; what is new is *which* backbones and
*which* sequences are built, and the deterministic stiffening edits in
`scripts/rigidify.py`.

Why rigidify, and what is being rigidified
------------------------------------------
Three batches have stalled at 2-of-3 on the Asp6 read-out, and batch-06's
diagnosis is that the contact is **intermittent, not weak** - when `bb54-g2`
reads the carboxylate it reads it at 2.5-3.5 A. Measured here on bb54's own
RFdiffusion backbone (`rigidify.donor_loop_metrics`), the donor anchor is
already on a **helix** (shell secondary structure `HHHH`, 63 % helix overall),
so the mobility that drops the contact is *not* a floppy donor loop: it is the
inter-helix connectors and the termini letting the whole three-helix bundle
rock on the groove. Batch-07 therefore stiffens the scaffold **away from the
paratope**, and keeps the donor shell untouched:

* **Track A - sequence-space rigidification of bb54-g2 itself.** ProteinMPNN is
  re-run on bb54's saved backbone at low temperature (rigid, low-entropy
  sequences), the same deterministic donor graft
  (`run_batch04.donor_anchored_variant`) is applied so the donor geometry that
  works is preserved byte-for-byte, and then `rigidify.rigidifying_edits` adds
  Gly->Ala on loops with negative phi, Pro in the proline phi window, and
  optionally one Cys-Cys staple - never within the 6 nearest Asp6-shell
  positions.
* **Track B - rigid new scaffolds carrying the same graft.** RFdiffusion at
  more diffusion steps (30 vs 15) with the hotspot on C6, at bb54-g2's paratope
  size (46 aa) and two buttressed sizes (52/58 aa) that add scaffold, not
  paratope - batch-06 measured that *shrinking* the binder does not re-centre
  the interface, so nothing here is smaller than bb54-g2. On top of batch-06's
  peptide-centred design-stage filter, a backbone must also be rigid where it
  matters: `--min-helix-frac` of the binder in regular secondary structure,
  the donor anchor on a helix/strand or in a loop no longer than
  `--max-donor-loop`, and no loop longer than `--max-longest-loop`.

Scoring (step 3) is the batch-06 gate, unchanged and stated before it is
applied: the Asp6 direct read-out (binder Arg/Lys/His/Ser/Thr/Tyr heavy atom
within 3.5 A of Asp6 OD1/OD2 **and** `asp6_buried_frac` > 0) in **3 of 3**
independent draws, taken as `diffusion_samples: 3` in a single call, then the
WT arm with `batch06_panel.WT_CRITERION` **verbatim**, then apo and decoy.
`src/design_gates.py` is untouched; no ipTM gap and no MM-GBSA ddG is computed.
"""
from __future__ import annotations

import argparse, json, os, sys, time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any, Dict, List

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(REPO / "scripts"))

import batch06_panel as P  # noqa: E402
import batch07_settings as step1  # noqa: E402
import pose_filter  # noqa: E402
import rigidify  # noqa: E402
import run_batch04 as batch  # noqa: E402
import run_batch06  # noqa: E402  (its peptide-centred design-stage filter)
import run_pipeline  # noqa: E402
import sampling_study as study  # noqa: E402
from nvidia_client import NvidiaApiError, NvidiaClient  # noqa: E402

OUT = REPO / "results" / "batch-07"
B06 = REPO / "results" / "batch-06"
P.OUT = OUT
batch.OUT = OUT
run_batch06.OUT = OUT

# The step-1 winner first, the runner-up second: the best Boltz-2 sampling
# setting is design-specific (sampling study, batch-06), so a design that
# misses 3-of-3 in the first cell is retried in the second and no further.
CONDITION_ORDER = ["r3np", "r3"]

PARENT = {"name": "B6-c6-46_w2_bb54_donor-anchored-g2",
          "sequence": "MMLPDDAIAELSGTPEARELARRLGGDERVAINYRAGRRGRVGVEA",
          "backbone_pdb": B06 / "all_structures" / "B6-c6-46_w2_bb54.pdb"}

# Track B strategies. bb54-g2 is 46 aa; nothing here is smaller.
STRATEGIES: List[Dict[str, Any]] = [
    {"name": "B7-rigid-46", "binder_length": 46, "hotspots": ["C6"],
     "sampling_temp": 0.1,
     "rationale": "bb54-g2's own paratope size, re-diffused at 30 steps and "
                  "filtered for a donor anchor on a regular element"},
    {"name": "B7-rigid-52", "binder_length": 52, "hotspots": ["C6"],
     "sampling_temp": 0.1,
     "rationale": "bb54-g2's paratope with ~6 residues of extra scaffold to "
                  "buttress the donor-carrying helix - more bundle, not more "
                  "interface"},
    {"name": "B7-rigid-58", "binder_length": 58, "hotspots": ["C6"],
     "sampling_temp": 0.15,
     "rationale": "upper buttressing band; the batch-04/05 length that gave "
                  "the project's tightest pose (bb50-g1, 1.58/3.78 A)"},
]


# --- design-stage rigidity filter ------------------------------------------

_PEPTIDE_CENTRED_OK = run_batch06.peptide_centred_design_pose_ok
RIGIDITY = {"min_ss_frac": 0.55, "max_donor_loop": 3, "max_longest_loop": 8}
_CURRENT_PDB: Dict[int, str] = {}


def rigid_design_pose_ok(metrics, **kwargs):
    """batch-06's peptide-centred filter plus the batch-07 rigidity terms."""
    ok, reason = _PEPTIDE_CENTRED_OK(metrics, **kwargs)
    reasons = [reason] if reason else []
    rig = metrics.get("rigidity") or {}
    if rig:
        ss_frac = 1.0 - float(rig.get("ss_loop_frac") or 0.0)
        if ss_frac < RIGIDITY["min_ss_frac"]:
            reasons.append(f"only {ss_frac:.2f} of the binder in regular "
                           f"secondary structure (need >= {RIGIDITY['min_ss_frac']:.2f})")
        if int(rig.get("donor_loop_length") or 0) > RIGIDITY["max_donor_loop"]:
            reasons.append(f"donor anchor sits in a {rig['donor_loop_length']}-residue "
                           f"loop (cap {RIGIDITY['max_donor_loop']})")
        if int(rig.get("longest_loop") or 0) > RIGIDITY["max_longest_loop"]:
            reasons.append(f"longest loop {rig['longest_loop']} residues "
                           f"(cap {RIGIDITY['max_longest_loop']})")
    return (not reasons), "; ".join(r for r in reasons if r)


def backbone_and_grafts(client, strategy, index, target_pdb, args) -> Dict[str, Any]:
    """RFdiffusion backbone -> rigidity metrics -> filter -> MPNN -> grafts."""
    name = f"{strategy['name']}_bb{index:02d}"
    try:
        pdb = client.rfdiffusion_generate(
            input_pdb=target_pdb,
            contigs=run_pipeline.binder_contig(int(strategy["binder_length"])),
            hotspot_res=list(strategy["hotspots"]),
            diffusion_steps=args.diffusion_steps,
            random_seed=args.seed + index * 17)
    except NvidiaApiError as exc:
        print(f"[rfdiffusion] {name} FAILED: {exc}", file=sys.stderr)
        return {"backbone": name, "ok": False, "candidates": []}
    (OUT / "all_structures" / f"{name}.pdb").write_text(pdb)

    try:
        metrics = pose_filter.design_pose_metrics(pdb)
    except Exception as exc:  # noqa: BLE001
        print(f"[pose] {name} unreadable: {exc}", file=sys.stderr)
        return {"backbone": name, "ok": False, "candidates": []}
    shell = batch.pos6_shell_positions(pdb)
    rig = rigidify.donor_loop_metrics(pdb, shell)
    metrics["rigidity"] = rig
    ok, reason = rigid_design_pose_ok(
        metrics, min_pos6_contacts=args.min_pos6_contacts,
        max_pos6_distance=args.max_pos6_distance)
    record = {"backbone": name, "strategy": strategy["name"], "ok": ok,
              "reason": reason, "pose": metrics, "rigidity": rig,
              "pos6_shell": shell[:8], "candidates": []}
    if not ok:
        print(f"[pose] drop {name}: {reason}")
        return record
    try:
        designed = client.proteinmpnn_design(
            input_pdb=pdb, input_pdb_chains=[run_pipeline.RFD_BINDER_CHAIN],
            num_seq_per_target=args.seqs_per_backbone,
            sampling_temp=[float(strategy["sampling_temp"])],
            use_soluble_model=True, random_seed=args.seed + index)
    except NvidiaApiError as exc:
        print(f"[proteinmpnn] {name} FAILED: {exc}", file=sys.stderr)
        return record
    record["candidates"] = build_variants(
        name, strategy, pdb, shell, designed["sequences"], rig, args,
        track="B")
    print(f"[design] {name}: kept, {len(record['candidates'])} variants "
          f"(ss_loop {rig['ss_loop_frac']}, donor_loop {rig['donor_loop_length']}, "
          f"anchor {rig['anchor_ss']})", flush=True)
    return record


def build_variants(backbone_name, strategy, pdb, shell, mpnn_sequences, rig,
                   args, track: str) -> List[Dict[str, Any]]:
    """Donor graft (control) + rigidified variants, per selected MPNN parent."""
    samples = []
    for entry in mpnn_sequences:
        sequence = entry["sequence"].split("/")[0]
        if run_pipeline.sequence_filter_failures(sequence):
            continue
        samples.append({"sequence": sequence, "mpnn_score": entry.get("score"),
                        "donor_bias_score": batch.donor_bias_score(sequence, shell)})
    if not samples:
        return []
    samples.sort(key=lambda s: (-s["donor_bias_score"],
                                s["mpnn_score"] if s["mpnn_score"] is not None else 9))
    out, seen = [], set()
    for index, parent in enumerate(samples[: args.parents_per_backbone], start=1):
        grafted = batch.donor_anchored_variant(parent["sequence"], shell)
        variants = [("graft", grafted, [])]
        rigid_np, notes_np = rigidify.rigidifying_edits(
            grafted, pdb, shell, do_proline=True, do_staple=False)
        variants.append((f"rigid-gp", rigid_np, notes_np))
        rigid_ss, notes_ss = rigidify.rigidifying_edits(
            grafted, pdb, shell, do_proline=True, do_staple=True)
        if rigid_ss != rigid_np:
            variants.append((f"rigid-gps", rigid_ss, notes_ss))
        for label, sequence, notes in variants:
            if sequence in seen or run_pipeline.sequence_filter_failures(sequence):
                continue
            seen.add(sequence)
            out.append({
                "strategy": strategy["name"], "backbone": backbone_name,
                "binder_length": len(sequence),
                "variant": f"{label}-p{index}", "track": track,
                "sequence": sequence, "parent_sequence": parent["sequence"],
                "graft_sequence": grafted,
                "mpnn_score": parent["mpnn_score"],
                "donor_bias_score": batch.donor_bias_score(sequence, shell),
                "rigidifying_edits": notes,
                "rigidity": rig,
                "design_pos6_contacts": "", "design_pos6_distance": ""})
    return out


# --- Track A: rigidify bb54-g2 on its own backbone --------------------------

def track_a_candidates(client, args) -> List[Dict[str, Any]]:
    pdb = PARENT["backbone_pdb"].read_text()
    (OUT / "all_structures" / "bb54_parent.pdb").write_text(pdb)
    shell = batch.pos6_shell_positions(pdb)
    rig = rigidify.donor_loop_metrics(pdb, shell)
    print(f"[track-A] bb54 backbone rigidity: {json.dumps(rig)[:400]}")
    out: List[Dict[str, Any]] = []

    # A0: bb54-g2's own sequence, rigidified - the most direct test.
    for label, staple in (("rigid-gp", False), ("rigid-gps", True)):
        sequence, notes = rigidify.rigidifying_edits(
            PARENT["sequence"], pdb, shell, do_proline=True, do_staple=staple)
        if sequence == PARENT["sequence"] or run_pipeline.sequence_filter_failures(sequence):
            continue
        if any(c["sequence"] == sequence for c in out):
            continue
        out.append({"strategy": "B7-A-bb54", "backbone": "bb54",
                    "binder_length": len(sequence), "variant": f"{label}-self",
                    "track": "A", "sequence": sequence,
                    "parent_sequence": PARENT["sequence"],
                    "graft_sequence": PARENT["sequence"],
                    "mpnn_score": "", "rigidity": rig,
                    "donor_bias_score": batch.donor_bias_score(sequence, shell),
                    "rigidifying_edits": notes,
                    "design_pos6_contacts": "", "design_pos6_distance": ""})

    # A1: fresh low-temperature MPNN sequences on the same backbone, grafted
    #     and rigidified - same donor geometry, stiffer sequence background.
    designed = client.proteinmpnn_design(
        input_pdb=pdb, input_pdb_chains=[run_pipeline.RFD_BINDER_CHAIN],
        num_seq_per_target=args.seqs_per_backbone,
        sampling_temp=[args.track_a_temp], use_soluble_model=True,
        random_seed=args.seed + 7)
    strategy = {"name": "B7-A-bb54"}
    out.extend(build_variants("bb54", strategy, pdb, shell,
                              designed["sequences"], rig, args, track="A"))
    return out


# --- scoring: 3 draws per candidate in one call -----------------------------

def score_candidates(client, candidates, states, args) -> List[Dict[str, Any]]:
    study.CONDITIONS.update(step1.SWEEP)
    summary = []
    calls = 0
    for candidate in candidates:
        if calls >= args.max_boltz_calls:
            print("[stop] Boltz-2 call cap reached"); break
        name = f"{candidate['backbone']}_{candidate['variant']}"
        design = {"name": name, "sequence": candidate["sequence"],
                  "backbone": candidate["backbone"],
                  "variant": candidate["variant"],
                  "binder_length": candidate["binder_length"],
                  "strategy": candidate.get("strategy", ""),
                  "mpnn_score": candidate.get("mpnn_score", ""),
                  "donor_bias_score": candidate.get("donor_bias_score", ""),
                  "design_pos6_contacts": "", "design_pos6_distance": ""}
        record = {"design": name, "track": candidate.get("track"),
                  "sequence": candidate["sequence"],
                  "binder_length": candidate["binder_length"],
                  "rigidifying_edits": candidate.get("rigidifying_edits"),
                  "rigidity": candidate.get("rigidity"),
                  "conditions": {}}
        chosen = None
        print(f"\n[score] === {name} ({candidate['binder_length']} aa, "
              f"track {candidate.get('track')}) ===", flush=True)
        for condition in CONDITION_ORDER[: args.conditions_per_design]:
            if calls >= args.max_boltz_calls:
                break
            already = P.done_keys()
            calls += P.run_draws(client, design, "g12d", condition, states,
                                 args.draws, already)
            rows = P.arm_rows(P.load_rows(), name, "g12d", condition)
            passes = sum(1 for r in rows
                         if str(r.get("asp6_gate_pass")).lower() == "true")
            record["conditions"][condition] = {"gate_passes": passes,
                                               "draws": len(rows)}
            print(f"  -> {condition}: Asp6 gate {passes}/{len(rows)}", flush=True)
            if passes >= args.draws and len(rows) >= args.draws:
                chosen = condition
                break
        record["stabilising_condition"] = chosen
        record["gate_3of3"] = bool(chosen)
        summary.append(record)
        append_summary(record)
    return summary


def append_summary(record: Dict[str, Any]) -> None:
    with (OUT / "screen_summary.jsonl").open("a") as handle:
        handle.write(json.dumps(record, default=str) + "\n")
        handle.flush()
        os.fsync(handle.fileno())


def main(argv=None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--backbones-per-strategy", type=int, default=16)
    parser.add_argument("--seqs-per-backbone", type=int, default=16)
    parser.add_argument("--parents-per-backbone", type=int, default=1)
    parser.add_argument("--diffusion-steps", type=int, default=30)
    parser.add_argument("--design-workers", type=int, default=4)
    parser.add_argument("--seed", type=int, default=7070)
    parser.add_argument("--timeout", type=float, default=900.0)
    parser.add_argument("--min-pos6-contacts", type=int, default=3)
    parser.add_argument("--max-pos6-distance", type=float, default=4.5)
    parser.add_argument("--min-peptide-contact-frac", type=float, default=0.25)
    parser.add_argument("--max-hla-contacts", type=int, default=90)
    parser.add_argument("--track-a-temp", type=float, default=0.1)
    parser.add_argument("--draws", type=int, default=3)
    parser.add_argument("--conditions-per-design", type=int, default=2)
    parser.add_argument("--max-candidates", type=int, default=24)
    parser.add_argument("--max-boltz-calls", type=int, default=40)
    parser.add_argument("--design-only", action="store_true")
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args(argv)

    run_batch06.PEPTIDE_CONTACT_FRAC_MIN = args.min_peptide_contact_frac
    run_batch06.MAX_HLA_CONTACTS = args.max_hla_contacts

    for sub in ("all_structures", "structures"):
        (OUT / sub).mkdir(parents=True, exist_ok=True)
    started = time.time()
    target_dir = REPO / "data" / "targets"
    states = run_pipeline.load_target_states(target_dir)
    target_pdb = run_pipeline.trimmed_design_target(target_dir / "g12d_complex.pdb")
    client = NvidiaClient(timeout=args.timeout)

    candidates_path = OUT / "candidates.json"
    poses_path = OUT / "backbone_poses.json"
    if args.resume and candidates_path.exists():
        candidates = json.loads(candidates_path.read_text())
        backbone_records = json.loads(poses_path.read_text()) if poses_path.exists() else []
        print(f"[design] resuming: {len(candidates)} candidates")
    else:
        candidates = track_a_candidates(client, args)
        print(f"[track-A] {len(candidates)} candidates")
        backbone_records: List[Dict[str, Any]] = []
        jobs = [(s, i) for s in STRATEGIES for i in range(args.backbones_per_strategy)]
        with ThreadPoolExecutor(max_workers=args.design_workers) as pool:
            futures = [pool.submit(backbone_and_grafts,
                                   NvidiaClient(timeout=args.timeout),
                                   s, i, target_pdb, args) for s, i in jobs]
            for future in as_completed(futures):
                record = future.result()
                backbone_records.append(record)
                candidates.extend(record.get("candidates") or [])
                print(f"[design] backbones={len(backbone_records)}/{len(jobs)} "
                      f"kept={sum(1 for r in backbone_records if r.get('ok'))} "
                      f"candidates={len(candidates)} "
                      f"({time.time() - started:.0f}s)", flush=True)
        seen, unique = set(), []
        for candidate in candidates:
            if candidate["sequence"] in seen:
                continue
            seen.add(candidate["sequence"])
            unique.append(candidate)
        candidates = unique
        candidates_path.write_text(json.dumps(candidates, indent=2, default=str))
        poses_path.write_text(json.dumps(backbone_records, indent=2, default=str))

    print(f"[design] {len(candidates)} unique candidates "
          f"({time.time() - started:.0f}s)")
    if args.design_only:
        return 0

    # Track A first (it is the direct test of the batch-07 hypothesis), then
    # Track B ordered by donor content, as batches 04-06 ordered their queues.
    candidates.sort(key=lambda c: (0 if c.get("track") == "A" else 1,
                                   -float(c.get("donor_bias_score") or 0)))
    summary = score_candidates(client, candidates[: args.max_candidates],
                               states, args)
    stats = {"elapsed_seconds": round(time.time() - started, 1),
             "candidates": len(candidates),
             "scored": len(summary),
             "gate_3of3": [r["design"] for r in summary if r["gate_3of3"]],
             "condition_order": CONDITION_ORDER,
             "rigidity_filter": RIGIDITY,
             "strategies": STRATEGIES,
             "constraint": batch.POCKET,
             "draws_per_call": args.draws}
    (OUT / "screen_stats.json").write_text(json.dumps(stats, indent=2, default=str))
    print(f"\n[batch-07] {len(summary)} designs scored, "
          f"{len(stats['gate_3of3'])} clear the Asp6 gate 3 of 3")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
