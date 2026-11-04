#!/usr/bin/env python3
"""Batch-06: a much smaller, peptide-centred paratope, with the WT arm in the gate.

This is **not** a new pipeline. Stages 1-3 are `scripts/run_batch04.py` and
`scripts/run_batch05.py` code, imported and re-used; this driver changes only
what `results/selectivity-panel/SUMMARY.md` said had to change:

1. **Much shorter binders.** 34 / 40 / 46 aa instead of the 62-72 aa used in
   batches 04/05, whose interfaces came out at 600-1100 A^2 dominated by HLA
   groove contact. A binder that is physically too small to drape onto the
   alpha1/alpha2 helices has to pay for its interface with the peptide.
2. **A peptide-centred design-stage pose filter.** On top of batch-05's filter
   (>= 10 binder heavy atoms on the peptide, >= 3 peptide residues, position-6
   contacts) the backbone must put at least `--min-peptide-contact-frac` of its
   target contact on the peptide itself, and must come within
   `--max-pos6-distance` of position 6 with at least `--min-pos6-contacts`
   atoms. An HLA-draped backbone costs zero ProteinMPNN and zero Boltz-2 calls.
3. **The optimised quantity is the ratio, not the area.**
   `buried_area_peptide` vs `buried_area_hla` is recorded for every scored
   complex and is what survivors are ranked on; `peptide_area_share` =
   buried_area_peptide / buried_area_total.
4. **A stabilising sampling setting from the start.** Screening runs at
   `recycling_steps` 8, the setting that collapsed the pose spread for two of
   three designs in the sampling study. The per-design choice between `recyc8`
   and `combo` is made in stage 4 (`scripts/batch06_panel.py`).

Nothing is ranked on ipTM, `specificity_gap` or `ddG_sel`; all three are
retired. `src/design_gates.py` is untouched - its thresholds are reported, not
applied as a filter here.
"""
from __future__ import annotations

import argparse, json, sys, time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any, Dict, List, Tuple

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(REPO / "scripts"))

import pose_filter  # noqa: E402
import run_batch04 as batch  # noqa: E402
import run_batch05 as b5  # noqa: E402
import run_pipeline  # noqa: E402
from nvidia_client import NvidiaApiError, NvidiaClient  # noqa: E402

OUT = REPO / "results" / "batch-06"
batch.OUT = OUT     # every batch-04 helper writes into results/batch-06
b5.OUT = OUT        # ... and so does batch-05's graft helper

STRATEGIES: List[Dict[str, Any]] = [
    {"name": "B6-c6-34", "binder_length": 34, "hotspots": ["C6"],
     "sampling_temp": 0.15,
     "rationale": "roughly half the batch-05 length - too small to reach the "
                  "alpha1/alpha2 rim and the peptide at the same time"},
    {"name": "B6-c6-40", "binder_length": 40, "hotspots": ["C6"],
     "sampling_temp": 0.15,
     "rationale": "smallest length that reliably folds as a mini-binder while "
                  "still being peptide-sized"},
    {"name": "B6-c6-46", "binder_length": 46, "hotspots": ["C6"],
     "sampling_temp": 0.2,
     "rationale": "upper end of the shrunken band, higher MPNN temperature "
                  "for topological diversity among survivors"},
]

CSV_COLUMNS = batch.CSV_COLUMNS + [
    "peptide_area_share", "buried_area_ratio_pep_hla", "design_peptide_contact_frac",
]


# --- stage 1: peptide-centred design-stage pose filter ----------------------

_ORIGINAL_DESIGN_POSE_OK = pose_filter.design_pose_ok
PEPTIDE_CONTACT_FRAC_MIN = 0.25     # overridden from argv in main()
MAX_HLA_CONTACTS = 90               # cap on groove contact of the backbone


def peptide_centred_design_pose_ok(metrics, **kwargs) -> Tuple[bool, str]:
    """batch-05's filter plus an explicit cap on groove contact.

    The extra terms are the whole point of batch-06: a backbone that reaches
    position 6 but spends most of its contact on the HLA helices produced
    exactly the 600-1100 A^2 groove-dominated interfaces the off-target panel
    found to be non-selective, so it is dropped before it costs a sequence.
    """
    ok, reason = _ORIGINAL_DESIGN_POSE_OK(metrics, **kwargs)
    reasons = [reason] if reason else []
    frac = float(metrics.get("peptide_contact_frac") or 0.0)
    if frac < PEPTIDE_CONTACT_FRAC_MIN:
        reasons.append(
            f"only {frac:.2f} of target contact is on the peptide "
            f"(need >= {PEPTIDE_CONTACT_FRAC_MIN:.2f})")
    hla = float(metrics.get("hla_contacts") or 0.0)
    if hla > MAX_HLA_CONTACTS:
        reasons.append(
            f"{hla:.0f} binder heavy atoms on the HLA groove "
            f"(cap {MAX_HLA_CONTACTS})")
    return (not reasons), "; ".join(r for r in reasons if r)


pose_filter.design_pose_ok = peptide_centred_design_pose_ok


# --- stage 3: screening call, at the stabilising setting --------------------

def patched_boltz2(client, polymers, args):
    """batch.boltz2_constrained with the sampling fields taken from `args`.

    Identical body otherwise - same polymers, same step_scale, same single
    diffusion sample, same validated flat POCKET constraint
    {"binder":"D","contacts":[{"id":"C","residue_index":6}],"max_distance":6.0}.
    """
    from nvidia_client import BOLTZ2_PATH
    body = {"polymers": [dict(p) for p in polymers],
            "recycling_steps": int(args.recycling_steps),
            "sampling_steps": int(args.sampling_steps),
            "diffusion_samples": 1, "step_scale": 1.638,
            "without_potentials": bool(getattr(args, "without_potentials", False)),
            "constraints": batch.POCKET}
    return client.post(BOLTZ2_PATH, body)


batch.boltz2_constrained = patched_boltz2


def area_columns(row: Dict[str, Any]) -> None:
    """buried_area_peptide vs buried_area_hla - the batch-06 objective."""
    total = row.get("buried_area_total")
    pep = row.get("buried_area_peptide")
    hla = row.get("buried_area_hla")
    try:
        total, pep = float(total), float(pep)
        row["peptide_area_share"] = round(pep / total, 3) if total > 0 else None
    except (TypeError, ValueError):
        row["peptide_area_share"] = None
    try:
        hla = float(hla)
        row["buried_area_ratio_pep_hla"] = round(float(pep) / hla, 3) if hla > 0 else None
    except (TypeError, ValueError):
        row["buried_area_ratio_pep_hla"] = None


def main(argv=None) -> int:
    global PEPTIDE_CONTACT_FRAC_MIN, MAX_HLA_CONTACTS
    parser = argparse.ArgumentParser()
    parser.add_argument("--backbones-per-strategy", type=int, default=32)
    parser.add_argument("--seqs-per-backbone", type=int, default=16)
    parser.add_argument("--seqs-scored-per-backbone", type=int, default=2)
    parser.add_argument("--diffusion-steps", type=int, default=15)
    parser.add_argument("--recycling-steps", type=int, default=8)
    parser.add_argument("--sampling-steps", type=int, default=50)
    parser.add_argument("--without-potentials", action="store_true")
    parser.add_argument("--design-workers", type=int, default=4)
    parser.add_argument("--seed", type=int, default=6060)
    parser.add_argument("--timeout", type=float, default=600.0)
    parser.add_argument("--min-pos6-contacts", type=int, default=3)
    parser.add_argument("--max-pos6-distance", type=float, default=4.5)
    parser.add_argument("--min-peptide-contact-frac", type=float, default=0.25)
    parser.add_argument("--max-hla-contacts", type=int, default=90)
    parser.add_argument("--max-boltz-calls", type=int, default=60)
    parser.add_argument("--score-time-budget", type=float, default=5400.0)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--tag", default="")
    args = parser.parse_args(argv)

    PEPTIDE_CONTACT_FRAC_MIN = args.min_peptide_contact_frac
    MAX_HLA_CONTACTS = args.max_hla_contacts

    started = time.time()
    for sub in ("all_structures", "scored_structures", "structures"):
        (OUT / sub).mkdir(parents=True, exist_ok=True)

    target_dir = REPO / "data" / "targets"
    target_states = run_pipeline.load_target_states(target_dir)
    target_pdb = run_pipeline.trimmed_design_target(target_dir / "g12d_complex.pdb")

    candidates_path = OUT / f"candidates{args.tag}.json"
    poses_path = OUT / f"backbone_poses{args.tag}.json"
    backbone_records: List[Dict[str, Any]] = []

    if args.resume and candidates_path.exists():
        candidates = json.loads(candidates_path.read_text())
        backbone_records = json.loads(poses_path.read_text()) if poses_path.exists() else []
        print(f"[design] resuming: {len(candidates)} candidates from {candidates_path}")
    else:
        candidates = []
        jobs = [(s, i) for s in STRATEGIES for i in range(args.backbones_per_strategy)]
        with ThreadPoolExecutor(max_workers=args.design_workers) as pool:
            futures = [pool.submit(batch.backbone_and_sequences,
                                   NvidiaClient(timeout=args.timeout),
                                   s, i, target_pdb, args) for s, i in jobs]
            for future in as_completed(futures):
                record = future.result()
                backbone_records.append(record)
                record["grafts"] = b5.grafts_for_backbone(record)
                candidates.extend(record["grafts"])
                print(f"[design] backbones={len(backbone_records)}/{len(jobs)} "
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
        poses_path.write_text(json.dumps(backbone_records, indent=2, default=str))

    pose_ok = sum(1 for r in backbone_records if r.get("ok"))
    print(f"[design] stage 1-2 done in {time.time() - started:.0f}s: "
          f"{len(backbone_records)} backbones, {pose_ok} reach position 6 and "
          f"stay peptide-centred, {len(candidates)} grafts")

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
    frac_by_bb = {r["backbone"]: (r.get("pose") or {}).get("peptide_contact_frac")
                  for r in backbone_records if r.get("backbone")}
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
        frac = frac_by_bb.get(candidate["backbone"])
        row["design_peptide_contact_frac"] = round(float(frac), 3) if frac is not None else ""
        area_columns(row)
        # unique per candidate: {backbone}_{variant} collided in batch-04 when a
        # backbone had two variants; variant labels are unique per candidate here.
        name = f"{candidate['backbone']}_{candidate['variant']}"
        (OUT / "scored_structures" / f"{name}_g12d.cif").write_text(result["structure"])
        batch.append_checkpoint({"row": row,
                                 "structure_file": f"scored_structures/{name}_g12d.cif",
                                 "boltz_kwargs": {"recycling_steps": args.recycling_steps,
                                                  "sampling_steps": args.sampling_steps,
                                                  "diffusion_samples": 1,
                                                  "without_potentials": bool(args.without_potentials)},
                                 "scored_at": time.strftime("%Y-%m-%dT%H:%M:%S%z")})
        print(f"[boltz2 {calls}] {scored}/{len(todo)} {name} "
              f"area={row.get('buried_area_total')} "
              f"pep={row.get('buried_area_peptide')} hla={row.get('buried_area_hla')} "
              f"share={row.get('peptide_area_share')} "
              f"pep_bf={row.get('peptide_buried_frac')} "
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
             "strategies": STRATEGIES,
             "design_pose_filter": {
                 "min_pos6_contacts": args.min_pos6_contacts,
                 "max_pos6_distance": args.max_pos6_distance,
                 "min_peptide_contact_frac": args.min_peptide_contact_frac,
                 "max_hla_contacts": args.max_hla_contacts},
             "constraint": batch.POCKET,
             "boltz_kwargs": {"recycling_steps": args.recycling_steps,
                              "sampling_steps": args.sampling_steps,
                              "diffusion_samples": 1,
                              "without_potentials": bool(args.without_potentials)},
             "gate": {"donor_residues": sorted(batch.GATE_DONOR_RESIDUES),
                      "donor_distance": batch.GATE_DONOR_DISTANCE,
                      "asp6_buried_frac": "> 0",
                      "replication": "3 of 3 independent constrained draws "
                                     "(scripts/batch06_panel.py)",
                      "wt_arm": "inside the gate (scripts/batch06_panel.py)"}}
    suffix = args.tag or ("_resumed" if (OUT / "run_stats.json").exists() else "")
    (OUT / f"run_stats{suffix}.json").write_text(json.dumps(stats, indent=2, default=str))
    print(f"[batch-06] {scored} scored this segment, {calls} Boltz-2 calls, "
          f"{stats['elapsed_seconds']:.0f}s {stopped}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
