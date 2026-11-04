#!/usr/bin/env python3
"""Batch-08 stage 1-2: steer the paratope by SHAPE, cutting alpha1/alpha2 contact.

This is **not** a new pipeline.  Stages 1-2 are `scripts/run_batch04.py` +
`scripts/run_batch05.py` code, imported and re-used exactly as batch-06 and
batch-07 re-used them; stage 3 onwards is `scripts/batch06_panel.py` driven by
`scripts/batch08_panel.py`.  Only the *selection* changes, and only in the two
places where interface shape is decided:

1. **Shape-steered backbone generation and selection.**  Batch-06 proved that
   shrinking a binder scales the whole interface down (peptide share of buried
   area 0.209 -> 0.213, nil) and batch-07 proved stiffening it changes nothing;
   the groove is still ~75 % of bb54-g2's interface.  So the lever here is
   topology: hotspot sets that pull the paratope *along the peptide axis*
   rather than across the groove, and a design-stage filter that measures, on
   the RFdiffusion backbone itself, how much of the contact lands on the
   alpha1 (HLA 57-84) and alpha2 (138-180) helices and from what direction the
   binder comes down (`src/groove_geometry.py`).  Binder length is held at
   bb54-g2's size or a little above (46-54 aa): never smaller - that is a
   measured dead end - and never rigidified.

2. **Explicit selection against alpha-helix contact in sequence design.**
   ProteinMPNN is over-sampled as before and the donor-biased selection of
   batch-04 is kept, but the selection score now carries a penalty for large
   hydrophobics at *helix-facing* backbone positions, and every donor graft
   gets a paired **de-wetted twin** (`-dw`): at helix-facing positions only,
   Phe/Trp/Tyr/Leu/Ile/Met/Val -> Ser, with the six positions nearest Asp6
   protected so the donor chemistry that works is untouched by construction.
   Graft and twin sit on the *same backbone*, so the effect of de-wetting the
   helix face is measurable rather than confounded.

Thresholds of the design-stage shape filter, fixed here before the batch ran
(they are design knobs, not result gates - nothing in `src/design_gates.py` is
touched and no scoring gate is loosened anywhere):

    peptide_contact_frac   >= 0.30     NEW, tighter than batch-06's 0.25
    helix_contact_frac     <= 0.30     NEW - the actual lever
    approach_height        >=  2.0 A   NEW - binder centroid at or above the
                                       peptide, i.e. coming down on the groove
                                       rather than hanging off its side
    abs_lateral_offset     <=  5.0 A   NEW - not displaced onto a helix wall
    over_peptide_frac      >= 0.55     NEW - CAs within 7 A of the groove axis

These numbers were calibrated on the 97 backbones batch-07 had already
generated (prior data, before a single batch-08 structure was scored): of the
14 that cleared the base pose filter there, the median helix contact fraction
was 0.23 and the median approach height was -6.7 A - most RFdiffusion
backbones hang off the side of the groove rather than coming down on it - and
exactly one of the 97 would have cleared the five thresholds above.  A ~7 %
yield is the design target; the reserve list below covers the downside.
  plus batch-06's base filter: >= 10 peptide contacts, >= 3 peptide residues,
  >= 3 position-6 contacts, closest approach to position 6 <= 4.5 A.

Backbones that clear the base filter but miss a shape threshold are kept on a
ranked reserve list (lowest `helix_contact_frac` first), so a batch is never
left with nothing to score; which list a design came from is recorded.

Retired metrics appear nowhere: no ipTM `specificity_gap`, no paired MM-GBSA
ddG.  ipTM is recorded downstream as a pose-quality read-out only.
"""
from __future__ import annotations

import argparse, json, sys, time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any, Dict, List, Tuple

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(REPO / "scripts"))

import complex_analysis  # noqa: E402
import groove_geometry as gg  # noqa: E402
import pose_filter  # noqa: E402
import run_batch04 as batch  # noqa: E402
import run_batch05 as b5  # noqa: E402
import run_pipeline  # noqa: E402
from nvidia_client import NvidiaApiError, NvidiaClient  # noqa: E402

OUT = REPO / "results" / "batch-08"
batch.OUT = OUT
b5.OUT = OUT

# --- the shape thresholds, stated before the run ---------------------------
SHAPE = {
    "min_peptide_contact_frac": 0.30,
    "max_helix_contact_frac": 0.30,
    "min_approach_height": 2.0,
    "max_abs_lateral_offset": 5.0,
    "min_over_peptide_frac": 0.55,
}

STRATEGIES: List[Dict[str, Any]] = [
    {"name": "B8-apex-46", "binder_length": 46, "hotspots": ["C6"],
     "sampling_temp": 0.2,
     "rationale": "bb54-g2's size and hotspot, but kept only when the backbone "
                  "comes down on position 6 from above the groove axis"},
    {"name": "B8-span-50", "binder_length": 50, "hotspots": ["C4", "C6", "C8"],
     "sampling_temp": 0.2,
     "rationale": "hotspots spread along the peptide long axis, so the "
                  "paratope has to lie along the peptide rather than across "
                  "the groove and onto the alpha1/alpha2 walls"},
    {"name": "B8-crest-46", "binder_length": 46, "hotspots": ["C5", "C6", "C7"],
     "sampling_temp": 0.2,
     "rationale": "the solvent-exposed crest of the bulge only - the three "
                  "residues whose side chains point out of the groove"},
    {"name": "B8-apex-54", "binder_length": 54, "hotspots": ["C6"],
     "sampling_temp": 0.25,
     "rationale": "slightly larger than bb54-g2 (never smaller: shrinking is a "
                  "measured dead end) at a higher MPNN temperature, for "
                  "topological diversity among shape survivors"},
]

HELIX_APOLAR = set("FWYLIMV")
DEWET_TO = "S"
PROTECTED_SHELL = 6          # positions nearest Asp6 are never edited


# --- stage 1: the shape filter ---------------------------------------------

_BASE_DESIGN_POSE_OK = pose_filter.design_pose_ok


def shape_reasons(shape: Dict[str, Any]) -> List[str]:
    out: List[str] = []
    if shape["peptide_contact_frac"] < SHAPE["min_peptide_contact_frac"]:
        out.append(f"peptide contact fraction {shape['peptide_contact_frac']:.2f} "
                   f"< {SHAPE['min_peptide_contact_frac']}")
    if shape["helix_contact_frac"] > SHAPE["max_helix_contact_frac"]:
        out.append(f"helix contact fraction {shape['helix_contact_frac']:.2f} "
                   f"> {SHAPE['max_helix_contact_frac']} (alpha1 "
                   f"{shape['alpha1_contacts']}, alpha2 {shape['alpha2_contacts']})")
    if shape["approach_height"] < SHAPE["min_approach_height"]:
        out.append(f"binder centroid only {shape['approach_height']:.1f} A above "
                   f"the peptide (need >= {SHAPE['min_approach_height']})")
    if shape["abs_lateral_offset"] > SHAPE["max_abs_lateral_offset"]:
        out.append(f"lateral offset {shape['abs_lateral_offset']:.1f} A towards a "
                   f"helix wall (max {SHAPE['max_abs_lateral_offset']})")
    if shape["over_peptide_frac"] < SHAPE["min_over_peptide_frac"]:
        out.append(f"only {shape['over_peptide_frac']:.2f} of binder CAs sit over "
                   f"the groove axis (need >= {SHAPE['min_over_peptide_frac']})")
    return out


# --- stage 2: de-wetting the helix face ------------------------------------

def dewetted_variant(sequence: str, helix_positions: List[Dict[str, Any]],
                     shell: List[Dict[str, Any]]) -> Tuple[str, List[int]]:
    """Large hydrophobics at helix-facing positions -> Ser, Asp6 shell protected."""
    protected = {e["position"] for e in shell[:PROTECTED_SHELL]}
    residues = list(sequence)
    edited: List[int] = []
    for entry in helix_positions:
        position = entry["position"]
        if position in protected or position > len(residues):
            continue
        if residues[position - 1] in HELIX_APOLAR:
            residues[position - 1] = DEWET_TO
            edited.append(position)
    return "".join(residues), edited


def helix_apolar_penalty(sequence: str, helix_positions: List[Dict[str, Any]],
                         shell: List[Dict[str, Any]]) -> float:
    """How much greasy surface this sequence points at the alpha1/alpha2 walls."""
    protected = {e["position"] for e in shell[:PROTECTED_SHELL]}
    score = 0.0
    for entry in helix_positions:
        position = entry["position"]
        if position in protected or position > len(sequence):
            continue
        if sequence[position - 1] in HELIX_APOLAR:
            score += max(0.0, 1.0 - (entry["helix_distance"] - 3.0) / 6.0)
    return round(score, 3)


# --- stage 1 + 2 driver (batch-04's, with the shape read-out added) --------

def backbone_job(client: NvidiaClient, strategy, index, target_pdb, args) -> Dict[str, Any]:
    name = f"{strategy['name']}{args.tag}_bb{index:02d}"
    try:
        pdb = client.rfdiffusion_generate(
            input_pdb=target_pdb,
            contigs=run_pipeline.binder_contig(int(strategy["binder_length"])),
            hotspot_res=list(strategy["hotspots"]),
            diffusion_steps=args.diffusion_steps,
            random_seed=args.seed + index * 17)
    except NvidiaApiError as exc:
        print(f"[rfdiffusion] {name} FAILED: {exc}", file=sys.stderr)
        return {"backbone": name, "strategy": strategy["name"], "ok": False,
                "reason": f"rfdiffusion: {exc}", "candidates": []}

    (OUT / "all_structures" / f"{name}.pdb").write_text(pdb)
    try:
        metrics = pose_filter.design_pose_metrics(pdb)
        shape = gg.design_shape_metrics(pdb)
    except complex_analysis.StructureAnalysisError as exc:
        print(f"[pose] {name} unreadable: {exc}", file=sys.stderr)
        return {"backbone": name, "strategy": strategy["name"], "ok": False,
                "reason": f"unreadable: {exc}", "candidates": []}

    base_ok, base_reason = _BASE_DESIGN_POSE_OK(
        metrics, min_pos6_contacts=args.min_pos6_contacts,
        max_pos6_distance=args.max_pos6_distance)
    reasons = shape_reasons(shape)
    record = {"backbone": name, "strategy": strategy["name"],
              "binder_length": strategy["binder_length"],
              "base_ok": base_ok, "base_reason": base_reason,
              "shape": shape, "shape_reasons": reasons,
              "ok": bool(base_ok and not reasons),
              "pose": metrics, "candidates": []}
    if not base_ok:
        print(f"[pose] drop {name}: {base_reason}")
        return record
    print(f"[shape] {name} pep_frac={shape['peptide_contact_frac']:.2f} "
          f"helix_frac={shape['helix_contact_frac']:.2f} "
          f"h={shape['approach_height']:.1f} lat={shape['abs_lateral_offset']:.1f} "
          f"over={shape['over_peptide_frac']:.2f} "
          f"-> {'SHAPE-OK' if record['ok'] else 'reserve: ' + '; '.join(reasons)}")
    return record


def sequences_for(client: NvidiaClient, record: Dict[str, Any], strategy,
                  index: int, args) -> Dict[str, Any]:
    """ProteinMPNN + donor graft + de-wetted twin for one kept backbone."""
    name = record["backbone"]
    pdb = (OUT / "all_structures" / f"{name}.pdb").read_text()
    shell = batch.pos6_shell_positions(pdb)
    helix_positions = gg.helix_facing_positions(pdb)
    record["pos6_shell"] = shell[:8]
    record["helix_facing_positions"] = helix_positions[:12]
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
        donor = batch.donor_bias_score(sequence, shell)
        greasy = helix_apolar_penalty(sequence, helix_positions, shell)
        samples.append({"sequence": sequence, "mpnn_score": entry.get("score"),
                        "donor_bias_score": donor,
                        "helix_apolar_penalty": greasy,
                        "shape_selection_score": round(donor - 0.5 * greasy, 3)})
    if not samples:
        return record
    # donor content first, helix greasiness penalised, MPNN score breaks ties
    samples.sort(key=lambda s: (-s["shape_selection_score"],
                                s["mpnn_score"] if s["mpnn_score"] is not None else 9))
    picked = samples[: args.seqs_scored_per_backbone]
    for sample in picked:
        sample.update({"strategy": strategy["name"], "backbone": name,
                       "binder_length": strategy["binder_length"],
                       "variant": "mpnn-selected",
                       "design_pos6_contacts": record["pose"]["pos6_contacts"],
                       "design_pos6_distance": record["pose"]["pos6_distance"]})

    candidates: List[Dict[str, Any]] = []
    for rank, parent in enumerate(picked, start=1):
        grafted = batch.donor_anchored_variant(parent["sequence"], shell)
        if run_pipeline.sequence_filter_failures(grafted):
            continue
        graft = {**parent, "sequence": grafted,
                 "parent_sequence": parent["sequence"],
                 "donor_bias_score": batch.donor_bias_score(grafted, shell),
                 "helix_apolar_penalty": helix_apolar_penalty(
                     grafted, helix_positions, shell),
                 "variant": f"donor-anchored-g{rank}",
                 "shape_list": record["list"],
                 "shape": record["shape"]}
        candidates.append(graft)
        dewet, edited = dewetted_variant(grafted, helix_positions, shell)
        if (dewet != grafted and not run_pipeline.sequence_filter_failures(dewet)):
            candidates.append({**graft, "sequence": dewet,
                               "parent_sequence": grafted,
                               "donor_bias_score": batch.donor_bias_score(dewet, shell),
                               "helix_apolar_penalty": helix_apolar_penalty(
                                   dewet, helix_positions, shell),
                               "dewet_positions": edited,
                               "variant": f"donor-anchored-g{rank}-dw"})
    record["candidates"] = candidates
    return record


def main(argv=None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--backbones-per-strategy", type=int, default=48)
    parser.add_argument("--seqs-per-backbone", type=int, default=16)
    parser.add_argument("--seqs-scored-per-backbone", type=int, default=2)
    parser.add_argument("--diffusion-steps", type=int, default=15)
    parser.add_argument("--design-workers", type=int, default=4)
    parser.add_argument("--seed", type=int, default=8080)
    parser.add_argument("--timeout", type=float, default=600.0)
    parser.add_argument("--min-pos6-contacts", type=int, default=3)
    parser.add_argument("--max-pos6-distance", type=float, default=4.5)
    parser.add_argument("--reserve", type=int, default=8,
                        help="backbones kept by shape *rank* when they clear the "
                             "base filter but miss a shape threshold")
    parser.add_argument("--tag", default="")
    parser.add_argument("--from-existing", default="",
                        help="re-run stage 2 only, re-using the backbones named "
                             "in this backbone_poses JSON (no new RFdiffusion "
                             "calls); use to widen the reserve list")
    args = parser.parse_args(argv)

    started = time.time()
    for sub in ("all_structures", "structures"):
        (OUT / sub).mkdir(parents=True, exist_ok=True)

    target_dir = REPO / "data" / "targets"
    target_pdb = run_pipeline.trimmed_design_target(target_dir / "g12d_complex.pdb")

    if args.from_existing:
        records = json.loads(Path(args.from_existing).read_text())
        print(f"[design] re-using {len(records)} existing backbones from "
              f"{args.from_existing}; no RFdiffusion calls")
    else:
        jobs = [(s, i) for s in STRATEGIES for i in range(args.backbones_per_strategy)]
        records = []
        with ThreadPoolExecutor(max_workers=args.design_workers) as pool:
            futures = {pool.submit(backbone_job, NvidiaClient(timeout=args.timeout),
                                   s, i, target_pdb, args): (s, i) for s, i in jobs}
            for future in as_completed(futures):
                records.append(future.result())
    print(f"[design] {len(records)} backbones in {time.time() - started:.0f}s: "
          f"{sum(1 for r in records if r.get('base_ok'))} clear the base pose "
          f"filter, {sum(1 for r in records if r.get('ok'))} also clear every "
          f"shape threshold")

    shape_ok = [r for r in records if r.get("ok")]
    for record in shape_ok:
        record["list"] = "shape"
    reserve = sorted((r for r in records if r.get("base_ok") and not r.get("ok")),
                     key=lambda r: (r["shape"]["helix_contact_frac"],
                                    -r["shape"]["over_peptide_frac"]))[: args.reserve]
    for record in reserve:
        record["list"] = "reserve"
    kept = shape_ok + reserve
    print(f"[design] {len(shape_ok)} shape-filtered + {len(reserve)} reserve "
          f"backbones go to ProteinMPNN")

    strategy_by_name = {s["name"]: s for s in STRATEGIES}
    with ThreadPoolExecutor(max_workers=args.design_workers) as pool:
        futures = [pool.submit(sequences_for, NvidiaClient(timeout=args.timeout),
                               record, strategy_by_name[record["strategy"]],
                               index, args)
                   for index, record in enumerate(kept)]
        for future in as_completed(futures):
            future.result()

    candidates: List[Dict[str, Any]] = []
    seen = set()
    for record in kept:
        for candidate in record.get("candidates") or []:
            if candidate["sequence"] in seen:
                continue
            seen.add(candidate["sequence"])
            candidates.append(candidate)

    (OUT / f"candidates{args.tag}.json").write_text(json.dumps(candidates, indent=2))
    (OUT / f"backbone_poses{args.tag}.json").write_text(json.dumps(
        [{k: v for k, v in r.items() if k != "candidates"} for r in records],
        indent=2, default=str))
    stats = {"elapsed_seconds": round(time.time() - started, 1),
             "backbones_generated": len(records),
             "backbones_base_ok": sum(1 for r in records if r.get("base_ok")),
             "backbones_shape_ok": len(shape_ok),
             "backbones_reserve": len(reserve),
             "candidates": len(candidates),
             "shape_thresholds": SHAPE, "strategies": STRATEGIES,
             "diffusion_steps": args.diffusion_steps, "seed": args.seed}
    (OUT / f"design_stats{args.tag}.json").write_text(json.dumps(stats, indent=2, default=str))
    print(f"[batch-08] {len(candidates)} candidates written to "
          f"{OUT / f'candidates{args.tag}.json'} in {time.time() - started:.0f}s")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
