#!/usr/bin/env python3
"""Analyse the Boltz-2 sampling-settings study against the batch-05 baseline.

Pose spread is measured exactly as `scripts/replicate_rmsd_batch04.py` does it
and that module's code is re-used: superpose on the pMHC frame only (HLA heavy
chain A + peptide C C-alphas), then measure binder (chain D) C-alpha RMSD
*without* re-fitting the binder. The frame RMSD is reported alongside so a
spread can be told apart from a fitting artefact.

Baseline = the 3 batch-04/batch-05 replicate draws at recycling_steps 3,
sampling_steps 50, diffusion_samples 1, without_potentials False. The
screening draw is excluded so every condition is a like-for-like 3-draw set.
"""
from __future__ import annotations

import itertools, json, statistics, sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(REPO / "scripts"))

import replicate_rmsd_batch04 as rmsdmod  # noqa: E402
import sampling_study as study  # noqa: E402

OUT = REPO / "results" / "sampling-study"
GATE_KEYS = ("asp6_min_donor_distance", "pos6_sidechain_contacts",
             "peptide_buried_frac", "asp6_buried_frac",
             "net_interface_charge", "iptm_g12d")


def num(value):
    if value in ("", None):
        return None
    return float(value)


def spread(values):
    vals = [v for v in values if v is not None]
    if not vals:
        return None
    return {"values": vals, "mean": round(statistics.fmean(vals), 3),
            "sd": round(statistics.pstdev(vals), 3) if len(vals) > 1 else 0.0}


def pose_spread(paths):
    models = [rmsdmod.model_of(p) for p in paths]
    rows = []
    for (i, a), (j, b) in itertools.combinations(list(enumerate(models)), 2):
        value, n, frame = rmsdmod.binder_rmsd(a, b)
        rows.append({"pair": f"r{i}~r{j}", "binder_ca_rmsd": round(value, 2),
                     "frame_ca_rmsd": frame, "n_ca": n})
    vals = [r["binder_ca_rmsd"] for r in rows]
    return {"pairwise": rows,
            "mean": round(statistics.fmean(vals), 2),
            "max": round(max(vals), 2),
            "sd": round(statistics.pstdev(vals), 2) if len(vals) > 1 else 0.0,
            "frame_max": round(max(r["frame_ca_rmsd"] for r in rows), 2)}


def main() -> int:
    subjects = {s["key"]: s for s in study.load_subjects()}
    checkpoint = [json.loads(l) for l in
                  (OUT / "checkpoint.jsonl").read_text().splitlines() if l.strip()]

    records = []
    for key, subject in subjects.items():
        source = subject["source"]
        stem = f"{subject['backbone']}_{subject['variant']}"
        base_dir = REPO / "results" / source / "replicate_structures"
        base_paths = sorted(base_dir.glob(f"{stem}_r*.cif"))
        base_runs = next(r for r in json.loads(
            (REPO / "results" / source / "replicates.json").read_text())
            if r["backbone"] == subject["backbone"])["runs"]

        conditions = {"baseline": {
            "boltz_kwargs": {"recycling_steps": 3, "sampling_steps": 50,
                             "diffusion_samples": 1, "without_potentials": False},
            "draws": len(base_runs),
            "gate_passes": sum(1 for r in base_runs if r["gate_pass"]),
            "per_draw": [{"gate_pass": bool(r["gate_pass"]),
                          **{k: num(r.get(k)) for k in GATE_KEYS}}
                         for r in base_runs],
            **{k: spread([num(r.get(k)) for r in base_runs]) for k in GATE_KEYS},
            "pose": pose_spread(base_paths) if len(base_paths) > 1 else None,
        }}

        for condition in study.CONDITIONS:
            rows = sorted([c for c in checkpoint
                           if c["design"] == key and c["condition"] == condition],
                          key=lambda c: c["draw"])
            if not rows:
                continue
            paths = [OUT / r["structure_file"] for r in rows]
            recs = [r["row"] for r in rows]
            conditions[condition] = {
                "boltz_kwargs": rows[0]["boltz_kwargs"],
                "draws": len(recs),
                "gate_passes": sum(1 for r in recs if r["gate_pass"]),
                "per_draw": [{"gate_pass": bool(r["gate_pass"]),
                              **{k: num(r.get(k)) for k in GATE_KEYS}}
                             for r in recs],
                **{k: spread([num(r.get(k)) for r in recs]) for k in GATE_KEYS},
                "pose": pose_spread(paths) if len(paths) > 1 else None,
            }
        records.append({"design": key, "backbone": subject["backbone"],
                        "variant": subject["variant"], "source": source,
                        "conditions": conditions})

    (OUT / "analysis.json").write_text(json.dumps(records, indent=2))

    # console table
    header = (f"{'design':<10} {'condition':<9} {'gate/3':>6} {'donor A (per draw)':<26} "
              f"{'pos6_sc':<12} {'pep_bf':>7} {'charge':<12} {'ipTM':>6} "
              f"{'RMSD mean/max':>15} {'frame':>6}")
    print(header)
    print("-" * len(header))
    for rec in records:
        for name, c in rec["conditions"].items():
            donors = ", ".join(f"{d['asp6_min_donor_distance']:.2f}"
                               if d["asp6_min_donor_distance"] else "-"
                               for d in c["per_draw"])
            pos6 = ", ".join(str(int(d["pos6_sidechain_contacts"] or 0))
                             for d in c["per_draw"])
            charge = ", ".join(f"{d['net_interface_charge']:+.1f}"
                               for d in c["per_draw"])
            pose = c["pose"]
            print(f"{rec['design']:<10} {name:<9} "
                  f"{c['gate_passes']:>2}/{c['draws']:<3} {donors:<26} {pos6:<12} "
                  f"{(c['peptide_buried_frac'] or {}).get('mean', 0):>7.3f} "
                  f"{charge:<12} {(c['iptm_g12d'] or {}).get('mean', 0):>6.3f} "
                  f"{(pose['mean'] if pose else 0):>7.2f}/{(pose['max'] if pose else 0):<7.2f} "
                  f"{(pose['frame_max'] if pose else 0):>6.2f}")
        print()
    print("wrote", OUT / "analysis.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
