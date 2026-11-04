#!/usr/bin/env python3
"""Batch-06 analysis: paratope size/centring, gate outcomes, pose spread.

Re-uses `scripts/replicate_rmsd_batch04.binder_rmsd` for the pose spread (binder
C-alpha RMSD after superposing on the pMHC frame only, binder not re-fitted) and
`src/design_gates.evaluate` to *report* - never to apply - the unchanged hard
gates, including the `peptide_buried_frac >= 0.40` floor that rejects the real
7OW6 TCR. Nothing here computes an ipTM gap or a ddG.
"""
from __future__ import annotations

import itertools, json, statistics, sys
from pathlib import Path
from typing import Any, Dict, List, Optional

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(REPO / "scripts"))

import design_gates  # noqa: E402
import replicate_rmsd_batch04 as rmsd_mod  # noqa: E402

OUT = REPO / "results" / "batch-06"
B05 = REPO / "results" / "batch-05"


def num(value) -> Optional[float]:
    try:
        if value is None or value == "":
            return None
        return float(value)
    except (TypeError, ValueError):
        return None


def describe(values: List[Optional[float]]) -> Dict[str, Any]:
    vals = [v for v in values if v is not None]
    if not vals:
        return {"n": 0}
    return {"n": len(vals), "min": round(min(vals), 3),
            "median": round(statistics.median(vals), 3),
            "mean": round(statistics.mean(vals), 3),
            "max": round(max(vals), 3)}


def load_rows(path: Path) -> List[Dict[str, Any]]:
    if not path.exists():
        return []
    return [json.loads(l) for l in path.read_text().splitlines() if l.strip()]


def pose_spread(design: str, arm: str, condition: str,
                panel_rows: List[Dict[str, Any]]) -> Dict[str, Any]:
    members = [(r["draw"], OUT / r["structure_file"]) for r in panel_rows
               if r["design"] == design and r["arm"] == arm
               and r["condition"] == condition]
    members = [(d, p) for d, p in sorted(members) if p.exists()]
    if len(members) < 2:
        return {"n_poses": len(members)}
    models = {d: rmsd_mod.model_of(p) for d, p in members}
    rows, frames = [], []
    # Boltz-2 relabels chains sequentially, so in the *apo* arm (no peptide
    # chain C) the binder comes back as chain C, not D, and the superposition
    # frame is HLA chain A alone - exactly as recorded in the off-target panel
    # (PR #11). Stated rather than silently substituted.
    saved = (rmsd_mod.HLA_CHAIN, rmsd_mod.PEPTIDE_CHAIN, rmsd_mod.BINDER_CHAIN)
    if arm == "apo":
        rmsd_mod.PEPTIDE_CHAIN = rmsd_mod.HLA_CHAIN   # frame = HLA chain A only
        rmsd_mod.BINDER_CHAIN = "C"
    try:
        for (da, ma), (db, mb) in itertools.combinations(models.items(), 2):
            value, n_ca, frame = rmsd_mod.binder_rmsd(ma, mb)
            rows.append({"pair": f"r{da}~r{db}", "binder_ca_rmsd": round(value, 2),
                         "n_ca": n_ca, "frame_ca_rmsd": frame})
            frames.append(frame)
    except Exception as exc:  # noqa: BLE001
        return {"n_poses": len(members), "error": f"{type(exc).__name__}: {exc}"}
    finally:
        rmsd_mod.HLA_CHAIN, rmsd_mod.PEPTIDE_CHAIN, rmsd_mod.BINDER_CHAIN = saved
    values = [r["binder_ca_rmsd"] for r in rows]
    return {"n_poses": len(members),
            "frame": "HLA chain A only (apo: no peptide)" if arm == "apo"
                     else "pMHC (HLA chain A + peptide chain C)",
            "pairwise": rows,
            "mean": round(statistics.mean(values), 2),
            "max": round(max(values), 2),
            "frame_max": round(max(frames), 2)}


def main() -> int:
    screening = [r["row"] for r in load_rows(OUT / "checkpoint.jsonl")]
    panel_rows = load_rows(OUT / "panel_checkpoint.jsonl")
    b05 = [r["row"] for r in load_rows(B05 / "checkpoint.jsonl")]

    def dist(rows, key):
        return describe([num(r.get(key)) for r in rows])

    out: Dict[str, Any] = {
        "screening": {
            "n": len(screening),
            "buried_area_total": dist(screening, "buried_area_total"),
            "buried_area_peptide": dist(screening, "buried_area_peptide"),
            "buried_area_hla": dist(screening, "buried_area_hla"),
            "peptide_area_share": dist(screening, "peptide_buried_frac"),
            "asp6_gate_passes": sum(1 for r in screening
                                    if str(r.get("gate_pass")).lower() == "true"),
            "asp6_buried_frac_gt0": sum(
                1 for r in screening if (num(r.get("asp6_buried_frac")) or 0) > 0),
            "pos6_sidechain_contact_ge1": sum(
                1 for r in screening if (num(r.get("pos6_sidechain_contacts")) or 0) >= 1),
            "net_interface_charge": dist(screening, "net_interface_charge"),
            "peptide_buried_frac_ge_0.40": sum(
                1 for r in screening if (num(r.get("peptide_buried_frac")) or 0) >= 0.40),
            "peptide_buried_frac_ge_0.18": sum(
                1 for r in screening if (num(r.get("peptide_buried_frac")) or 0) >= 0.18),
        },
        "batch05_reference": {
            "n": len(b05),
            "buried_area_total": dist(b05, "buried_area_total"),
            "buried_area_peptide": dist(b05, "buried_area_peptide"),
            "buried_area_hla": dist(b05, "buried_area_hla"),
            "peptide_area_share": dist(b05, "peptide_buried_frac"),
        },
        "binder_length_breakdown": {},
        "panel": {},
    }

    for length in sorted({r.get("binder_length") for r in screening}):
        subset = [r for r in screening if r.get("binder_length") == length]
        out["binder_length_breakdown"][str(length)] = {
            "n": len(subset),
            "buried_area_total": dist(subset, "buried_area_total"),
            "peptide_area_share": dist(subset, "peptide_buried_frac"),
            "asp6_gate_passes": sum(1 for r in subset
                                    if str(r.get("gate_pass")).lower() == "true"),
        }

    # --- unchanged hard gates, reported not applied ---
    gate_report = []
    for row in screening:
        failures = design_gates.evaluate(row)
        gate_report.append({"design": f"{row.get('backbone')}_{row.get('variant')}",
                            "peptide_buried_frac": num(row.get("peptide_buried_frac")),
                            "net_interface_charge": num(row.get("net_interface_charge")),
                            "hard_gate_failures": failures})
    out["unchanged_hard_gates_reported"] = gate_report

    # --- panel arms ---
    if panel_rows:
        designs = sorted({r["design"] for r in panel_rows})
        for design in designs:
            drows = [r for r in panel_rows if r["design"] == design]
            entry: Dict[str, Any] = {}
            for condition in sorted({r["condition"] for r in drows}):
                for arm in sorted({r["arm"] for r in drows
                                   if r["condition"] == condition}):
                    sel = [r for r in drows if r["condition"] == condition
                           and r["arm"] == arm]
                    rows = [r["row"] for r in sel]
                    entry[f"{arm}@{condition}"] = {
                        "draws": len(sel),
                        "distinct_structures": len({r["structure_sha1"] for r in sel}),
                        "asp6_gate_passes": sum(
                            1 for r in rows
                            if str(r.get("asp6_gate_pass")).lower() == "true"),
                        "donor_distance": [r.get("pos6_min_donor_distance") for r in rows],
                        "pos6_contacts": [r.get("pos6_contacts") for r in rows],
                        "pos6_sidechain_contacts": [r.get("pos6_sidechain_contacts") for r in rows],
                        "buried_area_total": [r.get("buried_area_total") for r in rows],
                        "buried_area_peptide": [r.get("buried_area_peptide") for r in rows],
                        "buried_area_hla": [r.get("buried_area_hla") for r in rows],
                        "peptide_buried_frac": [r.get("peptide_buried_frac") for r in rows],
                        "net_interface_charge": [r.get("net_interface_charge") for r in rows],
                        "iptm_pose_readout": [r.get("iptm") for r in rows],
                        "pose_spread": pose_spread(design, arm, condition, panel_rows),
                    }
            out["panel"][design] = entry

    (OUT / "analysis.json").write_text(json.dumps(out, indent=2, default=str))
    print(json.dumps({k: v for k, v in out.items()
                      if k not in ("unchanged_hard_gates_reported",)},
                     indent=2, default=str)[:6000])
    print("wrote", OUT / "analysis.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
