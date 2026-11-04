#!/usr/bin/env python3
"""Batch-08 analysis: where the interface actually lands, and on what.

`scripts/batch06_analyze.py` is imported and re-used for the pose spread
(binder C-alpha RMSD in the pMHC frame, binder **not** re-fitted, frame RMSD
reported as the artefact floor, apo framed on HLA chain A alone) and for the
reported-not-applied `src/design_gates.evaluate` read-out.  What is new is the
partition: every scored structure's buried area is split peptide / HLA alpha-1
(57-84) / HLA alpha-2 (138-180) / rest of HLA / b2m with
`src/groove_geometry.py`, so "cut groove contact" is measured rather than
asserted.

Three comparators are computed on structures that already exist, at zero API
cost:
  * batch-06's peptide share of buried area (0.213 median) - the shrink arm;
  * bb54-g2's own batch-07 draws at r3np - the design batch-08 must beat;
  * the **real 7OW6 TCR in the crystal** - the only validated G12D-selective
    binder in this system, partitioned the same way.
No ipTM gap, no ddG, anywhere.
"""
from __future__ import annotations

import json, statistics, sys
from pathlib import Path
from typing import Any, Dict, List, Optional

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(REPO / "scripts"))

import batch06_analyze as a6  # noqa: E402
import complex_analysis as ca  # noqa: E402
import design_gates  # noqa: E402
import groove_geometry as gg  # noqa: E402

OUT = REPO / "results" / "batch-08"
B06 = REPO / "results" / "batch-06"
B07 = REPO / "results" / "batch-07"
a6.OUT = OUT

num, describe, load_rows = a6.num, a6.describe, a6.load_rows


def partition_rows(rows: List[Dict[str, Any]], root: Path) -> List[Dict[str, Any]]:
    """Buried-area partition for every checkpointed structure that still exists."""
    out = []
    for record in rows:
        path = root / record["structure_file"]
        if not path.exists():
            continue
        try:
            part = gg.buried_area_partition(path.read_text(), record["sequence"],
                                            arm=record["arm"])
            approach = (gg.predicted_approach(path.read_text(), record["sequence"],
                                              arm=record["arm"])
                        if record["arm"] != "apo" else {})
        except Exception as exc:  # noqa: BLE001
            out.append({**{k: record[k] for k in ("design", "arm", "condition", "draw")},
                        "error": f"{type(exc).__name__}: {exc}"})
            continue
        out.append({**{k: record[k] for k in ("design", "arm", "condition", "draw")},
                    "variant": record.get("variant"), **part, **approach})
    return out


def share_summary(parts: List[Dict[str, Any]], arm: str = "g12d") -> Dict[str, Any]:
    sel = [p for p in parts if p.get("arm") == arm and "error" not in p]
    return {"n": len(sel),
            "peptide_area_share": describe([num(p.get("peptide_area_share")) for p in sel]),
            "helix_area_share": describe([num(p.get("helix_area_share")) for p in sel]),
            "buried_area_total": describe([num(p.get("buried_area_total")) for p in sel]),
            "buried_area_peptide": describe([num(p.get("buried_area_peptide")) for p in sel]),
            "buried_area_hla_alpha1": describe([num(p.get("buried_area_hla_alpha1")) for p in sel]),
            "buried_area_hla_alpha2": describe([num(p.get("buried_area_hla_alpha2")) for p in sel]),
            "buried_area_hla_other": describe([num(p.get("buried_area_hla_other")) for p in sel]),
            "approach_height": describe([num(p.get("approach_height")) for p in sel]),
            "over_peptide_frac": describe([num(p.get("over_peptide_frac")) for p in sel])}


def crystal_tcr_reference() -> Dict[str, Any]:
    """The real 7OW6 TCR (chains D+E), partitioned identically."""
    text = (REPO / "data" / "raw" / "7OW6.pdb").read_text()

    def model():
        return ca._clean_model(next(iter(ca.parse_structure(text, fmt="pdb"))))

    bound = ca._atom_sasa(model())
    free_model = model()
    for cid in ("D", "E"):
        free_model.detach_child(cid)
    free = ca._atom_sasa(free_model)
    parts = {"peptide": 0.0, "hla_alpha1": 0.0, "hla_alpha2": 0.0,
             "hla_other": 0.0, "b2m": 0.0}
    total = 0.0
    for (cid, resseq, atom), value in free.items():
        delta = value - bound.get((cid, resseq, atom), 0.0)
        if delta <= 0:
            continue
        total += delta
        if cid == "C":
            parts["peptide"] += delta
        elif cid == "A":
            parts[gg.segment_of(resseq)] += delta
        elif cid == "B":
            parts["b2m"] += delta
    helix = parts["hla_alpha1"] + parts["hla_alpha2"]
    return {"source": "7OW6 crystal, TCR chains D+E as the binder",
            **{f"buried_area_{k}": round(v, 1) for k, v in parts.items()},
            "buried_area_total": round(total, 1),
            "peptide_area_share": round(parts["peptide"] / total, 3) if total else None,
            "helix_area_share": round(helix / total, 3) if total else None}


def main() -> int:
    panel_rows = load_rows(OUT / "panel_checkpoint.jsonl")
    parts = partition_rows(panel_rows, OUT)
    b07_rows = load_rows(B07 / "panel_checkpoint.jsonl")
    bb54 = [r for r in b07_rows
            if r["design"].endswith("bb54_donor-anchored-g2") and r["condition"] == "r3np"]
    bb54_parts = partition_rows(bb54, B07)
    b06_screen = [r["row"] for r in load_rows(B06 / "checkpoint.jsonl")]

    out: Dict[str, Any] = {
        "comparators": {
            "batch06_screening_peptide_area_share": describe(
                [num(r.get("peptide_area_share")) for r in b06_screen]),
            "bb54_g2_r3np": {arm: share_summary(bb54_parts, arm)
                             for arm in ("g12d", "wt", "apo", "decoy")},
            "crystal_7ow6_tcr": crystal_tcr_reference(),
        },
        "batch08": {arm: share_summary(parts, arm)
                    for arm in sorted({p.get("arm") for p in parts})},
        "per_structure_partition": parts,
        "designs": {},
    }

    designs = sorted({r["design"] for r in panel_rows})
    for design in designs:
        drows = [r for r in panel_rows if r["design"] == design]
        entry: Dict[str, Any] = {}
        for condition in sorted({r["condition"] for r in drows}):
            for arm in sorted({r["arm"] for r in drows if r["condition"] == condition}):
                sel = [r for r in drows if r["condition"] == condition and r["arm"] == arm]
                rows = [r["row"] for r in sel]
                dparts = [p for p in parts if p["design"] == design
                          and p["arm"] == arm and p["condition"] == condition]
                entry[f"{arm}@{condition}"] = {
                    "draws": len(sel),
                    "distinct_structures": len({r["structure_sha1"] for r in sel}),
                    "asp6_gate_passes": sum(
                        1 for r in rows if str(r.get("asp6_gate_pass")).lower() == "true"),
                    "donor_distance": [r.get("pos6_min_donor_distance") for r in rows],
                    "pos6_contacts": [r.get("pos6_contacts") for r in rows],
                    "pos6_sidechain_contacts": [r.get("pos6_sidechain_contacts") for r in rows],
                    "buried_area_total": [r.get("buried_area_total") for r in rows],
                    "buried_area_peptide": [r.get("buried_area_peptide") for r in rows],
                    "buried_area_hla": [r.get("buried_area_hla") for r in rows],
                    "peptide_buried_frac": [r.get("peptide_buried_frac") for r in rows],
                    "peptide_area_share": [p.get("peptide_area_share") for p in dparts],
                    "helix_area_share": [p.get("helix_area_share") for p in dparts],
                    "buried_area_hla_alpha1": [p.get("buried_area_hla_alpha1") for p in dparts],
                    "buried_area_hla_alpha2": [p.get("buried_area_hla_alpha2") for p in dparts],
                    "iptm_pose_readout": [r.get("iptm") for r in rows],
                    "pose_spread": a6.pose_spread(design, arm, condition, panel_rows),
                }
        entry["hard_gates_reported_not_applied"] = [
            {"draw": r["draw"],
             "peptide_buried_frac": num(r["row"].get("peptide_buried_frac")),
             "net_interface_charge": num(r["row"].get("net_interface_charge")),
             "failures": design_gates.evaluate(r["row"])}
            for r in drows if r["arm"] == "g12d"]
        out["designs"][design] = entry

    g12d = [p for p in parts if p.get("arm") == "g12d" and "error" not in p]
    g12d_rows = [r["row"] for r in panel_rows if r["arm"] == "g12d"]
    # a degenerate structure can bury nothing, leaving the share undefined (None);
    # those draws are excluded from the medians and counted separately.
    pep_shares = [p["peptide_area_share"] for p in g12d
                  if p.get("peptide_area_share") is not None]
    helix_shares = [p["helix_area_share"] for p in g12d
                    if p.get("helix_area_share") is not None]
    out["headline"] = {
        "g12d_draws": len(g12d),
        "g12d_draws_with_defined_share": len(pep_shares),
        "peptide_area_share_median": (round(statistics.median(pep_shares), 3)
                                      if pep_shares else None),
        "helix_area_share_median": (round(statistics.median(helix_shares), 3)
                                    if helix_shares else None),
        "peptide_buried_frac_max": max(
            [num(r.get("peptide_buried_frac")) or 0 for r in g12d_rows], default=None),
        "peptide_buried_frac_ge_0.40": sum(
            1 for r in g12d_rows if (num(r.get("peptide_buried_frac")) or 0) >= 0.40),
        "peptide_buried_frac_ge_0.18": sum(
            1 for r in g12d_rows if (num(r.get("peptide_buried_frac")) or 0) >= 0.18),
        "asp6_gate_passes": sum(
            1 for r in g12d_rows if str(r.get("asp6_gate_pass")).lower() == "true"),
    }

    (OUT / "analysis.json").write_text(json.dumps(out, indent=2, default=str))
    print(json.dumps({k: out[k] for k in ("headline", "comparators", "batch08")},
                     indent=2, default=str)[:7000])
    print("wrote", OUT / "analysis.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
