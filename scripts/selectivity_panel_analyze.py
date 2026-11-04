#!/usr/bin/env python3
"""Analyse the replication confirmation and the off-target panel.

Pose spread is measured exactly as `scripts/replicate_rmsd_batch04.py` does it,
and that module's code is re-used: superpose on the **pMHC frame only** (HLA
heavy chain A + peptide chain C C-alphas) and measure the binder (chain D)
C-alpha RMSD *without* re-fitting the binder, with the frame RMSD reported
alongside as the fitting-artefact floor. For the **apo** arm there is no
peptide chain, so the frame is HLA chain A alone - stated, not silently
substituted.

No ipTM gap and no MM-GBSA ddG are computed anywhere: both failed the 7OW6
positive control and are retired. ipTM is reported per draw as a pose-quality
read-out only.
"""
from __future__ import annotations

import hashlib, itertools, json, statistics, sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(REPO / "scripts"))

import numpy as np  # noqa: E402
import replicate_rmsd_batch04 as rm  # noqa: E402

OUT = REPO / "results" / "selectivity-panel"
STUDY = REPO / "results" / "sampling-study"
SETTINGS = {"bb50-g1": "recyc8", "bb25": "recyc8", "bb135-g1": "combo"}


def binder_chain_of(model, binder_length: int) -> str:
    """The binder's chain id in this model.

    Boltz-2 relabels chains sequentially, so in the **apo** arm - where there
    is no peptide - the binder comes back as chain C rather than D. The binder
    is identified by length, never by a hardcoded id.
    """
    hits = [c.id for c in model
            if len([r for r in c if rm.ca.is_polymer_residue(r)]) == binder_length]
    if len(hits) != 1:
        raise ValueError(f"cannot identify the binder chain among {hits}")
    return hits[0]


def binder_rmsd_frame(ref, mob, chains, ref_binder="D", mob_binder="D"):
    """rm.binder_rmsd, with the fitting frame made explicit (apo has no C)."""
    fit_m, fit_r = [], []
    for cid in chains:
        m, r = rm.pairs(rm.ca_list(mob, cid), rm.ca_list(ref, cid))
        fit_m.append(m); fit_r.append(r)
    M, R_ = np.vstack(fit_m), np.vstack(fit_r)
    rot, t = rm.kabsch(M, R_)
    frame = float(np.sqrt((((M @ rot.T + t) - R_) ** 2).sum(1).mean()))
    bm, br = rm.pairs(rm.ca_list(mob, mob_binder), rm.ca_list(ref, ref_binder))
    moved = bm @ rot.T + t
    return float(np.sqrt(((moved - br) ** 2).sum(1).mean())), round(frame, 2)


def digest(path: Path) -> str:
    return hashlib.md5(path.read_bytes()).hexdigest()[:8]


def stat(values):
    vals = [v for v in values if v is not None]
    if not vals:
        return None
    return {"mean": round(statistics.fmean(vals), 3),
            "sd": round(statistics.pstdev(vals), 3) if len(vals) > 1 else 0.0,
            "min": round(min(vals), 3), "max": round(max(vals), 3), "n": len(vals)}


def load():
    rows = [json.loads(l) for l in
            (OUT / "checkpoint.jsonl").read_text().splitlines() if l.strip()]
    study = [json.loads(l) for l in
             (STUDY / "checkpoint.jsonl").read_text().splitlines() if l.strip()]
    return rows, study


def g12d_draws(design, rows, study):
    """The sampling-study draws plus this run's, in draw order."""
    out = []
    for r in study:
        if r["design"] == design and r["condition"] == SETTINGS[design]:
            out.append({"draw": r["draw"], "row": r["row"],
                        "path": STUDY / r["structure_file"], "source": "sampling-study"})
    for r in rows:
        if r["design"] == design and r["arm"] == "g12d":
            out.append({"draw": r["draw"], "row": r["row"],
                        "path": OUT / r["structure_file"], "source": "this run"})
    return sorted(out, key=lambda d: d["draw"])


def main() -> int:
    rows, study = load()
    report = {"designs": {}}
    for design, condition in SETTINGS.items():
        draws = g12d_draws(design, rows, study)
        per_draw = []
        for d in draws:
            r = d["row"]
            per_draw.append({
                "draw": d["draw"], "source": d["source"],
                "gate_pass": bool(r.get("gate_pass", r.get("asp6_gate_pass"))),
                "asp6_min_donor_distance": r.get("asp6_min_donor_distance") or None,
                "pos6_sidechain_contacts": r.get("pos6_sidechain_contacts"),
                "pos6_contacts": r.get("pos6_contacts"),
                "peptide_buried_frac": r.get("peptide_buried_frac"),
                "asp6_buried_frac": r.get("asp6_buried_frac"),
                "net_interface_charge": r.get("net_interface_charge"),
                "iptm": r.get("iptm_g12d", r.get("iptm")),
                "structure_md5": digest(d["path"]),
            })
        models = [rm.model_of(d["path"]) for d in draws]
        pairwise = []
        for (i, a), (j, b) in itertools.combinations(list(enumerate(models)), 2):
            value, frame = binder_rmsd_frame(a, b, ("A", "C"))
            pairwise.append({"pair": f"r{draws[i]['draw']}~r{draws[j]['draw']}",
                             "binder_ca_rmsd": round(value, 2), "frame_ca_rmsd": frame})
        vals = [p["binder_ca_rmsd"] for p in pairwise]
        entry = {
            "condition": condition,
            "n_draws": len(draws),
            "n_distinct_structures": len({p["structure_md5"] for p in per_draw}),
            "gate_passes": sum(1 for p in per_draw if p["gate_pass"]),
            "per_draw": per_draw,
            "pose": {"pairwise": pairwise,
                     "mean": round(statistics.fmean(vals), 2),
                     "max": round(max(vals), 2),
                     "frame_max": round(max(p["frame_ca_rmsd"] for p in pairwise), 2)},
            "summary": {k: stat([p[k] for p in per_draw]) for k in
                        ("asp6_min_donor_distance", "pos6_sidechain_contacts",
                         "peptide_buried_frac", "net_interface_charge", "iptm")},
            "arms": {},
        }

        # --- off-target arms ---
        ref = models[0]
        for arm in ("wt", "apo", "decoy"):
            arm_rows = sorted([r for r in rows
                               if r["design"] == design and r["arm"] == arm],
                              key=lambda r: r["draw"])
            if not arm_rows:
                continue
            arm_draws = []
            for r in arm_rows:
                row = r["row"]
                path = OUT / r["structure_file"]
                model = rm.model_of(path)
                frame_chains = ("A",) if arm == "apo" else ("A", "C")
                mob_binder = binder_chain_of(model, len(r["sequence"]))
                rmsd_vs_g12d, frame = binder_rmsd_frame(
                    ref, model, frame_chains, "D", mob_binder)
                arm_draws.append({
                    "draw": r["draw"], "constrained": r["constrained"],
                    "pos6_residue": row.get("pos6_residue"),
                    "pos6_contacts": row.get("pos6_contacts"),
                    "pos6_sidechain_contacts": row.get("pos6_sidechain_contacts"),
                    "pos6_min_donor_distance": row.get("pos6_min_donor_distance"),
                    "pos6_donor_partner": row.get("pos6_donor_partner", ""),
                    "peptide_buried_frac": row.get("peptide_buried_frac"),
                    "buried_area_total": row.get("buried_area_total"),
                    "buried_area_peptide": row.get("buried_area_peptide"),
                    "buried_area_by_chain": row.get("buried_area_by_chain"),
                    "contacts_by_chain": row.get("contacts_by_chain"),
                    "net_interface_charge": row.get("net_interface_charge"),
                    "iptm": row.get("iptm"),
                    "donor_contacts": row.get("donor_contacts"),
                    "binder_rmsd_vs_g12d_r0": round(rmsd_vs_g12d, 2),
                    "frame_ca_rmsd": frame,
                    "frame_chains": "+".join(frame_chains),
                    "structure_md5": digest(path),
                })
            # internal pose spread within the arm
            arm_models = [rm.model_of(OUT / r["structure_file"]) for r in arm_rows]
            internal = []
            for (i, a), (j, b) in itertools.combinations(list(enumerate(arm_models)), 2):
                binder_id = binder_chain_of(a, len(arm_rows[i]["sequence"]))
                value, frame = binder_rmsd_frame(
                    a, b, ("A",) if arm == "apo" else ("A", "C"),
                    binder_id, binder_chain_of(b, len(arm_rows[j]["sequence"])))
                internal.append({"pair": f"r{i}~r{j}",
                                 "binder_ca_rmsd": round(value, 2),
                                 "frame_ca_rmsd": frame})
            entry["arms"][arm] = {
                "n_draws": len(arm_draws),
                "n_distinct_structures": len({d["structure_md5"] for d in arm_draws}),
                "per_draw": arm_draws,
                "internal_pose": {
                    "pairwise": internal,
                    "mean": round(statistics.fmean(
                        [p["binder_ca_rmsd"] for p in internal]), 2) if internal else None,
                    "max": round(max(p["binder_ca_rmsd"] for p in internal), 2)
                    if internal else None},
                "summary": {k: stat([d[k] for d in arm_draws]) for k in
                            ("pos6_contacts", "pos6_sidechain_contacts",
                             "pos6_min_donor_distance", "peptide_buried_frac",
                             "buried_area_total", "net_interface_charge", "iptm",
                             "binder_rmsd_vs_g12d_r0")},
            }
        report["designs"][design] = entry

    (OUT / "analysis.json").write_text(json.dumps(report, indent=2))

    # --- console ---
    for design, e in report["designs"].items():
        print(f"\n=== {design}  ({e['condition']})  G12D n={e['n_draws']} "
              f"({e['n_distinct_structures']} distinct structures) "
              f"gate {e['gate_passes']}/{e['n_draws']}")
        print(f"{'draw':>4} {'gate':>5} {'donor':>6} {'pos6_sc':>7} {'pos6':>5} "
              f"{'pep_bf':>7} {'asp6_bf':>7} {'charge':>7} {'ipTM':>6}  md5")
        for p in e["per_draw"]:
            print(f"{p['draw']:>4} {str(p['gate_pass']):>5} "
                  f"{(p['asp6_min_donor_distance'] or '-'):>6} "
                  f"{str(p['pos6_sidechain_contacts']):>7} {str(p['pos6_contacts']):>5} "
                  f"{str(p['peptide_buried_frac']):>7} {str(p['asp6_buried_frac'])[:5]:>7} "
                  f"{str(p['net_interface_charge']):>7} {str(p['iptm']):>6}  "
                  f"{p['structure_md5']} [{p['source']}]")
        print(f"  binder RMSD mean/max {e['pose']['mean']}/{e['pose']['max']} A "
              f"(frame <= {e['pose']['frame_max']} A)")
        for arm, a in e["arms"].items():
            s = a["summary"]
            print(f"  -- {arm}: n={a['n_draws']} ({a['n_distinct_structures']} distinct) "
                  f"pos6_res={a['per_draw'][0]['pos6_residue']} "
                  f"pos6={s['pos6_contacts']['mean'] if s['pos6_contacts'] else '-'} "
                  f"pos6_sc={s['pos6_sidechain_contacts']['mean'] if s['pos6_sidechain_contacts'] else '-'} "
                  f"donor={s['pos6_min_donor_distance']['mean'] if s['pos6_min_donor_distance'] else 'none'} "
                  f"pep_bf={s['peptide_buried_frac']['mean'] if s['peptide_buried_frac'] else '-'} "
                  f"area={s['buried_area_total']['mean']} "
                  f"charge={s['net_interface_charge']['mean']} "
                  f"ipTM={s['iptm']['mean']} "
                  f"RMSD_vs_g12d={s['binder_rmsd_vs_g12d_r0']['mean']} A")
            for d in a["per_draw"]:
                print(f"       r{d['draw']} pos6={d['pos6_contacts']} "
                      f"sc={d['pos6_sidechain_contacts']} "
                      f"donor={d['pos6_min_donor_distance']} {d['pos6_donor_partner']} "
                      f"pep_bf={d['peptide_buried_frac']} "
                      f"areas={d['buried_area_by_chain']} "
                      f"rmsd={d['binder_rmsd_vs_g12d_r0']} (frame {d['frame_ca_rmsd']}, "
                      f"{d['frame_chains']}) donors={d['donor_contacts']}")
    print("\nwrote", OUT / "analysis.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
