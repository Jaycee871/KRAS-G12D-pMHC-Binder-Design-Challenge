#!/usr/bin/env python3
"""Batch-05 ranking, with replication folded into the gate.

`--stage screen` is batch-04's `scripts/rank_batch04.py` re-used verbatim with
its output directory pointed at `results/batch-05`: it applies the single-draw
Asp6 donor gate to `checkpoint.jsonl` and writes `gate_results.json` (the list
`scripts/replicate_batch05.py` then re-runs) plus `designs.csv`.

`--stage final` applies the batch-05 gate proper: a design passes only if it
clears the Asp6 direct read-out gate (binder Arg/Lys/His/Ser/Thr/Tyr donor
heavy atom within 3.5 A of Asp6 OD1/OD2, and `asp6_buried_frac` > 0) in **3 of
3** independent constrained draws. Survivors are ranked on measured geometry
only, in the order fixed after both `specificity_gap` and `ddG_sel` failed the
7OW6 control:

    1. peptide_buried_frac   (0.40 floor REPORTED, never applied here - the
                              recalibration to 0.18 is the user's call)
    2. contacts to the Asp6 **side chain** (backbone contacts excluded)
    3. net interface charge
    4. replicate pose stability (binder C-alpha RMSD spread, pMHC frame)
    5. topological diversity (70 % sequence-identity clusters)

ipTM is reported as a pose-quality read-out only: never a gap, never a
tie-breaker. No gate is loosened to manufacture survivors.
"""
from __future__ import annotations

import argparse, csv, json, statistics, sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(REPO / "scripts"))

import rank_batch04, run_batch04  # noqa: E402

OUT = REPO / "results" / "batch-05"
PEPTIDE_BURIED_FLOOR = 0.40          # reported, not applied - user's call
REQUIRED_GATE_PASSES = 3


def mean_of(spread, key="mean"):
    return None if not spread else spread.get(key)


def final() -> int:
    rows = {json.loads(l)["row"]["sequence"]: json.loads(l)["row"]
            for l in (OUT / "checkpoint.jsonl").read_text().splitlines() if l.strip()}
    replicates = json.loads((OUT / "replicates.json").read_text())
    rmsd = {r["backbone"] + "|" + r["variant"]: r
            for r in json.loads((OUT / "replicate_rmsd.json").read_text())}

    records = []
    for rep in replicates:
        key = rep["backbone"] + "|" + rep["variant"]
        screen = rows.get(rep["sequence"], {})
        geometry = rmsd.get(key, {})
        replicate_runs = rep["runs"]
        passes = sum(1 for r in replicate_runs if r["gate_pass"])
        # pose spread over the 3 replicate draws only (pair labels r0/r1/r2)
        rep_pairs = [p for p in geometry.get("pairwise", [])
                     if "original" not in p["pair"]]
        records.append({
            "backbone": rep["backbone"], "variant": rep["variant"],
            "sequence": rep["sequence"],
            "binder_length": screen.get("binder_length", ""),
            "strategy": screen.get("strategy", ""),
            "draws": len(replicate_runs),
            "gate_passes": passes,
            "gate_3of3": passes >= REQUIRED_GATE_PASSES and len(replicate_runs) >= 3,
            "screen_gate_pass": bool(screen.get("gate_pass")),
            "screen_asp6_min_donor_distance": screen.get("asp6_min_donor_distance", ""),
            "screen_asp6_donor_partner": screen.get("asp6_donor_partner", ""),
            "asp6_min_donor_distance_draws": [r.get("asp6_min_donor_distance", "")
                                              for r in replicate_runs],
            "asp6_donor_partner_draws": [r.get("asp6_donor_partner", "")
                                         for r in replicate_runs],
            "asp6_min_donor_distance": rep.get("asp6_min_donor_distance"),
            "asp6_buried_frac": rep.get("asp6_buried_frac"),
            "peptide_buried_frac": rep.get("peptide_buried_frac"),
            "pos6_sidechain_contacts": rep.get("pos6_sidechain_contacts"),
            "net_interface_charge": rep.get("net_interface_charge"),
            "iptm_g12d": rep.get("iptm_g12d"),
            "binder_ca_rmsd_mean_all": geometry.get("binder_ca_rmsd_mean"),
            "binder_ca_rmsd_max_all": geometry.get("binder_ca_rmsd_max"),
            "binder_ca_rmsd_mean_replicates": (
                round(statistics.fmean([p["binder_ca_rmsd"] for p in rep_pairs]), 2)
                if rep_pairs else None),
            "binder_ca_rmsd_max_replicates": (
                max(p["binder_ca_rmsd"] for p in rep_pairs) if rep_pairs else None),
            "pairwise_rmsd": geometry.get("pairwise", []),
        })

    survivors = [r for r in records if r["gate_3of3"]]
    survivors.sort(key=lambda r: (
        -(mean_of(r["peptide_buried_frac"]) or 0),
        -(mean_of(r["pos6_sidechain_contacts"]) or 0),
        -(mean_of(r["net_interface_charge"]) or 0),
        (r["binder_ca_rmsd_mean_replicates"] if r["binder_ca_rmsd_mean_replicates"]
         is not None else 1e9),
    ))
    clusters = []
    for row in survivors:
        for cluster in clusters:
            if rank_batch04.identity(cluster[0]["sequence"], row["sequence"]) >= 0.70:
                cluster.append(row)
                break
        else:
            clusters.append([row])
    for index, cluster in enumerate(clusters, start=1):
        for member in cluster:
            member["cluster"] = index
    for rank, row in enumerate(survivors, start=1):
        row["rank"] = rank
        row["name"] = f"KRASG12D_B05_{rank:03d}_{row['backbone']}"
        row["meets_peptide_floor_0.40"] = (
            (mean_of(row["peptide_buried_frac"]) or 0) >= PEPTIDE_BURIED_FLOOR)

    summary = {
        "screening_gate_passers_replicated": len(records),
        "survivors_3of3": len(survivors),
        "clusters": len(clusters),
        "required_gate_passes": REQUIRED_GATE_PASSES,
        "peptide_buried_floor_reported_not_applied": PEPTIDE_BURIED_FLOOR,
        "survivors_meeting_peptide_floor": sum(
            1 for r in survivors if r["meets_peptide_floor_0.40"]),
    }
    (OUT / "final_ranking.json").write_text(json.dumps(
        {"summary": summary, "survivors": survivors, "all_replicated": records},
        indent=2, default=str))

    columns = ["rank", "name", "backbone", "variant", "cluster", "sequence",
               "gate_passes", "draws", "peptide_buried_frac_mean",
               "pos6_sidechain_contacts_mean", "asp6_min_donor_distance_mean",
               "asp6_buried_frac_mean", "net_interface_charge_mean",
               "iptm_g12d_mean", "binder_ca_rmsd_mean_replicates",
               "binder_ca_rmsd_max_replicates", "meets_peptide_floor_0.40"]
    with (OUT / "survivors.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        for row in survivors:
            writer.writerow({
                **{c: row.get(c, "") for c in columns},
                "peptide_buried_frac_mean": mean_of(row["peptide_buried_frac"]),
                "pos6_sidechain_contacts_mean": mean_of(row["pos6_sidechain_contacts"]),
                "asp6_min_donor_distance_mean": mean_of(row["asp6_min_donor_distance"]),
                "asp6_buried_frac_mean": mean_of(row["asp6_buried_frac"]),
                "net_interface_charge_mean": mean_of(row["net_interface_charge"]),
                "iptm_g12d_mean": mean_of(row["iptm_g12d"]),
            })
    print(json.dumps(summary, indent=2))
    for row in survivors:
        print(f"#{row['rank']} {row['backbone']} ({row['variant']}) "
              f"cluster={row['cluster']} passes={row['gate_passes']}/{row['draws']} "
              f"pep={mean_of(row['peptide_buried_frac'])} "
              f"pos6_sc={mean_of(row['pos6_sidechain_contacts'])} "
              f"donor={mean_of(row['asp6_min_donor_distance'])} "
              f"asp6_buried={mean_of(row['asp6_buried_frac'])} "
              f"charge={mean_of(row['net_interface_charge'])} "
              f"rmsd={row['binder_ca_rmsd_mean_replicates']}/"
              f"{row['binder_ca_rmsd_max_replicates']} A "
              f"iptm={mean_of(row['iptm_g12d'])}")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--stage", choices=("screen", "final"), default="screen")
    args = parser.parse_args()
    if args.stage == "screen":
        rank_batch04.OUT = OUT
        run_batch04.OUT = OUT
        return rank_batch04.main()
    return final()


if __name__ == "__main__":
    raise SystemExit(main())
