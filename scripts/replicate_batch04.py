#!/usr/bin/env python3
"""Replicate the constrained G12D prediction for batch-04's gate survivors.

Any quantity used to rank designs gets >= 3 identical calls and a reported
spread (results/batch-04/SUMMARY.md, and recommendation 2 of
results/batch-04/SUMMARY.md's predecessor). This re-runs the *identical*
constrained Boltz-2 request N times per survivor and reports the spread of the
gate quantities, so a survivor that only passes on one sampling draw is
visible as such.
"""
from __future__ import annotations

import argparse, json, statistics, sys, time
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(REPO / "scripts"))

import run_batch04 as batch  # noqa: E402
import run_pipeline  # noqa: E402
from nvidia_client import NvidiaClient  # noqa: E402

OUT = REPO / "results" / "batch-04"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--replicates", type=int, default=3)
    parser.add_argument("--timeout", type=float, default=150.0)
    parser.add_argument("--recycling-steps", type=int, default=3)
    parser.add_argument("--sampling-steps", type=int, default=50)
    args = parser.parse_args()

    survivors = json.loads((OUT / "gate_results.json").read_text())["survivors"]
    states = run_pipeline.load_target_states(REPO / "data" / "targets")
    client = NvidiaClient(timeout=args.timeout)
    records = []
    for survivor in survivors:
        runs = []
        for replicate in range(args.replicates):
            candidate = {k: survivor.get(k, "") for k in
                         ("sequence", "strategy", "backbone", "variant",
                          "binder_length", "mpnn_score", "donor_bias_score",
                          "design_pos6_contacts", "design_pos6_distance")}
            try:
                result = batch.score_candidate(client, candidate, states, args)
            except Exception as exc:  # noqa: BLE001
                print(f"[replicate] {survivor['backbone']} r{replicate} failed: "
                      f"{type(exc).__name__}: {exc}", file=sys.stderr)
                continue
            row = result["row"]
            runs.append(row)
            (OUT / "replicate_structures").mkdir(exist_ok=True)
            (OUT / "replicate_structures" /
             f"{survivor['backbone']}_{survivor['variant']}_r{replicate}.cif"
             ).write_text(result["structure"])
            print(f"[replicate] {survivor['backbone']} r{replicate} "
                  f"gate={row['gate_pass']} donor={row['asp6_min_donor_distance'] or '-'} "
                  f"pep={row['peptide_buried_frac']} pos6_sc={row['pos6_sidechain_contacts']} "
                  f"charge={row['net_interface_charge']} iptm={row['iptm_g12d']}",
                  flush=True)
            time.sleep(1)

        def spread(key):
            values = [float(r[key]) for r in runs
                      if r.get(key) not in ("", None)]
            if not values:
                return None
            return {"values": values,
                    "mean": round(statistics.fmean(values), 3),
                    "sd": round(statistics.pstdev(values), 3) if len(values) > 1 else 0.0}

        records.append({
            "backbone": survivor["backbone"], "variant": survivor["variant"],
            "sequence": survivor["sequence"],
            "replicates": len(runs),
            "gate_pass_count": sum(1 for r in runs if r["gate_pass"]),
            "asp6_min_donor_distance": spread("asp6_min_donor_distance"),
            "peptide_buried_frac": spread("peptide_buried_frac"),
            "pos6_sidechain_contacts": spread("pos6_sidechain_contacts"),
            "asp6_buried_frac": spread("asp6_buried_frac"),
            "net_interface_charge": spread("net_interface_charge"),
            "iptm_g12d": spread("iptm_g12d"),
            "runs": runs,
        })
        (OUT / "replicates.json").write_text(json.dumps(records, indent=2,
                                                        default=str))
    print(json.dumps([{k: v for k, v in r.items() if k != "runs"}
                      for r in records], indent=2, default=str))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
