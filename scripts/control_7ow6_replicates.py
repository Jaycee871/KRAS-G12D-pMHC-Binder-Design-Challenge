#!/usr/bin/env python3
"""How much of the 7OW6 control's specificity gap is noise?

The single-sample control gave ipTM(G12D) 0.927 vs ipTM(WT) 0.931 - a gap of
-0.004 for a *known selective* binder. This script measures the run-to-run and
sample-to-sample spread of that number so the gap can be compared with its own
noise floor: two independent Boltz-2 calls per state, five diffusion samples
each, same pocket constraint, G12D / WT / decoy.
"""
from __future__ import annotations
import json, sys, time, statistics
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))
import complex_analysis as ca            # noqa: E402
import run_pipeline                      # noqa: E402
from nvidia_client import NvidiaClient, BOLTZ2_PATH   # noqa: E402
from control_7ow6_constrained import tcr_sequences, binder_vs_target_iptm  # noqa: E402

OUT = REPO / "results" / "validation"
CKPT = OUT / "7ow6_replicates.jsonl"
REPLICATES = 2
SAMPLES = 5

def main() -> int:
    alpha, beta = tcr_sequences(REPO / "data" / "validation" / "7ow6.cif")
    states = run_pipeline.load_target_states(REPO / "data" / "targets")
    contacts = [{"id": "C", "residue_index": i} for i in range(4, 9)]
    pocket = [{"binder": "D", "contacts": contacts, "max_distance": 6.0},
              {"binder": "E", "contacts": contacts, "max_distance": 6.0}]
    client = NvidiaClient()
    done = set()
    if CKPT.exists():
        done = {(json.loads(l)["state"], json.loads(l)["replicate"]) for l in CKPT.open()}
    series = {}
    for state in ("g12d", "wt", "decoy"):
        for rep in range(REPLICATES):
            if (state, rep) in done:
                continue
            polymers = [dict(p) for p in states[state]] + [
                {"id": "D", "molecule_type": "protein", "sequence": alpha},
                {"id": "E", "molecule_type": "protein", "sequence": beta}]
            body = {"polymers": polymers, "recycling_steps": 3, "sampling_steps": 50,
                    "diffusion_samples": SAMPLES, "step_scale": 1.638,
                    "without_potentials": False, "constraints": pocket}
            t0 = time.time()
            payload = client.post(BOLTZ2_PATH, body)
            bi = {len(polymers) - 2, len(polymers) - 1}
            n = len(payload.get("pair_chains_iptm_scores") or [])
            values = [binder_vs_target_iptm(payload, bi, sample=i) for i in range(n)]
            # ipTM of each binder chain against the *peptide* chain alone (index 2)
            pep = []
            for i in range(n):
                m = payload["pair_chains_iptm_scores"][i]
                pep.append(max(float(m[str(b)]["2"]) for b in bi))
            rec = {"state": state, "replicate": rep, "iptm_vs_target": values,
                   "iptm_vs_peptide": pep,
                   "iptm_scores": payload.get("iptm_scores"),
                   "seconds": round(time.time() - t0, 1)}
            with CKPT.open("a") as fh:
                fh.write(json.dumps(rec, default=str) + "\n")
            print(f"[rep] {state} r{rep} {values} {rec['seconds']}s", flush=True)
            time.sleep(8)
    for line in CKPT.open():
        r = json.loads(line)
        series.setdefault(r["state"], {"t": [], "p": []})
        series[r["state"]]["t"] += r["iptm_vs_target"]
        series[r["state"]]["p"] += r["iptm_vs_peptide"]
    summary = {}
    for state, v in series.items():
        summary[state] = {
            "n": len(v["t"]),
            "iptm_vs_target_mean": round(statistics.mean(v["t"]), 4),
            "iptm_vs_target_sd": round(statistics.pstdev(v["t"]), 4),
            "iptm_vs_target_min": round(min(v["t"]), 4),
            "iptm_vs_target_max": round(max(v["t"]), 4),
            "iptm_vs_peptide_mean": round(statistics.mean(v["p"]), 4),
            "iptm_vs_peptide_sd": round(statistics.pstdev(v["p"]), 4),
        }
    (OUT / "7ow6_replicates_summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
