#!/usr/bin/env python3
"""Is a batch-03 specificity gap reproducible? Re-run the top three.

The 7OW6 control showed that Boltz-2's ipTM is reproducible to ~0.001 on a
*confident* complex but cannot tell G12D from WT for a known selective binder.
This asks the companion question in the design regime: re-run the full
four-state constrained panel for the three best batch-03 designs and compare
the gap with the one already recorded in
`results/constrained-rescore/checkpoint.jsonl`.
"""
from __future__ import annotations
import json, sys, time
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))
import run_pipeline                                   # noqa: E402
from nvidia_client import NvidiaClient, BOLTZ2_PATH   # noqa: E402

TARGETS = ("S1-pos6-68_bb18", "S1-pos6-68_bb19", "S1-pos6-68_bb13")
OUT = REPO / "results" / "constrained-rescore"
CKPT = OUT / "replicate_top3.jsonl"

def main() -> int:
    first = {}
    for line in (OUT / "checkpoint.jsonl").open():
        r = json.loads(line)
        first[r["backbone"]] = r
    states = run_pipeline.load_target_states(REPO / "data" / "targets")
    contacts = [{"id": "C", "residue_index": i} for i in range(4, 9)]
    pocket = [{"binder": "D", "contacts": contacts, "max_distance": 6.0}]
    client = NvidiaClient()
    done = set()
    if CKPT.exists():
        done = {json.loads(l)["backbone"] for l in CKPT.open()}
    for name in TARGETS:
        if name in done:
            continue
        seq = first[name]["sequence"]
        rec = {"backbone": name, "first_pass_gap": first[name].get("specificity_gap")}
        for state in ("g12d", "wt", "apo", "decoy"):
            polymers = [dict(p) for p in states[state]] + [
                {"id": "D", "molecule_type": "protein", "sequence": seq}]
            body = {"polymers": polymers, "recycling_steps": 3, "sampling_steps": 50,
                    "diffusion_samples": 1, "step_scale": 1.638,
                    "without_potentials": False}
            if state != "apo":
                body["constraints"] = pocket
            payload = client.post(BOLTZ2_PATH, body)
            rec[f"iptm_{state}"] = run_pipeline.binder_interface_iptm(
                payload, len(polymers) - 1)
            time.sleep(8)
        offs = [rec[f"iptm_{s}"] for s in ("wt", "apo", "decoy")
                if rec.get(f"iptm_{s}") is not None]
        rec["replicate_gap"] = round(rec["iptm_g12d"] - max(offs), 4)
        with CKPT.open("a") as fh:
            fh.write(json.dumps(rec, default=str) + "\n")
        print(f"[rep3] {name} first={rec['first_pass_gap']} "
              f"replicate={rec['replicate_gap']} g12d={rec['iptm_g12d']:.3f}", flush=True)
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
