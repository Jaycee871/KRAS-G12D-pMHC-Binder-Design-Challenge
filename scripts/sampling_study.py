#!/usr/bin/env python3
"""Boltz-2 sampling-settings study: can a setting collapse the pose spread?

Batches 04 and 05 converged on one diagnosis: the chemistry is solved, the
*pose* is not. Every Asp6 gate failure among the near-misses coincides with the
binder flipping between two docked modes, and the pMHC reference frame is
stable (0.08-0.33 A), so the 4-20 A binder RMSD spreads are real.

This script re-scores five *existing* designs - it generates no new designs and
consumes zero Track-1 design slots - under alternative Boltz-2 sampling
settings and measures whether the replicate spread collapses:

    baseline  recycling_steps 3,  sampling_steps 50, without_potentials False
    nopot     recycling_steps 3,  sampling_steps 50, without_potentials True
    recyc8    recycling_steps 8,  sampling_steps 50, without_potentials False
    combo     recycling_steps 8,  sampling_steps 50, without_potentials True

Everything else is held identical to the batch-05 run, including the validated
flat pocket constraint
`{"binder":"D","contacts":[{"id":"C","residue_index":6}],"max_distance":6.0}`
(nested or [chain,index] pair forms return HTTP 422), the strictly serial
>= 8 s Boltz-2 pacing from `src/nvidia_client.py`, and the scoring/analysis
path: `scripts/run_batch04.score_candidate`, which is what batches 04 and 05
both called. The only thing overridden is the request body's sampling fields.

Baseline numbers are *not* re-run: the batch-04/batch-05 replicate draws are
the baseline, read from their `replicates.json` / `replicate_structures/`.

Checkpointed per (design, condition, draw) with an fsync, as the batch scripts
do, so an interruption loses nothing; re-running resumes.
"""
from __future__ import annotations

import argparse, json, os, sys, time
from pathlib import Path
from typing import Any, Dict, List

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(REPO / "scripts"))

import run_batch04 as batch  # noqa: E402
import run_pipeline  # noqa: E402
from nvidia_client import BOLTZ2_PATH, NvidiaClient  # noqa: E402

OUT = REPO / "results" / "sampling-study"

# (key, source batch dir, label used in the write-up)
SUBJECTS = [
    ("bb135-g1", "batch-05", "B5-c6-62_bb135"),
    ("bb173-g1", "batch-05", "B5-c6-62_bb173"),
    ("bb90-g2",  "batch-05", "B5-c6-62_bb90"),
    ("bb50-g1",  "batch-05", "B5-c6-62_bb50"),
    ("bb25",     "batch-04", "B4-c6-62_bb25"),
]

CONDITIONS = {
    # name: (recycling_steps, sampling_steps, without_potentials)
    "nopot":  (3, 50, True),
    "recyc8": (8, 50, False),
    "combo":  (8, 50, True),
}
BASELINE = (3, 50, False)


def load_subjects() -> List[Dict[str, Any]]:
    out = []
    for key, source, backbone in SUBJECTS:
        survivors = json.loads(
            (REPO / "results" / source / "gate_results.json").read_text())["survivors"]
        match = next(s for s in survivors if s["backbone"] == backbone)
        out.append({"key": key, "source": source, **match})
    return out


def patched_boltz2(client, polymers, args):
    """batch.boltz2_constrained with the sampling fields taken from `args`.

    Identical body otherwise - same polymers, same step_scale, same single
    diffusion sample, same validated POCKET constraint.
    """
    body = {"polymers": [dict(p) for p in polymers],
            "recycling_steps": int(args.recycling_steps),
            "sampling_steps": int(args.sampling_steps),
            "diffusion_samples": 1, "step_scale": 1.638,
            "without_potentials": bool(args.without_potentials),
            "constraints": batch.POCKET}
    return client.post(BOLTZ2_PATH, body)


batch.boltz2_constrained = patched_boltz2


def append_checkpoint(record: Dict[str, Any]) -> None:
    path = OUT / "checkpoint.jsonl"
    with path.open("a") as handle:
        handle.write(json.dumps(record, default=str) + "\n")
        handle.flush()
        os.fsync(handle.fileno())


def done_keys() -> set:
    path = OUT / "checkpoint.jsonl"
    if not path.exists():
        return set()
    return {(json.loads(l)["design"], json.loads(l)["condition"], json.loads(l)["draw"])
            for l in path.read_text().splitlines() if l.strip()}


def main(argv=None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--conditions", default="nopot,recyc8")
    parser.add_argument("--draws", type=int, default=3)
    parser.add_argument("--designs", default="",
                        help="comma-separated subject keys; default all five")
    parser.add_argument("--timeout", type=float, default=600.0)
    parser.add_argument("--max-boltz-calls", type=int, default=40)
    args = parser.parse_args(argv)

    (OUT / "structures").mkdir(parents=True, exist_ok=True)
    subjects = load_subjects()
    if args.designs:
        wanted = set(args.designs.split(","))
        subjects = [s for s in subjects if s["key"] in wanted]
    conditions = [c.strip() for c in args.conditions.split(",") if c.strip()]

    states = run_pipeline.load_target_states(REPO / "data" / "targets")
    client = NvidiaClient(timeout=args.timeout)
    already = done_keys()
    calls = 0
    started = time.time()

    for condition in conditions:
        recyc, samp, nopot = CONDITIONS[condition]
        for subject in subjects:
            for draw in range(args.draws):
                if (subject["key"], condition, draw) in already:
                    print(f"[skip] {subject['key']} {condition} r{draw} already done")
                    continue
                if calls >= args.max_boltz_calls:
                    print(f"[stop] call cap {args.max_boltz_calls} reached")
                    return 0
                call_args = argparse.Namespace(
                    recycling_steps=recyc, sampling_steps=samp,
                    without_potentials=nopot)
                candidate = {k: subject.get(k, "") for k in
                             ("sequence", "strategy", "backbone", "variant",
                              "binder_length", "mpnn_score", "donor_bias_score",
                              "design_pos6_contacts", "design_pos6_distance")}
                try:
                    result = batch.score_candidate(client, candidate, states, call_args)
                except Exception as exc:  # noqa: BLE001
                    print(f"[fail] {subject['key']} {condition} r{draw}: "
                          f"{type(exc).__name__}: {exc}", file=sys.stderr)
                    calls += 1
                    continue
                calls += 1
                stem = f"{subject['key']}_{condition}_r{draw}"
                (OUT / "structures" / f"{stem}.cif").write_text(result["structure"])
                row = result["row"]
                append_checkpoint({
                    "design": subject["key"], "condition": condition, "draw": draw,
                    "boltz_kwargs": {"recycling_steps": recyc,
                                     "sampling_steps": samp,
                                     "diffusion_samples": 1,
                                     "without_potentials": nopot},
                    "constraint": batch.POCKET,
                    "row": row, "structure_file": f"structures/{stem}.cif",
                    "scored_at": time.strftime("%Y-%m-%dT%H:%M:%S%z")})
                print(f"[boltz2 {calls}] {stem} GATE={row.get('gate_pass')} "
                      f"donor={row.get('asp6_min_donor_distance') or '-'} "
                      f"pos6_sc={row.get('pos6_sidechain_contacts')} "
                      f"pep={row.get('peptide_buried_frac')} "
                      f"asp6bf={row.get('asp6_buried_frac')} "
                      f"charge={row.get('net_interface_charge')} "
                      f"iptm={row.get('iptm_g12d')} "
                      f"({time.time() - started:.0f}s)", flush=True)

    print(f"[sampling-study] {calls} Boltz-2 calls in {time.time() - started:.0f}s")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
