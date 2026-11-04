#!/usr/bin/env python3
"""Batch-09 matched regional-anchor counter-screen.

Why this exists
---------------
The historical off-target panel is asymmetric:

* G12D / WT / decoy are constrained to peptide chain C position 6.
* apo has no peptide chain C, so it is unconstrained.

That makes the apo arm scientifically useful as a warning but not a
like-for-like comparison. Batch-09 fixes only that asymmetry.

All four arms use the SAME HLA-groove regional anchor:

    HLA-A*11:01 residues A69 and A155, max_distance 10 A.

Those residues were fixed before any batch-09 result was seen. They are on the
two groove helices around the central peptide-facing surface and are contacted
by the experimentally characterised JDIa41b1 TCR (Poole et al. 2022). The
constraint is deliberately regional rather than mutation-specific: it asks the
sampler to place the binder over the same groove patch in every state.

This script DOES NOT replace the direct Asp6 experiment and changes no
historical gate. It answers one narrow question:

    When the docking prior is matched, does peptide-free HLA still attract the
    existing lead binders more strongly / bury more surface than loaded pHLA?

Subjects are fixed from batch-09's control audit:
    bb50-g1, bb54-g2, bb25

One Boltz-2 call per design x arm with diffusion_samples=3, so the default run
is 12 calls / 36 internally distinct diffusion trajectories.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
import time
from pathlib import Path
from typing import Any, Dict, List

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(REPO / "scripts"))

import batch07_settings as b7  # noqa: E402
import run_pipeline  # noqa: E402
import sampling_study as study  # noqa: E402
import selectivity_panel as panel  # noqa: E402
from nvidia_client import BOLTZ2_PATH, NvidiaClient  # noqa: E402

OUT = REPO / "results" / "batch-09"

ARMS = ("g12d", "wt", "apo", "decoy")
SUBJECT_KEYS = ("bb50-g1", "bb54-g2", "bb25")

# Fixed from the positive-control literature before any batch-09 result.
# A69 is on the alpha-1 helix; A155 is on the alpha-2 helix. The two points
# flank the central peptide-facing groove surface contacted by the known TCR.
MATCHED_GROOVE_CONTACTS = (
    {"id": run_pipeline.HLA_CHAIN, "residue_index": 69},
    {"id": run_pipeline.HLA_CHAIN, "residue_index": 155},
)
MATCHED_GROOVE_MAX_DISTANCE = 10.0

# Preserve each design's previously identified sampling regime, but use it
# identically across all four arms for that design.
SETTINGS = {
    "bb50-g1": (8, 50, False),   # recyc8
    "bb25": (8, 50, False),      # recyc8
    "bb54-g2": (3, 50, True),    # r3np / nopot
}


def matched_constraint() -> List[Dict[str, Any]]:
    """Return a fresh copy of the fixed HLA regional anchor."""
    return [{
        "binder": run_pipeline.BINDER_CHAIN,
        "contacts": [dict(x) for x in MATCHED_GROOVE_CONTACTS],
        "max_distance": MATCHED_GROOVE_MAX_DISTANCE,
    }]


def load_subjects() -> Dict[str, Dict[str, Any]]:
    """Load the exact historical sequences; generate nothing new."""
    subjects = {s["key"]: s for s in study.load_subjects()
                if s["key"] in {"bb50-g1", "bb25"}}
    bb54 = b7.subject_design()
    subjects["bb54-g2"] = {
        "key": "bb54-g2",
        "sequence": bb54["sequence"],
        "backbone": bb54.get("backbone", ""),
        "variant": bb54.get("variant", ""),
        "binder_length": bb54.get("binder_length", 46),
        "source_name": bb54["name"],
    }
    missing = [k for k in SUBJECT_KEYS if k not in subjects]
    if missing:
        raise RuntimeError(f"could not load batch-09 subject(s): {missing}")
    return subjects


def call_boltz(client: NvidiaClient, polymers, design_key: str, samples: int):
    recyc, sampling_steps, nopot = SETTINGS[design_key]
    body = {
        "polymers": [dict(p) for p in polymers],
        "recycling_steps": recyc,
        "sampling_steps": sampling_steps,
        "diffusion_samples": int(samples),
        "step_scale": 1.638,
        "without_potentials": nopot,
        # IMPORTANT: unlike every earlier apo arm, this is present for apo too.
        "constraints": matched_constraint(),
    }
    return client.post(BOLTZ2_PATH, body), body


def checkpoint_path() -> Path:
    return OUT / "matched_anchor_checkpoint.jsonl"


def append_checkpoint(record: Dict[str, Any]) -> None:
    with checkpoint_path().open("a") as handle:
        handle.write(json.dumps(record, default=str) + "\n")
        handle.flush()
        os.fsync(handle.fileno())


def done_call_keys() -> set:
    """A call is complete when draw 0 has already been checkpointed."""
    if not checkpoint_path().exists():
        return set()
    out = set()
    for line in checkpoint_path().read_text().splitlines():
        if not line.strip():
            continue
        row = json.loads(line)
        if row.get("draw") == 0:
            out.add((row["design"], row["arm"]))
    return out


def sample_payload(payload: Dict[str, Any], index: int) -> Dict[str, Any]:
    def one(key):
        values = payload.get(key) or []
        return [values[index]] if index < len(values) else []
    return {
        "structures": [payload["structures"][index]],
        "pair_chains_iptm_scores": one("pair_chains_iptm_scores"),
        "iptm_scores": one("iptm_scores"),
        "complex_ipde_scores": one("complex_ipde_scores"),
        "complex_plddt_scores": one("complex_plddt_scores"),
    }


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--designs", default=",".join(SUBJECT_KEYS))
    ap.add_argument("--arms", default=",".join(ARMS))
    ap.add_argument("--draws", type=int, default=3)
    ap.add_argument("--max-boltz-calls", type=int, default=12)
    ap.add_argument("--timeout", type=float, default=900.0)
    ap.add_argument("--max-retries", type=int, default=8)
    args = ap.parse_args(argv)

    if args.draws < 1:
        raise ValueError("--draws must be >= 1")

    wanted = [x.strip() for x in args.designs.split(",") if x.strip()]
    arms = [x.strip() for x in args.arms.split(",") if x.strip()]
    bad_designs = sorted(set(wanted) - set(SUBJECT_KEYS))
    bad_arms = sorted(set(arms) - set(ARMS))
    if bad_designs:
        raise ValueError(f"unknown design(s): {bad_designs}")
    if bad_arms:
        raise ValueError(f"unknown arm(s): {bad_arms}")

    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "structures").mkdir(parents=True, exist_ok=True)

    subjects = load_subjects()
    states = run_pipeline.load_target_states(REPO / "data" / "targets")
    client = NvidiaClient(timeout=args.timeout, max_retries=args.max_retries)
    done = done_call_keys()

    calls = 0
    started = time.time()
    for key in wanted:
        subject = subjects[key]
        sequence = subject["sequence"]
        for arm in arms:
            if (key, arm) in done:
                print(f"[skip] {key} {arm} already complete")
                continue
            if calls >= args.max_boltz_calls:
                print(f"[stop] call cap {args.max_boltz_calls} reached")
                return 0

            polymers = list(states[arm]) + [{
                "id": run_pipeline.BINDER_CHAIN,
                "molecule_type": "protein",
                "sequence": sequence,
            }]
            payload, body = call_boltz(client, polymers, key, args.draws)
            calls += 1
            structures = payload.get("structures") or []
            hashes = [hashlib.sha1(s["structure"].encode()).hexdigest()
                      for s in structures]

            for draw, structure in enumerate(structures):
                text = structure["structure"]
                stem = f"{key}_{arm}_matched_r{draw}"
                (OUT / "structures" / f"{stem}.cif").write_text(text)

                shim = sample_payload(payload, draw)
                row: Dict[str, Any] = {
                    "iptm": run_pipeline._round(
                        run_pipeline.binder_interface_iptm(
                            shim, len(polymers) - 1), 4),
                    "pde_interface": run_pipeline._round(
                        run_pipeline._first(shim, "complex_ipde_scores"), 2),
                }
                try:
                    row.update(panel.panel_metrics(text, sequence, arm))
                    row["analysis_error"] = ""
                except Exception as exc:  # noqa: BLE001
                    row["analysis_error"] = f"{type(exc).__name__}: {exc}"

                append_checkpoint({
                    "design": key,
                    "source_name": subject.get("source_name", key),
                    "arm": arm,
                    "draw": draw,
                    "sequence": sequence,
                    "sampling": {
                        "recycling_steps": body["recycling_steps"],
                        "sampling_steps": body["sampling_steps"],
                        "without_potentials": body["without_potentials"],
                        "diffusion_samples": body["diffusion_samples"],
                        "step_scale": body["step_scale"],
                    },
                    "constraint": body["constraints"],
                    "constraint_kind": "matched_hla_regional_anchor",
                    "structure_sha1": hashes[draw],
                    "distinct_structures_in_call": len(set(hashes)),
                    "row": row,
                    "structure_file": f"structures/{stem}.cif",
                    "scored_at": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
                })

            print(
                f"[boltz2 {calls}] {key} {arm}: "
                f"{len(structures)} draws, distinct={len(set(hashes))} "
                f"({time.time() - started:.0f}s)",
                flush=True,
            )

    meta = {
        "purpose": "matched regional-anchor counter-screen; diagnostic only",
        "designs": wanted,
        "arms": arms,
        "draws_per_call": args.draws,
        "boltz2_calls_this_run": calls,
        "matched_constraint": matched_constraint(),
        "settings": {k: SETTINGS[k] for k in wanted},
        "historical_gates_changed": False,
        "track1_slots_spent": 0,
        "elapsed_seconds": round(time.time() - started, 1),
    }
    (OUT / "matched_anchor_run.json").write_text(
        json.dumps(meta, indent=2, default=str)
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
