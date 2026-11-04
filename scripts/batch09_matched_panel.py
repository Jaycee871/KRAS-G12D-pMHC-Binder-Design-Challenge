#!/usr/bin/env python3
"""Batch-09: matched groove-centred counter-screen before any new backbone batch.

Why this experiment exists
--------------------------
All peptide-containing panel arms in batches 06-08 were predicted with the
same peptide-C6 pocket constraint, while peptide-free HLA (apo) necessarily ran
unconstrained.  The resulting apo burial was extremely large, but that
comparison is not like-for-like: the docking prior changes together with the
biological state.

Challenge 5 itself identifies HLA-A*11:01 Gln70 and Arg114 as the groove
residues flanking KRAS peptide position 6 in 7OW6.  Those two HLA positions are
present in G12D, WT, apo and decoy targets, so batch-09 uses them as one
*identical* docking prior in every arm:

    binder D -> HLA A70/A114, max_distance 6 A

This does not weaken the apo objective and does not change any existing gate.
It is a diagnostic asking whether the severe apo signal survives when the
counter-screen is matched to the same physical region of the groove.

No RFdiffusion or ProteinMPNN generation occurs here.  The candidate set is
fixed before model calls and reuses informative designs already measured by
the project.

Usage
-----
    python scripts/batch09_matched_panel.py --dry-run
    python scripts/batch09_matched_panel.py
    python scripts/batch09_matched_panel.py --designs bb54-g2,bb135-g1

The live path requires NVIDIA_API_KEY through src/nvidia_client.py.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import statistics
import sys
import time
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Sequence

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(REPO / "scripts"))

import batch06_panel as b6  # noqa: E402
import batch07_settings as b7s  # noqa: E402
import run_pipeline  # noqa: E402
import sampling_study as study  # noqa: E402
import selectivity_panel as panel  # noqa: E402
from nvidia_client import BOLTZ2_PATH, NvidiaClient  # noqa: E402

OUT = REPO / "results" / "batch-09"

# Same physical pocket in all four states.  The 6 A value is deliberately the
# same distance already used by the project's validated peptide-C6 constraint;
# only the invariant anchor residues change.
GROOVE_POCKET = [{
    "binder": run_pipeline.BINDER_CHAIN,
    "contacts": [
        {"id": run_pipeline.HLA_CHAIN, "residue_index": 70},
        {"id": run_pipeline.HLA_CHAIN, "residue_index": 114},
    ],
    "max_distance": 6.0,
}]

CONDITION = "nopot"
RECYCLING_STEPS, SAMPLING_STEPS, WITHOUT_POTENTIALS = study.CONDITIONS[CONDITION]
ARMS = ("g12d", "wt", "apo", "decoy")

# Fixed before any batch-09 calls.  Four historical mechanisms plus one
# peptide-centred batch-08 control; none was selected from batch-09 outcomes.
DEFAULT_DESIGNS = (
    "bb54-g2",       # strongest pooled WT evidence / current parent
    "bb135-g1",      # tight Arg-Asp readout; 2/3 at baseline, 3/3 in nopot
    "bb173-g1",      # high Asp6 burial; intermittent second pose
    "bb90-g2",       # tightest pose spread among batch-05 screening passers
    "b8-apex46-g2dw",# peptide-centred batch-08 apo diagnostic control
)


def _batch08_control() -> Dict[str, Any]:
    data = json.loads((REPO / "results" / "batch-08" / "apo_probe.json").read_text())
    name = "B8-apex-46_bb46_donor-anchored-g2-dw"
    row = next(d for d in data["designs"] if d["design"] == name)
    return {
        "key": "b8-apex46-g2dw",
        "name": name,
        "sequence": row["sequence"],
        "source": "batch-08/apo_probe",
        "reason": "peptide-centred control from the predeclared batch-08 apo probe",
    }


def candidate_registry() -> Dict[str, Dict[str, Any]]:
    """Informative existing designs; deterministic and independent of batch-09."""
    historical = {s["key"]: s for s in study.load_subjects()}

    bb54 = b7s.subject_design()
    out: Dict[str, Dict[str, Any]] = {
        "bb54-g2": {
            "key": "bb54-g2",
            "name": bb54["name"],
            "sequence": bb54["sequence"],
            "source": "batch-06/07/08 parent",
            "reason": "current parent; reproducible Asp6 reader under r3np",
        },
        "bb135-g1": {
            "key": "bb135-g1",
            "name": "bb135-g1",
            "sequence": historical["bb135-g1"]["sequence"],
            "source": "batch-05/sampling-study",
            "reason": "tight Arg-Asp readout; nopot stabilised 3/3",
        },
        "bb173-g1": {
            "key": "bb173-g1",
            "name": "bb173-g1",
            "sequence": historical["bb173-g1"]["sequence"],
            "source": "batch-05/sampling-study",
            "reason": "high Asp6 burial with intermittent alternate pose",
        },
        "bb90-g2": {
            "key": "bb90-g2",
            "name": "bb90-g2",
            "sequence": historical["bb90-g2"]["sequence"],
            "source": "batch-05/sampling-study",
            "reason": "tightest historical pose spread among the batch-05 gate screeners",
        },
    }
    b8 = _batch08_control()
    out[b8["key"]] = b8
    return out


def selected_candidates(keys: Sequence[str]) -> List[Dict[str, Any]]:
    registry = candidate_registry()
    missing = [k for k in keys if k not in registry]
    if missing:
        raise ValueError(
            "unknown design key(s): " + ", ".join(missing)
            + "; choices: " + ", ".join(registry)
        )
    rows = [registry[k] for k in keys]
    if len({r["sequence"] for r in rows}) != len(rows):
        raise ValueError("candidate panel contains duplicate sequences")
    return rows


def build_body(polymers: Sequence[Mapping[str, Any]], draws: int) -> Dict[str, Any]:
    """One matched Boltz-2 request body, identical for every biological arm."""
    return {
        "polymers": [dict(p) for p in polymers],
        "recycling_steps": RECYCLING_STEPS,
        "sampling_steps": SAMPLING_STEPS,
        "diffusion_samples": int(draws),
        "step_scale": 1.638,
        "without_potentials": WITHOUT_POTENTIALS,
        "constraints": GROOVE_POCKET,
    }


def target_anchor_residues(states: Mapping[str, Sequence[Mapping[str, str]]]) -> Dict[str, str]:
    """Return the HLA amino acids at the two invariant groove anchors."""
    found: Dict[str, str] = {}
    for arm in ARMS:
        hla = next(p for p in states[arm] if p["id"] == run_pipeline.HLA_CHAIN)
        seq = str(hla["sequence"])
        found[arm] = seq[69] + seq[113]  # 1-based positions 70 and 114
    return found


def validate_matched_target(states: Mapping[str, Sequence[Mapping[str, str]]]) -> None:
    anchors = target_anchor_residues(states)
    bad = {arm: aa for arm, aa in anchors.items() if aa != "QR"}
    if bad:
        raise ValueError(
            "HLA groove anchors do not match 7OW6 Q70/R114: "
            + json.dumps(bad, sort_keys=True)
        )


def checkpoint_path() -> Path:
    return OUT / "matched_checkpoint.jsonl"


def done_keys() -> set:
    path = checkpoint_path()
    if not path.exists():
        return set()
    out = set()
    for line in path.read_text().splitlines():
        if line.strip():
            row = json.loads(line)
            out.add((row["design_key"], row["arm"], row["draw"]))
    return out


def append_checkpoint(record: Mapping[str, Any]) -> None:
    with checkpoint_path().open("a") as fh:
        fh.write(json.dumps(dict(record), default=str) + "\n")
        fh.flush()
        os.fsync(fh.fileno())


def load_checkpoint() -> List[Dict[str, Any]]:
    path = checkpoint_path()
    if not path.exists():
        return []
    return [json.loads(line) for line in path.read_text().splitlines() if line.strip()]


def sample_payload(payload: Mapping[str, Any], index: int) -> Dict[str, Any]:
    def pick(key: str):
        values = payload.get(key) or []
        return [values[index]] if index < len(values) else []

    return {
        "structures": [payload["structures"][index]],
        "pair_chains_iptm_scores": pick("pair_chains_iptm_scores"),
        "iptm_scores": pick("iptm_scores"),
        "complex_ipde_scores": pick("complex_ipde_scores"),
    }


def _hla_buried(metrics: Mapping[str, Any]) -> float | None:
    by_chain = metrics.get("buried_area_by_chain") or {}
    hla = by_chain.get(run_pipeline.HLA_CHAIN)
    return None if hla is None else float(hla)


def score_structure(
    text: str,
    sequence: str,
    arm: str,
    shim: Mapping[str, Any],
    binder_index: int,
) -> Dict[str, Any]:
    metrics = panel.panel_metrics(text, sequence, arm)
    metrics["buried_area_hla"] = _hla_buried(metrics)
    metrics["iptm"] = run_pipeline._round(
        run_pipeline.binder_interface_iptm(shim, binder_index), 4)
    metrics["pde_interface"] = run_pipeline._round(
        run_pipeline._first(shim, "complex_ipde_scores"), 2)
    total = metrics.get("buried_area_total")
    pep = metrics.get("buried_area_peptide")
    metrics["peptide_area_share"] = (
        round(float(pep) / float(total), 3)
        if pep is not None and total not in (None, 0, 0.0) else None
    )
    return metrics


def run_arm(
    client,
    design: Mapping[str, Any],
    arm: str,
    states: Mapping[str, Sequence[Mapping[str, str]]],
    draws: int,
    already: set,
) -> int:
    if (design["key"], arm, 0) in already:
        print(f"[skip] {design['key']} {arm} already checkpointed")
        return 0

    polymers = list(states[arm]) + [{
        "id": run_pipeline.BINDER_CHAIN,
        "molecule_type": "protein",
        "sequence": design["sequence"],
    }]
    body = build_body(polymers, draws)
    payload = client.post(BOLTZ2_PATH, body)
    structures = payload.get("structures") or []
    if len(structures) < draws:
        raise RuntimeError(
            f"{design['key']} {arm}: requested {draws} structures, got {len(structures)}"
        )

    binder_index = len(polymers) - 1
    hashes = [hashlib.sha1(s["structure"].encode()).hexdigest() for s in structures]
    for index, structure in enumerate(structures[:draws]):
        text = structure["structure"]
        stem = f"{design['key']}_{arm}_{CONDITION}_r{index}"
        (OUT / "structures" / f"{stem}.cif").write_text(text)
        shim = sample_payload(payload, index)
        try:
            row = score_structure(text, design["sequence"], arm, shim, binder_index)
            error = ""
        except Exception as exc:  # keep the run checkpointable
            row = {}
            error = f"{type(exc).__name__}: {exc}"

        append_checkpoint({
            "design_key": design["key"],
            "design_name": design["name"],
            "source": design["source"],
            "sequence": design["sequence"],
            "arm": arm,
            "condition": CONDITION,
            "draw": index,
            "constraint": GROOVE_POCKET,
            "boltz_kwargs": {
                "recycling_steps": RECYCLING_STEPS,
                "sampling_steps": SAMPLING_STEPS,
                "diffusion_samples": draws,
                "without_potentials": WITHOUT_POTENTIALS,
                "step_scale": 1.638,
            },
            "structure_sha1": hashes[index],
            "distinct_structures_in_call": len(set(hashes)),
            "row": row,
            "analysis_error": error,
            "structure_file": f"structures/{stem}.cif",
            "scored_at": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        })
        print(
            f"  [{stem}] ipTM={row.get('iptm')} total={row.get('buried_area_total')} "
            f"hla={row.get('buried_area_hla')} pep={row.get('buried_area_peptide')} "
            f"pos6={row.get('pos6_contacts')} gate={row.get('asp6_gate_pass')} "
            f"error={error or '-'}",
            flush=True,
        )
    return 1


def _rows(records: Iterable[Mapping[str, Any]], key: str, arm: str) -> List[Dict[str, Any]]:
    return [
        dict(r.get("row") or {})
        for r in records
        if r.get("design_key") == key and r.get("arm") == arm
        and not r.get("analysis_error")
    ]


def _median(rows: Sequence[Mapping[str, Any]], key: str) -> float | None:
    values = [float(r[key]) for r in rows if r.get(key) not in (None, "")]
    return statistics.median(values) if values else None


def summarise(keys: Sequence[str]) -> Dict[str, Any]:
    records = load_checkpoint()
    designs = []
    for key in keys:
        arms = {arm: _rows(records, key, arm) for arm in ARMS}
        g12d, wt, apo, decoy = (arms[a] for a in ARMS)
        wt_result = b6.wt_verdict(g12d, wt) if g12d and wt else None

        g_total = _median(g12d, "buried_area_total")
        a_total = _median(apo, "buried_area_total")
        d_total = _median(decoy, "buried_area_total")
        designs.append({
            "design_key": key,
            "draws": {arm: len(arms[arm]) for arm in ARMS},
            "g12d_asp6_gate_passes": sum(bool(r.get("asp6_gate_pass")) for r in g12d),
            "wt": wt_result,
            "median_buried_area": {
                arm: _median(arms[arm], "buried_area_total") for arm in ARMS
            },
            "median_hla_buried_area": {
                arm: _median(arms[arm], "buried_area_hla") for arm in ARMS
            },
            "apo_to_g12d_total_area_ratio": (
                round(a_total / g_total, 3)
                if a_total is not None and g_total not in (None, 0.0) else None
            ),
            "decoy_to_g12d_total_area_ratio": (
                round(d_total / g_total, 3)
                if d_total is not None and g_total not in (None, 0.0) else None
            ),
            "note": (
                "apo/decoy area ratios are descriptive diagnostics only; "
                "batch-09 pre-registers no new pass threshold"
            ),
        })

    return {
        "experiment": "batch-09 matched HLA-groove counter-screen",
        "condition": CONDITION,
        "matched_constraint": GROOVE_POCKET,
        "invariant_hla_anchors": {"70": "Q", "114": "R"},
        "arms": list(ARMS),
        "candidate_keys_predeclared": list(keys),
        "existing_wt_criterion_reused_unchanged": b6.WT_CRITERION,
        "new_gate_or_threshold_added": False,
        "designs": designs,
    }


class DryRunClient:
    """Record request bodies without model calls; useful for pre-flight checks."""

    def __init__(self) -> None:
        self.bodies: List[Dict[str, Any]] = []

    def post(self, path: str, body: Mapping[str, Any]) -> Dict[str, Any]:
        self.bodies.append({"path": path, "body": dict(body)})
        return {"structures": []}


def dry_run(designs, states, draws: int) -> int:
    validate_matched_target(states)
    expected = json.dumps(GROOVE_POCKET, sort_keys=True)
    requests = []
    for design in designs:
        for arm in ARMS:
            polymers = list(states[arm]) + [{
                "id": run_pipeline.BINDER_CHAIN,
                "molecule_type": "protein",
                "sequence": design["sequence"],
            }]
            body = build_body(polymers, draws)
            assert json.dumps(body["constraints"], sort_keys=True) == expected
            requests.append({
                "design": design["key"],
                "arm": arm,
                "chains": [p["id"] for p in body["polymers"]],
                "constraint": body["constraints"],
            })
    print(json.dumps({
        "dry_run": True,
        "anchors": target_anchor_residues(states),
        "condition": {
            "recycling_steps": RECYCLING_STEPS,
            "sampling_steps": SAMPLING_STEPS,
            "without_potentials": WITHOUT_POTENTIALS,
            "diffusion_samples": draws,
        },
        "n_model_calls": len(requests),
        "requests": requests,
    }, indent=2))
    return 0


def main(argv=None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--designs", default=",".join(DEFAULT_DESIGNS))
    parser.add_argument("--draws", type=int, default=3)
    parser.add_argument("--timeout", type=float, default=900.0)
    parser.add_argument("--max-boltz-calls", type=int, default=20)
    parser.add_argument("--max-retries", type=int, default=8)
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args(argv)

    if args.draws < 1:
        parser.error("--draws must be >= 1")

    keys = [k.strip() for k in args.designs.split(",") if k.strip()]
    designs = selected_candidates(keys)
    states = run_pipeline.load_target_states(REPO / "data" / "targets")
    validate_matched_target(states)

    if args.dry_run:
        return dry_run(designs, states, args.draws)

    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "structures").mkdir(parents=True, exist_ok=True)
    already = done_keys()
    client = NvidiaClient(timeout=args.timeout, max_retries=args.max_retries)
    calls, started = 0, time.time()

    for design in designs:
        print(f"\n[batch-09] === {design['key']} :: {design['reason']} ===", flush=True)
        for arm in ARMS:
            if calls >= args.max_boltz_calls:
                print("[stop] model-call cap reached", flush=True)
                break
            calls += run_arm(client, design, arm, states, args.draws, already)
            already = done_keys()
        if calls >= args.max_boltz_calls:
            break

    result = summarise(keys)
    result["boltz2_calls_this_segment"] = calls
    result["elapsed_seconds"] = round(time.time() - started, 1)
    (OUT / "matched_summary.json").write_text(
        json.dumps(result, indent=2, default=str) + "\n"
    )
    print(f"\n[batch-09] {calls} model calls -> {OUT / 'matched_summary.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
