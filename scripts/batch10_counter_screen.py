#!/usr/bin/env python3
"""Batch-10: matched-anchor counter-screen of review ranks 4-11.

The project already has 11 direct-readout review candidates, but ranks 4-11
lack a matched apo counter-screen. Before generating any new protein, spend
information-efficient Boltz-2 calls on the official third objective.

Phase screen:
    G12D + apo, 3 diffusion samples per call, 8 designs = 16 calls.

Phase confirm:
    WT + decoy only for designs that keep >=2/3 corrected Asp6 readout on G12D
    AND have median apo total buried area < median G12D total buried area.

The phase-screen advance rule is a predeclared *triage heuristic*, not a claim
that buried area is an experimental binding-affinity metric.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import statistics
import sys
import time
from pathlib import Path
from typing import Any, Dict, Iterable, List

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(REPO / "scripts"))

import batch09_matched_anchor as matched  # noqa: E402
import run_pipeline  # noqa: E402
import selectivity_panel as panel  # noqa: E402
from nvidia_client import BOLTZ2_PATH, NvidiaClient  # noqa: E402

OUT = REPO / "results" / "batch-10"
SHORTLIST = REPO / "results" / "readout-correction" / "submission-review" / "review_shortlist.csv"
SCREEN_ARMS = ("g12d", "apo")
CONFIRM_ARMS = ("wt", "decoy")
ALL_ARMS = SCREEN_ARMS + CONFIRM_ARMS
DRAWS = 3

# Sampling regime is fixed before batch-10 results.
# Each design keeps the best previously observed *non-matched* G12D regime
# where one exists. B4 bb26 had no later settings sweep, so it keeps the
# original batch-04 3/50/potentials setting.
SETTINGS = {
    "KRAS_REVIEW_04_8827a211": (8, 50, False),  # bb36: recyc8, 2/3
    "KRAS_REVIEW_05_84aaeb67": (8, 50, False),  # bb173: recyc8, 2/3
    "KRAS_REVIEW_06_7cc39d25": (3, 50, False),  # bb54 graft p3: r3, 2/3
    "KRAS_REVIEW_07_610c10d2": (8, 50, False),  # bb90: recyc8, 1/3
    "KRAS_REVIEW_08_c4b05d07": (3, 50, False),  # bb26: batch-04 default
    "KRAS_REVIEW_09_9f4561b9": (8, 50, False),  # bb38: recyc8 standard
    "KRAS_REVIEW_10_08ffd5a8": (3, 50, False),  # bb54 graft p1-v3: r3
    "KRAS_REVIEW_11_71561413": (3, 50, False),  # bb54 graft p2: r3
}


def load_subjects() -> Dict[str, Dict[str, Any]]:
    with SHORTLIST.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    subjects = {}
    for row in rows:
        rank = int(row["evidence_rank"])
        if rank < 4 or rank > 11:
            continue
        key = row["name"]
        if key not in SETTINGS:
            raise ValueError(f"batch-10 settings missing for {key}")
        sequence = row["sequence"]
        digest = hashlib.sha256(sequence.encode()).hexdigest()
        if digest != row["sequence_sha256"]:
            raise ValueError(f"shortlist sequence hash mismatch for {key}")
        subjects[key] = {
            "key": key,
            "rank": rank,
            "sequence": sequence,
            "sequence_sha256": digest,
            "family": row["family"],
            "aliases": row["aliases"],
            "tier": row["tier"],
        }
    if set(subjects) != set(SETTINGS):
        raise ValueError("review ranks 4-11 changed; update batch-10 explicitly")
    return dict(sorted(subjects.items(), key=lambda kv: kv[1]["rank"]))


def request_body(polymers, design_key: str, samples: int = DRAWS) -> dict:
    recyc, sampling_steps, nopot = SETTINGS[design_key]
    return {
        "polymers": [dict(p) for p in polymers],
        "recycling_steps": recyc,
        "sampling_steps": sampling_steps,
        "diffusion_samples": int(samples),
        "step_scale": 1.638,
        "without_potentials": nopot,
        "constraints": matched.matched_constraint(),
    }


def checkpoint_path() -> Path:
    return OUT / "counter_screen_checkpoint.jsonl"


def checkpoint_records() -> List[dict]:
    path = checkpoint_path()
    if not path.exists():
        return []
    records, seen = [], set()
    for line_no, line in enumerate(path.read_text().splitlines(), 1):
        if not line.strip():
            continue
        try:
            record = json.loads(line)
            key = (record["design"], record["arm"], record["draw"])
        except (ValueError, KeyError, TypeError) as exc:
            raise ValueError(f"invalid checkpoint line {line_no}") from exc
        if key in seen:
            raise ValueError(f"duplicate checkpoint sample: {key}")
        seen.add(key)
        records.append(record)
    return records


def append_checkpoint(record: dict) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    path = checkpoint_path()
    with path.open("a") as handle:
        handle.write(json.dumps(record, default=str) + "\n")
        handle.flush()
        os.fsync(handle.fileno())


def validate_records(records: Iterable[dict], subjects: Dict[str, dict],
                     draws: int = DRAWS) -> None:
    for record in records:
        key, arm, draw = record["design"], record["arm"], record["draw"]
        if key not in subjects or arm not in ALL_ARMS:
            raise ValueError(f"unknown checkpoint design/arm: {key} {arm}")
        if type(draw) is not int or draw not in range(draws):
            raise ValueError(f"unexpected checkpoint draw: {key} {arm} {draw}")
        recyc, steps, nopot = SETTINGS[key]
        expected = {
            "recycling_steps": recyc,
            "sampling_steps": steps,
            "without_potentials": nopot,
            "diffusion_samples": draws,
            "step_scale": 1.638,
        }
        if (record.get("sampling") != expected
                or record.get("constraint") != matched.matched_constraint()
                or record.get("sequence") != subjects[key]["sequence"]):
            raise ValueError(f"checkpoint protocol mismatch: {key} {arm} {draw}")
        relative = Path(record.get("structure_file", ""))
        path = (OUT / relative).resolve()
        if not path.is_relative_to((OUT / "structures").resolve()):
            raise ValueError(f"invalid structure path: {relative}")
        if (not path.is_file()
                or hashlib.sha1(path.read_bytes()).hexdigest() != record.get("structure_sha1")):
            raise ValueError(f"missing or changed checkpoint structure: {relative}")


def completed_call_keys(records: Iterable[dict], draws: int = DRAWS) -> set:
    counts: Dict[tuple, set] = {}
    for record in records:
        counts.setdefault((record["design"], record["arm"]), set()).add(record["draw"])
    return {key for key, ds in counts.items() if ds == set(range(draws))}


def _median(values):
    values = [float(x) for x in values if isinstance(x, (int, float))]
    return statistics.median(values) if values else None


def screen_summary(records: Iterable[dict], subjects: Dict[str, dict]) -> Dict[str, dict]:
    records = list(records)
    out = {}
    for key, subject in subjects.items():
        g = [r for r in records if r["design"] == key and r["arm"] == "g12d"]
        a = [r for r in records if r["design"] == key and r["arm"] == "apo"]
        complete = len(g) == DRAWS and len(a) == DRAWS
        g_area = _median([r["row"].get("buried_area_total") for r in g])
        a_area = _median([r["row"].get("buried_area_total") for r in a])
        passes = sum(bool(r["row"].get("asp6_gate_pass")) for r in g)
        ratio = (a_area / g_area) if g_area and a_area is not None else None
        advance = bool(complete and passes >= 2 and ratio is not None and ratio < 1.0)
        out[key] = {
            "rank": subject["rank"],
            "family": subject["family"],
            "complete": complete,
            "g12d_Asp6_gate_passes": passes,
            "g12d_median_total_buried_area": g_area,
            "apo_median_total_buried_area": a_area,
            "apo_to_g12d_total_area_ratio": ratio,
            "advance_to_wt_decoy": advance,
            "advance_rule": "G12D direct Asp6 >=2/3 AND apo/G12D median total area <1.0",
        }
    return out


def selected_for_phase(phase: str, records: List[dict],
                       subjects: Dict[str, dict]) -> tuple[List[str], tuple[str, ...]]:
    if phase == "screen":
        return list(subjects), SCREEN_ARMS
    summary = screen_summary(records, subjects)
    incomplete = [k for k, v in summary.items() if not v["complete"]]
    if incomplete:
        raise ValueError(
            "phase confirm requires complete G12D+apo screen first: " + ", ".join(incomplete)
        )
    return [k for k, v in summary.items() if v["advance_to_wt_decoy"]], CONFIRM_ARMS


def atomic_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w") as handle:
        json.dump(value, handle, indent=2, default=str)
        handle.flush()
        os.fsync(handle.fileno())
    temporary.replace(path)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--phase", choices=("screen", "confirm"), default="screen")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--draws", type=int, default=DRAWS)
    ap.add_argument("--max-boltz-calls", type=int, default=16)
    ap.add_argument("--timeout", type=float, default=900.0)
    ap.add_argument("--max-retries", type=int, default=8)
    args = ap.parse_args(argv)
    if args.draws != DRAWS:
        raise ValueError("batch-10 protocol is frozen at --draws 3")

    subjects = load_subjects()
    states = run_pipeline.load_target_states(REPO / "data" / "targets")
    records = checkpoint_records()
    validate_records(records, subjects, args.draws)
    wanted, arms = selected_for_phase(args.phase, records, subjects)
    done = completed_call_keys(records, args.draws)

    planned = [(key, arm) for key in wanted for arm in arms if (key, arm) not in done]
    if args.dry_run:
        requests = []
        for key, arm in planned:
            polymers = list(states[arm]) + [{
                "id": run_pipeline.BINDER_CHAIN,
                "molecule_type": "protein",
                "sequence": subjects[key]["sequence"],
            }]
            requests.append({"design": key, "arm": arm,
                             "body": request_body(polymers, key, args.draws)})
        print(json.dumps({
            "dry_run": True,
            "phase": args.phase,
            "model_calls": 0,
            "planned_calls": len(requests),
            "advance_summary": screen_summary(records, subjects),
            "requests": requests,
        }, indent=2))
        return 0

    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "structures").mkdir(parents=True, exist_ok=True)
    (OUT / "responses").mkdir(parents=True, exist_ok=True)

    client = None
    calls = 0
    started = time.time()
    status = "running"
    try:
        for key in wanted:
            sequence = subjects[key]["sequence"]
            for arm in arms:
                if (key, arm) in done:
                    print(f"[skip] {key} {arm} complete")
                    continue
                polymers = list(states[arm]) + [{
                    "id": run_pipeline.BINDER_CHAIN,
                    "molecule_type": "protein",
                    "sequence": sequence,
                }]
                body = request_body(polymers, key, args.draws)
                cache = OUT / "responses" / f"{key}_{arm}.json"
                existing = {r["draw"] for r in records
                            if r["design"] == key and r["arm"] == arm}
                if cache.exists():
                    saved = json.loads(cache.read_text())
                    if saved.get("request") != body:
                        raise ValueError(f"cached request mismatch: {key} {arm}")
                    payload = saved["payload"]
                else:
                    if existing:
                        raise ValueError(f"partial checkpoint without response cache: {key} {arm}")
                    if calls >= args.max_boltz_calls:
                        status = "call_cap_reached"
                        break
                    if client is None:
                        client = NvidiaClient(timeout=args.timeout, max_retries=args.max_retries)
                    calls += 1
                    payload = client.post(BOLTZ2_PATH, body)
                    atomic_json(cache, {"request": body, "payload": payload})

                structures = payload.get("structures") or []
                if len(structures) != args.draws:
                    raise ValueError(f"{key} {arm}: expected {args.draws} structures")
                hashes = [hashlib.sha1(s["structure"].encode()).hexdigest()
                          for s in structures]
                if len(set(hashes)) != args.draws:
                    raise ValueError(f"{key} {arm}: duplicate diffusion structures")

                for draw, structure in enumerate(structures):
                    if draw in existing:
                        continue
                    text = structure["structure"]
                    stem = f"{key}_{arm}_matched_r{draw}"
                    relative = Path("structures") / f"{stem}.cif"
                    (OUT / relative).write_text(text)
                    shim = matched.sample_payload(payload, draw)
                    row: Dict[str, Any] = {
                        "iptm": run_pipeline._round(
                            run_pipeline.binder_interface_iptm(shim, len(polymers) - 1), 4),
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
                        "evidence_rank": subjects[key]["rank"],
                        "family": subjects[key]["family"],
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
                        "structure_file": str(relative),
                        "scored_at": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
                    })
                records = checkpoint_records()
                done = completed_call_keys(records, args.draws)
                print(f"[boltz2 {calls}] {key} {arm}: {args.draws} distinct draws", flush=True)
            if status == "call_cap_reached":
                break
        if status == "running":
            status = "completed"
    except BaseException:
        status = "failed"
        raise
    finally:
        records = checkpoint_records()
        atomic_json(OUT / "counter_screen_run.json", {
            "status": status,
            "phase": args.phase,
            "boltz2_calls_this_run": calls,
            "draws_per_call": args.draws,
            "subjects": list(subjects),
            "selected_designs_this_phase": wanted,
            "arms_this_phase": arms,
            "matched_constraint": matched.matched_constraint(),
            "settings": {k: SETTINGS[k] for k in subjects},
            "screen_advance_rule": (
                "G12D direct Asp6 >=2/3 AND apo/G12D median total buried area <1.0"
            ),
            "screen_summary": screen_summary(records, subjects),
            "historical_gates_changed": False,
            "track1_slots_spent": 0,
            "elapsed_seconds": round(time.time() - started, 1),
        })
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
