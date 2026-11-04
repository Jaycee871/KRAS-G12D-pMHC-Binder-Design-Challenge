#!/usr/bin/env python3
"""Batch-08 stage 3-6: the gate, unchanged, on shape-steered designs.

`scripts/batch06_panel.py` is imported and re-used **verbatim** - the draws,
the checkpointing, the Asp6 read-out gate, the WT-loss criterion
(`batch06_panel.WT_CRITERION`, `batch06_panel.wt_verdict`), the apo handling.
Nothing is retuned here and nothing may be: the criterion was written before
batch-06's results and is applied to batch-08's results unchanged.

Two deliberate differences from batch-06, both about *sampling*, not gates:

* The sampling cell is `nopot` = recycling_steps 3, sampling_steps 50,
  without_potentials True - the cell batch-07 named **r3np**, the only cell in
  which bb54-g2 reached the Asp6 gate 3-of-3 twice (8/9 pooled).  Every result
  below is reported with its cell, because gate reproducibility and pose
  tightness are anti-correlated in this system.
* There is no separate 1-draw screening stage.  Screening *is* the 3-draw set:
  one Boltz-2 call with `diffusion_samples: 3`, which is both the screen and
  the 3-of-3 gate count, so no design is advanced on a single draw.

The WT arm is run only for designs at 3-of-3; apo (unconstrained, **not**
like-for-like) and decoy only for designs that clear the WT criterion.
"""
from __future__ import annotations

import argparse, json, sys, time
from pathlib import Path
from typing import Any, Dict, List

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(REPO / "scripts"))

import batch06_panel as p6  # noqa: E402
import run_batch04 as batch  # noqa: E402
import run_pipeline  # noqa: E402

OUT = REPO / "results" / "batch-08"
p6.OUT = OUT
batch.OUT = OUT

CONDITION = "nopot"   # = r3np: recycling 3, sampling 50, without_potentials


def rank_candidates(candidates: List[Dict[str, Any]], limit: int) -> List[Dict[str, Any]]:
    """Shape list first, then least helix contact, keeping -dw twins with parents."""
    def key(c):
        shape = c.get("shape") or {}
        return (0 if c.get("shape_list") == "shape" else 1,
                shape.get("helix_contact_frac", 1.0),
                -shape.get("peptide_contact_frac", 0.0),
                -float(c.get("donor_bias_score") or 0))
    by_backbone: Dict[str, List[Dict[str, Any]]] = {}
    for c in candidates:
        by_backbone.setdefault(c["backbone"], []).append(c)
    ordered: List[Dict[str, Any]] = []
    for backbone in sorted(by_backbone, key=lambda b: key(by_backbone[b][0])):
        group = sorted(by_backbone[backbone], key=lambda c: c["variant"])
        ordered.extend(group)
    return ordered[:limit]


def main(argv=None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-designs", type=int, default=36)
    parser.add_argument("--draws", type=int, default=3)
    parser.add_argument("--timeout", type=float, default=900.0)
    parser.add_argument("--max-boltz-calls", type=int, default=64)
    parser.add_argument("--candidates", default="candidates.json")
    parser.add_argument("--max-retries", type=int, default=8,
                        help="the hosted Boltz-2 NIM returns bursts of HTTP 429 "
                             "under load; the throttle is the primary defence "
                             "and this is the backstop")
    args = parser.parse_args(argv)

    (OUT / "structures").mkdir(parents=True, exist_ok=True)
    states = run_pipeline.load_target_states(REPO / "data" / "targets")
    client = p6.NvidiaClient(timeout=args.timeout, max_retries=args.max_retries)

    candidates = json.loads((OUT / args.candidates).read_text())
    advanced = rank_candidates(candidates, args.max_designs)
    designs = [{"name": f"{c['backbone']}_{c['variant']}", **c} for c in advanced]
    print(f"[panel] {len(candidates)} candidates, screening {len(designs)} at "
          f"{CONDITION} (r3np) with diffusion_samples={args.draws}", flush=True)

    calls, started = 0, time.time()
    already = p6.done_keys()
    summary: List[Dict[str, Any]] = []

    for design in designs:
        name = design["name"]
        shape = design.get("shape") or {}
        print(f"\n[panel] === {name} ({design.get('shape_list')} list, "
              f"helix_frac {shape.get('helix_contact_frac')}, "
              f"pep_frac {shape.get('peptide_contact_frac')}) ===", flush=True)
        record: Dict[str, Any] = {
            "design": name, "backbone": design.get("backbone"),
            "variant": design.get("variant"), "strategy": design.get("strategy"),
            "sequence": design["sequence"], "binder_length": design.get("binder_length"),
            "shape_list": design.get("shape_list"), "design_shape": shape,
            "donor_bias_score": design.get("donor_bias_score"),
            "helix_apolar_penalty": design.get("helix_apolar_penalty"),
            "dewet_positions": design.get("dewet_positions"),
            "conditions": {}}
        if calls >= args.max_boltz_calls:
            print("[stop] Boltz-2 call cap reached"); break
        try:
            calls += p6.run_draws(client, design, "g12d", CONDITION, states,
                                  args.draws, already)
        except Exception as exc:  # noqa: BLE001  - a 429 burst must not kill the run
            print(f"  !! {name} g12d draws failed: {type(exc).__name__}: {exc}",
                  flush=True)
            record["error"] = f"{type(exc).__name__}: {exc}"
            summary.append(record)
            time.sleep(60)
            continue
        rows = p6.arm_rows(p6.load_rows(), name, "g12d", CONDITION)
        passes = sum(1 for r in rows if str(r.get("asp6_gate_pass")).lower() == "true")
        record["conditions"][CONDITION] = {"gate_passes": passes, "draws": len(rows),
                                           "rows": rows}
        print(f"  -> {CONDITION}: Asp6 gate {passes}/{len(rows)}", flush=True)
        record["stabilising_condition"] = CONDITION if passes >= args.draws else None
        record["gate_3of3"] = bool(passes >= args.draws and len(rows) >= args.draws)
        if not record["gate_3of3"]:
            summary.append(record); continue

        if calls < args.max_boltz_calls:
            try:
                calls += p6.run_draws(client, design, "wt", CONDITION, states,
                                      args.draws, already)
            except Exception as exc:  # noqa: BLE001
                print(f"  !! {name} wt draws failed: {exc}", flush=True)
                record["wt_error"] = str(exc); summary.append(record); continue
        rows = p6.load_rows()
        verdict = p6.wt_verdict(p6.arm_rows(rows, name, "g12d", CONDITION),
                                p6.arm_rows(rows, name, "wt", CONDITION))
        record["wt"] = verdict
        print(f"  -> WT arm: peptide-area loss {verdict['peptide_area_loss']}, "
              f"pos6-contact loss {verdict['pos6_contact_loss']}, "
              f"PASS={verdict['wt_gate_pass']}", flush=True)
        if not verdict["wt_gate_pass"]:
            summary.append(record); continue

        for arm in ("apo", "decoy"):
            if calls >= args.max_boltz_calls:
                break
            try:
                calls += p6.run_draws(client, design, arm, CONDITION, states,
                                      args.draws, already)
            except Exception as exc:  # noqa: BLE001
                print(f"  !! {name} {arm} draws failed: {exc}", flush=True)
        rows = p6.load_rows()
        record["apo"] = p6.arm_rows(rows, name, "apo", CONDITION)
        record["decoy"] = p6.arm_rows(rows, name, "decoy", CONDITION)
        summary.append(record)

    out = {"wt_criterion": p6.WT_CRITERION, "condition": CONDITION,
           "condition_alias": "r3np (recycling_steps 3, sampling_steps 50, "
                              "without_potentials True)",
           "draws_per_arm": args.draws,
           "boltz2_calls_this_segment": calls,
           "elapsed_seconds": round(time.time() - started, 1),
           "candidates_available": len(candidates),
           "screened": [d["name"] for d in designs],
           "designs": summary}
    (OUT / "panel_results.json").write_text(json.dumps(out, indent=2, default=str))
    print(f"\n[batch-08 panel] {calls} Boltz-2 calls in {time.time() - started:.0f}s "
          f"-> {OUT / 'panel_results.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
