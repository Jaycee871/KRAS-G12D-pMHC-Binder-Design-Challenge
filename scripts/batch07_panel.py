#!/usr/bin/env python3
"""Batch-07 stage 3: the WT / apo / decoy panel, inside the gate.

Only designs that cleared the Asp6 direct read-out gate in **3 of 3**
independent draws (one Boltz-2 call, `diffusion_samples: 3`) get here:

  * `bb54_rigid-gp-p2` - rigidified redesign on bb54's backbone, 3/3 under
    `r3` (recycling_steps 3, sampling_steps 50, potentials on);
  * `B6-c6-46_w2_bb54_donor-anchored-g2` (bb54-g2, the batch-06 parent) - 3/3
    under `r3np` (recycling_steps 3, without_potentials), reproduced over two
    independent calls in step 1.

The WT-loss criterion is `scripts/batch06_panel.WT_CRITERION`, re-used
**verbatim** and not retuned: median `buried_area_peptide` >= 30 % lower on
VVVGAGGVGK and median position-6 heavy-atom contacts (all atoms, backbone
included) >= 30 % lower, with strict separation recorded alongside as the
stronger read. The Asp6 side-chain contact necessarily vanishing on WT is
excluded - it is true by construction. No ipTM gap, no MM-GBSA ddG.

The apo arm has no chain C to constrain against, so it is run **unconstrained**
exactly as `CONTROL_7OW6.md` did, and is therefore NOT like-for-like with the
three constrained arms. Every call is strictly serial through
`src/nvidia_client.py`; nothing in `src/design_gates.py` is touched.
"""
from __future__ import annotations

import argparse, json, sys, time
from pathlib import Path
from typing import Any, Dict, List

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(REPO / "scripts"))

import batch06_panel as P  # noqa: E402
import batch07_settings as S  # noqa: E402  (registers the r* sampling cells)
import run_batch04 as batch  # noqa: E402
import run_pipeline  # noqa: E402
from nvidia_client import NvidiaClient  # noqa: E402

OUT = REPO / "results" / "batch-07"
P.OUT = OUT
batch.OUT = OUT


def screened_designs() -> List[Dict[str, Any]]:
    rows = [json.loads(l) for l in
            (OUT / "screen_summary.jsonl").read_text().splitlines() if l.strip()]
    out = []
    for r in rows:
        if r.get("gate_3of3"):
            out.append({"name": r["design"], "sequence": r["sequence"],
                        "backbone": r.get("design", "").split("_")[0],
                        "variant": r["design"],
                        "binder_length": r.get("binder_length"),
                        "condition": r["stabilising_condition"],
                        "track": r.get("track"),
                        "rigidifying_edits": r.get("rigidifying_edits", []),
                        "rigidity": r.get("rigidity", {})})
    return out


def parent_design() -> Dict[str, Any]:
    d = S.subject_design()
    d["condition"] = "r3np"       # the reproducible 3-of-3 cell from step 1
    d["track"] = "parent"
    return d


def main(argv=None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--draws", type=int, default=3)
    parser.add_argument("--timeout", type=float, default=900.0)
    parser.add_argument("--max-boltz-calls", type=int, default=12)
    args = parser.parse_args(argv)

    (OUT / "structures").mkdir(parents=True, exist_ok=True)
    states = run_pipeline.load_target_states(REPO / "data" / "targets")
    client = NvidiaClient(timeout=args.timeout)

    designs = screened_designs() + [parent_design()]
    print(f"[panel] {len(designs)} designs at 3-of-3 on the Asp6 gate: "
          + ", ".join(f"{d['name']} ({d['condition']})" for d in designs))

    calls, started = 0, time.time()
    already = P.done_keys()
    summary: List[Dict[str, Any]] = []

    for design in designs:
        name, cond = design["name"], design["condition"]
        print(f"\n[panel] === {name} @ {cond} ===", flush=True)
        record: Dict[str, Any] = {k: design.get(k) for k in
                                  ("name", "sequence", "binder_length", "track",
                                   "rigidifying_edits", "rigidity")}
        record["condition"] = cond
        rows = P.load_rows()
        g12d = P.arm_rows(rows, name, "g12d", cond)
        record["g12d_gate_passes"] = sum(
            1 for r in g12d if str(r.get("asp6_gate_pass")).lower() == "true")
        record["g12d_draws"] = len(g12d)
        record["g12d_rows"] = g12d

        if calls < args.max_boltz_calls:
            calls += P.run_draws(client, design, "wt", cond, states, args.draws, already)
        rows = P.load_rows()
        verdict = P.wt_verdict(P.arm_rows(rows, name, "g12d", cond),
                               P.arm_rows(rows, name, "wt", cond))
        record["wt"] = verdict
        record["wt_rows"] = P.arm_rows(rows, name, "wt", cond)
        print(f"  -> WT arm: peptide-area loss {verdict['peptide_area_loss']}, "
              f"pos6-contact loss {verdict['pos6_contact_loss']}, "
              f"PASS={verdict['wt_gate_pass']}")
        if not verdict["wt_gate_pass"]:
            summary.append(record)
            continue
        for arm in ("apo", "decoy"):
            if calls >= args.max_boltz_calls:
                break
            calls += P.run_draws(client, design, arm, cond, states, args.draws, already)
        rows = P.load_rows()
        record["apo_rows"] = P.arm_rows(rows, name, "apo", cond)
        record["decoy_rows"] = P.arm_rows(rows, name, "decoy", cond)
        summary.append(record)

    out = {"wt_criterion": P.WT_CRITERION,
           "draws_per_arm": args.draws,
           "boltz2_calls_this_segment": calls,
           "elapsed_seconds": round(time.time() - started, 1),
           "designs": summary}
    (OUT / "panel_results.json").write_text(json.dumps(out, indent=2, default=str))
    print(f"\n[batch-07 panel] {calls} Boltz-2 calls in {time.time()-started:.0f}s")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
