#!/usr/bin/env python3
"""Batch-06 stage 4-6: replication AND the WT arm, both *inside* the gate.

Batch-05 folded replication into the gate. The off-target panel (PR #11) showed
that is not enough: all three designs that cleared it still bound the wild-type
peptide. So batch-06 folds the WT arm in too. A design passes only if it clears
**all** of:

  (a) the Asp6 direct read-out gate - a binder Arg/Lys/His/Ser/Thr/Tyr donor
      heavy atom within 3.5 A of Asp6 OD1/OD2 **and** `asp6_buried_frac` > 0
      (`scripts/run_batch04.gate_donor_distance`, thresholds untouched),
  (b) in **>= 3 of 3** independent constrained draws on the G12D state,
  (c) and **demonstrably loses engagement on the WT arm** - see WT_CRITERION
      below, judged on measured geometry only.

The Asp6 *side-chain* contact necessarily vanishing on WT (Gly has no side
chain) is true by construction and carries **no** information: it is excluded
from the WT criterion. No ipTM gap is computed anywhere; `specificity_gap` and
the paired MM-GBSA `ddG_sel` are retired.

Draws and the byte-identical-response artefact
----------------------------------------------
Repeat Boltz-2 calls with byte-identical request bodies sometimes return
byte-identical structures, so "n = 10" in PR #11 was really 6-8 distinct
samples. Here each set of 3 draws is taken as `diffusion_samples: 3` in a
*single* call: three independent diffusion trajectories from one request, which
cannot be duplicates of one another by construction. Every structure is hashed
and `distinct_structures` is reported per design/arm/condition.
Everything else is held identical to batch-04/05: same polymers, same
`step_scale` 1.638, same validated flat pocket constraint
`{"binder":"D","contacts":[{"id":"C","residue_index":6}],"max_distance":6.0}`,
same strictly-serial >= 8 s pacing, same analysis code.

The apo arm has no chain C to point the pocket constraint at, so it is run
**unconstrained**, exactly as `scripts/control_7ow6_constrained.py` did. It is
therefore **not** a like-for-like comparison with the three constrained arms.
"""
from __future__ import annotations

import argparse, hashlib, json, os, statistics, sys, time
from pathlib import Path
from typing import Any, Dict, List, Optional

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(REPO / "scripts"))

import run_batch04 as batch  # noqa: E402
import run_pipeline  # noqa: E402
import sampling_study as study  # noqa: E402
import selectivity_panel as panel  # noqa: E402
from nvidia_client import BOLTZ2_PATH, NvidiaClient  # noqa: E402

OUT = REPO / "results" / "batch-06"
batch.OUT = OUT

# Per-design stabilising settings stay in play: `recyc8` first (it collapsed
# two of three designs in the sampling study), `combo` as the fallback, because
# the best setting is design-specific and the two are not additive.
CONDITION_ORDER = ["recyc8", "combo"]

# --- the WT-loss criterion, stated before it is applied ---------------------
WT_AREA_LOSS_MIN = 0.30     # >= 30 % less peptide buried area on WT than G12D
WT_POS6_LOSS_MIN = 0.30     # >= 30 % fewer position-6 heavy-atom contacts
WT_CRITERION = (
    "A design clears the WT arm when, over 3 independent draws per state under "
    "the same sampling setting and the same pocket constraint, BOTH of: "
    f"(1) median buried_area_peptide on VVVGAGGVGK is at least {WT_AREA_LOSS_MIN:.0%} "
    "lower than on VVVGADGVGK, and (2) median position-6 heavy-atom contacts "
    f"(all atoms, backbone included) is at least {WT_POS6_LOSS_MIN:.0%} lower. "
    "Strict separation (max over WT draws < min over G12D draws, for both "
    "quantities) is recorded alongside as the stronger read. The disappearance "
    "of Asp6 side-chain contacts on WT is excluded: glycine has no side chain, "
    "so that loss is true by construction and carries no information. No ipTM "
    "quantity enters this criterion.")


def boltz_kwargs(condition: str) -> Dict[str, Any]:
    recyc, samp, nopot = study.CONDITIONS[condition]
    return {"recycling_steps": recyc, "sampling_steps": samp,
            "without_potentials": nopot}


def call_boltz(client: NvidiaClient, polymers, condition: str, constrained: bool,
               samples: int) -> Dict[str, Any]:
    body = {"polymers": [dict(p) for p in polymers], "step_scale": 1.638,
            "diffusion_samples": int(samples), **boltz_kwargs(condition)}
    if constrained:
        body["constraints"] = batch.POCKET
    return client.post(BOLTZ2_PATH, body)


def sample_payload(payload: Dict[str, Any], index: int) -> Dict[str, Any]:
    """One diffusion sample of a multi-sample payload, shaped like a 1-sample one.

    `pair_chains_iptm_scores` is only returned for sample 0 by the hosted NIM,
    so for later samples `binder_interface_iptm` falls back to the global ipTM.
    That is acceptable here and only here because ipTM is a pose-quality
    read-out and is never ranked on.
    """
    def pick(key):
        values = payload.get(key) or []
        return [values[index]] if index < len(values) else []
    return {"structures": [payload["structures"][index]],
            "pair_chains_iptm_scores": pick("pair_chains_iptm_scores"),
            "iptm_scores": pick("iptm_scores"),
            "ptm_scores": pick("ptm_scores"),
            "complex_ipde_scores": pick("complex_ipde_scores"),
            "complex_plddt_scores": pick("complex_plddt_scores")}


# --- checkpointing ----------------------------------------------------------

def checkpoint_path() -> Path:
    return OUT / "panel_checkpoint.jsonl"


def append_checkpoint(record: Dict[str, Any]) -> None:
    with checkpoint_path().open("a") as handle:
        handle.write(json.dumps(record, default=str) + "\n")
        handle.flush()
        os.fsync(handle.fileno())


def done_keys() -> set:
    path = checkpoint_path()
    if not path.exists():
        return set()
    out = set()
    for line in path.read_text().splitlines():
        if line.strip():
            r = json.loads(line)
            out.add((r["design"], r["arm"], r["condition"], r["draw"]))
    return out


def load_rows() -> List[Dict[str, Any]]:
    path = checkpoint_path()
    if not path.exists():
        return []
    return [json.loads(l) for l in path.read_text().splitlines() if l.strip()]


def arm_rows(rows, design, arm, condition) -> List[Dict[str, Any]]:
    return [r["row"] for r in rows if r["design"] == design and r["arm"] == arm
            and r["condition"] == condition]


# --- one set of draws -------------------------------------------------------

def run_draws(client, design: Dict[str, Any], arm: str, condition: str,
              states, draws: int, already: set) -> int:
    """One Boltz-2 call -> `draws` independent samples, checkpointed per sample."""
    key0 = (design["name"], arm, condition, 0)
    if key0 in already:
        print(f"[skip] {design['name']} {arm} {condition} already done")
        return 0
    sequence = design["sequence"]
    constrained = arm != "apo"
    polymers = list(states[arm]) + [
        {"id": run_pipeline.BINDER_CHAIN, "molecule_type": "protein",
         "sequence": sequence}]
    payload = call_boltz(client, polymers, condition, constrained, draws)
    structures = payload.get("structures") or []
    hashes = [hashlib.sha1(s["structure"].encode()).hexdigest() for s in structures]
    for index, structure in enumerate(structures):
        text = structure["structure"]
        stem = f"{design['name']}_{arm}_{condition}_r{index}"
        (OUT / "structures" / f"{stem}.cif").write_text(text)
        shim = sample_payload(payload, index)
        if arm == "g12d":
            # batch-04's scoring path verbatim, on this sample's structure
            saved = batch.boltz2_constrained
            batch.boltz2_constrained = lambda c, p, a, _s=shim: _s
            try:
                call_args = argparse.Namespace(**boltz_kwargs(condition))
                candidate = {k: design.get(k, "") for k in
                             ("sequence", "strategy", "backbone", "variant",
                              "binder_length", "mpnn_score", "donor_bias_score",
                              "design_pos6_contacts", "design_pos6_distance")}
                result = batch.score_candidate(client, candidate, states, call_args)
            finally:
                batch.boltz2_constrained = saved
            row = dict(result["row"])
            row["iptm"] = row.get("iptm_g12d")
            row["pos6_min_donor_distance"] = row.get("asp6_min_donor_distance")
            row["asp6_gate_pass"] = row.get("gate_pass")
            row["buried_area_hla"] = row.get("buried_area_hla")
        else:
            row = {"iptm": run_pipeline._round(
                       run_pipeline.binder_interface_iptm(shim, len(polymers) - 1), 4),
                   "pde_interface": run_pipeline._round(
                       run_pipeline._first(shim, "complex_ipde_scores"), 2)}
            try:
                metrics = panel.panel_metrics(text, sequence, arm)
                by_chain = metrics.get("buried_area_by_chain") or {}
                hla = max((v for k, v in by_chain.items()
                           if k not in (metrics.get("binder_chain"),
                                        metrics.get("peptide_chain"))), default=None)
                metrics["buried_area_hla"] = hla
                row.update(metrics)
                row["analysis_error"] = ""
            except Exception as exc:  # noqa: BLE001
                row["analysis_error"] = f"{type(exc).__name__}: {exc}"
        total = row.get("buried_area_total")
        pep = row.get("buried_area_peptide")
        try:
            row["peptide_area_share"] = round(float(pep) / float(total), 3) if total else None
        except (TypeError, ValueError):
            row["peptide_area_share"] = None
        append_checkpoint({
            "design": design["name"], "backbone": design.get("backbone", ""),
            "variant": design.get("variant", ""), "arm": arm,
            "condition": condition, "draw": index,
            "boltz_kwargs": {**boltz_kwargs(condition), "diffusion_samples": draws},
            "constraint": batch.POCKET if constrained else None,
            "constrained": constrained, "sequence": sequence,
            "structure_sha1": hashes[index],
            "distinct_structures_in_call": len(set(hashes)),
            "row": row, "structure_file": f"structures/{stem}.cif",
            "scored_at": time.strftime("%Y-%m-%dT%H:%M:%S%z")})
        print(f"  [{stem}] gate={row.get('asp6_gate_pass')} "
              f"donor={row.get('pos6_min_donor_distance')} "
              f"pos6={row.get('pos6_contacts')} pos6_sc={row.get('pos6_sidechain_contacts')} "
              f"pep_area={row.get('buried_area_peptide')} hla={row.get('buried_area_hla')} "
              f"share={row.get('peptide_area_share')} "
              f"pep_bf={row.get('peptide_buried_frac')}", flush=True)
    return 1


def median(values) -> Optional[float]:
    values = [float(v) for v in values if v is not None and v != ""]
    return statistics.median(values) if values else None


def wt_verdict(g12d_rows, wt_rows) -> Dict[str, Any]:
    g_area = [r.get("buried_area_peptide") for r in g12d_rows]
    w_area = [r.get("buried_area_peptide") for r in wt_rows]
    g_pos6 = [r.get("pos6_contacts") for r in g12d_rows]
    w_pos6 = [r.get("pos6_contacts") for r in wt_rows]
    ga, wa, gp, wp = median(g_area), median(w_area), median(g_pos6), median(w_pos6)
    area_loss = (1 - wa / ga) if (ga and wa is not None and ga > 0) else None
    pos6_loss = (1 - wp / gp) if (gp and wp is not None and gp > 0) else None
    def clean(xs):
        return [float(x) for x in xs if x is not None and x != ""]
    sep_area = (clean(w_area) and clean(g_area)
                and max(clean(w_area)) < min(clean(g_area)))
    sep_pos6 = (clean(w_pos6) and clean(g_pos6)
                and max(clean(w_pos6)) < min(clean(g_pos6)))
    passes = bool(area_loss is not None and pos6_loss is not None
                  and area_loss >= WT_AREA_LOSS_MIN
                  and pos6_loss >= WT_POS6_LOSS_MIN)
    return {"g12d_median_peptide_area": ga, "wt_median_peptide_area": wa,
            "peptide_area_loss": None if area_loss is None else round(area_loss, 3),
            "g12d_median_pos6_contacts": gp, "wt_median_pos6_contacts": wp,
            "pos6_contact_loss": None if pos6_loss is None else round(pos6_loss, 3),
            "strict_separation_area": bool(sep_area),
            "strict_separation_pos6": bool(sep_pos6),
            "wt_gate_pass": passes,
            "criterion": WT_CRITERION}


def main(argv=None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-designs", type=int, default=8)
    parser.add_argument("--draws", type=int, default=3)
    parser.add_argument("--timeout", type=float, default=900.0)
    parser.add_argument("--max-boltz-calls", type=int, default=60)
    args = parser.parse_args(argv)

    (OUT / "structures").mkdir(parents=True, exist_ok=True)
    states = run_pipeline.load_target_states(REPO / "data" / "targets")
    client = NvidiaClient(timeout=args.timeout)

    screened = [json.loads(l)["row"] for l in
                (OUT / "checkpoint.jsonl").read_text().splitlines() if l.strip()]
    passed = [r for r in screened if str(r.get("gate_pass")).lower() == "true"]
    # rank the screening survivors on the batch-06 objective: peptide share of
    # buried area, i.e. how peptide-centred the paratope actually came out.
    passed.sort(key=lambda r: -float(r.get("peptide_area_share") or 0))
    advanced = passed[: args.max_designs]
    designs = [{"name": f"{r['backbone']}_{r['variant']}", **r} for r in advanced]
    print(f"[panel] {len(screened)} screened, {len(passed)} cleared the Asp6 "
          f"read-out gate on the screening draw, advancing {len(designs)}")

    calls, started = 0, time.time()
    already = done_keys()
    summary: List[Dict[str, Any]] = []

    for design in designs:
        name = design["name"]
        print(f"\n[panel] === {name} (screening share "
              f"{design.get('peptide_area_share')}) ===", flush=True)
        record: Dict[str, Any] = {"design": name,
                                  "backbone": design.get("backbone"),
                                  "variant": design.get("variant"),
                                  "sequence": design["sequence"],
                                  "binder_length": design.get("binder_length"),
                                  "screening": {k: design.get(k) for k in
                                                ("buried_area_total", "buried_area_peptide",
                                                 "buried_area_hla", "peptide_area_share",
                                                 "peptide_buried_frac", "asp6_buried_frac",
                                                 "asp6_min_donor_distance",
                                                 "net_interface_charge", "iptm_g12d")},
                                  "conditions": {}}
        chosen = None
        for condition in CONDITION_ORDER:
            if calls >= args.max_boltz_calls:
                print("[stop] Boltz-2 call cap reached"); break
            calls += run_draws(client, design, "g12d", condition, states,
                               args.draws, already)
            rows = arm_rows(load_rows(), name, "g12d", condition)
            passes = sum(1 for r in rows
                         if str(r.get("asp6_gate_pass")).lower() == "true")
            record["conditions"][condition] = {
                "gate_passes": passes, "draws": len(rows),
                "rows": rows}
            print(f"  -> {condition}: Asp6 gate {passes}/{len(rows)}")
            if passes >= args.draws and len(rows) >= args.draws:
                chosen = condition
                break
        record["stabilising_condition"] = chosen
        record["gate_3of3"] = bool(chosen)
        if not chosen:
            summary.append(record)
            continue

        # --- the WT arm, inside the gate ---
        if calls < args.max_boltz_calls:
            calls += run_draws(client, design, "wt", chosen, states, args.draws, already)
        rows = load_rows()
        verdict = wt_verdict(arm_rows(rows, name, "g12d", chosen),
                             arm_rows(rows, name, "wt", chosen))
        record["wt"] = verdict
        print(f"  -> WT arm: peptide-area loss {verdict['peptide_area_loss']}, "
              f"pos6-contact loss {verdict['pos6_contact_loss']}, "
              f"PASS={verdict['wt_gate_pass']}")
        if not verdict["wt_gate_pass"]:
            summary.append(record)
            continue

        # --- apo and decoy, only for designs that clear the WT arm ---
        for arm in ("apo", "decoy"):
            if calls >= args.max_boltz_calls:
                break
            calls += run_draws(client, design, arm, chosen, states, args.draws, already)
        rows = load_rows()
        record["apo"] = arm_rows(rows, name, "apo", chosen)
        record["decoy"] = arm_rows(rows, name, "decoy", chosen)
        summary.append(record)

    out = {"wt_criterion": WT_CRITERION,
           "condition_order": CONDITION_ORDER,
           "draws_per_arm": args.draws,
           "boltz2_calls_this_segment": calls,
           "elapsed_seconds": round(time.time() - started, 1),
           "screened": len(screened), "screening_gate_passes": len(passed),
           "advanced": [d["name"] for d in designs],
           "designs": summary}
    (OUT / "panel_results.json").write_text(json.dumps(out, indent=2, default=str))
    print(f"\n[batch-06 panel] {calls} Boltz-2 calls in "
          f"{time.time() - started:.0f}s -> {OUT / 'panel_results.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
