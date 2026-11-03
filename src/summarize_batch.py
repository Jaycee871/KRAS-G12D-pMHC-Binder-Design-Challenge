#!/usr/bin/env python3
"""Turn a `run_batch.py` results directory into the batch summary markdown.

Reads `designs.csv` + `run_stats.json` from the batch directory and writes
`SUMMARY.md`: per-gate survival, per-strategy comparison, the leaderboard by
specificity gap, **where on the target each binder actually landed**, and what
to change in the next batch - using the thresholds in `src/design_gates.py`
(and therefore `docs/PIPELINE.md` section 6).

The "where it landed" section exists because batch-01's all-zero buried
fractions turned out to be real rather than a bug: the binders were docking
onto the HLA / beta-2m face instead of onto the peptide, and the old summary
had no column that would have said so.
"""
from __future__ import annotations

import argparse
import csv
import json
import statistics
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List

sys.path.insert(0, str(Path(__file__).resolve().parent))
import design_gates  # noqa: E402

REPO_ROOT = Path(__file__).resolve().parents[1]
VALIDATION_JSON = REPO_ROOT / "results" / "validation" / "7ow6_buried_area.json"

def _b(row, key):
    raw = row.get(key, "")
    if raw is None or (isinstance(raw, str) and not raw.strip()):
        return None
    return str(raw).strip().lower() in ("true", "1", "yes")


# Hard gates (geometry and chemistry). The ipTM thresholds are ranking terms
# since batch-02 and are reported separately.
GATES = [
    ("peptide_buried_frac >= 0.40",
     lambda r: (_f(r, "peptide_buried_frac") or 0.0)
     >= design_gates.PEPTIDE_BURIED_FRACTION_MIN),
    ("Asp6 H-bond / salt bridge (or >= 7 pos6 contacts)",
     lambda r: not design_gates.asp6_readout_failures(r)),
    ("pos6 contacts >= 3",
     lambda r: (_f(r, "pos6_contacts") or 0.0)
     >= design_gates.MIN_POS6_HEAVY_ATOM_CONTACTS),
    ("net interface charge in -2..+4",
     lambda r: _f(r, "net_interface_charge") is not None
     and design_gates.NET_INTERFACE_CHARGE_MIN
     <= _f(r, "net_interface_charge")
     <= design_gates.NET_INTERFACE_CHARGE_MAX),
    ("PAE/PDE-interface(G12D) < 10 A",
     lambda r: (_f(r, "pae_interface_g12d") if _f(r, "pae_interface_g12d")
                is not None else 99) < design_gates.PAE_INTERFACE_MAX),
    ("pose filter (binder on the peptide)",
     lambda r: _b(r, "pose_ok") is True),
]

# Ranking targets: reported, never vetoing.
RANKING_TARGETS = [
    ("ipTM(G12D) >= 0.75",
     lambda r: _f(r, "iptm_g12d") is not None
     and _f(r, "iptm_g12d") >= design_gates.ON_TARGET_IPTM_MIN),
    ("max(WT, apo, decoy) <= 0.50",
     lambda r: max(_f(r, c) or 0.0 for c in
                   ("iptm_wt", "iptm_apo", "iptm_decoy"))
     <= design_gates.OFF_TARGET_IPTM_MAX),
    ("specificity_gap >= 0.25",
     lambda r: (_f(r, "specificity_gap") or -99)
     >= design_gates.SPECIFICITY_GAP_MIN),
    ("asp6_buried_frac >= 0.75 (preferred)",
     lambda r: (_f(r, "asp6_buried_frac") or 0.0)
     >= design_gates.ASP6_BURIED_FRACTION_PREFERRED),
]


def _f(row, key):
    raw = row.get(key, "")
    if raw is None or (isinstance(raw, str) and not raw.strip()):
        return None
    try:
        return float(raw)
    except ValueError:
        return None


def _mean(values):
    values = [v for v in values if v is not None]
    return statistics.fmean(values) if values else float("nan")


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dir", required=True)
    parser.add_argument("--title", default="Batch 01 - first real design batch")
    parser.add_argument("--validation", default=str(VALIDATION_JSON))
    args = parser.parse_args(argv)
    out_dir = Path(args.dir)

    rows: List[Dict[str, str]] = list(csv.DictReader(
        (out_dir / "designs.csv").open()))
    stats = json.loads((out_dir / "run_stats.json").read_text())
    n = max(len(rows), 1)

    lines: List[str] = []
    w = lines.append
    w(f"# {args.title}")
    w("")
    w(f"Run with `src/run_batch.py` against the hosted NVIDIA BioNeMo NIMs "
      f"(RFdiffusion 2.3, ProteinMPNN, Boltz-2), {stats['elapsed_seconds']:.0f} s "
      f"wall clock, {stats['total_calls']} API calls "
      f"({', '.join(f'{k} {v}' for k, v in sorted(stats['call_counts'].items()))}).")
    w("")
    w("## What was produced")
    w("")
    w("| quantity | count |")
    w("| --- | --- |")
    backbones = stats.get("backbones") or len(
        {r.get("backbone") for r in rows if r.get("backbone")})
    w(f"| RFdiffusion backbones | "
      f"{stats.get('backbones_generated', backbones)} |")
    w(f"| candidates offered for scoring | {stats['raw_sequences']} |")
    if stats.get("backbones_generated") is not None:
        w(f"| RFdiffusion backbones passing the design pose filter | "
          f"{stats.get('backbones_pose_ok')} |")
    w(f"| candidates scored (>= the G12D state) | {stats['scored']} |")
    if stats.get("pose_filter_survivors") is not None:
        w(f"| candidates surviving the Boltz-2 pose filter | "
          f"{stats['pose_filter_survivors']} |")
    w(f"| candidates with a complete four-state score | "
      f"{stats.get('scored_four_state', stats['scored'])} |")
    w(f"| Boltz-2 predictions | {stats['call_counts'].get('boltz2', 0)} |")
    w(f"| designs clearing **every** section 6 hard gate | {stats['passing']} |")
    if stats.get("analysis_errors") is not None:
        w(f"| designs whose interface analysis raised an error | "
          f"{stats['analysis_errors']} |")
    w("")
    retries = stats.get("retry_statuses") or {}
    w(f"Retryable HTTP statuses hit: "
      f"{', '.join(f'{k} x{v}' for k, v in retries.items()) if retries else 'none'}.")
    throttle = stats.get("throttle") or {}
    boltz = (throttle.get("endpoints") or {}).get(
        "/v1/biology/mit/boltz2/predict") or {}
    if boltz:
        w("")
        w(f"Boltz-2 throttle: strictly serial, {boltz.get('min_interval')} s "
          f"minimum interval between requests, {boltz.get('grants')} slots "
          f"granted, {boltz.get('waited_seconds')} s spent waiting on the "
          f"throttle (the alternative, measured in "
          f"`results/calibration/`, is 31 x HTTP 429 and 6 lost calls out of "
          f"12 at 6 workers).")
    w("")

    w("## Strategies")
    w("")
    w("| strategy | binder length | hotspots | MPNN T | rationale |")
    w("| --- | --- | --- | --- | --- |")
    for s in stats["strategies"]:
        w(f"| `{s['name']}` | {s['binder_length']} | "
          f"{' '.join(s['hotspots'])} | {s['sampling_temp']} | {s['rationale']} |")
    w("")

    w("## Gate survival (of the scored candidates)")
    w("")
    w("| gate | passing | % |")
    w("| --- | --- | --- |")
    for label, test in GATES:
        passing = sum(1 for r in rows if test(r))
        w(f"| {label} | {passing}/{len(rows)} | {100 * passing / n:.0f} % |")
    w(f"| **all hard gates** | {stats['passing']}/{len(rows)} | "
      f"{100 * stats['passing'] / n:.0f} % |")
    w("")
    w("Ranking targets (ipTM; demoted from hard gates in batch-02 because the "
      "reported ipTM is the strongest binder-to-*any*-chain interface and so "
      "says nothing about where the binder sits):")
    w("")
    w("| ranking target | meeting it | % |")
    w("| --- | --- | --- |")
    for label, test in RANKING_TARGETS:
        meeting = sum(1 for r in rows if test(r))
        w(f"| {label} | {meeting}/{len(rows)} | {100 * meeting / n:.0f} % |")
    w("")

    w("## Per-strategy results")
    w("")
    by_strategy = defaultdict(list)
    for r in rows:
        by_strategy[r["strategy"]].append(r)
    w("| strategy | scored | mean ipTM(G12D) | best ipTM(G12D) | mean gap | "
      "best gap | mean peptide buried frac | mean buried area on peptide | "
      "designs touching the peptide | passing all gates |")
    w("| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |")
    for name in sorted(by_strategy):
        group = by_strategy[name]
        gaps = [g for g in (_f(r, "specificity_gap") for r in group)
                if g is not None]
        touching = sum(1 for r in group
                       if (_f(r, "buried_area_peptide") or 0.0) > 0)
        w(f"| `{name}` | {len(group)} | "
          f"{_mean([_f(r, 'iptm_g12d') for r in group]):.3f} | "
          f"{max([_f(r, 'iptm_g12d') or 0 for r in group]):.3f} | "
          f"{_mean(gaps):+.3f} | {max(gaps or [float('nan')]):+.3f} | "
          f"{_mean([_f(r, 'peptide_buried_frac') for r in group]):.2f} | "
          f"{_mean([_f(r, 'buried_area_peptide') for r in group]):.0f} A^2 | "
          f"{touching}/{len(group)} | "
          f"{sum(1 for r in group if not r['gate_failures'])} |")
    w("")

    # --- where the binders landed ------------------------------------------
    touching = [r for r in rows if (_f(r, "buried_area_peptide") or 0.0) > 0]
    pos6 = [r for r in rows if (_f(r, "pos6_contacts") or 0.0) >= 3]
    no_contact = [r for r in rows if (_f(r, "buried_area_total") or 0.0) <= 0]

    w("## Where the binders actually landed")
    w("")
    w("Buried area is the change in solvent-accessible surface on the "
      "*target* side: the same predicted coordinates with the binder chain "
      "present and with it deleted. A design can score a respectable ipTM by "
      "gripping the HLA alpha-3 face or beta-2m and never touch the peptide, "
      "which is what most of this batch does, so the split matters far more "
      "than the total.")
    w("")
    w("| quantity | value |")
    w("| --- | --- |")
    w(f"| designs burying **any** peptide surface | {len(touching)}/{len(rows)} |")
    w(f"| designs with >= 3 heavy-atom contacts to peptide position 6 | "
      f"{len(pos6)}/{len(rows)} |")
    w(f"| designs that touch the target nowhere at all | "
      f"{len(no_contact)}/{len(rows)} |")
    w(f"| mean buried area, total | "
      f"{_mean([_f(r, 'buried_area_total') for r in rows]):.0f} A^2 |")
    w(f"| mean buried area on the peptide | "
      f"{_mean([_f(r, 'buried_area_peptide') for r in rows]):.0f} A^2 |")
    w(f"| mean buried area on the HLA heavy chain | "
      f"{_mean([_f(r, 'buried_area_hla') for r in rows]):.0f} A^2 |")
    w(f"| best peptide-buried fraction seen | "
      f"{max([_f(r, 'peptide_buried_frac') or 0.0 for r in rows] or [0]):.2f} |")
    w("")

    validation = Path(args.validation)
    if validation.exists():
        cases = json.loads(validation.read_text())["cases"]
        groove = cases["hla_groove_as_binder"]
        tcr = cases["tcr_as_binder"]
        w("### Positive control: PDB 7OW6")
        w("")
        w("The same analysis code run on the real HLA-A*11:01 / beta-2m / "
          "VVVGADGVGK complex (`scripts/validate_7ow6.py`), treating a known "
          "binding partner as the binder:")
        w("")
        w("| reference binder | peptide buried frac | Asp6 side-chain buried "
          "frac | pos6 contacts |")
        w("| --- | --- | --- | --- |")
        w(f"| HLA groove (chain A) | {groove['peptide_buried_frac']:.2f} | "
          f"**{groove['asp6_buried_frac']:.2f}** "
          f"({groove['asp6_sasa_free']:.0f} -> "
          f"{groove['asp6_sasa_bound']:.0f} A^2) "
          f"| {groove['pos6_contacts']:.0f} |")
        w(f"| TCR (chains D+E) | {tcr['peptide_buried_frac']:.2f} | "
          f"{tcr['asp6_buried_frac']:.2f} "
          f"({tcr['asp6_sasa_free']:.1f} -> {tcr['asp6_sasa_bound']:.1f} A^2) "
          f"| {tcr['pos6_contacts']:.0f} |")
        w("")
        w("The second row is the uncomfortable one: in the real complex the "
          "G12D aspartate is **already ~90 % buried by the groove itself** "
          f"(only {tcr['asp6_sasa_free']:.1f} A^2 of side chain left exposed), "
          "and the selective TCR of Poole et al. 2022 buries none of what "
          "remains while still making "
          f"{tcr['pos6_contacts']:.0f} heavy-atom contacts with position 6. "
          "Measured this way the `asp6_buried_frac >= 0.60` gate was "
          "unreachable by construction, and batch-02 replaced it with the "
          "Asp6 hydrogen-bond / salt-bridge read-out criterion.")
        w("")

    if stats.get("pose_filter"):
        pf = stats["pose_filter"]
        survivors = stats.get("pose_filter_survivors", 0)
        w("## Pose filter")
        w("")
        w("Two filters, both new in batch-02. The first reads the "
          "*RFdiffusion* complex and rejects a backbone whose binder is not "
          "on the peptide before any sequence is designed or any Boltz-2 call "
          "is spent. The second reads the Boltz-2 **G12D** prediction and "
          "stops there when the binder has been re-docked away from the "
          "peptide, saving the WT, apo and decoy calls.")
        w("")
        w("| quantity | value |")
        w("| --- | --- |")
        w(f"| enabled | {pf['enabled']} |")
        w(f"| design filter: peptide contacts / residues / pos6 contacts / "
          f"pos6 distance | >= {pf['design_min_peptide_contacts']} / >= "
          f"{pf['design_min_peptide_residues']} / >= "
          f"{pf['design_min_pos6_contacts']} / <= "
          f"{pf['design_max_pos6_distance']} A |")
        w(f"| backbones generated | {stats.get('backbones_generated')} |")
        w(f"| backbones passing the design filter | "
          f"{stats.get('backbones_pose_ok')} |")
        w(f"| predicted-pose filter: peptide area / pos6 contacts | >= "
          f"{pf['predicted_min_peptide_area']} A^2 / >= "
          f"{pf['predicted_min_pos6_contacts']} |")
        w(f"| candidates surviving the predicted-pose filter | "
          f"{survivors}/{len(rows)} |")
        w(f"| Boltz-2 calls saved by it | "
          f"{3 * (len(rows) - survivors)} |")
        w("")

    w("## Leaderboard - top 10 by specificity gap")
    w("")
    w("| # | name | strategy | ipTM G12D | WT | apo | decoy | gap | "
      "peptide buried | buried on peptide | pos6 contacts | Asp6 H-bond | "
      "salt bridge | asp6 buried | iface charge | gates failed |")
    w("| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | "
      "--- | --- | --- | --- | --- |")
    for i, r in enumerate(rows[:10], start=1):
        w(f"| {i} | `{r['name']}` | `{r['strategy']}` | {r['iptm_g12d']} | "
          f"{r['iptm_wt']} | {r['iptm_apo']} | {r['iptm_decoy']} | "
          f"{r['specificity_gap']} | {r['peptide_buried_frac']} | "
          f"{r.get('buried_area_peptide', '')} | {r['pos6_contacts']} | "
          f"{r.get('asp6_hbond', '')} | {r.get('asp6_salt_bridge', '')} | "
          f"{r['asp6_buried_frac']} | {r['net_interface_charge']} | "
          f"{r['n_gate_failures']} |")
    w("")

    w("### Most common gate failures")
    w("")
    counts = defaultdict(int)
    for r in rows:
        for failure in (r["gate_failures"] or "").split("; "):
            if failure:
                parts = failure.split(" ")
                counts[" ".join(parts[:2])] += 1
    for key, count in sorted(counts.items(), key=lambda kv: -kv[1]):
        w(f"* `{key}` - {count} designs")
    w("")

    # --- recommendations ----------------------------------------------------
    best_gap = max([_f(r, "specificity_gap") or -99 for r in rows] or [-99])
    best_iptm = max([_f(r, "iptm_g12d") or 0.0 for r in rows] or [0.0])
    best_strategy = None
    if by_strategy:
        best_strategy = max(
            by_strategy,
            key=lambda k: max(
                [_f(r, "specificity_gap") or -99 for r in by_strategy[k]]),
        )

    survivors = stats.get("pose_filter_survivors", 0)
    hbond = sum(1 for r in rows if _b(r, "asp6_hbond"))
    four_state = stats.get("scored_four_state", len(rows))

    w("## What to change next")
    w("")
    w(f"1. **Pose, again.** {len(touching)}/{len(rows)} scored designs bury "
      f"any peptide surface and {survivors}/{len(rows)} survived the "
      f"predicted-pose filter; {four_state} were worth a full four-state "
      f"panel. The filter is doing its job as a budget saver, but it cannot "
      f"create a pose Boltz-2 does not predict: the next step is to score the "
      f"*designed* complex directly (a predictor that accepts a template or a "
      f"docked pose), rather than asking a co-folder to rediscover it from "
      f"sequence.")
    w(f"2. **Asp6 read-out.** {hbond}/{len(rows)} designs hydrogen bond or "
      f"salt bridge to the Asp6 carboxylate. Note what the 7OW6 control says: "
      f"HLA Arg114 already caps that carboxylate, so a binder has to compete "
      f"with the groove for it - which is why the Poole TCR reads the "
      f"aspartate indirectly and why a dense position-6 contact footprint "
      f"substitutes for the bond in the gate.")
    w(f"3. **ipTM stays a ranking term.** Best ipTM(G12D) this batch was "
      f"{best_iptm:.3f} and the best specificity gap {best_gap:+.3f}. Until a "
      f"design clears the geometry gates there is nothing for an ipTM "
      f"threshold to filter.")
    if best_strategy:
        w(f"4. **Best strategy so far: `{best_strategy}`.** Concentrate "
          f"backbones there, and keep spending the budget on *more "
          f"backbones*, not more sequences per backbone - same-backbone "
          f"sequences give correlated poses.")
    w(f"5. **Keep the serial Boltz-2 throttle.** The calibrated sustainable "
      f"rate is ~8 requests/minute from one key. With the pose filter an "
      f"off-peptide candidate now costs one call instead of four, so batch "
      f"size should be planned as "
      f"`candidates + 3 x (candidates surviving the pose filter)`.")
    w("")

    (out_dir / "SUMMARY.md").write_text("\n".join(lines) + "\n")
    print(f"wrote {out_dir / 'SUMMARY.md'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
