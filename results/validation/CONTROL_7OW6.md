# Positive control: the real 7OW6 TCR under the pocket-constrained protocol

**Date:** 2026-11-03 (Asia/Taipei)
**Scripts:** `scripts/control_7ow6_constrained.py`, `scripts/control_7ow6_replicates.py`,
`scripts/control_7ow6_rmsd.py`, `scripts/replicate_top3_gaps.py`,
`scripts/g12d_repeat_spread.py`
**Raw data:** `results/validation/7ow6_constrained_control.json`,
`7ow6_constrained_checkpoint.jsonl`, `7ow6_replicates.jsonl`,
`7ow6_control_rmsd.json`, `results/constrained-rescore/replicate_top3.jsonl`,
`results/constrained-rescore/g12d_repeat_spread.json`

## The question

`results/constrained-rescore/FINDINGS.md` left one item open before any new
batch was allowed to run: the pocket-constrained Boltz-2 protocol had never
been run on a binder known to be right. 7OW6 is that binder - the selective
TCR of Poole et al. 2022, which discriminates KRAS G12D `VVVGADGVGK` on
HLA-A\*11:01 from the wild-type peptide. The control asks two things:

1. **Pose** - does the constrained protocol put a known binder where it really
   sits?
2. **Discrimination** - does the protocol's ranking quantity, the ipTM
   specificity gap, come out clearly positive for a known selective binder?

## Method

Exactly the protocol a design goes through, with one unavoidable difference:
this binder is two chains.

* Polymers: the four target states from `data/targets` (chain A HLA-A\*11:01
  heavy, B β2m, C peptide) plus the 7OW6 TCR α chain as **D** (181 aa) and the
  β chain as **E** (245 aa), sequences taken from the deposited structure.
* Constraint: the flat NIM pocket term, one entry per binder chain -
  `{"binder": "D"|"E", "contacts": [{"id": "C", "residue_index": 4..8}],
  "max_distance": 6.0}` - applied to **g12d, wt and decoy** so the gap is a
  like-for-like comparison. `apo` has no peptide and is run unconstrained.
  This is the same constraint batch-03 was rescored under.
* Score: the largest ipTM between a *binder* chain and a *target* chain, the
  binder-internal α/β pair excluded. For a one-chain design this is exactly
  `run_pipeline.binder_interface_iptm`.
* Interface analysis: `complex_analysis.interface_metrics` on the predicted
  G12D complex with **E merged into D**, which is the convention
  `results/validation/7ow6_buried_area.json` (`tcr_as_binder`) already uses, so
  every number below is comparable with the crystal row.

Five predictions: g12d/wt/decoy constrained, apo unconstrained, plus an
unconstrained g12d as the A/B comparison. Each took 10-27 s.

## Result 1 - the pose: PASS, and quantitatively

| quantity | crystal 7OW6 (TCR as binder) | predicted, **constrained** | predicted, unconstrained |
| --- | --- | --- | --- |
| peptide buried fraction | 0.207 | **0.224** | 0.220 |
| buried area, total | 813 Å² | **825 Å²** | 828 Å² |
| buried area on the peptide | 169 Å² | **184 Å²** | 182 Å² |
| buried area on the HLA | 645 Å² | **640 Å²** | 646 Å² |
| contacts to the peptide | 33 | **31** | 31 |
| position-6 contacts | 7 | **7** | 8 |
| Asp6 buried fraction | 0.00 | **0.00** | 0.00 |
| Asp6 H-bond / salt bridge | no | **no** | no |
| net interface charge | −1.9 | **−2.9** | −2.9 |

C-alpha RMSD to the crystal after superposing on the pMHC frame only
(HLA + peptide, 286 residues, frame RMSD 0.85 Å) - i.e. the binder is *not*
re-fitted, so this measures whether it docks where the real one does:

| | constrained | unconstrained |
| --- | --- | --- |
| TCR α variable domain (115 res) | **1.29 Å** | 1.81 Å |
| TCR β variable domain (115 res) | **1.51 Å** | 1.78 Å |
| TCR α, full chain | 2.74 Å | 3.45 Å |
| TCR β, full chain | 2.56 Å | 2.68 Å |

The constrained protocol reproduces the native TCR pose to ~1.3-1.5 Å over the
variable domains and reproduces every buried-area number within 10 % of the
crystal. The geometry half of the pipeline is validated.

**Caveat worth stating:** the constraint *barely changed anything* here
(1.3 Å vs 1.8 Å). Boltz-2 already docks a real TCR correctly. So this control
proves the constraint does not **distort** a correct pose; it does not by
itself prove the constraint **rescues** a wrong one. The batch-03 rescore
(26/28 on the peptide vs 8/28) is the evidence for the rescue.

## Result 2 - the discrimination: FAIL

| state | ipTM (binder ↔ target) | ipTM (binder ↔ peptide chain only) |
| --- | --- | --- |
| **G12D** `VVVGADGVGK` | **0.9268 ± 0.0010** | 0.8966 ± 0.0009 |
| WT `VVVGAGGVGK` | **0.9325 ± 0.0001** | 0.9001 ± 0.0012 |
| decoy `AVFDRKSDAK` | 0.9136 ± 0.0005 | 0.8847 ± 0.0010 |
| apo (no peptide) | 0.7844 | — |

**specificity_gap = score(G12D) − max(WT, apo, decoy) = 0.9268 − 0.9325 =
−0.006.**

The known selective TCR gets a *negative* gap. Boltz-2 scores it as binding the
wild-type peptide very slightly **better** than the G12D peptide, and the
unrelated decoy peptide only 0.013 worse than the on-target one. Restricting
the ipTM to the binder↔peptide chain pair does not help: 0.8966 vs 0.9001, same
sign. The spread over repeat calls is 0.001, so the −0.006 is reproducible and
far outside the noise, not a coin flip.

The one thing ipTM does separate is **peptide present vs absent** (0.93 vs
0.78). That is a real but trivial signal: it says the TCR needs a peptide in
the groove, not that it needs *this* peptide.

## Result 3 - the gap is not even reproducible in the design regime

Prompted by the above, the full four-state constrained panel was re-run for the
three best batch-03 designs, byte-identical requests:

| design | gap, first pass | gap, replicate | ipTM(G12D) first → replicate |
| --- | --- | --- | --- |
| `S1-pos6-68_bb18` | **+0.298** | **−0.300** | 0.708 → 0.440 |
| `S1-pos6-68_bb19` | **+0.221** | **−0.041** | 0.680 → 0.417 |
| `S1-pos6-68_bb13` | +0.183 | **+0.151** | 0.758 → 0.758 |

And three further repeats of the *same* constrained G12D request:

| design | ipTM(G12D) over repeats | SD |
| --- | --- | --- |
| `S1-pos6-68_bb18` | 0.486, 0.619, 0.524 | 0.056 |
| `S1-pos6-68_bb13` | 0.277, 0.758, 0.758 | 0.227 |

For a *confident* complex (the real TCR) the constrained prediction is
reproducible to 0.001. For a *design* it is reproducible to 0.06-0.23 - one to
two orders of magnitude larger than any specificity gap this project has ever
measured. Likely mechanism: with `constraints` and `without_potentials: false`
the steering potentials make the diffusion stochastic, and a weak binder's
energy landscape has several docked modes that the sampling flips between; a
real TCR's funnel is deep enough that it never flips. (Consistent with this,
the *unconstrained* apo state reproduced to four decimal places.)

## Verdict

**The control fails**, on the half that matters for ranking.

* Geometry: PASS - near-native pose, crystal-matching buried areas.
* Specificity gap: **FAIL** - the one binder we know to be selective scores
  −0.006, and design-regime gaps of ±0.3 are not reproducible across identical
  calls.

Per the plan agreed before this run, **no new design batch was launched.**
Batch-04 would have spent ~14 hours of Boltz-2 calls ranking designs by a
quantity that, measured on the only ground truth available, carries no
mutation-specificity signal at all and in the design regime is dominated by
sampling noise. `results/batch-04/SUMMARY.md` records that decision.

## What the real TCR scores on every one of our hard gates

Measured on the constrained G12D prediction (crystal value in brackets where it
differs):

| hard gate (PIPELINE.md §6) | threshold | real 7OW6 TCR | verdict |
| --- | --- | --- | --- |
| peptide share of buried area | ≥ 0.40 | **0.224** (crystal 0.207) | **FAIL** |
| Asp6 read-out (H-bond/salt bridge **or** ≥ 7 pos-6 contacts) | — | no bond, **7** contacts | pass, on the Poole escape clause exactly |
| position-6 heavy-atom contacts | ≥ 3 | 7 | pass |
| net interface charge | −2.0 … +4.0 | **−2.9** (crystal −1.9) | **FAIL** |
| PDE-interface(G12D) | < 10 Å | 0.40 | pass |
| predicted-pose filter (≥ 25 Å² on peptide, ≥ 1 pos-6 contact) | — | 184 Å², 7 | pass |
| *(ranking)* ipTM(G12D) | ≥ 0.75 | 0.927 | meets |
| *(ranking)* max(WT, apo, decoy) | ≤ 0.50 | 0.933 | misses |
| *(ranking)* specificity gap | ≥ 0.25 | −0.006 | misses |
| *(ranking)* asp6_buried_frac | ≥ 0.75 | 0.00 | misses |

**The real, experimentally validated, selective binder fails two of our six
hard gates and three of our four ranking targets.** Any gate set a known-good
binder cannot pass is a gate set that is measuring the wrong thing.

## Proposed gate recalibration (evidence, not applied)

Step 2 of the plan was conditional on the control passing, so
`src/design_gates.py` has **not** been changed. The evidence for the changes is
recorded here so the recalibration is a one-line review rather than a guess:

| threshold | now | proposed | evidence |
| --- | --- | --- | --- |
| `PEPTIDE_BURIED_FRACTION_MIN` | 0.40 | **0.18** | The real TCR measures 0.207 (crystal) / 0.224 (predicted). 0.40 is unreachable by a known selective binder because the groove, not the binder, owns most of the peptide surface. 0.18 ≈ 0.85 × the control, i.e. "at least 85 % of what a real selective TCR achieves". |
| `NET_INTERFACE_CHARGE_MIN` | −2.0 | **−3.5** | The real TCR measures −1.9 (crystal) / −2.9 (predicted). The Xie et al. caution against excess negative charge is sound, but a floor that a validated TCR trips at −2.9 is too tight; −3.5 keeps the caution while admitting the control. |
| `SPECIFICITY_GAP_MIN` = 0.25 as a *ranking* term | 0.25 | **withdraw** | The control's gap is −0.006 and design-regime gaps swing ±0.3 between identical calls. Ranking on it ranks noise. |
| `ON_TARGET_IPTM_MIN` 0.75 | 0.75 | keep, but as a *reproducibility-checked* term | The control sits at 0.927; designs that reach 0.75 **in every one of ≥ 3 repeat calls** are a meaningful shortlist. A single-call 0.75 is not. |
| Asp6 read-out gate | — | **keep unchanged** | The control passes it exactly on the Poole escape clause (7 contacts, no bond). This gate was already calibrated on 7OW6 and the control confirms it. |
| pose filter | — | **keep unchanged** | The control passes comfortably (184 Å², 7 contacts). |

Under the proposed thresholds, **one** batch-03 design clears every hard gate:
`S1-pos6-68_bb13` (peptide buried 0.223, 12 position-6 contacts, Asp6 hydrogen
bond *and* salt bridge, interface charge −1.0, PDE 0.4, pose OK) - and it is
also the only one of the top three whose score vector reproduces.

## What to do instead of another ipTM batch

1. **Stop ranking on co-folder ipTM.** It is validated here as a pose/plausibility
   score and invalidated as a specificity score. Keep it as a sanity filter
   (`is this thing even a complex`), not as the ranking column.
2. **Replicate before you believe a number.** Any score used for ranking needs
   ≥ 3 identical calls and a reported spread. Three repeats cost 24 s; a wrong
   leaderboard costs a submission slot.
3. **Score the mutation directly, not the complex globally.** The quantity the
   competition wants is ΔΔG(G12D → WT) for a *fixed* pose. That means a paired
   calculation on one structure - relax the predicted G12D complex, mutate
   position 6 to Gly in place, re-score - rather than two independent
   co-folding runs whose difference is swamped by docking noise. FoldX/Rosetta
   ddG or an MM-GBSA-style rescoring on the box would do this for free, with no
   API budget at all.
4. **Use the geometry that the control *does* validate.** Position-6 contact
   density, an Asp6 hydrogen bond or salt bridge, and the peptide share of
   buried area are reproducible, physically meaningful, and the real TCR gives
   us calibrated values for all three. A geometry-only ranking is defensible
   today; an ipTM-gap ranking is not.
