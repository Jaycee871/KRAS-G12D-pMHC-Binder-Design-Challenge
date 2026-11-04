# Sampling-settings study — can a Boltz-2 setting collapse the pose spread?

**Date:** 2026-11-03 (Asia/Taipei).
**Scripts:** `scripts/sampling_study.py` (re-uses `scripts/run_batch04.py`'s
`score_candidate` verbatim — same constraint, same analysis, same gate) and
`scripts/sampling_study_analyze.py` (re-uses
`scripts/replicate_rmsd_batch04.py`'s pMHC-frame RMSD).
**Raw data:** `checkpoint.jsonl` (39 rows, authoritative), `structures/`
(39 mmCIF), `analysis.json`.
**Boltz-2 calls:** **39**. Strictly serial, ≥ 8 s interval, one process.
**Track-1 design slots consumed: 0 of 40.** No design was generated: every
sequence scored here already existed in batch-04/batch-05. Nothing submitted,
nothing pushed to `main`, no PR merged.

## The question

Batches 04 and 05 both ended with the same diagnosis: **the chemistry is
solved, the pose is not.** The designs make the shortest Asp6 donor contacts in
the project (2.03 Å) but only intermittently, and in every batch-05 near-miss
the gate failure and a flip between two docked modes are the *same event*. The
pMHC reference frame is stable (0.08–0.33 Å), so the 4–20 Å binder RMSD spreads
are real, not fitting artefacts.

Recommendation 1 of `results/batch-05/SUMMARY.md` was therefore to test whether
changing Boltz-2's **sampling settings** collapses the replicate spread. This
is that experiment.

## Design

Five existing designs — the four batch-05 near-misses plus batch-04's `bb25`,
the project's strongest candidate:

| key | design | source | baseline gate |
| --- | --- | --- | --- |
| `bb135-g1` | `B5-c6-62_bb135` donor-anchored-g1 | batch-05 | 2 / 3 |
| `bb173-g1` | `B5-c6-62_bb173` donor-anchored-g1 | batch-05 | 2 / 3 |
| `bb90-g2` | `B5-c6-62_bb90` donor-anchored-g2 | batch-05 | 1 / 3 |
| `bb50-g1` | `B5-c6-62_bb50` donor-anchored-g1 | batch-05 | 1 / 3 |
| `bb25` | `B4-c6-62_bb25` donor-anchored | batch-04 | 2 / 3 |

Four conditions, 3 independent draws each:

| condition | `recycling_steps` | `sampling_steps` | `without_potentials` | calls |
| --- | --- | --- | --- | --- |
| **baseline** (batch-05) | 3 | 50 | false | 0 — re-used from batch-04/05 replicates |
| **nopot** | 3 | 50 | **true** | 15 |
| **recyc8** | **8** | 50 | false | 15 |
| **combo** | **8** | 50 | **true** | 9 (`bb135-g1`, `bb50-g1`, `bb25` — the three that reached 3/3 under a single setting) |

Everything else is identical to the batch-05 run: same polymers, same
`diffusion_samples: 1`, same `step_scale: 1.638`, and the validated flat pocket
constraint `{"binder":"D","contacts":[{"id":"C","residue_index":6}],"max_distance":6.0}`
(nested or `[chain, index]` pair forms return HTTP 422). The baseline is **not**
re-run — the batch-04/batch-05 replicate draws *are* the baseline, so the
comparison costs nothing. The screening draw is excluded everywhere, so every
cell is a like-for-like 3-draw set.

Gate, unchanged: a binder Arg/Lys/His/Ser/Thr/Tyr heavy atom within 3.5 Å of
Asp6 OD1/OD2 **and** `asp6_buried_frac > 0`. **No gate or threshold was
loosened or modified**; `src/design_gates.py` is untouched.

Pose spread: binder (chain D) C-alpha RMSD after superposing on the **pMHC
frame only** (HLA chain A + peptide chain C C-alphas), binder **not** re-fitted;
all 3 pairwise comparisons per condition. `frame` below is the largest pMHC
frame RMSD in that set — the fitting-artefact floor.

## Results

| design | condition | gate / 3 | `asp6_min_donor_distance` per draw (Å) | Asp6 **side-chain** contacts per draw | `peptide_buried_frac` mean | net interface charge per draw | ipTM* mean | **binder RMSD mean / max (Å)** | frame (Å) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `bb135-g1` | baseline | 2 / 3 | 2.03, —, 2.03 | 5, 1, 5 | 0.333 | 0.0, +1.0, 0.0 | 0.451 | 4.54 / 5.81 | 0.33 |
| `bb135-g1` | **nopot** | **3 / 3** | 2.55, 2.69, 2.70 | 5, 4, 5 | 0.331 | 0.0, 0.0, 0.0 | 0.500 | **2.67 / 3.34** | 0.17 |
| `bb135-g1` | recyc8 | 2 / 3 | —, 2.66, 2.67 | 1, 5, 5 | 0.294 | 0.0, 0.0, 0.0 | 0.521 | 5.94 / 8.78 | 0.22 |
| `bb135-g1` | **combo** | **3 / 3** | 2.81, 2.42, 2.59 | 5, 5, 5 | 0.334 | 0.0, 0.0, 0.0 | 0.469 | **2.04 / 2.97** | 0.25 |
| `bb173-g1` | baseline | 2 / 3 | 2.48, 2.67, — | 6, 6, 6 | 0.297 | −1.0, +1.0, 0.0 | 0.531 | 8.78 / 11.79 | 0.25 |
| `bb173-g1` | nopot | 1 / 3 | —, —, 2.76 | 11, 1, 6 | 0.282 | +1.0, 0.0, 0.0 | 0.402 | 9.06 / 12.72 | 0.26 |
| `bb173-g1` | recyc8 | 2 / 3 | 2.09, 3.20, — | 8, 8, 1 | 0.272 | +1.0, +4.0, 0.0 | 0.377 | 8.15 / 10.05 | 0.24 |
| `bb90-g2` | baseline | 1 / 3 | —, —, 3.37 | 4, 6, 3 | 0.241 | 0.0, +2.0, +2.0 | 0.267 | 3.28 / 3.73 | 0.09 |
| `bb90-g2` | nopot | **0 / 3** | —, —, — | 3, 3, 3 | 0.234 | +2.0, +1.0, +1.0 | 0.240 | 2.51 / 3.76 | 0.11 |
| `bb90-g2` | recyc8 | 1 / 3 | —, 3.39, — | 5, 3, 3 | 0.225 | 0.0, +1.0, +1.0 | 0.277 | 14.71 / 21.09 | 0.39 |
| `bb50-g1` | baseline | 1 / 3 | —, —, 2.98 | 8, 6, 5 | 0.225 | −1.0, −1.0, −1.0 | 0.427 | 14.44 / 20.37 | 0.33 |
| `bb50-g1` | **nopot** | **3 / 3** | 2.69, 2.47, 2.65 | 4, 5, 5 | 0.172 | 0.0, 0.0, −1.0 | 0.263 | **1.38 / 1.65** | 0.10 |
| `bb50-g1` | **recyc8** | **3 / 3** | 2.77, 2.73, 2.77 | 7, 7, 7 | 0.219 | −2.0, **−3.0**, −2.0 | 0.490 | **0.70 / 1.05** | 0.06 |
| `bb50-g1` | **combo** | **3 / 3** | 2.79, 3.11, 2.77 | 7, 4, 6 | 0.226 | −2.0, −1.0, −2.0 | 0.384 | **2.95 / 3.76** | 0.11 |
| `bb25` | baseline | 2 / 3 | 2.72, —, 3.25 | 15, 10, 15 | 0.226 | +1.0, −2.0, +1.0 | 0.427 | 7.27 / 10.36 | 0.36 |
| `bb25` | nopot | 2 / 3 | 2.22, —, 2.70 | 13, 8, 16 | 0.227 | +1.0, −1.0, +1.0 | 0.404 | 4.09 / 5.67 | 0.39 |
| `bb25` | **recyc8** | **3 / 3** | 2.12, 2.52, 2.26 | 15, 12, 13 | 0.281 | +1.0, +2.0, +1.0 | 0.409 | **1.96 / 2.54** | 0.10 |
| `bb25` | combo | 2 / 3 | 2.26, 2.34, — | 13, 13, 1 | 0.257 | +1.0, +1.0, +1.0 | 0.431 | 7.40 / 10.83 | 0.30 |

\* ipTM is a **pose-quality read-out only** — never a gap, never a tie-breaker.
The ipTM `specificity_gap` and the paired same-pose MM-GBSA ddG both failed the
7OW6 positive control and are retired; neither is used anywhere above.

## Verdict 1 — yes, the spread collapses, and it is design-specific

**This is a positive result, and the first one the pose problem has produced.**

| design | baseline mean / max (Å) | best condition | best mean / max (Å) | collapse? |
| --- | --- | --- | --- | --- |
| `bb50-g1` | 14.44 / 20.37 | **recyc8** | **0.70 / 1.05** | **yes — 20× ; two unrelated docked modes become one** |
| `bb25` | 7.27 / 10.36 | **recyc8** | **1.96 / 2.54** | **yes — 3.7×** |
| `bb135-g1` | 4.54 / 5.81 | **combo** | **2.04 / 2.97** | **yes — 2.2×** |
| `bb90-g2` | 3.28 / 3.73 | nopot | 2.51 / 3.76 | marginal (already tight) |
| `bb173-g1` | 8.78 / 11.79 | recyc8 | 8.15 / 10.05 | **no — spread persists in every condition** |

Three of five designs drop out of the two-mode regime entirely. `bb50-g1` is
the extreme case: at baseline its four poses sat in two docked modes ~20 Å
apart; at `recycling_steps` 8 all three draws agree to **0.70 Å mean / 1.05 Å
max**, with a frame RMSD of 0.06 Å, and it holds **7 Asp6 side-chain contacts in
every draw** with the donor distance pinned at 2.73–2.77 Å. That is a single
pose, not a coin flip.

**Neither setting wins outright — they differ by design.** `recyc8` is best for
`bb50-g1` and `bb25`; `without_potentials: true` is best for `bb135-g1` and is
actively *worse* for `bb173-g1` and `bb90-g2`. The combination is not additive:
it rescued `bb135-g1` (2.04 Å, 3/3) but *re-opened* the spread for `bb25`
(7.40 Å, back to 2/3). So "more recycling" is not a universal fix, and
`without_potentials` is not either; what the experiment establishes is that
**the spread is a property of the sampling settings, not an irreducible
property of these designs** — which is exactly what was in doubt.

A caution that the data itself supplies: a collapsed spread is **a more
reproducible prediction, not evidence of a better binder**. Nothing here
validates Boltz-2 against experiment, and the same settings that pin
`bb50-g1`'s pose also lowered its ipTM in one condition and raised it in
another. `bb90-g2` remains the standing counter-example in the other direction:
its pose was already stable at 3.28 Å and it still reads Asp6 on 0–1 draws of 3.
**A stable pose is necessary, not sufficient.**

## Verdict 2 — three designs now clear the Asp6 gate 3 of 3

| design | conditions reaching 3/3 | draws passed |
| --- | --- | --- |
| **`bb50-g1`** | **nopot, recyc8 and combo — all three** | **9 / 9** |
| **`bb135-g1`** | **nopot and combo** | **6 / 6** |
| **`bb25`** | **recyc8** | **3 / 3** |
| `bb173-g1` | none (best 2/3) | — |
| `bb90-g2` | none (best 1/3) | — |

Batch-05 ended with **zero** designs clearing the gate in 3 of 3. Under
alternative sampling settings, three do. `bb50-g1` is the strongest result: it
passes 9 draws out of 9 across three different sampling settings, i.e. its
read-out is not an artefact of one setting.

**What this does and does not mean.** The gate here is the Asp6 read-out gate
replicated 3 of 3, the same gate batch-05 used. It is *not* the full hard-gate
set of `docs/PIPELINE.md` section 6. Evaluated against that full set with
`src/design_gates.evaluate()` and thresholds untouched:

* **every** design × condition still fails `peptide_buried_frac ≥ 0.40` — the
  best cell in this study is `bb135-g1` at 0.33, and the floor also rejects the
  real 7OW6 TCR (0.224 predicted / 0.207 crystal). The recalibration to 0.18
  proposed in `results/validation/CONTROL_7OW6.md` remains **unapplied and is
  the user's call.** Under 0.18, `bb135-g1`, `bb25` and `bb50-g1`/recyc8+combo
  would clear it; `bb50-g1`/**nopot** would *still* fail at 0.172.
* `bb50-g1` under **recyc8** breaches the net interface charge floor on one
  draw (−3.0 against a −2.0 minimum) and sits exactly on it in the other two.
  Its `combo` and `nopot` draws sit at −2.0/−1.0/0.0 and do not.

Nothing above was adjusted to produce survivors. These are reported values.

## Honest limitations

1. **n = 3 draws per cell, 5 designs.** `bb50-g1` at 0.70 Å across three draws
   is a strong signal; `bb135-g1` at 2.04 vs 4.54 Å is suggestive, not
   conclusive. A confirmation run at 5–10 draws on `bb50-g1` + `bb25` would
   settle it for ~30 more calls and zero design slots.
2. **No off-target panel was run.** Every call here is the G12D state only, as
   in batches 04 and 05. These designs have **not** been shown to discriminate
   G12D from WT, apo or the decoy under the new settings — and the project has
   no validated ranking quantity for that (both candidates failed the 7OW6
   control). Selectivity remains unmeasured.
3. **Only one sampling axis was moved per condition**, at one value each
   (`recycling_steps` 8). `sampling_steps`, `step_scale` and
   `diffusion_samples` were deliberately held at the batch-05 values so the
   comparison stays clean.
4. **A collapsed spread is self-consistency, not accuracy.**

## Recommendation

1. **Confirm `bb50-g1` and `bb25` at 5–10 draws** under their best condition
   (`recyc8`) before either is treated as a survivor. ~30 calls, 0 design slots.
2. **Run the off-target panel** (WT / apo / decoy) on the 3-of-3 designs under
   the settings that stabilised them. Until that is done, no design in this
   project has been shown to be *selective*, only to *read the aspartate*.
3. **Adopt the per-design setting sweep as standard practice** for future
   batches' replication gate: the sampling setting is now a measured variable,
   and running the gate at a single setting under-counts designs that are stable
   at another.
4. **Settle the `peptide_buried_frac` floor** (0.40 vs the control-calibrated
   0.18) — it is now the binding constraint on every design in this study, and
   it should not be decided under deadline pressure.
5. **`bb173-g1` is the one genuine null** — no setting helps it. Its failure
   mode is the design, not the sampling; the rigidifying redesign proposed in
   `results/batch-05/SUMMARY.md` recommendation 2 is still the right move.

**Submission:** unchanged from batch-05 — this study does not make the call.
Track-1 stands at **0 of 40 design slots used**; nothing has been submitted, no
PR merged, nothing pushed to `main`. The window opens 07:30 Taipei on Nov 4 and
the deadline is Nov 8 (AoE). The submission decision is the user's alone.
