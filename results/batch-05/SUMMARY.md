# Batch-05 — strategy `B5-c6-62` at volume, with replication folded into the gate

**Date:** 2026-11-03 (Asia/Taipei).
**Scripts:** `scripts/run_batch05.py`, `scripts/rank_batch05.py`,
`scripts/replicate_batch05.py` (which re-use `run_batch04.py`,
`rank_batch04.py`, `replicate_batch04.py`, `replicate_rmsd_batch04.py`
verbatim with their output directory pointed at `results/batch-05`).
**Raw data:** `checkpoint.jsonl` (authoritative, 60 rows), `candidates.json`,
`backbone_poses.json`, `run_stats.json`, `gate_results.json`, `designs.csv`,
`replicates.json`, `replicate_rmsd.json`, `final_ranking.json`,
`survivors.csv`, `scored_structures/`, `replicate_structures/`.
**Boltz-2 calls:** 60 screening (one constrained G12D call per candidate) + 12
replicate calls = **72**. Strictly serial throughout, ≥ 8 s interval.
**Track-1 design slots consumed: 0 of 40.** Nothing submitted.

## Headline

**Zero designs pass the batch-05 gate.** Four of sixty cleared the Asp6 direct
read-out gate on the screening draw; none of the four cleared it in 3 of 3
independent constrained draws. This is an honest negative and no gate was
loosened to avoid it.

## What was run

This is not a new pipeline. It is batch-04's code with the three changes
`results/batch-04/SUMMARY.md` recommended:

1. **One strategy, at volume.** `B5-c6-62` is `B4-c6-62` verbatim — 62 aa
   binder, hotspot `C6` alone, ProteinMPNN `sampling_temp` 0.15 — the only
   batch-04 strategy that produced a gate survivor whose pose replicated.
2. **A deterministic donor-anchored graft for every backbone** and every MPNN
   parent selected on it, rather than one graft off the single best sample.
   MPNN-selected sequences are not scored on their own (0 of 38 passed in
   batch-04).
3. **Replication inside the gate**, not after it.

Scoring unchanged: one constrained Boltz-2 call per candidate with the
validated flat pocket term
`{"binder":"D","contacts":[{"id":"C","residue_index":6}],"max_distance":6.0}`.

| stage | count |
| --- | --- |
| RFdiffusion backbones generated | **200** |
| backbones whose pose reaches peptide position 6 (≥ 2 contacts, ≤ 5.0 Å) | **46** |
| unique donor-anchored grafts built | **60** (32 `-g1`, 28 `-g2`, over 32 backbones) |
| complexes scored with constrained Boltz-2 | **60** |
| screening-draw Asp6 gate passes | **4** |
| **gate passes in 3 of 3 replicate draws** | **0** |

Run time 2801 s for design + screening; the 7200 s scoring budget was never
reached and no call cap was hit, so the batch is complete, not truncated.

## Result 1 — the chemistry held up at volume

The donor graft keeps doing what batch-04 showed it does:

| quantity | batches 02+03 (88) | batch-04 (54) | **batch-05 (60)** |
| --- | --- | --- | --- |
| any donor within 3.5 Å of Asp6 OD1/OD2 | 0 | 2 | **4** |
| `asp6_buried_frac` > 0 | 5 / 88 | 23 / 54 | **39 / 60** |
| max `asp6_buried_frac` | 0.538 | 0.888 | **0.966** |
| designs with ≥ 1 Asp6 **side-chain** contact | — | 24 / 54 | **36 / 60** |
| max `peptide_buried_frac` | — | 0.353 | 0.341 |
| ipTM median / max | — | 0.348 / 0.708 | 0.410 / 0.845 |

39 of 60 bury the aspartate at all and 36 make a side-chain contact to it —
the highest rates this project has reached. Scaling the strategy scaled the
chemistry roughly linearly: 4 donor-in-range designs from 60 candidates versus
2 from 54.

## Result 2 — interface charge stayed fixed

| | median net interface charge | mean | min / max | negative : ≥ 0 |
| --- | --- | --- | --- | --- |
| batch-02 | −1.0 | — | — | — |
| batch-03 | −1.0 | — | — | — |
| batch-04 (54) | **0.0** | −0.5 | −4.0 / +3.0 | 23 : 31 |
| **batch-05 (60)** | **0.0** | **+0.39** | −5.0 / +5.1 | **22 : 38** |

The donor bias's charge correction replicates at volume. 38 of 60 designs sit
at or above neutral, against a −1.0 median in both pre-donor-bias batches, and
the mean has moved further positive (+0.39 vs −0.5) now that every backbone
gets a graft. **This is the one unambiguous success of batch-05.**

## Result 3 — the gate, and why nothing survives it

Four designs cleared the Asp6 donor gate on the screening draw:

| screening rank | design | `asp6_min_donor_distance` | partner | `asp6_buried_frac` | `peptide_buried_frac` | Asp6 side-chain contacts | net charge | ipTM* |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | `B5-c6-62_bb135` (`-g1`) | **2.80 Å** | ARG23:NH1 | 0.500 | 0.322 | 4 | 0.0 | 0.437 |
| 2 | `B5-c6-62_bb173` (`-g1`) | **3.40 Å** | ARG17:NH2 | 0.747 | 0.286 | 6 | −1.0 | 0.400 |
| 3 | `B5-c6-62_bb90` (`-g2`) | **3.14 Å** | SER22:OG | 0.320 | 0.240 | 3 | +1.0 | 0.253 |
| 4 | `B5-c6-62_bb50` (`-g1`) | **2.96 Å** | ARG44:NH2 | 0.363 | 0.188 | 5 | 0.0 | 0.294 |

\* ipTM is a pose-quality read-out only — never a gap, never a tie-breaker.

Each was then re-run in **3 independent constrained draws** (the screening call
is reported as draw 0 and does not count towards the 3):

| design | **gate passes / 3** | `asp6_min_donor_distance` over draws | `asp6_buried_frac` (mean ± SD) | `peptide_buried_frac` (mean ± SD) | Asp6 side-chain contacts | net charge | ipTM (mean ± SD) |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `bb135` `-g1` | **2 / 3** | 2.03, —, 2.03 | 0.327 ± 0.117 | **0.333** ± 0.005 | 5, 1, 5 | 0.0, +1.0, 0.0 | 0.451 ± 0.058 |
| `bb173` `-g1` | **2 / 3** | 2.48, 2.67, — | **0.626** ± 0.206 | 0.297 ± 0.019 | **6, 6, 6** | −1.0, +1.0, 0.0 | 0.531 ± 0.224 |
| `bb90` `-g2` | **1 / 3** | —, —, 3.37 | 0.524 ± 0.096 | 0.241 ± 0.001 | 4, 6, 3 | 0.0, +2.0, +2.0 | 0.267 ± 0.024 |
| `bb50` `-g1` | **1 / 3** | —, —, 2.98 | 0.449 ± 0.151 | 0.225 ± 0.044 | 8, 6, 5 | −1.0, −1.0, −1.0 | 0.427 ± 0.070 |

**`survivors_3of3` = 0.** The gate requires 3 of 3; the best designs manage 2.

Worth recording about the two near-misses: when `bb135` and `bb173` *do* read
the carboxylate they read it **tightly** — 2.03 Å twice with zero spread for
`bb135`, 2.48/2.67 Å for `bb173`, all shorter than either batch-04 survivor
(2.94 / 2.70 Å). `bb173` also holds **6 Asp6 side-chain contacts in every
single draw** (SD 0.0). The failure is not a weak contact; it is an
intermittent one.

## Result 4 — replicate pose stability

Binder C-alpha RMSD measured in the pMHC frame (HLA + peptide superposed,
binder **not** re-fitted). Frame RMSD 0.08–0.33 Å throughout, so these numbers
are real. "Replicates only" excludes the screening pose.

| design | pairwise mean / max, all 4 poses | pairwise mean / max, 3 replicates only | pattern |
| --- | --- | --- | --- |
| `bb90` `-g2` | **3.66 Å / 4.11 Å** | **3.28 Å / 3.73 Å** | tightest spread in the batch; every pose within ~4 Å |
| `bb135` `-g1` | 3.65 Å / 5.81 Å | 4.54 Å / 5.81 Å | two modes ~5.6–5.8 Å apart; the off-mode draw is the gate failure |
| `bb173` `-g1` | 7.29 Å / 11.79 Å | 8.78 Å / 11.79 Å | r0/r1 agree to 2.89 Å and both pass; r2 flips 11.7 Å away and fails |
| `bb50` `-g1` | **14.26 Å / 20.37 Å** | 14.44 Å / 20.37 Å | r0/r1 agree to 2.68 Å, the screening pose and r2 sit ~20 Å away — two unrelated docked modes |

The mechanism is exactly the one `CONTROL_7OW6.md` predicted and batch-04's
`bb26` demonstrated: with `constraints` and `without_potentials: false`, a weak
binder's landscape holds several docked modes and the sampling flips between
them. In every one of the four cases the gate failure and the pose flip are the
**same event** — `bb135` r1, `bb173` r2 and `bb50` r0/r1 each lose the
carboxylate contact precisely on the draws that land in the other mode.
`bb90` is the instructive exception: its pose is stable to 3.3 Å yet it still
only reads Asp6 on one draw in three, i.e. a stable pose is necessary but not
sufficient.

No batch-05 design reproduces as well as batch-04's `bb25` (5.67 Å mean, 3 of 4
passes), which therefore remains the project's best design.

## The `peptide_buried_frac` floor — reported, deliberately not applied

**No design in batch-05 reaches the 0.40 floor.** Batch max is 0.341
(`bb54 -g1`); the median is 0.209.

| design | `peptide_buried_frac` (replicate mean) | clears 0.40? | clears the proposed 0.18? |
| --- | --- | --- | --- |
| `bb135` `-g1` | 0.333 | **no** | yes |
| `bb173` `-g1` | 0.297 | **no** | yes |
| `bb90` `-g2` | 0.241 | **no** | yes |
| `bb50` `-g1` | 0.225 | **no** | yes |
| *7OW6 real TCR (reference)* | *0.224 predicted / 0.207 crystal* | ***no*** | *yes* |

Across all 60 screened designs: **0 clear 0.40**, 42 clear 0.18.

The floor rejects the real, experimentally validated 7OW6 TCR itself. The
recalibration to 0.18 proposed in `CONTROL_7OW6.md` remains **unapplied and is
the user's call**; `src/design_gates.py` is unchanged. Note that this threshold
is not what decided batch-05 — the 3-of-3 replication gate is, and recalibrating
the floor to 0.18 would still leave **zero** survivors.

## Verdict

* **The design half keeps improving.** 200 backbones, 46 reaching position 6,
  60 grafts, 39 burying Asp6, 36 contacting its side chain, charge median held
  at 0.0 with the mean now positive, and the two best donor contacts in project
  history at 2.03 Å. Scaling `B5-c6-62` did what it was expected to do.
* **The scoring half still cannot hold a pose.** Every gate failure in this
  batch is a sampling flip between docked modes, not a chemistry failure. The
  3-of-3 gate is correctly rejecting designs whose read-out is a coin flip —
  this is the gate working, not the gate being too strict.
* **0 of 60 designs are submittable.** Batch-04's `bb25` (3 of 4 draws, 5.67 Å
  mean RMSD) is still the only design in the project that reads the G12D
  aspartate reproducibly, and batch-05 did not beat it.

## Recommendation

**Do not submit when the window opens at 07:30 on Nov 4.** The Track-1 cap is
40 designs with **0 spent** and the deadline is Nov 8 (AoE) — four days of
runway. There is no scarcity argument for spending a slot on a design whose
carboxylate contact appears in 2 draws out of 3, and a 1-design submission
built on batch-04's `bb25` alone would waste the ≤ 20-per-submission slot.

The submission decision is the user's alone. Nothing has been submitted, no PR
merged, nothing pushed to `main`.

**Highest-value next steps, in order:**

1. **Attack the pose flip directly — it is now the binding constraint.** Three
   batches of evidence say the chemistry is solvable and the pose is not.
   Concretely: test `without_potentials: true` and/or more `recycling_steps` on
   the four near-misses and on `bb25`, and measure whether the RMSD spread
   collapses. This costs ~15 calls and no design slots, and it is the only
   experiment that could convert existing designs into survivors.
2. **Rescue `bb173` and `bb135` by sequence, not by threshold.** Both hold a
   genuine short Arg···Asp6 bond in the mode they reach 2 draws in 3;
   `bb173` keeps 6 Asp6 side-chain contacts in *every* draw. A small
   rigidifying redesign around the graft (proline/disulfide placement, loop
   shortening) targets the actual failure mode.
3. **Settle the `peptide_buried_frac` floor** — 0.40 or the control-calibrated
   0.18. It did not change this batch's outcome, but it will block any future
   survivor and it should not be decided under deadline pressure.
4. **Keep `bb25` as the held reference design to beat.**
