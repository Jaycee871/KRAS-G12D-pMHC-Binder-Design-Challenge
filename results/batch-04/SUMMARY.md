# Batch-04 - not run, and why

**Date:** 2026-11-03 (Asia/Taipei). **Boltz-2 calls spent: 20** (the positive
control and its replicates), **0** on new designs. No RFdiffusion, no
ProteinMPNN.

Batch-04 was planned as a fresh constrained design batch built on
`S1-pos6-68` / `S1-pos6-65`, with ProteinMPNN biased toward hydrogen-bond
donors near peptide position 6. It was gated on one thing agreed in advance:
the pocket-constrained Boltz-2 protocol had to reproduce known-correct
behaviour on the real 7OW6 TCR first. **It did not**, so the batch was not
launched. Full control write-up and raw numbers:
[`results/validation/CONTROL_7OW6.md`](../validation/CONTROL_7OW6.md).

## The control, in one table

| | result |
| --- | --- |
| Pose, TCR α/β variable domains, C-alpha RMSD to crystal in the pMHC frame | **1.29 Å / 1.51 Å** — near-native |
| Buried area on the peptide, predicted vs crystal | 184 Å² vs 169 Å² |
| Peptide buried fraction, predicted vs crystal | 0.224 vs 0.207 |
| Position-6 contacts, predicted vs crystal | 7 vs 7 |
| ipTM(G12D) | 0.9268 ± 0.0010 |
| ipTM(WT) | **0.9325 ± 0.0001** |
| ipTM(decoy `AVFDRKSDAK`) | 0.9136 ± 0.0005 |
| ipTM(apo) | 0.7844 |
| **specificity_gap of the known selective TCR** | **−0.006** |

The geometry half passes beautifully. The ranking half fails: the one binder in
the world we *know* is G12D-selective is scored by Boltz-2 as marginally
preferring the wild-type peptide, reproducibly (SD 0.001), and as barely
distinguishing an unrelated HLA-A\*11:01 decoy peptide.

## And the gap does not reproduce on designs either

Identical four-state constrained panels, re-run for batch-03's best three:

| design | gap, first pass | gap, replicate | ipTM(G12D) first → replicate |
| --- | --- | --- | --- |
| `S1-pos6-68_bb18` | +0.298 | **−0.300** | 0.708 → 0.440 |
| `S1-pos6-68_bb19` | +0.221 | **−0.041** | 0.680 → 0.417 |
| `S1-pos6-68_bb13` | +0.183 | +0.151 | 0.758 → 0.758 |

Three further repeats of the same constrained G12D request gave ipTM SD 0.056
(`bb18`: 0.486 / 0.619 / 0.524) and 0.227 (`bb13`: 0.277 / 0.758 / 0.758).
Sampling noise on a design is 0.06-0.23; the best specificity gap this project
has ever recorded is 0.30. **The leaderboards of batches 01-03 are, to a large
extent, rankings of noise.**

## Why that means no batch-04 tonight

A batch-04 sized to finish by 06:30 would have been roughly 60-80 candidates
and ~200 Boltz-2 calls, all of them ranked by `specificity_gap`. Spending 14
hours to produce another ranking by a quantity that

* carries no mutation-specificity signal on the only ground truth available, and
* is not reproducible across identical calls in the regime where our designs live,

would have produced a longer leaderboard with the same information content as
batch-03's: none. The honest move is to fix the read-out first. Batches 01, 02
and 03 each failed on the pose; batch-04 was stopped because the *score* is
the thing that is broken.

## What the real TCR scores on our hard gates

| hard gate | threshold | real 7OW6 TCR (constrained prediction) | verdict |
| --- | --- | --- | --- |
| peptide share of buried area | ≥ 0.40 | 0.224 (crystal 0.207) | **FAIL** |
| Asp6 read-out | H-bond/salt bridge or ≥ 7 pos-6 contacts | no bond, 7 contacts | pass (exactly on the escape clause) |
| position-6 contacts | ≥ 3 | 7 | pass |
| net interface charge | −2.0 … +4.0 | −2.9 (crystal −1.9) | **FAIL** |
| PDE-interface(G12D) | < 10 Å | 0.40 | pass |
| predicted-pose filter | ≥ 25 Å², ≥ 1 pos-6 contact | 184 Å², 7 | pass |

A real, experimentally validated selective binder fails two of our six hard
gates. Proposed recalibration with the evidence for each number is in
`CONTROL_7OW6.md`; it is **not applied to `src/design_gates.py`**, because
step 2 of the plan was explicitly conditional on the control passing and this
is the user's call to make. The two substantive proposals are
`PEPTIDE_BURIED_FRACTION_MIN` 0.40 → 0.18 and `NET_INTERFACE_CHARGE_MIN`
−2.0 → −3.5, both anchored on the control rather than chosen to let designs
through.

## The one design that survives everything

Under the control-calibrated thresholds, exactly one batch-03 design clears
every hard gate, and it is also the only one of the top three whose score
vector reproduces:

**`S1-pos6-68_bb13`** (68 aa, strategy `S1-pos6-68`, constrained rescore)

| quantity | value |
| --- | --- |
| ipTM(G12D) | 0.758 (0.758 on replicate; 0.277 / 0.758 / 0.758 over three further repeats) |
| ipTM(WT) | 0.575 → 0.398 on replicate |
| ipTM(apo) | 0.259 |
| ipTM(decoy) | 0.490 → 0.607 on replicate |
| specificity gap | +0.183 first pass, +0.151 replicate |
| peptide buried fraction | 0.223 (control: 0.224) |
| buried area on the peptide | 235 Å² (control: 184 Å²) |
| buried area, total | 1055 Å² |
| position-6 contacts | **12** (control: 7) |
| Asp6 hydrogen bond | **yes** |
| Asp6 salt bridge | **yes** |
| net interface charge | −1.0 |
| PDE-interface | < 10 Å |
| pose filter | pass |

It is the only design in three batches that both reads the Asp6 carboxylate
chemically (H-bond *and* salt bridge, which the real TCR does not even do) and
buries as large a share of the peptide as the real TCR does, with a gap that
survived a replicate. That is a weak but genuine recommendation, resting on
geometry rather than on the discredited ipTM gap.

Runners-up on geometry alone: `S1-pos6-65_bb03` (peptide buried 0.190, 8
position-6 contacts, Asp6 H-bond and salt bridge, charge 0.0, gap +0.045).
`bb18` and `bb19` should be dropped despite their headline gaps: both have
weak position-6 footprints (4 and 5 contacts), no Asp6 bond, and their gaps
inverted on replicate.

## What to try next, in order

1. **Replace the ranking quantity with a paired, same-pose ΔΔG.** Take the
   predicted G12D complex, mutate peptide position 6 D→G *in that structure*,
   repack locally and re-score. One structure, one perturbation, no re-docking
   - which removes the docking noise that is currently drowning the signal.
   FoldX, Rosetta ddG or an MM-GBSA rescore all run on the box for free and
   cost zero API budget. This is the single highest-value change available.
2. **Make replication mandatory.** Any score that ranks designs gets ≥ 3
   identical calls and a reported SD. Budget `4 states × 3 repeats` per
   finalist, not per candidate: triage on geometry, replicate only the top ~10.
3. **Recalibrate the two mis-set gates** (`peptide_buried_frac` 0.40 → 0.18,
   charge floor −2.0 → −3.5), since a known-good binder fails both.
4. **Then, and only then, run a design batch.** With a reproducible score and
   calibrated gates, the `S1-pos6-68` strategy plus a fixed Arg anchor at the
   position-6 contact and donor-biased ProteinMPNN is still the right design
   plan - it produced `bb13`. It is the scoring, not the designing, that was
   blocking.
