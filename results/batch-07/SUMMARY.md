# Batch-07 - sampling sweep on bb54-g2, and a rigidified redesign around its donor graft

Branch `results/batch-07-rigidified-graft`, stacked on `results/batch-06-peptide-centred`.
Target: KRAS G12D decamer **VVVGADGVGK** on HLA-A\*11:01 (7OW6). WT comparator
**VVVGAGGVGK**, decoy **AVFDRKSDAK**, apo = peptide-free HLA + b2m.

Nothing in `src/design_gates.py` was touched. No gate or threshold was loosened,
retuned or re-derived. The WT-loss criterion is `scripts/batch06_panel.WT_CRITERION`
re-used **verbatim**, fixed before any batch-07 result was seen. No ipTM
`specificity_gap` and no paired MM-GBSA `ddG` was computed anywhere - both stay retired;
ipTM appears only as a pose-quality read-out. **72 Boltz-2 calls** (216 structures),
every one strictly serial at >= 8 s spacing. **0 of the 40 Track-1 design slots used**
(nothing submitted).

---

## Step 1 - the per-design sampling sweep on bb54-g2 (17 calls)

`B6-c6-46_w2_bb54_donor-anchored-g2`, 46 aa, one Boltz-2 call per cell with
`diffusion_samples: 3` (three independent trajectories, so they cannot be byte-identical
to one another). Same flat pocket constraint
`{"binder":"D","contacts":[{"id":"C","residue_index":6}],"max_distance":6.0}`,
`step_scale` 1.638, scoring path `run_batch04.score_candidate` unchanged.

| cell | recycling | sampling | without_potentials | Asp6 gate | distinct | binder CA RMSD mean/max (frame floor) | median pep area | median pbf |
|---|---|---|---|---|---|---|---|---|
| r1 | 1 | 50 | no | **3/3** | 3/3 | 11.13 / 13.98 (0.55) | 221.9 | 0.303 |
| r1_rep | 1 | 50 | no | 2/3 | 3/3 | 13.29 / 15.78 (0.36) | 180.4 | 0.228 |
| r1np | 1 | 50 | yes | 1/3 | 3/3 | 16.17 / 17.53 (0.35) | 163.2 | 0.193 |
| r3 | 3 | 50 | no | **3/3** | 3/3 | 12.95 / 16.84 (0.36) | 184.5 | 0.269 |
| r3_rep | 3 | 50 | no | 2/3 | 3/3 | 16.41 / 22.16 (0.33) | 210.8 | 0.300 |
| **r3np** | 3 | 50 | **yes** | **3/3** | 3/3 | 13.93 / 16.97 (0.39) | 180.0 | 0.291 |
| **r3np_rep** | 3 | 50 | **yes** | **3/3** | 3/3 | 13.02 / 17.37 (0.33) | 186.1 | 0.286 |
| r5 | 5 | 50 | no | 1/3 | 3/3 | 13.50 / 17.38 (0.38) | 177.8 | 0.267 |
| r5np | 5 | 50 | yes | 2/3 | 3/3 | 12.84 / 16.68 (0.40) | 177.5 | 0.283 |
| r6 | 6 | 50 | no | 1/3 | 3/3 | 3.51 / 4.59 (0.43) | 276.9 | 0.307 |
| r6np | 6 | 50 | yes | 2/3 | 3/3 | 12.23 / 16.69 (0.38) | 189.5 | 0.239 |
| r8 (= batch-06 `recyc8`) | 8 | 50 | no | 2/3 | 3/3 | 7.17 / 7.76 (0.29) | 288.4 | 0.319 |
| r8np (= batch-06 `combo`) | 8 | 50 | yes | 0/3 | 3/3 | 3.87 / 5.05 (0.44) | 278.5 | 0.304 |
| r10 | 10 | 50 | no | 2/3 | 3/3 | 4.93 / 6.99 (0.36) | 245.5 | 0.300 |
| r10np | 10 | 50 | yes | 2/3 | 3/3 | **2.25 / 2.81** (0.28) | 274.2 | 0.302 |
| r8s100 | 8 | 100 | no | 2/3 | 3/3 | 12.67 / 18.44 (0.41) | 278.2 | 0.320 |
| r10s100np | 10 | 100 | yes | 0/3 | 3/3 | 13.69 / 19.80 (0.48) | 214.1 | 0.273 |

**Result 1 - bb54-g2 does reach 3 of 3, and the setting that does it is the opposite of
the batch-06 intuition.** `r3np` (recycling_steps **3** + `without_potentials`) is the only
cell that is 3-of-3 **twice**, on two independent calls (6/6 draws). `r1` and `r3` each hit
3/3 once and fell back to 2/3 on an independent repeat, so those are lucky draws, not
settings. More recycling (8-10) and longer trajectories do **not** help the gate: `r8np` and
`r10s100np` are 0/3. Across three further independent `r3np` calls run later in this batch
the pooled count is **8 of 9 draws** gated.

**Result 2 - the gate and the pose are anti-correlated.** The tightest pose in the whole
sweep is `r10np` at 2.25/2.81 A (frame floor 0.28) - and it gates only 2/3. The cells that
gate 3/3 are loose: `r3np` sits at 13.93/16.97 A. Pose tightness and donor-contact
reproducibility are being bought from different parts of the sampler. The pose under the
3-of-3 setting is **not** tight; bb50-g1 (1.58/3.78 A at n = 10) remains the pose reference.

## Step 2 - rigidifying redesign around the donor graft

First measurement, before designing: bb54-g2's donor anchor is **already on a helix**
(anchor residue 35, `shell_ss` = `HHHH`, 63 % helix, donor-loop length 0, longest loop 6).
There is no floppy donor loop to stiffen. The redesign therefore targeted the *connectors
and termini* that let the donor-bearing helix rock, never the donor shell itself: the 6
nearest Asp6-shell positions are protected from every edit, so the donor geometry that
works is unchanged by construction (`scripts/rigidify.py`).

Deterministic stiffening edits: Gly -> Ala on negative-phi loop positions, Pro at loop
positions with phi in (-95, -45), optional Cys-Cys staple on CB-CB pairs 3.3-6.5 A apart and
>= 5 residues apart. **No staple candidate existed on bb54's backbone** at either the strict
(3.3-5.2 A) or the relaxed window, so no `rigid-gps` variant was ever built - reported, not
worked around. Positive-phi glycines were kept (they are there for a reason).

Fresh scaffolds came from the existing machinery: RFdiffusion `/v1/biology/ipd/rfdiffusion/generate`
+ ProteinMPNN `/v1/biology/ipd/proteinmpnn/predict`, three strategies at bb54-g2's paratope
size or larger (never smaller) - B7-rigid-46 / -52 / -58, hotspot C6 - filtered by batch-06's
peptide-centred pose filter **plus** a rigidity filter (`ss_frac >= 0.55`, donor loop <= 3,
longest loop <= 8).

- Wave 1 at `diffusion_steps 30`: **0 of 96 backbones kept** (closest position-6 approach
  4.66 A, median 18.4 A). The diffusion-steps change, not the filter, was the cause.
- Waves 2-3 back at `diffusion_steps 15` (the batch-04/05/06 value): 4 backbones kept of 192.
  One kept backbone (`B7-rigid-52_bb00`) yielded 0 usable sequences.

**21 candidates scored**, each as a paired set wherever possible - a `graft` control and the
`rigid-gp` rigidified twin on the **same backbone with the same donor graft** - so the
rigidification effect is measurable rather than confounded. Track A (bb54's own backbone,
including bb54-g2 self-rigidified as `bb54_rigid-gp-self`): 15. Track B (fresh rigid
scaffolds): 6. 42 Boltz-2 calls (`r3np` first, `r3` only when `r3np` missed 3/3).

## Step 3 - gate and panel

**Asp6 gate, 3 of 3, among the 21 new designs: 1.**

| design | track | best cell | Asp6 gate | notes |
|---|---|---|---|---|
| `bb54_rigid-gp-p2` | A | r3 | **3/3** | repeat at the same cell: **1/3**; `r3np` 1/3; `recyc8` 0/3 |
| `bb54_graft-p3` | A | r3 | 2/3 | unrigidified control |
| `bb54_rigid-gp-p3` | A | r3np | 2/3 | |
| `bb54_rigid-gp-self` | A | r3np / r3 | 1/3 | bb54-g2's own sequence + 3 Pro + 1 Gly->Ala |
| all 6 Track-B designs | B | both | **0/3** | fresh rigid scaffolds never engaged Asp6 |
| remaining 11 | A | both | 0-1/3 | |

**Result 3 - rigidification did not make the contact constitutive.** The single 3-of-3
(`bb54_rigid-gp-p2`) did not reproduce: an independent repeat under the identical cell gave
1/3. Its paired unrigidified control behaved the same way, so the pass is sampling noise, not
rigidification. Self-rigidifying bb54-g2 itself made it **worse** (1/3 vs the parent's 3/3 at
r3np) - the prolines that stiffen the connectors also cost the donor helix its register.
Track B is a clean 0/6: fresh backbones filtered for rigidity simply do not find Asp6.

### WT arm - `batch06_panel.WT_CRITERION`, verbatim

Criterion (fixed in advance, not retuned): median `buried_area_peptide` >= 30 % lower on WT
**and** median position-6 heavy-atom contacts (all atoms, backbone included) >= 30 % lower,
over 3 draws per state in the **same** cell. Asp6 side-chain loss on WT is excluded - Gly has
no side chain, so it is true by construction and carries no information.

| design | cell | G12D gate | pep-area loss | pos6-contact loss | verdict |
|---|---|---|---|---|---|
| `bb54_rigid-gp-p2` | r3 | 3/3 | -17.5 % | **+66.7 % (gains contacts)** | FAIL |
| `bb54_rigid-gp-p2` | r3 repeat | 1/3 | -4.3 % | +40 % | FAIL |
| `bb54_rigid-gp-p2` | recyc8 | 0/3 | -21.2 % | +100 % | FAIL |
| bb54-g2 (parent) | r3np call 1 | 3/3 | -27.9 % | +25 % | FAIL |
| bb54-g2 | r3np call 2 | - | -33.8 % | -50.0 % | PASS |
| bb54-g2 | r3np call 3 | 3/3 | -20.3 % | -42.1 % | FAIL (area) |
| bb54-g2 | r3np call 4 | 2/3 | -40.7 % | -58.3 % | PASS |
| bb54-g2 pooled, 9 G12D vs 12 WT draws, distinct calls only (context, not the criterion) | r3np | 8/9 | **-33.4 %** | **-52.4 %** | PASS, no strict separation |

**Result 4 - no design clears the gate and the WT arm in the same cell with any
reproducibility.** The new rigidified designs fail the WT arm outright, and several *gain*
position-6 contacts on WT. The parent bb54-g2 passes the WT criterion in 2 of 4 independent
3-draw calls at the very cell where its gate is 3-of-3: at n = 3 the verdict flips call to
call. Pooled over 9 vs 12 distinct draws it passes on both medians (-33.4 % / -52.4 %) but
with **no strict separation** on either quantity. Batch-06's much larger margins
(-59.2 % / -73.1 %, measured under `recyc8`) are not reproduced at n > 3; they were a
three-draw high excursion. This is the honest read: **bb54-g2's WT discrimination is real but
marginal, roughly a 30-50 % median effect with fully overlapping distributions.**

### Apo and decoy (bb54-g2, r3np, the only design reaching the WT arm)

Decoy **AVFDRKSDAK** on the same HLA: peptide buried area 18.5 / 84.5 / 111.7 A^2
(median 84.5) vs 208.3 A^2 on G12D - **-59 %**; position-6 contacts 0-7; no Asp6 gate pass in
any draw. The decoy arm holds.

Apo (peptide-free HLA + b2m) - run **unconstrained**, as `CONTROL_7OW6.md` did, because there
is no chain C to constrain against, and therefore **not like-for-like**: HLA buried area
**1102 / 1212 / 1239 A^2**, against 379-758 A^2 (median 613) of HLA buried in the G12D
complexes. **Result 5 - the apo failure persists, unchanged.** The binder buries roughly
twice as much empty groove as it does peptide-loaded groove. Rigidification did nothing to
it, which is expected: the problem is interface *shape*, not flexibility, and this batch
changed no interface shape (the donor shell was protected by design).

### Peptide share and the unapplied 0.40 floor

Across all 183 constrained G12D draws in this batch: `peptide_buried_frac` median **0.229**,
max **0.396**. **0 of 183 reach the 0.40 hard floor.** 137 of 183 clear the written-up but
**deliberately unapplied** 0.18 recalibration. For reference the floor rejects the real 7OW6
TCR itself (0.224 predicted, 0.207 crystal). Nothing was changed in `src/design_gates.py`;
these are values only, and whether 0.40 or 0.18 applies remains the user's call.

Buried area, bb54-g2 at r3np (medians): peptide 208.3 A^2 vs HLA 613 A^2 - the HLA groove is
still ~75 % of the interface, unchanged from batch-06's ~79 %. Consistent with batch-06's
finding that groove contact has to be cut by **shape**, not by size or stiffness.

## Pose spreads

Binder CA RMSD in the pMHC frame (HLA chain A + peptide chain C), binder **not** re-fitted,
frame RMSD reported as the artefact floor; in the apo arm Boltz-2 relabels the binder to
chain C and there is no peptide, so the frame is HLA chain A alone (stated, not substituted).
Frame floors are 0.21-0.55 A throughout, so the spreads below are real, not artefact.

| design / arm / cell | mean | max | frame floor |
|---|---|---|---|
| bb54-g2 g12d r10np (tightest in sweep, gate 2/3) | 2.25 | 2.81 | 0.28 |
| bb54-g2 g12d r8np | 3.87 | 5.05 | 0.44 |
| bb54-g2 g12d r3np (the 3-of-3 cell) | 13.93 | 16.97 | 0.39 |
| bb54-g2 wt r3np | 17.98 | 22.78 | 0.23 |
| bb54-g2 apo r3np (HLA-only frame) | 9.53 | 13.42 | 0.36 |
| bb54-g2 decoy r3np | 27.42 | 31.89 | 0.25 |
| bb54_rigid-gp-p2 g12d r3 (its 3-of-3 cell) | 13.74 | 16.26 | 0.26 |

## Artefacts and honesty notes

- One pair of independent calls with byte-identical request bodies (`r3np_rep` and
  `r3np_rep3`) returned **byte-identical structures** - the known Boltz-2 artefact. All
  pooled statistics above exclude the duplicate call; 213 of 216 structures are distinct, and
  every 3-draw call was internally 3-distinct.
- Wave-1 (`diffusion_steps 30`) artifacts are preserved as `*_w1_ds30.json` rather than
  deleted.
- No staple was ever formable on bb54; the relaxed 6.5 A CB-CB window in `scripts/rigidify.py`
  is documented in-file as a hypothesis, and it produced nothing, so it changed no result.

## Verdict

Three straight batches stalled at 2-of-3; batch-07 shows that was a **sampling-settings**
artefact, not a design property - bb54-g2 is 8/9 on the Asp6 read-out under `r3np`. But the
two things that would make it submittable do not come together: the WT margin is marginal and
flips at n = 3, there is no strict separation, the apo failure is untouched, and
`peptide_buried_frac` never reaches 0.40 (best 0.396 of any draw in the batch). Rigidification
is a measured dead end, like shrinking was: 1 of 21 designs hit 3-of-3 once and did not
repeat, and self-rigidifying the parent degraded it.

Design slots used: **0 of 40**.

## Files

- `panel_checkpoint.jsonl` - every scored complex, one line per draw, fsync'd as it was written.
- `screen_summary.jsonl`, `screen_stats.json` - the 21 candidates and their per-cell gate counts.
- `analysis.json` - per-cell medians, distinct-structure counts and pose spreads.
- `cross_condition.json`, `depth_results.json`, `pooled_wt.json` - the WT-reproducibility probes.
- `structures/` - the 18 curated complexes (bb54-g2's four arms at `r3np`, `bb54_rigid-gp-p2`'s
  G12D and WT at `r3`). The other 198 are in the untracked `all_structures/`, exactly as the
  repo's bulk-dump convention requires; `panel_checkpoint.jsonl` still names them under
  `structures/`, which is where they were written before curation.
