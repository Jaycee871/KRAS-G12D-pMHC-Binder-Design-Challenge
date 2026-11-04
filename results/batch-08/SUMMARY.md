# Batch-08 - steering the interface by SHAPE, to cut HLA groove contact

Branch `results/batch-08-shape-steered`, stacked on `results/batch-07-rigidified-graft`.
Target: KRAS G12D decamer **VVVGADGVGK** on HLA-A\*11:01 (7OW6). WT comparator
**VVVGAGGVGK**, decoy **AVFDRKSDAK**, apo = peptide-free HLA + b2m.

Nothing in `src/design_gates.py` was touched. No gate or threshold was loosened,
retuned or re-derived anywhere in this batch. The WT-loss criterion is
`scripts/batch06_panel.WT_CRITERION` re-used **verbatim**, fixed before batch-06's
results and applied here unchanged. No ipTM `specificity_gap` and no paired MM-GBSA
`ddG` was computed anywhere - both stay retired; ipTM appears only as a pose-quality
read-out. **43 Boltz-2 calls** (129 structures, every one distinct), strictly serial at
>= 8 s spacing, plus 192 RFdiffusion and 24 ProteinMPNN calls.
**0 of the 40 Track-1 design slots used** (nothing submitted).

## Headline

**Null result, and a costly one for the incumbent.** Shape steering did not re-centre
the interface: the peptide share of buried area over 114 on-target draws has median
**0.195**, against batch-06's 0.213 and bb54-g2's ~0.25-0.29. The HLA alpha-1/alpha-2
helices still take a median **0.744** of the buried area. Not one of the 36 new designs
reached the Asp6 gate 3-of-3. The de-wetting twin test, run paired on bb54-g2 itself,
*destroyed* the gate (3/3 -> 0/3) without moving the helix share. And pooling the three
new bb54-g2 draws into the r3np record drops its WT margin from -33.4 % to **-29.4 %**,
which is **below** the 30 % criterion - bb54-g2 no longer clears the WT arm pooled.

**The apo failure is not broken.** On the two most peptide-centred batch-08 designs the
binder buries **1110-1287 A^2** on peptide-free HLA against **542-758 A^2** on the G12D
complex - apo burial is still ~1.8x on-target burial, the same sign and a worse ratio
than bb54-g2's (apo median 1211.5 vs on-target 842.1 A^2). Groove dependence is intact.

---

## Step 1 - shape-steered backbone generation (192 RFdiffusion calls)

Four strategies x 48 backbones, `diffusion_steps` 15 (30 is the value that gave 0 of 96),
seed 8080, run by `scripts/run_batch08.py`. Binder length held at bb54-g2's 46 aa or a
little above - never smaller, since shrinking is a measured dead end - and never
rigidified.

| strategy | length | hotspots | rationale | base-filter OK |
|---|---|---|---|---|
| B8-apex-46 | 46 | C6 | bb54-g2's size and hotspot, kept only when the backbone comes down on position 6 from above | 4/48 |
| B8-span-50 | 50 | C4, C6, C8 | hotspots spread along the peptide long axis, so the paratope lies *along* the peptide rather than across the groove | 6/48 |
| B8-crest-46 | 46 | C5, C6, C7 | the solvent-exposed crest of the bulge only | 4/48 |
| B8-apex-54 | 54 | C6 | slightly larger than bb54-g2 at MPNN temp 0.25, for topological diversity | 2/48 |

`src/groove_geometry.py` is the new measurement module and contains no gate. It provides:

* `buried_area_partition` - the same Shrake-Rupley SASA difference
  `complex_analysis.interface_metrics` uses (1.4 A probe, binder present vs binder chain
  detached), but split **peptide / HLA alpha-1 (57-84) / HLA alpha-2 (138-180) / rest of
  the heavy chain / b2m**, so `buried_area_total` is directly comparable with every
  earlier batch while "cut groove contact" becomes measurable instead of asserted.
  Earlier batches only split peptide vs HLA chain, which cannot tell "lying along the
  helices" apart from "resting on the groove rim".
* `design_shape_metrics` / `contact_partition` - the same split on heavy-atom contact
  counts, computed on the **RFdiffusion backbone** (binder chain A, concatenated target
  chain B) before any Boltz-2 call is paid for.
* `groove_frame` / `approach_metrics` - a local orthonormal frame on the groove: peptide
  long axis (first PC of its CAs), `up` = groove normal from the beta-sheet floor
  centroid towards the peptide, `lateral` = up x axis, helix wall to helix wall. The
  binder's CA cloud is projected into it to give `approach_height` (mean height above
  the peptide centroid), `abs_lateral_offset` (displacement towards a wall) and
  `over_peptide_frac` (share of CAs within 7 A of the groove axis).
* `helix_facing_positions` - binder positions whose closest heavy-atom approach to
  alpha-1/alpha-2 is under 6 A and closer than their approach to the peptide. Coarse by
  design: it is used to *select against* helix contact in sequence design, not to score.

The five shape thresholds were fixed **before the batch ran**, calibrated on the 97
backbones batch-07 had already generated: `peptide_contact_frac >= 0.30`,
`helix_contact_frac <= 0.30`, `approach_height >= 2.0 A`, `abs_lateral_offset <= 5.0 A`,
`over_peptide_frac >= 0.55`, on top of batch-06's base pose filter.

**Result 1 - RFdiffusion essentially cannot be asked for this shape.** Of 192 backbones,
16 cleared the base pose filter and **0 cleared all five shape thresholds**. The binding
constraint is geometric, not chemical: median `approach_height` over all 192 backbones is
**-20.4 A** and over the 16 base-OK ones **-7.0 A**, i.e. the binder centroid sits
*below* the peptide, hanging off the side of the groove. `approach_height` alone killed
14 of the 16; helix contact killed 5, peptide contact 3, lateral offset 2. The ~7 % yield
the thresholds were calibrated for did not materialise at all. Median
`helix_contact_frac` among base-OK backbones was 0.222 and median
`peptide_contact_frac` 0.338 - the *contact* numbers were fine; the *approach* was not.

Everything scored below therefore came from the ranked **reserve list** (lowest
`helix_contact_frac` first), widened in a second wave from 8 to 16 backbones with
`--from-existing` (no new RFdiffusion calls): 22 then 37 candidates, 13 backbones.
**Every design in this batch is a reserve-list design.** That is recorded on each row.

## Step 2 - selecting against helix contact in sequence design (24 ProteinMPNN calls)

Two changes, both in selection only. The MPNN selection score became
`donor_bias_score - 0.5 * helix_apolar_penalty`, where the penalty sums distance-weighted
F/W/Y/L/I/M/V at helix-facing positions. And every donor graft got a paired **de-wetted
twin** (`-dw`): at helix-facing positions only, F/W/Y/L/I/M/V -> Ser, with the six
positions nearest Asp6 protected so the donor chemistry is untouched by construction.
Graft and twin sit on the **same backbone**, so de-wetting is measurable, not confounded.

## Step 3 - the gate, unchanged, on 36 designs (36 calls)

`scripts/batch08_panel.py` imports `scripts/batch06_panel.py` verbatim. Cell **r3np**
(`nopot`: recycling_steps 3, sampling_steps 50, without_potentials True) - the cell in
which bb54-g2 reaches the gate 3-of-3 reproducibly. One call per design/arm with
`diffusion_samples: 3`, so screening *is* the 3-draw gate count and no design advances on
one draw. Same validated flat pocket constraint
`{"binder":"D","contacts":[{"id":"C","residue_index":6}],"max_distance":6.0}`,
`step_scale` 1.638. **All 129 structures across the batch were byte-distinct**
(`distinct_structures_in_call` = 3 for all 43 calls); no duplicate-response artefact.

**Result 2 - 0 of 36 designs reached the Asp6 gate 3-of-3.** Seven reached it 1-of-3, 29
never. Ten gate passes in 114 on-target draws (8.8 %). Because the gate is the entry
condition, the WT, apo and decoy arms were never triggered by the normal order.

### Buried area partitioned, on-target (114 draws, medians)

| quantity | batch-08 | bb54-g2 @ r3np (b07) | batch-06 screen | 7OW6 crystal TCR |
|---|---|---|---|---|
| peptide share of buried area | **0.195** | 0.291 | 0.213 | 0.207 |
| helix (alpha1+alpha2) share | **0.744** | 0.704 | - | 0.793 |
| buried area total (A^2) | 871.7 | 842.1 | - | 813.1 |
| buried area peptide (A^2) | 181.8 | 180.0 | - | 168.6 |
| buried area alpha-1 (A^2) | 307.6 | 249.8 | - | 477.2 |
| buried area alpha-2 (A^2) | 302.9 | 372.2 | - | 167.3 |
| buried area rest of HLA (A^2) | 14.5 | 3.2 | - | 0.0 |

**Result 3 - shape steering did not re-centre the interface.** The peptide-share
distribution (n=113 with a defined share) runs 0.0 - 0.48 with median **0.195**, i.e.
*below* batch-06's 0.213 and well below bb54-g2's 0.29. Nothing in the distribution
approaches a peptide-dominated interface. By strategy the medians are B8-crest-46 0.223,
B8-apex-46 0.213, B8-span-50 0.157, B8-apex-54 0.148 - the "lie along the peptide axis"
strategy was the *worst*. Note the comparator that matters most: the real, validated,
G12D-selective 7OW6 TCR itself sits at peptide share **0.207** and helix share 0.793. A
groove-dominated interface is what a working binder on this target looks like; three
batches of trying to move that number have now moved it by less than its draw-to-draw
noise.

### The de-wetted twins, paired on the same backbone (13 pairs)

**Result 4 - de-wetting the helix face does nothing to the helix share.** Median change
graft -> twin: peptide share **+0.002**, helix share **+0.017** (means -0.003 and -0.023).
Removing the greasy side chains that point at the walls does not reduce how much wall the
binder buries - the contact is driven by backbone placement, not by side-chain chemistry.
Five twins picked up a single gate pass their parent lacked and one lost its interface
entirely; none reproduced.

## Step 4 - the controlled de-wetting test on bb54-g2 itself (3 calls)

The shape-steered backbones never reached the gate, so the hypothesis was still untested
on a design that works. `scripts/batch08_dewet.py` runs the controlled version: take
bb54-g2, change **nothing** except the chemistry it points at the helices. Helix-facing
positions are read off bb54-g2's own predicted r3np poses, the six positions nearest Asp6
are protected, and the large hydrophobics among the rest become serine - six edits,
M1S M2S L3S V30S V42S V44S. Fold, length and the whole donor shell untouched. Parent and
twin run side by side, same cell, same segment, paired.

| design | Asp6 gate | median peptide area | median peptide share | median helix share |
|---|---|---|---|---|
| B8-bb54-g2-parent | **3/3** | 186.1 | 0.286 | 0.714 |
| B8-bb54-g2-dewet | **0/3** | 204.3 | 0.269 | 0.721 |

**Result 5 - de-wetting destroys the donor contact and leaves the groove untouched.** The
twin lost the Asp6 gate in all three draws (`asp6_buried_frac` 0.0, no donor within
3.5 A) while its helix share moved by +0.007. Six serines outside the protected shell
were enough to stop the binder engaging Asp6 at all. This is a clean negative: the lever
does not work, and it is not a near miss.

## Step 5 - the WT arm, pooled, criterion verbatim (1 call)

The parent reached 3-of-3, so the WT arm ran. On its own 3+3 draws the batch-08 segment
gives peptide-area loss **0.252** and pos6-contact loss 0.421 - **FAIL**, because the area
leg is under 0.30. Pooling every r3np draw of that sequence across batch-07 and batch-08
(`scripts/batch08_pooled_wt.py`, `results/batch-08/pooled_wt.json`), **n = 12 G12D vs 15
WT**, Asp6 gate 11/12:

| pooled quantity | value |
|---|---|
| median peptide area, G12D / WT | 197.2 / 139.2 A^2 |
| peptide-area loss | **0.294** (need >= 0.30) |
| median pos6 contacts, G12D / WT | 20.0 / 10.0 |
| pos6-contact loss | 0.500 (need >= 0.30) |
| strict separation (area / pos6) | false / false |
| **wt_gate_pass** | **false** |

**Result 6 - bb54-g2's WT margin does not survive more draws.** It was -59.2 % at n=3
(batch-06), -33.4 % at n=9/12 (batch-07), and **-29.4 % at n=12/15** - a monotone decay
towards the criterion and now just under it. The pos6 leg holds comfortably at -50 %; it
is the buried-area leg that fails. Nothing was retuned to produce this: it is the
batch-06 criterion applied verbatim to more data. The honest reading is that bb54-g2's
WT discrimination is real on position-6 contacts and not established on buried area.

## Step 6 - the apo and decoy diagnostic (4 calls)

Because nothing cleared the gate, the normal order never reached apo, and batch-08 would
have carried no apo arm at all. `scripts/batch08_apo_probe.py` runs one read on the two
most peptide-centred designs, by a rule stated before the calls (highest median
on-target peptide share). It changes no gate, threshold or verdict and cannot promote
anything - both designs had already failed at 1-of-3. apo is **unconstrained** (there is
no chain C to point the pocket constraint at), exactly as `CONTROL_7OW6.md` did, so it is
**not like-for-like** with the constrained arms; the apo binder is relabelled to chain C
by Boltz-2, so framing uses HLA chain A only.

| design | arm | buried area total (A^2) | peptide share | helix share |
|---|---|---|---|---|
| B8-apex-46_bb46-g2-dw | g12d | 659.3 (595-758) | 0.320 | 0.680 |
| B8-apex-46_bb46-g2-dw | **apo** | **1224.0 (1189-1287)** | 0 | 0.889 |
| B8-apex-46_bb46-g2-dw | decoy | 244.0 (224-318) | 0.299 | 0.701 |
| B8-crest-46_bb46-g2-dw | g12d | 543.8 (542-703) | 0.299 | 0.701 |
| B8-crest-46_bb46-g2-dw | **apo** | **1179.3 (1110-1253)** | 0 | 0.874 |
| B8-crest-46_bb46-g2-dw | decoy | 400.8 (367-571) | 0.253 | 0.681 |

**Result 7 - the apo failure is NOT broken. Apo burial does not drop below on-target
burial; it is roughly double it.** Median apo total across the six draws is **1206.7 A^2**
against **627.3 A^2** on target for the same two designs - a ratio of 1.9, worse than
bb54-g2's 1.4 (1211.5 / 842.1). The two designs chosen were the *most* peptide-centred in
the batch, so this is the best case, not the worst. On apo, 87.9 % of the buried area is
alpha-1/alpha-2 helix. These binders are HLA-A\*11:01 groove binders that tolerate a
peptide, not peptide binders.

The decoy arm is the one encouraging number: burial collapses to 244-401 A^2, **-63 % and
-26 %** against on-target, so the designs are not indiscriminate. But a decoy drop with no
apo drop is exactly the signature of groove binding that is sterically disturbed by a
different peptide, not of peptide recognition.

## Pose spread, peptide_buried_frac, and the hard gates

Binder C-alpha RMSD is reported in the pMHC frame with the binder **not** re-fitted;
the frame RMSD is the artefact floor. Across all 38 designs on target: mean RMSD median
**19.25 A** (range 1.85 - 34.65), per-design max median 23.48 A, **frame floor median
0.30 A** (0.18 - 0.40). Only 3 of 38 designs held a mean under 3 A. Reserve-list
backbones are much less reproducible than bb50-g1's 1.58/3.78 A reference; the tightest
here were B8-crest-46_bb46-g2 at 1.85/2.58 (floor 0.19) and B8-apex-46_bb46-g2 at
2.98/3.77 (floor 0.24) - and both are 0/3 on the gate, consistent with the established
anti-correlation between pose tightness and gate reproducibility.

`src/design_gates.py` was **not** modified and the 0.18 recalibration stays **UNAPPLIED**.
Reported, not applied, over the 114 on-target draws:

* `peptide_buried_frac >= 0.40`: **1 of 114 draws**, max **0.48**.
* `peptide_buried_frac >= 0.18` (the unapplied recalibration): 60 of 114.
* The single draw at 0.48 is `B8-crest-46_bb08_donor-anchored-g1` draw 0, and it is the
  only draw in the batch with **zero** `design_gates.evaluate` failures. It is also a
  near-collapsed interface - total buried area 386.2 A^2 - that clears the Asp6 route by
  contact count (10 pos6 contacts) with **no donor within 3.5 A**, so it fails the
  project's Asp6 direct read-out gate, and its two sibling draws are 0.091 and 0.068. It
  is noise, not a survivor.
* For reference the real 7OW6 TCR scores 0.224 predicted / 0.207 crystal against this
  same 0.40 floor, which is why the floor is not being applied.

## Cost

| stage | calls |
|---|---|
| RFdiffusion | 192 |
| ProteinMPNN | 24 |
| Boltz-2 panel, 36 designs x g12d | 36 |
| Boltz-2 de-wetting test (2 designs g12d + 1 WT) | 3 |
| Boltz-2 apo/decoy diagnostic | 4 |
| **Boltz-2 total** | **43** (129 structures, all distinct) |

38 designs scored. **0 of the 40 Track-1 design slots used.** Nothing submitted.

## What this batch establishes

1. **Shape steering fails the same way size and rigidity did.** Peptide share median
   0.195 vs 0.213 (batch-06) and ~0.29 (bb54-g2). Three independent levers - smaller,
   stiffer, re-shaped - have now each moved the peptide share by less than the
   draw-to-draw noise. The groove share is not a property of the binder; it is a property
   of the target surface.
2. **RFdiffusion will not produce a top-down paratope on this groove.** 0 of 192, with
   median approach height -20.4 A. The shape thresholds are not the problem; the
   generator's prior is. Any further attempt at this needs a different generator or
   explicit backbone placement, not more sampling.
3. **De-wetting the helix face is a dead end, cleanly measured both ways.** Nil effect on
   helix share across 13 paired twins; total loss of the Asp6 contact when applied to the
   one design that works.
4. **The apo failure stands, and it is the central result of this project.** Every design
   measured, in every batch, buries about twice as much on peptide-free HLA as on the
   G12D complex.
5. **bb54-g2's WT margin has decayed below the criterion at n=12/15.** Its claim to be
   G12D-selective now rests on position-6 contacts alone.

**Recommendation: do not submit.** There is no design in this project that clears the
gate as written - Asp6 3-of-3 *and* a pooled WT loss meeting the verbatim criterion - and
the one design that came closest has just fallen below the WT criterion as more draws
accumulated. A null result honestly reported is worth more than a submission of designs
that are demonstrably HLA-A\*11:01 groove binders rather than G12D neoantigen binders.
