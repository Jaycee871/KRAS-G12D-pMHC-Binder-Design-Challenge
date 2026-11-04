# Replication confirmation and the first off-target selectivity panel

**Date:** 2026-11-03 (Asia/Taipei).
**Scripts:** `scripts/selectivity_panel.py` (step 1 calls
`scripts/run_batch04.score_candidate` verbatim, with the sampling fields
patched exactly as `scripts/sampling_study.patched_boltz2` does; step 2 adds a
generic per-arm measurement built from the same `src/complex_analysis.py`
primitives) and `scripts/selectivity_panel_analyze.py` (pose spread re-uses
`scripts/replicate_rmsd_batch04.binder_rmsd`).
**Raw data:** `checkpoint.jsonl` (57 rows, authoritative), `structures/`
(57 mmCIF), `analysis.json`.
**Boltz-2 calls:** **57** — 21 in step 1, 36 in step 2. Strictly serial,
≥ 8 s interval, one process, nothing else running.
**Track-1 design slots consumed: 0 of 40.** No design was generated: all three
sequences already existed in batch-04/batch-05. Re-scoring costs API calls
only. Nothing submitted, nothing pushed to `main`, no PR merged.
**No gate or threshold was changed.** `src/design_gates.py` is untouched.

## The two questions

1. The sampling study (PR #10) collapsed the replicate pose spread of three
   designs under design-specific Boltz-2 settings — but at n = 3 draws per
   cell, which is too thin to call them survivors. **Does the collapse hold at
   n = 10?**
2. No design in this project has ever been shown to be **selective**, only to
   read the aspartate. **Does the Asp6 read-out survive the off-target panel?**

Selectivity is judged here on **measured geometry** and on **reproducibility
across draws**. The ipTM `specificity_gap` and the paired same-pose MM-GBSA
ddG both failed the 7OW6 positive control (`results/validation/CONTROL_7OW6.md`,
`results/ddg-scoring/FINDINGS.md`) and are retired: neither is computed, ranked
on, revived or tuned anywhere in this run. ipTM appears below as a
**pose-quality read-out only** — never a gap, never a tie-breaker.

## Method

* **Step 1 — G12D, n = 10 per design.** Each design under its own stabilising
  setting: `bb50-g1` and `bb25` at `recycling_steps` 8 (`recyc8`), `bb135-g1`
  at `without_potentials: true` + `recycling_steps` 8 (`combo`). The three
  draws already on disk from the sampling study are re-used, not re-run, so
  only 7 new draws per design were paid for. Same polymers, same
  `diffusion_samples: 1`, same `step_scale: 1.638`, same validated flat pocket
  constraint `{"binder":"D","contacts":[{"id":"C","residue_index":6}],"max_distance":6.0}`.
* **Step 2 — the full four-state panel**, 4 draws per off-target arm, under the
  same per-design setting:

  | arm | polymers | constraint |
  | --- | --- | --- |
  | **g12d** `VVVGADGVGK` | HLA-A\*11:01 + β2m + peptide | pocket term on C6 (step 1, n = 10) |
  | **wt** `VVVGAGGVGK` | same | pocket term on C6, identical |
  | **decoy** `AVFDRKSDAK` | same | pocket term on C6, identical |
  | **apo** | HLA-A\*11:01 + β2m only | **dropped — run unconstrained** |

  **What was done for the apo arm, and why.** The constraint names peptide
  chain `C` residue 6; in the apo state there is no chain C, so the term cannot
  be expressed and the call is made **unconstrained**. This is exactly what
  `scripts/control_7ow6_constrained.py` did for this arm in the 7OW6 control.
  The consequence is stated plainly in the verdict: the apo arm is **not**
  a like-for-like comparison with the three constrained arms, because nothing
  is pulling the binder onto the groove face in it.
* **Off-target measurement.** `complex_analysis.interface_metrics` cannot run
  unchanged on the off-target arms — WT position 6 is a glycine with no side
  chain, and apo has no peptide chain, both of which it (correctly) refuses.
  The arms are therefore measured with the same primitives that function uses,
  generalised over the peptide: position-6 **side-chain** contacts for any
  residue, closest gate-qualifying binder donor (Arg/Lys/His/Ser/Thr/Tyr) to
  the position-6 side chain, and donor approach to **every** acidic peptide
  side chain, so a binder that has merely learned to like a carboxylate shows
  up on the decoy's Asp4/Asp8. The generic code was validated against a
  sampling-study G12D draw and reproduces `interface_metrics` exactly
  (7 side-chain contacts, donor 2.77 Å, `peptide_buried_frac` 0.216,
  charge −2.0 — identical to the published row).
* **Pose.** Binder (chain D) C-alpha RMSD after superposing on the **pMHC
  frame only** (HLA chain A + peptide chain C), binder **not** re-fitted; the
  frame RMSD is reported alongside as the fitting-artefact floor. For apo the
  frame is HLA chain A alone (there is no peptide), stated rather than
  silently substituted. Boltz-2 relabels chains sequentially, so in the apo
  arm the binder comes back as chain **C**, not D; it is identified by length.

### One measurement caveat that affects how n should be read

Repeat calls with a byte-identical request body sometimes return a
**byte-identical structure**. Over the 10 G12D draws the number of *distinct*
structures is 6 (`bb50-g1`), 7 (`bb25`) and 8 (`bb135-g1`); on the
**unconstrained apo** arm all 4 draws are identical for two of three designs.
So the effective number of independent samples is lower than the nominal n —
somewhere between 4 and 8 per design rather than 10. This does not invalidate
anything below (every distinct draw is a genuine independent sample, and the
duplicates are consistent with the control's observation that unconstrained
calls reproduce to four decimal places), but **n = 10 should be read as
"6–8 independent draws plus repeats"**, and it is why no claim here rests on a
single draw.

## Step 1 — does the collapse hold at n = 10?

### `bb50-g1` under `recyc8` — **yes, 10 / 10**

| draw | gate | donor (Å) | Asp6 **side-chain** contacts | pos6 contacts (all) | `peptide_buried_frac` | `asp6_buried_frac` | net charge | ipTM\* |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 0 † | pass | 2.77 | 7 | 13 | 0.216 | 0.724 | −2.0 | 0.468 |
| 1 † | pass | 2.73 | 7 | 16 | 0.224 | 0.706 | **−3.0** | 0.532 |
| 2 † | pass | 2.77 | 7 | 13 | 0.216 | 0.724 | −2.0 | 0.468 |
| 3 | pass | 2.71 | 7 | 17 | 0.221 | 0.764 | −2.0 | 0.467 |
| 4 | pass | 2.79 | 7 | 14 | 0.218 | 0.708 | **−3.0** | 0.324 |
| 5 | pass | 2.71 | 7 | 17 | 0.221 | 0.764 | −2.0 | 0.467 |
| 6 | pass | 2.71 | 7 | 17 | 0.221 | 0.764 | −2.0 | 0.467 |
| 7 | pass | 2.53 | 5 | 13 | 0.218 | 0.863 | −2.0 | 0.447 |
| 8 | pass | 2.53 | 5 | 13 | 0.218 | 0.863 | −2.0 | 0.447 |
| 9 | pass | 3.11 | 4 | 9 | 0.239 | 0.619 | −1.0 | 0.375 |

† the three sampling-study draws, re-used. \* pose-quality read-out only.

**Binder C-alpha RMSD, all 45 pairs: mean 1.58 Å, max 3.78 Å, frame ≤ 0.12 Å.**
The n = 3 figure was 0.70 / 1.05 Å. The collapse **holds** — there is one
docked mode, not two — but it is **looser than n = 3 advertised**: draw 9
(donor 3.11 Å, 4 side-chain contacts) sits ~3.8 Å off the main cluster. Every
draw still reads the carboxylate. Donor distance 2.53–3.11 Å, side-chain
contacts 4–7 (7 in six of ten draws).

### `bb25` under `recyc8` — **9 / 10, and the second mode is still there**

| draw | gate | donor (Å) | Asp6 side-chain contacts | pos6 contacts | `peptide_buried_frac` | `asp6_buried_frac` | net charge | ipTM\* |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 0 † | pass | 2.12 | 15 | 19 | 0.284 | 0.852 | +1.0 | 0.367 |
| 1 † | pass | 2.52 | 12 | 16 | 0.253 | 0.876 | +2.0 | 0.459 |
| 2 † | pass | 2.26 | 13 | 21 | 0.307 | 0.837 | +1.0 | 0.399 |
| 3 | pass | 2.26 | 13 | 21 | 0.307 | 0.837 | +1.0 | 0.399 |
| 4 | pass | 1.94 | 8 | 14 | 0.305 | 0.703 | +1.0 | 0.378 |
| 5 | pass | 1.94 | 8 | 14 | 0.305 | 0.703 | +1.0 | 0.378 |
| 6 | pass | 2.26 | 13 | 21 | 0.307 | 0.837 | +1.0 | 0.399 |
| 7 | pass | 2.88 | 16 | 21 | 0.281 | 0.788 | +1.0 | 0.437 |
| **8** | **FAIL** | — | **0** | **1** | 0.206 | 0.086 | −1.0 | 0.456 |
| 9 | pass | 2.67 | 16 | 21 | 0.307 | 0.728 | +1.0 | 0.398 |

**Binder C-alpha RMSD: mean 3.69 Å, max 13.94 Å, frame ≤ 0.49 Å.** The n = 3
figure was 1.96 / 2.54 Å. **The collapse does not fully hold.** Nine draws sit
in one tight mode with the shortest donor contacts in the project
(1.94–2.88 Å, 8–16 side-chain contacts); one draw in ten flips ~14 Å away and
loses the aspartate entirely. `recyc8` suppressed the flip rate — batch-04 saw
it 1 draw in 4, batch-05's baseline 7.27 Å mean — but did not eliminate it.
At n = 3 we were simply lucky.

### `bb135-g1` under `combo` — **yes, 10 / 10**

| draw | gate | donor (Å) | Asp6 side-chain contacts | pos6 contacts | `peptide_buried_frac` | `asp6_buried_frac` | net charge | ipTM\* |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 0 † | pass | 2.81 | 5 | 11 | 0.347 | 0.339 | 0.0 | 0.384 |
| 1 † | pass | 2.42 | 5 | 11 | 0.325 | 0.274 | 0.0 | 0.562 |
| 2 † | pass | 2.59 | 5 | 10 | 0.331 | 0.516 | 0.0 | 0.463 |
| 3 | pass | 2.69 | 5 | 10 | 0.303 | 0.326 | 0.0 | 0.450 |
| 4 | pass | 2.59 | 5 | 10 | 0.331 | 0.516 | 0.0 | 0.463 |
| 5 | pass | 2.61 | 5 | 12 | 0.314 | 0.396 | 0.0 | 0.626 |
| 6 | pass | 2.88 | 4 | 9 | 0.282 | 0.256 | 0.0 | 0.426 |
| 7 | pass | 2.61 | 5 | 12 | 0.314 | 0.396 | 0.0 | 0.626 |
| 8 | pass | 2.71 | 5 | 12 | 0.318 | 0.440 | 0.0 | 0.449 |
| 9 | pass | 2.67 | 5 | 12 | 0.316 | 0.458 | 0.0 | 0.443 |

**Binder C-alpha RMSD: mean 2.17 Å, max 4.79 Å, frame ≤ 0.25 Å** (n = 3 gave
2.04 / 2.97 Å). The collapse **holds**, and this is the most monotonous design
in the set: 4–5 side-chain contacts and net interface charge **0.0 in every
single draw**, donor pinned at 2.42–2.88 Å.

### Step 1 verdict

| design | setting | gate at n = 10 | donor (Å), range | Asp6 side-chain contacts | RMSD mean / max (frame) | collapse holds? |
| --- | --- | --- | --- | --- | --- | --- |
| `bb50-g1` | recyc8 | **10 / 10** | 2.53–3.11 | 4–7 | 1.58 / 3.78 Å (≤ 0.12) | **yes**, looser than n = 3 (0.70 / 1.05) |
| `bb135-g1` | combo | **10 / 10** | 2.42–2.88 | 4–5 | 2.17 / 4.79 Å (≤ 0.25) | **yes**, essentially unchanged |
| `bb25` | recyc8 | **9 / 10** | 1.94–2.88 | 0–16 | 3.69 / 13.94 Å (≤ 0.49) | **partly** — one 14 Å mode flip in ten |

Stated plainly: **the collapse is real and it holds at higher n for two of the
three designs; for `bb25` the two-mode behaviour survives at roughly a
1-in-10 rate and the n = 3 number understated the spread.** All three were
carried into step 2 — `bb25` because 9/10 is still a defensible survivor, with
the caveat recorded.

### The unchanged hard-gate set, for the record

`src/design_gates.evaluate()` with thresholds untouched, over the 10 G12D
draws each:

* `peptide_buried_frac ≥ 0.40` — **failed by all 30 draws**. Best single draw
  in this run is `bb135-g1` at 0.347. The floor also rejects the real 7OW6 TCR
  (0.224 predicted / 0.207 crystal); the 0.18 recalibration written up in
  `CONTROL_7OW6.md` remains **unapplied and is the user's call**.
* net interface charge ≥ −2.0 — `bb50-g1` breaches it on **2 of 10** draws
  (−3.0); `bb25` and `bb135-g1` never do.
* `bb25` draw 8 additionally fails the position-6 contact and Asp6 read-out
  gates, as the table shows.

Nothing above was adjusted to produce survivors. These are reported values.

## Step 2 — the off-target panel

Four draws per arm, per design. **Read the WT column carefully**: peptide
position 6 in the WT state is a *glycine*, which has no side chain, so "Asp6
side-chain contacts = 0" there is **true by construction and carries no
information**. The informative WT quantities are whether the binder still
engages position 6 at all (backbone contacts), how much peptide it still
buries, and where its pose goes.

### `bb50-g1` (recyc8)

| arm | pos6 contacts (all heavy atoms) | pos6 **side-chain** contacts | donor to pos6 side chain | `peptide_buried_frac` | buried area total / on peptide | net charge | ipTM\* | binder RMSD vs G12D draw 0 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| **g12d** (n = 10) | 9–17 | **4–7** | **2.53–3.11 Å, 10/10** | 0.216–0.239 | 1076 / 238 Å² | −1 … −3 | 0.32–0.53 | — (internal spread 1.58 Å) |
| **wt** (n = 4) | 12–19 | 0 (Gly) | none | 0.191–0.242 | 895 / 185 Å² | 0.0 … −2.0 | 0.36–0.56 | **15.1–17.0 Å** |
| **apo** (n = 4) | — | — | — | — | **1454 / —** | −2.0 | 0.36 | 12.3 Å |
| **decoy** (n = 4) | **0–2** | 0–2 | none ≤ 3.5 Å | **0.108–0.132** | **501 / 62 Å²** | 0 … −2 | 0.26–0.41 | 22.1–26.4 Å |

* **Decoy: rejected, clearly.** The peptide footprint collapses from 238 Å² to
  62 Å² and position-6 contacts from 9–17 to 0–2 in every draw. One draw does
  make an **Arg44:NH2···Asp8 salt bridge at 2.95 Å** on the decoy's *other*
  aspartate — carboxylate promiscuity is visible, but it does not buy a
  peptide-centred pose.
* **WT: not rejected.** The binder still makes 12–19 heavy-atom contacts with
  the WT glycine at position 6 and still buries 185 Å² of the peptide against
  238 Å² on G12D. It moves to a different pose (15–17 Å) but does not let go.
* **Apo: binds harder, not less.** 1454 Å² buried on empty HLA, *more* than the
  1076 Å² it buries on the G12D complex. Caveat: this arm is unconstrained, so
  the binder is free to find its best HLA site; it is not like-for-like.

### `bb25` (recyc8)

| arm | pos6 contacts | pos6 side-chain contacts | donor to pos6 side chain | `peptide_buried_frac` | buried area total / peptide | net charge | ipTM\* | RMSD vs G12D draw 0 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| **g12d** (n = 10) | 1–21 | **0–16** (8–16 on 9/10) | **1.94–2.88 Å, 9/10** | 0.206–0.307 | 667 / 190 Å² | +1 … −1 | 0.37–0.46 | — (internal spread 3.69 Å) |
| **wt** (n = 4) | **14–17** | 0 (Gly) | none | **0.280–0.290** | 685 / 196 Å² | −1 … −2 | 0.38–0.45 | 10.6–11.4 Å |
| **apo** (n = 4) | — | — | — | — | **1254 / —** | −2.0 | 0.44 | 25.3 Å |
| **decoy** (n = 4) | 0–23 | **0–23** | 3.27 Å (1 draw, to Lys6 **CB**) | 0.065–0.391 | 431 / 118 Å² | 0 … −2 | 0.26–0.52 | 18.7–31.6 Å |

* **WT: not rejected, and the footprint is if anything larger.** 196 Å² of WT
  peptide buried against 190 Å² on G12D, 14–17 position-6 contacts in all four
  draws, `peptide_buried_frac` 0.287 vs 0.29. There is no measured loss of
  engagement when the aspartate is removed.
* **Decoy: ambiguous.** Two draws bury 129–208 Å² of the decoy peptide with
  10–23 contacts at position 6 (a lysine); two collapse. No carboxylate bond to
  Asp4/Asp8 was formed in any draw.
* **Apo: 1254 Å² on empty HLA**, nearly double its G12D total. Same
  unconstrained caveat.

### `bb135-g1` (combo)

| arm | pos6 contacts | pos6 side-chain contacts | donor to pos6 side chain | `peptide_buried_frac` | buried area total / peptide | net charge | ipTM\* | RMSD vs G12D draw 0 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| **g12d** (n = 10) | 9–12 | **4–5** | **2.42–2.88 Å, 10/10** | 0.282–0.347 | 627 / 200 Å² | **0.0** | 0.38–0.63 | — (internal spread 2.17 Å) |
| **wt** (n = 4) | **13–24** | 0 (Gly) | none | **0.303–0.328** | 603 / 192 Å² | 0 … −1 | **0.52–0.72** | **4.8–8.0 Å** |
| **apo** (n = 4) | — | — | — | — | **1101 / —** | −3.0 | 0.46 | 15.2 Å |
| **decoy** (n = 4) | **18–28** | **18–28** | 4.10–5.17 Å (Ser11 OG → Lys6) | 0.225–0.284 | 757 / 189 Å² | 0 … −2 | 0.30–0.39 | 17.0–17.4 Å |

* **WT: the clearest negative result in the panel.** The binder keeps
  essentially the same pose (4.8–8.0 Å, against an internal G12D spread of
  2.17 Å), buries 192 Å² of WT peptide against 200 Å² on G12D, makes *more*
  contacts with position 6 (13–24 vs 9–12), and its ipTM is **higher on WT than
  on G12D in every draw** (0.52–0.72 vs 0.38–0.63). It loses the Asp6 hydrogen
  bond because the acceptor no longer exists, and nothing else changes.
* **Decoy: it reads the decoy's aspartate.** Arg14:NH2 forms a salt bridge to
  the decoy **Asp8** carboxylate in 3 of 4 draws at **2.70, 3.26 and 3.89 Å**,
  while burying 189 Å² of the decoy peptide — the same order as its on-target
  footprint. This is exactly the failure mode `docs/PIPELINE.md` section 1
  warns the decoy arm exists to catch: a paratope that has learned to like *a*
  carboxylate on the peptide face, not *this* carboxylate at *this* position.
* **Apo: 1101 Å² on empty HLA**, against 627 Å² total on the G12D complex.

## Verdict — is anything selective?

**No. On geometric evidence, none of the three designs is selective for G12D,
and one of them is positively anti-selective.** This is a clean negative
result and nothing was tuned to avoid it.

| design | G12D read-out | WT | apo | decoy | selective? |
| --- | --- | --- | --- | --- | --- |
| `bb50-g1` | **10/10**, donor 2.53–3.11 Å, 4–7 side-chain contacts | keeps 185 Å² of peptide, 12–19 pos6 contacts | buries **1454 Å²** on empty HLA (> its G12D total) | **collapses** (62 Å², 0–2 contacts) | **no** — best of three, decoy-clean, but WT- and apo-blind |
| `bb25` | 9/10, donor 1.94–2.88 Å, 8–16 contacts | **196 Å² vs 190 Å² on G12D** — no loss at all | 1254 Å² on empty HLA | ambiguous (2 of 4 draws engage) | **no** |
| `bb135-g1` | 10/10, donor 2.42–2.88 Å, charge 0.0 every draw | **same pose, same burial, more contacts, higher ipTM than G12D** | 1101 Å² on empty HLA | **salt-bridges the decoy's Asp8 in 3 of 4 draws** | **no — anti-selective** |

Three things are now measured that were not before:

1. **The Asp6 read-out is reproducible.** Two designs make the contact in 10
   draws out of 10, one in 9 of 10, with donor distances of 1.9–3.1 Å. The
   chemistry result from batches 04/05 and the pose result from the sampling
   study both survive at higher n (with the `bb25` caveat above).
2. **Reading the aspartate is not the same as needing it.** Removing the
   aspartate (WT) costs these binders the hydrogen bond and **nothing else** —
   not peptide burial, not position-6 engagement, not ipTM. The interface is
   carried by the peptide backbone and the groove rim, which are identical in
   both states. That is the single most important number in this document.
3. **Peptide-dependence is not demonstrated either.** All three bury more
   surface on **peptide-free** HLA than on the G12D complex. The apo arm is
   unconstrained and so is not a like-for-like comparison — it is suggestive,
   not conclusive — but it points the wrong way, and the real 7OW6 TCR's one
   genuine ipTM signal in the control was precisely peptide-present vs absent
   (0.93 vs 0.78), which none of these designs reproduces.

**What would have counted as selectivity and did not happen:** the Asp6
side-chain contact persisting on G12D while the *peptide footprint itself*
shrank on WT, apo and decoy. Only the decoy arm of `bb50-g1` shows that
pattern, and one arm out of three is not selectivity.

**Honest limits of this panel.** 4 draws per off-target arm, with duplicate
structures reducing the effective count (2–4 distinct per arm). The apo arm is
unconstrained, as it must be. ipTM is a pose-quality read-out and the ipTM
differences quoted above are *descriptions of the predictions*, not a
specificity gap — no gap was computed. And a Boltz-2 co-folding prediction is
not an experiment: nothing here is validated against binding data.

## Recommendation

1. **There is no submittable set today.** Not one design, not three. A
   submission needs evidence of mutation specificity, and this panel is the
   first direct measurement of it in the project: the answer is negative. Do
   **not** submit when the window opens at 07:30 Taipei on Nov 4. The deadline
   is Nov 8 (AoE) and **0 of 40 Track-1 design slots are spent**, so there is
   no scarcity argument for spending one on a design now measured to bind the
   wild-type peptide as well as the mutant. **The submission decision is the
   user's alone.**
2. **The design objective has to change, and this panel says how.** Every one
   of these binders is dominated by HLA and peptide-backbone contact
   (`peptide_buried_frac` 0.21–0.35, 600–1100 Å² total, most of it on the
   groove). A design that needs the aspartate has to pay for the interface
   *with* the aspartate — i.e. a much smaller, peptide-centred paratope whose
   total buried area is low and whose Asp6 contribution is a large share of it,
   rather than a 1000 Å² interface with one good salt bridge in it. Screening
   future batches on `buried_area_peptide / buried_area_total` **plus the WT
   arm**, from the start, is the concrete change.
3. **Put the WT arm inside the gate, not after it.** Batch-05 folded
   replication into the gate and it worked. The same move is now needed for
   selectivity: a design should have to lose peptide footprint on WT before it
   is ranked at all. At 4 draws × 2 arms that is 8 calls per survivor — cheap
   compared with what was spent getting here.
4. **`bb50-g1` is the one worth keeping as a reference**: it is the only design
   whose footprint collapses on an off-target peptide, and its pose is the
   tightest in the project (1.58 Å mean over 10 draws). Its problems are the
   apo arm and a net interface charge that breaches −2.0 on 2 of 10 draws.
5. **Settle the `peptide_buried_frac` floor** (0.40 vs the control-calibrated
   0.18) — still unapplied, still the user's call, and still the binding
   constraint on every design in this project.

**Nothing has been submitted, no PR merged, nothing pushed to `main`.
Track-1 stands at 0 of 40 design slots used.**
