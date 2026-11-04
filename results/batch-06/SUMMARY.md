# Batch-06 — a peptide-centred paratope, with the WT arm inside the gate

**Date:** 2026-11-03 → 2026-11-04 (Asia/Taipei).
**Branch:** `results/batch-06-peptide-centred`. PR against `results/selectivity-panel`.
**Scripts:** `scripts/run_batch06.py` (stages 1–3), `scripts/batch06_panel.py`
(stages 4–6: replication, WT arm, apo/decoy), `scripts/batch06_wt_diagnostic.py`
(a clearly-labelled **diagnostic**, see below), `scripts/batch06_analyze.py`.
**Raw data:** `checkpoint.jsonl` (41 rows, screening, authoritative),
`panel_checkpoint.jsonl` (30 rows), `analysis.json`, `panel_results.json`,
`wt_diagnostic.json`, `structures/` (30 mmCIF).
**Boltz-2 calls: 51** — 41 screening (1 per candidate) + 6 panel + 2 WT
diagnostic + 2 apo/decoy diagnostic. Strictly serial, ≥ 8 s interval, one
process at a time, never two Boltz-2 processes.
**Track-1 design slots consumed: 0 of 40.** Nothing has been submitted.
**No gate or threshold was changed.** `src/design_gates.py` is untouched.

> This run was **resumed** after the previous executor was interrupted by the
> platform. Stages 1–3 (both waves) were already complete on disk and were
> **not** recomputed: all 41 screened candidates were reconciled against
> `candidates.json` + `candidates_w2.json` and **0 remained to score**. Stage 4
> onwards is new work. The WT-loss criterion below was **already defined by the
> previous executor** in `scripts/batch06_panel.py` and is **kept verbatim** —
> it was not retuned after seeing any result.

## The question

The off-target panel (PR #11) found the project's central failure: **no design
is selective.** All three panel designs buried *more* area on apo peptide-free
HLA (1101–1454 Å²) than on the G12D complex (627–1076 Å²), and none let go of
the wild-type peptide. The diagnosis was that interfaces of 600–1100 Å²
dominated by HLA-groove contact make the aspartate a rounding error:
**reading the aspartate is not the same as needing it.**

Batch-06 tests the obvious remedy — make the paratope much smaller and
peptide-centred — and puts the WT arm **inside** the gate so that "needs it" is
what is actually measured.

## What was changed

1. **Much shorter binders:** 34 / 40 / 46 aa, against 62–72 aa in batches 04/05.
2. **A peptide-centred design-stage pose filter:** on top of batch-05's filter,
   a backbone must put ≥ 0.25 of its target contact on the peptide and may put
   no more than 90 heavy atoms on the HLA groove. An HLA-draped backbone is
   dropped before it costs a ProteinMPNN sequence or a Boltz-2 call.
3. **The optimised quantity is the ratio**, `buried_area_peptide` vs
   `buried_area_hla`, not raw buried area.
4. **Stabilising sampling settings stay in play per design**: screening at
   `recycling_steps` 8, and the panel tries `recyc8` then `combo`
   (`without_potentials` + `recyc8`) per design, because the best setting is
   design-specific and the two are not additive.

Unchanged: the validated flat pocket constraint
`{"binder":"D","contacts":[{"id":"C","residue_index":6}],"max_distance":6.0}`,
`step_scale` 1.638, and batch-04's scoring code called verbatim.

**Retired metrics appear nowhere.** No ipTM `specificity_gap` and no paired
same-pose MM-GBSA ddG is computed, ranked on, revived or tuned. ipTM is
recorded as a pose-quality read-out only.

## The gate, stated before it was applied

A design passes batch-06 only if it clears **all** of:

* **(a) the Asp6 direct read-out gate** — a binder Arg/Lys/His/Ser/Thr/Tyr
  heavy atom within 3.5 Å of Asp6 OD1/OD2 **and** `asp6_buried_frac` > 0;
* **(b) in ≥ 3 of 3 independent constrained draws** on the G12D state;
* **(c) the WT arm**, verbatim from `batch06_panel.WT_CRITERION`:

  > A design clears the WT arm when, over 3 independent draws per state under
  > the same sampling setting and the same pocket constraint, BOTH of:
  > (1) median `buried_area_peptide` on VVVGAGGVGK is at least **30 %** lower
  > than on VVVGADGVGK, and (2) median position-6 heavy-atom contacts (all
  > atoms, backbone included) is at least **30 %** lower. Strict separation
  > (max over WT draws < min over G12D draws, for both quantities) is recorded
  > alongside as the stronger read. **The disappearance of Asp6 side-chain
  > contacts on WT is excluded: glycine has no side chain, so that loss is true
  > by construction and carries no information.** No ipTM quantity enters this
  > criterion.

### The byte-identical-response artefact is fixed

PR #11's "n = 10" was really 6–8 distinct samples, because repeat calls with
byte-identical bodies return byte-identical structures. Here each set of 3
draws is taken as `diffusion_samples: 3` in a **single** call — three
independent diffusion trajectories that cannot be duplicates of one another by
construction. Every structure is SHA-1 hashed: **3 of 3 distinct in all 10
call-sets.** The artefact is gone.

## Stage 1–3 — design and screening

| stage | count |
| --- | --- |
| RFdiffusion backbones generated (2 waves: 96 + 168) | **264** |
| backbones reaching peptide position 6 **and** staying peptide-centred | **28** (14 + 14) |
| unique donor-anchored grafts built | **41** |
| complexes scored with constrained Boltz-2 | **41** |
| screening-draw Asp6 gate passes | **3** |
| **gate passes in 3 of 3 draws** | **0** |

Yield by length — the short binders almost never reach the aspartate:

| binder length | scored | median total area (Å²) | median peptide share | Asp6 gate passes |
| --- | --- | --- | --- | --- |
| 34 aa | 4 | 736 | 0.205 | **0** |
| 40 aa | 7 | 830 | 0.179 | **0** |
| 46 aa | 30 | 730 | 0.217 | **3** |

The peptide-centred design filter is brutal: 28 of 264 backbones survived it
(10.6 %), against 46 of 200 (23 %) for batch-05's looser filter, and the two
shortest strategies produced **no** gate pass at all.

## Result 1 — the paratope shrank, but it did **not** re-centre

This is the headline, and it is negative.

| quantity (median over scored designs) | batch-05 (n = 60) | **batch-06 (n = 41)** | change |
| --- | --- | --- | --- |
| `buried_area_total` | 960.1 Å² | **760.8 Å²** | **−21 %** |
| `buried_area_hla` | 772.9 Å² | **591.3 Å²** | **−24 %** |
| `buried_area_peptide` | 201.2 Å² | **158.2 Å²** | **−21 %** |
| **peptide share of buried area** | **0.209** | **0.213** | **+0.004 — nil** |
| max peptide share | 0.341 | 0.350 | +0.009 |

**Halving the binder shrank every part of the interface by the same ~21 %, and
moved the peptide/HLA ratio not at all.** The groove contact was not traded for
peptide contact; it was simply scaled down along with everything else. Batch-06
designs are still groove-dominated — the HLA still takes ~79 % of the buried
area, exactly as in batch-05.

That is a clean refutation of the batch-06 hypothesis as stated. A smaller
binder is not automatically a more peptide-centred one: the peptide is a
10-residue strip lying *inside* the groove, so any binder that touches it at a
sensible angle also touches the α1/α2 helices, and the ratio is set by the
geometry of the target, not by the size of the binder.

Interface charge stayed fixed where the donor bias put it: median
**0.0**, mean +0.08, range −4.0 to +5.0 (batch-05: median 0.0, mean +0.39).
Asp6 was buried at all in 21 of 41, with ≥ 1 position-6 side-chain contact in 18.

## Result 2 — the 3-of-3 Asp6 gate: **0 survivors**

The 3 screening-draw passes were advanced to the panel, each tried at `recyc8`
and then `combo`:

| design | screening share | `recyc8` | `combo` | best |
| --- | --- | --- | --- | --- |
| `bb38 -g1` | 0.321 | **0 / 3** | **0 / 3** | 0 / 3 |
| `bb54 -g2` | 0.273 | **2 / 3** | 1 / 3 | **2 / 3** |
| `bb36 -g2` | 0.146 | **2 / 3** | 0 / 3 | **2 / 3** |

**`survivors_3of3` = 0**, so by the stated gate **no design was entitled to a WT
arm**, and none is submittable. This is the third consecutive batch in which the
best designs manage 2 of 3: the read-out is still an intermittent contact, not a
weak one.

`bb38 -g1` is worth a note as the *most* peptide-centred design in the project
(share 0.32, HLA area only 290–491 Å², the smallest interfaces here) — and it
lost the aspartate entirely in 6 of 6 panel draws. Shrinking the interface far
enough to re-centre it also shrank it below the point where it reads anything.

### Replicate pose spread (binder Cα RMSD, pMHC frame, binder not re-fitted)

| design | arm @ condition | mean / max (Å) | frame (Å) |
| --- | --- | --- | --- |
| `bb54 -g2` | g12d @ combo | **3.74 / 5.10** | ≤ 0.46 |
| `bb54 -g2` | decoy @ recyc8 | 4.97 / 6.36 | ≤ 0.18 |
| `bb54 -g2` | g12d @ recyc8 | 7.26 / 7.96 | ≤ 0.29 |
| `bb54 -g2` | apo @ recyc8 | 9.90 / 13.47 | ≤ 0.41 |
| `bb54 -g2` | wt @ recyc8 | 17.70 / 24.00 | ≤ 0.33 |
| `bb38 -g1` | g12d @ recyc8 | 8.05 / 9.20 | ≤ 0.30 |
| `bb38 -g1` | g12d @ combo | 8.30 / 10.52 | ≤ 0.29 |
| `bb36 -g2` | g12d @ recyc8 | **20.81 / 26.44** | ≤ 0.49 |
| `bb36 -g2` | g12d @ combo | 20.62 / 26.44 | ≤ 0.47 |

Frame RMSD ≤ 0.49 Å throughout, so these are real. No batch-06 design matches
`bb50-g1`'s 1.58 / 3.78 Å at n = 10 — **`bb50-g1` remains the tightest pose in
the project and the reference design.** `bb36 -g2` at ~21 Å mean is two
unrelated docked modes and should be read as noise, which is exactly how its WT
diagnostic came out.

## Result 3 — the WT arm, run as an explicit diagnostic

Since 0 designs cleared (a)+(b), the WT arm was run **as a diagnostic** on the
two 2-of-3 near-misses (`scripts/batch06_wt_diagnostic.py`). **This loosens
nothing**: the criterion is the same one stated above, and the result is
reported as a diagnostic, never as a pass. It was run because it answers the
batch-06 scientific question directly — did re-centring buy *any* WT
discrimination? — for 2 calls and 0 design slots.

| design | median peptide area G12D → WT | **area loss** | median pos6 contacts G12D → WT | **pos6 loss** | would meet criterion? |
| --- | --- | --- | --- | --- | --- |
| **`bb54 -g2`** @ recyc8 | 287.1 → 117.0 Å² | **−59.2 %** | 26 → 7 | **−73.1 %** | **yes (both ≥ 30 %)** |
| `bb36 -g2` @ recyc8 | 178.3 → 175.9 Å² | −1.3 % | 12 → 17 | **+41.7 %** (anti-selective) | no |

**`bb54 -g2` is the first design in this project to show genuine wild-type
discrimination on measured geometry.** It loses 59 % of its peptide buried area
and 73 % of its position-6 contacts on VVVGAGGVGK, and the loss is not the
by-construction glycine artefact — position-6 contacts here are **all** heavy
atoms, backbone included, and the Asp6 side-chain term is excluded. On the G12D
state it holds 15/19/24 position-6 side-chain contacts; on WT, 0/0/0 while
still making 3–16 backbone-level contacts, i.e. it is still *near* position 6
and has stopped *gripping* it.

Honest caveats, which are why this is a diagnostic and not a survivor:
* it clears the Asp6 read-out in only **2 of 3** draws;
* the WT separation is **not strict** (WT draw r0 at 205.2 Å² overlaps the G12D
  range 250.4–320.2 Å²); the medians clear 30 %, the distributions touch;
* its WT pose spread is 17.7 Å mean — part of the "loss" is the binder simply
  falling off into a different mode, which is a real loss of engagement but a
  noisy one;
* its net interface charge reaches **+5.0** in one G12D draw, above the
  unchanged `NET_INTERFACE_CHARGE_MAX` of +4.0.

## Result 4 — apo and decoy, for `bb54 -g2` (also diagnostic)

| arm | buried area on HLA (3 draws) | peptide area | verdict |
| --- | --- | --- | --- |
| g12d @ recyc8 | 531.9 / 554.5 / 793.1 Å² | 250–320 Å² | on-target |
| **decoy** `AVFDRKSDAK` | 346.0 / 386.5 / 396.2 Å² | **33.2 / 34.1 / 47.6 Å²** | **clean rejection** |
| **apo** (unconstrained) | **1116.4 / 1198.2 / 1220.1 Å²** | n/a (no peptide) | **fails — buries more than on-target** |

* **Decoy: a clean rejection**, and a strong one — peptide buried area collapses
  by ~87 % (median 287 → 34 Å²) and position-6 side-chain contacts go 15–24 → 0–3.
  The binder has not merely learned to like a carboxylate; the decoy's Asp4/Asp8
  do not rescue it.
* **Apo: the PR #11 failure is unchanged.** `bb54 -g2` buries ~1116–1220 Å² on
  peptide-free HLA, *more* than the 532–793 Å² it buries on the G12D complex.
  **Stated plainly: the apo arm has no chain C, so the pocket constraint cannot
  be expressed and the call is made unconstrained, exactly as
  `CONTROL_7OW6.md` did. It is therefore NOT a like-for-like comparison with the
  three constrained arms** — nothing is pulling the binder onto the groove face
  in it. But the number is the same number PR #11 reported for every design, and
  the peptide-centred paratope did not improve it.

So even the single most encouraging design in batch-06 is selective against the
*wrong peptide* and against *no peptide at the position-6 level*, but still
drapes happily onto an empty groove.

## The `peptide_buried_frac` floor — reported, deliberately not applied

**0 of 41 screened designs reach the 0.40 floor**; 26 of 41 clear the proposed
0.18. Batch max 0.350, median 0.213.

| design | `peptide_buried_frac` (G12D panel draws) | clears 0.40? | clears 0.18? |
| --- | --- | --- | --- |
| `bb54 -g2` | 0.266 / 0.311 / 0.376 | **no** | yes |
| `bb38 -g1` | 0.264 / 0.299 / 0.248 | **no** | yes |
| `bb36 -g2` | 0.154 / 0.309 / 0.197 | **no** | partly |
| *7OW6 real TCR (reference)* | *0.224 predicted / 0.207 crystal* | ***no*** | *yes* |

The floor still rejects the real, experimentally validated 7OW6 TCR. The
recalibration to 0.18 remains **unapplied and is the user's call**;
`src/design_gates.py` is unchanged. As in batch-05, it did not decide this
batch — the 3-of-3 replication gate did, and recalibrating to 0.18 would still
leave **zero** survivors.

## Verdict

* **The batch-06 hypothesis is refuted as stated.** Shorter binders shrank the
  interface 21 % uniformly and moved the peptide/HLA ratio by 0.004. You cannot
  re-centre this paratope by shrinking it, because the peptide sits inside the
  groove and the ratio is a property of the target's geometry. The designs that
  *did* re-centre (`bb38 -g1`, share 0.32) stopped reading the aspartate
  altogether.
* **0 of 41 designs clear the gate.** This is a null result and it is reported
  as one.
* **One genuinely new finding:** `bb54 -g2` is the first design in this project
  to lose engagement on the wild-type peptide on measured geometry (−59 % area,
  −73 % position-6 contacts) *and* to reject the decoy cleanly (−87 % area).
  That is a real advance on PR #11, where no design let go of WT at all. It is
  still not submittable: 2 of 3 on the Asp6 read-out, non-strict WT separation,
  and unchanged apo over-burial.
* **`bb50-g1` remains the project's reference design** on pose tightness
  (1.58 / 3.78 Å at n = 10); nothing in batch-06 approaches it.

## Recommendation

**Do not submit on the window that opened 07:30 Taipei on Nov 4.** There are
zero gate-clearing designs; a submission today would be a 0-design or
a knowingly-failing-design submission. The Track-1 cap is 40 designs with
**0 spent** and the deadline is **Nov 8 (AoE)** — roughly four days of runway,
so there is no scarcity argument for spending a slot now.

**Aim for a submission window on Nov 6–7 Taipei**, which leaves a full day of
margin before Nov 8 AoE, and spend the intervening time on the one thread that
actually moved:

1. **Push `bb54 -g2` from 2-of-3 to 3-of-3.** It is the only design that has
   ever shown WT discrimination here. Its failure is the same intermittent-pose
   failure as every other batch, and its g12d pose is tightest under `combo`
   (3.74 Å) while its gate passes come under `recyc8` — those two facts point at
   a per-design setting that has not been tried (e.g. `combo` at more draws, or
   `recyc8` with `step_scale` moved). Cheap: ~10–15 calls, 0 design slots.
2. **A rigidifying redesign around `bb54 -g2`'s graft** (loop shortening,
   proline/disulfide placement) — attack the pose flip in sequence space, which
   is what batch-05 recommended and what nothing has yet done.
3. **Do not pursue "shrink the binder" further.** It is now measured and it does
   not work. If the groove contact is to be cut, it has to be cut by *shape*
   (a binder that approaches end-on along the peptide axis), not by *size*.
4. **Settle the `peptide_buried_frac` floor** (0.40 vs the control-calibrated
   0.18) outside deadline pressure. It has not decided any batch yet, but it
   will block the first real survivor.

**If a submission is made anyway**, the only defensible set is **1 design:
`bb54 -g2`**, submitted explicitly as a 2-of-3 read-out with documented WT and
decoy discrimination and documented apo failure. My recommendation is **not** to
spend the slot on it yet — item 1 above could convert it into a 3-of-3 design
within a day, and the same slot then buys far more.

**The submission decision is the user's alone.** Nothing has been submitted to
Proteinbase, nothing pushed to `main`, no PR merged.

## Reproducing

```
python3 scripts/run_batch06.py --resume                 # stages 1-3 (complete)
python3 scripts/batch06_panel.py --max-designs 8 --draws 3
python3 scripts/batch06_wt_diagnostic.py                # diagnostic only
python3 scripts/batch06_analyze.py
```

`logs/` is **not committed**: `*.log` is already covered by the repository's
existing `.gitignore` policy and has been since batch-01, the logs are raw
stdout duplicating `checkpoint.jsonl` / `panel_checkpoint.jsonl` (which are
authoritative), and keeping API-calling transcripts out of git is the safer
default. No `.gitignore` change was needed.
