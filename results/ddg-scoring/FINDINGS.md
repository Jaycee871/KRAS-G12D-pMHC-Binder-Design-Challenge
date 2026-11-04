# Paired same-pose MM-GBSA ddG: the 7OW6 control, and why the metric is dead

**Date:** 2026-11-03 (Asia/Taipei)
**Module:** `src/ddg_scoring.py` · **Tests:** `tests/test_ddg_scoring.py` (10 passed)
**Scripts:** `scripts/ddg_control_7ow6.py`, `scripts/ddg_control_diagnostics.py`,
`scripts/ddg_groove_control.py`, `scripts/ddg_rescore_designs.py` (**written, never run on designs**)
**Raw data:** `results/ddg-scoring/control_7ow6_ddg.json`,
`control_diagnostics.jsonl`, `groove_control.jsonl`

## Verdict, up front

**The paired same-pose MM-GBSA ddG fails its 7OW6 control and is not a
validated ranking quantity. No design has been or will be ranked on it.**

The positive control comes out with the **wrong sign** under every one of the
four parameter corners tested (+1.27, +2.32, +7.13, +9.60 kcal/mol, where
negative is required). The replicate spread (SD 0.8–8.7 kcal/mol) is larger
than the effect being measured. Both candidate remedies — tighter minimisation
and a raised solute dielectric — were tested and **both made it worse**. There
is no defensible fix, so none was adopted, and the metric was not tuned until
the control passed. That is the mistake `specificity_gap` made and it is not
being repeated.

This is the **second** ranking quantity to die on the same control. See
[Where this leaves ranking](#where-this-leaves-ranking).

## Method chosen, and why

`results/validation/CONTROL_7OW6.md` killed

```
specificity_gap = ipTM(G12D) - max(ipTM(WT), ipTM(apo), ipTM(decoy))
```

because each ipTM comes from an **independent** Boltz-2 co-folding run, and
byte-identical re-runs move ipTM(G12D) by 0.06–0.23 — far more than one side
chain is worth. The diagnosis was *run-to-run variance swamping a one-side-chain
signal*, so the replacement was built to cancel that variance by construction:

```
ddG_sel = dG_bind(G12D pose) - dG_bind(same pose, peptide position 6 -> Gly)
negative favours G12D
```

* **One pose, two states.** Both terms are evaluated on the *same* predicted
  complex. The docking geometry, the HLA groove, β2m and the other nine peptide
  side chains are shared atom-for-atom and cancel exactly in the subtraction.
* **Energy model:** Amber ff14SB (`amber14-all.xml`) + GBn2 implicit solvent
  (`implicit/gbn2.xml`, ACE non-polar surface term), `NoCutoff`.
  Single-trajectory MM-GBSA, so `dG_bind = E_complex − E_binder − E_target`
  carries interaction and desolvation but no conformational entropy — defensible
  *only* for a difference of two dG_bind values that share a pose.
* **Local truncation** to residues with a heavy atom within 16 Å of peptide
  position 6; the kept-residue set is computed once on the G12D coordinates and
  reused for the WT state, so the two states are provably identical in
  composition (enforced in code and unit-tested).
* **Relaxation:** restrained L-BFGS; heavy atoms beyond 8 Å of position 6
  harmonically restrained (10 kcal/mol/Å²), everything inside free.
* **Replicates:** Gaussian jitter (σ = 0.01 nm) on the free atoms, applied
  *identically to both states of a pair*. Replicate 0 is unjittered.
* **Asp→Gly is a pure deletion** — no rotamer is guessed, no side chain is
  rebuilt. This was chosen deliberately: it removes the one step that would have
  introduced a modelling free parameter.

### What was rejected, and why

| rejected | reason |
| --- | --- |
| ipTM `specificity_gap` | wrong sign on 7OW6, and not reproducible run-to-run (`CONTROL_7OW6.md`) |
| independent re-docking of the WT complex | reintroduces exactly the run-to-run variance the paired form exists to cancel |
| full-complex MM-GBSA, no truncation | ~10× the cost for terms that cancel; 16 Å truncation is provably state-identical |
| side-chain rebuilding / repacking of the WT state | Asp→Gly is a deletion; rebuilding would add a rotamer guess to a 1 kcal/mol measurement |
| raising `soluteDielectric` to 4 | tested — see diagnostics; it inflates `dG_bind` ~4× and makes the control worse |
| conformational-entropy (normal-mode) term | not affordable, and the paired form was supposed to make it cancel |

## Control panel 1 — the 7OW6 TCR

`results/ddg-scoring/control_7ow6_ddg.json`, 5 replicates for the two positives,
3 for the negatives, default knobs (eps_in = 1, 500 L-BFGS iterations).
**Negative is the required sign for a G12D-selective binder.**

| system | kind | ddG_sel | SD | replicates | n_res | n_frag | expected | result |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `crystal_7ow6` — deposited crystal pose of the experimentally G12D-selective Poole TCR, the ground truth | positive | **+1.27** | 0.81 | +0.41, +2.06, +1.89, +0.40, +1.59 | 153 | 14 | clearly negative | **WRONG SIGN** |
| `predicted_7ow6` — pocket-constrained Boltz-2 prediction of the same complex | positive | **+0.81** | 1.90 | +1.87, −2.49, +1.15, +1.25, +2.29 | 151 | 15 | clearly negative | **WRONG SIGN, SD > signal** |
| `crystal_offsite_p3` — same crystal, paired mutation moved to peptide position 3 | negative | −0.83 | 1.58 | −0.82, +0.75, −2.41 | 153 | 17 | ≈ 0 | indistinguishable from zero *and* from the positive |
| `decoy_7ow6` — same TCR on the unrelated `AVFDRKSDAK` decoy, paired at position 6 | negative | **−3.53** | **0.12** | −3.40, −3.58, −3.63 | 156 | 12 | ≈ 0 | **strongest and most confident "G12D-favouring" number in the panel** |
| `offpeptide_design` — a batch-01 design docked away from the peptide | negative | — | — | — | — | — | raise | correctly raised `DdgError: no binder residue within 16 Å of peptide position 6` |

The panel is **inverted**: the two positives are positive, the strongest
negative value belongs to a decoy that cannot have a G12D preference, and the
positive and the off-site negative are not separable given their spreads
(+1.27 ± 0.81 vs −0.83 ± 1.58; SEM 0.36 and 0.91).

Only the error-handling row behaved: a binder that does not reach the site
*raises* rather than quietly returning 0, which is the one thing batch-01 could
not distinguish.

## Control panel 2 — the groove, i.e. can the energy function see a carboxylate at all?

A failing control is ambiguous: the metric could be ill-posed *for this epitope*,
or the MM-GBSA energy function could be blind to this chemistry. A contact
census on the deposited crystal separates them. Computed here from
`data/raw/7OW6.pdb`, peptide `VVVGADGVGK`, site residue ASP 6:

| probe → Asp6 | nearest heavy atom | atoms < 4 Å | atoms < 5 Å |
| --- | --- | --- | --- |
| TCR (chains D+E) → Asp6 **side chain** (CB/CG/OD1/OD2) | **5.73 Å** | **0** | **0** |
| TCR (chains D+E) → Asp6 **backbone** (N/CA/C/O) | 3.37 Å | 1 | 7 |
| HLA chain A → Asp6 **side chain** | **2.95 Å** (Arg114 NH2) | 6 | 16 |
| HLA chain A → Asp6 backbone | 3.02 Å | 2 | 9 |

**The TCR never touches the Asp6 side chain.** Every "position-6 contact" the
pipeline counts for 7OW6 is to the peptide *backbone*, which both states of a
pair share atom-for-atom and which therefore cancels exactly. This is the
structural restatement of the caution already in `docs/LITERATURE.md`: Poole et
al. 2022 report *few discernible differences in direct contacts* and attribute
the selectivity to **indirect electrostatic** effects.

So the groove was used as a positive control *inside the same crystal*: take HLA
chain A as the "binder" and the peptide as the target, where the Arg114–Asp6
salt bridge is real, and run the identical calculation
(`scripts/ddg_groove_control.py`, `groove_control.jsonl`, 3 replicates each):

| system | eps_in | ddG_sel | SD | replicates |
| --- | --- | --- | --- | --- |
| `groove_asp6_eps1` — Arg114 NH2 2.95 Å from Asp6 OD2 | 1 | **−7.42** | 0.75 | −7.83, −7.88, −6.56 |
| `groove_asp6_eps4` — same, with the buried-charge "remedy" | 4 | −46.32 | 1.95 | −44.82, −48.52, −45.61 |
| `groove_offsite_p3_eps1` — **negative**: Val3, a buried apolar anchor, no salt bridge | 1 | **−7.32** | 4.46 | −7.81, −2.63, −11.52 |

Two things fall out of this, and the second is the one that kills the metric.

1. **The energy function is not blind.** Given a real 2.95 Å salt bridge it
   returns −7.4 ± 0.75 kcal/mol, right sign, well clear of its own noise.
2. **It cannot tell that signal from plain side-chain burial.** The matched
   negative — Val3, an apolar groove anchor with no salt bridge, no charge and
   no hydrogen bond — returns **−7.32**, the same magnitude. The quantity being
   computed is not "is this carboxylate specifically read", it is "how much
   binding energy does whatever side chain sits at this position contribute by
   being buried at all".

## Control panel 3 — the diagnostics

`scripts/ddg_control_diagnostics.py --replicates 3 --workers 4`,
`control_diagnostics.jsonl`. Two separable hypotheses for the inverted panel,
run as a 2×2 on the ground-truth crystal plus the decoy at the best corner:

* **A — buried-charge over-penalty.** GB at eps_in = 1 over-pays desolvation of
  an interfacial charge; here the charge *is* the thing being measured. Remedy:
  eps_in = 4.
* **B — unconverged minimisation.** `dG_bind` ≈ −75 and the ddG ≈ 1: a 1.5 %
  difference of two large numbers, at 500 iterations / 1 kJ/mol/nm tolerance on
  a ~2500-atom system. Remedy: 5000 iterations / 0.05.

| run | eps_in | iterations | ddG_sel | SD | replicates | dG_bind (G12D) |
| --- | --- | --- | --- | --- | --- | --- |
| **baseline** `crystal_7ow6` | 1 | 500 | **+1.27** | 0.81 | 5 reps, +0.40 … +2.06 | ≈ −74 |
| `crystal_eps1_tight` — B alone | 1 | 5000 | **+2.32** | 1.22 | +1.19, +3.61, +2.15 | ≈ −73 |
| `crystal_eps4_loose` — A alone | 4 | 500 | **+9.60** | 7.29 | +11.79, +1.47, +15.56 | ≈ −304 |
| `crystal_eps4_tight` — A+B, the best case for the metric | 4 | 5000 | **+7.13** | 4.44 | +12.20, +3.95, +5.23 | ≈ −307 |
| `decoy_eps4_tight` — discrimination at the best corner | 4 | 5000 | −2.29 | **8.70** | −9.86, +7.22, −4.24 | ≈ −305 |

**Both hypotheses are rejected.**

* **B fails.** A 10× tighter minimisation did not restore the sign and did not
  reduce the scatter — it moved the control *further* wrong (+1.27 → +2.32) and
  raised the SD (0.81 → 1.22). The spread is therefore **not** numerical noise
  from stopping short of a stationary point; the jittered replicates are
  relaxing into genuinely different local minima, and the restrained-minimisation
  basin is not unique. More compute does not converge this.
* **A fails, badly.** eps_in = 4 inflates `dG_bind` from ≈ −74 to ≈ −305
  kcal/mol — a ~4× magnification, and the same artefact appears in the groove
  control (−7.42 → −46.32). The reason is mechanical: in OpenMM's GB
  implementation `soluteDielectric` rescales the **polar desolvation** term by
  (1/ε_in − 1/ε_out) but leaves the direct Coulomb interaction in
  `NonbondedForce` computed at ε = 1. Raising ε_in therefore removes most of the
  desolvation penalty while leaving the full unscreened attraction, which is not
  a physical correction at all. It also *tripled* the error on the control
  (SD 0.81 → 7.29). This knob is not a remedy and must not be used.
* **Discrimination at the best corner is gone entirely.** `decoy_eps4_tight`
  returns −2.29 ± **8.70** over three replicates spanning −9.86 to +7.22. At the
  corner most favourable to the metric, a single replicate carries no
  information whatsoever.

Across the 2×2 the ground-truth crystal spans **+1.27 to +9.60 kcal/mol** — an
8.3 kcal/mol range produced entirely by two numerical knobs that have nothing to
do with the biology, for an effect that should be on the order of 1 kcal/mol.

## Diagnosis — why it failed

Three causes, in descending order of importance. They are independent, and the
first alone is fatal.

**1. The metric measures side-chain burial, not specific recognition.**
`ddG_sel` as constructed is "the binding contribution of whatever side chain sits
at the queried peptide position". Deleting any buried side chain to Gly removes
vdW contact and non-polar burial, so the quantity is systematically negative for
any buried residue regardless of whether the binder discriminates. The evidence
is direct and threefold: the decoy's Lys6→Gly gives −3.53 ± 0.12 (the tightest
number in the panel) because Lys6 is a large buried side chain, not because a
decoy peptide has a G12D preference; the groove's apolar Val3→Gly gives −7.32,
matching the real Arg114–Asp6 salt bridge's −7.42; and the only genuinely
specific term in the whole problem — electrostatic complementarity to a
carboxylate — never separates from that background. **A quantity that scores a
Lys, a Val and an Asp alike is not a selectivity metric.** The decoy row is also
mis-specified in a way worth stating plainly: with no G12D in the decoy peptide,
`pos6→Gly` asks a different physical question (Lys6 vs Gly6) than it does in the
KRAS context, so that row is weaker evidence of failure than it first appears.
It is not weaker evidence of the burial confound, which is the real problem.

**2. The 7OW6 positive control is ill-posed for any single-pose contact-based
energy.** No TCR atom comes within 5 Å of the Asp6 side chain (nearest 5.73 Å).
What `ddG_sel` is actually computing for this complex is the GB self-desolvation
of a **solvent-exposed** carboxylate — slightly more costly in the complexed
state than free, hence the consistent small **positive** value — plus noise.
There is no direct interaction term for it to find, because there is no direct
interaction. Poole's TCR is selective by indirect electrostatics acting at
>5.7 Å, and a single minimised structure in implicit solvent has no term for
that: it contains no reorganisation free energy, no explicit water
restructuring, no ionic screening, and no conformational entropy. This also
means the project currently has **no positive control that any direct-contact
metric could pass**, which is a problem in its own right.

**3. The numerical floor sits above the signal.** The required effect is of
order 1 kcal/mol; the method's replicate SD is 0.8–8.7 kcal/mol and its
parameter sensitivity is 8.3 kcal/mol. The pairing cancelled the *Boltz-2*
run-to-run variance successfully — that part worked — but it replaced it with
minimisation-basin variance of the same order. Tightening convergence did not
help, which says the basins are real, not numerical. Contributing: the local
truncation is not identical across systems (`n_fragments` 12–17, `n_residues`
151–156), so the dangling-terminus and fragment-interaction error differs
between the systems being compared, though this is a second-order effect next to
causes 1 and 2.

## No fix was applied

The brief allowed one defensible fix to be tested and re-controlled. The two
candidates the diagnostics were designed to test were both tested and both
rejected by the data above. No third fix is indicated by the diagnostics, and
continuing to turn knobs until the control passes is precisely how
`specificity_gap` survived as long as it did. **The metric is retired.**
`scripts/ddg_rescore_designs.py` is committed but has **never been run against a
design**, and must not be until a control passes.

## Where this leaves ranking

Two ranking quantities have now failed the same 7OW6 control:

| quantity | control result | status |
| --- | --- | --- |
| `specificity_gap` = ipTM(G12D) − max(ipTM(WT), ipTM(apo), ipTM(decoy)) | −0.006 on a known selective binder; re-run spread 0.06–0.23 flips a +0.30 gap | **not validated** |
| `ddG_sel` paired same-pose MM-GBSA | +1.27 ± 0.81 (wrong sign) baseline; wrong sign at all four corners; decoy −3.53 ± 0.12 | **not validated** |

The honest reading is that **the project has no validated continuous ranking
quantity, and probably cannot obtain one before 2026-11-08**. The shared cause
is not bad luck with two methods; it is that the thing being asked for — a
one-side-chain, partly *indirect* electrostatic preference — is below the
resolution of both a co-folding confidence score and a single-structure implicit-
solvent energy.

### What to rank by instead — recommendation

Stop trying to score selectivity and start *requiring the structure that
produces it*. The groove control is the usable result here: when a design makes
a real 2.95 Å salt bridge to the Asp6 carboxylate, the energy function resolves
it at −7.4 ± 0.75 kcal/mol, comfortably above its own noise. 7OW6 fails the
control because the TCR does not do that; a design that **does** do it is in the
regime where the physics is measurable.

1. **Hard gate, not a score — direct Asp6 read-out.** A binder donor (Arg/Lys/
   His/Ser/Thr/Tyr) within 3.5 Å of Asp6 OD1/OD2, with `asp6_buried_frac > 0`.
   This is the batch-02 Asp6 read-out gate; it should become the primary filter.
   Note 7OW6 itself would fail this gate (0.00 buried, no H-bond) — that is
   intentional. We cannot reproduce indirect selectivity, so we should design for
   the direct mechanism we *can* verify.
2. **Rank the survivors on geometry and chemistry only**, all of which are
   measured on a structure rather than inferred from an energy:
   `peptide_buried_frac` (≥ 0.40), count of contacts to the Asp6 **side chain**
   (not the backbone — 7OW6 shows that distinction is everything), net interface
   charge kept out of strongly negative territory (Xie et al. 2025), then
   topological diversity.
3. **Demote ipTM to a pose-quality filter.** Use ipTM(G12D) as a threshold on
   whether a prediction is confident at all. Never as a gap, never as a
   tie-breaker.
4. **Use pose reproducibility as a real ranking signal.** Re-run each candidate's
   constrained prediction N times and rank on C-alpha RMSD spread of the binder.
   This is cheap, it is measurable, and it is exactly the axis on which
   `specificity_gap` was shown to be worthless — a design whose pose is stable
   across re-runs is a better bet than one whose pose is not, independent of any
   energy.
5. **Use `ddG_sel` only as a binary sanity check, clearly labelled unvalidated:**
   for designs that *do* form a direct salt bridge, a value more negative than
   about −3 kcal/mol is consistent with the groove control. It is not a ranker
   and it must not appear in a sort key or in a submission column.
6. **Spend the remaining submission windows on diversity, not on a fine
   ranking.** With no validated ranker, a fine ordering of the top 20 is noise.
   Submitting structurally distinct designs that all satisfy the hard gate is
   strictly better than submitting 20 variations ordered by a quantity that
   failed its control.

### What I would try next, given time

* **A positive control that a direct-contact metric could pass.** The single
  most valuable next step, and currently missing. Find a pMHC–TCR or
  pMHC–binder pair whose published selectivity rests on a *direct* side-chain
  contact to the discriminating residue, and control against that. Until such a
  control exists, no energetic metric in this repo can be validated at all, and
  that includes any future one.
* **Ensemble-averaged MM-GBSA**, 10–50 snapshots from a short restrained MD per
  state rather than one minimised structure. This attacks cause 3 directly: the
  diagnostics show the single-structure result depends on which local basin the
  minimiser lands in, and averaging over a basin ensemble is the standard
  remedy. It would not touch causes 1 or 2.
* **Alchemical free energy (TI or FEP) for Asp6→Gly in explicit solvent**,
  bound and free legs. This is the only approach in reach that contains the
  terms cause 2 identifies as missing — explicit water reorganisation and ionic
  screening of an indirect electrostatic effect. It is also far beyond the
  compute budget before the deadline, and should be scoped as post-competition
  work rather than attempted now.
* **Decompose rather than aggregate.** Report the per-term breakdown
  (Coulomb / GB polar / non-polar surface) of the Asp6 contribution instead of a
  single ddG. The burial confound in cause 1 is visible in the decomposition —
  a Val and an Asp do not look alike term-by-term even when their totals match —
  so a decomposed read-out might discriminate where the total does not. Cheap
  to add to `src/ddg_scoring.py`; worth doing before any further ensemble work.
