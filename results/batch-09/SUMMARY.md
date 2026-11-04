# Batch-09 — control-calibrated selectivity audit

**Date:** 2026-11-04 (Asia/Taipei)  
**Audit compute spent:** 0 RFdiffusion, 0 ProteinMPNN, 0 Boltz-2 calls.
**Subsequent matched panel:** 12 Boltz-2 calls / 36 structures; see
[`MATCHED_ANCHOR_RESULTS.md`](MATCHED_ANCHOR_RESULTS.md).
**Track-1 slots spent:** 0.  
**Purpose:** audit the *decision rule* before generating any more binders.

## Headline

The batch-06/07/08 WT geometry-loss criterion is **retired as a hard
selectivity gate for future decisions**.

That criterion required both median peptide buried area and median position-6
heavy-atom contacts to fall by at least 30% on WT. It was pre-registered before
batch-06 and was therefore applied honestly in batches 06–08, but it was never
validated against the known positive control.

The known G12D-selective JDIa41b1 TCR is the counterexample.

Poole et al. (Nature Communications 2022, DOI
10.1038/s41467-022-32811-1) report:

- G12D affinity: **743 ± 18 pM**;
- WT affinity: **3.0 µM**;
- affinity window: **>4000-fold**;
- the TCR variable-domain pose on G12D vs WT is almost identical
  (**0.29 Å C-alpha RMSD**);
- there are few discernible differences in the direct TCR-pHLA contacts;
- selectivity is driven mainly by **indirect electrostatic interactions**.

The structural mechanism is explicit: peptide D6 forms a salt bridge with
HLA-A*11 R114, is stabilised by HLA D116, and hydrogen-bonds HLA Q70. Those
interactions are absent for WT G6. The buried D6 changes the electrostatic
surface seen by the TCR even though the TCR does not need to geometrically
"let go" of WT.

Therefore a rule that requires a selective binder to lose 30% of its WT
footprint is not control-calibrated. A real >4000-fold selective binder need
not satisfy that geometry pattern.

**Historical results are not rewritten.** The batch-06/07/08 files remain
unchanged and continue to report the criterion they actually used. This audit
only changes how future work interprets that criterion.

---

## What remains valid

The project already established two other important control results.

### 1. Boltz-2 pose geometry is useful

`results/validation/CONTROL_7OW6.md` showed that constrained Boltz-2 places
the real 7OW6 TCR close to the crystal pose and reproduces buried-area geometry
well. So Boltz-2 remains useful for:

- whether a binder reaches the intended epitope;
- whether it reproducibly adopts one docked mode;
- direct donor/Asp6 geometry;
- peptide/HLA interface partitioning;
- obvious decoy cross-reactivity.

### 2. Boltz-2 confidence is not a selectivity score

The real selective TCR gives essentially no useful G12D-vs-WT ipTM gap, and
design-regime repeats fluctuate more than the one-residue signal. The ipTM
specificity gap remains retired.

### 3. Paired single-pose MM-GBSA is not a selectivity score

`results/ddg-scoring/FINDINGS.md` showed the paired Asp6->Gly MM-GBSA metric
has the wrong sign on the positive control and cannot separate specific
electrostatics from generic side-chain burial. It remains retired.

---

## New evidence policy

For the remainder of the challenge, no quantity becomes a hard selectivity gate
unless it survives an appropriate positive control.

### Hard evidence

1. **Reproducible on-target epitope engagement.**
   A design should repeatedly reach the intended pMHC surface rather than flip
   between unrelated poses.

2. **A mechanistically interpretable G12D read-out.**
   For this de novo campaign the intentionally chosen mechanism is a direct
   donor to the Asp6 carboxylate, because that chemistry is measurable even
   though the natural TCR uses an indirect electrostatic mechanism.

3. **No demonstrated decoy-carboxylate promiscuity.**
   A binder that simply salt-bridges another Asp/Glu on an unrelated
   HLA-A*11 peptide is not a G12D reader.

4. **No demonstrated catastrophic peptide-independent behaviour.**
   Apo remains an important safety diagnostic, but the historical apo arm is
   not like-for-like because G12D/WT are peptide-constrained and apo is
   unconstrained. Until a matched anchor experiment is run, apo over-burial is
   evidence of concern, not a calibrated quantitative selectivity threshold.

### Diagnostic / ranking evidence, not hard gates

- magnitude of WT peptide buried-area loss;
- magnitude of WT position-6 contact loss;
- WT pose displacement;
- ipTM differences;
- the historical `peptide_buried_frac >= 0.40` floor.

The 0.40 peptide-share floor is especially unsuitable as an absolute gate:
the experimentally validated 7OW6 TCR is ~0.21–0.22 by the same analysis.

---

## Existing-candidate re-audit

### 1. `bb50-g1` — **priority A**

Why it stays:

- direct Asp6 read-out: **10/10**;
- donor distance: **2.53–3.11 Å**;
- Asp6 side-chain contacts: **4–7**;
- best replicated pose in the project:
  **1.58 Å mean / 3.78 Å max** binder C-alpha RMSD;
- decoy is clean: peptide footprint collapses and position-6 engagement falls
  to 0–2 contacts.

Why it is not called selective:

- the WT structure did not geometrically collapse;
- the historical apo diagnostic buried ~1454 Å² on peptide-free HLA;
- that apo comparison is not matched to the peptide-constrained on-target
  prior.

**Decision:** strongest existing sequence for the next matched counter-screen.

### 2. `bb54-g2` — **priority A-/B+**

Why it stays:

- under r3np, the direct Asp6 read-out became highly reproducible
  (pooled **11/12** by batch-08);
- it repeatedly rejects the unrelated decoy;
- it is the only later-batch design that showed a substantial WT geometry
  change in some calls.

Why its old "29.4% vs 30%" failure is no longer decisive:

- the 30% WT geometry-loss criterion itself is not positive-control validated;
- a real >4000-fold selective TCR can bind WT in almost the same pose.

Why it remains below bb50-g1:

- its G12D pose under the gate-favouring setting is much less stable;
- the apo diagnostic remains poor;
- direct-readout reproducibility is slightly below bb50-g1.

**Decision:** keep as the second mechanism/topology family.

### 3. `bb25` — **priority B**

Why it stays:

- direct Asp6 read-out: **9/10**;
- extremely short donor contacts (down to ~1.94 Å);
- sequence/topology diversity relative to bb50-g1 and bb54-g2.

Why it is lower priority:

- one clear mode flip in ten;
- decoy behaviour is mixed;
- historical apo burial is high.

**Decision:** retain as a diversity hedge, not as the lead.

### 4. `bb135-g1` — **deprioritise**

Despite a 10/10 direct Asp6 read-out, it forms a salt bridge to the unrelated
decoy peptide's **Asp8 in 3/4 draws** while retaining a large peptide
footprint. That is measured carboxylate promiscuity.

**Decision:** do not advance unless redesigned specifically to remove the
decoy-Asp8 route.

### 5. `bb54_rigid-gp-p2` and batch-08 shape-steered designs — **retire**

The rigidified design's 3/3 event did not reproduce. Batch-08 produced
0/36 new designs with a reproducible 3/3 Asp6 read-out. Shrinking, generic
rigidification, shape steering and helix de-wetting are therefore still
measured dead ends.

---

## Predeclared next experiment — matched apo / WT regional anchor

**Execution update:** this panel is now complete. All three retain greater
apo buried area under the matched regional anchor. The prospective rationale
and interpretation branch below are preserved; measured outcomes are in
`MATCHED_ANCHOR_RESULTS.md`.

The highest-value new computation is **not another backbone-generation batch**.

Run a small diagnostic on `bb50-g1`, `bb54-g2`, and `bb25` in which
G12D, WT and apo are all given the **same HLA-groove regional anchor**, instead
of constraining G12D/WT to peptide C6 while leaving apo unconstrained.

This does not replace the peptide-specific Asp6 experiment. It asks one narrow
question:

> When the docking prior is matched, does peptide-free HLA still attract the
> binder more strongly / bury more surface than peptide-loaded HLA?

The regional anchor must be fixed from the known pHLA geometry before looking
at the three design results and applied identically to all arms. No result from
this diagnostic is allowed to retune an existing historical gate.

If the matched apo failure remains, stop trying to rescue these scaffolds and
prepare a diversity-ranked Track-1 submission from the best existing evidence.
If it disappears, the old apo result was confounded by unequal docking priors
and the top three deserve a final focused rescore.

---

## Submission implication

Track 1 screens the top 20 ranked designs. Fine ordering by a metric that failed
its positive control is false precision.

The final 20–40 should therefore be ordered by **evidence tiers**:

1. reproducible direct Asp6 engagement;
2. clean unrelated-peptide behaviour;
3. matched peptide-dependence evidence once available;
4. pose reproducibility;
5. topology/sequence diversity;
6. on-target confidence only after the above.

No design is being submitted in this batch.
