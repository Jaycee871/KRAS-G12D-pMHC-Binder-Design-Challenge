# Submission plan — current as of 2026-11-05

**Official deadline:** Sunday 8 Nov 2026, 23:59 Anywhere on Earth (UTC−12),
which is **Monday 9 Nov 2026, 19:59 Asia/Taipei**.

Official challenge page:
https://proteinbase.com/competitions/anthropic-adaptyv-2026/challenges/pmhc

## Current project state

Batches 01–08 have been run and analysed. Batch-09's control audit and
matched-anchor panel are complete (12 calls / 36 samples). Batch-10 has now
completed the matched G12D/apo counter-screen of review ranks 4–11
(16 calls / 48 samples).

- No sequence has been submitted to Proteinbase.
- The batch-09 positive-control audit retires the 30% WT geometry-loss rule as
  a future hard selectivity gate; historical batch reports remain unchanged.
- `bb50-g1` is the lead, followed by `bb54-g2` and diversity hedge `bb25`.
- Historical apo measurements used a different docking prior from loaded pMHC;
  the completed matched regional-anchor panel still shows apo over-burial for
  all three candidates (1.34–1.85× G12D area). This remains a structural warning,
  not an affinity measurement.
- Batch-10 produced **0/8 Phase-1 advancers**. Corrected G12D Asp6 readout was
  at most 1/3 for every rank 4–11 design, and every matched apo/G12D buried-area
  ratio was >1 (1.258–2.348). By the predeclared adaptive rule, Phase 2 uses
  **0 WT/decoy calls**.
- No candidate's biological selectivity is established by these predictions.

See `results/batch-09/SUMMARY.md` and
`results/batch-09/MATCHED_ANCHOR_RESULTS.md` for the lead-panel audit, and
`results/batch-10/SUMMARY.md` for the lower-ranked matched apo counter-screen.

## Current official submission rules

The live FAQ now says:

- submit a **CSV ordered by preferred rank**;
- minimum columns: `name`, `sequence`, `molecule_class`;
- `molecule_class` is one of
  `protein`, `nanobody`, `scfv`, `fab_kappa`, `fab_lambda`;
- single-chain proteins are 10–250 aa;
- Fab is submitted as `{VH}:{VL}`;
- Track 1: submit **at least 20 and at most 40** designs;
- Tracks 2/3: submit at most **20** designs;
- Track 1: the **top 20 ranked designs that pass filtering** are screened;
- all designs must be unique, de novo and zero-shot. An existing binder may be
  used for calibration/evaluation, but not as the starting binder to modify.

The current FAQ does **not** document the old repo assumptions of “one
submission every 24 h”, “≤20 designs per file for everyone”, a mandatory 5 MB
CSV limit, or a hard “novelty ≥3/4” score. Those old assumptions are retired
unless the live submission UI explicitly reintroduces one of them.

## What the internal gates mean

The repo's gates are **research diagnostics**, not competition rules. They are
deliberately stricter than “submit anything that might bind” and should not be
silently retuned after seeing results. The control-based retirement of the WT
geometry-loss hard gate is documented explicitly in the batch-09 audit.

Two important consequences:

1. The old `peptide_buried_frac >= 0.40` threshold is not an appropriate
   absolute veto. The experimentally validated 7OW6 TCR itself is around
   0.21–0.22 by the same analysis, and the official challenge recommends the
   peptide **and the groove around it** as the epitope.
2. A failed internal gate is evidence to lower a design's rank, not permission
   to rewrite history. If a near-miss is submitted, its failure mode should be
   documented honestly in the metadata.

## Legacy-pool decision — 5 Nov

The fallback audit of batches 01–03 is now complete with **zero new model calls**.
Across 163 checkpoint rows and 144 unique historical sequences, **no sequence**
has a recorded Asp6 carboxylate-oxygen donor within 3.5 Å together with positive
Asp6 side-chain burial. Seven sequences show burial without the required oxygen
readout and are not rescued by contact count alone.

This closes the “assess older excluded batches” branch of the decision tree.
The current 11-design review pool therefore remains 9 short of the Track-1
minimum. Any expansion now requires a genuinely new mechanism; silently
reintroducing legacy failures is not permitted. See
`results/legacy-readout-audit/README.md`.

## Plan from here

### 4 Nov — infrastructure and rule repair

- Align submission tooling with the live FAQ.
- Repair stale `single_chain` molecule classes to `protein`.
- Restore CI without installing the optional OpenMM/MM-GBSA stack into the
  lightweight test environment.
- Do not spend experimental design slots.

### 5 Nov — matched counter-screen diagnostic

Before generating another large RFdiffusion batch, test the strongest and most
diverse existing candidates under a **matched docking/counter-screen setup**.
The present apo comparison is not like-for-like: G12D/WT use a peptide-position
constraint while apo has no peptide and therefore runs unconstrained.

The purpose of this diagnostic is not to weaken the apo objective. It is to ask
whether the current computational counter-screen is measuring peptide
dependence or an artefact of different docking priors.

The merged, frozen panel tests `bb50-g1`, `bb54-g2` and `bb25`, each on
G12D, WT, apo and decoy. All arms use HLA A69/A155 anchors; each design uses its
preselected sampling regime across all arms. This is 12 calls / 36 samples.
The exact settings and recovery procedure are in `docs/BATCH09_RUNBOOK.md`.

No gate is changed by this diagnostic.

### 6 Nov — mechanism decision and candidate-pool freeze

Batch-10 closed the remaining matched-apo evidence gap for ranks 4–11 and
produced no advancer. Do **not** run its WT/decoy Phase 2 and do **not** repeat
the already-refuted “smaller”, “rigidified”, generic shape-steering or helix
de-wetting strategies.

A final generation experiment is justified only if it changes the mechanism:
the binder must be geometrically dependent on the peptide itself at design
time, for example a predeclared bidentate peptide motif combining direct Asp6
carboxylate recognition with a second independent peptide contact while
explicitly disfavoring HLA-only docking. If that protocol cannot be implemented
cleanly, freeze generation and build the ranked near-miss package instead.

### 7 Nov — build the ranked submission package

For Track 1, prepare 20–40 unique sequences, with the strongest 20 first.
For Track 2/3, prepare at most 20.

Rank by a composite evidence table rather than a single ipTM:

1. reproducible direct/indirect read-out of peptide position 6;
2. WT diagnostics, interpreted under the batch-09 control audit;
3. peptide-free HLA and unrelated-peptide counter-screens;
4. pose reproducibility;
5. sequence/topology diversity and design provenance;
6. sequence liabilities and manufacturability;
7. on-target affinity/confidence only after the specificity evidence above.

Include methods/provenance and links to the relevant batch summaries.

### 8 Nov (Taipei) — submit with margin

Prefer submitting by **18:00 Asia/Taipei on 8 Nov**, more than a day before the
AoE cutoff. The official cutoff corresponds to 19:59 Taipei on 9 Nov, but that
extra day is contingency, not working time.

After submission, verify the accepted design count and saved rank order in the
Proteinbase UI and preserve the confirmation.

## Pre-flight checklist

```bash
python -m pytest -q
python src/make_submission.py <scored>.csv -o submissions/final.csv -n <N> \
  --rank-column evidence_rank --ascending
python src/validate_submission.py submissions/final.csv
```

- [ ] Reviewed `evidence_rank` is explicit (1 = best); legacy ipTM gaps do not reorder it.
- [ ] CSV is ordered best-first.
- [ ] `name`, `sequence`, `molecule_class` are present.
- [ ] molecule classes use the official vocabulary.
- [ ] 10–250 aa per single protein chain; canonical amino acids only.
- [ ] Fab only: `{VH}:{VL}`, with `fab_kappa` or `fab_lambda`.
- [ ] unique sequences and unique names.
- [ ] Track 1: 20–40 submitted; Track 2/3: ≤20.
- [ ] designs are de novo/zero-shot; known binders were calibration controls,
      not starting templates.
- [ ] no API key, token, private path or secret in the CSV or metadata.
- [ ] methods state negative/null results as well as positive evidence.

## Submission log

| date (Taipei) | track | designs | accepted? | notes |
| --- | --- | ---: | --- | --- |
| — | — | 0 | — | nothing submitted as of 2026-11-04 |
