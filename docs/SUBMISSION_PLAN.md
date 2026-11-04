# Submission plan — current as of 2026-11-04

**Official deadline:** Sunday 8 Nov 2026, 23:59 Anywhere on Earth (UTC−12),
which is **Monday 9 Nov 2026, 19:59 Asia/Taipei**.

Official challenge page:
https://proteinbase.com/competitions/anthropic-adaptyv-2026/challenges/pmhc

## Current project state

This is no longer a scaffold. Batches 01–08 have been run and analysed.

- RFdiffusion and ProteinMPNN have produced multiple de novo design batches.
- Hosted Boltz-2 has been used for G12D, WT and off-target diagnostics.
- No sequence has been submitted to Proteinbase yet.
- No candidate currently clears the project's full internal
  **Asp6 reproducibility + pooled WT-discrimination** gate.
- The central unresolved failure is peptide-free HLA: candidates tested on apo
  HLA still show substantial groove binding.
- `bb54-g2` remains an informative near-miss, but its pooled WT peptide-area
  loss fell to 29.4%, just below the pre-registered 30% criterion in batch-08.

See `results/batch-08/SUMMARY.md` for the latest controlled result.

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
silently retuned after seeing results.

Two important consequences:

1. The old `peptide_buried_frac >= 0.40` threshold is not an appropriate
   absolute veto. The experimentally validated 7OW6 TCR itself is around
   0.21–0.22 by the same analysis, and the official challenge recommends the
   peptide **and the groove around it** as the epitope.
2. A failed internal gate is evidence to lower a design's rank, not permission
   to rewrite history. If a near-miss is submitted, its failure mode should be
   documented honestly in the metadata.

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

Prioritise:
- `bb54-g2` (best mutation-specific evidence, but pooled WT near-miss);
- `bb50-g1` (tightest pose reference);
- the strongest diverse 3-of-3/2-of-3 Asp6 readers from earlier batches;
- one or two peptide-centred batch-08 controls.

No gate is changed by this diagnostic.

### 6 Nov — freeze the candidate pool

If the matched diagnostic produces a clear improvement, make a small targeted
batch-09 around that mechanism only. Do **not** repeat the already-refuted
“smaller”, “rigidified”, generic shape-steering, or helix de-wetting strategies.

If it does not improve the evidence, stop generating backbones and rank the
best existing near-misses. Additional RFdiffusion sampling without a new
mechanism is unlikely to add information before the deadline.

### 7 Nov — build the ranked submission package

For Track 1, prepare 20–40 unique sequences, with the strongest 20 first.
For Track 2/3, prepare at most 20.

Rank by a composite evidence table rather than a single ipTM:

1. reproducible direct/indirect read-out of peptide position 6;
2. pooled WT discrimination;
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
python src/make_submission.py <scored>.csv -o submissions/final.csv -n <N>
python src/validate_submission.py submissions/final.csv
```

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
