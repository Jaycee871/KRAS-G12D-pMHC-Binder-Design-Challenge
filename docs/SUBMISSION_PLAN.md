# Submission plan to the Nov 8 deadline

**Deadline: Sunday 8 Nov 2026, 23:59 Anywhere on Earth** = 9 Nov 11:59 UTC
= **9 Nov 19:59 Asia/Taipei**. Plan to be finished a full day early; AoE is a
safety margin, not a schedule.

**Status: no designs generated, nothing submitted. All batches below are
planned, none are done.**

## The two limits that shape the schedule

1. **Cadence:** one submission per 24 hours, ≤ 20 designs per submission.
   From 2 Nov that leaves about **6 submission slots** before the deadline.
2. **Per-participant design caps** (these are totals, *not* per submission):
   * Track 1 — up to **40 designs**
   * Track 2 — up to **20 designs**
   * Track 3 — up to **20 designs**

The caps bind harder than the cadence. Six slots × 20 designs is 120 designs,
but only 80 can ever be spent. So batches must be **small and deliberate**:
roughly 6–10 designs each, keeping reserve capacity for the last two days when
the modelling is best. Never burn 20 Track-1 slots on batch 1 — half the budget
would be gone before a single scoring run has been sanity-checked.

Proposed budget: Track 1 → 8 + 8 + 12 + 12 = 40. Track 2 → 6 + 14 = 20.
Track 3 → 6 + 14 = 20.

## Day by day

### Mon 3 Nov — Batch 0: the safety batch (small, early, deliberate)
* Submit **4–6 designs**, the best available from a first fast
  RFdiffusion → ProteinMPNN → AF2-monomer → AF-Multimer pass.
* Purpose is **not** score. It is to prove the whole path end to end: file
  parses, columns accepted, novelty check (≥ 3/4) passes, no length or
  alphabet rejection, and the 24-hour clock behaves as expected.
* Run `python src/validate_submission.py` before uploading, every time.
* Record the exact upload timestamp (Taipei + UTC) below — the next slot opens
  24 h after it.

### Tue 4 Nov — Batch 1: first real Track 1 batch (8 designs)
* Full pipeline with the hard gates from `docs/PIPELINE.md`
  (ipTM G12D ≥ 0.75, WT and apo ≤ 0.50, gap ≥ 0.25).
* Cluster at 70 % identity, submit cluster representatives only — topological
  spread matters more than squeezing in near-duplicates.
* Carry any novelty-flagged design back for redesign rather than resubmitting.

### Wed 5 Nov — Batch 2: scaffold diversification (8 Track 1 + 6 Track 2)
* Add BindCraft-derived backbones as an independent topology prior.
* Introduce the antibody formats: nanobody and scFv written `{VH}:{VL}`.
* Explicitly include two designs that over-weight peptide position 6 contacts
  even at some cost in predicted affinity — specificity ranks above KD.

### Thu 6 Nov — Batch 3: specificity hardening (12 Track 1 + 6 Track 3)
* Re-score every surviving candidate against **all three** targets with the
  second predictor (Boltz if Multimer was primary, or vice versa) and keep only
  designs where both predictors agree on the gap.
* Drop anything whose buried interface is < 40 % on the peptide chain, however
  good its G12D ipTM — those are groove-rim binders and will hit WT and apo.

### Fri 7 Nov — Batch 4: best-of (12 Track 1 + 14 Track 2)
* Spend the remaining Track 1 budget on the highest `specificity_gap` designs,
  plus mutational micro-variants (1–3 substitutions) around the two best
  interfaces to hedge prediction error.
* Full Track 2 top-up.

### Sat 8 Nov — Batch 5: final top-up (14 Track 3) — **submit by 18:00 Taipei**
* Last remaining capacity only. Submit the day before the AoE cut-off, not
  into it.
* After uploading: confirm the submission appears on the competition page,
  screenshot the confirmation, and fill in the log below.

### Sun 9 Nov — buffer, no new work
* Reserved for a failed upload, a site outage, or a novelty rejection on
  batch 5. If everything went in on Sat, do nothing.

## Pre-flight checklist (run before every upload)

```bash
python -m pytest -q
python src/make_submission.py <scored>.csv -o submissions/batchN.csv -n <N>
python src/validate_submission.py submissions/batchN.csv
```

* [ ] ≤ 20 designs, each chain 10–250 aa, canonical amino acids only
* [ ] scFv/Fab written `{VH}:{VL}`, exactly one colon
* [ ] file ≤ 5 MB, `name` and `sequence` columns present, names unique
* [ ] remaining track budget checked against the caps above
* [ ] ≥ 24 h since the previous upload
* [ ] novelty expected ≥ 3/4; designs are de novo, not grafted from 7OW6's TCR
* [ ] no API key anywhere in the file or in the repo

Submit page:
<https://proteinbase.com/competitions/anthropic-adaptyv-2026/rounds/5/submit>

## Submission log

Fill in as batches go out. Nothing here yet.

| batch | date (Taipei) | upload time | track | designs | cumulative T1 / T2 / T3 | accepted? | notes |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 0 | — | — | — | — | 0 / 0 / 0 | — | not submitted |
