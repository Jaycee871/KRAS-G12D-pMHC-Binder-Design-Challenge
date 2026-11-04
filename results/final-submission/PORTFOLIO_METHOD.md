# Track-1 portfolio ranking v2

This ranking is a **submission-allocation policy**, not a replacement for the
scientific evidence hierarchy in `results/submission-review/`.

The previous fallback-20 file preserved the review order, but that produced a
poor wet-lab portfolio: six of twenty rows came from the bb54 family, two
sequences carried four cysteines, and eight rows in the middle had already
failed the batch-10 matched G12D/apo advance rule.

The v2 portfolio keeps those scientific results unchanged and asks a different
question:

> Given only the already-generated sequence inventory, how should twenty
> reserved Track-1 wet-lab slots be diversified under severe model uncertainty?

## Frozen submission-level policy

### Ranks 1–3: exploitation

Keep the three corrected matched-panel leads:

1. bb50-g1
2. bb25
3. bb54-g2

They are not claimed to be selective. All three retain their apo-HLA warnings.
They remain at the front because they contain the strongest direct-Asp6 evidence
in the frozen inventory.

### Ranks 4–20: exploration

Use one sequence from each of seventeen additional backbone families.

Eligibility is deliberately conservative at the sequence/submission level:

- exclusion reason must be only `no recorded successful direct Asp6 readout`;
- 40–100 aa;
- at most two cysteines;
- must pass the repository's frozen sequence-liability filter;
- must not reuse a family represented by the batch-10 failed lower-rank panel;
- one sequence per family in the top 20.

Among eligible candidates, selection is greedy and deterministic:

1. maximize the minimum normalized edit distance to sequences already selected;
2. prefer a `donor-anchored` candidate;
3. prefer zero cysteines, then fewer cysteines;
4. break ties lexicographically.

This deliberately chooses **unknown but decorrelated** candidates over exact
sequences with stronger evidence of the wrong failure mode.

### Ranks 21–30: reserve

Ten additional sequences are submitted as reserve rows. Track 1 screens the
highest-ranked designs that pass competition filtering, so these rows provide
backups if a higher-ranked design is removed by novelty or other hard filters.

Reserve rows may include known batch-10 failures or retired families, but:

- every sequence still passes the repository sequence-liability filter;
- every sequence has at most two cysteines;
- no family appears more than twice in the full 30-row file.

Direct-readout evidence is prioritized within the reserve tail.

## What this policy does not claim

- Sequence diversity is not proof of structural/topological independence.
  Backbone-family labels are used as the best available provenance proxy.
- The exploration rows do not have positive direct-Asp6 evidence.
- The top three still carry unresolved apo-HLA risk.
- This policy does not change any batch gate or reinterpret a failed model run
  as a success.

The purpose is narrower: reduce correlated experimental failure while keeping
the evidence record intact.

## Reproduction

```bash
python scripts/build_portfolio_submission.py
python src/validate_submission.py results/final-submission/track1_portfolio30.csv
python -m pytest -q
```

The authoritative upload candidate is:

`results/final-submission/track1_portfolio30.csv`

The old `track1_fallback20.csv` is retained only as an audit artifact and
should not be uploaded if portfolio v2 is accepted.
