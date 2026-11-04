# Track-1 fallback-20 package — review before upload

**Date:** 2026-11-05 (Asia/Taipei)

This directory contains a competition-format 20-design package assembled only
after the predeclared batch-11 stop rule froze further backbone generation.

## Files

- `track1_fallback20.csv` — the three-column upload file, ranked top to bottom.
- `fallback20_audit.csv` — the same 20 sequences plus provenance, evidence tier
  and known negative evidence. This file is for researcher review, not upload.
- `METHODS.md` — concise ranking/design-method note.

## Important limitation

This is a **fallback submission package, not a validated selective-binder set**.

Ranks 1–11 are the existing retained review pool. Ranks 12–20 are explicitly
labelled `fallback-only`: they had previously been excluded from the curated
11-sequence review set for measured or family-level failure reasons. They are
included only because Track 1 asks for at least 20 designs and because the
predeclared generation stop rule forbids inventing a new post-hoc rescue batch
after batch-11.

No sequence has been altered to fill the count, no excluded candidate is
silently promoted, and the upload file preserves the exact rank in
`fallback20_audit.csv`.

## Ranking policy

Ranks 1–3 retain the corrected matched-panel lead order. Ranks 4–11 retain the
existing review order even though batch-10 subsequently found 0/8 Phase-1
advancers; preserving that order avoids post-hoc metric shopping.

The fallback-only tail is ordered by the strength of recorded direct-Asp6
evidence first, then by family/topology diversity. Known negative evidence is
never erased. In particular:

- the bb135 fallback has strong corrected direct-Asp6 readout but measured
  unrelated-peptide Asp8 promiscuity;
- rigidified bb54 variants were retired after replication failure;
- batch-08 shape/de-wetting variants were retired after the family failed to
  produce a reproducible 3/3 direct-Asp6 design.

Because exactly 20 designs are present, all 20 are within Track 1's screened
top-20 window if they pass competition filtering.

## Validation

Run:

```bash
python src/validate_submission.py results/final-submission/track1_fallback20.csv
python -m pytest -q
```

The package is intentionally not uploaded automatically. Researcher review of
the final ranked sequences is the last step before submission.
