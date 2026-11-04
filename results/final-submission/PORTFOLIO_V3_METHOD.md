# Track-1 portfolio ranking v3 — balanced evidence × diversity

Portfolio v3 keeps the scientific evidence table frozen and changes only how
Track-1 wet-lab slots are allocated.

## Why v3 exists

Portfolio v2 solved three real problems in the original fallback-20 package:
family concentration, high cysteine burden and known batch-10 failures occupying
many top-20 slots. But v2 over-corrected: ranks 4–20 became seventeen mostly
under-tested exploration sequences with no positive oxygen-specific Asp6
readout. Only three of the top twenty retained direct-readout evidence.

v3 therefore uses three buckets.

### 1–3 — exploitation leads

The three corrected matched-panel leads remain first:

1. bb50-g1
2. bb25
3. bb54-g2

Their unresolved apo-HLA warnings remain authoritative.

### 4–8 — evidence-bearing, high-risk hedges

Five additional families receive one slot each because they contain some
positive direct-readout evidence despite important negative evidence:

- **bb135-g1** — oxygen-specific G12D Asp6 7/7, the strongest single direct
  readout in the inventory. It also shows unrelated-decoy Asp8 donor contact in
  3/4 draws and strong apo-HLA burial in the earlier asymmetric panel. It has
  four cysteines. Two independent G12D predictions place C20-C53 and C27-C42
  as recurring pair-like SG-SG contacts, so four cysteines are treated as a
  developability risk flag rather than an automatic rejection. No disulfide was
  explicitly constrained.
- **bb36-g2** — oxygen-specific G12D 2/6 in batch-06; batch-10 matched G12D
  1/3 and apo/G12D median total-area ratio 1.258. This is the single deliberate
  batch-10 failed sequence retained in the top 20 as an evidence-risk hedge.
- **B8-apex-54_bb06-g2**, **B8-apex-46_bb08-g1-dw** and
  **B8-span-50_bb28-g1-dw** — each recorded one corrected direct-Asp6 success
  in three G12D draws. Their batch-08 families were retired after poor
  replication and matched apo counter-screens were not run.

These five are not promoted to "validated" status. Their negative evidence is
the reason they are hedges rather than leads.

### 9–20 — family-diverse exploration

Twelve additional families are selected from the frozen candidate pool using
the same deterministic max-min sequence-diversity rule as v2.

Eligibility:

- exclusion reason is only `no recorded successful direct Asp6 readout`;
- 40–100 aa;
- at most two cysteines;
- passes the repository sequence-liability filter;
- no family duplication with ranks 1–8;
- no additional batch-10 failed family.

These rows are weak-negative / under-tested candidates, not positive hits.

### 21–30 — reserve

Ten reserve sequences provide backups for competition hard filters. Historical
exclusion reasons remain visible and no family may occur more than twice in the
full 30-row file.

## Frozen portfolio properties

- top 20: 20 distinct family labels;
- top 20: 8 designs with some direct-readout evidence;
- top 20 allocation: 3 leads + 5 evidence-risk hedges + 12 exploration;
- bb135 is the only >2-Cys exception;
- no batch-10 exact failure returns to the top 20 except bb36, which is
  explicitly predeclared as the single evidence-risk hedge;
- no biological threshold or historical batch result is changed.

## Files

Candidate upload:
`track1_portfolio30_v3.csv`

Audit:
`portfolio30_v3_audit.csv`

Summary:
`portfolio30_v3_summary.json`

Rebuild:

```bash
python scripts/build_portfolio_v3.py
python src/validate_submission.py results/final-submission/track1_portfolio30_v3.csv
python -m pytest -q
```
