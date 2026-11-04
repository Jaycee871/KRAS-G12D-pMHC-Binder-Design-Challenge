# Legacy direct-readout audit — 2026-11-05

This is the predeclared fallback from the submission plan: before inventing a
new mechanism or spending more model calls, assess batches 01–03 under the
current direct Asp6 readout.

**Result: zero additional direct-readout candidates.**

| Source | checkpoint rows | unique sequences | oxygen + burial passes | burial-only sequences |
| --- | ---: | ---: | ---: | ---: |
| batch-01 | 75 | 56 | **0** | 2 |
| batch-02 | 60 | 60 | **0** | 4 |
| batch-03 | 28 | 28 | **0** | 1 |
| **Total** | **163** | **144 across all three** | **0** | **7** |

The seven burial-only sequences are not rescued. They have recorded positive
Asp6 side-chain burial but no qualifying carboxylate-oxygen donor within the
3.5 Å direct-readout cutoff. Peptide contact counts alone are not substituted
for the oxygen-specific readout.

## What this does and does not establish

The audit uses the historical `g12d_metrics` already stored in each checkpoint.
Those metrics record the Asp6 carboxylate donor distance and Asp6 side-chain
burial separately. The current two-leg definition requires both:

1. donor-to-Asp6-carboxylate distance ≤ 3.5 Å; and
2. positive Asp6 side-chain burial.

No new Boltz-2, RFdiffusion or ProteinMPNN call is made. This is deliberately
labelled **recorded-unverified**, not coordinate-remeasured: it does not turn
legacy checkpoint fields into new structural evidence, and it does not infer
missing counter-screens.

The result closes the “assess older excluded batches” branch of
`docs/SUBMISSION_PLAN.md`. Batches 01–03 cannot fill the Track-1 shortfall
under the same direct-readout criterion used to build the current review pool.

## Reproduce

```bash
python scripts/audit_legacy_readout.py --out /tmp/legacy_readout_audit.json
python -m pytest -q tests/test_legacy_readout_audit.py
```

The script records SHA-256 hashes of all three checkpoint inputs in its JSON
output so later changes cannot be confused with this audit.
