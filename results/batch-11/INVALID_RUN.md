# Invalid infrastructure run — 2026-11-05

The first attempted execution of batch-11 did **not** produce a scientific
result.

All 24 RFdiffusion requests were rejected before a model backbone was returned:

```
NvidiaAuthError: POST /v1/biology/ipd/rfdiffusion/generate rejected the API key
(HTTP 401): Invalid JWT serialization: Missing dot delimiter(s).
```

Therefore:

- RFdiffusion model outputs: **0**
- ProteinMPNN calls on valid backbones: **0**
- candidate sequences: **0**
- G12D/apo Boltz-2 calls: **0**
- WT/decoy Boltz-2 calls: **0**
- Track-1 submission slots: **0**

The earlier generated `design_run.json` incorrectly labelled this as
`completed` because per-backbone API exceptions were caught and treated like
ordinary rejected backbones. That bookkeeping bug has been fixed: batch-11 now
fails closed on generation-service errors and uses its first planned
RFdiffusion request as a fail-fast infrastructure preflight.

**Do not cite the first run as "0/24 peptide-only backbones passed".**
The peptide-only hypothesis remains untested until RFdiffusion authentication
succeeds and actual model backbones are returned.

The original error-bearing `backbone_poses.json` is retained as audit
evidence rather than deleted or rewritten.
