# Existing candidate pool — review only, 2026-11-05

**Correction notice (2026-11-04):** the original readout flags below are historical.
The panel acceptor search admitted side-chain carbons and omitted the burial leg.
Use [`../readout-correction/README.md`](../readout-correction/README.md) for the
oxygen-specific two-leg remeasurement, current batch-09 counts and measured RMSDs.
Original checkpoint/analysis files and raw responses are preserved unchanged.


**11 retained review candidates; 9 short of the Track-1 minimum recorded in the
project plan. No final submission package and no Proteinbase upload. No new
model calls, sequences or backbones. No experimentally validated selective
candidate.** The format preview passes the CSV schema validator, which does
not enforce the Track-1 minimum count; that does not make it submission-ready.

## Scope and evidence

This inventory covers 213 distinct sequences recorded downstream in batches
04–09, the sampling study, the selectivity panel and the original batch-04/05
repeat tables. Batches 01–03 preceded the direct-readout selection and are not
included; this is not a claim to have exhausted every sequence ever generated.
Full input file SHA-256 hashes and per-record locations are in
`evidence_pool.json`. Every exported sequence is unchanged from those records.
The known TCR control is outside this generation pool. Competition novelty and
eligibility still require verification.

| Disposition | Sequences |
| --- | ---: |
| Retained for review | 11 |
| No recorded successful direct Asp6 readout | 147 |
| Retired batch-08 shape/de-wetting family | 37 |
| Retired rigidified family | 16 |
| bb135 family deprioritized after g1 decoy-Asp8 promiscuity | 2 |

The bb135 variant exclusion is a conservative **family curation decision**;
only the g1 sequence has the cited measured decoy-Asp8 result. It does not
assert that every family variant has been measured to be promiscuous.

## Review order

| Rank | Recorded candidate | Evidence and limitation |
| ---: | --- | --- |
| 1 | `bb50-g1` | Corrected matched Asp6 2/3; unresolved apo over-burial |
| 2 | `bb25` | Matched Asp6 2/3; apo and decoy concerns |
| 3 | `bb54-g2` | Matched Asp6 1/3; several older sampling successes; apo concern |
| 4 | `B6-c6-46_w2_bb36_donor-anchored-g2` | Batch-10 matched G12D 1/3; apo/G12D area 1.258 — no advance |
| 5 | `B5-c6-62_bb173_donor-anchored-g1` | Batch-10 matched G12D 1/3; apo/G12D area 1.372 — no advance |
| 6 | `bb54_graft-p3` | Batch-10 matched G12D 0/3; apo/G12D area 1.776 — no advance |
| 7 | `B5-c6-62_bb90_donor-anchored-g2` | Batch-10 matched G12D 1/3; apo/G12D area 1.750 — no advance |
| 8 | `B4-c67-65_bb26_donor-anchored` | Batch-10 matched G12D 0/3; apo/G12D area 1.671 — no advance |
| 9 | `B6-c6-46_w2_bb38_donor-anchored-g1` | Batch-10 matched G12D 0/3; apo/G12D area 2.180 — no advance |
| 10 | `bb54_graft-p1-v3` | Batch-10 matched G12D 0/3; apo/G12D area 2.348 — no advance |
| 11 | `bb54_graft-p2` | Batch-10 matched G12D 0/3; apo/G12D area 1.799 — no advance |

This is a transparent review heuristic, not a new biological gate or an
estimated probability of binding. The first three are the matched-panel leads,
ordered using the batch-09 evidence. Others are separated into multi-sample
mixed evidence (at least two distinct successful coordinates within a recorded
protocol) and limited-readout tiers. Within each tier, prefer a generation
backbone family not already selected, then the lowest maximum normalized edit
similarity to earlier sequences, with SHA-256 as a deterministic tie-break.
Generation backbone labels are provenance, not measured structural clusters;
normalized edit similarity is sequence distance, not topology or epitope
similarity. No ipTM gap, MM-GBSA score or historical WT-loss threshold ranks
this list. Lower ranks can have severe known failures and remain review entries.

Per-protocol outcomes are separated by source file, target state, constraint,
sampling settings and condition. A repeated coordinate hash contributes once
to that protocol's distinct-structure count. Rows lacking retrievable coordinates
or a recorded hash remain explicitly unhashed; they do not become verified
replicates. For example, bb50-g1's selectivity-panel G12D group has seven records
but four distinct structures. Distinct coordinates still do not prove statistical
independence. No pooled cross-protocol success rate is used for ranking.

## Files and reproduction

- `review_shortlist.csv`: all 11 ranks, unchanged sequences, aliases, family,
  sequence hashes and per-sequence warnings.
- `candidate_audit.csv`: all 213 sequence identities and curation reasons.
- `evidence_pool.json`: source hashes, record locations and all protocol counts,
  including failures, missing hashes and untested counter-screen arms.
- `format_preview_only.csv`: the same 11 in the three-column Proteinbase schema;
  it is below the Track-1 count requirement and is not a final submission.

```bash
python scripts/build_evidence_pool.py
python src/validate_submission.py results/submission-review/format_preview_only.csv
python -m pytest -q
```

## Next decision

Batch-10 has now directly counter-screened ranks 4–11 under the same matched
HLA regional anchor used in batch-09. **0/8 advance**: every apo/G12D median
total-area ratio is >1 and none reproduces corrected direct Asp6 recognition
in at least 2/3 G12D draws. Phase-2 WT/decoy confirmation is therefore zero
calls by the frozen adaptive rule; see `../batch-10/SUMMARY.md`.

Batch-11 tested a genuinely new peptide-only generation mechanism. Offline
side-chain restoration found 9/24 bidentate Asp6+Val8 peptide-contact backbones,
but 0/24 fit after HLA was restored (41–154 severe clash atoms, frozen limit
<=2). The predeclared stop rule therefore freezes generation.

A separate 20-sequence **fallback-only submission package** now lives in
`../final-submission/`. Its first 11 rows preserve this retained review pool;
rows 12–20 are explicitly documented previously excluded near-misses. They are
not silently reclassified as validated or retained candidates.
