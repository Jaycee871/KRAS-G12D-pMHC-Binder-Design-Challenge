# Batch-09 matched-anchor results — 2026-11-04

**Correction notice (2026-11-04):** the original readout flags below are historical.
The panel acceptor search admitted side-chain carbons and omitted the burial leg.
Use [`../readout-correction/README.md`](../readout-correction/README.md) for the
oxygen-specific two-leg remeasurement, current batch-09 counts and measured RMSDs.
Original checkpoint/analysis files and raw responses are preserved unchanged.


## Conclusion

Matching the HLA regional docking prior does **not** remove the apo over-burial
pattern in these three designs. Every candidate has greater median buried
surface on peptide-free HLA than on the G12D complex. This is a computational
warning of peptide-independent groove engagement; buried area is **not** an
affinity measurement and does not establish experimental binding or selectivity.

`bb50-g1` remains the best supported of these three existing candidates. It
retains direct Asp6 engagement in 3/3 draws and the lowest unrelated-peptide
footprint, but it still fails to resolve the apo concern. The regional-anchor
readout is 2/3 for `bb25` and 1/3 for `bb54-g2`. These new draws must not be
pooled with the earlier peptide-C6-constrained draws as if the protocol were
unchanged.

## Frozen protocol and execution

- Source commit: `75176b2053690109b3b3039c5634ac1ac3083508`.
- Actions run: https://github.com/Jaycee871/KRAS-G12D-pMHC-Binder-Design-Challenge/actions/runs/41828512409
- Completed 2026-11-04 at 20:24 Asia/Taipei.
- 12 hosted Boltz-2 requests; 36 structures, all SHA-1 hashes distinct;
  3 samples per candidate/state; no scoring errors; ~150 seconds in runner.
- HLA A69/A155 regional anchor, maximum distance 10 Å, identical in all arms.
- `bb50-g1` and `bb25`: 8 recycling / 50 sampling steps, potentials enabled.
  `bb54-g2`: 3 recycling / 50 sampling steps, without potentials.
- Only within-design state comparisons hold the sampling settings constant.
- No new backbones or sequences; no Proteinbase submission; no new hard gate.

## Measured comparison

Buried areas below are medians in Å²; ratios are ratios of the respective
medians, not affinities. Each state has n=3.

| Candidate | G12D total area | Apo total area | Apo/G12D total | G12D Asp6 readout | WT/G12D peptide area | Decoy/G12D peptide area |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `bb50-g1` | 1037.3 | 1387.0 | 1.337 | 3/3 | 0.763 | 0.371 |
| `bb54-g2` | 725.7 | 1200.6 | 1.654 | 1/3 | 0.829 | 0.825 |
| `bb25` | 694.2 | 1284.3 | 1.850 | 2/3 | 0.969 | 0.865 |

`bb50-g1` has a 23.7% lower WT peptide footprint and 62.9% lower decoy peptide
footprint relative to G12D under this setup. Those are geometric diagnostics,
not proof of WT or decoy rejection. In particular, the decoy footprint does
not vanish and a low peptide footprint alone does not exclude HLA binding.
The previously retired 30% WT geometry-loss gate stays retired.

The regional anchor encourages groove engagement in every state, including
apo. Thus the persistence of apo over-burial rules out **unequal regional
anchors as the sole explanation for this pattern in this experiment**; it
does not demonstrate that a binder prefers apo in solution. This small,
constrained panel does not provide a free-energy difference, a dissociation
constant or positive-control-calibrated apo selectivity threshold.

## Decision

Follow the predeclared branch in `SUMMARY.md`: stop attempts to rescue these
scaffolds through further unguided backbone generation. Prepare a diversity-
ranked package from the existing pool, retaining `bb50-g1` as the strongest
lead and labelling its unresolved peptide dependence explicitly. Reassess
`bb54-g2`'s rank given its 1/3 readout under the matched regional prior.
No candidate is promoted to “validated selective” on this result.

## Evidence and reproduction

- `matched_anchor_checkpoint.jsonl`: all 36 sample-level geometry/confidence records.
- `matched_anchor_run.json`: model-call count, exact settings and completion status.
- `matched_anchor_analysis.json`: complete-panel summaries and diagnostic ratios.
- `structures/`: all 36 coordinates, hash-checked against the checkpoint.
- `matched_anchor_responses.tar.gz` and `.sha256`: full model responses and
  exact request bodies, archived without authentication headers.
- Original Actions artifact: `batch09-matched-anchor-41828512409-1`,
  ID `13322180421`; ZIP SHA-256
  `885e70c1d410f3dbcff38859d6c2ce09417c48a420c7bb3a3ecccdc7fded296c`.

After restoring the coordinate files, run `python scripts/batch09_analyze.py`.
For response-level recovery, extract the response archive in `results/batch-09/`
and follow `docs/BATCH09_RUNBOOK.md`. No additional model calls are needed to
reanalyze the completed panel. All 36 coordinate-derived metric dictionaries
were also independently recomputed locally against the checkpoint.
