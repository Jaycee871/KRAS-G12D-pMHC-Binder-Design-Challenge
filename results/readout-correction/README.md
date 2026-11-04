# Carboxylate readout correction — 2026-11-04

**This additive correction supersedes the archived panel readout flags for the
coordinates actually remeasured. It spends zero model calls. Historical
checkpoints, original analysis JSONs, model responses and structures are not
overwritten. Missing coordinates remain unverified.**

## Defect and fixed definition

`scripts/selectivity_panel.py:panel_metrics` formerly searched every position-6
side-chain heavy atom, including ASP CB/CG, for a qualifying binder donor. A
carbon contact could therefore pass the carboxylate readout. It also omitted
`asp6_buried_frac > 0`, the second leg of the original batch-04 gate.

The fixed implementation searches ASP OD1/OD2 (GLU OE1/OE2 for descriptive
acidic-residue diagnostics), requires a chemically matched donor residue/atom,
uses the unrounded distance for the unchanged 3.5 Å cutoff, and requires positive
Asp6 side-chain burial. Carbon contacts remain legitimate contact-census entries;
they are not acceptors. WT Gly6 and peptide-free HLA cannot pass the Asp6 gate.
The version label is `carboxylate-oxygen-and-burial-v1`.

The real counterexample is `bb50-g1_g12d_matched_r1.cif`: ARG44 NH1 lies 3.306 Å
from ASP6 CB, but 4.199 Å from OD1 and 4.519 Å from OD2. Its burial fraction is
positive (~0.577), so oxygen specificity, rather than the missing burial check,
is what changes this draw from pass to fail.

## Coverage measured in this correction

Of **468 checkpoint rows**, **169** have retained coordinates and were
remeasured. **299** lack coordinates in the clone and remain unverified.
Exactly **one archived gate flag changes** among the 169 remeasurements:
`bb50-g1`, batch-09 G12D draw 1, true → false. Restoring the missing burial
leg changes none of these 169 flags by itself; the observed flip is the carbon
acceptor defect. This does not clear the 299 unavailable records.

## Corrected matched-panel result

All nine G12D coordinate-derived flags agree with the original two-leg
`complex_analysis.interface_metrics` + `run_batch04.gate_donor_distance` path.

| Candidate | Archived flag count | Corrected full gate | G12D binder Cα RMSD mean / max (Å) |
| --- | ---: | ---: | ---: |
| bb50-g1 | 3/3 | **2/3** | 7.32 / 9.97 |
| bb54-g2 | 1/3 | **1/3** | 5.34 / 6.59 |
| bb25 | 2/3 | **2/3** | 15.96 / 20.53 |

RMSDs are all three pairwise comparisons among three draws, fitted to HLA A +
peptide C; the binder is measured without refitting it. Apo uses HLA A only.
The corrected analysis reports all four states, exact pairs, fitting chains,
frame RMSDs and coordinate coverage. These are within-arm pose spreads, not
accuracy against an experimental binder structure. Earlier selectivity-panel
RMSDs (such as bb50-g1 1.58/3.78 Å) describe another constraint protocol and do
not characterize these matched-anchor poses.

Total/peptide buried-area ratios are unchanged: the corrected code uses the
same area primitives and only repairs the readout chemistry. Apo over-burial
remains a structural warning, not an affinity measurement. No candidate gains
validated selectivity from this correction.

## Scope and limitations

The shared function is called directly for every batch-09 arm. In the earlier
selectivity panel and batch-06/07/08 panels, G12D scoring generally used the
original batch-04 path while off-target scoring used `panel_metrics`. Thus a
shared function does not imply that all historical G12D pass counts used the
bug. `audit.json` records actual remeasurement coverage and changed flags rather
than assuming that every historical row is affected.

Rows are grouped by checkpoint source, design, biological state and condition.
Coordinates are checked against stored SHA-1 where present; a newly computed
SHA-256 is recorded for every remeasurement. A missing recorded hash does not
become proof of original-response identity. Missing coordinates and their
recorded hashes are listed separately. No failure, absence or unchanged value
is inferred for an unavailable sample.

## Files and reproduction

- `audit.json`: input/code hashes, per-row corrected metrics, original flags,
  changed flags, per-group coverage and unavailable coordinates.
- `batch09_corrected_analysis.json`: corrected full-gate counts and newly
  measured pose spreads, with a hash of the historical analysis it supersedes.
- `analysis_environment.json`: versions used for this offline reanalysis;
  hosted model builds remain unknown.
- `submission-review/`: separately regenerated review inventory using the audit
  overlay and corrected matched analysis. It is still not a final submission.

```bash
python scripts/audit_carboxylate_readout.py
python scripts/build_evidence_pool.py \
  --readout-audit results/readout-correction/audit.json \
  --matched-analysis results/readout-correction/batch09_corrected_analysis.json
python src/validate_submission.py results/readout-correction/submission-review/format_preview_only.csv
python -m pytest -q
```

The candidate builder checks the audit's source hashes and keeps remeasured
flags separate from unverified archived flags. Its historical mode remains
available for reproducing the original inventory. The unchanged heuristic lead
order is not evidence that the old 3/3 claim remains valid.

`batch09_analyze.py` now rejects unversioned legacy readout rows by default.
To reproduce the old recorded summary explicitly, run
`python scripts/batch09_analyze.py --historical-recorded-metrics`; it writes
`historical_recorded_analysis_rebuilt.json`, not the original analysis file.

## Regression coverage

Tests include a carbon-near/oxygen-far counterexample; each ASP/GLU carboxylate
oxygen; a 3.504 Å distance that rounds to 3.50 but must fail the 3.5 Å cutoff;
positive burial required even with a valid oxygen contact; and all nine real
G12D samples compared against the original gate implementation. Audit tests
check missing structures, source preservation and refusal of changed coordinate
hashes. Version checks prevent a corrected analyzer silently consuming legacy
flags. Chemistry correction does not validate hydrogen-bond angles, protonation
states or experimental affinity.
