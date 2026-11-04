# Batch-10 — matched G12D / apo counter-screen of review ranks 4–11

**Date:** 2026-11-05 (Asia/Taipei)  
**Status:** complete.  
**Compute:** 16 Boltz-2 calls / 48 structures.  
**Track-1 slots spent:** 0.  
**Historical gates changed:** no.

## Question

Ranks 4–11 of the corrected submission-review pool had some prior direct-Asp6
evidence but lacked a like-for-like peptide-free HLA counter-screen. Batch-10
therefore tested whether any of those eight sequences was less apo-prone than
the three headline leads before spending calls on WT and unrelated-peptide
confirmation.

Every G12D and apo arm used the same HLA regional anchor as batch-09
(HLA-A*11:01 A69 + A155 within 10 Å). Each design kept its frozen,
pre-selected sampling regime across both arms.

The predeclared Phase-1 advance heuristic was:

1. corrected oxygen-specific Asp6 readout on G12D in at least 2 of 3 draws; and
2. median apo total buried area / median G12D total buried area < 1.0.

Buried area is used here only as a directional triage diagnostic, not as a
calibrated affinity predictor.

## Result

**0 of 8 designs advance. Phase 2 is therefore 0 calls.**

| review rank | candidate | G12D corrected Asp6 | G12D median area | apo median area | apo/G12D | advance |
| ---: | --- | ---: | ---: | ---: | ---: | --- |
| 4 | B6-c6-46_w2_bb36_donor-anchored-g2 | 1/3 | 1018.1 | 1280.9 | 1.258 | no |
| 5 | B5-c6-62_bb173_donor-anchored-g1 | 1/3 | 797.6 | 1094.5 | 1.372 | no |
| 6 | bb54_graft-p3 | 0/3 | 665.1 | 1181.3 | 1.776 | no |
| 7 | B5-c6-62_bb90_donor-anchored-g2 | 1/3 | 761.2 | 1332.4 | 1.750 | no |
| 8 | B4-c67-65_bb26_donor-anchored | 0/3 | 939.3 | 1569.5 | 1.671 | no |
| 9 | B6-c6-46_w2_bb38_donor-anchored-g1 | 0/3 | 572.0 | 1247.0 | 2.180 | no |
| 10 | bb54_graft-p1-v3 | 0/3 | 545.9 | 1281.6 | 2.348 | no |
| 11 | bb54_graft-p2 | 0/3 | 666.4 | 1198.9 | 1.799 | no |

Two observations are especially strong.

First, the corrected direct-Asp6 chemistry did not replicate under the matched
regional-anchor prior: the best designs reached only 1/3, and five of eight were
0/3.

Second, the apo direction is uniform. Every design buries more total interface
area on peptide-free HLA than on G12D pMHC under the same regional docking
prior. The smallest ratio is 1.258 and the largest is 2.348.

This extends batch-09's result: the apo warning is not confined to bb50-g1,
bb25 or bb54-g2, and it is not explained by the old unequal docking priors.

## Decision

Do not run WT/decoy Phase 2 for these eight sequences. The adaptive protocol
predeclared zero confirm calls when Phase 1 produced no advancers.

Do not rescue these families by changing the threshold after observing the
result. In particular, another generic sampling sweep, smaller scaffold,
rigidification, shape-steering or helix de-wetting pass would repeat already
measured dead ends.

If another generation experiment is run before the deadline, it must introduce
a genuinely different mechanism whose design-stage geometry depends on the
peptide itself rather than merely on the HLA groove. A defensible example is a
predeclared bidentate peptide mechanism: direct Asp6 carboxylate recognition
plus a second independent peptide contact, with HLA-only docking explicitly
disfavoured before Boltz-2 ranking.

## Evidence

- `counter_screen_checkpoint.jsonl`: 48 scored structures.
- `counter_screen_run.json`: frozen protocol, settings and Phase-1 summary.
- `structures/`: all 48 coordinates.
- `counter_screen_responses.tar.gz`: raw cached model responses.
- `counter_screen_responses.sha256`: archive integrity hash.

The response archive SHA-256 is
`0c4ae470efc18a6f506845551c4c7855fe9d48ae0f1d692704253967a57cd69f`.
