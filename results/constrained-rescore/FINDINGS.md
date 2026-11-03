# Pocket-constrained Boltz-2 rescore of batch-03 (2026-11-03)

## Why
Batch-03 finished with 0/28 designs clearing every hard gate; 20/28 were thrown
out because Boltz-2, re-docking from sequence, put the binder on the HLA face
instead of on the peptide. The agreed plan-B was to stop asking the co-folder to
rediscover the designed pose.

## The unblock
The hosted NVIDIA Boltz-2 NIM (`/v1/biology/mit/boltz2/predict`) accepts a
pocket constraint. The schema is *flat* (not nested under "pocket") and the
contact items must be objects keyed `id` / `residue_index`:

```json
"constraints": [{"binder": "D",
                 "contacts": [{"id": "C", "residue_index": 6}],
                 "max_distance": 6.0}]
```

Nested `{"pocket": {...}}`, `[chain, index]` pairs and `chain_id` keys are all
rejected (HTTP 422).

## Result: all 28 batch-03 sequences rescored under the constraint
Same constraint applied to g12d, wt and decoy so the specificity gap is a like
for like comparison; apo has no peptide and was run unconstrained.

| quantity | batch-03 (unconstrained) | constrained |
| --- | --- | --- |
| binder lands on the peptide | 8/28 | 26/28 |
| Asp6 H-bond / salt bridge | 0/28 | 3/28 |
| ipTM(G12D) >= 0.75 | 0/28 | 2/28 |
| specificity gap > 0 | 1/28 | 6/28 |
| specificity gap >= 0.25 | 0/28 | 1/28 |
| peptide-buried fraction >= 0.40 | 0/28 | 0/28 (max 0.32) |

Top rows (gap, ipTM G12D, peptide buried frac, pos6 contacts, Asp6 H-bond):
* `S1-pos6-68_bb18`  +0.298  0.708  0.068  4   no
* `S1-pos6-68_bb19`  +0.221  0.680  0.182  5   no
* `S1-pos6-68_bb13`  +0.183  0.758  0.223  12  yes
* `S1-pos6-65_bb03`  +0.045  0.754  0.190  8   yes

## Open items
1. The `peptide_buried_frac >= 0.40` hard gate is still unmet by everything, but
   the 7OW6 positive control (the selective Poole et al. TCR) scores 0.21 by the
   same measure, i.e. the threshold is stricter than a real selective binder.
   Decide whether to recalibrate it to the control.
2. ~~The constrained protocol has not yet been run on 7OW6 itself.~~ **Done,
   2026-11-03, and it is a split result** - see
   `results/validation/CONTROL_7OW6.md`. The pose is reproduced to 1.3-1.5 A
   C-alpha RMSD over the TCR variable domains (crystal-matching buried areas),
   but the specificity gap of the known selective TCR is **-0.006**
   (ipTM G12D 0.9268 +/- 0.0010 vs WT 0.9325 +/- 0.0001), and re-running the
   identical constrained panel on the three best designs below flips
   `bb18` from +0.298 to -0.300 and `bb19` from +0.221 to -0.041. The gaps in
   the table above are therefore **not trustworthy**; only `bb13` reproduced
   (+0.183 -> +0.151). Batch-04 was not launched because of this.
3. Raw data: `results/constrained-rescore/checkpoint.jsonl` and the per-design
   constrained G12D mmCIFs in the same folder. Driver: `/workspace/constrained_rescore.py`.
