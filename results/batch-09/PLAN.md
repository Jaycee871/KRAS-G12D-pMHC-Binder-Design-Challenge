# Batch-09 — matched HLA-groove counter-screen

**Status: experiment defined and pre-registered; no batch-09 model result is
claimed in this directory yet.**

## Question

Does the severe peptide-free HLA signal observed in earlier batches survive
when G12D, WT, apo and unrelated-peptide controls are all docked with the same
physical prior?

Earlier panels were biologically sensible but computationally asymmetric:

- G12D / WT / decoy: pocket-constrained to peptide position 6;
- apo: unconstrained, because peptide chain C does not exist.

Therefore an apo-vs-G12D difference mixes two changes: peptide removal **and**
removal of the docking constraint.

## Matched prior

The Challenge 5 description of 7OW6 states that KRAS peptide position 6 lies
against HLA-A*11:01 **Gln70 and Arg114**.  The exact target sequence stored in
this repository contains Q at HLA position 70 and R at 114 in every arm.

Batch-09 therefore uses one constraint for every arm:

```json
[
  {
    "binder": "D",
    "contacts": [
      {"id": "A", "residue_index": 70},
      {"id": "A", "residue_index": 114}
    ],
    "max_distance": 6.0
  }
]
```

The 6 Å distance is inherited unchanged from the already-used peptide-C6
pocket.  This experiment changes the invariant anchor, not the allowed
distance.

## Sampling cell

One cell only:

- recycling steps: 3
- sampling steps: 50
- without potentials: true
- diffusion samples: 3

This is the existing `nopot` / `r3np` cell.  There is no parameter sweep in
batch-09.

## Candidate panel fixed before calls

| key | reason |
| --- | --- |
| `bb54-g2` | current parent; strongest accumulated WT-discrimination evidence and reproducible Asp6 read-out under r3np |
| `bb135-g1` | tight Arg–Asp read-out; stabilised under no-potentials sampling |
| `bb173-g1` | high Asp6 burial but an intermittent alternate pose |
| `bb90-g2` | tightest historical pose spread among the batch-05 gate screeners |
| `b8-apex46-g2dw` | peptide-centred batch-08 control already measured in the old unconstrained apo diagnostic |

Each design is evaluated against:

1. G12D pMHC;
2. WT pMHC;
3. peptide-free HLA + B2M;
4. unrelated A*11:01 peptide decoy.

Five designs × four states = **20 Boltz-2 calls**, each returning three
diffusion samples.

## Decision discipline

This diagnostic adds **no new pass threshold**.

- The existing 30% WT peptide-area + 30% position-6-contact loss criterion is
  reused unchanged where applicable.
- Apo and decoy total/HLA buried-area ratios are reported descriptively.
- ipTM remains a pose/confidence read-out, not a selectivity gap.
- A result may falsify the old apo interpretation; it may not retroactively
  alter batch-06/07/08 outcomes.
- No new RFdiffusion or ProteinMPNN generation is allowed before this diagnostic
  is interpreted.

## Interpretation

Three useful outcomes are possible:

**A. Apo remains strongly bound under the matched groove prior.**  
Then the failure is real: the current sequence family is fundamentally
HLA-groove-seeking.  Stop local optimisation around those families and either
freeze a low-confidence submission pool or introduce a genuinely new
peptide-dependent topology.

**B. Apo binding collapses under the matched prior while G12D remains stable.**  
Then earlier apo measurements were substantially confounded by the unequal
docking prior.  A small, mechanism-specific batch-09b may be justified; the
apo gate itself is not weakened.

**C. All arms collapse or become highly unstable.**  
Then the A70/A114 prior is not informative enough for these designs.  Report
the diagnostic as inconclusive and do not use it to rescue candidates.

## Reproduce

Pre-flight only, no model calls:

```bash
python scripts/batch09_matched_panel.py --dry-run
python -m pytest -q tests/test_batch09_matched_panel.py
```

Live run when the NVIDIA credential is available in the execution environment:

```bash
python scripts/batch09_matched_panel.py
```

Outputs are checkpointed per design / arm / draw to
`results/batch-09/matched_checkpoint.jsonl`, with structures under
`results/batch-09/structures/`.  The aggregate diagnostic is written to
`matched_summary.json`.
