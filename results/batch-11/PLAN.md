# Batch-11 — peptide-only generation with HLA re-insertion

**Status:** predeclared protocol; no model calls at the time this file is created.  
**Purpose:** test a genuinely different mechanism after batches 06–10 repeatedly
produced HLA-dominant / apo-prone interfaces.

## First-principles intervention

Earlier batches generated a binder while the HLA groove was present and then
tried to persuade or filter the binder into using more peptide surface. Batch-08
already tested multi-residue peptide hotspots (including C4/C6/C8 and C5/C6/C7)
plus explicit groove-shape filters, so simply adding more hotspots is not a new
mechanism.

Batch-11 removes the HLA heavy chain and B2M from the **generative target**.

RFdiffusion receives only the crystallographic KRAS G12D 10-mer in its bound
7OW6 conformation. ProteinMPNN likewise sees only that peptide plus the
generated binder backbone. HLA therefore cannot contribute a favourable design
surface or sequence signal during either generation step.

The complete HLA is re-inserted geometrically before a sequence is accepted,
and complete pMHC / apo states return for Boltz-2 screening.

## Design-stage protocol

Three binder lengths are sampled: 40, 46 and 52 aa. Every RFdiffusion request
uses peptide hotspots C5, C6 and C8. The mechanism is **not** defined by the
hotspot request; each returned backbone must pass an explicit residue-level
test:

1. Asp6 side chain: at least one binder heavy-atom contact within 5 A and
   closest side-chain approach <= 4.5 A.
2. Val8 side chain: the same requirement.
3. The two contacts must involve at least two distinct binder residues.
4. The peptide coordinates returned by RFdiffusion are aligned back to 7OW6.
   After HLA residues A1–180 are restored, at most two binder heavy atoms may
   lie within 2.2 A of any HLA heavy atom.

There is deliberately **no batch-08 peptide/helix contact-fraction, approach
height, lateral-offset or de-wetting gate**. The experimental variable is the
absence of HLA during generation, not another retuning of those thresholds.

For each surviving backbone, soluble ProteinMPNN samples are ranked by a
predeclared paired-contact sequence score: donor chemistry around Asp6 plus
hydrophobic/aromatic chemistry around Val8. One deterministic paired-anchor
variant is allowed: nearest Asp6-facing binder position -> Arg and nearest
distinct Val8-facing position -> Tyr. Standard sequence-liability filters still
apply.

## Phase 1 — unbiased G12D / apo screen

At most 12 sequence-diverse candidates enter Phase 1, preferring one sequence
per generated backbone before a second sequence from the same backbone.

Both G12D pMHC and peptide-free HLA+B2M are predicted **without a docking
constraint**, with identical Boltz-2 sampling settings:

- recycling steps: 8
- sampling steps: 50
- diffusion samples per call: 3
- potentials: on

This is at most 24 Boltz-2 calls / 72 structures.

A candidate advances only if all three predeclared conditions hold:

1. corrected oxygen-specific Asp6 readout on G12D >= 2/3 draws;
2. Val8 side-chain contact on G12D >= 2/3 draws; and
3. median apo total buried area / median G12D total buried area < 1.0.

Condition 3 is a directional triage heuristic, **not** a calibrated affinity
claim.

## Phase 2 — WT / unrelated-peptide confirmation

Only Phase-1 advancers receive unconstrained WT and decoy predictions, three
draws per arm. If Phase 1 has zero advancers, Phase 2 spends zero calls.

WT and decoy are recorded as diagnostics; no new hard selectivity threshold is
invented from them. Decoy carboxylate donor contacts remain a strong negative
signal because the unrelated peptide contains acidic side chains.

## Stop rules

- No historical threshold is changed after seeing batch-11.
- No failed batch-06/07/08 scaffold is silently rescued.
- If peptide-only generation produces no Phase-1 advancer, generation is frozen
  and the project moves to an explicitly labelled fallback submission package.
- Track-1 slots are not spent by this computation; only final submitted
  sequences consume the competition's design list.
