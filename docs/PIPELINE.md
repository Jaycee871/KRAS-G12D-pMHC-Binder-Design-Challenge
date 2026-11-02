# Design pipeline

Planned pipeline for Challenge 5. **Nothing in this document has been run yet** —
the GPU stages below are not implemented in this repository; only target
preparation, submission building and validation are (see `src/`).

## 0. The central difficulty

The three objectives are not three separate problems, they are one problem:

| state | peptide | required outcome |
| --- | --- | --- |
| G12D pMHC | `VVVGADGVGK` | **bind** |
| WT pMHC | `VVVGAGGVGK` | **do not bind** |
| peptide-free HLA + B2M | — | **do not bind** |

The G12D and WT complexes differ by **one side chain**: peptide position 6,
aspartate vs glycine. Everything else — the entire HLA-A\*11:01 α1/α2 groove
rim, B2M, the peptide backbone, the other nine side chains — is identical.

**Therefore: a binder whose interface sits on the groove rim or the helices
will bind all three targets and score zero.** The paratope has to be centred on
the solvent-exposed face of the peptide, directly over position 6, and it has
to read the Asp carboxylate specifically — a hydrogen-bond donor / salt-bridge
partner (Arg, Lys, His, Ser/Thr/Tyr pairs) buried against it, so that the
glycine version loses contact *and* pays a desolvation penalty rather than
simply losing one weak contact. Peptide-contact area should dominate the
interface; HLA contact should be the minimum needed to position the binder
(roughly: aim for ≥ 40 % of buried surface on the peptide chain).

The apo (peptide-free) requirement follows from the same design rule — a binder
that needs the peptide's backbone and position-6 side chain to complete its
interface has nothing to grip on empty HLA.

## 1. Target preparation (implemented — `src/fetch_target.py`)

Downloads RCSB **7OW6**, strips the TCR chains, waters and heteroatoms, and
writes three target states to `data/targets/`:

* `g12d_complex.pdb` / `.fasta` — HLA-A\*11:01 heavy chain (chain A, 276 aa) +
  B2M (chain B, 99 aa) + peptide `VVVGADGVGK` (chain C, 10 aa)
* `wt_complex_template.pdb` + `wt_complex.fasta` — same complex with the
  peptide sequence changed at position 6 (D→G). The PDB is a *template*
  carrying G12D coordinates; the WT side chain must be rebuilt before scoring,
  and the FASTA is the authoritative input for complex prediction.
* `apo_hla_b2m.pdb` / `.fasta` — peptide-free HLA + B2M

Hotspot definition for the design stage: peptide chain C residues **4–8**,
weighted on residue 6 (Asp), plus the small number of α1/α2 rim residues within
8 Å of that patch. Define the epitope from the peptide outward, never from the
helices inward.

## 2. Backbone generation — RFdiffusion / BindCraft

* RFdiffusion binder design (`scaffoldguided` or free hallucination) against
  the `g12d_complex` target with the hotspot residues above.
* Length sweep 55–120 aa — small enough to sit over the peptide without
  draping onto the helices, large enough to fold. Keep to the 10–250 aa
  submission window regardless.
* 3,000–6,000 backbones per length band; discard any whose predicted contact
  map does not place ≥ 3 residues within 5 Å of peptide position 6.
* BindCraft run in parallel on the same hotspot as an independent prior — its
  AF2-in-the-loop hallucination tends to produce different topologies from
  RFdiffusion, which matters for the novelty requirement.
* Topology preference: compact mini-binders (three-helix bundle, ferredoxin-like
  αβ) presenting a concave or loop-rich face. For the antibody-format tracks,
  nanobody scaffolds with an extended CDR3 that can insert into the groove.

## 3. Sequence design — ProteinMPNN

* **8–16 sequences per backbone**, sampling temperature 0.1–0.2.
* Fixed-chain mode with the target chains held fixed; design only the binder.
* Bias against free cysteine (no unpaired Cys) and against long hydrophobic
  surface patches; soft bias toward Arg/Lys/His/Tyr at positions whose Cβ points
  at peptide position 6.
* Filters: net charge −6…+6, no N-glycosylation sequon (N-X-S/T), no more than
  4 consecutive identical residues, no Met at position 1 issues for expression.

## 4. Monomer refolding filter — AlphaFold2 (single sequence / monomer)

Refold each designed sequence on its own and keep only:

* **pLDDT ≥ 85** (mean over the binder), and
* **RMSD < 2 Å** (Cα) to the design backbone it came from.

This removes the bulk of ProteinMPNN output cheaply and is the single most
predictive cheap filter for expression and folding. Expect ~5–15 % survival.

## 5. Complex prediction — run three times per design

For every surviving design, run complex prediction **three times**, against:

1. the **G12D** complex (HLA + B2M + `VVVGADGVGK`),
2. the **WT** complex (HLA + B2M + `VVVGAGGVGK`),
3. **peptide-free** HLA + B2M.

Use AlphaFold-Multimer and/or Boltz; where compute allows, run both and require
agreement. Record per run:

* `ipTM` (and pTM), the primary interface confidence
* **PAE-interface**: mean predicted aligned error across the binder↔target
  residue pairs (lower is better; < 10 Å is the usual cut)
* buried surface area, split into peptide-buried vs HLA-buried
* number of binder heavy atoms within 5 Å of peptide position 6

## 6. Scoring and ranking — the specificity gap

Primary ranking quantity, mirroring the competition's own ordering
(mutation specificity → peptide-dependent binding → affinity):

```
specificity_gap = score(G12D) - max( score(WT), score(apo) )
```

where `score` is ipTM (or a normalised ipTM / PAE-interface composite). Rank
descending on `specificity_gap`; this column is what `src/make_submission.py`
sorts by out of the box.

Hard gates applied before ranking:

* `ipTM(G12D) ≥ 0.75`
* `ipTM(WT) ≤ 0.50` **and** `ipTM(apo) ≤ 0.50`
* `specificity_gap ≥ 0.25`
* ≥ 40 % of buried interface area on the peptide chain
* ≥ 3 binder heavy atoms within 5 Å of peptide position 6
* PAE-interface(G12D) < 10 Å

Tie-breakers, in order: larger `specificity_gap`, then lower PAE-interface on
G12D, then higher monomer pLDDT, then topological diversity (do not submit 20
near-identical sequences — cluster at 70 % identity and take cluster
representatives).

## 7. Novelty

Submissions need a novelty score **≥ 3/4**, checked automatically. Designs must
be de novo. Before each batch: BLAST / foldseek every candidate against the PDB
and against known TCR/antibody germlines, and drop anything that looks like a
retrieved natural binder or a lightly mutated known TCR from the 7OW6 entry
itself. De novo hallucinated scaffolds normally clear this easily; grafted CDR
loops from a natural TCR do not.

## 8. Submission

`src/make_submission.py` takes the scored table, ranks it, truncates to ≤ 20,
and validates the result against every documented rule before writing. Then
upload at the submit page and log the batch in `docs/SUBMISSION_PLAN.md`.
