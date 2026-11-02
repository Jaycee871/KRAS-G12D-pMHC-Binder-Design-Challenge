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
| decoy pMHC | `AVFDRKSDAK` | **do not bind** |

The G12D and WT complexes differ by **one side chain**: peptide position 6,
aspartate vs glycine. Everything else — the entire HLA-A\*11:01 α1/α2 groove
rim, B2M, the peptide backbone, the other nine side chains — is identical.

**Therefore: a binder whose interface sits on the groove rim or the helices
will bind all four target states and score zero.** The paratope has to be centred on
the solvent-exposed face of the peptide, directly over position 6, and it has
to read the Asp carboxylate specifically — a hydrogen-bond donor / salt-bridge
partner (Arg, Lys, His, Ser/Thr/Tyr pairs) buried against it, so that the
glycine version loses contact *and* pays a desolvation penalty rather than
simply losing one weak contact. Peptide-contact area should dominate the
interface; HLA contact should be the minimum needed to position the binder
(roughly: aim for ≥ 40 % of buried surface on the peptide chain).

The apo (peptide-free) requirement follows from the same design rule — a binder
that needs the peptide's backbone and position-6 side chain to complete its
interface has nothing to grip on empty HLA. The decoy row is the fourth state
added after the literature review (sources and rationale in section 1): an unrelated, well
characterised HLA-A\*11:01 10-mer. It catches the failure the first three states
miss — a binder that reads the groove rim plus generic peptide backbone and
clears the WT gate only by luck.

Two cautions from the literature review ([`LITERATURE.md`](LITERATURE.md)) bear
directly on how the interface is judged:

* **Contact counting alone can mislead.** In Poole et al. 2022 — the one paper
  on exactly this pair — the selective TCR shows *few discernible differences in
  direct contacts* between the WT and G12D complexes; the selectivity is
  attributed to **indirect electrostatic** effects of the Asp side chain. So a
  design that looks non-selective by contact count may still be selective
  energetically, and vice versa. Treat the contact count as a necessary filter,
  not as evidence of specificity, and let the buried-Asp fraction (section 6)
  and the four-way ipTM comparison carry the judgement.
* **Excess negative interface charge hurts.** Xie et al. 2025 found that the
  refinement round which raised their experimental success rate also **reduced
  excessive negative charge** at the interface. That matters here more than
  usual, because the residue we are discriminating on *is itself an aspartate*:
  the temptation is to decorate the paratope with carboxylates, which both
  repels the target and costs solubility. Present complementary **positive
  charge locally** — Arg/Lys/His, or Ser/Thr/Tyr donors — directly over position
  6, while keeping the binder's **net interface charge** from running strongly
  negative. Both are checked as a gate in section 6.

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

The counter-screen decoy is kept alongside them as
`decoy_a1101_complex.fasta` — the same HLA-A\*11:01 heavy chain and B2M with
the peptide replaced by `AVFDRKSDAK`, the EBV EBNA-3B (EBNA-4) 399–408 "AVF"
decamer, one of the two immunodominant HLA-A\*11:01 epitopes in A11-positive
donors and a standard A\*11:01 reference binder
([de Campos-Lima et al. 1994, *J Exp Med* 179:1297,
10.1084/jem.179.4.1297](https://doi.org/10.1084/jem.179.4.1297);
[Midgley et al. 2003, *J Virol* 77:11517,
10.1128/JVI.77.21.11517-11530.2003](https://doi.org/10.1128/JVI.77.21.11517-11530.2003)).
It was chosen because it is a 10-mer — same length and therefore the same
groove register and bulge geometry as `VVVGADGVGK` — is unrelated in sequence,
and carries its own acidic side chains (Asp4, Asp8), so it also catches a
paratope that has merely learned to like a carboxylate on the peptide face
rather than *this* carboxylate at *this* position. There is no PDB entry for it
on A\*11:01, so the complex is predicted rather than templated; where a
template is wanted instead, `QVPLRPMTYK` (HIV-1 Nef 73–82, PDB 1QVO) is a
crystallographically characterised A\*11:01 decamer and is a drop-in
alternative.

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
  at peptide position 6 — complementary donors for the Asp carboxylate, and
  deliberately *not* more carboxylates (Xie et al. 2025; section 0).
* Filters: net charge −6…+6, net *interface* charge −2…+4 (section 6), no
  N-glycosylation sequon (N-X-S/T), no more than 4 consecutive identical
  residues, no Met at position 1 issues for expression.

## 4. Monomer refolding filter — AlphaFold2 (single sequence / monomer)

Refold each designed sequence on its own and keep only:

* **pLDDT ≥ 85** (mean over the binder), and
* **RMSD < 2 Å** (Cα) to the design backbone it came from.

This removes the bulk of ProteinMPNN output cheaply and is the single most
predictive cheap filter for expression and folding. Expect ~5–15 % survival.

## 5. Complex prediction — run four times per design

For every surviving design, run complex prediction **four times**, against:

1. the **G12D** complex (HLA + B2M + `VVVGADGVGK`),
2. the **WT** complex (HLA + B2M + `VVVGAGGVGK`),
3. **peptide-free** HLA + B2M,
4. the **decoy** complex (HLA + B2M + `AVFDRKSDAK`).

Use AlphaFold-Multimer and/or Boltz; where compute allows, run both and require
agreement. Record per run:

* `ipTM` (and pTM), the primary interface confidence
* **PAE-interface**: mean predicted aligned error across the binder↔target
  residue pairs (lower is better; < 10 Å is the usual cut)
* buried surface area, split into peptide-buried vs HLA-buried
* number of binder heavy atoms within 5 Å of peptide position 6
* **buried fraction of the mutated Asp at peptide position 6** (definition in
  section 6)
* net charge of the binder residues forming the interface (section 6)

The fourth run is one extra prediction per design — the cheapest available
insurance against the dominant failure mode, an allele binder that passes a
three-state panel. Run it on every design that survives section 4; there is no
point ranking a design before it is known not to be a generic A\*11:01 binder.

## 6. Scoring and ranking — the specificity gap

Primary ranking quantity, mirroring the competition's own ordering
(mutation specificity → peptide-dependent binding → affinity):

```
specificity_gap = score(G12D) - max( score(WT), score(apo), score(decoy) )
```

where `score` is ipTM (or a normalised ipTM / PAE-interface composite). Rank
descending on `specificity_gap`; this column is what `src/make_submission.py`
sorts by out of the box.

### Buried fraction of the mutated Asp (peptide position 6)

A second design objective, alongside the specificity gap. Define

```
asp6_buried_frac = ( SASA_free - SASA_complex ) / SASA_free
```

where `SASA_free` is the solvent-accessible surface area of the **Asp6 side
chain** (peptide chain C residue 6, side-chain atoms CB/CG/OD1/OD2) in the
**free pMHC** — HLA + B2M + peptide, binder absent — and `SASA_complex` is the
same atoms' SASA in the **binder-bound complex**. Compute both with the same
probe radius (1.4 Å) and the same algorithm (FreeSASA / Shrake–Rupley) on the
predicted G12D complex, so the two numbers differ only by the binder.

* **Gate: `asp6_buried_frac ≥ 0.60`**, i.e. the binder must bury at least 60 %
  of the aspartate side chain it is supposed to be reading.
* **Preferred: ≥ 0.75**, used as a ranking term behind the specificity gap.

The reasoning is energetic, not geometric. A binder that genuinely *buries* the
carboxylate makes the glycine variant pay a real desolvation and cavity penalty;
a binder that merely touches position 6 loses one weak contact and keeps
binding. Where compute allows, confirm with Xie et al.'s reverse-mutation score
— rescore the top designs with the D→G substitution in place and keep those with
the largest predicted ΔΔG.

Hard gates applied before ranking:

* `ipTM(G12D) ≥ 0.75`
* `max( ipTM(WT), ipTM(apo), ipTM(decoy) ) ≤ 0.50`
* `specificity_gap ≥ 0.25`
* ≥ 40 % of buried interface area on the peptide chain
* `asp6_buried_frac ≥ 0.60`
* ≥ 3 binder heavy atoms within 5 Å of peptide position 6 (necessary, not
  sufficient — see the Poole caution in section 0)
* net interface charge in **−2 … +4**, with at least one Arg/Lys/His or
  Ser/Thr/Tyr donor within 4 Å of the Asp6 carboxylate
* PAE-interface(G12D) < 10 Å

These thresholds are mirrored in code in [`src/design_gates.py`](../src/design_gates.py)
so the document and the tooling cannot drift apart; `evaluate()` returns the
list of gates a scored design fails.

Tie-breakers, in order: larger `specificity_gap`, then higher
`asp6_buried_frac`, then lower PAE-interface on G12D, then higher monomer
pLDDT, then topological diversity (do not submit 20 near-identical sequences —
cluster at 70 % identity and take cluster representatives).

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
