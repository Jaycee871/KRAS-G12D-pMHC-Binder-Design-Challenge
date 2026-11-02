# Literature review — mutation-specific pMHC binder design

**Scope.** What is published on designing protein binders that recognise a
peptide–MHC complex carrying a single point mutation, and what that implies for
Challenge 5 (a binder to KRAS G12D `VVVGADGVGK` on HLA-A\*11:01, PDB 7OW6, that
must not bind wild-type `VVVGAGGVGK` or peptide-free HLA + B2M).

**Source of the search.** Undermind literature search, query
*"Mutation specific pMHC binder design"*, run 2 Nov 2026 — 115 papers screened,
88 full texts read:
<https://app.undermind.ai/projects/362994e3-20f1-4bb8-9984-9a088a13c27f>

---

## 1. The gap: nothing published hits this exact target

**No published de novo binder exists for KRAS G12D presented on
HLA-A\*11:01.** The literature brackets the target from two sides but never
lands on it.

* **Poole et al. 2022, *Nature Communications*** — *Therapeutic high affinity T
  cell receptor targeting a KRAS G12D cancer neoantigen*
  ([10.1038/s41467-022-32811-1](https://doi.org/10.1038/s41467-022-32811-1)).
  This is the one paper covering exactly our pair: the decamer `VVVGADGVGK`
  presented by HLA-A\*11:01, discriminated against wild-type `VVVGAGGVGK`. But
  it is an **affinity-enhanced natural TCR** (the JDI TCR, affinity-matured by
  CDR randomisation), not a designed binder. Two results from it matter here.
  First, the parental TCR bound the G12D complex at K_D ≈ 63 µM with no
  measurable binding to the WT complex, and affinity maturation reached
  picomolar affinity while *widening* rather than collapsing the specificity
  window — specificity and affinity are not necessarily in tension on this
  target. Second, and more instructive: crystal structures of the TCR on WT
  versus G12D complexes show **few discernible differences in direct contacts**;
  thermodynamic analysis and MD attribute the selectivity to *indirect
  electrostatic* effects of the Asp side chain. Discrimination of this epitope
  is therefore not purely a matter of counting hydrogen bonds to the aspartate.
  Notably, their phage libraries were **depleted on the WT complex before
  positive selection** on G12D — counter-selection built into the campaign, not
  applied afterwards.

* **Ahn et al. 2025 (bioRxiv) / Ahn et al. 2026 (*Molecular Therapy*)** —
  *Discovery of TCR-like antibodies to the KRAS G12D neoantigen via in
  silico–in vitro workflow*
  ([10.1101/2025.08.31.673313](https://doi.org/10.1101/2025.08.31.673313);
  [10.1016/j.ymthe.2026.05.032](https://doi.org/10.1016/j.ymthe.2026.05.032)).
  A TCR-mimic antibody against KRAS G12D with nanomolar affinity and no
  detectable cross-reactivity to wild-type KRAS peptide — but presented on
  **HLA-C\*08:02**, a different allele with a different groove and a different
  peptide register. The workflow (in silico design followed by in vitro
  selection) transfers; the molecule does not.

So the target is open. The transferable knowledge has to come from the general
pMHC binder-design literature, which matured sharply in 2025.

## 2. The 2025 state of the art in de novo pMHC binder design

* **Liu et al. 2025, *Science*** — *Design of high-specificity binders for
  peptide–MHC-I complexes*
  ([10.1126/science.adv0185](https://doi.org/10.1126/science.adv0185)).
  The reference pipeline. RFdiffusion is targeted at the **outward-facing
  peptide residues** to generate backbones that arc *above* the groove and make
  extensive peptide contact while **avoiding extensive contact with the
  ubiquitously expressed MHC**; ProteinMPNN designs the sequence; AF2 filters
  folding and binding. Specific binders were obtained for **11 pMHC targets**,
  from experimental or predicted structures, and designs for eight targets
  conferred peptide-specific T cell activation as CARs. Two of their selection
  steps are directly reusable: they used ProteinMPNN to *predict the peptide
  sequence given the binder* and kept designs where the on-target peptide
  scored above off-target sequences, and they re-predicted each design against
  closely related proteome peptides and kept only those predicted to bind the
  on-target considerably more confidently.

* **Johansen et al. 2025, *Science*** — *De novo-designed pMHC binders
  facilitate T cell–mediated cytotoxicity toward cancer cells*
  ([10.1126/science.adv0422](https://doi.org/10.1126/science.adv0422)).
  Minibinders against NY-ESO-1 `SLLMWITQC`/HLA-A\*02:01 built with
  RFdiffusion + ProteinMPNN + AlphaFold2, with the key addition of
  **in silico cross-panning**: every surviving design was scored by AF2 against
  a deliberately constructed panel of single- and double-mutant peptide
  variants plus a proteome-derived decoy peptide on the same allele, and
  candidates were ranked by the *number of undesired predicted binding events*.
  Molecular dynamics improved the predictability of in vitro success. The
  method generalised to a second neoantigen pMHC of unknown structure
  (`RVTDESILSY`/HLA-A\*01:01).

* **Householder et al. 2025, *Science*** — *De novo design and structure of a
  peptide-centric TCR mimic binding module*
  ([10.1126/science.adv3813](https://doi.org/10.1126/science.adv3813)).
  A de novo α-helical TCR mimic for NY-ESO-1 on HLA-A\*02 with K_D = 9.5 nM and
  no detectable affinity for the same MHC presenting a different peptide. The
  2.05 Å structure shows a rigid, TCR-like docking mode with an **unusual degree
  of focus on the up-facing peptide side chains** — structural confirmation
  that the peptide-centric interface geometry is what buys specificity. A
  structure-informed in silico screen over 14,363 HLA-A\*02 peptides correctly
  predicted the two genuine off-targets, so off-target prediction from the
  designed interface is tractable.

* **Motmaen et al. 2025, bioRxiv** — *Targeting peptide–MHC complexes with
  designed T cell receptors and antibodies* (ADAPT)
  ([10.1101/2025.11.19.689381](https://doi.org/10.1101/2025.11.19.689381)).
  Extends structure-based design from minibinders to the **immunoglobulin
  formats** — designed TCRs and TCR-mimic antibodies — via fine-tuned AF2
  docking, ProteinMPNN CDR design, redocking and genetic-algorithm refinement.
  Binders were obtained for 6 of 9 TCR targets and for three antibody targets,
  with two antibody–pMHC complexes validated by **cryo-EM at 2.6 Å and 2.5 Å**,
  confirming that the designed interfaces are accurate at atomic level.

* **Xie et al. 2025, bioRxiv** — *Design of TCR-mimicking binders for pHLA with
  high potency*
  ([10.1101/2025.10.03.680404](https://doi.org/10.1101/2025.10.03.680404)).
  The **closest analogue of our task**: a single-point-mutation discrimination
  problem, p53 R175H (`HMTEVVRHC`) versus wild type (`HMTEVVRRC`) on
  HLA-A\*02:01. BindCraft generated backbones with the mutation site set as the
  design **hotspot**; MD-guided electrostatic filtering and free-energy
  perturbation — scoring the ΔΔG of reverting the mutant residue to wild type —
  selected for mutant-selective interfaces. The result was nanomolar binders
  (K_D 2.0 and 91.6 nM) with **>10× selectivity over wild type**. A later
  solMPNN-based optimisation round, designed to work on targets *without*
  electrostatic cues (motivated by KRAS G12V on HLA-A\*03:01), raised the
  experimental success rate from 2/5 to 5/5, with most candidates showing
  minimal or no wild-type binding, and reached K_D 1.6–126.1 nM. Sequence
  analysis found that refinement **reduced excessive negative charge** at the
  interface and balanced hydrophobic and polar residues around the mutated
  position.

* **Brasas et al. 2026, bioRxiv** — *Specificity-driven protein binder design
  with Odin-Multi*
  ([10.64898/2026.09.08.749745](https://doi.org/10.64898/2026.09.08.749745)).
  Makes specificity an **optimisation objective inside the design loop**: one
  sequence is jointly optimised against multiple on-target complexes while
  being actively *repelled* from specified off-target structures. The paper
  explicitly addresses single-residue pMHC discrimination. This is the
  architectural answer to the problem that filtering alone throws away most of
  the design budget.

* **Du et al. 2024, *Nature Biotechnology*** — *Targeting peptide antigens
  using a multiallelic MHC I-binding system*
  ([10.1038/s41587-024-02505-8](https://doi.org/10.1038/s41587-024-02505-8)).
  TRACeR-I: a small helical-bundle scaffold carrying a single flexible
  antigen-recognition loop that docks the peptide in a parallel orientation
  suited to the bulged MHC-I conformation, with the scaffold reused across
  alleles and only the loop re-diversified. Peptide-specific binders spanning
  roughly 22 nM to 3.7 µM were obtained across HLA-A and HLA-B targets. It
  demonstrates that **one engineered loop over the peptide**, rather than a
  large grafted interface, can carry peptide specificity — a useful scaffold
  prior for antibody/nanobody-format tracks.

## 3. Synthesis

The literature converges on one principle: **specificity must be an explicit
design objective, not a post-hoc filter.** Every successful campaign above
encoded discrimination into the generation or optimisation step — hotspot
placement on the mutated residue, negative design against off-targets,
specificity-aware scoring — and used filtering only to confirm it. Campaigns
that design for affinity and screen for specificity afterwards spend their
budget on groove-rim binders that bind everything on the allele.

Three recurring levers:

1. **Peptide-centric interfaces.** Bury the mutated residue, not the MHC
   surface. Liu, Householder and Johansen all steer the backbone over the
   outward-facing peptide side chains and deliberately minimise helix contact;
   Householder's crystal structure shows this is what the selective molecules
   actually do.
2. **Counter-design / cross-panning during design.** Score every candidate
   against the wild-type peptide *and* unrelated peptides on the same allele
   while designs are still being selected or optimised — Johansen's in silico
   cross-panning, Liu's inverse-folding peptide-identity check, Poole's
   wild-type depletion of the phage library, Odin-Multi's explicit off-target
   repulsion term.
3. **Mutation-aware sequence optimisation.** Designate the mutation as the
   hotspot and optimise the residues contacting it specifically — Xie's
   MD/FEP reverse-mutation scoring and solMPNN interface refinement, which also
   warns that **excess negative charge at the interface hurts**: relevant for
   us, because our mutated residue is itself an aspartate, so the paratope
   should present complementary donors rather than more carboxylates.

On validation, the papers agree on the panel: test **mutant pMHC vs wild-type
pMHC vs unrelated peptides on the same allele vs peptide-free MHC**. Anything
less cannot distinguish a true neoantigen binder from an allele binder. Poole's
result is the cautionary note — selectivity there arose from indirect
electrostatics with almost no change in direct contacts, so a design that looks
non-selective by contact counting may still be selective energetically, and
vice versa.

## 4. Implications for this challenge

The repository's existing design rules in [`docs/PIPELINE.md`](PIPELINE.md)
already sit on the right side of this literature. Section 0 states the
peptide-centric principle explicitly (paratope centred over peptide position 6,
≥ 40 % of buried surface on the peptide chain, hotspots defined "from the
peptide outward, never from the helices inward") — the same geometry Liu,
Johansen and Householder converged on. The hard gates in section 6 —
`ipTM(G12D) ≥ 0.75`, `ipTM(WT) ≤ 0.50`, `ipTM(apo) ≤ 0.50`,
`specificity_gap ≥ 0.25` — **already implement counter-screening**, and the
three-way prediction in section 5 is structurally the same comparison Johansen
performed as in silico cross-panning. Two concrete changes would bring the
pipeline in line with the strongest results above.

1. **Score the buried fraction of the mutated Asp (peptide position 6) as a
   design objective, not only as a contact count.** The pipeline currently
   requires ≥ 3 binder heavy atoms within 5 Å of position 6 and ≥ 40 % buried
   area on the peptide chain. Add a per-design quantity — the fraction of the
   Asp6 side-chain solvent-accessible surface buried by the binder — and use it
   both as a gate and as a ranking term, in the spirit of Xie's hotspot-centred
   design and Householder's up-facing-side-chain focus. A binder that buries
   the carboxylate pays a real desolvation penalty on the glycine variant; a
   binder that merely touches position 6 does not. Where compute allows,
   approximate Xie's reverse-mutation ΔΔG by rescoring the top designs with the
   D→G substitution in place. Carry over Xie's sequence lesson: bias the
   interface toward Arg/Lys/His/Ser/Thr/Tyr donors around position 6 (section 3
   already does) and keep net interface negative charge low.

2. **Add an unrelated same-allele decoy peptide to the counter-screen panel.**
   The current panel is G12D vs WT vs apo. That catches apo binders and
   position-6-blind binders, but it does not catch a design that binds the
   HLA-A\*11:01 groove rim plus generic peptide backbone — which would pass the
   WT gate only by luck. Following Johansen and Householder, add at least one
   well-characterised, unrelated HLA-A\*11:01 10-mer (and ideally a 9-mer) as a
   fourth prediction target, and extend the gate to
   `max(ipTM(WT), ipTM(apo), ipTM(decoy)) ≤ 0.50`, redefining
   `specificity_gap = score(G12D) − max(score(WT), score(apo), score(decoy))`.
   This costs one extra prediction run per design and is the cheapest available
   insurance against the dominant failure mode.

A third item, if time permits before the 8 Nov deadline: the design stage
currently generates against G12D only and filters afterwards. Odin-Multi's
formulation — optimising one sequence with explicit repulsion from the WT and
apo structures — would move specificity into the generator, and is the single
change most likely to raise the yield of gate-passing designs per GPU hour.
Given the submission caps (40/20/20 designs) this matters: the binding
constraint is design *quality per slot*, not throughput.
