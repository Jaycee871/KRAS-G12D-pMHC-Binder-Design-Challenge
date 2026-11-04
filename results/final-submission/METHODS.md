# Methods note for the KRAS G12D pMHC Track-1 package

The campaign used RFdiffusion/ProteinMPNN generation followed by repeated
Boltz-2 structural screening against KRAS G12D pMHC, wild-type KRAS pMHC,
peptide-free HLA-A*11:01/B2M and an unrelated-peptide counter-screen where
available.

The design hypothesis was direct recognition of the mutant Asp6 carboxylate.
During audit, the original direct-readout implementation was found to admit a
side-chain carbon as an acceptor in some poses. All decision-relevant evidence
was therefore remeasured with an oxygen-specific Asp6 OD1/OD2 criterion plus
positive Asp6 burial, while the original files were preserved.

A positive-control audit also showed that large WT pose/footprint loss is not a
valid universal specificity requirement for this target, so historical WT-loss,
ipTM-gap and single-pose MM-GBSA criteria were not used as final hard gates.

The final diagnostic batches tested peptide dependence directly. A matched
regional-anchor counter-screen showed greater apo-HLA buried area than G12D
pMHC for all three lead sequences and all eight lower-ranked retained review
sequences. A final peptide-only RFdiffusion experiment produced 9/24 backbones
with the predeclared Asp6+Val8 peptide-contact geometry, but 0/24 fit after the
crystallographic HLA groove was restored; severe HLA clashes ranged from 41 to
154 atoms under a <=2-atom predeclared threshold. Generation was therefore
stopped rather than retuning the threshold after observing the result.

The 20-design file is consequently an evidence-tiered fallback package. The
top 11 preserve the curated review order; the remaining nine are explicitly
documented fallback-only near-misses selected from previously excluded families
to satisfy the Track-1 minimum without fabricating new evidence.
