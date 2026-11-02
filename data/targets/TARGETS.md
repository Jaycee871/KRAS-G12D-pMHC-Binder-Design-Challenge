# Challenge 5 target states

Source structure: RCSB PDB `7OW6` (TCR chains stripped; waters and heteroatoms removed).

- HLA-A*11:01 heavy chain: chain `A`, 276 aa
- beta-2-microglobulin: chain `B`, 99 aa
- peptide: chain `C`, 10 aa

| state | peptide | bind? | files |
| --- | --- | --- | --- |
| G12D complex | `VVVGADGVGK` | YES | `g12d_complex.pdb`, `g12d_complex.fasta` |
| WT complex | `VVVGAGGVGK` | NO | `wt_complex_template.pdb`, `wt_complex.fasta` |
| peptide-free HLA+B2M | - | NO | `apo_hla_b2m.pdb`, `apo_hla_b2m.fasta` |
| decoy complex (unrelated A*11:01 10-mer) | `AVFDRKSDAK` | NO | `decoy_a1101_complex.fasta` |

The two peptides differ only at position 6 (`D` vs `G`). The binder interface must cover that position.

## Decoy state

`decoy_a1101_complex.fasta` is the counter-screen decoy added after the literature review (see `docs/LITERATURE.md` section 4 and `docs/PIPELINE.md` sections 1, 5 and 6). It reuses the HLA-A*11:01 heavy chain and beta-2-microglobulin sequences byte-for-byte from `g12d_complex.fasta`; only the chain `C` peptide differs.

Peptide `AVFDRKSDAK` - EBV EBNA-3B (EBNA-4) residues 399-408, the "AVF" epitope: one of the two immunodominant HLA-A*11:01-restricted EBV epitopes in A11-positive donors, alongside the nonamer `IVTDFSVIK` (EBNA-3B 416-424). Sources: de Campos-Lima et al. 1994, J Exp Med 179:1297 (doi:10.1084/jem.179.4.1297); Midgley et al. 2003, J Virol 77:11517 (doi:10.1128/JVI.77.21.11517-11530.2003).

Chosen because it is a 10-mer - the same length, groove register and bulge geometry as `VVVGADGVGK` - is unrelated in sequence, and carries its own acidic side chains (Asp4, Asp8), so it also flags a paratope that merely likes a carboxylate on the peptide face rather than this one at position 6. There is no deposited A*11:01 structure for it, so the complex is predicted rather than templated; `QVPLRPMTYK` (HIV-1 Nef 73-82, PDB 1QVO) is a crystallographically characterised A*11:01 decamer alternative if a template is needed.

The purpose of the decoy is negative: a design that scores well on it is binding the HLA-A*11:01 groove rim plus generic peptide backbone, not the G12D neoantigen. Gate: `max(ipTM(WT), ipTM(apo), ipTM(decoy)) <= 0.50`.

Note: `wt_complex_template.pdb` carries G12D coordinates and is a modelling template only - rebuild position 6 before scoring. `wt_complex.fasta` is authoritative.
