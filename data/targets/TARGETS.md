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

The two peptides differ only at position 6 (`D` vs `G`). The binder interface must cover that position.

Note: `wt_complex_template.pdb` carries G12D coordinates and is a modelling template only - rebuild position 6 before scoring. `wt_complex.fasta` is authoritative.
