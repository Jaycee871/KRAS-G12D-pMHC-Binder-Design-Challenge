#!/usr/bin/env python3
"""Fetch PDB 7OW6 from RCSB and build the three Challenge 5 target states.

Targets produced under data/targets/:
  1. g12d_complex      - HLA-A*11:01 heavy chain + B2M + KRAS G12D peptide VVVGADGVGK
  2. wt_complex        - same, with peptide position 6 D->G (wild-type VVVGAGGVGK)
  3. apo_hla           - peptide-free HLA heavy chain + B2M

The TCR chains present in the crystal structure (D/E) are stripped: the design
target is the pMHC surface, not the ternary TCR complex.

Usage:
    python3 src/fetch_target.py [--outdir data/targets] [--pdb-id 7OW6]
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import requests
from Bio.PDB import PDBIO, PDBParser, Select
from Bio.PDB.Polypeptide import is_aa
from Bio.SeqUtils import seq1

RCSB_PDB_URL = "https://files.rcsb.org/download/{pdb_id}.pdb"

# Expected biology for 7OW6 (used for validation, not for fabricating data).
EXPECTED_PEPTIDE = "VVVGADGVGK"
WT_PEPTIDE = "VVVGAGGVGK"
MUTATED_POSITION = 6  # 1-based position in the 10-mer that differs (D vs G)


class ChainSelect(Select):
    """Write out only the requested chains, protein residues only."""

    def __init__(self, chain_ids):
        self.chain_ids = set(chain_ids)

    def accept_chain(self, chain):
        return 1 if chain.id in self.chain_ids else 0

    def accept_residue(self, residue):
        # Drop waters / heteroatoms so downstream tools get a clean target.
        return 1 if is_aa(residue, standard=True) else 0


def download_pdb(pdb_id: str, dest: Path) -> Path:
    """Download a PDB file from RCSB (cached on disk)."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0:
        print(f"[fetch] using cached {dest}")
        return dest
    url = RCSB_PDB_URL.format(pdb_id=pdb_id.upper())
    print(f"[fetch] GET {url}")
    resp = requests.get(url, timeout=60)
    resp.raise_for_status()
    dest.write_bytes(resp.content)
    print(f"[fetch] wrote {dest} ({dest.stat().st_size} bytes)")
    return dest


def chain_sequence(chain) -> str:
    """One-letter sequence of the standard amino acids in a chain."""
    return "".join(
        seq1(res.get_resname()) for res in chain if is_aa(res, standard=True)
    )


def classify_chains(structure):
    """Identify the HLA heavy chain, beta-2-microglobulin and the peptide.

    Classification is by sequence length and sequence content rather than by
    hard-coded chain IDs, so the script still works if RCSB re-letters chains.
    """
    model = next(iter(structure))
    seqs = {ch.id: chain_sequence(ch) for ch in model}
    seqs = {cid: s for cid, s in seqs.items() if s}

    peptide_id = None
    # The peptide is the short chain (<= 15 aa).
    short = sorted(
        (cid for cid, s in seqs.items() if len(s) <= 15),
        key=lambda cid: len(seqs[cid]),
    )
    if short:
        peptide_id = short[0]

    # B2M is ~99 aa; HLA class I heavy chain ectodomain is ~270-290 aa.
    b2m_id = None
    heavy_id = None
    for cid, s in seqs.items():
        if cid == peptide_id:
            continue
        if 85 <= len(s) <= 120 and b2m_id is None:
            b2m_id = cid
    heavy_candidates = [
        cid
        for cid, s in seqs.items()
        if cid not in {peptide_id, b2m_id} and 250 <= len(s) <= 320
    ]
    # The HLA heavy chain starts with the canonical GSHSMRYF motif.
    for cid in heavy_candidates:
        if seqs[cid].startswith("GSHSMRY") or "GSHSMRY" in seqs[cid][:10]:
            heavy_id = cid
            break
    if heavy_id is None and heavy_candidates:
        heavy_id = heavy_candidates[0]

    missing = [n for n, v in
               (("HLA heavy chain", heavy_id), ("B2M", b2m_id), ("peptide", peptide_id))
               if v is None]
    if missing:
        raise RuntimeError(
            f"could not identify: {', '.join(missing)}; chains found: "
            + ", ".join(f"{c}({len(s)}aa)" for c, s in seqs.items())
        )
    return heavy_id, b2m_id, peptide_id, seqs


def mutate_peptide(peptide: str, position: int, new_aa: str) -> str:
    """Return peptide with 1-based `position` replaced by `new_aa`."""
    if not 1 <= position <= len(peptide):
        raise ValueError(f"position {position} out of range for {len(peptide)}-mer")
    return peptide[: position - 1] + new_aa + peptide[position:]


def write_pdb(structure, chain_ids, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    io = PDBIO()
    io.set_structure(structure)
    io.save(str(path), ChainSelect(chain_ids))
    print(f"[write] {path}")


def write_fasta(records, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as fh:
        for header, seq in records:
            fh.write(f">{header}\n")
            for i in range(0, len(seq), 60):
                fh.write(seq[i : i + 60] + "\n")
    print(f"[write] {path}")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--pdb-id", default="7OW6")
    ap.add_argument("--outdir", default="data/targets", type=Path)
    ap.add_argument("--raw", default="data/raw", type=Path,
                    help="where the downloaded PDB is cached")
    args = ap.parse_args(argv)

    raw_pdb = download_pdb(args.pdb_id, args.raw / f"{args.pdb_id.upper()}.pdb")

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(args.pdb_id.upper(), str(raw_pdb))
    heavy_id, b2m_id, peptide_id, seqs = classify_chains(structure)

    heavy, b2m, peptide = seqs[heavy_id], seqs[b2m_id], seqs[peptide_id]
    print(f"[id] HLA heavy chain = {heavy_id} ({len(heavy)} aa)")
    print(f"[id] beta-2-microglobulin = {b2m_id} ({len(b2m)} aa)")
    print(f"[id] peptide = {peptide_id} ({len(peptide)} aa): {peptide}")

    if peptide != EXPECTED_PEPTIDE:
        print(
            f"[warn] extracted peptide {peptide} != expected {EXPECTED_PEPTIDE}; "
            "check the structure before using these targets",
            file=sys.stderr,
        )
    wt_peptide = mutate_peptide(peptide, MUTATED_POSITION, "G")
    if peptide == EXPECTED_PEPTIDE and wt_peptide != WT_PEPTIDE:
        raise RuntimeError(f"WT peptide construction failed: got {wt_peptide}")
    print(f"[id] wild-type counterpart (pos {MUTATED_POSITION} D->G): {wt_peptide}")

    out = args.outdir

    # Target 1: the G12D pMHC complex (the thing we want to bind).
    write_pdb(structure, [heavy_id, b2m_id, peptide_id], out / "g12d_complex.pdb")
    write_fasta(
        [
            (f"g12d_complex|HLA-A*11:01_heavy|chain_{heavy_id}", heavy),
            (f"g12d_complex|B2M|chain_{b2m_id}", b2m),
            (f"g12d_complex|KRAS_G12D_peptide|chain_{peptide_id}", peptide),
        ],
        out / "g12d_complex.fasta",
    )

    # Target 2: the wild-type pMHC complex (must NOT be bound).
    # Only the peptide sequence differs; the PDB is emitted with the G12D
    # coordinates as a modelling template - the position-6 side chain must be
    # rebuilt (ProteinMPNN/AF/Boltz) before use. The FASTA carries the real WT
    # sequence and is the authoritative input for complex prediction.
    write_pdb(structure, [heavy_id, b2m_id, peptide_id],
              out / "wt_complex_template.pdb")
    write_fasta(
        [
            (f"wt_complex|HLA-A*11:01_heavy|chain_{heavy_id}", heavy),
            (f"wt_complex|B2M|chain_{b2m_id}", b2m),
            (f"wt_complex|KRAS_WT_peptide|chain_{peptide_id}", wt_peptide),
        ],
        out / "wt_complex.fasta",
    )

    # Target 3: peptide-free HLA + B2M (must NOT be bound).
    write_pdb(structure, [heavy_id, b2m_id], out / "apo_hla_b2m.pdb")
    write_fasta(
        [
            (f"apo_hla_b2m|HLA-A*11:01_heavy|chain_{heavy_id}", heavy),
            (f"apo_hla_b2m|B2M|chain_{b2m_id}", b2m),
        ],
        out / "apo_hla_b2m.fasta",
    )

    summary = out / "TARGETS.md"
    summary.write_text(
        "# Challenge 5 target states\n\n"
        f"Source structure: RCSB PDB `{args.pdb_id.upper()}` "
        f"(TCR chains stripped; waters and heteroatoms removed).\n\n"
        f"- HLA-A*11:01 heavy chain: chain `{heavy_id}`, {len(heavy)} aa\n"
        f"- beta-2-microglobulin: chain `{b2m_id}`, {len(b2m)} aa\n"
        f"- peptide: chain `{peptide_id}`, {len(peptide)} aa\n\n"
        "| state | peptide | bind? | files |\n"
        "| --- | --- | --- | --- |\n"
        f"| G12D complex | `{peptide}` | YES | `g12d_complex.pdb`, `g12d_complex.fasta` |\n"
        f"| WT complex | `{wt_peptide}` | NO | `wt_complex_template.pdb`, `wt_complex.fasta` |\n"
        "| peptide-free HLA+B2M | - | NO | `apo_hla_b2m.pdb`, `apo_hla_b2m.fasta` |\n\n"
        f"The two peptides differ only at position {MUTATED_POSITION} "
        f"(`{peptide[MUTATED_POSITION-1]}` vs `{wt_peptide[MUTATED_POSITION-1]}`). "
        "The binder interface must cover that position.\n\n"
        "Note: `wt_complex_template.pdb` carries G12D coordinates and is a "
        "modelling template only - rebuild position "
        f"{MUTATED_POSITION} before scoring. `wt_complex.fasta` is authoritative.\n"
    )
    print(f"[write] {summary}")
    print("[done] 3 target states written to", out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
