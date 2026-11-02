#!/usr/bin/env python3
"""Run the Challenge 5 design pipeline end to end against NVIDIA BioNeMo NIMs.

Chains the three hosted NIMs (see `src/nvidia_client.py`) into the pipeline of
`docs/PIPELINE.md`:

    1. RFdiffusion   - binder backbones scaffolded onto the peptide-facing
                       surface of 7OW6, hotspots on peptide chain C 4-8 and
                       weighted on the G12D aspartate at position 6.
    2. ProteinMPNN   - sequences for the designed chain only, target fixed,
                       with the section 3 sequence filters applied.
    3. Boltz-2       - every surviving sequence predicted against all four
                       target states (G12D, WT, apo, decoy).
    4. design_gates  - the section 6 hard gates, written to a scored CSV with
                       the same columns as examples/scored_designs_example.csv.

Usage:
    python3 src/run_pipeline.py --dry-run                 # no API calls at all
    python3 src/run_pipeline.py --backbones 4 --seqs-per-backbone 8 --out scored.csv

`--dry-run` exercises every code path with deterministic stub responses and
never constructs a client, so it needs neither an API key nor quota.

The NVIDIA API key is read from the environment (`NVIDIA_API_KEY`) by
`nvidia_client`; this script never takes it as an argument and never logs it.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import random
import re
import sys
import time
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

import complex_analysis
import design_gates
from nvidia_client import NvidiaApiError, NvidiaClient, parse_mpnn_mfasta

REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_TARGET_DIR = REPO_ROOT / "data" / "targets"

# --- target geometry (7OW6, as written by src/fetch_target.py) --------------
HLA_CHAIN = "A"
B2M_CHAIN = "B"
PEPTIDE_CHAIN = "C"
BINDER_CHAIN = "D"           # the chain id we give the binder in Boltz-2 input
ALPHA12_LAST_RESIDUE = 180   # a1/a2 groove domain of the heavy chain
MUTATED_PEPTIDE_POSITION = 6  # 1-based, the G12D aspartate
HOTSPOT_RESIDUES = [f"{PEPTIDE_CHAIN}{i}" for i in range(4, 9)]  # C4..C8

# In RFdiffusion NIM output the designed chain is A and the kept target
# segments are concatenated into chain B, in contig order.
RFD_BINDER_CHAIN = "A"
RFD_TARGET_CHAIN = "B"

# Interface definitions live in complex_analysis so the analysis module is the
# single source of truth; re-exported here for callers and tests.
CONTACT_CUTOFF = complex_analysis.CONTACT_CUTOFF
INTERFACE_CUTOFF = complex_analysis.INTERFACE_CUTOFF
SASA_PROBE_RADIUS = complex_analysis.SASA_PROBE_RADIUS
ASP_SIDECHAIN_ATOMS = ("CB", "CG", "OD1", "OD2")

STATES = ("g12d", "wt", "apo", "decoy")
STATE_FASTA = {
    "g12d": "g12d_complex.fasta",
    "wt": "wt_complex.fasta",
    "apo": "apo_hla_b2m.fasta",
    "decoy": "decoy_a1101_complex.fasta",
}

CSV_COLUMNS = [
    "name",
    "sequence",
    "molecule_class",
    "iptm_g12d",
    "iptm_wt",
    "iptm_apo",
    "iptm_decoy",
    "specificity_gap",
    "peptide_buried_frac",
    "asp6_buried_frac",
    "net_interface_charge",
    "pos6_contacts",
    "pae_interface_g12d",
]

# Formal charge per residue at pH 7.4; His counted at ~10 % protonation.
RESIDUE_CHARGE = complex_analysis.RESIDUE_CHARGE

N_GLYC_SEQUON = re.compile(r"N[^P][ST]")
RUN_OF_FIVE = re.compile(r"(.)\1{4,}")


# ---------------------------------------------------------------------------
# target preparation
# ---------------------------------------------------------------------------

def read_fasta(path: Path) -> List[str]:
    """Sequences from a (multi-record) FASTA, in file order."""
    sequences: List[str] = []
    for line in path.read_text().splitlines():
        line = line.strip()
        if line.startswith(">"):
            sequences.append("")
        elif line and sequences:
            sequences[-1] += line
    return sequences


def load_target_states(target_dir: Path) -> Dict[str, List[Dict[str, str]]]:
    """Polymer lists (minus the binder) for each of the four target states."""
    states: Dict[str, List[Dict[str, str]]] = {}
    for state, filename in STATE_FASTA.items():
        path = target_dir / filename
        if not path.exists():
            raise FileNotFoundError(
                f"{path} is missing - run `python3 src/fetch_target.py` first"
            )
        sequences = [s for s in read_fasta(path) if s]
        if state == "apo":
            if len(sequences) != 2:
                raise ValueError(f"{path}: expected HLA + B2M, got {len(sequences)}")
            chains = [HLA_CHAIN, B2M_CHAIN]
        else:
            if len(sequences) != 3:
                raise ValueError(
                    f"{path}: expected HLA + B2M + peptide, got {len(sequences)}"
                )
            chains = [HLA_CHAIN, B2M_CHAIN, PEPTIDE_CHAIN]
        states[state] = [
            {"id": cid, "molecule_type": "protein", "sequence": seq}
            for cid, seq in zip(chains, sequences)
        ]
    return states


def trimmed_design_target(target_pdb: Path) -> str:
    """PDB text of the design target: a1/a2 groove + peptide only.

    B2M and the a3 domain are dropped before RFdiffusion: they are 40 A from the
    epitope, cost diffusion time, and give the model surface it must not use.
    """
    from Bio.PDB import PDBIO, PDBParser, Select

    structure = PDBParser(QUIET=True).get_structure("target", str(target_pdb))

    class GrooveAndPeptide(Select):
        def accept_chain(self, chain):
            return 1 if chain.id in (HLA_CHAIN, PEPTIDE_CHAIN) else 0

        def accept_residue(self, residue):
            chain_id = residue.get_parent().id
            if chain_id == PEPTIDE_CHAIN:
                return 1
            return 1 if residue.id[1] <= ALPHA12_LAST_RESIDUE else 0

    import io as _io

    writer = PDBIO()
    writer.set_structure(structure)
    handle = _io.StringIO()
    writer.save(handle, GrooveAndPeptide())
    return handle.getvalue()


def binder_contig(binder_length: int, peptide_length: int = 10) -> str:
    """Contig string for binder design against the trimmed groove + peptide."""
    return (
        f"{HLA_CHAIN}1-{ALPHA12_LAST_RESIDUE}/0 "
        f"{PEPTIDE_CHAIN}1-{peptide_length}/0 "
        f"{binder_length}-{binder_length}"
    )


# ---------------------------------------------------------------------------
# sequence filters (PIPELINE.md section 3)
# ---------------------------------------------------------------------------

def sequence_filter_failures(sequence: str) -> List[str]:
    """Section 3 sequence filters; empty list means the sequence is acceptable."""
    failures: List[str] = []
    charge = sum(
        1.0 if aa in "KR" else -1.0 if aa in "DE" else 0.1 if aa == "H" else 0.0
        for aa in sequence
    )
    if not -6.0 <= charge <= 6.0:
        failures.append(f"net charge {charge:+.1f} outside -6..+6")
    if sequence.count("C") % 2 == 1:
        failures.append("odd number of cysteines (unpaired Cys)")
    if N_GLYC_SEQUON.search(sequence):
        failures.append("N-glycosylation sequon (N-X-S/T)")
    if RUN_OF_FIVE.search(sequence):
        failures.append("more than 4 consecutive identical residues")
    return failures


# ---------------------------------------------------------------------------
# structure metrics from the predicted complex
# ---------------------------------------------------------------------------

def _parse_mmcif(text: str):
    from Bio.PDB import MMCIFParser
    import io as _io

    return MMCIFParser(QUIET=True).get_structure(
        "complex", _io.StringIO(text)
    )


def _heavy_atoms(entity) -> List:
    return [a for a in entity.get_atoms() if a.element != "H"]


def complex_metrics(
    mmcif_text: str,
    binder_chain: str = BINDER_CHAIN,
    peptide_chain: Optional[str] = None,
    peptide_position: int = MUTATED_PEPTIDE_POSITION,
    binder_sequence: Optional[str] = None,
) -> Dict[str, float]:
    """Interface geometry of a predicted binder/pMHC complex.

    Thin wrapper over :func:`complex_analysis.interface_metrics`, which does
    the real work: it maps the chains by *content* (sequence and length) rather
    than by the hard-coded ids this function used to assume, and raises instead
    of silently returning 0.0 when a chain, a residue or a side chain is
    missing. See that module's docstring for why - batch-01's all-zero buried
    fractions were indistinguishable from a mapping bug.

    `binder_chain` is only a *hint*, used when `binder_sequence` is not given;
    `peptide_chain` is accepted for backwards compatibility and, when passed,
    is checked against the chain the mapping actually identified.

    Returns `peptide_buried_frac`, `asp6_buried_frac`, `pos6_contacts` and
    `net_interface_charge` as defined in `docs/PIPELINE.md` section 6, plus the
    per-chain buried areas, contact counts and the chain map behind them. SASA
    is Shrake-Rupley at a 1.4 A probe, computed twice on the *same predicted
    coordinates*: once with the binder present and once with the binder chain
    removed, so the two numbers differ only by the binder.
    """
    metrics = complex_analysis.interface_metrics(
        mmcif_text,
        binder_sequence=binder_sequence,
        binder_chain_hint=binder_chain,
        peptide_position=peptide_position,
    )
    if peptide_chain is not None and metrics["peptide_chain"] != peptide_chain:
        raise complex_analysis.ChainMappingError(
            f"peptide was expected in chain {peptide_chain!r} but the sequence "
            f"mapping puts it in chain {metrics['peptide_chain']!r} "
            f"(chains: {metrics['chain_map']['chains']})"
        )
    return metrics


def _count_contacts(model, binder_chain: str, residue) -> int:
    """Binder heavy atoms within CONTACT_CUTOFF of `residue`'s heavy atoms."""
    return complex_analysis.count_contacts(model, binder_chain, residue)


def _interface_charge(model, binder_chain: str) -> float:
    """Net formal charge of the binder residues forming the interface."""
    return complex_analysis.interface_charge(model, binder_chain)


# ---------------------------------------------------------------------------
# Boltz-2 score extraction
# ---------------------------------------------------------------------------

def binder_interface_iptm(
    payload: Mapping[str, Any], binder_index: int, sample: int = 0
) -> Optional[float]:
    """Binder-vs-target ipTM from `pair_chains_iptm_scores`.

    The top-level `iptm_scores` of a pMHC + binder prediction is dominated by
    the three *native* interfaces (HLA/B2M/peptide) and barely moves when the
    binder misses, so it is useless as a binding score here. The per-chain-pair
    matrix is used instead: the score for a state is the **largest** ipTM
    between the binder chain and any target chain, i.e. the strongest interface
    the model predicts the binder to make with that state.
    """
    pairs = payload.get("pair_chains_iptm_scores")
    if not pairs or sample >= len(pairs):
        values = payload.get("iptm_scores") or []
        return float(values[sample]) if sample < len(values) else None
    matrix = pairs[sample]
    row = matrix.get(str(binder_index), matrix.get(binder_index, {}))
    best: Optional[float] = None
    for key, value in row.items():
        if str(key) == str(binder_index):
            continue  # the binder's own intra-chain pTM
        value = float(value)
        best = value if best is None else max(best, value)
    return best


def _first(payload: Mapping[str, Any], key: str) -> Optional[float]:
    values = payload.get(key) or []
    return float(values[0]) if values else None


# ---------------------------------------------------------------------------
# dry-run stubs
# ---------------------------------------------------------------------------

class DryRunClient:
    """Deterministic stand-in for `NvidiaClient`; makes no network calls.

    Shapes its return values exactly like the real endpoints (schemas
    documented in `nvidia_client`), so `--dry-run` exercises the same parsing,
    metric and gate code as a live run.
    """

    def __init__(self, seed: int = 0) -> None:
        self._random = random.Random(seed)
        self.calls: List[str] = []

    def rfdiffusion_generate(
        self, input_pdb: str, contigs: str, hotspot_res=None,
        diffusion_steps: int = 15, random_seed: Optional[int] = None,
    ) -> str:
        self.calls.append("rfdiffusion")
        length = int(contigs.split()[-1].split("-")[0])
        lines = []
        serial = 1
        for index in range(length):
            for name, element in (("N", "N"), ("CA", "C"), ("C", "C"), ("O", "O")):
                lines.append(
                    f"ATOM  {serial:>5} {name:<4}GLY {RFD_BINDER_CHAIN}"
                    f"{index + 1:>4}    "
                    f"{index * 3.8:>8.3f}{0.0:>8.3f}{0.0:>8.3f}  1.00  0.00"
                )
                serial += 1
        lines.append("TER")
        return "\n".join(lines) + "\n"

    def proteinmpnn_design(
        self, input_pdb: str, input_pdb_chains=None, num_seq_per_target: int = 8,
        sampling_temp=(0.1,), ca_only: bool = False,
        use_soluble_model: bool = True, random_seed: Optional[int] = None,
    ) -> Dict[str, Any]:
        self.calls.append("proteinmpnn")
        length = sum(
            1 for line in input_pdb.splitlines()
            if line.startswith("ATOM") and line[12:16].strip() == "CA"
        ) or 60
        alphabet = "AEKLRSTVYWQNGIDFHMP"
        records = [">input, designed_chains=['A']", "G" * length]
        for sample in range(1, num_seq_per_target + 1):
            sequence = "".join(self._random.choice(alphabet) for _ in range(length))
            records.append(
                f">T={sampling_temp[0]}, sample={sample}, "
                f"score={0.9 + sample / 100:.4f}, seq_recovery=0.3"
            )
            records.append(sequence)
        mfasta = "\n".join(records) + "\n"
        return {
            "mfasta": mfasta,
            "sequences": parse_mpnn_mfasta(mfasta),
            "scores": [0.9] * num_seq_per_target,
        }

    def boltz2_predict(self, polymers, **kwargs) -> Dict[str, Any]:
        self.calls.append("boltz2")
        count = len(polymers)
        binder_index = count - 1
        peptides = [
            str(p["sequence"]) for p in polymers if len(str(p["sequence"])) <= 15
        ]
        on_target = any(p == "VVVGADGVGK" for p in peptides)
        matrix = {
            str(i): {
                str(j): (
                    0.95 if i != binder_index and j != binder_index
                    else (0.80 if on_target else 0.30)
                )
                for j in range(count)
            }
            for i in range(count)
        }
        structure = _stub_mmcif(polymers)
        return {
            "structures": [
                {"structure": structure, "format": "mmcif", "name": "stub",
                 "source": "dry-run"}
            ],
            "iptm_scores": [0.79],
            "ptm_scores": [0.80],
            "confidence_scores": [0.70],
            "complex_plddt_scores": [0.88],
            "complex_iplddt_scores": [0.85],
            "complex_pde_scores": [3.0],
            "complex_ipde_scores": [4.3],
            "pair_chains_iptm_scores": [matrix],
            "metrics": {"total_time_seconds": 0.0},
            "pae": None,
            "pde": None,
        }


def _stub_mmcif(polymers: Sequence[Mapping[str, Any]]) -> str:
    """A minimal but genuinely parseable mmCIF: one CB-bearing residue trace.

    Coordinates are a straight line per chain with the chains stacked 6 A apart,
    which is enough for the SASA / contact code to run and produce numbers; the
    numbers themselves are meaningless, which is the point of a dry run.
    """
    header = [
        "data_stub",
        "#",
        "loop_",
        "_atom_site.group_PDB",
        "_atom_site.id",
        "_atom_site.type_symbol",
        "_atom_site.label_atom_id",
        "_atom_site.label_alt_id",
        "_atom_site.label_comp_id",
        "_atom_site.label_asym_id",
        "_atom_site.label_seq_id",
        "_atom_site.pdbx_PDB_ins_code",
        "_atom_site.Cartn_x",
        "_atom_site.Cartn_y",
        "_atom_site.Cartn_z",
        "_atom_site.occupancy",
        "_atom_site.B_iso_or_equiv",
        "_atom_site.auth_seq_id",
        "_atom_site.auth_asym_id",
        "_atom_site.pdbx_PDB_model_num",
    ]
    three = {
        "A": "ALA", "C": "CYS", "D": "ASP", "E": "GLU", "F": "PHE", "G": "GLY",
        "H": "HIS", "I": "ILE", "K": "LYS", "L": "LEU", "M": "MET", "N": "ASN",
        "P": "PRO", "Q": "GLN", "R": "ARG", "S": "SER", "T": "THR", "V": "VAL",
        "W": "TRP", "Y": "TYR", "X": "GLY",
    }
    rows: List[str] = []
    serial = 1
    for chain_index, polymer in enumerate(polymers):
        chain_id = str(polymer.get("id") or chr(ord("A") + chain_index))
        for residue_index, aa in enumerate(str(polymer["sequence"])):
            comp = three.get(aa.upper(), "GLY")
            atoms = [("N", "N"), ("CA", "C"), ("C", "C"), ("O", "O")]
            if comp == "ASP":
                atoms += [("CB", "C"), ("CG", "C"), ("OD1", "O"), ("OD2", "O")]
            elif comp != "GLY":
                atoms.append(("CB", "C"))
            for offset, (name, element) in enumerate(atoms):
                x = residue_index * 3.8 + offset * 0.4
                y = chain_index * 6.0
                z = offset * 0.3
                rows.append(
                    f"ATOM {serial} {element} {name} . {comp} {chain_id} "
                    f"{residue_index + 1} ? {x:.3f} {y:.3f} {z:.3f} 1.00 50.00 "
                    f"{residue_index + 1} {chain_id} 1"
                )
                serial += 1
    return "\n".join(header + rows) + "\n#\n"


# ---------------------------------------------------------------------------
# pipeline
# ---------------------------------------------------------------------------

def design_name(sequence: str, index: int) -> str:
    digest = hashlib.sha1(sequence.encode()).hexdigest()[:6]
    return f"KRASG12D_NIM_{index:03d}_{digest}"


def score_sequence(
    client,
    sequence: str,
    target_states: Mapping[str, List[Dict[str, str]]],
    boltz_kwargs: Mapping[str, Any],
    states: Sequence[str] = STATES,
) -> Tuple[Dict[str, Any], Dict[str, Any]]:
    """Predict one binder against every target state and score the interface.

    Returns `(row, detail)`: `row` has the scored-CSV columns, `detail` keeps
    the raw per-state confidences for the sidecar JSON.
    """
    row: Dict[str, Any] = {"sequence": sequence, "molecule_class": "single_chain"}
    detail: Dict[str, Any] = {"sequence": sequence, "states": {}}

    for state in states:
        polymers = list(target_states[state]) + [
            {"id": BINDER_CHAIN, "molecule_type": "protein", "sequence": sequence}
        ]
        payload = client.boltz2_predict(polymers, **boltz_kwargs)
        binder_index = len(polymers) - 1
        iptm = binder_interface_iptm(payload, binder_index)
        row[design_gates.IPTM_COLUMNS[state]] = (
            round(iptm, 4) if iptm is not None else ""
        )
        detail["states"][state] = {
            "binder_interface_iptm": iptm,
            "iptm_scores": payload.get("iptm_scores"),
            "complex_plddt": _first(payload, "complex_plddt_scores"),
            "complex_iplddt": _first(payload, "complex_iplddt_scores"),
            "complex_ipde": _first(payload, "complex_ipde_scores"),
            "pair_chains_iptm_scores": payload.get("pair_chains_iptm_scores"),
        }
        if state == "g12d":
            structure = payload["structures"][0]["structure"]
            try:
                metrics = complex_metrics(structure)
            except Exception as exc:  # noqa: BLE001 - never lose a whole run to one design
                detail["metrics_error"] = f"{type(exc).__name__}: {exc}"
                metrics = {}
            detail["g12d_metrics"] = metrics
            row["peptide_buried_frac"] = _round(metrics.get("peptide_buried_frac"))
            row["asp6_buried_frac"] = _round(metrics.get("asp6_buried_frac"))
            row["net_interface_charge"] = _round(metrics.get("net_interface_charge"), 1)
            contacts = metrics.get("pos6_contacts")
            row["pos6_contacts"] = int(contacts) if contacts is not None else ""
            # The hosted Boltz-2 does not return per-residue PAE (see
            # nvidia_client.boltz2_predict), so the interface *PDE* stands in.
            row["pae_interface_g12d"] = _round(
                detail["states"]["g12d"]["complex_ipde"], 2
            )

    gap = design_gates.specificity_gap(row)
    row["specificity_gap"] = _round(gap)
    row["gate_failures"] = design_gates.evaluate(row)
    return row, detail


def _round(value: Optional[float], digits: int = 3):
    return round(float(value), digits) if value is not None else ""


def run(args: argparse.Namespace) -> int:
    target_dir = Path(args.target_dir)
    target_states = load_target_states(target_dir)

    if args.dry_run:
        client: Any = DryRunClient(seed=args.seed)
        target_pdb_text = "DRY-RUN"
        print("[pipeline] dry run: no API calls, stub responses")
    else:
        target_pdb = target_dir / "g12d_complex.pdb"
        if not target_pdb.exists():
            raise FileNotFoundError(
                f"{target_pdb} is missing - run `python3 src/fetch_target.py` "
                f"(the PDB files are git-ignored)"
            )
        target_pdb_text = trimmed_design_target(target_pdb)
        client = NvidiaClient(timeout=args.timeout)
        print(
            f"[pipeline] live run against NVIDIA BioNeMo NIMs; "
            f"target {target_pdb.name} trimmed to the a1/a2 groove + peptide "
            f"({len(target_pdb_text)} bytes)"
        )

    boltz_kwargs = {
        "recycling_steps": args.recycling_steps,
        "sampling_steps": args.sampling_steps,
        "diffusion_samples": 1,
        "without_potentials": True,
    }

    lengths = [int(x) for x in args.binder_lengths.split(",")]
    rows: List[Dict[str, Any]] = []
    details: List[Dict[str, Any]] = []
    candidates: List[str] = []
    seen: set = set()
    started = time.time()

    for backbone_index in range(args.backbones):
        length = lengths[backbone_index % len(lengths)]
        contigs = binder_contig(length)
        print(
            f"[rfdiffusion] backbone {backbone_index + 1}/{args.backbones}, "
            f"length {length}, hotspots {','.join(HOTSPOT_RESIDUES)}"
        )
        try:
            backbone_pdb = client.rfdiffusion_generate(
                input_pdb=target_pdb_text,
                contigs=contigs,
                hotspot_res=HOTSPOT_RESIDUES,
                diffusion_steps=args.diffusion_steps,
                random_seed=args.seed + backbone_index,
            )
        except NvidiaApiError as exc:
            print(f"[rfdiffusion] FAILED: {exc}", file=sys.stderr)
            continue

        if args.save_backbones:
            out_dir = Path(args.save_backbones)
            out_dir.mkdir(parents=True, exist_ok=True)
            (out_dir / f"backbone_{backbone_index:03d}.pdb").write_text(backbone_pdb)

        try:
            designed = client.proteinmpnn_design(
                input_pdb=backbone_pdb,
                input_pdb_chains=[RFD_BINDER_CHAIN],
                num_seq_per_target=args.seqs_per_backbone,
                sampling_temp=[args.sampling_temp],
                use_soluble_model=True,
                random_seed=args.seed + backbone_index,
            )
        except NvidiaApiError as exc:
            print(f"[proteinmpnn] FAILED: {exc}", file=sys.stderr)
            continue

        for record in designed["sequences"]:
            sequence = record["sequence"].split("/")[0]  # designed chain only
            failures = sequence_filter_failures(sequence)
            if failures:
                print(f"[filter] drop {sequence[:20]}...: {'; '.join(failures)}")
                continue
            if sequence in seen:
                continue
            seen.add(sequence)
            candidates.append(sequence)

    if args.max_candidates:
        candidates = candidates[: args.max_candidates]
    print(f"[pipeline] {len(candidates)} candidate sequences to score "
          f"({len(candidates) * len(STATES)} Boltz-2 calls)")

    for index, sequence in enumerate(candidates, start=1):
        print(f"[boltz2] scoring {index}/{len(candidates)} against "
              f"{', '.join(STATES)}")
        try:
            row, detail = score_sequence(
                client, sequence, target_states, boltz_kwargs
            )
        except NvidiaApiError as exc:
            print(f"[boltz2] FAILED for candidate {index}: {exc}", file=sys.stderr)
            continue
        row["name"] = design_name(sequence, index)
        detail["name"] = row["name"]
        rows.append(row)
        details.append(detail)

    rows.sort(
        key=lambda r: (r.get("specificity_gap") or -99, r.get("asp6_buried_frac") or 0),
        reverse=True,
    )

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    write_scored_csv(rows, out_path)
    sidecar = out_path.with_suffix(".details.json")
    sidecar.write_text(json.dumps(details, indent=2, default=str))

    passing = [r for r in rows if not r["gate_failures"]]
    print(
        f"[pipeline] wrote {len(rows)} scored designs to {out_path} "
        f"({len(passing)} clear every hard gate), details in {sidecar}, "
        f"{time.time() - started:.1f}s elapsed"
    )
    for row in rows:
        if row["gate_failures"]:
            print(f"  - {row['name']}: {'; '.join(row['gate_failures'])}")
    return 0


def write_scored_csv(rows: Iterable[Mapping[str, Any]], path: Path) -> None:
    """Write the scored table with exactly the example CSV's columns."""
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=CSV_COLUMNS, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row.get(column, "") for column in CSV_COLUMNS})


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--dry-run", action="store_true",
                        help="stub every NIM response; no network, no quota, no key")
    parser.add_argument("--backbones", type=int, default=2,
                        help="RFdiffusion backbones to generate (default 2)")
    parser.add_argument("--binder-lengths", default="70,85",
                        help="comma-separated binder lengths to cycle through")
    parser.add_argument("--seqs-per-backbone", type=int, default=4,
                        help="ProteinMPNN sequences per backbone (1-100)")
    parser.add_argument("--sampling-temp", type=float, default=0.15,
                        help="ProteinMPNN sampling temperature (0.1-0.3)")
    parser.add_argument("--diffusion-steps", type=int, default=15,
                        help="RFdiffusion denoising steps (1-50)")
    parser.add_argument("--recycling-steps", type=int, default=3,
                        help="Boltz-2 recycling steps (1-10)")
    parser.add_argument("--sampling-steps", type=int, default=50,
                        help="Boltz-2 sampling steps (10-1000)")
    parser.add_argument("--max-candidates", type=int, default=0,
                        help="cap the sequences sent to Boltz-2 (0 = no cap)")
    parser.add_argument("--target-dir", default=str(DEFAULT_TARGET_DIR))
    parser.add_argument("--out", default="designs/scored_designs.csv")
    parser.add_argument("--save-backbones", default="",
                        help="directory to write RFdiffusion backbones to")
    parser.add_argument("--timeout", type=float, default=600.0,
                        help="per-request HTTP timeout in seconds")
    parser.add_argument("--seed", type=int, default=1)
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    return run(args)


if __name__ == "__main__":
    raise SystemExit(main())
