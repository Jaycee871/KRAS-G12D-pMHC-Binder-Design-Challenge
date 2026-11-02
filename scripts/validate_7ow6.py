#!/usr/bin/env python3
"""Validate the buried-area code against a real complex: PDB 7OW6.

Problem 2 of the batch-01 post-mortem: every design reported
`peptide_buried_frac = 0.00` and `asp6_buried_frac = 0.00`. Before trusting any
replacement, the analysis has to reproduce a *known* answer. 7OW6 is the
HLA-A*11:01 / β2m / VVVGADGVGK (KRAS G12D) complex bound by the selective TCR
of Poole et al. 2022 - the same target this project designs against - so it is
the natural positive control:

* chain A  HLA-A*11:01 heavy chain
* chain B  β2-microglobulin
* chain C  the 10-mer peptide VVVGADGVGK (position 6 = the G12D aspartate)
* chains D/E  the TCR α and β chains

Two checks are run, each treating a real binding partner as "the binder":

1. **HLA groove as the binder.** The peptide sits in the groove, so burying it
   should be nearly total and the mapping must find position 6 as an ASP.
2. **TCR (chains D+E) as the binder.** This is the geometry the designs are
   meant to imitate: a protein landing on the solvent-exposed face of the
   peptide. The peptide share of the buried area and the Asp6 burial should
   both be substantial.

    python3 scripts/validate_7ow6.py            # downloads 7OW6 if needed
"""
from __future__ import annotations

import argparse
import io
import json
import sys
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import complex_analysis as ca  # noqa: E402

REPO_ROOT = Path(__file__).resolve().parents[1]
PDB_ID = "7OW6"
RCSB_URL = "https://files.rcsb.org/download/{pdb_id}.cif"
DEFAULT_CACHE = REPO_ROOT / "data" / "validation" / "7ow6.cif"


def fetch(pdb_id: str, cache: Path) -> str:
    if cache.exists():
        return cache.read_text()
    cache.parent.mkdir(parents=True, exist_ok=True)
    url = RCSB_URL.format(pdb_id=pdb_id.upper())
    print(f"[7ow6] downloading {url}")
    with urllib.request.urlopen(url, timeout=120) as response:
        text = response.read().decode()
    cache.write_text(text)
    return text


def subset(cif_text: str, keep_chains, merge_into=None) -> str:
    """mmCIF text containing only `keep_chains`, optionally merged to one id.

    Merging is how a two-chain binder (the TCR αβ heterodimer) is presented to
    `interface_metrics`, which treats a single chain as the binder.
    """
    from Bio.PDB import MMCIFIO, MMCIFParser, Select

    structure = MMCIFParser(QUIET=True).get_structure("x", io.StringIO(cif_text))
    model = next(iter(structure))
    if merge_into:
        source, destination = merge_into
        src_chain = model[source]
        dst_chain = model[destination]
        offset = max(r.id[1] for r in dst_chain) + 100
        for residue in list(src_chain):
            src_chain.detach_child(residue.id)
            residue.id = (residue.id[0], residue.id[1] + offset, residue.id[2])
            dst_chain.add(residue)
        model.detach_child(source)
        keep_chains = [c for c in keep_chains if c != source]

    class Keep(Select):
        def accept_chain(self, chain):
            return 1 if chain.id in keep_chains else 0

        def accept_residue(self, residue):
            return 1 if ca.is_polymer_residue(residue) else 0

        def accept_atom(self, atom):
            return 1 if atom.element != "H" and atom.get_altloc() in (" ", "A") else 0

    writer = MMCIFIO()
    writer.set_structure(structure)
    handle = io.StringIO()
    writer.save(handle, Keep())
    return handle.getvalue()


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--cache", default=str(DEFAULT_CACHE))
    parser.add_argument("--out", default=str(REPO_ROOT / "results" / "validation"
                                             / "7ow6_buried_area.json"))
    args = parser.parse_args(argv)

    cif = fetch(PDB_ID, Path(args.cache))
    report = {"pdb_id": PDB_ID, "cases": {}}

    print("[7ow6] case 1: HLA heavy chain A as the 'binder' "
          "(peptide in the groove)")
    groove = ca.interface_metrics(
        subset(cif, ["A", "B", "C"]), binder_chain_hint="A",
    )
    report["cases"]["hla_groove_as_binder"] = groove
    print(json.dumps(_summary(groove), indent=2))

    print("[7ow6] case 2: TCR (chains D+E, merged) as the binder "
          "(the geometry the designs imitate)")
    tcr = ca.interface_metrics(
        subset(cif, ["A", "B", "C", "D", "E"], merge_into=("E", "D")),
        binder_chain_hint="D",
    )
    report["cases"]["tcr_as_binder"] = tcr
    print(json.dumps(_summary(tcr), indent=2))

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2, default=str))
    print(f"[7ow6] wrote {out}")

    # Sanity assertions: a true complex must not look like batch-01 did.
    problems = []
    if groove["asp6_buried_frac"] < 0.6:
        problems.append(f"groove Asp6 burial {groove['asp6_buried_frac']:.2f} < 0.60")
    if groove["peptide_buried_frac"] < 0.2:
        problems.append("groove peptide share of buried area implausibly low")
    if tcr["pos6_contacts"] < 3:
        problems.append("TCR makes fewer than 3 contacts with position 6")
    if problems:
        print("[7ow6] VALIDATION FAILED: " + "; ".join(problems), file=sys.stderr)
        return 1
    print("[7ow6] validation OK")
    return 0


def _summary(metrics):
    return {
        "chains": {k: v["length"] for k, v in metrics["chain_map"]["chains"].items()},
        "mapping": {k: metrics["chain_map"][k]
                    for k in ("binder", "peptide", "hla", "b2m")},
        "pos6_residue": f"{metrics['pos6_residue']}{metrics['pos6_resseq']}",
        "peptide_buried_frac": round(metrics["peptide_buried_frac"], 3),
        "asp6_buried_frac": round(metrics["asp6_buried_frac"], 3),
        "asp6_sasa_free_A2": round(metrics["asp6_sasa_free"], 1),
        "asp6_sasa_bound_A2": round(metrics["asp6_sasa_bound"], 1),
        "buried_area_total_A2": round(metrics["buried_area_total"], 1),
        "buried_area_peptide_A2": round(metrics["buried_area_peptide"], 1),
        "buried_area_by_chain": metrics["buried_area_by_chain"],
        "pos6_contacts": metrics["pos6_contacts"],
        "contacts_by_chain": metrics["contacts_by_chain"],
        "net_interface_charge": round(metrics["net_interface_charge"], 1),
        "notes": metrics["notes"],
    }


if __name__ == "__main__":
    raise SystemExit(main())
