#!/usr/bin/env python3
"""Positive control: the real 7OW6 TCR run through the *constrained* protocol.

Before trusting the pocket-constrained Boltz-2 protocol (see
`results/constrained-rescore/FINDINGS.md`) on a new design batch, run it on a
binder that is known to be correct and selective: the Poole et al. 2022 TCR of
PDB 7OW6, which discriminates KRAS G12D VVVGADGVGK / HLA-A*11:01 from the WT
peptide.

Protocol, identical to a design except that the binder is two chains:

* polymers  = the four target states from `data/targets` (chains A/B/C) plus
  the TCR alpha chain D and beta chain E taken from 7OW6;
* constraint = the flat NIM pocket term, one entry per binder chain, contacts
  on peptide C4-C8, max_distance 6.0 A, applied to every peptide-bearing state
  (g12d, wt, decoy) so the specificity gap is a like-for-like comparison; apo
  has no peptide and is run unconstrained;
* score      = max ipTM between a *binder* chain and a *target* chain. The
  binder-internal alpha/beta pair is excluded: for a one-chain design this is
  exactly `run_pipeline.binder_interface_iptm`;
* analysis   = `complex_analysis.interface_metrics` on the predicted G12D
  complex with D and E merged into one chain, which is the convention
  `results/validation/7ow6_buried_area.json` ("tcr_as_binder") already uses.

An unconstrained G12D prediction is run alongside as the A/B comparison.
"""
from __future__ import annotations

import io
import json
import sys
import time
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

import complex_analysis as ca            # noqa: E402
import design_gates                      # noqa: E402
import pose_filter                       # noqa: E402
import run_pipeline                      # noqa: E402
from nvidia_client import NvidiaClient, BOLTZ2_PATH   # noqa: E402

OUT = REPO / "results" / "validation"
OUT.mkdir(parents=True, exist_ok=True)
CKPT = OUT / "7ow6_constrained_checkpoint.jsonl"

TCR_ALPHA_CHAIN = "D"
TCR_BETA_CHAIN = "E"
MIN_INTERVAL = 8.0


def tcr_sequences(cif_path: Path):
    from Bio.PDB import MMCIFParser
    model = next(iter(MMCIFParser(QUIET=True).get_structure(
        "x", io.StringIO(cif_path.read_text()))))
    out = {}
    for chain in model:
        residues = [r for r in chain if ca.is_polymer_residue(r)]
        out[chain.id] = "".join(ca.one_letter(r) for r in residues)
    return out[TCR_ALPHA_CHAIN], out[TCR_BETA_CHAIN]


def merge_chains(cif_text: str, source: str, destination: str) -> str:
    """mmCIF text with chain `source` folded into `destination` (renumbered)."""
    from Bio.PDB import MMCIFIO, MMCIFParser
    structure = MMCIFParser(QUIET=True).get_structure("x", io.StringIO(cif_text))
    model = next(iter(structure))
    src, dst = model[source], model[destination]
    offset = max(r.id[1] for r in dst) + 100
    for residue in list(src):
        src.detach_child(residue.id)
        residue.id = (residue.id[0], residue.id[1] + offset, residue.id[2])
        dst.add(residue)
    model.detach_child(source)
    writer = MMCIFIO()
    writer.set_structure(structure)
    handle = io.StringIO()
    writer.save(handle)
    return handle.getvalue()


def binder_vs_target_iptm(payload, binder_indices, sample: int = 0):
    """max ipTM over (binder chain, target chain) pairs; binder-internal excluded."""
    pairs = payload.get("pair_chains_iptm_scores")
    if not pairs or sample >= len(pairs):
        return None
    matrix = pairs[sample]
    best = None
    for bi in binder_indices:
        row = matrix.get(str(bi), matrix.get(bi, {}))
        for key, value in row.items():
            if int(key) in binder_indices:
                continue
            value = float(value)
            best = value if best is None else max(best, value)
    return best


def main() -> int:
    alpha, beta = tcr_sequences(REPO / "data" / "validation" / "7ow6.cif")
    print(f"[control] TCR alpha {len(alpha)} aa, beta {len(beta)} aa")

    states = run_pipeline.load_target_states(REPO / "data" / "targets")
    contacts = [{"id": run_pipeline.PEPTIDE_CHAIN, "residue_index": i}
                for i in range(4, 9)]
    pocket = [{"binder": TCR_ALPHA_CHAIN, "contacts": contacts, "max_distance": 6.0},
              {"binder": TCR_BETA_CHAIN, "contacts": contacts, "max_distance": 6.0}]

    client = NvidiaClient()
    done = {}
    if CKPT.exists():
        for line in CKPT.open():
            rec = json.loads(line)
            done[rec["run"]] = rec

    def predict(state, constrained):
        polymers = [dict(p) for p in states[state]] + [
            {"id": TCR_ALPHA_CHAIN, "molecule_type": "protein", "sequence": alpha},
            {"id": TCR_BETA_CHAIN, "molecule_type": "protein", "sequence": beta},
        ]
        body = {"polymers": polymers, "recycling_steps": 3, "sampling_steps": 50,
                "diffusion_samples": 1, "step_scale": 1.638,
                "without_potentials": False}
        if constrained:
            body["constraints"] = pocket
        payload = client.post(BOLTZ2_PATH, body)
        binder_indices = {len(polymers) - 2, len(polymers) - 1}
        return payload, binder_indices

    runs = [("g12d_constrained", "g12d", True),
            ("wt_constrained", "wt", True),
            ("decoy_constrained", "decoy", True),
            ("apo_unconstrained", "apo", False),
            ("g12d_unconstrained", "g12d", False)]

    for name, state, constrained in runs:
        if name in done:
            print(f"[control] {name}: cached")
            continue
        t0 = time.time()
        rec = {"run": name, "state": state, "constrained": constrained}
        try:
            payload, binder_indices = predict(state, constrained)
            rec["iptm_binder_vs_target"] = binder_vs_target_iptm(payload, binder_indices)
            rec["iptm_scores"] = payload.get("iptm_scores")
            rec["pair_chains_iptm_scores"] = payload.get("pair_chains_iptm_scores")
            rec["complex_plddt"] = (payload.get("complex_plddt_scores") or [None])[0]
            rec["complex_iplddt"] = (payload.get("complex_iplddt_scores") or [None])[0]
            rec["complex_ipde"] = (payload.get("complex_ipde_scores") or [None])[0]
            structure = payload["structures"][0]["structure"]
            (OUT / f"7ow6_control_{name}.cif").write_text(structure)
            if state != "apo":
                merged = merge_chains(structure, TCR_BETA_CHAIN, TCR_ALPHA_CHAIN)
                metrics = ca.interface_metrics(
                    merged, binder_chain_hint=TCR_ALPHA_CHAIN,
                    expected_residue=None)
                rec["metrics"] = {k: v for k, v in metrics.items()
                                  if k not in ("chain_map",)}
        except Exception as exc:  # noqa: BLE001
            rec["error"] = f"{type(exc).__name__}: {exc}"[:600]
        rec["seconds"] = round(time.time() - t0, 1)
        with CKPT.open("a") as fh:
            fh.write(json.dumps(rec, default=str) + "\n")
        done[name] = rec
        print(f"[control] {name}: iptm={rec.get('iptm_binder_vs_target')} "
              f"{rec.get('seconds')}s err={rec.get('error', '')}", flush=True)
        time.sleep(MIN_INTERVAL)

    # ---- assemble the gate read-out on the constrained G12D run -----------
    g = done.get("g12d_constrained", {})
    metrics = g.get("metrics") or {}
    row = {
        "name": "7OW6_native_TCR",
        "iptm_g12d": g.get("iptm_binder_vs_target"),
        "iptm_wt": done.get("wt_constrained", {}).get("iptm_binder_vs_target"),
        "iptm_apo": done.get("apo_unconstrained", {}).get("iptm_binder_vs_target"),
        "iptm_decoy": done.get("decoy_constrained", {}).get("iptm_binder_vs_target"),
        "peptide_buried_frac": metrics.get("peptide_buried_frac"),
        "asp6_buried_frac": metrics.get("asp6_buried_frac"),
        "buried_area_peptide": metrics.get("buried_area_peptide"),
        "buried_area_total": metrics.get("buried_area_total"),
        "pos6_contacts": metrics.get("pos6_contacts"),
        "asp6_hbond": metrics.get("asp6_hbond"),
        "asp6_salt_bridge": metrics.get("asp6_salt_bridge"),
        "net_interface_charge": metrics.get("net_interface_charge"),
        "pae_interface_g12d": g.get("complex_ipde"),
    }
    row["specificity_gap"] = design_gates.specificity_gap(
        {k: v for k, v in row.items() if v is not None})
    row["gate_failures"] = design_gates.evaluate(
        {k: v for k, v in row.items() if v is not None})
    pose_ok, pose_reason = pose_filter.predicted_pose_ok(metrics)
    row["pose_filter_ok"] = pose_ok
    row["pose_filter_reason"] = pose_reason

    report = {
        "generated": time.strftime("%Y-%m-%d %H:%M:%S %z"),
        "runs": {k: {kk: vv for kk, vv in v.items()
                     if kk not in ("pair_chains_iptm_scores", "metrics")}
                 for k, v in done.items()},
        "g12d_constrained_metrics": metrics,
        "g12d_unconstrained_metrics": done.get("g12d_unconstrained", {}).get("metrics"),
        "gate_readout": row,
    }
    path = OUT / "7ow6_constrained_control.json"
    path.write_text(json.dumps(report, indent=2, default=str))
    print(json.dumps(row, indent=2, default=str))
    print(f"[control] wrote {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
