#!/usr/bin/env python3
"""Replication confirmation + off-target selectivity panel.

Two steps, one script, both re-using code that already exists:

**Step 1 - confirm replication.** The sampling-settings study
(`results/sampling-study/SUMMARY.md`) collapsed the replicate pose spread of
three designs under design-specific Boltz-2 settings, but at only n = 3 draws
per cell. This step takes each of them to n = 10 on the **G12D** state under
its own stabilising setting, re-using `scripts/run_batch04.score_candidate`
verbatim (same polymers, same validated flat pocket constraint, same analysis,
same Asp6 read-out gate) with the sampling fields patched exactly as
`scripts/sampling_study.patched_boltz2` does. The 3 draws already in
`results/sampling-study/checkpoint.jsonl` are re-used, not re-run.

**Step 2 - off-target panel.** For every design confirmed in step 1, the same
binder is scored against the full four-state panel of
`results/validation/CONTROL_7OW6.md` under its own stabilising setting:

    g12d   VVVGADGVGK on HLA-A*11:01   (from step 1, n = 10)
    wt     VVVGAGGVGK                  constrained, as g12d
    decoy  AVFDRKSDAK                  constrained, as g12d
    apo    peptide-free HLA + B2M      unconstrained - there is no chain C to
                                       point a pocket constraint at, so the
                                       constraint is dropped, which is what
                                       `scripts/control_7ow6_constrained.py`
                                       did for this arm too.

How selectivity is judged here
------------------------------
On **measured geometry** and on **reproducibility across draws** only. The ipTM
`specificity_gap` and the paired same-pose MM-GBSA ddG both failed the 7OW6
positive control and are retired: neither is computed, ranked on or revived
anywhere in this script. ipTM is recorded as a pose-quality read-out.

The off-target arms cannot go through `complex_analysis.interface_metrics`
unchanged - WT position 6 is a glycine with no side chain, and apo has no
peptide chain at all, both of which that function (correctly) refuses. So the
arms are measured with the same primitives that function uses
(`_atom_sasa`, `count_contacts`, `interface_charge`, `asp6_hbond_partners`
logic) generalised over the peptide: side-chain contacts at position 6 for any
residue, and donor approach to **every** peptide side-chain acceptor, so a
binder that has merely learned to like a carboxylate shows up on the decoy's
Asp4/Asp8.

Checkpointed per draw with an fsync; re-running resumes and re-spends nothing.
No design is generated: Track-1 design slots consumed = 0.
"""
from __future__ import annotations

import argparse, json, os, sys, time
from pathlib import Path
from typing import Any, Dict, List, Optional

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(REPO / "scripts"))

import complex_analysis as ca  # noqa: E402
import run_batch04 as batch  # noqa: E402
import run_pipeline  # noqa: E402
import sampling_study as study  # noqa: E402
from nvidia_client import BOLTZ2_PATH, NvidiaClient  # noqa: E402

OUT = REPO / "results" / "selectivity-panel"
STUDY = REPO / "results" / "sampling-study"

# design -> the condition that stabilised it in the sampling study
SETTINGS = {
    "bb50-g1": "recyc8",
    "bb25": "recyc8",
    "bb135-g1": "combo",
}
READOUT_VERSION = "carboxylate-oxygen-and-burial-v1"
ACCEPTOR_ATOMS = {"ASP": ("OD1", "OD2"), "GLU": ("OE1", "OE2")}


def boltz_kwargs(condition: str) -> Dict[str, Any]:
    recyc, samp, nopot = study.CONDITIONS[condition]
    return {"recycling_steps": recyc, "sampling_steps": samp,
            "diffusion_samples": 1, "without_potentials": nopot}


def call_boltz(client: NvidiaClient, polymers, condition: str,
               constrained: bool) -> Dict[str, Any]:
    """One Boltz-2 call, identical to batch-04's except for the two fields."""
    body = {"polymers": [dict(p) for p in polymers], "step_scale": 1.638,
            **boltz_kwargs(condition)}
    if constrained:
        body["constraints"] = batch.POCKET
    return client.post(BOLTZ2_PATH, body)


# --- generic panel-arm measurement ------------------------------------------

def nearest_carboxylate_donor(residue, search, cutoff=6.0):
    """Only ASP OD1/OD2 or GLU OE1/OE2 acceptors; decide on unrounded distance."""
    names = ACCEPTOR_ATOMS.get(residue.get_resname().strip().upper(), ())
    donor_atoms = {**ca.CHARGED_DONORS, **ca.HYDROXYL_DONORS}
    best = None
    for atom in residue:
        if atom.get_id() not in names:
            continue
        for neighbour in search.search(atom.coord, cutoff):
            donor_res = neighbour.get_parent()
            rname = donor_res.get_resname().strip().upper()
            aname = neighbour.get_id().strip().upper()
            if rname not in batch.GATE_DONOR_RESIDUES or aname not in donor_atoms.get(rname, ()):
                continue
            distance = float(((neighbour.coord - atom.coord) ** 2).sum() ** 0.5)
            if best is None or distance < best["distance"]:
                best = {"donor_residue": rname, "donor_resseq": donor_res.id[1],
                        "donor_atom": aname, "acceptor_atom": atom.get_id(),
                        "distance": distance}
    return best


def panel_metrics(text: str, binder_sequence: str, arm: str) -> Dict[str, Any]:
    """Interface geometry of one panel arm (works for apo and for Gly6)."""
    from Bio.PDB import NeighborSearch

    structure = ca.parse_structure(text, fmt="mmcif")
    model = ca._clean_model(next(iter(structure)))
    chain_map = ca.map_complex_chains(
        model, binder_sequence=binder_sequence,
        binder_chain_hint=run_pipeline.BINDER_CHAIN,
        require_peptide=(arm != "apo"))
    binder_chain = chain_map.binder
    peptide_chain = chain_map.peptide

    bound = ca._atom_sasa(model)
    free_model = ca._clean_model(next(iter(ca.parse_structure(text, fmt="mmcif"))))
    free_model.detach_child(binder_chain)
    free = ca._atom_sasa(free_model)

    buried_total, buried_by_chain = 0.0, {}
    for key, free_value in free.items():
        delta = free_value - bound.get(key, 0.0)
        if delta <= 0:
            continue
        buried_total += delta
        buried_by_chain[key[0]] = buried_by_chain.get(key[0], 0.0) + delta

    def f(value):
        return None if value is None else float(value)

    out: Dict[str, Any] = {
        "arm": arm,
        "binder_chain": binder_chain,
        "peptide_chain": peptide_chain,
        "buried_area_total": f(round(buried_total, 1)),
        "buried_area_by_chain": {k: f(round(v, 1)) for k, v in buried_by_chain.items()},
        "net_interface_charge": f(round(ca.interface_charge(model, binder_chain), 1)),
        "contacts_by_chain": ca.contacts_by_chain(model, binder_chain),
    }
    if peptide_chain is None:
        out.update({"readout_version": READOUT_VERSION, "asp6_buried_frac": None,
                    "asp6_donor_distance_pass": False, "asp6_gate_pass": False,
                    "peptide_buried_frac": None, "pos6_residue": None,
                    "pos6_contacts": None, "pos6_sidechain_contacts": None,
                    "pos6_min_donor_distance": None, "donor_contacts": []})
        return out

    buried_peptide = buried_by_chain.get(peptide_chain, 0.0)
    out["buried_area_peptide"] = f(round(buried_peptide, 1))
    out["peptide_buried_frac"] = f(
        round(buried_peptide / buried_total, 3) if buried_total > 0 else 0.0)

    residues = [r for r in model[peptide_chain] if ca.is_polymer_residue(r)]
    pos6 = residues[ca.MUTATED_PEPTIDE_POSITION - 1]
    out["pos6_residue"] = pos6.get_resname().strip().upper()
    out["pos6_contacts"] = int(ca.count_contacts(model, binder_chain, pos6))

    binder_atoms = ca.heavy_atoms(model[binder_chain])
    search = NeighborSearch(binder_atoms)

    # side-chain-only contacts at position 6 (backbone contacts do not
    # discriminate G12D from WT - the backbone is identical in both)
    sc_names = set(ca.sidechain_atom_names(pos6))
    sc_atoms = [a for a in pos6 if a.get_id() in sc_names and a.element != "H"]
    hits = set()
    for atom in sc_atoms:
        for neighbour in search.search(atom.coord, ca.CONTACT_CUTOFF):
            hits.add(id(neighbour))
    out["pos6_sidechain_contacts"] = len(hits) if sc_atoms else 0

    # Carbon contacts remain in the contact census, never in the readout gate.
    best = nearest_carboxylate_donor(pos6, search)
    out["readout_version"] = READOUT_VERSION
    out["pos6_min_donor_distance"] = round(best["distance"], 2) if best else None
    out["pos6_donor_partner"] = (
        f"{best['donor_residue']}{best['donor_resseq']}:{best['donor_atom']}"
        f"->{best['acceptor_atom']}@{best['distance']:.2f}A" if best else "")
    out["asp6_buried_frac"] = None
    if out["pos6_residue"] == "ASP":
        free_sc = sum(free.get((peptide_chain, pos6.id[1], n), 0.0) for n in sc_names)
        bound_sc = sum(bound.get((peptide_chain, pos6.id[1], n), 0.0) for n in sc_names)
        if free_sc <= 0:
            raise ca.StructureAnalysisError("Asp6 has no accessible side-chain area in binder-free reference")
        out["asp6_buried_frac"] = (free_sc - min(bound_sc, free_sc)) / free_sc
    out["asp6_donor_distance_pass"] = bool(
        best and best["distance"] <= batch.GATE_DONOR_DISTANCE
        and out["pos6_residue"] == "ASP")
    out["asp6_gate_pass"] = bool(out["asp6_donor_distance_pass"]
                                 and out["asp6_buried_frac"] is not None
                                 and out["asp6_buried_frac"] > 0)

    # donor approach to *every* acidic peptide side chain - carboxylate
    # promiscuity shows up here (the decoy carries Asp4 and Asp8)
    donor_contacts = []
    for index, residue in enumerate(residues, start=1):
        closest = nearest_carboxylate_donor(residue, search, cutoff=4.0)
        if closest:
            donor_contacts.append({
                "peptide_position": index,
                "peptide_residue": residue.get_resname().strip().upper(),
                "donor": f"{closest['donor_residue']}{closest['donor_resseq']}:{closest['donor_atom']}",
                "acceptor_atom": closest["acceptor_atom"],
                "distance": round(closest["distance"], 2),
            })
    out["donor_contacts"] = donor_contacts
    return out


# --- checkpointing ----------------------------------------------------------

def checkpoint_path() -> Path:
    return OUT / "checkpoint.jsonl"


def append_checkpoint(record: Dict[str, Any]) -> None:
    with checkpoint_path().open("a") as handle:
        handle.write(json.dumps(record, default=str) + "\n")
        handle.flush()
        os.fsync(handle.fileno())


def done_keys() -> set:
    path = checkpoint_path()
    if not path.exists():
        return set()
    out = set()
    for line in path.read_text().splitlines():
        if line.strip():
            r = json.loads(line)
            out.add((r["design"], r["arm"], r["draw"]))
    return out


def study_draws(design: str, condition: str) -> int:
    """G12D draws already paid for in the sampling study."""
    rows = [json.loads(l) for l in
            (STUDY / "checkpoint.jsonl").read_text().splitlines() if l.strip()]
    return sum(1 for r in rows
               if r["design"] == design and r["condition"] == condition)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--step", choices=["1", "2"], required=True)
    parser.add_argument("--designs", default="bb50-g1,bb25,bb135-g1")
    parser.add_argument("--target-draws", type=int, default=10,
                        help="step 1: total G12D draws per design, counting "
                             "the sampling-study ones")
    parser.add_argument("--arm-draws", type=int, default=3,
                        help="step 2: draws per off-target arm")
    parser.add_argument("--arms", default="wt,apo,decoy")
    parser.add_argument("--max-boltz-calls", type=int, default=60)
    parser.add_argument("--timeout", type=float, default=600.0)
    args = parser.parse_args(argv)

    (OUT / "structures").mkdir(parents=True, exist_ok=True)
    subjects = {s["key"]: s for s in study.load_subjects()}
    designs = [d.strip() for d in args.designs.split(",") if d.strip()]
    states = run_pipeline.load_target_states(REPO / "data" / "targets")
    client = NvidiaClient(timeout=args.timeout)
    already = done_keys()
    calls, started = 0, time.time()

    jobs: List[Dict[str, Any]] = []
    if args.step == "1":
        for design in designs:
            condition = SETTINGS[design]
            have = study_draws(design, condition)
            for draw in range(have, args.target_draws):
                jobs.append({"design": design, "arm": "g12d", "draw": draw,
                             "condition": condition})
    else:
        for design in designs:
            condition = SETTINGS[design]
            for arm in [a.strip() for a in args.arms.split(",") if a.strip()]:
                for draw in range(args.arm_draws):
                    jobs.append({"design": design, "arm": arm, "draw": draw,
                                 "condition": condition})

    for job in jobs:
        key = (job["design"], job["arm"], job["draw"])
        if key in already:
            print(f"[skip] {key} already done")
            continue
        if calls >= args.max_boltz_calls:
            print(f"[stop] call cap {args.max_boltz_calls} reached")
            break
        subject = subjects[job["design"]]
        sequence = subject["sequence"]
        arm, condition = job["arm"], job["condition"]
        constrained = arm != "apo"
        stem = f"{job['design']}_{arm}_r{job['draw']}"
        if arm == "g12d":
            # Step 1 goes through batch-04's `score_candidate` verbatim - the
            # same call the sampling study made - so the extra draws are
            # directly comparable with the three already on disk.
            call_args = argparse.Namespace(**{
                k: v for k, v in boltz_kwargs(condition).items()
                if k in ("recycling_steps", "sampling_steps")},
                without_potentials=boltz_kwargs(condition)["without_potentials"])
            candidate = {k: subject.get(k, "") for k in
                         ("sequence", "strategy", "backbone", "variant",
                          "binder_length", "mpnn_score", "donor_bias_score",
                          "design_pos6_contacts", "design_pos6_distance")}
            try:
                result = batch.score_candidate(client, candidate, states, call_args)
            except Exception as exc:  # noqa: BLE001
                print(f"[fail] {key}: {type(exc).__name__}: {exc}", file=sys.stderr)
                calls += 1
                continue
            calls += 1
            text = result["structure"]
            row = dict(result["row"])
            row["iptm"] = row.get("iptm_g12d")
            row["pos6_min_donor_distance"] = row.get("asp6_min_donor_distance")
            row["asp6_gate_pass"] = row.get("gate_pass")
            (OUT / "structures" / f"{stem}.cif").write_text(text)
        else:
            polymers = list(states[arm]) + [
                {"id": run_pipeline.BINDER_CHAIN, "molecule_type": "protein",
                 "sequence": sequence}]
            try:
                payload = call_boltz(client, polymers, condition, constrained)
            except Exception as exc:  # noqa: BLE001
                print(f"[fail] {key}: {type(exc).__name__}: {exc}", file=sys.stderr)
                calls += 1
                continue
            calls += 1
            text = payload["structures"][0]["structure"]
            (OUT / "structures" / f"{stem}.cif").write_text(text)
            row = {
                "iptm": run_pipeline._round(
                    run_pipeline.binder_interface_iptm(payload, len(polymers) - 1), 4),
                "pde_interface": run_pipeline._round(
                    run_pipeline._first(payload, "complex_ipde_scores"), 2),
            }
            try:
                row.update(panel_metrics(text, sequence, arm))
                row["analysis_error"] = ""
            except Exception as exc:  # noqa: BLE001
                row["analysis_error"] = f"{type(exc).__name__}: {exc}"
        append_checkpoint({
            "design": job["design"], "backbone": subject["backbone"],
            "variant": subject["variant"], "arm": arm, "draw": job["draw"],
            "condition": condition, "boltz_kwargs": boltz_kwargs(condition),
            "constraint": batch.POCKET if constrained else None,
            "constrained": constrained, "sequence": sequence,
            "row": row, "structure_file": f"structures/{stem}.cif",
            "scored_at": time.strftime("%Y-%m-%dT%H:%M:%S%z")})
        print(f"[boltz2 {calls}] {stem} gate={row.get('asp6_gate_pass')} "
              f"donor={row.get('pos6_min_donor_distance')} "
              f"pos6_sc={row.get('pos6_sidechain_contacts')} "
              f"pep_bf={row.get('peptide_buried_frac')} "
              f"charge={row.get('net_interface_charge')} "
              f"iptm={row.get('iptm')} ({time.time() - started:.0f}s)", flush=True)

    print(f"[selectivity-panel] step {args.step}: {calls} Boltz-2 calls in "
          f"{time.time() - started:.0f}s")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
