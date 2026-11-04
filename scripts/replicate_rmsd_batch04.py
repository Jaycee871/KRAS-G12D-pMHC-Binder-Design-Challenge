#!/usr/bin/env python3
"""Replicate pose stability for batch-04 gate survivors: binder C-alpha RMSD.

The 7OW6 post-mortem (results/validation/CONTROL_7OW6.md) left exactly one
ranking axis standing: whether a design's *pose* reproduces across identical
constrained calls. ipTM gaps and paired MM-GBSA ddG are retired.

Method, identical in spirit to scripts/control_7ow6_rmsd.py: superpose each
replicate on the first one using the **pMHC frame only** (HLA heavy chain A +
peptide C C-alphas), then measure the binder (chain D) C-alpha RMSD *without*
re-fitting it. That isolates "does the binder dock in the same place" from
"did the binder fold the same way".

Reported per survivor: all pairwise binder RMSDs across the ensemble (the
original scored pose plus every replicate), with mean/max/SD.
"""
from __future__ import annotations
import io, itertools, json, sys
from pathlib import Path
import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))
import complex_analysis as ca  # noqa: E402
from Bio.PDB import MMCIFParser  # noqa: E402

OUT = REPO / "results" / "batch-04"
HLA_CHAIN, PEPTIDE_CHAIN, BINDER_CHAIN = "A", "C", "D"


def model_of(path: Path):
    return next(iter(MMCIFParser(QUIET=True).get_structure(
        "x", io.StringIO(path.read_text()))))


def ca_list(model, chain_id):
    return np.array([r["CA"].coord for r in model[chain_id]
                     if ca.is_polymer_residue(r) and "CA" in r], dtype=float)


def pairs(a, b):
    n = min(len(a), len(b))
    return a[:n], b[:n]


def kabsch(P, Q):
    pc, qc = P.mean(0), Q.mean(0)
    H = (P - pc).T @ (Q - qc)
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    return R, qc - R @ pc


def binder_rmsd(ref, mob):
    """Binder C-alpha RMSD of `mob` vs `ref` in the pMHC frame."""
    fit_m, fit_r = [], []
    for cid in (HLA_CHAIN, PEPTIDE_CHAIN):
        m, r = pairs(ca_list(mob, cid), ca_list(ref, cid))
        fit_m.append(m); fit_r.append(r)
    M, R_ = np.vstack(fit_m), np.vstack(fit_r)
    rot, t = kabsch(M, R_)
    frame = float(np.sqrt((((M @ rot.T + t) - R_) ** 2).sum(1).mean()))
    bm, br = pairs(ca_list(mob, BINDER_CHAIN), ca_list(ref, BINDER_CHAIN))
    moved = bm @ rot.T + t
    return (float(np.sqrt(((moved - br) ** 2).sum(1).mean())), len(bm),
            round(frame, 2))


def main() -> int:
    survivors = json.loads((OUT / "gate_results.json").read_text())["survivors"]
    records = []
    for s in survivors:
        stem = f"{s['backbone']}_{s['variant']}"
        members = [("original", OUT / "scored_structures" / f"{stem}_g12d.cif")]
        for path in sorted((OUT / "replicate_structures").glob(f"{stem}_r*.cif")):
            members.append((path.stem.rsplit("_", 1)[-1], path))
        members = [(label, p) for label, p in members if p.exists()]
        models = {label: model_of(p) for label, p in members}
        rows = []
        for (la, ma), (lb, mb) in itertools.combinations(models.items(), 2):
            rmsd, n, frame = binder_rmsd(ma, mb)
            rows.append({"pair": f"{la}~{lb}", "binder_ca_rmsd": round(rmsd, 2),
                         "n_ca": n, "frame_ca_rmsd": frame})
        values = [r["binder_ca_rmsd"] for r in rows]
        records.append({
            "backbone": s["backbone"], "variant": s["variant"],
            "ensemble": [l for l, _ in members],
            "pairwise": rows,
            "binder_ca_rmsd_mean": round(float(np.mean(values)), 2) if values else None,
            "binder_ca_rmsd_max": round(float(np.max(values)), 2) if values else None,
            "binder_ca_rmsd_sd": round(float(np.std(values)), 2) if values else None,
        })
        print(f"{s['backbone']} ({s['variant']}): "
              f"mean {records[-1]['binder_ca_rmsd_mean']} A, "
              f"max {records[-1]['binder_ca_rmsd_max']} A over "
              f"{len(members)} poses")
        for r in rows:
            print(f"   {r['pair']:>18}  binder {r['binder_ca_rmsd']:>6} A "
                  f"(frame {r['frame_ca_rmsd']} A, {r['n_ca']} CA)")
    (OUT / "replicate_rmsd.json").write_text(json.dumps(records, indent=2))
    print("wrote", OUT / "replicate_rmsd.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
