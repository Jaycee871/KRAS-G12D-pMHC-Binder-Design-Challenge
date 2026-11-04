#!/usr/bin/env python3
"""Pose accuracy of the constrained 7OW6 control: RMSD to the crystal.

Superpose the predicted complex on the crystal using the pMHC only (HLA heavy
chain + peptide C-alphas, which is the frame the binder has to be correct in),
then report the C-alpha RMSD of the TCR chains *without* re-fitting them. That
is the number that says whether the predicted binder sits where the real one
does, independently of any buried-area definition.
"""
from __future__ import annotations
import io, json, sys
from pathlib import Path
import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))
import complex_analysis as ca  # noqa: E402
from Bio.PDB import MMCIFParser  # noqa: E402

def model_of(path: Path):
    return next(iter(MMCIFParser(QUIET=True).get_structure("x", io.StringIO(path.read_text()))))

def ca_list(model, chain_id):
    out = []
    for r in model[chain_id]:
        if ca.is_polymer_residue(r) and "CA" in r:
            out.append(r["CA"].coord)
    return np.array(out, dtype=float)

def pairs(a, b):
    n = min(len(a), len(b))
    return a[:n], b[:n]

def kabsch(P, Q):
    """Rotation+translation taking P onto Q."""
    pc, qc = P.mean(0), Q.mean(0)
    H = (P - pc).T @ (Q - qc)
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    D = np.diag([1, 1, d])
    R = Vt.T @ D @ U.T
    return R, qc - R @ pc

def main() -> int:
    xtal = model_of(REPO / "data" / "validation" / "7ow6.cif")
    out = {}
    for name in ("g12d_constrained", "g12d_unconstrained"):
        pred = model_of(REPO / "results" / "validation" / f"7ow6_control_{name}.cif")
        # frame: HLA heavy chain + peptide
        fit_p, fit_x = [], []
        for cid in ("A", "C"):
            p, x = pairs(ca_list(pred, cid), ca_list(xtal, cid))
            fit_p.append(p); fit_x.append(x)
        P = np.vstack(fit_p); X = np.vstack(fit_x)
        R, t = kabsch(P, X)
        frame_rmsd = float(np.sqrt((((P @ R.T + t) - X) ** 2).sum(1).mean()))
        res = {"frame_residues": len(P), "frame_ca_rmsd": round(frame_rmsd, 2)}
        for cid, label in (("D", "tcr_alpha"), ("E", "tcr_beta")):
            p, x = pairs(ca_list(pred, cid), ca_list(xtal, cid))
            moved = p @ R.T + t
            res[f"{label}_n"] = len(p)
            res[f"{label}_ca_rmsd_in_pmhc_frame"] = round(
                float(np.sqrt(((moved - x) ** 2).sum(1).mean())), 2)
            # variable-domain only (first 115 residues) - the part that binds
            k = min(115, len(p))
            res[f"{label}_v_domain_ca_rmsd"] = round(
                float(np.sqrt(((moved[:k] - x[:k]) ** 2).sum(1).mean())), 2)
        out[name] = res
        print(name, json.dumps(res))
    path = REPO / "results" / "validation" / "7ow6_control_rmsd.json"
    path.write_text(json.dumps(out, indent=2))
    print("wrote", path)
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
