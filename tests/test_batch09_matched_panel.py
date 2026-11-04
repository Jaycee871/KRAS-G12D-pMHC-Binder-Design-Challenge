"""Pre-flight tests for the batch-09 matched groove counter-screen."""
from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

import batch09_matched_panel as b9
import run_pipeline


def states():
    return run_pipeline.load_target_states(ROOT / "data" / "targets")


def test_groove_constraint_uses_invariant_hla_positions_not_peptide():
    pocket = b9.GROOVE_POCKET
    assert len(pocket) == 1
    assert pocket[0]["binder"] == run_pipeline.BINDER_CHAIN
    assert pocket[0]["max_distance"] == 6.0
    assert pocket[0]["contacts"] == [
        {"id": run_pipeline.HLA_CHAIN, "residue_index": 70},
        {"id": run_pipeline.HLA_CHAIN, "residue_index": 114},
    ]
    assert all(c["id"] != run_pipeline.PEPTIDE_CHAIN for c in pocket[0]["contacts"])


def test_hla_anchor_identity_is_q70_r114_in_every_arm():
    s = states()
    assert b9.target_anchor_residues(s) == {
        "g12d": "QR",
        "wt": "QR",
        "apo": "QR",
        "decoy": "QR",
    }
    b9.validate_matched_target(s)


def test_same_constraint_is_applied_to_apo_and_peptide_arms():
    s = states()
    design = b9.selected_candidates(["bb54-g2"])[0]
    constraints = []
    for arm in b9.ARMS:
        polymers = list(s[arm]) + [{
            "id": run_pipeline.BINDER_CHAIN,
            "molecule_type": "protein",
            "sequence": design["sequence"],
        }]
        body = b9.build_body(polymers, draws=3)
        constraints.append(body["constraints"])
        assert body["diffusion_samples"] == 3
        assert body["recycling_steps"] == 3
        assert body["sampling_steps"] == 50
        assert body["without_potentials"] is True
    assert all(c == b9.GROOVE_POCKET for c in constraints)


def test_candidate_panel_is_fixed_unique_and_contains_expected_mechanisms():
    rows = b9.selected_candidates(b9.DEFAULT_DESIGNS)
    assert [r["key"] for r in rows] == list(b9.DEFAULT_DESIGNS)
    assert len(rows) == 5
    assert len({r["sequence"] for r in rows}) == len(rows)
    assert {"bb54-g2", "bb135-g1", "bb173-g1", "bb90-g2", "b8-apex46-g2dw"} == {
        r["key"] for r in rows
    }


def test_dry_run_request_count_is_designs_times_four_arms(capsys):
    s = states()
    rows = b9.selected_candidates(["bb54-g2", "bb135-g1"])
    assert b9.dry_run(rows, s, draws=3) == 0
    out = capsys.readouterr().out
    assert '"n_model_calls": 8' in out
    assert '"apo"' in out
    assert '"residue_index": 70' in out
    assert '"residue_index": 114' in out
