"""Unit tests for the batch runner's checkpointing (batch-01 problem 1).

Batch-01 scored 50 of 56 candidates and then lost all of them, because
`designs.csv` was only written after the last future resolved. These tests pin
the replacement: every finished candidate is on disk immediately, a resumed run
skips what is already there, and the CSV can be rebuilt from the checkpoint
without touching the network.

Nothing here makes an API call.
"""
from __future__ import annotations

import argparse
import csv
import json

import pytest

import run_batch


def make_result(sequence, gap=0.5, strategy="S1-core-short", structure="data_x\n",
                **row_extra):
    row = {
        "sequence": sequence,
        "strategy": strategy,
        "backbone": f"{strategy}_bb00",
        "molecule_class": "single_chain",
        "iptm_g12d": 0.80,
        "iptm_wt": 0.80 - gap,
        "iptm_apo": 0.10,
        "iptm_decoy": 0.10,
        "specificity_gap": gap,
        "peptide_buried_frac": 0.55,
        "asp6_buried_frac": 0.70,
        "pos6_contacts": 9,
        "net_interface_charge": 1.0,
        "analysis_error": "",
        "gate_failures": "",
        "n_gate_failures": 0,
    }
    row.update(row_extra)
    return {
        "row": row,
        "detail": {"sequence": sequence, "states": {}, "strategy": strategy},
        "structure": structure,
    }


def fake_args(out_dir, **kwargs):
    defaults = dict(
        top_structures=5, recycling_steps=3, sampling_steps=50,
        diffusion_steps=15, seqs_per_backbone=8, seqs_scored_per_backbone=4,
        workers=1, resume=False, out_dir=str(out_dir),
    )
    defaults.update(kwargs)
    return argparse.Namespace(**defaults)


# --- checkpoint round trip --------------------------------------------------

def test_a_scored_candidate_is_on_disk_immediately(tmp_path):
    run_batch.append_checkpoint(tmp_path, make_result("AAAK"))
    path = run_batch.checkpoint_path(tmp_path)
    assert path.exists()
    record = json.loads(path.read_text().splitlines()[0])
    assert record["sequence"] == "AAAK"
    assert record["row"]["specificity_gap"] == 0.5
    assert (tmp_path / record["structure_file"]).read_text() == "data_x\n"


def test_checkpoint_appends_one_line_per_candidate(tmp_path):
    for index in range(3):
        run_batch.append_checkpoint(tmp_path, make_result(f"SEQ{index}"))
    lines = run_batch.checkpoint_path(tmp_path).read_text().splitlines()
    assert len(lines) == 3
    assert set(run_batch.load_checkpoint(tmp_path)) == {"SEQ0", "SEQ1", "SEQ2"}


def test_a_rescored_candidate_replaces_the_earlier_line(tmp_path):
    run_batch.append_checkpoint(tmp_path, make_result("SEQ", gap=0.1))
    run_batch.append_checkpoint(tmp_path, make_result("SEQ", gap=0.9))
    done = run_batch.load_checkpoint(tmp_path)
    assert len(done) == 1
    assert done["SEQ"]["row"]["specificity_gap"] == 0.9


def test_a_truncated_last_line_does_not_lose_the_rest(tmp_path):
    """An interrupt mid-write must cost one candidate, not the whole file."""
    run_batch.append_checkpoint(tmp_path, make_result("GOOD"))
    with run_batch.checkpoint_path(tmp_path).open("a") as handle:
        handle.write('{"sequence": "TRUNC", "row": {"seq')
    done = run_batch.load_checkpoint(tmp_path)
    assert set(done) == {"GOOD"}


def test_results_from_checkpoint_restores_the_structures(tmp_path):
    run_batch.append_checkpoint(tmp_path, make_result("SEQ", structure="cif!"))
    results = run_batch.results_from_checkpoint(tmp_path)
    assert len(results) == 1
    assert results[0]["structure"] == "cif!"
    assert results[0]["row"]["sequence"] == "SEQ"


def test_missing_checkpoint_is_an_empty_dict(tmp_path):
    assert run_batch.load_checkpoint(tmp_path) == {}
    assert run_batch.results_from_checkpoint(tmp_path) == []


def test_structure_names_are_stable_and_sequence_derived():
    assert run_batch.structure_name("ABC") == run_batch.structure_name("ABC")
    assert run_batch.structure_name("ABC") != run_batch.structure_name("ABD")


# --- rebuilding the batch artifacts ----------------------------------------

def test_finalize_writes_the_csv_sorted_by_specificity_gap(tmp_path):
    results = [make_result("LOW", gap=0.1), make_result("HIGH", gap=0.9)]
    stats = run_batch.finalize(tmp_path, results, fake_args(tmp_path), 0.0,
                               backbones=2, candidates=2,
                               strategies=run_batch.STRATEGIES)
    rows = list(csv.DictReader((tmp_path / "designs.csv").open()))
    assert [r["sequence"] for r in rows] == ["HIGH", "LOW"]
    assert rows[0]["name"].startswith("KRASG12D_NIM_001_")
    assert stats["scored"] == 2
    assert stats["passing"] == 2
    assert json.loads((tmp_path / "run_stats.json").read_text())["scored"] == 2


def test_finalize_counts_analysis_errors_as_failures(tmp_path):
    broken = make_result("BROKEN", analysis_error="ChainMappingError: boom",
                         gate_failures="interface analysis failed",
                         n_gate_failures=1)
    stats = run_batch.finalize(tmp_path, [broken], fake_args(tmp_path), 0.0,
                               backbones=0, candidates=1,
                               strategies=run_batch.STRATEGIES)
    assert stats["analysis_errors"] == 1
    assert stats["passing"] == 0


def test_finalize_writes_only_the_top_structures(tmp_path):
    results = [make_result(f"SEQ{i}", gap=i / 10) for i in range(6)]
    run_batch.finalize(tmp_path, results, fake_args(tmp_path, top_structures=2),
                       0.0, backbones=0, candidates=6,
                       strategies=run_batch.STRATEGIES)
    written = sorted(p.name for p in (tmp_path / "structures").glob("*.cif"))
    assert len(written) == 2
    assert written[0].startswith("01_")


def test_rebuild_regenerates_the_csv_without_any_api_call(tmp_path):
    for index in range(4):
        run_batch.append_checkpoint(tmp_path, make_result(f"SEQ{index}",
                                                          gap=index / 10))
    exit_code = run_batch.main(["--rebuild", "--out-dir", str(tmp_path)])
    assert exit_code == 0
    rows = list(csv.DictReader((tmp_path / "designs.csv").open()))
    assert len(rows) == 4
    assert rows[0]["sequence"] == "SEQ3"      # highest gap first


# --- resuming ---------------------------------------------------------------

def test_score_all_skips_candidates_already_checkpointed(tmp_path, monkeypatch):
    run_batch.append_checkpoint(tmp_path, make_result("DONE"))
    scored = []

    def fake_score(candidate, target_states, boltz_kwargs, args):
        scored.append(candidate["sequence"])
        return make_result(candidate["sequence"])

    monkeypatch.setattr(run_batch, "score_candidate", fake_score)
    candidates = [{"sequence": "DONE", "strategy": "S1-core-short",
                   "backbone": "S1-core-short_bb00"},
                  {"sequence": "NEW", "strategy": "S1-core-short",
                   "backbone": "S1-core-short_bb01"}]
    run_batch.score_all(candidates, tmp_path, {}, {},
                        fake_args(tmp_path, resume=True), 0.0, backbones=0,
                        strategies=run_batch.STRATEGIES)
    assert scored == ["NEW"]
    assert set(run_batch.load_checkpoint(tmp_path)) == {"DONE", "NEW"}


def test_score_all_without_resume_rescores_everything(tmp_path, monkeypatch):
    run_batch.append_checkpoint(tmp_path, make_result("DONE"))
    scored = []

    def fake_score(candidate, target_states, boltz_kwargs, args):
        scored.append(candidate["sequence"])
        return make_result(candidate["sequence"])

    monkeypatch.setattr(run_batch, "score_candidate", fake_score)
    run_batch.score_all([{"sequence": "DONE", "strategy": "S1-core-short",
                          "backbone": "bb"}], tmp_path, {}, {},
                        fake_args(tmp_path, resume=False), 0.0, backbones=0,
                        strategies=run_batch.STRATEGIES)
    assert scored == ["DONE"]


def test_a_dropped_candidate_does_not_reach_the_checkpoint(tmp_path, monkeypatch):
    monkeypatch.setattr(run_batch, "score_candidate",
                        lambda *a, **k: None)      # simulates an API failure
    run_batch.score_all([{"sequence": "LOST", "strategy": "S1-core-short",
                          "backbone": "bb"}], tmp_path, {}, {},
                        fake_args(tmp_path), 0.0, backbones=0,
                        strategies=run_batch.STRATEGIES)
    assert run_batch.load_checkpoint(tmp_path) == {}


def test_analysis_columns_are_in_the_csv_schema():
    for column in ("buried_area_peptide", "buried_area_hla", "contacts_peptide",
                   "contacts_hla", "analysis_error"):
        assert column in run_batch.BATCH_CSV_COLUMNS
