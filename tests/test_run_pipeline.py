"""Unit tests for src/run_pipeline.py. No live API calls: the pipeline is driven
either by `DryRunClient` or by a recording stub, and the Boltz-2 score parsing
is checked against the real response shape observed on 2026-11-02.
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest

import design_gates
import run_pipeline as rp

EXAMPLE_CSV = Path(__file__).resolve().parents[1] / "examples" / "scored_designs_example.csv"
TARGET_DIR = Path(__file__).resolve().parents[1] / "data" / "targets"


# --- target handling -------------------------------------------------------

def test_read_fasta_joins_wrapped_lines(tmp_path):
    path = tmp_path / "x.fasta"
    path.write_text(">a|one\nACDE\nFGHI\n>b|two\nKLMN\n")
    assert rp.read_fasta(path) == ["ACDEFGHI", "KLMN"]


def test_load_target_states_gives_the_four_panel_states():
    states = rp.load_target_states(TARGET_DIR)
    assert set(states) == set(rp.STATES)
    assert [p["id"] for p in states["g12d"]] == ["A", "B", "C"]
    assert states["g12d"][2]["sequence"] == "VVVGADGVGK"
    assert states["wt"][2]["sequence"] == "VVVGAGGVGK"
    assert states["decoy"][2]["sequence"] == design_gates.DECOY_PEPTIDE
    assert [p["id"] for p in states["apo"]] == ["A", "B"]      # no peptide chain
    assert all(p["molecule_type"] == "protein"
               for state in states.values() for p in state)


def test_load_target_states_reports_a_missing_file(tmp_path):
    with pytest.raises(FileNotFoundError, match="fetch_target"):
        rp.load_target_states(tmp_path)


def test_binder_contig_keeps_groove_and_peptide_then_builds_the_binder():
    assert rp.binder_contig(70) == "A1-180/0 C1-10/0 70-70"
    assert rp.binder_contig(85, peptide_length=9).endswith("85-85")


def test_hotspots_cover_the_mutated_position():
    assert "C6" in rp.HOTSPOT_RESIDUES
    assert rp.HOTSPOT_RESIDUES == ["C4", "C5", "C6", "C7", "C8"]


# --- sequence filters ------------------------------------------------------

@pytest.mark.parametrize(
    "sequence,expected",
    [
        ("AEKLRSTVYWQNGIDFHM", []),
        ("A" * 5 + "EK", ["more than 4 consecutive identical residues"]),
        ("KKKKRRRRKKKK", None),           # charge, checked loosely below
        ("AANASAA", ["N-glycosylation sequon (N-X-S/T)"]),
        ("ACAAAA", ["odd number of cysteines (unpaired Cys)"]),
    ],
)
def test_sequence_filter_failures(sequence, expected):
    failures = rp.sequence_filter_failures(sequence)
    if expected is None:
        assert any("net charge" in f for f in failures)
    else:
        assert failures == expected


def test_paired_cysteines_pass():
    assert rp.sequence_filter_failures("ACDEFGHIKLMNPQRSTVCY") == []


# --- Boltz-2 score extraction ---------------------------------------------

REAL_SHAPE = {
    "pair_chains_iptm_scores": [
        {
            "0": {"0": 0.963, "1": 0.952, "2": 0.899, "3": 0.105},
            "1": {"0": 0.963, "1": 0.988, "2": 0.856, "3": 0.093},
            "2": {"0": 0.958, "1": 0.936, "2": 0.963, "3": 0.162},
            "3": {"0": 0.147, "1": 0.143, "2": 0.140, "3": 0.72},
        }
    ],
    "iptm_scores": [0.792],
}


def test_binder_interface_iptm_ignores_the_native_interfaces():
    # The global ipTM is 0.79 only because HLA/B2M/peptide fit each other; the
    # binder (chain index 3) barely touches anything.
    assert rp.binder_interface_iptm(REAL_SHAPE, binder_index=3) == pytest.approx(0.147)


def test_binder_interface_iptm_skips_the_diagonal():
    payload = {"pair_chains_iptm_scores": [{"0": {"0": 0.9, "1": 0.2},
                                           "1": {"0": 0.2, "1": 0.99}}]}
    assert rp.binder_interface_iptm(payload, binder_index=1) == pytest.approx(0.2)


def test_binder_interface_iptm_falls_back_to_the_global_score():
    assert rp.binder_interface_iptm({"iptm_scores": [0.42]}, 3) == pytest.approx(0.42)
    assert rp.binder_interface_iptm({}, 3) is None


def test_first_reads_the_leading_sample():
    assert rp._first({"complex_ipde_scores": [4.3, 9.9]}, "complex_ipde_scores") == 4.3
    assert rp._first({}, "complex_ipde_scores") is None


# --- structure metrics -----------------------------------------------------

def _two_chain_mmcif():
    return rp._stub_mmcif(
        [
            {"id": "C", "molecule_type": "protein", "sequence": "VVVGADGVGK"},
            {"id": "D", "molecule_type": "protein", "sequence": "KKDEAAAA"},
        ]
    )


def test_complex_metrics_runs_on_a_parseable_structure():
    metrics = rp.complex_metrics(_two_chain_mmcif())
    assert 0.0 <= metrics["peptide_buried_frac"] <= 1.0
    assert 0.0 <= metrics["asp6_buried_frac"] <= 1.0
    assert metrics["pos6_contacts"] >= 0
    assert "net_interface_charge" in metrics


def test_complex_metrics_requires_the_binder_chain():
    with pytest.raises(ValueError, match="binder chain"):
        rp.complex_metrics(_two_chain_mmcif(), binder_chain="Z")


def test_interface_charge_counts_only_contacting_residues():
    # Two chains 6 A apart in the stub: with a 5 A cutoff nothing is in contact,
    # so the interface charge is zero even though the binder carries +2/-2.
    structure = rp._parse_mmcif(_two_chain_mmcif())
    model = next(iter(structure))
    assert rp._interface_charge(model, "D") == pytest.approx(0.0)


def test_residue_charges_match_the_documented_convention():
    assert rp.RESIDUE_CHARGE["ARG"] == 1.0
    assert rp.RESIDUE_CHARGE["ASP"] == -1.0
    assert rp.RESIDUE_CHARGE["HIS"] == pytest.approx(0.1)


# --- scoring one sequence --------------------------------------------------

class RecordingClient(rp.DryRunClient):
    """DryRunClient that also records the polymer sets it was asked to predict."""

    def __init__(self):
        super().__init__(seed=0)
        self.polymer_sets = []

    def boltz2_predict(self, polymers, **kwargs):
        self.polymer_sets.append([p["sequence"] for p in polymers])
        return super().boltz2_predict(polymers, **kwargs)


def test_score_sequence_runs_all_four_states_with_the_binder_appended():
    client = RecordingClient()
    states = rp.load_target_states(TARGET_DIR)
    binder = "KKDEAAAAHHRR"
    row, detail = rp.score_sequence(client, binder, states, {})
    assert len(client.polymer_sets) == 4
    assert all(polymers[-1] == binder for polymers in client.polymer_sets)
    assert {"VVVGADGVGK", "VVVGAGGVGK", design_gates.DECOY_PEPTIDE} <= {
        s for polymers in client.polymer_sets for s in polymers
    }
    assert len(client.polymer_sets[2]) == 3          # apo: HLA + B2M + binder
    for state in rp.STATES:
        assert row[design_gates.IPTM_COLUMNS[state]] != ""
        assert state in detail["states"]
    assert row["specificity_gap"] == pytest.approx(
        row["iptm_g12d"] - max(row["iptm_wt"], row["iptm_apo"], row["iptm_decoy"])
    )
    assert row["molecule_class"] == "single_chain"


def test_score_sequence_applies_the_gates():
    client = RecordingClient()
    states = rp.load_target_states(TARGET_DIR)
    row, _ = rp.score_sequence(client, "KKDEAAAAHHRR", states, {})
    assert row["gate_failures"] == design_gates.evaluate(row)


def test_score_sequence_survives_an_unparseable_structure(monkeypatch):
    client = RecordingClient()
    monkeypatch.setattr(
        rp, "complex_metrics", lambda *a, **k: (_ for _ in ()).throw(ValueError("bad"))
    )
    states = rp.load_target_states(TARGET_DIR)
    row, detail = rp.score_sequence(client, "KKDEAAAA", states, {})
    assert "ValueError: bad" in detail["metrics_error"]
    assert row["iptm_g12d"] != ""          # the ipTM panel still came through


# --- CSV output ------------------------------------------------------------

def test_csv_columns_match_the_example_exactly():
    with EXAMPLE_CSV.open() as handle:
        header = next(csv.reader(handle))
    assert rp.CSV_COLUMNS == header


def test_write_scored_csv_drops_extra_keys(tmp_path):
    out = tmp_path / "scored.csv"
    rp.write_scored_csv(
        [{"name": "D1", "sequence": "ACDE", "iptm_g12d": 0.9,
          "gate_failures": ["nope"]}],
        out,
    )
    rows = list(csv.DictReader(out.open()))
    assert list(rows[0]) == rp.CSV_COLUMNS
    assert rows[0]["name"] == "D1"
    assert rows[0]["iptm_apo"] == ""


def test_design_name_is_deterministic_and_submission_safe():
    name = rp.design_name("ACDEFG", 3)
    assert name == rp.design_name("ACDEFG", 3)
    assert name.startswith("KRASG12D_NIM_003_")
    assert name.replace("_", "").isalnum()


# --- end to end, dry run ---------------------------------------------------

def test_dry_run_end_to_end_writes_a_scored_csv(tmp_path, monkeypatch):
    monkeypatch.delenv("NVIDIA_API_KEY", raising=False)   # must not be needed
    out = tmp_path / "scored.csv"
    exit_code = rp.main(
        [
            "--dry-run",
            "--backbones", "1",
            "--seqs-per-backbone", "3",
            "--out", str(out),
            "--target-dir", str(TARGET_DIR),
        ]
    )
    assert exit_code == 0
    rows = list(csv.DictReader(out.open()))
    assert rows, "dry run produced no designs"
    assert list(rows[0]) == rp.CSV_COLUMNS
    for row in rows:
        assert float(row["iptm_g12d"]) > float(row["iptm_wt"])
        assert row["sequence"] and row["name"].startswith("KRASG12D_NIM_")
    details = json.loads(out.with_suffix(".details.json").read_text())
    assert len(details) == len(rows)
    assert set(details[0]["states"]) == set(rp.STATES)


def test_dry_run_makes_no_client_and_no_calls(monkeypatch):
    def explode(*args, **kwargs):
        raise AssertionError("--dry-run must not construct a live client")

    monkeypatch.setattr(rp, "NvidiaClient", explode)
    args = rp.build_parser().parse_args(
        ["--dry-run", "--backbones", "1", "--seqs-per-backbone", "1",
         "--out", "/tmp/_dryrun_test.csv", "--target-dir", str(TARGET_DIR)]
    )
    assert rp.run(args) == 0


def test_dry_run_client_shapes_match_the_real_endpoints():
    client = rp.DryRunClient(seed=3)
    backbone = client.rfdiffusion_generate("pdb", "A1-180/0 C1-10/0 70-70")
    assert backbone.count(" CA ") == 70
    designed = client.proteinmpnn_design(backbone, input_pdb_chains=["A"],
                                         num_seq_per_target=2)
    assert len(designed["sequences"]) == 2
    assert all(len(s["sequence"]) == 70 for s in designed["sequences"])
    payload = client.boltz2_predict(
        [{"id": "C", "molecule_type": "protein", "sequence": "VVVGADGVGK"},
         {"id": "D", "molecule_type": "protein", "sequence": "ACDEF"}]
    )
    assert payload["structures"][0]["format"] == "mmcif"
    assert client.calls == ["rfdiffusion", "proteinmpnn", "boltz2"]


def test_live_run_requires_the_target_pdb(tmp_path, monkeypatch):
    monkeypatch.setenv("NVIDIA_API_KEY", "test-key")
    for name in rp.STATE_FASTA.values():
        (tmp_path / name).write_text(
            ">a\nACDE\n>b\nKLMN\n" + (">c\nVVVGADGVGK\n" if "apo" not in name else "")
        )
    args = rp.build_parser().parse_args(
        ["--target-dir", str(tmp_path), "--out", str(tmp_path / "o.csv")]
    )
    with pytest.raises(FileNotFoundError, match="fetch_target"):
        rp.run(args)
