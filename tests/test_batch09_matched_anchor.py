from __future__ import annotations

import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "batch09_matched_anchor.py"

spec = importlib.util.spec_from_file_location("batch09_matched_anchor", SCRIPT)
b9 = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(b9)


def test_matched_anchor_is_fixed_to_same_hla_region():
    constraint = b9.matched_constraint()
    assert constraint == [{
        "binder": "D",
        "contacts": [
            {"id": "A", "residue_index": 69},
            {"id": "A", "residue_index": 155},
        ],
        "max_distance": 10.0,
    }]


def test_all_four_arms_use_same_constraint_body():
    class FakeClient:
        def __init__(self):
            self.body = None

        def post(self, path, body):
            self.body = body
            return {"structures": []}

    fake = FakeClient()
    polymers = [
        {"id": "A", "molecule_type": "protein", "sequence": "A" * 276},
        {"id": "B", "molecule_type": "protein", "sequence": "A" * 100},
        {"id": "D", "molecule_type": "protein", "sequence": "A" * 46},
    ]
    _payload, body = b9.call_boltz(fake, polymers, "bb50-g1", 3)
    assert body["constraints"] == b9.matched_constraint()
    assert body["diffusion_samples"] == 3
    assert body["recycling_steps"] == 8
    assert body["without_potentials"] is False


def test_subjects_and_settings_are_frozen():
    assert b9.SUBJECT_KEYS == ("bb50-g1", "bb54-g2", "bb25")
    assert b9.SETTINGS["bb50-g1"] == (8, 50, False)
    assert b9.SETTINGS["bb25"] == (8, 50, False)
    assert b9.SETTINGS["bb54-g2"] == (3, 50, True)


import json
import pytest


@pytest.fixture
def sandbox_run(monkeypatch, tmp_path):
    monkeypatch.setattr(b9, "OUT", tmp_path / "batch-09")
    monkeypatch.delenv("NVIDIA_API_KEY", raising=False)
    monkeypatch.setattr(b9.panel, "panel_metrics", lambda *args: {
        "buried_area_total": 100.0, "buried_area_peptide": 20.0,
        "peptide_buried_frac": 0.2, "pos6_contacts": 2,
    })
    return b9.OUT


def payload(n=3):
    return {"structures": [{"structure": f"data_test_{i}\n"} for i in range(n)]}


def install_client(monkeypatch, response):
    calls = []

    class FakeClient:
        def __init__(self, **kwargs):
            pass

        def post(self, path, body):
            calls.append(body)
            return response

    monkeypatch.setattr(b9, "NvidiaClient", FakeClient)
    return calls


def test_preflight_needs_no_key_and_writes_no_results(sandbox_run, monkeypatch, capsys):
    def forbidden(**kwargs):
        raise AssertionError("dry-run must not construct a client")

    monkeypatch.setattr(b9, "NvidiaClient", forbidden)
    assert b9.main(["--dry-run"]) == 0
    plan = json.loads(capsys.readouterr().out)
    assert plan["planned_calls"] == 12
    assert plan["boltz2_calls_this_run"] == 0
    assert len(plan["requests"]) == 12
    assert not sandbox_run.exists()


def test_resume_recovers_unscored_draws_without_another_api_call(sandbox_run, monkeypatch):
    calls = install_client(monkeypatch, payload())
    append = b9.append_checkpoint

    def interrupt_after_draw_zero(record):
        if record["draw"] == 1:
            raise RuntimeError("simulated interruption")
        append(record)

    monkeypatch.setattr(b9, "append_checkpoint", interrupt_after_draw_zero)
    selection = ["--designs", "bb50-g1", "--arms", "g12d"]
    with pytest.raises(RuntimeError, match="simulated interruption"):
        b9.main(selection)
    assert b9.done_call_keys() == set()
    assert len(b9.checkpoint_records()) == 1
    assert json.loads((sandbox_run / "matched_anchor_run.json").read_text())["status"] == "failed"

    monkeypatch.setattr(b9, "append_checkpoint", append)
    assert b9.main(selection) == 0
    assert len(calls) == 1
    assert len(b9.checkpoint_records()) == 3
    assert b9.done_call_keys() == {("bb50-g1", "g12d")}
    meta = json.loads((sandbox_run / "matched_anchor_run.json").read_text())
    assert meta["status"] == "completed"
    assert meta["boltz2_calls_this_run"] == 0
    assert b9.main(selection) == 0
    assert len(calls) == 1
    assert len(b9.checkpoint_records()) == 3


@pytest.mark.parametrize("response", [payload(0), payload(2), payload(4),
                                      {"structures": payload(1)["structures"] * 3}])
def test_bad_model_response_is_preserved_and_not_marked_complete(sandbox_run, monkeypatch, response):
    calls = install_client(monkeypatch, response)
    with pytest.raises(ValueError, match="structures"):
        b9.main(["--designs", "bb50-g1", "--arms", "g12d"])
    assert len(calls) == 1
    assert b9.done_call_keys() == set()
    assert (sandbox_run / "responses" / "bb50-g1_g12d.json").is_file()
    assert json.loads((sandbox_run / "matched_anchor_run.json").read_text())["status"] == "failed"


def test_call_cap_preserves_run_status_without_key(sandbox_run):
    assert b9.main(["--max-boltz-calls", "0"]) == 0
    meta = json.loads((sandbox_run / "matched_anchor_run.json").read_text())
    assert meta["status"] == "call_cap_reached"
    assert meta["completed_calls"] == 0
    assert meta["boltz2_calls_this_run"] == 0


def test_resume_rejects_missing_structure_and_mixed_protocol(sandbox_run, monkeypatch):
    install_client(monkeypatch, payload())
    b9.main(["--designs", "bb50-g1", "--arms", "g12d"])
    records = b9.checkpoint_records()
    records[0]["sampling"]["recycling_steps"] = 3
    with pytest.raises(ValueError, match="protocol mismatch"):
        b9.validate_records(records, b9.load_subjects(), 3)
    records = b9.checkpoint_records()
    (sandbox_run / records[0]["structure_file"]).unlink()
    with pytest.raises(ValueError, match="missing or changed"):
        b9.main(["--dry-run"])


def test_partial_legacy_call_requires_response_recovery(sandbox_run, monkeypatch):
    install_client(monkeypatch, payload())
    b9.main(["--designs", "bb50-g1", "--arms", "g12d"])
    first = b9.checkpoint_path().read_text().splitlines()[0]
    b9.checkpoint_path().write_text(first + "\n")
    (sandbox_run / "responses" / "bb50-g1_g12d.json").unlink()
    with pytest.raises(ValueError, match="recover the original response"):
        b9.main(["--designs", "bb50-g1", "--arms", "g12d"])


def test_duplicate_checkpoint_draw_is_rejected(sandbox_run, monkeypatch):
    install_client(monkeypatch, payload())
    b9.main(["--designs", "bb50-g1", "--arms", "g12d"])
    records = b9.checkpoint_records()
    b9.append_checkpoint(records[0])
    with pytest.raises(ValueError, match="duplicate checkpoint"):
        b9.done_call_keys()


def test_full_cached_panel_round_trip_to_analysis(sandbox_run, monkeypatch):
    calls = install_client(monkeypatch, payload())
    assert b9.main([]) == 0
    assert len(calls) == 12
    assert len(b9.checkpoint_records()) == 36
    spec = importlib.util.spec_from_file_location("batch09_roundtrip_analysis", ROOT / "scripts/batch09_analyze.py")
    analyzer = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(analyzer)
    monkeypatch.setattr(analyzer, "OUT", sandbox_run)
    monkeypatch.setattr(analyzer, "CHECKPOINT", b9.checkpoint_path())
    assert analyzer.main() == 0
    report = json.loads((sandbox_run / "matched_anchor_analysis.json").read_text())
    assert report["panel_complete"] is True
    assert set(report["designs"]) == set(b9.SUBJECT_KEYS)
    assert b9.main([]) == 0
    assert len(calls) == 12
