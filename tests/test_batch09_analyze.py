from __future__ import annotations

import copy
import importlib.util
import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("batch09_analyze", ROOT / "scripts/batch09_analyze.py")
analyzer = importlib.util.module_from_spec(spec)
spec.loader.exec_module(analyzer)


@pytest.fixture
def panel(monkeypatch):
    # Protocol/structure validation is covered by the runner tests.
    monkeypatch.setattr(analyzer.experiment, "validate_records", lambda *args: None)
    return [{"design": key, "arm": arm, "draw": draw,
             "row": {"buried_area_total": 100, "buried_area_peptide": 20,
                     "peptide_buried_frac": 0.2, "pos6_contacts": 2,
                     "analysis_error": "", "readout_version": analyzer.experiment.panel.READOUT_VERSION,
                     "asp6_buried_frac": 0.3}}
            for key in analyzer.experiment.SUBJECT_KEYS
            for arm in analyzer.experiment.ARMS for draw in range(3)]


def test_complete_panel_can_be_summarized(panel, monkeypatch, tmp_path):
    monkeypatch.setattr(analyzer, "OUT", tmp_path)
    checkpoint = tmp_path / "matched_anchor_checkpoint.jsonl"
    checkpoint.write_text("\n".join(json.dumps(r) for r in panel))
    monkeypatch.setattr(analyzer, "CHECKPOINT", checkpoint)
    assert analyzer.main() == 0
    report = json.loads((tmp_path / "matched_anchor_analysis.json").read_text())
    assert report["panel_complete"] is True
    assert report["new_hard_gate"] is False
    assert all(arm["n"] == 3 for d in report["designs"].values() for arm in d["arms"].values())


@pytest.mark.parametrize("failure", ["missing", "duplicate", "analysis_error", "nan", "missing_metric"])
def test_incomplete_or_invalid_panel_is_rejected(panel, failure):
    data = copy.deepcopy(panel)
    if failure == "missing":
        data.pop()
    elif failure == "duplicate":
        data.append(data[0])
    elif failure == "analysis_error":
        data[0]["row"]["analysis_error"] = "scoring failed"
    elif failure == "nan":
        data[0]["row"]["buried_area_total"] = float("nan")
    else:
        del data[0]["row"]["buried_area_total"]
    with pytest.raises(ValueError):
        analyzer.validate_panel(data)


def test_legacy_flags_require_explicit_historical_mode(panel):
    for r in panel:
        r['row'].pop('readout_version')
    with pytest.raises(ValueError, match='legacy readout'):
        analyzer.validate_panel(panel)
    analyzer.validate_panel(panel, historical=True)


def test_historical_mode_does_not_overwrite_original_analysis(panel,monkeypatch,tmp_path):
    for r in panel:
        r['row'].pop('readout_version')
    checkpoint=tmp_path/'matched_anchor_checkpoint.jsonl'
    checkpoint.write_text('\n'.join(json.dumps(r) for r in panel))
    monkeypatch.setattr(analyzer,'OUT',tmp_path)
    monkeypatch.setattr(analyzer,'CHECKPOINT',checkpoint)
    original=tmp_path/'matched_anchor_analysis.json'
    original.write_text('original evidence sentinel')
    assert analyzer.main(['--historical-recorded-metrics'])==0
    assert original.read_text()=='original evidence sentinel'
    rebuilt=json.loads((tmp_path/'historical_recorded_analysis_rebuilt.json').read_text())
    assert rebuilt['historical_recorded_metrics'] is True
    assert rebuilt['readout_version']=='legacy-recorded-flags'
