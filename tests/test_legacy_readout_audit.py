import importlib.util
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location(
    "legacy_readout_audit", ROOT / "scripts/audit_legacy_readout.py"
)
audit_module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(audit_module)


def test_legacy_batches_do_not_add_direct_readout_candidates():
    report = audit_module.audit(ROOT)
    assert report["model_calls"] == 0
    assert report["readout_status"] == "recorded-unverified"
    assert report["totals"]["checkpoint_records"] == 163
    assert report["totals"]["unique_sequences"] == 144
    assert report["totals"]["two_leg_sequences"] == 0
    assert report["totals"]["burial_only_sequences"] == 7
    assert report["two_leg_candidates"] == []


def test_legacy_source_breakdown_is_frozen():
    report = audit_module.audit(ROOT)
    observed = {
        path: (
            info["records"],
            info["unique_sequences"],
            info["two_leg_sequences"],
            info["oxygen_only_sequences"],
            info["burial_only_sequences"],
        )
        for path, info in report["sources"].items()
    }
    assert observed == {
        "results/batch-01/checkpoint.jsonl": (75, 56, 0, 0, 2),
        "results/batch-02/checkpoint.jsonl": (60, 60, 0, 0, 4),
        "results/batch-03/checkpoint.jsonl": (28, 28, 0, 0, 1),
    }
