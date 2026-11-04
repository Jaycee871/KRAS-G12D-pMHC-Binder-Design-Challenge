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
