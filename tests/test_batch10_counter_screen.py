import importlib.util
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location(
    "batch10_counter_screen", ROOT / "scripts/batch10_counter_screen.py"
)
b10 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(b10)


def test_batch10_freezes_review_ranks_4_to_11_and_no_generation():
    subjects = b10.load_subjects()
    assert [v["rank"] for v in subjects.values()] == list(range(4, 12))
    assert len(subjects) == 8
    assert set(subjects) == set(b10.SETTINGS)
    assert b10.SCREEN_ARMS == ("g12d", "apo")
    assert b10.CONFIRM_ARMS == ("wt", "decoy")


def test_phase_screen_is_exactly_16_calls_when_empty():
    subjects = b10.load_subjects()
    wanted, arms = b10.selected_for_phase("screen", [], subjects)
    assert len(wanted) * len(arms) == 16


def test_advance_rule_is_directional_and_requires_direct_readout():
    subjects = b10.load_subjects()
    key = next(iter(subjects))

    def rows(gate_passes, g_area=800.0, apo_area=600.0):
        out = []
        for draw in range(3):
            out.append({
                "design": key, "arm": "g12d", "draw": draw,
                "row": {"asp6_gate_pass": draw < gate_passes,
                        "buried_area_total": g_area},
            })
            out.append({
                "design": key, "arm": "apo", "draw": draw,
                "row": {"asp6_gate_pass": False,
                        "buried_area_total": apo_area},
            })
        return out

    summary = b10.screen_summary(rows(2), {key: subjects[key]})[key]
    assert summary["advance_to_wt_decoy"] is True
    assert summary["apo_to_g12d_total_area_ratio"] == 0.75

    no_readout = b10.screen_summary(rows(1), {key: subjects[key]})[key]
    assert no_readout["advance_to_wt_decoy"] is False

    apo_dominant = b10.screen_summary(
        rows(3, g_area=800.0, apo_area=900.0), {key: subjects[key]}
    )[key]
    assert apo_dominant["advance_to_wt_decoy"] is False
