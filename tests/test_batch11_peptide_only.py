import importlib.util
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "batch11_peptide_only", ROOT / "scripts" / "batch11_peptide_only.py"
)
b11 = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(b11)


def test_peptide_only_contig_is_not_hla_contig():
    assert b11.peptide_only_contig(46, "C") == "C1-10/0 46-46"
    assert "A1-180" not in b11.peptide_only_contig(46, "C")


def test_design_backbone_assessment_requires_two_residue_span_and_no_clash():
    pos6 = {
        "sidechain_contact_atoms": 2,
        "sidechain_min_distance": 3.4,
        "binder_residues": [7],
    }
    pos8 = {
        "sidechain_contact_atoms": 1,
        "sidechain_min_distance": 4.1,
        "binder_residues": [19],
    }
    ok, reasons = b11.design_backbone_assessment(
        pos6, pos8, {"severe_clash_atoms": 0}
    )
    assert ok
    assert reasons == []

    bad, reasons = b11.design_backbone_assessment(
        pos6,
        {**pos8, "binder_residues": [7]},
        {"severe_clash_atoms": 0},
    )
    assert not bad
    assert any("two binder residues" in reason for reason in reasons)


def test_design_backbone_assessment_rejects_reinserted_hla_clash():
    metric = {
        "sidechain_contact_atoms": 1,
        "sidechain_min_distance": 4.0,
        "binder_residues": [1],
    }
    pos8 = {**metric, "binder_residues": [2]}
    ok, reasons = b11.design_backbone_assessment(
        metric, pos8, {"severe_clash_atoms": 3}
    )
    assert not ok
    assert any("severe HLA clash" in reason for reason in reasons)


def test_paired_anchor_variant_uses_distinct_positions():
    sequence = "A" * 12
    shell6 = [{"position": 3, "distance": 3.0}]
    shell8 = [
        {"position": 3, "distance": 3.0},
        {"position": 8, "distance": 3.5},
    ]
    variant = b11.paired_anchor_variant(sequence, shell6, shell8)
    assert variant[2] == "R"
    assert variant[7] == "Y"
    assert variant.count("R") == 1
    assert variant.count("Y") == 1


def test_screen_selection_prefers_backbone_diversity():
    candidates = [
        {"sequence": "A" * 10, "backbone": "bb1", "paired_score": 10},
        {"sequence": "C" * 10, "backbone": "bb1", "paired_score": 9},
        {"sequence": "D" * 10, "backbone": "bb2", "paired_score": 8},
    ]
    selected = b11.select_screen_candidates(candidates, 2)
    assert [c["backbone"] for c in selected] == ["bb1", "bb2"]


def _row(design, arm, draw, area, asp=False, p8=0):
    return {
        "design": design,
        "arm": arm,
        "draw": draw,
        "row": {
            "buried_area_total": area,
            "asp6_gate_pass": asp,
            "pos8_sidechain_contacts": p8,
        },
    }


def test_screen_summary_uses_frozen_three_leg_rule():
    sequence = "ACDEFGHIKL"
    design = b11._design_key(sequence)
    selected = [{
        "sequence": sequence,
        "backbone": "bb1",
        "variant": "mpnn",
        "paired_score": 3.2,
    }]
    rows = []
    rows += [
        _row(design, "g12d", 0, 1000, True, 2),
        _row(design, "g12d", 1, 900, True, 1),
        _row(design, "g12d", 2, 1100, False, 0),
    ]
    rows += [
        _row(design, "apo", 0, 700),
        _row(design, "apo", 1, 800),
        _row(design, "apo", 2, 900),
    ]
    summary = b11.summarize_screen(rows, selected)[design]
    assert summary["g12d_Asp6_gate_passes"] == 2
    assert summary["g12d_Val8_contact_draws"] == 2
    assert summary["apo_to_g12d_total_area_ratio"] == 0.8
    assert summary["advance_to_wt_decoy"] is True

    rows[4]["row"]["buried_area_total"] = 1200
    rows[5]["row"]["buried_area_total"] = 1300
    summary = b11.summarize_screen(rows, selected)[design]
    assert summary["apo_to_g12d_total_area_ratio"] > 1
    assert summary["advance_to_wt_decoy"] is False


def test_generation_service_failure_is_infrastructure_not_geometry():
    records = [
        {
            "backbone": "bb0",
            "stage": "rfdiffusion",
            "error_type": "NvidiaAuthError",
            "error": "NvidiaAuthError: rejected",
        },
        {
            "backbone": "bb1",
            "ok": False,
            "stage": "geometry",
            "error_type": "StructureAnalysisError",
            "error": "geometry parse issue",
        },
        {
            "backbone": "bb2",
            "ok": True,
            "mpnn_error_type": "NvidiaAuthError",
            "mpnn_error": "NvidiaAuthError: rejected",
        },
    ]
    failures = b11._design_infrastructure_failures(records)
    assert [(x["backbone"], x["service"]) for x in failures] == [
        ("bb0", "rfdiffusion"),
        ("bb2", "proteinmpnn"),
    ]
