"""Pure-function tests for fetch_target (no network, no structure parsing)."""
import pytest

import fetch_target as ft


def test_wt_peptide_construction():
    assert ft.mutate_peptide(ft.EXPECTED_PEPTIDE, ft.MUTATED_POSITION, "G") == ft.WT_PEPTIDE


def test_peptides_differ_at_exactly_one_position():
    diffs = [
        i for i, (a, b) in enumerate(zip(ft.EXPECTED_PEPTIDE, ft.WT_PEPTIDE), start=1)
        if a != b
    ]
    assert diffs == [ft.MUTATED_POSITION]
    assert len(ft.EXPECTED_PEPTIDE) == len(ft.WT_PEPTIDE) == 10


def test_mutate_peptide_out_of_range():
    with pytest.raises(ValueError):
        ft.mutate_peptide("VVVGADGVGK", 0, "G")
    with pytest.raises(ValueError):
        ft.mutate_peptide("VVVGADGVGK", 11, "G")


def test_write_fasta_wraps_at_60(tmp_path):
    out = tmp_path / "x.fasta"
    ft.write_fasta([("hdr", "A" * 130)], out)
    lines = out.read_text().splitlines()
    assert lines[0] == ">hdr"
    assert [len(x) for x in lines[1:]] == [60, 60, 10]
