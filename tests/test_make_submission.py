import csv

import pytest

import make_submission as ms
import validate_submission as vs

GOOD = "MKQLEDKVEELLSKNYHLENEVARLKKLVGER"


def write_scored(path, rows, header=("name", "sequence", "specificity_gap")):
    with path.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(header)
        for row in rows:
            w.writerow(row)
    return path


def read_rows(path):
    with path.open(newline="") as fh:
        return list(csv.DictReader(fh))


def test_builds_valid_submission(tmp_path):
    src = write_scored(tmp_path / "scored.csv",
                       [(f"d{i}", GOOD, i / 10) for i in range(30)])
    out = ms.build(src, tmp_path / "sub.csv", top_n=20)
    rows = read_rows(out)
    assert len(rows) == 20
    assert vs.validate_file(out) == []


def test_ranks_by_specificity_gap_descending(tmp_path):
    src = write_scored(tmp_path / "scored.csv", [
        ("low", GOOD, 0.1), ("high", GOOD, 0.9), ("mid", GOOD, 0.5),
    ])
    out = ms.build(src, tmp_path / "sub.csv", top_n=3)
    assert [r["name"] for r in read_rows(out)] == ["high", "mid", "low"]


def test_top_n_truncates(tmp_path):
    src = write_scored(tmp_path / "scored.csv", [
        ("a", GOOD, 0.1), ("b", GOOD, 0.9), ("c", GOOD, 0.5),
    ])
    out = ms.build(src, tmp_path / "sub.csv", top_n=2)
    assert [r["name"] for r in read_rows(out)] == ["b", "c"]


def test_ascending_rank_for_kd(tmp_path):
    src = write_scored(tmp_path / "scored.csv",
                       [("a", GOOD, 5.0), ("b", GOOD, 1.0)],
                       header=("name", "sequence", "kd_nM"))
    out = ms.build(src, tmp_path / "sub.csv", top_n=2,
                   rank_column="kd_nm", ascending=True)
    assert [r["name"] for r in read_rows(out)] == ["b", "a"]


def test_falls_back_to_score_column(tmp_path):
    src = write_scored(tmp_path / "scored.csv",
                       [("a", GOOD, 0.2), ("b", GOOD, 0.8)],
                       header=("name", "sequence", "score"))
    out = ms.build(src, tmp_path / "sub.csv", top_n=2)
    assert [r["name"] for r in read_rows(out)] == ["b", "a"]


def test_preserves_order_without_score_columns(tmp_path):
    src = write_scored(tmp_path / "scored.csv",
                       [("first", GOOD), ("second", GOOD)],
                       header=("name", "sequence"))
    out = ms.build(src, tmp_path / "sub.csv", top_n=2)
    assert [r["name"] for r in read_rows(out)] == ["first", "second"]


def test_unscored_rows_sort_last(tmp_path):
    src = write_scored(tmp_path / "scored.csv",
                       [("blank", GOOD, ""), ("scored", GOOD, 0.3)])
    out = ms.build(src, tmp_path / "sub.csv", top_n=2)
    assert [r["name"] for r in read_rows(out)] == ["scored", "blank"]


def test_molecule_class_passed_through(tmp_path):
    src = write_scored(tmp_path / "scored.csv",
                       [("s1", f"{GOOD}:{GOOD}", "scFv")],
                       header=("name", "sequence", "molecule_class"))
    out = ms.build(src, tmp_path / "sub.csv", top_n=1)
    rows = read_rows(out)
    assert rows[0]["molecule_class"] == "scFv"
    assert vs.validate_file(out) == []


def test_molecule_class_column_omitted_when_unused(tmp_path):
    src = write_scored(tmp_path / "scored.csv", [("a", GOOD, 0.5)])
    out = ms.build(src, tmp_path / "sub.csv", top_n=1)
    with out.open() as fh:
        assert fh.readline().strip() == "name,sequence"


def test_top_n_above_20_rejected(tmp_path):
    src = write_scored(tmp_path / "scored.csv", [("a", GOOD, 0.5)])
    with pytest.raises(ms.ValidationError, match="exceeds the per-submission limit"):
        ms.build(src, tmp_path / "sub.csv", top_n=21)


def test_top_n_zero_rejected(tmp_path):
    src = write_scored(tmp_path / "scored.csv", [("a", GOOD, 0.5)])
    with pytest.raises(ms.ValidationError, match=">= 1"):
        ms.build(src, tmp_path / "sub.csv", top_n=0)


def test_missing_sequence_column_rejected(tmp_path):
    src = tmp_path / "scored.csv"
    src.write_text("name,specificity_gap\na,0.5\n")
    with pytest.raises(ms.ValidationError, match="missing required column 'sequence'"):
        ms.build(src, tmp_path / "sub.csv")


def test_missing_name_column_rejected(tmp_path):
    src = tmp_path / "scored.csv"
    src.write_text("sequence,specificity_gap\n%s,0.5\n" % GOOD)
    with pytest.raises(ms.ValidationError, match="missing required column 'name'"):
        ms.build(src, tmp_path / "sub.csv")


def test_unknown_rank_column_rejected(tmp_path):
    src = write_scored(tmp_path / "scored.csv", [("a", GOOD, 0.5)])
    with pytest.raises(ms.ValidationError, match="rank column 'nope' not found"):
        ms.build(src, tmp_path / "sub.csv", rank_column="nope")


def test_empty_input_rejected(tmp_path):
    src = tmp_path / "scored.csv"
    src.write_text("name,sequence,specificity_gap\n")
    with pytest.raises(ms.ValidationError, match="contains no designs"):
        ms.build(src, tmp_path / "sub.csv")


def test_invalid_design_blocks_output(tmp_path):
    """A 9 aa design must abort the build and leave no output file behind."""
    src = write_scored(tmp_path / "scored.csv", [("tooshort", "AAAAAAAAA", 0.9)])
    out = tmp_path / "sub.csv"
    with pytest.raises(ms.ValidationError, match="would be invalid"):
        ms.build(src, out, top_n=1)
    assert not out.exists()


def test_bad_amino_acid_blocks_output(tmp_path):
    src = write_scored(tmp_path / "scored.csv", [("badaa", GOOD[:-1] + "X", 0.9)])
    out = tmp_path / "sub.csv"
    with pytest.raises(ms.ValidationError, match="invalid amino acid"):
        ms.build(src, out, top_n=1)
    assert not out.exists()


def test_malformed_vh_vl_blocks_output(tmp_path):
    src = write_scored(tmp_path / "scored.csv",
                       [("badscfv", f"{GOOD}:{GOOD}:{GOOD}", 0.9)])
    out = tmp_path / "sub.csv"
    with pytest.raises(ms.ValidationError, match="exactly one"):
        ms.build(src, out, top_n=1)
    assert not out.exists()


def test_duplicate_names_block_output(tmp_path):
    src = write_scored(tmp_path / "scored.csv",
                       [("same", GOOD, 0.9), ("same", GOOD, 0.8)])
    with pytest.raises(ms.ValidationError, match="duplicate name"):
        ms.build(src, tmp_path / "sub.csv", top_n=2)


def test_cli(tmp_path, capsys):
    src = write_scored(tmp_path / "scored.csv", [("a", GOOD, 0.5)])
    out = tmp_path / "sub.csv"
    assert ms.main([str(src), "-o", str(out), "-n", "1"]) == 0
    assert out.exists()
    bad = write_scored(tmp_path / "bad.csv", [("a", "AAAAAAAAA", 0.5)])
    assert ms.main([str(bad), "-o", str(tmp_path / "bad_out.csv")]) == 1
