import csv

import pytest

import validate_submission as vs


def write_csv(path, rows, header=("name", "sequence")):
    with path.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(header)
        for row in rows:
            w.writerow(row)
    return path


def seq(n, aa="A"):
    """A valid-alphabet sequence of length n."""
    return (aa * n)[:n]


GOOD = "MKQLEDKVEELLSKNYHLENEVARLKKLVGER"  # 32 aa, canonical letters


def test_valid_submission(tmp_path):
    p = write_csv(tmp_path / "s.csv", [(f"d{i}", GOOD) for i in range(5)])
    assert vs.validate_file(p) == []


def test_valid_at_exactly_20_designs(tmp_path):
    p = write_csv(tmp_path / "s.csv", [(f"d{i}", GOOD) for i in range(20)])
    assert vs.validate_file(p) == []


def test_21_designs_rejected(tmp_path):
    p = write_csv(tmp_path / "s.csv", [(f"d{i}", GOOD) for i in range(21)])
    errors = vs.validate_file(p)
    assert any("21 designs" in e and "limit of 20" in e for e in errors)


def test_chain_of_9_aa_rejected(tmp_path):
    p = write_csv(tmp_path / "s.csv", [("short", seq(9))])
    errors = vs.validate_file(p)
    assert any("9 aa" in e and "minimum of 10" in e for e in errors)


def test_chain_of_10_aa_accepted(tmp_path):
    p = write_csv(tmp_path / "s.csv", [("edge", seq(10))])
    assert vs.validate_file(p) == []


def test_chain_of_250_aa_accepted(tmp_path):
    p = write_csv(tmp_path / "s.csv", [("edge", seq(250))])
    assert vs.validate_file(p) == []


def test_chain_of_251_aa_rejected(tmp_path):
    p = write_csv(tmp_path / "s.csv", [("long", seq(251))])
    errors = vs.validate_file(p)
    assert any("251 aa" in e and "maximum of 250" in e for e in errors)


@pytest.mark.parametrize("bad_letter", ["X", "B", "Z", "U", "O", "1", "*"])
def test_bad_amino_acid_rejected(tmp_path, bad_letter):
    p = write_csv(tmp_path / "s.csv", [("bad", GOOD[:-1] + bad_letter)])
    errors = vs.validate_file(p)
    assert any("invalid amino acid" in e for e in errors)


def test_valid_vh_vl_pair(tmp_path):
    p = write_csv(tmp_path / "s.csv", [("scfv1", f"{GOOD}:{GOOD}")],
                  header=("name", "sequence", "molecule_class"))
    # rewrite with the third column populated
    with p.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(("name", "sequence", "molecule_class"))
        w.writerow(("scfv1", f"{GOOD}:{GOOD}", "scFv"))
    assert vs.validate_file(p) == []


def test_malformed_vh_vl_two_colons(tmp_path):
    p = write_csv(tmp_path / "s.csv", [("scfv_bad", f"{GOOD}:{GOOD}:{GOOD}")])
    errors = vs.validate_file(p)
    assert any("exactly one ':'" in e for e in errors)


def test_malformed_vh_vl_empty_chain(tmp_path):
    p = write_csv(tmp_path / "s.csv", [("scfv_bad", f"{GOOD}:")])
    errors = vs.validate_file(p)
    assert any("empty" in e for e in errors)


def test_short_vl_in_pair_rejected(tmp_path):
    p = write_csv(tmp_path / "s.csv", [("scfv_bad", f"{GOOD}:{seq(9)}")])
    errors = vs.validate_file(p)
    assert any("VL" in e and "minimum of 10" in e for e in errors)


def test_scfv_class_without_colon_rejected(tmp_path):
    p = tmp_path / "s.csv"
    with p.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(("name", "sequence", "molecule_class"))
        w.writerow(("x", GOOD, "scFv"))
    errors = vs.validate_file(p)
    assert any("{VH}:{VL}" in e for e in errors)


def test_missing_sequence_column(tmp_path):
    p = tmp_path / "s.csv"
    p.write_text("name,seq\nd1,AAAAAAAAAA\n")
    with pytest.raises(vs.ValidationError, match="missing required column"):
        vs.validate_file(p)


def test_missing_name_column(tmp_path):
    p = tmp_path / "s.csv"
    p.write_text("id,sequence\nd1,AAAAAAAAAA\n")
    with pytest.raises(vs.ValidationError, match="missing required column"):
        vs.validate_file(p)


def test_missing_name_value(tmp_path):
    p = write_csv(tmp_path / "s.csv", [("", GOOD)])
    errors = vs.validate_file(p)
    assert any("missing name" in e for e in errors)


def test_duplicate_names_rejected(tmp_path):
    p = write_csv(tmp_path / "s.csv", [("same", GOOD), ("same", GOOD)])
    errors = vs.validate_file(p)
    assert any("duplicate name" in e for e in errors)


def test_empty_file_raises(tmp_path):
    p = tmp_path / "s.csv"
    p.write_text("")
    with pytest.raises(vs.ValidationError, match="empty"):
        vs.validate_file(p)


def test_missing_file_raises(tmp_path):
    with pytest.raises(vs.ValidationError, match="file not found"):
        vs.validate_file(tmp_path / "nope.csv")


def test_unsupported_extension(tmp_path):
    p = tmp_path / "s.txt"
    p.write_text("name,sequence\nd1,AAAAAAAAAA\n")
    with pytest.raises(vs.ValidationError, match="unsupported file type"):
        vs.validate_file(p)


def test_fasta_valid(tmp_path):
    p = tmp_path / "s.fasta"
    p.write_text(f">d1\n{GOOD}\n>d2\n{GOOD[:16]}\n{GOOD[16:]}\n")
    assert vs.validate_file(p) == []


def test_fasta_short_chain_rejected(tmp_path):
    p = tmp_path / "s.fasta"
    p.write_text(f">d1\n{seq(9)}\n")
    errors = vs.validate_file(p)
    assert any("minimum of 10" in e for e in errors)


def test_fasta_too_many_records(tmp_path):
    p = tmp_path / "s.fasta"
    p.write_text("".join(f">d{i}\n{GOOD}\n" for i in range(21)))
    errors = vs.validate_file(p)
    assert any("limit of 20" in e for e in errors)


def test_fasta_without_header_raises(tmp_path):
    p = tmp_path / "s.fasta"
    p.write_text(f"{GOOD}\n")
    with pytest.raises(vs.ValidationError, match="before any"):
        vs.validate_file(p)


def test_oversize_file_rejected(tmp_path, monkeypatch):
    monkeypatch.setattr(vs, "MAX_FILE_BYTES", 100)
    p = write_csv(tmp_path / "s.csv", [(f"d{i}", GOOD) for i in range(5)])
    errors = vs.validate_file(p)
    assert any("above the 5 MB limit" in e for e in errors)


def test_cli_exit_codes(tmp_path, capsys):
    good = write_csv(tmp_path / "good.csv", [("d1", GOOD)])
    bad = write_csv(tmp_path / "bad.csv", [("d1", seq(9))])
    assert vs.main([str(good)]) == 0
    assert vs.main([str(bad)]) == 1
    assert vs.main([str(tmp_path / "missing.csv")]) == 1


def test_lowercase_sequence_accepted(tmp_path):
    p = write_csv(tmp_path / "s.csv", [("d1", GOOD.lower())])
    assert vs.validate_file(p) == []
