import pandas as pd
import pytest

import fixed_differences
from quasigraph import read_ivar_variants
from quasigraph.annotation import Annotator, Cds
from quasigraph.shared_fixed import difference_matrix, fixed_set, groups, shared_table

HEADER = "REGION\tPOS\tREF\tALT\tREF_DP\tREF_RV\tREF_QUAL\tALT_DP\tALT_RV\tALT_QUAL\tALT_FREQ\tTOTAL_DP\tPVAL\tPASS\n"


def table(*rows):
    return HEADER + "".join(f"r\t{p}\t{r}\t{a}\t0\t0\t0\t0\t0\t0\t{f}\t{d}\t0\tTRUE\n" for p, r, a, f, d in rows)


def test_fixed_set_boundaries_and_exclusions():
    v = read_ivar_variants(table((10, "A", "G", 0.5, 100),      # exactly at both limits: kept
                                 (20, "A", "G", 0.49, 1000),    # just below frequency
                                 (30, "A", "G", 0.9, 99),       # just below depth
                                 (40, "A", "+AG", 0.9, 1000),   # indel
                                 (50, "C", "T", 1.0, 5000)))
    assert fixed_set(v) == {(10, "A", "G"), (50, "C", "T")}
    assert fixed_set(v, min_depth=1000) == {(50, "C", "T")}


def test_shared_table_counts_carriers_and_missing_samples_and_sorts():
    sets = {"S1": {(10, "A", "G"), (50, "C", "T")}, "S2": {(50, "C", "T")}, "S3": {(50, "C", "T"), (20, "A", "C")}}
    t = shared_table(sets)
    assert list(t["position"]) == [50, 10, 20]                  # most shared first, then by position
    assert list(t["n_samples"]) == [3, 1, 1]
    r = t.set_index("position")
    assert r.loc[50, "carried_by"] == "S1,S2,S3" and r.loc[50, "not_in"] == ""
    assert r.loc[10, "carried_by"] == "S1" and r.loc[10, "not_in"] == "S2,S3"
    assert r.loc[20, "carried_by"] == "S3" and r.loc[20, "not_in"] == "S1,S2"
    assert (t["annotation"] == "").all()


def test_shared_table_with_no_changes_is_empty_with_columns():
    t = shared_table({"S1": set()})
    assert t.empty and "n_samples" in t.columns


def test_annotation_label_is_added():
    ref = "CCC" + "ATGGATTGG"
    ann = Annotator([Cds("c", "A", "", "+", ((4, 12),))], ref)
    t = shared_table({"S1": {(7, "G", "T"), (1, "C", "T")}}, ann)
    lab = dict(zip(t["position"], t["annotation"]))
    assert lab[7] == "A:D2Y:missense" and lab[1] == "noncoding"


def test_groups_identical_sets_largest_first():
    sets = {"S1": {(1, "A", "C")}, "S2": {(1, "A", "C")}, "S3": {(1, "A", "C"), (2, "A", "C")}, "S4": set()}
    g = groups(sets)
    assert g[0] == (("S1", "S2"), 1)
    assert (("S3",), 2) in g and ((("S4",), 0)) in g and len(g) == 3


def test_difference_matrix_is_symmetric_difference():
    sets = {"S1": {(1, "A", "C"), (2, "A", "C")}, "S2": {(1, "A", "C"), (3, "A", "C")}, "S3": {(1, "A", "C"), (2, "A", "C")}}
    m = difference_matrix(sets)
    assert m.loc["S1", "S2"] == 2 and m.loc["S2", "S1"] == 2 and m.loc["S1", "S3"] == 0
    assert all(m.loc[s, s] == 0 for s in sets)


def test_cli_end_to_end_and_errors(tmp_path, capsys):
    ref = "CCC" + "ATGGATTGG"
    (tmp_path / "ref.fa").write_text(">chrT\n" + ref + "\n")
    (tmp_path / "a.gff3").write_text("##gff-version 3\nchrT\tx\tCDS\t4\t12\t.\t+\t0\tID=c;gene=A\n")
    (tmp_path / "S1.variants.tsv").write_text(table((7, "G", "T", 0.9, 500), (1, "C", "T", 0.1, 500)))
    (tmp_path / "S2.variants.tsv").write_text(table((7, "G", "T", 0.8, 500), (12, "G", "A", 0.7, 500)))
    base = ["--variants", str(tmp_path / "S1.variants.tsv"), str(tmp_path / "S2.variants.tsv"), "--out", str(tmp_path / "o")]
    assert fixed_differences.main(base + ["--ref", str(tmp_path / "ref.fa"), "--gff", str(tmp_path / "a.gff3")]) == 0
    out = capsys.readouterr().out
    assert "2 fixed differences in 2 samples" in out and "A:D2Y:missense" in out and "A:W3*:stop_gained" in out
    t = pd.read_csv(tmp_path / "o.tsv", sep="\t")
    assert list(t["position"]) == [7, 12] and list(t["n_samples"]) == [2, 1]   # the 10% change at 1 is not fixed
    assert fixed_differences.main(base) == 0                                   # works without annotation
    assert fixed_differences.main(base + ["--gff", str(tmp_path / "a.gff3")]) == 2
    assert "--gff needs --ref" in capsys.readouterr().err
    (tmp_path / "e.gff3").write_text("##gff-version 3\n")
    assert fixed_differences.main(base + ["--ref", str(tmp_path / "ref.fa"), "--gff", str(tmp_path / "e.gff3")]) == 2
    assert "no CDS" in capsys.readouterr().err
    (tmp_path / "x").mkdir()
    (tmp_path / "x" / "S1.variants.tsv").write_text(table((7, "G", "T", 0.9, 500)))
    dup = ["--variants", str(tmp_path / "S1.variants.tsv"), str(tmp_path / "x" / "S1.variants.tsv")]
    assert fixed_differences.main(dup) == 2 and "two variants tables" in capsys.readouterr().err


def test_cli_min_depth_option_is_applied(tmp_path, capsys):
    (tmp_path / "S1.variants.tsv").write_text(table((7, "G", "T", 0.9, 500), (9, "A", "C", 0.9, 2000)))
    args = ["--variants", str(tmp_path / "S1.variants.tsv"), "--out", str(tmp_path / "o")]
    assert fixed_differences.main(args) == 0
    assert "2 fixed differences" in capsys.readouterr().out
    assert fixed_differences.main(args + ["--min-depth", "1000"]) == 0
    assert "1 fixed differences" in capsys.readouterr().out
