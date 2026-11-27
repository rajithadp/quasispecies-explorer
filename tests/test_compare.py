import pandas as pd
import pytest

import compare_samples
from quasigraph import read_ivar_variants
from quasigraph.compare import (
    compare_samples as compare,
    drop_classes,
    minor_variants,
    run_id_from_path,
    substitution_counts,
    summarize_sample,
)

HEADER = "REGION\tPOS\tREF\tALT\tREF_DP\tREF_RV\tREF_QUAL\tALT_DP\tALT_RV\tALT_QUAL\tALT_FREQ\tTOTAL_DP\tPVAL\tPASS\n"


def table(*rows):
    body = "".join(f"r\t{p}\t{r}\t{a}\t0\t0\t0\t0\t0\t0\t{f}\t{d}\t0\tTRUE\n" for p, r, a, f, d in rows)
    return HEADER + body


ROWS = [
    (100, "G", "T", 0.03, 1000),   # oxidative class, minor
    (200, "C", "A", 0.04, 1000),   # oxidative class, minor
    (300, "C", "T", 0.05, 1000),   # other class, minor
    (400, "G", "T", 0.99, 1000),   # fixed difference of the oxidative class: must be kept
    (500, "A", "G", 0.01, 1000),   # below the frequency threshold
    (600, "C", "A", 0.30, 50),     # too shallow
    (700, "A", "+AG", 0.20, 1000), # indel: not a substitution
]


def v():
    return read_ivar_variants(table(*ROWS))


def test_minor_variants_apply_frequency_depth_and_snv_filters():
    m = minor_variants(v())
    assert sorted(m["position"]) == [100, 200, 300]


def test_substitution_counts():
    c = substitution_counts(v())
    assert c == {"G>T": 1, "C>A": 1, "C>T": 1}


def test_drop_classes_keeps_fixed_differences_and_other_classes():
    kept = drop_classes(v())
    assert sorted(kept["position"]) == [300, 400, 500, 700]  # 600 is a low-depth C>A minor variant: dropped too


def test_summary_with_and_without_the_class():
    s = summarize_sample(v(), genome_length=1000)
    assert s["minor_variants"] == 3 and s["oxidative_class"] == 2
    assert s["oxidative_fraction"] == pytest.approx(2 / 3)
    assert s["variable_sites_all"] == 3 and s["variable_sites_without_GT_CA"] == 1
    assert s["fixed_differences"] == 1
    # pi is the sum of 1 - sum(p^2) over variable sites, divided by genome length
    het = lambda f: 2 * f * (1 - f)
    assert s["pi_all"] == pytest.approx((het(0.03) + het(0.04) + het(0.05)) / 1000)
    assert s["pi_without_GT_CA"] == pytest.approx(het(0.05) / 1000)


def test_sample_without_minor_variants_has_nan_fraction():
    s = summarize_sample(read_ivar_variants(table((10, "C", "T", 0.99, 1000))), genome_length=100)
    assert s["minor_variants"] == 0 and s["oxidative_fraction"] != s["oxidative_fraction"]


def test_run_id_from_path():
    assert run_id_from_path("samples/SRR12447392.variants.tsv") == "SRR12447392"


def test_compare_joins_metadata(tmp_path):
    (tmp_path / "A.variants.tsv").write_text(table(*ROWS))
    (tmp_path / "B.variants.tsv").write_text(table((300, "C", "T", 0.05, 1000)))
    meta = pd.DataFrame({"run_accession": ["A", "B", "C"], "country": ["x", "y", "z"]})
    out = compare([tmp_path / "A.variants.tsv", tmp_path / "B.variants.tsv"], 1000, metadata=meta)
    assert list(out["sample"]) == ["A", "B"] and list(out["country"]) == ["x", "y"]
    with pytest.raises(ValueError, match="run_accession"):
        compare([tmp_path / "A.variants.tsv"], 1000, metadata=pd.DataFrame({"id": ["A"]}), id_column="run_accession")


def test_cli_writes_table_and_figure(tmp_path, capsys):
    p = tmp_path / "A.variants.tsv"
    p.write_text(table(*ROWS))
    rc = compare_samples.main(["--variants", str(p), "--genome-length", "1000", "--out", str(tmp_path / "out")])
    assert rc == 0
    assert (tmp_path / "out.tsv").exists() and (tmp_path / "out.png").stat().st_size > 0
    assert "67%" in capsys.readouterr().out


def test_threshold_boundaries_are_inclusive_where_documented():
    # exactly min_freq and exactly min_depth count as minor variants; exactly 0.5 is a fixed difference
    rows = [(10, "C", "T", 0.02, 100), (20, "C", "T", 0.5, 1000), (30, "C", "T", 0.0199, 1000), (40, "C", "T", 0.05, 99)]
    m = minor_variants(read_ivar_variants(table(*rows)), min_freq=0.02, min_depth=100)
    assert list(m["position"]) == [10]


def test_oxidative_variant_at_exactly_half_is_kept_as_fixed_difference():
    vv = read_ivar_variants(table((10, "G", "T", 0.5, 1000), (20, "G", "T", 0.4999, 1000)))
    assert list(drop_classes(vv)["position"]) == [10]
    assert summarize_sample(vv, 1000)["fixed_differences"] == 1
