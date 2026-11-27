import csv
from pathlib import Path

import pandas as pd
import pytest

from quasigraph import frequency_spectrum, read_ivar_variants, site_summary, site_table
from quasigraph.pipeline_utils import describe_cliquesnv_error, read_cliquesnv_json
from quasigraph.plot import plot_site_positions, plot_site_spectrum

REAL = Path(__file__).parent / "data" / "real_variants.tsv"
HEADER = "REGION\tPOS\tREF\tALT\tREF_DP\tREF_RV\tREF_QUAL\tALT_DP\tALT_RV\tALT_QUAL\tALT_FREQ\tTOTAL_DP\tPVAL\tPASS\n"


def _row(pos, ref, alt, freq, depth=1000, passed="TRUE"):
    return f"r\t{pos}\t{ref}\t{alt}\t0\t0\t0\t0\t0\t0\t{freq}\t{depth}\t0\t{passed}\n"


def _table(*rows):
    return read_ivar_variants(HEADER + "".join(rows))


# ---------- parsing ----------

def test_only_pass_rows_are_read_and_indels_flagged():
    v = _table(_row(5, "A", "T", 0.3), _row(9, "C", "G", 0.4, passed="FALSE"), _row(12, "C", "+AG", 0.2))
    assert list(v["position"]) == [5, 12]
    assert list(v["is_indel"]) == [False, True]


def test_non_ivar_table_rejected():
    with pytest.raises(ValueError, match="missing column"):
        read_ivar_variants("a\tb\n1\t2\n")


# ---------- per-site statistics ----------

def test_single_alt_at_30_percent():
    t = site_table(_table(_row(5, "A", "T", 0.3)))
    r = t.iloc[0]
    assert r["minor_freq"] == pytest.approx(0.3)
    assert r["entropy_bits"] == pytest.approx(0.8812909, abs=1e-6)  # binary entropy of 0.3
    assert r["heterozygosity"] == pytest.approx(2 * 0.3 * 0.7)
    assert not r["fixed_difference"]


def test_near_fixed_alt_is_a_fixed_difference_with_tiny_minor_freq():
    t = site_table(_table(_row(5, "C", "T", 0.999)))
    r = t.iloc[0]
    assert r["fixed_difference"]
    assert r["minor_freq"] == pytest.approx(0.001)


def test_two_alts_at_one_position_are_combined_not_double_counted():
    t = site_table(_table(_row(7, "A", "T", 0.2), _row(7, "A", "G", 0.1)))
    assert len(t) == 1
    r = t.iloc[0]
    assert r["alt_freq"] == pytest.approx(0.3)
    assert r["minor_freq"] == pytest.approx(0.3)  # reference 0.7 is the major allele
    assert r["heterozygosity"] == pytest.approx(1 - (0.7**2 + 0.2**2 + 0.1**2))


def test_depth_and_snv_filters():
    v = _table(_row(5, "A", "T", 0.3, depth=40), _row(9, "A", "G", 0.3, depth=500), _row(12, "C", "+AG", 0.3, depth=500))
    assert list(site_table(v, min_depth=100)["position"]) == [9]
    assert list(site_table(v, min_depth=100, snv_only=False)["position"]) == [9, 12]


# ---------- summary ----------

def test_summary_at_explicit_threshold():
    v = _table(_row(10, "A", "T", 0.30), _row(20, "A", "G", 0.015), _row(30, "C", "T", 1.0))
    t = site_table(v)
    s = site_summary(t, genome_length=1000, min_freq=0.02)
    assert s["variable_sites"] == 1 and s["fixed_differences"] == 1
    assert s["nucleotide_diversity_estimate"] == pytest.approx(2 * 0.3 * 0.7 / 1000)
    assert site_summary(t, 1000, min_freq=0.01)["variable_sites"] == 2


def test_genome_shorter_than_last_variant_is_rejected():
    with pytest.raises(ValueError):
        site_summary(site_table(_table(_row(500, "A", "T", 0.3))), genome_length=100)


def test_spectrum_bins_and_exclusion_of_fixed_differences():
    t = site_table(_table(_row(1, "A", "T", 0.015), _row(2, "A", "T", 0.03), _row(3, "A", "T", 0.10),
                          _row(4, "A", "T", 0.30), _row(5, "A", "T", 1.0)))
    sp = frequency_spectrum(t)
    assert dict(zip(sp["minor_allele_frequency"], sp["sites"])) == {"<1%": 0, "1-2%": 1, "2-5%": 1, "5-20%": 1, "20-50%": 1}


def test_figures_render():
    t = site_table(_table(_row(100, "A", "T", 0.3), _row(500, "A", "G", 0.015), _row(900, "C", "T", 1.0)))
    assert plot_site_positions(t, 1000, 0.02) is not None
    assert plot_site_spectrum(frequency_spectrum(t)) is not None


# ---------- the real table that exposed the earlier problems ----------

def _independent_counts():
    """Counts straight from the file with the csv module, independent of the package code."""
    rows = [r for r in csv.DictReader(REAL.open(), delimiter="\t") if r["PASS"] == "TRUE"]
    snv = [r for r in rows if r["ALT"][0] not in "+-"]
    f = [float(r["ALT_FREQ"]) for r in snv]
    return len(rows), len(snv), sum(0.02 <= x < 0.5 for x in f), sum(x >= 0.5 for x in f)


def test_real_table_matches_independent_counts():
    n_pass, n_snv, n_variable, n_fixed = _independent_counts()
    v = read_ivar_variants(REAL.read_text())
    assert len(v) == n_pass == 260
    t = site_table(v, min_depth=0)
    assert len(t) == n_snv == 258
    s = site_summary(t, 29903, min_freq=0.02)
    assert s["variable_sites"] == n_variable == 18
    assert s["fixed_differences"] == n_fixed == 9


def test_real_table_default_depth_filter_removes_the_shallow_site():
    """Position 19452 has the highest minor frequency (10.9%) but only 46 reads; depth >= 100 removes it."""
    v = read_ivar_variants(REAL.read_text())
    assert 19452 in set(site_table(v, min_depth=0)["position"])
    assert 19452 not in set(site_table(v, min_depth=100)["position"])


# ---------- CliqueSNV failure detection ----------

FAILED = '{"version": "2.0.3", "settings": {}, "error": "Runtime error", "errorCode": 6, "foundHaplotypes": 1, "haplotypes": []}'
OK = '{"version": "2.0.3", "settings": {}, "error": "none", "foundHaplotypes": 4, "haplotypes": []}'


def test_cliquesnv_runtime_error_is_detected():
    info = read_cliquesnv_json(FAILED)
    assert info["failed"] and info["error_code"] == 6 and info["found"] == 1
    assert "runtime error" in describe_cliquesnv_error(info).lower()


def test_cliquesnv_success_is_not_flagged():
    info = read_cliquesnv_json(OK)
    assert not info["failed"] and info["found"] == 4


def test_cliquesnv_unknown_code_is_still_a_failure():
    info = read_cliquesnv_json('{"error": "Something new", "errorCode": 99, "foundHaplotypes": 0}')
    assert info["failed"]
    assert "unknown" in describe_cliquesnv_error(info)
