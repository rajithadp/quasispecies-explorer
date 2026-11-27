import numpy as np
import pytest

from quasigraph import site_entropy, summarize, top_variable_sites
from quasigraph.metrics import collapse_haplotypes


def test_two_equal_haplotypes_one_difference():
    seqs = ["AAAAAAAAAA", "AAAAATAAAA"]
    s = summarize(seqs, [0.5, 0.5])
    assert s.n_haplotypes == 2
    assert s.segregating_sites == 1
    assert s.haplotype_entropy_bits == pytest.approx(1.0)
    assert s.effective_haplotypes == pytest.approx(2.0)
    # 2 * 0.5 * 0.5 * 1 difference / 10 sites
    assert s.nucleotide_diversity == pytest.approx(0.05)
    assert s.mean_pairwise_distance == pytest.approx(0.5)


def test_site_entropy_values():
    ent = site_entropy(["AAAA", "AATA"], [0.5, 0.5])
    assert ent[2] == pytest.approx(1.0)
    assert ent[[0, 1, 3]].sum() == 0


def test_invariant_sites_have_exactly_zero_entropy_despite_float_error():
    # frequencies whose partial sums are not exactly 1.0 in floating point
    seqs = ["AAAA", "AAAT", "AAAC"]
    ent = site_entropy(seqs, [0.1, 0.2, 0.7])
    assert (ent[:3] == 0).all()
    assert ent[3] > 0
    assert len(top_variable_sites(seqs, [0.1, 0.2, 0.7])) == 1


def test_identical_sequences_have_no_diversity():
    s = summarize(["ACGT", "ACGT", "ACGT"])
    assert s.n_input_sequences == 3
    assert s.n_haplotypes == 1
    assert s.haplotype_entropy_bits == 0
    assert s.nucleotide_diversity == 0
    assert s.segregating_sites == 0


def test_duplicates_are_merged_and_frequencies_sum():
    haps, f = collapse_haplotypes(["AAA", "AAA", "AAT"], [1, 1, 2])
    assert haps == ["AAA", "AAT"] or haps == ["AAT", "AAA"]
    assert f.sum() == pytest.approx(1.0)
    assert sorted(f) == pytest.approx([0.5, 0.5])


def test_unaligned_input_rejected():
    with pytest.raises(ValueError):
        summarize(["AAAA", "AAA"])


def test_bad_frequencies_rejected():
    with pytest.raises(ValueError):
        summarize(["AA", "AT"], [1.0])
    with pytest.raises(ValueError):
        summarize(["AA", "AT"], [-1.0, 2.0])


def test_top_variable_sites_ranking():
    seqs = ["AAAA", "ATAA", "ATGA"]
    df = top_variable_sites(seqs, [0.4, 0.3, 0.3], k=5)
    assert list(df["position"])[0] in (2, 3)
    assert set(df["position"]) == {2, 3}
    assert (df["entropy_bits"] > 0).all()


def test_gap_counts_as_difference():
    s = summarize(["AAAA", "A-AA"], [0.5, 0.5])
    assert s.segregating_sites == 1
    assert s.mean_pairwise_distance == pytest.approx(0.5)
