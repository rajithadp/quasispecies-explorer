import pytest

from quasigraph.io import parse_haplotypes, to_fasta
from quasigraph.pipeline_utils import (
    Variant,
    build_synthetic_cloud,
    convert_cliquesnv_fasta,
    parse_ivar_variants,
    read_first_fasta_sequence,
)

CLIQUESNV = ">2_fr_0.25\nAATA\n>1_fr_0.70\nAAAA\n>3_fr_0.05\nAATT\n"

IVAR_TSV = (
    "REGION\tPOS\tREF\tALT\tREF_DP\tREF_RV\tREF_QUAL\tALT_DP\tALT_RV\tALT_QUAL\tALT_FREQ\tTOTAL_DP\tPVAL\tPASS\n"
    "ref\t3\tA\tT\t70\t30\t38\t30\t10\t37\t0.30\t100\t0\tTRUE\n"
    "ref\t5\tG\tC\t98\t40\t38\t2\t1\t37\t0.015\t100\t0.2\tTRUE\n"    # below min_freq
    "ref\t7\tT\tA\t60\t30\t38\t40\t20\t37\t0.40\t100\t0\tFALSE\n"   # failed PASS
    "ref\t9\tC\t+AG\t80\t30\t38\t20\t10\t37\t0.20\t100\t0\tTRUE\n"  # insertion
)


def test_cliquesnv_conversion_keeps_frequencies_sorted_and_renamed():
    names, seqs, freqs = convert_cliquesnv_fasta(CLIQUESNV)
    assert names == ["hap_1", "hap_2", "hap_3"]
    assert seqs == ["AAAA", "AATA", "AATT"]
    assert freqs == pytest.approx([0.70, 0.25, 0.05])


def test_converted_fasta_survives_round_trip():
    names, seqs, freqs = convert_cliquesnv_fasta(CLIQUESNV)
    hs = parse_haplotypes(to_fasta(names, seqs, freqs))
    assert hs.seqs == seqs
    assert list(hs.freqs) == pytest.approx(freqs, abs=1e-6)


def test_non_cliquesnv_fasta_is_rejected():
    with pytest.raises(ValueError):
        convert_cliquesnv_fasta(">a\nAAAA\n>b\nAATA\n")


def test_ivar_variants_filtering():
    v = parse_ivar_variants(IVAR_TSV, min_freq=0.02)
    assert [(x.pos, x.alt) for x in v] == [(3, "T"), (9, "+AG")]
    assert v[0].is_indel is False and v[1].is_indel is True


def test_reference_reader_takes_first_record_only():
    assert read_first_fasta_sequence(">a\nacg\nT\n>b\nGGGG\n") == "ACGT"
    with pytest.raises(ValueError):
        read_first_fasta_sequence("")


def test_synthetic_cloud_is_reproducible_and_skips_indels():
    ref = "A" * 50
    vs = [Variant(10, "A", "T", 0.3), Variant(20, "A", "+GG", 0.2)]
    n1, s1, skipped = build_synthetic_cloud(ref, vs, 40, seed=3)
    n2, s2, _ = build_synthetic_cloud(ref, vs, 40, seed=3)
    assert s1 == s2 and skipped == 1
    assert all(len(s) == 50 for s in s1)
    assert all(s[19] == "A" for s in s1)  # the indel was not applied


def test_synthetic_cloud_matches_site_frequency():
    ref = "A" * 30
    _, seqs, _ = build_synthetic_cloud(ref, [Variant(5, "A", "T", 0.25)], 4000, seed=1)
    observed = sum(s[4] == "T" for s in seqs) / len(seqs)
    assert observed == pytest.approx(0.25, abs=0.03)


def test_synthetic_cloud_two_alts_at_one_site_do_not_overwrite_each_other():
    ref = "A" * 10
    vs = [Variant(3, "A", "T", 0.3), Variant(3, "A", "G", 0.2)]
    _, seqs, _ = build_synthetic_cloud(ref, vs, 5000, seed=2)
    t = sum(s[2] == "T" for s in seqs) / 5000
    g = sum(s[2] == "G" for s in seqs) / 5000
    assert t == pytest.approx(0.3, abs=0.03) and g == pytest.approx(0.2, abs=0.03)


def test_variant_outside_reference_rejected():
    with pytest.raises(ValueError):
        build_synthetic_cloud("AAAA", [Variant(9, "A", "T", 0.5)], 3)


def test_synthetic_output_is_marked_synthetic():
    names, seqs, _ = build_synthetic_cloud("A" * 10, [Variant(2, "A", "T", 0.5)], 5, seed=0)
    assert parse_haplotypes(to_fasta(names, seqs, synthetic=True)).synthetic
