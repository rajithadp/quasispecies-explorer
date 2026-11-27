import numpy as np
import pytest

from quasigraph import (
    build_variation_graph,
    mutations_vs_consensus,
    parse_haplotypes,
    rare_haplotype_load,
    simulate_quasispecies,
    site_entropy,
    summarize,
    to_fasta,
    to_gfa,
)
from quasigraph.graph import END, START

# ---------- frequency parsing ----------

CLIQUESNV = ">1_fr_0.6\nAAAAAAAA\n>2_fr_0.3\nAATAAAAA\n>3_fr_0.1\nAATAAGAA\n"


def test_cliquesnv_header_frequencies_are_read():
    hs = parse_haplotypes(CLIQUESNV)
    assert hs.names == ["1_fr_0.6", "2_fr_0.3", "3_fr_0.1"]
    assert list(hs.freqs) == pytest.approx([0.6, 0.3, 0.1])
    assert not hs.weights_inferred and not hs.synthetic


def test_freq_equals_syntax_and_scientific_notation():
    hs = parse_haplotypes(">a freq=1e-3\nAA\n>b freq=0.999\nAT\n")
    assert list(hs.freqs) == pytest.approx([0.001, 0.999])


def test_no_frequencies_is_flagged_as_inferred():
    hs = parse_haplotypes(">a\nAA\n>b\nAT\n")
    assert hs.freqs is None and hs.weights_inferred


def test_mixed_frequencies_rejected():
    with pytest.raises(ValueError):
        parse_haplotypes(">a freq=0.5\nAA\n>b\nAT\n")


def test_synthetic_flag_round_trip():
    text = to_fasta(["h1", "h2"], ["AA", "AT"], [0.5, 0.5], synthetic=True)
    assert parse_haplotypes(text).synthetic is True
    assert parse_haplotypes(to_fasta(["h1"], ["AA"], [1.0])).synthetic is False


def test_weighted_metrics_differ_from_unweighted():
    """The old app treated every record equally; frequencies must change the answer."""
    seqs = ["AAAA", "AAAT"]
    equal = summarize(seqs)
    skewed = summarize(seqs, [0.99, 0.01])
    assert equal.nucleotide_diversity > skewed.nucleotide_diversity
    assert skewed.dominant_haplotype_freq == pytest.approx(0.99)


# ---------- explicit threshold (replaces the hidden 5% filter) ----------

def test_min_minor_freq_filters_only_when_asked():
    seqs = ["AAAA", "AAAT"]
    f = [0.97, 0.03]
    assert site_entropy(seqs, f)[3] > 0  # default: nothing hidden
    assert site_entropy(seqs, f, min_minor_freq=0.05)[3] == 0  # explicit filter
    assert site_entropy(seqs, f, min_minor_freq=0.02)[3] > 0


def test_filter_sums_minor_alleles_not_just_second_allele():
    # two minor alleles of 3% each (6% total) must survive a 5% filter
    seqs = ["AAAA", "AAAT", "AAAC"]
    ent = site_entropy(seqs, [0.94, 0.03, 0.03], min_minor_freq=0.05)
    assert ent[3] > 0


# ---------- rare haplotype load ----------

def test_rare_haplotype_load_both_readings():
    r = rare_haplotype_load([0.80, 0.12, 0.04, 0.03, 0.01], threshold=0.05)
    assert r["n_rare"] == 3
    assert r["mass"] == pytest.approx(0.08)
    assert r["count_fraction"] == pytest.approx(0.6)


def test_rhl_in_summary():
    s = summarize(["AAAA", "AAAT", "AATT"], [0.9, 0.06, 0.04], rhl_threshold=0.05)
    assert s.n_rare_haplotypes == 1
    assert s.rhl_mass == pytest.approx(0.04)


def test_pi_equals_mean_site_heterozygosity():
    """Weighted pairwise diversity must equal the mean over sites of 1 - sum(p^2)."""
    sim = simulate_quasispecies(length=80, n_haplotypes=7, seed=11)
    s = summarize(sim.seqs, sim.freqs)
    A = np.array([list(x) for x in sim.seqs])
    het = []
    for j in range(A.shape[1]):
        p = [sim.freqs[A[:, j] == b].sum() for b in set(A[:, j])]
        het.append(1 - sum(x * x for x in p))
    assert s.nucleotide_diversity == pytest.approx(np.mean(het))


def test_mutations_vs_consensus_is_weighted():
    seqs = ["AAAA", "AATA", "AATA"]
    # unweighted consensus would be AATA; with 80% weight on AAAA it is AAAA
    assert mutations_vs_consensus(seqs, [0.8, 0.1, 0.1]) == [0, 1, 1]
    assert mutations_vs_consensus(seqs) == [1, 0, 0]


# ---------- GFA export ----------

def _parse_gfa(text):
    segs, links, paths = {}, set(), {}
    for line in text.strip().splitlines():
        f = line.split("\t")
        if f[0] == "S":
            segs[f[1]] = f[2]
        elif f[0] == "L":
            links.add((f[1] + f[2], f[3] + f[4]))
        elif f[0] == "P":
            paths[f[1]] = f[2].split(",")
    return segs, links, paths


def test_gfa_paths_reconstruct_haplotypes_and_links_exist():
    sim = simulate_quasispecies(length=150, n_haplotypes=9, seed=6)
    vg = build_variation_graph(sim.seqs, sim.freqs)
    names = [f"h{i}" for i in range(len(vg.haplotypes))]
    segs, links, paths = _parse_gfa(to_gfa(vg, names))
    assert len(segs) == vg.stats()["nodes"]
    for i, name in enumerate(names):
        walk = paths[name]
        assert "".join(segs[s[:-1]] for s in walk) == vg.haplotypes[i]
        for a, b in zip(walk, walk[1:]):
            assert (a, b) in links


def test_gfa_name_count_checked():
    vg = build_variation_graph(["AA", "AT"])
    with pytest.raises(ValueError):
        to_gfa(vg, ["only_one"])


def test_gfa_header_and_default_names():
    vg = build_variation_graph(["AAAA", "AATA"], [0.7, 0.3])
    text = to_gfa(vg)
    assert text.startswith("H\tVN:Z:1.0")
    assert "P\thap1\t" in text and "P\thap2\t" in text
