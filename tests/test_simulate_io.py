import numpy as np
import pytest

from quasigraph import analyze, parse_fasta, simulate_quasispecies, to_fasta
from quasigraph.metrics import site_entropy


def test_simulation_is_reproducible():
    a = simulate_quasispecies(seed=5)
    b = simulate_quasispecies(seed=5)
    assert a.seqs == b.seqs
    assert np.allclose(a.freqs, b.freqs)


def test_simulation_shape():
    sim = simulate_quasispecies(length=200, n_haplotypes=10, seed=2)
    assert len(sim.seqs) == len(set(sim.seqs)) == 10
    assert all(len(s) == 200 for s in sim.seqs)
    assert sim.freqs.sum() == pytest.approx(1.0)


@pytest.mark.parametrize("seed", range(5))
def test_mutations_concentrate_on_hotspots(seed):
    sim = simulate_quasispecies(length=300, n_haplotypes=15, n_hotspots=8, seed=seed)
    variable = site_entropy(sim.seqs, sim.freqs) > 0
    is_hotspot = np.zeros(300, dtype=bool)
    is_hotspot[sim.hotspots] = True
    hotspot_rate = variable[is_hotspot].mean()  # share of hotspot sites that vary
    background_rate = variable[~is_hotspot].mean()  # share of other sites that vary
    assert hotspot_rate > 5 * background_rate


def test_fasta_round_trip_with_frequencies():
    sim = simulate_quasispecies(length=60, n_haplotypes=5, seed=1)
    names, seqs, freqs = parse_fasta(to_fasta(sim.names, sim.seqs, sim.freqs))
    assert seqs == sim.seqs
    assert np.allclose(freqs, sim.freqs, atol=1e-6)


def test_fasta_without_frequencies():
    names, seqs, freqs = parse_fasta(">a\nACGT\n>b\nACGA\n")
    assert freqs is None and seqs == ["ACGT", "ACGA"]


def test_fasta_mixed_frequencies_rejected():
    with pytest.raises(ValueError):
        parse_fasta(">a freq=0.5\nACGT\n>b\nACGA\n")


def test_analyze_returns_everything():
    sim = simulate_quasispecies(seed=7)
    res = analyze(sim.seqs, sim.freqs)
    assert set(res) == {"summary", "site_entropy", "top_sites", "graph"}
    assert res["summary"].n_haplotypes == 12
    assert len(res["site_entropy"]) == 300
