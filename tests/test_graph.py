import pytest

from quasigraph import build_variation_graph, simulate_quasispecies
from quasigraph.graph import END, START
from quasigraph.plot import draw_graph


def test_blocks_and_bubbles():
    seqs = ["AAAAAAAA", "AATAAAAA", "AATAAGAA"]
    vg = build_variation_graph(seqs, [0.5, 0.3, 0.2])
    # AA | T/A | AAA/AAG... columns: col 2 and col 5 variable, separate runs
    assert vg.n_bubbles == 2
    assert [b[2] for b in vg.blocks] == [False, True, False, True, False]


def test_adjacent_variable_columns_form_one_bubble():
    seqs = ["AAAA", "ATTA", "AGCA"]
    vg = build_variation_graph(seqs)
    assert vg.n_bubbles == 1
    alleles = {n[1] for n in vg.graph.nodes if n not in (START, END) and vg.graph.nodes[n]["variable"]}
    assert alleles == {"AA", "TT", "GC"}


def test_paths_reconstruct_original_haplotypes():
    sim = simulate_quasispecies(length=120, n_haplotypes=8, seed=3)
    vg = build_variation_graph(sim.seqs, sim.freqs)
    for i, h in enumerate(vg.haplotypes):
        assert vg.reconstruct(i) == h


def test_weights_are_conserved():
    sim = simulate_quasispecies(length=120, n_haplotypes=8, seed=4)
    vg = build_variation_graph(sim.seqs, sim.freqs)
    G = vg.graph
    assert sum(d["weight"] for _, _, d in G.out_edges(START, data=True)) == pytest.approx(1.0)
    assert sum(d["weight"] for _, _, d in G.in_edges(END, data=True)) == pytest.approx(1.0)
    for b in range(len(vg.blocks)):
        total = sum(G.nodes[n]["weight"] for n in G.nodes if n not in (START, END) and n[0] == b)
        assert total == pytest.approx(1.0)


def test_identical_haplotypes_make_a_single_path():
    vg = build_variation_graph(["ACGT", "ACGT"])
    assert vg.n_bubbles == 0
    assert vg.stats()["nodes"] == 1


def test_plot_runs():
    sim = simulate_quasispecies(length=100, n_haplotypes=6, seed=1)
    fig = draw_graph(build_variation_graph(sim.seqs, sim.freqs))
    assert fig is not None
