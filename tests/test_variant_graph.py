import random
from collections import Counter

import networkx as nx
import pytest

from quasigraph import read_ivar_variants
from quasigraph.linkage import Site
from quasigraph.variant_graph import build_variant_graph, to_gfa

HEADER = "REGION\tPOS\tREF\tALT\tREF_DP\tREF_RV\tREF_QUAL\tALT_DP\tALT_RV\tALT_QUAL\tALT_FREQ\tTOTAL_DP\tPVAL\tPASS\n"
random.seed(11)
REF = "".join(random.choice("ACGT") for _ in range(60))


def alt_of(pos, k=0):
    return [b for b in "ACGT" if b != REF[pos - 1]][k]


def table(*rows):
    """rows: (pos, alt, freq[, depth])"""
    body = ""
    for r in rows:
        pos, alt, freq = r[:3]
        depth = r[3] if len(r) > 3 else 1000
        ref = REF[pos - 1] if not alt.startswith(("+", "-")) else REF[pos - 1]
        body += f"r\t{pos}\t{ref}\t{alt}\t0\t0\t0\t0\t0\t0\t{freq}\t{depth}\t0\tTRUE\n"
    return read_ivar_variants(HEADER + body)


def parse_gfa(text):
    segs, links, paths = {}, [], {}
    for line in text.strip().split("\n"):
        f = line.split("\t")
        if f[0] == "S":
            segs[f[1]] = (f[2], dict(t.split(":", 2)[0::2] for t in f[3:]))
        elif f[0] == "L":
            links.append((f[1], f[3], dict(t.split(":", 2)[0::2] for t in f[6:])))
        elif f[0] == "P":
            paths[f[1]] = ([x[:-1] for x in f[2].split(",")], dict(t.split(":", 2)[0::2] for t in f[4:]))
    return segs, links, paths


def consensus_expected(variants):
    """Independent: reference with the most frequent allele at each position."""
    seq = list(REF)
    by = {}
    for r in variants.itertuples():
        by.setdefault(r.position, {})[r.alt] = r.freq
    for pos, alts in by.items():
        ref_f = 1 - sum(alts.values())
        best = max(list(alts.items()) + [(REF[pos - 1], ref_f)], key=lambda kv: kv[1])
        seq[pos - 1] = best[0]
    return "".join(seq)


V = table((10, alt_of(10), 0.03), (12, alt_of(12), 0.05), (12, alt_of(12, 1), 0.04), (30, alt_of(30), 0.95), (59, alt_of(59), 0.02))


def test_reference_path_spells_the_reference_and_consensus_path_the_major_alleles():
    vg = build_variant_graph(REF, V, region="chrT")
    assert vg.sequence(vg.paths["reference"]) == REF
    assert vg.sequence(vg.paths["consensus"]) == consensus_expected(V)
    assert consensus_expected(V)[29] == alt_of(30)          # fixed difference at 30 is the major allele


def test_nodes_edges_and_path_count_follow_the_structure():
    vg = build_variant_graph(REF, V)
    st = vg.stats()
    assert st["bubbles"] == 4 and st["fixed_difference_bubbles"] == 1 and st["minor_alleles"] == 4
    G = nx.DiGraph(vg.edges)
    assert nx.is_directed_acyclic_graph(G)
    sources = [n for n in G if G.in_degree(n) == 0]
    sinks = [n for n in G if G.out_degree(n) == 0]
    assert len(sources) == len(sinks) == 1
    # number of source-to-sink walks = product of the number of alleles per bubble: 2 * 3 * 2 * 2
    assert len(list(nx.all_simple_paths(G, sources[0], sinks[0]))) == 2 * 3 * 2 * 2
    # every walk spells a sequence of the reference length (SNV-only graph)
    for walk in nx.all_simple_paths(G, sources[0], sinks[0]):
        assert len(vg.sequence(walk)) == len(REF)


def test_allele_frequencies_and_roles():
    vg = build_variant_graph(REF, V)
    n = vg.nodes
    assert n[f"a10{REF[9]}"].freq == pytest.approx(0.97) and n[f"a10{REF[9]}"].role == "major"
    assert n[f"a10{alt_of(10)}"].freq == pytest.approx(0.03) and n[f"a10{alt_of(10)}"].role == "minor"
    assert n[f"a12{REF[11]}"].freq == pytest.approx(0.91)                       # 1 - 0.05 - 0.04
    assert n[f"a30{alt_of(30)}"].role == "major" and n[f"a30{REF[29]}"].role == "minor"
    assert n[f"a30{REF[29]}"].matches_reference and not n[f"a30{alt_of(30)}"].matches_reference


def test_thresholds_depth_indels_and_reference_mismatch_are_counted():
    rows = [(10, alt_of(10), 0.019), (20, alt_of(20), 0.5, 50), (25, "+AG", 0.2), (40, alt_of(40), 0.02, 100)]
    v = table(*rows)
    bad = read_ivar_variants(HEADER + f"r\t45\t{alt_of(45)}\t{alt_of(45, 1)}\t0\t0\t0\t0\t0\t0\t0.3\t500\t0\tTRUE\n")  # wrong REF base
    import pandas as pd
    vg = build_variant_graph(REF, pd.concat([v, bad], ignore_index=True), min_freq=0.02, min_depth=100)
    assert vg.stats()["bubbles"] == 1 and f"a40{alt_of(40)}" in vg.nodes          # exactly min_freq and min_depth are kept
    assert vg.skipped["indels"] == 1 and vg.skipped["below_min_depth"] == 1
    assert vg.skipped["reference_mismatch"] == 1 and vg.skipped["no_alt_at_min_freq"] == 1


def test_variants_at_the_first_and_last_base_and_adjacent_sites():
    v = table((1, alt_of(1), 0.1), (2, alt_of(2), 0.1), (60, alt_of(60), 0.1))
    vg = build_variant_graph(REF, v)
    assert vg.sequence(vg.paths["reference"]) == REF
    assert vg.columns[0] == [f"a1{REF[0]}", f"a1{alt_of(1)}"] and vg.columns[-1][0] == f"a60{REF[59]}"
    assert not any(n.startswith("b") and n.endswith("_0") for n in vg.nodes)       # no empty backbone nodes
    assert (f"a1{REF[0]}", f"a2{alt_of(2)}") in vg.edges                           # adjacent bubbles are joined directly


def test_no_variants_gives_a_single_backbone_node():
    vg = build_variant_graph(REF, read_ivar_variants(HEADER))
    assert list(vg.nodes) == ["b1_60"] and vg.edges == [] and vg.sequence(vg.paths["consensus"]) == REF


# ---------- read-supported evidence paths ----------

def sites(*pos):
    return [Site(p, REF[p - 1], alt_of(p)) for p in pos]


def test_evidence_paths_spell_the_chosen_alleles_and_carry_fragment_counts():
    v = table((10, alt_of(10), 0.03), (12, alt_of(12), 0.03), (14, alt_of(14), 0.03))
    ev = [(sites(10, 12, 14), Counter({"AAA": 8, "RRR": 300, "-RR": 3, "ARA": 1, "ROA": 2}))]
    vg = build_variant_graph(REF, v, evidence=ev)
    by = {e.name: e for e in vg.evidence}
    assert set(by) == {"linked_10_12_14_AAA", "linked_10_12_14_RRR", "linked_10_12_14_ARA"}   # partial / 'O' patterns skipped
    assert by["linked_10_12_14_AAA"].fragments == 8
    aaa = vg.sequence(by["linked_10_12_14_AAA"].nodes)
    assert aaa == alt_of(10) + REF[10] + alt_of(12) + REF[12] + alt_of(14)           # spans sites 10..14 only
    rrr = vg.sequence(by["linked_10_12_14_RRR"].nodes)
    assert rrr == REF[9:14]
    ara = vg.sequence(by["linked_10_12_14_ARA"].nodes)
    assert ara == alt_of(10) + REF[10] + REF[11] + REF[12] + alt_of(14)


def test_every_evidence_path_is_a_walk_in_the_graph():
    v = table((10, alt_of(10), 0.03), (12, alt_of(12), 0.03), (13, alt_of(13), 0.04))
    vg = build_variant_graph(REF, v, evidence=[(sites(10, 12, 13), Counter({"AAA": 5, "RRR": 9, "ARR": 2}))])
    edges = set(vg.edges)
    assert vg.evidence
    for e in vg.evidence:
        assert all((a, b) in edges for a, b in zip(e.nodes, e.nodes[1:]))


def test_group_with_an_unlisted_variable_site_between_is_not_placed():
    v = table((10, alt_of(10), 0.03), (12, alt_of(12), 0.03), (14, alt_of(14), 0.03))
    vg = build_variant_graph(REF, v, evidence=[(sites(10, 14), Counter({"AA": 4}))])
    assert vg.evidence == [] and len(vg.unplaced) == 1 and "[12]" in vg.unplaced[0][1]


def test_group_with_a_site_that_is_not_a_bubble_is_not_placed():
    v = table((10, alt_of(10), 0.03))
    vg = build_variant_graph(REF, v, evidence=[(sites(10, 12), Counter({"AA": 4}))])
    assert vg.evidence == [] and "not bubbles" in vg.unplaced[0][1]


# ---------- GFA ----------

def test_gfa_round_trip_with_an_independent_parser():
    v = table((10, alt_of(10), 0.03), (12, alt_of(12), 0.05), (12, alt_of(12, 1), 0.04), (30, alt_of(30), 0.95))
    vg = build_variant_graph(REF, v, region="chrT", evidence=[(sites(10, 12), Counter({"AA": 6, "RR": 40}))])
    segs, links, paths = parse_gfa(to_gfa(vg))
    assert len(segs) == len(vg.nodes) and len(links) == len(vg.edges)
    link_set = {(a, b) for a, b, _ in links}
    for name, (walk, tags) in paths.items():
        assert all((a, b) in link_set for a, b in zip(walk, walk[1:])), name
    spell = lambda name: "".join(segs[s][0] for s in paths[name][0])
    assert spell("reference") == REF and spell("consensus") == consensus_expected(v)
    assert paths["linked_10_12_AA"][1] == {"fc": "6"} and paths["linked_10_12_RR"][1] == {"fc": "40"}
    # segment tags: 0-based offset, length, and allele frequency/role
    tags = segs["a10" + alt_of(10)][1]
    assert tags["SO"] == "9" and tags["LN"] == "1" and float(tags["FR"]) == pytest.approx(0.03) and tags["rl"] == "minor"
    assert segs["a30" + REF[29]][1]["mr"] == "1" and segs["a30" + alt_of(30)][1]["mr"] == "0"
    # links into an allele carry its frequency, links into backbone carry none
    for a, b, t in links:
        assert ("FR" in t) == (b.startswith("a"))
        if "FR" in t:
            assert float(t["FR"]) == pytest.approx(float(segs[b][1]["FR"]))


def test_gfa_sequence_lengths_match_ln_and_offsets_tile_the_reference():
    v = table((10, alt_of(10), 0.03), (30, alt_of(30), 0.4))
    segs, _, paths = parse_gfa(to_gfa(build_variant_graph(REF, v)))
    pos = 0
    for s in paths["reference"][0]:
        seq, t = segs[s]
        assert int(t["SO"]) == pos and int(t["LN"]) == len(seq)
        pos += len(seq)
    assert pos == len(REF)


def test_random_variant_sets_always_give_consistent_graphs():
    rng = random.Random(3)
    for _ in range(40):
        positions = sorted(rng.sample(range(1, 61), rng.randint(0, 12)))
        rows = [(p, alt_of(p, rng.randint(0, 2)), round(rng.uniform(0.02, 0.99), 3)) for p in positions]
        v = table(*rows) if rows else read_ivar_variants(HEADER)
        # one alt per position here: frequencies up to 0.99 are allowed
        vg = build_variant_graph(REF, v)
        assert vg.sequence(vg.paths["reference"]) == REF
        if rows:
            assert vg.sequence(vg.paths["consensus"]) == consensus_expected(v)
        G = nx.DiGraph(vg.edges)
        if rows:
            assert nx.is_directed_acyclic_graph(G)
        segs, links, paths = parse_gfa(to_gfa(vg))
        assert {(a, b) for a, b, _ in links} == set(vg.edges)


# ---------- cases found by mutation testing ----------

def raw(pos, ref, alt, freq=0.3, depth=500):
    return f"r\t{pos}\t{ref}\t{alt}\t0\t0\t0\t0\t0\t0\t{freq}\t{depth}\t0\tTRUE\n"


def test_a_fifty_fifty_site_counts_as_a_fixed_difference_like_in_the_site_level_analysis():
    vg = build_variant_graph(REF, table((20, alt_of(20), 0.5)))
    assert vg.nodes[f"a20{alt_of(20)}"].role == "major" and vg.nodes[f"a20{REF[19]}"].role == "minor"
    assert vg.stats()["fixed_difference_bubbles"] == 1


def test_wrong_length_patterns_are_ignored():
    v = table((10, alt_of(10), 0.03), (12, alt_of(12), 0.03))
    vg = build_variant_graph(REF, v, evidence=[(sites(10, 12), Counter({"A": 3, "AAA": 2, "AA": 4}))])
    assert [e.name for e in vg.evidence] == ["linked_10_12_AA"]


def test_unusable_alleles_are_skipped_not_added():
    text = HEADER + raw(10, REF[9], "N") + raw(15, REF[14], alt_of(15) + alt_of(15, 1)) \
        + raw(20, REF[19], REF[19] + "C") + raw(61, "A", "C") + raw(0, REF[-1], "C" if REF[-1] != "C" else "A") + raw(30, REF[29], alt_of(30))
    vg = build_variant_graph(REF, read_ivar_variants(text))
    assert vg.stats()["bubbles"] == 1 and f"a30{alt_of(30)}" in vg.nodes
    assert vg.skipped["other_base"] + vg.skipped["indels"] >= 3 and vg.skipped["reference_mismatch"] == 2


# ---------- command line ----------

def test_cli_builds_graph_with_linkage_evidence_from_a_sam(tmp_path, capsys):
    import build_variant_graph
    (tmp_path / "ref.fa").write_text(">chrT desc\n" + REF + "\n")
    v = HEADER + raw(10, REF[9], alt_of(10), 0.1) + raw(12, REF[11], alt_of(12), 0.1)
    (tmp_path / "v.tsv").write_text(v)
    reads = ["@HD\tVN:1.6"]
    def rd(name, seq_alleles):
        seq = list(REF[0:40])
        for p, b in seq_alleles.items():
            seq[p - 1] = b
        return "\t".join([name, "0", "chrT", "1", "60", "40M", "*", "0", "0", "".join(seq), "I" * 40])
    for i in range(4):
        reads.append(rd(f"alt{i}", {10: alt_of(10), 12: alt_of(12)}))
    for i in range(6):
        reads.append(rd(f"ref{i}", {}))
    (tmp_path / "x.sam").write_text("\n".join(reads) + "\n")
    out = str(tmp_path / "g")
    rc = build_variant_graph.main(["--variants", str(tmp_path / "v.tsv"), "--ref", str(tmp_path / "ref.fa"),
                                   "--alignment", str(tmp_path / "x.sam"), "--sites", "10,12", "--out", out])
    assert rc == 0
    segs, links, paths = parse_gfa((tmp_path / "g.gfa").read_text())
    assert paths["linked_10_12_AA"][1] == {"fc": "4"} and paths["linked_10_12_RR"][1] == {"fc": "6"}
    assert segs["b1_9"][1]["SN"] == "chrT"
    assert "2 variable positions" in capsys.readouterr().out


def test_cli_errors_are_reported(tmp_path, capsys):
    import build_variant_graph
    (tmp_path / "ref.fa").write_text(">r\n" + REF + "\n")
    (tmp_path / "v.tsv").write_text(HEADER + raw(10, REF[9], alt_of(10), 0.1))
    base = ["--variants", str(tmp_path / "v.tsv"), "--ref", str(tmp_path / "ref.fa"), "--out", str(tmp_path / "g")]
    assert build_variant_graph.main(base) == 0 and (tmp_path / "g.gfa").exists()          # no groups: plain graph
    assert build_variant_graph.main(base + ["--auto"]) == 0                               # nothing to group: fine
    assert build_variant_graph.main(base + ["--sites", "10,11"]) == 2                     # 11 is not a variant; also no alignment
    assert "error" in capsys.readouterr().err
    (tmp_path / "v2.tsv").write_text(HEADER + raw(10, REF[9], alt_of(10), 0.1) + raw(12, REF[11], alt_of(12), 0.1))
    assert build_variant_graph.main(["--variants", str(tmp_path / "v2.tsv"), "--ref", str(tmp_path / "ref.fa"), "--sites", "10,12"]) == 2
    assert "--alignment" in capsys.readouterr().err


def test_gfa_marks_minor_gt_ca_alleles_as_artifact_prone():
    small = "AGCATGCATC"                      # 1-based: G at 2, C at 3, T at 5, G at 6, C at 7
    def row(pos, alt, freq):
        return raw(pos, small[pos - 1], alt, freq)
    v = read_ivar_variants(HEADER + row(2, "T", 0.03)        # G>T minor    -> flagged
                           + row(3, "A", 0.03)               # C>A minor    -> flagged
                           + row(5, "C", 0.03)               # T>C minor    -> not flagged
                           + row(6, "T", 0.6)                # G>T but major (fixed difference) -> not flagged
                           + row(7, "G", 0.03))              # C>G minor    -> not flagged
    segs, _, _ = parse_gfa(to_gfa(build_variant_graph(small, v)))
    flagged = sorted(k for k, (_, t) in segs.items() if t.get("ox") == "1")
    assert flagged == ["a2T", "a3A"]


# ---------- drawing ----------

def test_columns_in_region_and_layout():
    from quasigraph.variant_graph_plot import columns_in_region, layout_region
    v = table((10, alt_of(10), 0.03), (12, alt_of(12), 0.05), (12, alt_of(12, 1), 0.04), (30, alt_of(30), 0.2))
    vg = build_variant_graph(REF, v)
    # columns: b1_9, a10, b11_11, a12, b13_29, a30, b31_60
    assert columns_in_region(vg, 10, 12) == [1, 2, 3]
    assert columns_in_region(vg, 5, 5) == [0]                          # inside a backbone stretch
    assert columns_in_region(vg, 12, 30) == [3, 4, 5]                  # window edges inclusive
    assert columns_in_region(vg, 61, 70) == []
    with pytest.raises(ValueError):
        columns_in_region(vg, 12, 10)
    pos = layout_region(vg, 9, 13)
    assert [pos[n][0] for n in ("b1_9", f"a10{REF[9]}", "b11_11", f"a12{REF[11]}", "b13_29")] == [0, 1, 2, 3, 4]
    assert pos[f"a10{REF[9]}"][1] == 0 and pos[f"a10{alt_of(10)}"][1] == 1          # reference on the middle line
    # three alleles at position 12: the two alternatives sit on opposite sides, the more frequent one above
    assert pos[f"a12{alt_of(12)}"][1] == 1 and pos[f"a12{alt_of(12, 1)}"][1] == -1
    assert set(pos) == {n for ci in (0, 1, 2, 3, 4) for n in vg.columns[ci]}


def test_plot_draws_nodes_paths_and_legend():
    from quasigraph.variant_graph_plot import plot_variant_graph
    v = table((10, alt_of(10), 0.03), (12, alt_of(12), 0.03))
    vg = build_variant_graph(REF, v, evidence=[(sites(10, 12), Counter({"AA": 6, "RR": 40}))])
    fig = plot_variant_graph(vg, 8, 14)
    ax = fig.axes[0]
    labels = [t.get_text() for t in ax.texts]
    assert any(f"10{alt_of(10)}  3.0%" in t for t in labels) and any(f"12{REF[11]}  97.0%" in t for t in labels)
    legend = [t.get_text() for t in ax.get_legend().get_texts()]
    assert sorted(legend) == ["10_12_AA: 6 fragments", "10_12_RR: 40 fragments"]
    assert len(ax.collections) == 4                                    # one marker per allele node
    with pytest.raises(ValueError, match="no graph columns"):
        plot_variant_graph(vg, 100, 200)


def test_plot_marks_artifact_prone_alleles():
    from quasigraph.variant_graph_plot import plot_variant_graph
    small = "AGCATGCATC"
    v = read_ivar_variants(HEADER + raw(2, "G", "T", 0.03) + raw(5, "T", "C", 0.03))
    ax = plot_variant_graph(build_variant_graph(small, v), 1, 10).axes[0]
    labels = [t.get_text() for t in ax.texts]
    assert sum(t.endswith("ox") for t in labels) == 1 and any(t.startswith("2T") and t.endswith("ox") for t in labels)


def test_cli_plot_option(tmp_path, capsys):
    import build_variant_graph
    (tmp_path / "ref.fa").write_text(">chrT\n" + REF + "\n")
    (tmp_path / "v.tsv").write_text(HEADER + raw(10, REF[9], alt_of(10), 0.1))
    base = ["--variants", str(tmp_path / "v.tsv"), "--ref", str(tmp_path / "ref.fa"), "--out", str(tmp_path / "g")]
    assert build_variant_graph.main(base + ["--plot", "5-20"]) == 0
    assert (tmp_path / "g.png").stat().st_size > 1000
    assert build_variant_graph.main(base + ["--plot", "abc"]) == 2
    assert build_variant_graph.main(base + ["--plot", "500-600"]) == 2
    assert "no graph columns" in capsys.readouterr().err


def test_partial_window_draws_only_edges_inside_it_and_marker_area_follows_frequency():
    from quasigraph.variant_graph_plot import plot_variant_graph
    v = table((10, alt_of(10), 0.03), (12, alt_of(12), 0.3))
    vg = build_variant_graph(REF, v, evidence=[(sites(10, 12), Counter({"AA": 6, "RR": 40}))])
    ax = plot_variant_graph(vg, 12, 12).axes[0]                         # only the bubble at 12 (one column)
    sizes = sorted(float(c.get_sizes()[0]) for c in ax.collections)
    assert len(sizes) == 2 and sizes[0] < sizes[1]                      # 30% allele drawn smaller than the 70% one
    assert ax.get_legend() is None                                      # an evidence path with one visible node is not drawn
    ax2 = plot_variant_graph(vg, 10, 12).axes[0]                        # three columns: a10, b11_11, a12
    assert ax2.get_legend() is not None


def test_reference_paths_are_always_teal_and_alternative_paths_are_not():
    from quasigraph.variant_graph_plot import ALT_PATHS, REF_PATH, plot_variant_graph
    v = table((10, alt_of(10), 0.03), (12, alt_of(12), 0.03), (30, alt_of(30), 0.03), (32, alt_of(32), 0.03))
    ev = [(sites(10, 12), Counter({"AA": 5, "RR": 90})), (sites(30, 32), Counter({"AA": 4, "RR": 70, "AR": 2}))]
    ax = plot_variant_graph(build_variant_graph(REF, v, evidence=ev), 8, 34).axes[0]
    by_label = {l.get_label(): l.get_color() for l in ax.lines if l.get_label().endswith("fragments")}
    assert by_label["10_12_RR: 90 fragments"] == by_label["30_32_RR: 70 fragments"] == REF_PATH
    alt_colors = [c for l, c in by_label.items() if "_RR" not in l]
    assert len(alt_colors) == 3 and REF_PATH not in alt_colors and len(set(alt_colors)) == 3 and set(alt_colors) <= set(ALT_PATHS)


# ---------- annotation on the graph ----------

def _annotated(evidence=True):
    from quasigraph.annotation import Annotator, Cds
    from quasigraph.graph_annotation import annotate_graph
    v = table((10, alt_of(10), 0.1), (12, alt_of(12), 0.1), (30, alt_of(30), 0.95), (59, alt_of(59), 0.05))
    ev = [(sites(10, 12), Counter({"AA": 6, "RR": 40}))] if evidence else None
    vg = build_variant_graph(REF, v, region="chrT", evidence=ev)
    ann = Annotator([Cds("c", "A", "", "+", ((1, 60),))], REF)
    return annotate_graph(vg, ann), ann


def _aa(codon):
    from quasigraph.annotation import CODON_TABLE
    return CODON_TABLE[codon]


def test_alleles_carry_gene_amino_acid_and_effect_of_the_allele_alone():
    vg, _ = _annotated()
    c4 = REF[9:12]
    alone = _aa(alt_of(10) + c4[1:])
    expect = f"A:{_aa(c4)}4{alone}:" + ("synonymous" if alone == _aa(c4) else "stop_gained" if alone == "*" else "missense")
    if _aa(c4) == "*":
        pytest.skip("random reference puts a stop codon here")
    assert vg.nodes["a10" + alt_of(10)].annotation == expect


def test_reference_alleles_and_backbone_are_not_annotated():
    vg, _ = _annotated()
    ref_alleles = [n for n in vg.nodes.values() if n.kind == "allele" and n.matches_reference]
    assert ref_alleles and all(n.annotation == "" for n in ref_alleles)
    assert all(n.annotation == "" for n in vg.nodes.values() if n.kind != "allele")
    assert all(n.annotation for n in vg.nodes.values() if n.kind == "allele" and not n.matches_reference)


def test_joint_label_combines_linked_changes_in_one_codon():
    vg, _ = _annotated()
    (e,) = [x for x in vg.evidence if x.name.endswith("_AA")]
    assert [x.joint for x in vg.evidence if x.name.endswith("_RR")] == [""]
    c4 = REF[9:12]
    joint = alt_of(10) + c4[1] + alt_of(12)
    assert e.joint == f"A:{_aa(c4)}4{_aa(joint)}"


def test_joint_is_empty_for_a_path_without_alternative_alleles():
    vg, ann = _annotated()
    from quasigraph.graph_annotation import annotate_graph
    for e in vg.evidence:
        e.joint = "stale"
        e.nodes = [n for n in e.nodes if vg.nodes[n].matches_reference]
    annotate_graph(vg, ann)
    assert all(e.joint == "" for e in vg.evidence)


def test_noncoding_alleles_are_labelled_noncoding():
    from quasigraph.annotation import Annotator, Cds
    from quasigraph.graph_annotation import annotate_graph
    v = table((10, alt_of(10), 0.1), (59, alt_of(59), 0.1))
    vg = build_variant_graph(REF, v)
    annotate_graph(vg, Annotator([Cds("c", "A", "", "+", ((1, 30),))], REF))
    assert vg.nodes["a59" + alt_of(59)].annotation == "noncoding"
    assert vg.nodes["a10" + alt_of(10)].annotation.startswith("A:")


def test_annotating_with_a_different_reference_is_refused():
    from quasigraph.annotation import Annotator, Cds
    from quasigraph.graph_annotation import annotate_graph
    vg = build_variant_graph(REF, table((10, alt_of(10), 0.1)))
    other = REF[:-1] + ("A" if REF[-1] != "A" else "C")
    with pytest.raises(ValueError, match="reference"):
        annotate_graph(vg, Annotator([Cds("c", "A", "", "+", ((1, 30),))], other))


def test_frameshift_join_gives_both_codon_contexts():
    from quasigraph.annotation import Annotator, Cds
    from quasigraph.graph_annotation import annotate_graph
    vg = build_variant_graph(REF, table((13, alt_of(13), 0.1)))
    # two segments sharing base 13 (a -1 ribosomal frameshift join): 13 is the last base of one and first of the next
    annotate_graph(vg, Annotator([Cds("c", "F", "", "+", ((1, 13), (13, 30)))], REF))
    parts = vg.nodes["a13" + alt_of(13)].annotation.split(" | ")
    assert len(parts) == 2 and all(p.startswith("F:") for p in parts)


def test_gfa_has_an_and_jn_tags_only_where_annotated():
    vg, _ = _annotated()
    segs, _, paths = parse_gfa(to_gfa(vg))
    assert segs["a10" + alt_of(10)][1]["an"] == vg.nodes["a10" + alt_of(10)].annotation
    assert not any("an" in t for k, (_, t) in segs.items() if vg.nodes[k].annotation == "")
    aa = next(e for e in vg.evidence if e.name.endswith("_AA"))
    assert paths[aa.name][1]["jn"] == aa.joint and aa.joint
    assert "jn" not in paths[aa.name.replace("_AA", "_RR")][1] and "jn" not in paths["reference"][1]
    plain = build_variant_graph(REF, table((10, alt_of(10), 0.1)))
    assert "an:Z:" not in to_gfa(plain) and "jn:Z:" not in to_gfa(plain)


def test_plot_shows_annotation_under_allele_and_joint_in_legend():
    import matplotlib
    matplotlib.use("Agg")
    from quasigraph.variant_graph_plot import plot_variant_graph
    vg, _ = _annotated()
    ax = plot_variant_graph(vg, 5, 35).axes[0]
    texts = [t.get_text() for t in ax.texts]
    first = vg.nodes["a10" + alt_of(10)].annotation.split(":")
    assert any(t.startswith(f"10{alt_of(10)}") and t.endswith(f"\nA {first[1]}") for t in texts)
    assert not any("\n" in t for t in texts if t.startswith(f"10{REF[9]}"))   # reference allele: one line
    aa = next(e for e in vg.evidence if e.name.endswith("_AA"))
    legend = [t.get_text() for t in ax.get_legend().get_texts()]
    assert f"10_12_AA: 6 fragments ({aa.joint})" in legend and "10_12_RR: 40 fragments" in legend


def test_cli_gff_option_prints_joint_effect_and_writes_tags(tmp_path, capsys):
    import build_variant_graph
    (tmp_path / "ref.fa").write_text(">chrT\n" + REF + "\n")
    (tmp_path / "v.tsv").write_text(HEADER + raw(10, REF[9], alt_of(10), 0.1))
    (tmp_path / "a.gff3").write_text("##gff-version 3\nchrT\tx\tCDS\t1\t60\t.\t+\t0\tID=c;gene=A\n")
    base = ["--variants", str(tmp_path / "v.tsv"), "--ref", str(tmp_path / "ref.fa"), "--out", str(tmp_path / "g")]
    assert build_variant_graph.main(base + ["--gff", str(tmp_path / "a.gff3")]) == 0
    assert "an:Z:A:" in (tmp_path / "g.gfa").read_text()
    (tmp_path / "e.gff3").write_text("##gff-version 3\n")
    assert build_variant_graph.main(base + ["--gff", str(tmp_path / "e.gff3")]) == 2
    assert "no CDS" in capsys.readouterr().err


def test_joint_label_lists_changes_in_different_codons_separated_by_commas():
    from quasigraph.annotation import Annotator, Cds
    from quasigraph.graph_annotation import annotate_graph
    v = table((10, alt_of(10), 0.1), (30, alt_of(30), 0.1))
    ev = [([Site(10, REF[9], alt_of(10)), Site(30, REF[29], alt_of(30))], Counter({"AA": 6, "RR": 40}))]
    vg = build_variant_graph(REF, v, evidence=ev)
    annotate_graph(vg, Annotator([Cds("c", "A", "", "+", ((1, 60),))], REF))
    aa = next(e for e in vg.evidence if e.name.endswith("_AA"))
    parts = aa.joint.split(",")
    c4, c10 = REF[9:12], REF[27:30]
    assert parts == [f"A:{_aa(c4)}4{_aa(alt_of(10) + c4[1:])}", f"A:{_aa(c10)}10{_aa(c10[:2] + alt_of(30))}"]


def test_plot_second_line_for_annotation_without_amino_acid_shows_gene_only():
    import matplotlib
    matplotlib.use("Agg")
    from quasigraph.variant_graph_plot import plot_variant_graph
    vg = build_variant_graph(REF, table((10, alt_of(10), 0.1)))
    vg.nodes["a10" + alt_of(10)].annotation = "A:incomplete_codon"
    texts = [t.get_text() for t in plot_variant_graph(vg, 5, 15).axes[0].texts]
    assert any(t.startswith(f"10{alt_of(10)}") and t.endswith("\nA") for t in texts)
