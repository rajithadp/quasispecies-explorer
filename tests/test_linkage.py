from pathlib import Path

import pytest

import linkage_check
from quasigraph import read_ivar_variants, site_table
from quasigraph.linkage import (
    Site,
    fragment_summary,
    alignment_lines,
    analyze_group,
    auto_groups,
    bases_at_sites,
    best_snv_per_position,
    collect_fragments,
    fisher_exact_2x2,
    pair_verdict,
    parse_cigar,
    parse_site_spec,
)

REAL = Path(__file__).parent / "data" / "real_variants.tsv"
HEADER = "REGION\tPOS\tREF\tALT\tREF_DP\tREF_RV\tREF_QUAL\tALT_DP\tALT_RV\tALT_QUAL\tALT_FREQ\tTOTAL_DP\tPVAL\tPASS\n"


def sam(name, pos, seq, cigar=None, flag=0, mapq=60, rname="ref", qual=None, pnext=0, tlen=0):
    cigar = cigar or f"{len(seq)}M"
    qual = qual or "I" * len(seq)
    return "\t".join([name, str(flag), rname, str(pos), str(mapq), cigar, "=" if pnext else "*", str(pnext), str(tlen), seq, qual])


def read(name, start, length, alts, **kw):
    seq = ["A"] * length
    for p in alts:
        seq[p - start] = "T"
    return sam(name, start, "".join(seq), **kw)


def pair_reads(prefix, start, counts, pos1, pos2, length=100):
    """counts = (alt-alt, alt1-only, alt2-only, ref-ref) single reads covering both positions."""
    lines, k = [], 0
    for n, a1, a2 in zip(counts, (1, 1, 0, 0), (1, 0, 1, 0)):
        for _ in range(n):
            lines.append(read(f"{prefix}{k}", start, length, [p for p, a in ((pos1, a1), (pos2, a2)) if a]))
            k += 1
    return lines


# ---------- CIGAR-aware base extraction ----------

def test_parse_cigar():
    assert parse_cigar("5S10M2I3D20M") == [(5, "S"), (10, "M"), (2, "I"), (3, "D"), (20, "M")]


def test_insertion_shifts_read_index_not_reference():
    seq = "A" * 10 + "GG" + "C" * 10  # 10M 2I 10M starting at 100: C covers reference 110..119
    out = bases_at_sites(100, parse_cigar("10M2I10M"), seq, "*", [105, 112, 119])
    assert out == {105: "A", 112: "C", 119: "C"}


def test_deletion_positions_are_unobserved():
    seq = "A" * 10 + "T" * 10  # 10M 3D 10M: reference 110..112 deleted
    out = bases_at_sites(100, parse_cigar("10M3D10M"), seq, "*", [111, 113, 122])
    assert out == {113: "T", 122: "T"}


def test_soft_clip_is_skipped():
    seq = "GGGGG" + "A" * 10
    out = bases_at_sites(100, parse_cigar("5S10M"), seq, "*", [99, 100])
    assert out == {100: "A"}


def test_site_just_past_the_aligned_block_is_not_read_from_the_clipped_tail():
    seq = "A" * 10 + "GGGGG"  # 10M5S at 100: reference 100..109 aligned, 110 is NOT covered
    assert bases_at_sites(100, parse_cigar("10M5S"), seq, "*", [109, 110]) == {109: "A"}


def test_low_base_quality_is_dropped():
    seq = "AAAAT" + "A" * 5
    qual = "IIII#" + "I" * 5  # '#' = Q2 at the variant base
    assert bases_at_sites(100, parse_cigar("10M"), seq, qual, [104], min_baseq=20) == {}
    assert bases_at_sites(100, parse_cigar("10M"), seq, qual, [104], min_baseq=0) == {104: "T"}


# ---------- fragments ----------

def test_mates_of_one_pair_are_merged_across_a_long_gap():
    lines = [read("p1", 250, 100, [300]), read("p1", 650, 100, [700], flag=147)]
    frags = collect_fragments(lines, [300, 700])
    assert frags == {"p1": {300: "T", 700: "T"}}


def test_mates_that_disagree_become_ambiguous():
    lines = [read("p1", 100, 100, [150]), read("p1", 120, 100, [])]  # both cover 150, different bases
    assert collect_fragments(lines, [150]) == {"p1": {150: None}}


def test_unmapped_secondary_other_contig_and_low_mapq_are_skipped():
    lines = [
        read("a", 100, 100, [150], flag=4),
        read("b", 100, 100, [150], flag=256),
        read("c", 100, 100, [150], rname="other"),
        read("d", 100, 100, [150], mapq=5),
        read("e", 100, 100, [150]),
    ]
    frags = collect_fragments(lines, [150], region="ref", min_mapq=10)
    assert list(frags) == ["e"]


def test_header_lines_and_reads_not_overlapping_sites_are_ignored():
    lines = ["@SQ\tSN:ref\tLN:1000", read("far", 500, 100, []), read("near", 100, 100, [150])]
    assert list(collect_fragments(lines, [150])) == ["near"]


# ---------- statistics ----------

def test_fisher_matches_known_values():
    assert fisher_exact_2x2(8, 2, 1, 5) == pytest.approx(0.034965, abs=1e-6)
    assert fisher_exact_2x2(3, 1, 1, 3) == pytest.approx(0.485714, abs=1e-6)
    assert fisher_exact_2x2(0, 0, 0, 0) == 1.0
    assert fisher_exact_2x2(12, 12, 12, 12) == pytest.approx(1.0)


# ---------- verdicts on synthetic reads with known structure ----------

G = [Site(120, "A", "T"), Site(150, "A", "T")]


def _analyze(counts):
    frags = collect_fragments(pair_reads("r", 100, counts, 120, 150), [120, 150])
    return analyze_group(frags, G)


def test_linked_alleles_are_called_linked():
    res = _analyze((30, 0, 0, 970))
    p = res["pairs"][0]
    assert (p.aa, p.ar, p.ra, p.rr) == (30, 0, 0, 970)
    assert p.verdict.startswith("LINKED")
    assert p.fisher_p < 1e-10
    assert res["patterns"]["AA"] == 30 and res["patterns"]["RR"] == 970


def test_independent_alleles_are_not_called_linked():
    p = _analyze((3, 47, 47, 903))["pairs"][0]  # expected co-occurrences: 2.5
    assert p.verdict.startswith("consistent with independent")


def test_mutually_exclusive_needs_enough_expected_coincidences():
    assert _analyze((0, 200, 200, 600))["pairs"][0].verdict.startswith("MUTUALLY EXCLUSIVE")
    assert "too rare" in _analyze((0, 60, 60, 880))["pairs"][0].verdict


def test_nested_alleles():
    assert _analyze((40, 0, 20, 940))["pairs"][0].verdict.startswith("NESTED")


def test_no_shared_coverage_is_reported_not_guessed():
    frags = collect_fragments(read("x", 100, 100, [120]).split("\n") + read("y", 300, 100, [350]).split("\n"), [120, 350])
    assert "no fragment covers both" in analyze_group(frags, [Site(120, "A", "T"), Site(350, "A", "T")])["pairs"][0].verdict


def test_pair_verdict_directly():
    assert "no alternative allele" in pair_verdict(0, 0, 0, 100, 1.0)
    assert pair_verdict(0, 0, 0, 0, 1.0).startswith("no fragment covers")


def test_third_allele_is_not_counted_as_reference():
    lines = pair_reads("r", 100, (0, 0, 0, 100), 120, 150)
    seq = list("A" * 100)
    seq[20] = "G"  # an allele that is neither ref (A) nor alt (T) at 120
    lines.append(sam("odd", 100, "".join(seq)))
    res = analyze_group(collect_fragments(lines, [120, 150]), G)
    assert res["per_site"][120]["O"] == 1
    assert res["pairs"][0].n_both == 100  # the 'O' fragment is excluded from the pair table


def test_mates_across_a_gap_give_linkage():
    lines = []
    for i in range(40):
        lines += [read(f"m{i}", 250, 100, [300]), read(f"m{i}", 650, 100, [700], flag=147)]
    for i in range(460):
        lines += [read(f"w{i}", 250, 100, []), read(f"w{i}", 650, 100, [], flag=147)]
    res = analyze_group(collect_fragments(lines, [300, 700]), [Site(300, "A", "T"), Site(700, "A", "T")])
    assert res["pairs"][0].verdict.startswith("LINKED") and res["pairs"][0].aa == 40


# ---------- choosing sites ----------

def test_auto_groups_on_the_real_variants_table_find_the_clusters_seen_in_the_dashboard():
    v = read_ivar_variants(REAL.read_text())
    known = best_snv_per_position(v)
    groups = auto_groups(site_table(v, min_depth=100), known, min_freq=0.02, window=10)
    sets = [tuple(s.pos for s in g) for g in groups]
    assert (5743, 5744) in sets
    assert (26042, 26044, 26046) in sets


def test_site_spec_parsing():
    v = read_ivar_variants(REAL.read_text())
    known = best_snv_per_position(v)
    g = parse_site_spec("5743,5744", known)
    assert [s.pos for s in g] == [5743, 5744] and g[0].alt == known[5743].alt
    assert parse_site_spec("5743:A,5744", known)[0].alt == "A"
    with pytest.raises(ValueError):
        parse_site_spec("5743", known)  # needs two sites
    with pytest.raises(ValueError):
        parse_site_spec("5743,1", known)  # position 1 is not in the table


# ---------- command line, end to end ----------

def test_cli_end_to_end(tmp_path, capsys):
    sam_path = tmp_path / "t.sam"
    sam_path.write_text("@SQ\tSN:ref\tLN:1000\n" + "\n".join(pair_reads("r", 100, (30, 0, 0, 970), 120, 150)) + "\n")
    tsv = tmp_path / "v.tsv"
    row = lambda pos: f"ref\t{pos}\tA\tT\t900\t0\t0\t100\t0\t0\t0.03\t1000\t0\tTRUE\n"
    tsv.write_text(HEADER + row(120) + row(150))
    out = tmp_path / "res"
    rc = linkage_check.main(["--alignment", str(sam_path), "--variants", str(tsv), "--sites", "120,150", "--out", str(out)])
    text = capsys.readouterr().out
    assert rc == 0 and "LINKED" in text
    assert "AA  30" in text and "RR  970" in text
    pairs = (tmp_path / "res.pairs.tsv").read_text().splitlines()
    assert pairs[1].split("\t")[:7] == ["1", "120", "150", "30", "0", "0", "970"]
    assert (tmp_path / "res.patterns.tsv").exists()

    rc = linkage_check.main(["--alignment", str(sam_path), "--variants", str(tsv), "--auto", "--out", str(out)])
    assert rc == 0 and "Group 1" in capsys.readouterr().out  # 0.03 minor freq, 30 bp apart, depth 1000


def test_cli_with_nothing_to_test_explains_itself(tmp_path, capsys):
    sam_path = tmp_path / "t.sam"
    sam_path.write_text("")
    tsv = tmp_path / "v.tsv"
    tsv.write_text(HEADER + f"ref\t120\tA\tT\t900\t0\t0\t100\t0\t0\t0.03\t1000\t0\tTRUE\n")
    assert linkage_check.main(["--alignment", str(sam_path), "--variants", str(tsv), "--auto"]) == 2
    assert "Nothing to test" in capsys.readouterr().err


def test_bam_without_samtools_gives_a_clear_error(monkeypatch):
    monkeypatch.setenv("PATH", "")
    with pytest.raises(RuntimeError, match="samtools"):
        list(alignment_lines("reads.bam"))


# ---------- reads behind the fragments (artifact check) ----------

def _mated(name, start, alts, flag1=99, flag2=147, length=100, gap=0):
    """A proper pair whose mates overlap the same sites: R1 at start, R2 at start+gap."""
    s2 = start + gap
    tlen = (s2 + length) - start
    return [read(name, start, length, [a for a in alts if start <= a < start + length], flag=flag1, pnext=s2, tlen=tlen),
            read(name, s2, length, [a for a in alts if s2 <= a < s2 + length], flag=flag2, pnext=start, tlen=-tlen)]


def test_meta_records_fragment_extent_and_reads():
    meta = {}
    collect_fragments(_mated("p", 100, [150], gap=20), [150], meta=meta)
    assert meta["p"]["extent"] == (100, 120)
    assert sorted(r[0] for r in meta["p"]["reads"]) == [99, 147]


def test_meta_falls_back_to_read_extent_without_mate_information():
    meta = {}
    collect_fragments([read("x", 100, 100, [150])], [150], meta=meta)
    assert meta["x"]["extent"] == (100, 100)


def test_alt_fragments_seen_by_both_mates_are_distinguished_from_single_mate_reference_fragments():
    """The pattern that would make a read-level frequency twice the fragment-level one."""
    lines = []
    for i in range(7):  # alt fragments: both mates cover the sites
        lines += _mated(f"alt{i}", 100, [120, 150], gap=10)
    for i in range(200):  # reference fragments: only one mate reaches the sites
        lines += [read(f"ref{i}", 100, 100, [], flag=99, pnext=300, tlen=300)]
    meta = {}
    frags = collect_fragments(lines, [120, 150], meta=meta)
    fs = fragment_summary(frags, meta, G)
    assert fs["alt"]["fragments"] == 7 and fs["alt"]["both_mates_overlap"] == 7
    assert fs["alt"]["reads"] == 14 and fs["alt"]["reverse_reads"] == 7
    assert fs["ref"]["fragments"] == 200 and fs["ref"]["both_mates_overlap"] == 0
    assert fs["alt"]["distinct_extents"] == 1


def test_duplicate_flagged_reads_are_counted():
    lines = _mated("a", 100, [120, 150], flag1=99 | 1024, flag2=147 | 1024)
    meta = {}
    frags = collect_fragments(lines, [120, 150], meta=meta)
    assert fragment_summary(frags, meta, G)["alt"]["duplicate_flagged"] == 2


def test_cli_writes_the_reads_behind_alt_fragments(tmp_path, capsys):
    lines = []
    for i in range(5):
        lines += _mated(f"alt{i}", 100, [120, 150], gap=10)
    for i in range(95):
        lines += _mated(f"ref{i}", 100, [], gap=10)
    sam_path = tmp_path / "t.sam"
    sam_path.write_text("\n".join(lines) + "\n")
    tsv = tmp_path / "v.tsv"
    row = lambda pos: f"ref\t{pos}\tA\tT\t900\t0\t0\t100\t0\t0\t0.05\t1000\t0\tTRUE\n"
    tsv.write_text(HEADER + row(120) + row(150))
    out = tmp_path / "res"
    assert linkage_check.main(["--alignment", str(sam_path), "--variants", str(tsv), "--sites", "120,150", "--out", str(out)]) == 0
    text = capsys.readouterr().out
    assert "Alt-carrying fragments: n=5" in text and "both mates over the sites: 5/5" in text
    rows = (tmp_path / "res.alt_fragments.tsv").read_text().splitlines()
    assert len(rows) == 6 and rows[1].split("\t")[:3] == ["1", "alt0", "AA"]
    assert rows[1].split("\t")[5] == "2" and "99:100:60:100M" in rows[1]


def test_fragment_start_positions_distinguish_amplicons():
    lines = []
    for i in range(4):  # alt fragments from one start
        lines += _mated(f"a{i}", 100, [120, 150], gap=10)
    for i in range(3):  # alt fragments from another start (a second amplicon)
        lines += _mated(f"b{i}", 110, [120, 150], gap=10)
    for i in range(50):
        lines += _mated(f"r{i}", 100, [], gap=10)
    meta = {}
    frags = collect_fragments(lines, [120, 150], meta=meta)
    fs = fragment_summary(frags, meta, G)
    assert dict(fs["alt"]["starts"]) == {100: 4, 110: 3}
    assert dict(fs["ref"]["starts"]) == {100: 50}
    from quasigraph.linkage import format_report
    text = format_report([analyze_group(frags, G, meta)], len(frags))
    assert "alt-carrying:   100 x4, 110 x3" in text and "reference-only: 100 x50" in text


def test_fragments_where_mates_disagree_are_reported_not_hidden():
    lines = []
    # one fragment: mate 1 shows the alt allele at 120, mate 2 (same position) shows the reference
    lines += [read("d1", 100, 100, [120], flag=99, pnext=100, tlen=100), read("d1", 100, 100, [], flag=147, pnext=100, tlen=-100)]
    lines += _mated("ok", 100, [120, 150], gap=0)   # a normal alt-carrying fragment
    lines += _mated("ref", 100, [], gap=0)
    meta = {}
    frags = collect_fragments(lines, [120, 150], meta=meta)
    res = analyze_group(frags, G, meta)
    assert res["discordant"] == {120: {"fragments": 1, "alt_on_one_mate": 1}}
    assert res["per_site"][120]["A"] == 1  # the disagreeing fragment is not counted as alt
    q, pos, bases = res["discordant_rows"][0]
    assert (q, pos) == ("d1", 120) and dict(bases) == {99: "T", 147: "A"}
    from quasigraph.linkage import format_report, write_tables
    assert "120: 1 fragment(s), 1 with the alt allele on one mate only" in format_report([res], len(frags))


def test_discordant_table_is_written_only_when_there_are_disagreements(tmp_path):
    lines = _mated("a", 100, [120, 150], gap=0) + _mated("r", 100, [], gap=0)
    meta = {}
    frags = collect_fragments(lines, [120, 150], meta=meta)
    from quasigraph.linkage import write_tables
    paths = write_tables(str(tmp_path / "x"), [analyze_group(frags, G, meta)])
    assert not any(p.endswith("discordant.tsv") for p in paths)
    lines += [read("d", 100, 100, [120], flag=99, pnext=100, tlen=100), read("d", 100, 100, [], flag=147, pnext=100, tlen=-100)]
    meta = {}
    frags = collect_fragments(lines, [120, 150], meta=meta)
    paths = write_tables(str(tmp_path / "y"), [analyze_group(frags, G, meta)])
    text = Path([p for p in paths if p.endswith("discordant.tsv")][0]).read_text()
    assert "d\t120\t99:T;147:A" in text
