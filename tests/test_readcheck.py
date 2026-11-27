import random

import pytest

import read_check
from quasigraph.linkage import Site
from quasigraph.readcheck import alt_fasta, compare_group, format_comparison, profile_group, profile_read

random.seed(5)
REF = "".join(random.choice("ACGT") for _ in range(600))
G = [Site(300, REF[299], "T" if REF[299] != "T" else "A"), Site(302, REF[301], "G" if REF[301] != "G" else "C")]


def sam(name, flag, pos, seq, cigar=None, mapq=60, qual=None):
    cigar = cigar or f"{len(seq)}M"
    qual = qual or "I" * len(seq)
    return "\t".join([name, str(flag), "ref", str(pos), str(mapq), cigar, "=", "1", "0", seq, qual])


def read_over_sites(name="r", flag=0, start=250, length=100, alt=False, extra=(), cigar=None, mapq=60, qual=None):
    seq = list(REF[start - 1:start - 1 + length])
    if alt:
        for s in G:
            seq[s.pos - start] = s.alt
    for p, b in extra:            # extra mismatches at reference position p
        seq[p - start] = b
    return sam(name, flag, start, "".join(seq), cigar, mapq, qual)


def other(b):
    return "A" if b != "A" else "C"


def prof(line, exclude=()):
    return profile_read(line.split("\t"), REF, G, set(exclude) | {s.pos for s in G})


def test_alt_ref_and_other_classification():
    assert prof(read_over_sites(alt=True)).kind == "alt"
    assert prof(read_over_sites(alt=False)).kind == "ref"
    mixed = read_over_sites()
    seq = list(mixed.split("\t")[9]); seq[300 - 250] = G[0].alt
    f = mixed.split("\t"); f[9] = "".join(seq)
    assert profile_read(f, REF, G, {300, 302}).kind == "other"
    # a third base at a group site is 'O', which is not alt and not ref
    f2 = read_over_sites().split("\t"); s2 = list(f2[9]); s2[0 + 50] = [b for b in "ACGT" if b not in (G[0].ref, G[0].alt)][0]
    f2[9] = "".join(s2)
    assert profile_read(f2, REF, G, {300, 302}).states[300] == "O"


def test_read_not_covering_a_site_is_skipped():
    assert prof(read_over_sites(start=100, length=50)) is None


def test_unexplained_mismatches_exclude_group_sites_and_fixed_differences():
    p = prof(read_over_sites(alt=True, extra=[(270, other(REF[269])), (280, other(REF[279]))]))
    assert p.unexplained_mismatches == 2 and p.mismatch_positions == [270, 280] and p.examined_bases == 98
    q = prof(read_over_sites(alt=True, extra=[(270, other(REF[269])), (280, other(REF[279]))]), exclude={270})
    assert q.unexplained_mismatches == 1 and q.mismatch_positions == [280]


def test_low_quality_bases_and_n_are_not_counted():
    qual = list("I" * 100); qual[20] = "#"      # position 270, Q2
    seq_line = read_over_sites(alt=True, extra=[(270, other(REF[269])), (280, "N")], qual="".join(qual))
    p = prof(seq_line)
    assert p.unexplained_mismatches == 0


def test_soft_clip_insertion_and_deletion_are_counted_and_alignment_is_walked_correctly():
    # 20S 40M 2I 20M 3D 18M : the alt allele sites must still be read at the right query positions
    start = 280
    aligned1 = list(REF[start - 1:start - 1 + 40])                     # ref 280..319
    ins = "TT"
    aligned2 = list(REF[319:339])                                     # ref 320..339
    aligned3 = list(REF[342:360])                                     # skip 340..342 (3D), ref 343..360
    aligned1[300 - start] = G[0].alt
    aligned1[302 - start] = G[1].alt
    seq = "G" * 20 + "".join(aligned1) + ins + "".join(aligned2) + "".join(aligned3)
    line = sam("c", 0, start, seq, cigar="20S40M2I20M3D18M")
    p = prof(line)
    assert p.kind == "alt" and p.soft_clipped == 20 and p.insertions == 1 and p.deletions == 1
    assert p.unexplained_mismatches == 0 and p.examined_bases == 76   # 40+20+18 aligned bases minus the 2 group sites


def test_group_split_stats_and_flags_for_mismapped_alt_reads():
    lines = ["@HD\tVN:1.6"]
    for i in range(30):   # clean reference reads
        lines.append(read_over_sites(f"ref{i}", start=240 + i % 5, length=100))
    for i in range(6):    # alt reads with 3 extra mismatches each: look like another sequence
        extra = [(p, other(REF[p - 1])) for p in (260 + i, 270 + i, 280 + i)]
        lines.append(read_over_sites(f"alt{i}", start=250, alt=True, extra=extra))
    lines.append(read_over_sites("secondary", flag=256, alt=True))              # excluded by flag
    lines.append(read_over_sites("unmapped", flag=4, alt=True))
    profiles = profile_group(lines, REF, G, set())
    res = compare_group(profiles)
    assert res["alt"]["n"] == 6 and res["ref"]["n"] == 30
    assert res["alt"]["mismatches_mean"] == pytest.approx(3.0) and res["ref"]["mismatches_mean"] == 0
    assert any("mismatch" in f for f in res["flags"])
    assert "FLAGS" in format_comparison("g", res)


def test_clean_alt_reads_raise_no_flag_and_other_flags_work():
    lines = [read_over_sites(f"ref{i}", start=240 + i) for i in range(10)]
    lines += [read_over_sites(f"alt{i}", start=250, alt=True) for i in range(4)]
    res = compare_group(profile_group(lines, REF, G, set()))
    assert res["flags"] == [] and "nothing here suggests" in format_comparison("g", res)
    # heavy soft clipping of alt reads
    clipped = [read_over_sites(f"alt{i}", start=250, alt=True, cigar="40S60M") for i in range(4)]
    # clipping 40 leading bases moves the aligned part: place the sites inside the aligned block instead
    clipped = []
    for i in range(4):
        seq = list(REF[289:389]); seq[300 - 290] = G[0].alt; seq[302 - 290] = G[1].alt
        clipped.append(sam(f"alt{i}", 0, 290, "G" * 30 + "".join(seq), cigar="30S100M"))
    res2 = compare_group(profile_group(lines[:10] + clipped, REF, G, set()))
    assert any("soft-clipped" in f for f in res2["flags"])
    # low mapping quality
    lowq = [read_over_sites(f"alt{i}", start=250, alt=True, mapq=5) for i in range(4)]
    res3 = compare_group(profile_group(lines[:10] + lowq, REF, G, set()))
    assert any("mapping quality" in f for f in res3["flags"])


def test_min_mapq_filter_and_missing_classes():
    lines = [read_over_sites("a", alt=True, mapq=3), read_over_sites("b", alt=False, mapq=3)]
    assert profile_group(lines, REF, G, set(), min_mapq=10) == []
    only_alt = compare_group(profile_group([read_over_sites("a", alt=True)], REF, G, set()))
    assert "Not enough reads" in format_comparison("g", only_alt)


def test_fasta_keeps_the_sam_orientation_and_labels_strand_and_mate():
    fwd = prof(read_over_sites("f", flag=64, alt=True))
    rev = prof(read_over_sites("r", flag=16 | 128, alt=True))       # SAM already stores reverse reads in reference orientation
    ref_read = prof(read_over_sites("x", alt=False))
    recs = alt_fasta([fwd, rev, ref_read], "g").strip().split("\n")
    assert recs == [">g|f/1|pos250|fwd", fwd.seq, ">g|r/2|pos250|rev", rev.seq]
    assert fwd.seq[300 - 250] == G[0].alt and rev.seq[300 - 250] == G[0].alt       # alt allele in the same place in both


def test_cli_end_to_end(tmp_path, capsys):
    (tmp_path / "ref.fa").write_text(">ref\n" + REF + "\n")
    H = "REGION\tPOS\tREF\tALT\tREF_DP\tREF_RV\tREF_QUAL\tALT_DP\tALT_RV\tALT_QUAL\tALT_FREQ\tTOTAL_DP\tPVAL\tPASS\n"
    rows = "".join(f"ref\t{s.pos}\t{s.ref}\t{s.alt}\t0\t0\t0\t0\t0\t0\t0.03\t500\t0\tTRUE\n" for s in G)
    fixed_pos = 270
    rows += f"ref\t{fixed_pos}\t{REF[fixed_pos - 1]}\t{other(REF[fixed_pos - 1])}\t0\t0\t0\t0\t0\t0\t0.99\t500\t0\tTRUE\n"
    (tmp_path / "v.tsv").write_text(H + rows)
    lines = ["@HD\tVN:1.6"]
    for i in range(8):
        lines.append(read_over_sites(f"ref{i}", start=250, extra=[(270, other(REF[269]))]))
    for i in range(3):
        lines.append(read_over_sites(f"alt{i}", start=250, alt=True, extra=[(270, other(REF[269]))]))
    (tmp_path / "x.sam").write_text("\n".join(lines) + "\n")
    out = str(tmp_path / "o")
    rc = read_check.main(["--alignment", str(tmp_path / "x.sam"), "--variants", str(tmp_path / "v.tsv"), "--ref", str(tmp_path / "ref.fa"),
                          "--sites", "300,302", "--out", out])
    assert rc == 0
    text = capsys.readouterr().out
    assert "alt=3  ref=8" in text and "No difference" in text            # the fixed difference at 270 is not a mismatch
    assert (tmp_path / "o.alt_reads.fasta").read_text().count(">") == 3
    assert len((tmp_path / "o.alt_reads.tsv").read_text().strip().split("\n")) == 4
    assert read_check.main(["--alignment", str(tmp_path / "x.sam"), "--variants", str(tmp_path / "v.tsv"), "--ref", str(tmp_path / "ref.fa"),
                            "--sites", "301,302"]) == 2


# ---------- thresholds and counts, tested at their boundaries ----------

def _refs(n=10):
    return [read_over_sites(f"ref{i}", start=240 + i) for i in range(n)]


def _alt_with(k_mismatches, n=4, name="alt"):
    out = []
    for i in range(n):
        extra = [(260 + 10 * j + i, other(REF[260 + 10 * j + i - 1])) for j in range(k_mismatches)]
        out.append(read_over_sites(f"{name}{i}", start=250, alt=True, extra=extra))
    return out


def test_mismatch_flag_threshold_is_one_extra_mismatch_per_read():
    flagged = compare_group(profile_group(_refs() + _alt_with(1), REF, G, set()))
    assert any("mismatch" in f for f in flagged["flags"])                         # difference exactly 1.0
    three_of_four = _alt_with(1, 3) + _alt_with(0, 1, "clean")                    # mean 0.75
    assert compare_group(profile_group(_refs() + three_of_four, REF, G, set()))["flags"] == []


def test_share_with_two_or_more_mismatches():
    res = compare_group(profile_group(_refs() + _alt_with(2), REF, G, set()))
    assert res["alt"]["share_2plus_mismatches"] == 1.0
    res = compare_group(profile_group(_refs() + _alt_with(1), REF, G, set()))
    assert res["alt"]["share_2plus_mismatches"] == 0.0


def _clipped(n_clip, n=4, start=290):
    out = []
    for i in range(n):
        seq = list(REF[start - 1:start + 99]); seq[300 - start] = G[0].alt; seq[302 - start] = G[1].alt
        out.append(sam(f"c{n_clip}_{i}", 0, start, "G" * n_clip + "".join(seq), cigar=f"{n_clip}S100M"))
    return out


def test_soft_clip_counts_from_ten_bases_and_flag_from_a_quarter_of_reads():
    ten = compare_group(profile_group(_refs() + _clipped(10, 1) + _alt_with(0, 3), REF, G, set()))
    assert ten["alt"]["share_clipped_10plus"] == pytest.approx(0.25)
    assert any("soft-clipped" in f for f in ten["flags"])                        # exactly 0.25 difference
    nine = compare_group(profile_group(_refs() + _clipped(9, 1) + _alt_with(0, 3), REF, G, set()))
    assert nine["alt"]["share_clipped_10plus"] == 0 and not nine["flags"]
    one_of_five = compare_group(profile_group(_refs() + _clipped(10, 1) + _alt_with(0, 4), REF, G, set()))
    assert not any("soft-clipped" in f for f in one_of_five["flags"])             # 0.20 < 0.25


def test_mapping_quality_flag_threshold_is_ten_points():
    ref_reads = [read_over_sites(f"r{i}", start=240 + i, mapq=60) for i in range(10)]
    ten = [read_over_sites(f"a{i}", start=250, alt=True, mapq=50) for i in range(4)]
    nine = [read_over_sites(f"a{i}", start=250, alt=True, mapq=51) for i in range(4)]
    assert any("mapping quality" in f for f in compare_group(profile_group(ref_reads + ten, REF, G, set()))["flags"])
    assert compare_group(profile_group(ref_reads + nine, REF, G, set()))["flags"] == []


def test_other_reads_are_counted_separately():
    mixed = read_over_sites("m", start=250)
    f = mixed.split("\t"); s = list(f[9]); s[300 - 250] = G[0].alt; f[9] = "".join(s)
    res = compare_group(profile_group(_refs() + _alt_with(0, 2) + ["\t".join(f)], REF, G, set()))
    assert (res["alt"]["n"], res["ref"]["n"], res["n_other"]) == (2, 10, 1)


def test_note_when_clipping_cannot_discriminate():
    ref_reads = []
    for i in range(10):
        seq = list(REF[289:389])
        ref_reads.append(sam(f"r{i}", 0, 290, "G" * 30 + "".join(seq), cigar="30S100M"))
    alt = [sam(f"a{i}", 0, 290, "G" * 30 + REF[289:289 + 10] + G[0].alt + REF[300:301] + G[1].alt + REF[302:389], cigar="30S100M")
           for i in range(3)]
    res = compare_group(profile_group(ref_reads + alt, REF, G, set()))
    assert res["alt"]["n"] == 3 and res["ref"]["share_clipped_10plus"] == 1.0
    assert len(res["notes"]) == 1 and "untrimmed alignment" in format_comparison("g", res)
    # ordinary unclipped reads: no note
    assert compare_group(profile_group(_refs() + _alt_with(0, 3), REF, G, set()))["notes"] == []


def test_clipping_note_threshold_is_ninety_percent_of_reference_reads():
    def mix(n_clipped, n=10):
        reads = [sam(f"c{i}", 0, 290, "G" * 30 + REF[289:389], cigar="30S100M") for i in range(n_clipped)]
        reads += [read_over_sites(f"u{i}", start=250 + i) for i in range(n - n_clipped)]
        return reads + _alt_with(0, 3)
    assert len(compare_group(profile_group(mix(9), REF, G, set()))["notes"]) == 1     # 90%
    assert compare_group(profile_group(mix(8), REF, G, set()))["notes"] == []         # 80%
