import pytest

from quasigraph.annotation import CODON_TABLE, Annotator, Cds, annotate_haplotype, parse_gff

GFF = """##gff-version 3
chr\tx\tgene\t1\t12\t.\t+\t.\tID=g1
chr\tx\tCDS\t4\t12\t.\t+\t0\tID=c1;gene=A;product=protein%3B A
chr\tx\tCDS\t3\t8\t.\t-\t0\tID=c2;gene=B;product=minus
"""


def test_codon_table_spot_checks():
    assert CODON_TABLE["ATG"] == "M" and CODON_TABLE["TGG"] == "W"
    assert CODON_TABLE["TAA"] == CODON_TABLE["TAG"] == CODON_TABLE["TGA"] == "*"
    assert CODON_TABLE["GAT"] == "D" and CODON_TABLE["GGT"] == "G" and CODON_TABLE["CCT"] == "P" and CODON_TABLE["CTT"] == "L"
    assert len(CODON_TABLE) == 64 and sorted(set(CODON_TABLE.values())) == sorted(set("ACDEFGHIKLMNPQRSTVWY*"))


def test_parse_gff_reads_cds_only_and_decodes_text():
    c = parse_gff(GFF)
    assert [x.gene for x in c] == ["A", "B"] and c[0].product == "protein; A" and c[1].strand == "-"


def test_plus_strand_synonymous_missense_stop():
    # positions 4-12 = ATG GAT TGG ; a leading 3 bases of noncoding sequence
    ref = "CCC" + "ATGGATTGG"
    a = Annotator(parse_gff(GFF)[:1], ref)
    (r,) = a.annotate(7, "G", "T")          # GAT -> TAT : D->Y, codon 2 position 1
    assert (r["gene"], r["codon"], r["codon_position"], r["aa_change"], r["effect"]) == ("A", 2, 1, "D2Y", "missense")
    (r,) = a.annotate(9, "T", "C")          # GAT -> GAC : synonymous
    assert r["effect"] == "synonymous" and r["aa_change"] == "D2D"
    (r,) = a.annotate(12, "G", "A")         # TGG -> TGA : stop gained
    assert r["effect"] == "stop_gained" and r["aa_change"] == "W3*"
    (r,) = a.annotate(1, "C", "T")
    assert r["effect"] == "noncoding" and r["gene"] == ""


def test_stop_lost():
    a = Annotator([Cds("c", "A", "", "+", ((1, 6),))], "ATGTAA")
    (r,) = a.annotate(4, "T", "C")           # TAA -> CAA : stop lost
    assert r["effect"] == "missense" or r["effect"] == "stop_lost"
    assert r["effect"] == "stop_lost" and r["ref_aa"] == "*"


def test_minus_strand_uses_reverse_complement():
    # genomic 3-8 is CAT CAT on '+' ; minus strand reads reverse complement: ATG ATG
    ref = "NN" + "CATCAT" + "NN"
    ref = ref.replace("N", "A")
    a = Annotator(parse_gff(GFF)[1:], ref)
    # first codon (coding order) uses genomic 8,7,6 -> complement of T,A,C = A,T,G
    (r,) = a.annotate(8, "T", "C")           # coding base A -> G : ATG -> GTG (M->V), codon 1 position 1
    assert (r["ref_codon"], r["alt_codon"], r["aa_change"], r["codon_position"]) == ("ATG", "GTG", "M1V", 1)


def test_frameshift_join_reads_the_shared_base_twice():
    # two segments (1-6) and (6-11) like orf1ab: position 6 is read twice (as codon 2 pos 3 and codon 3 pos 1)
    c = Cds("c", "O", "", "+", ((1, 6), (6, 11)))
    ref = "ATGGCTTGGACG"[:11]
    a = Annotator([c], ref)
    rows = a.annotate(6, "T", "C")
    assert len(rows) == 2
    assert {(r["codon"], r["codon_position"]) for r in rows} == {(2, 3), (3, 1)}


def test_incomplete_trailing_codon_and_unknown_base():
    a = Annotator([Cds("c", "A", "", "+", ((1, 4),))], "ATGA")
    assert a.annotate(4, "A", "C")[0]["effect"] == "incomplete_codon"
    b = Annotator([Cds("c", "A", "", "+", ((1, 3),))], "ANG")
    assert b.annotate(3, "G", "T")[0]["effect"] == "unknown"


def test_indels_and_bad_input():
    a = Annotator([Cds("c", "A", "", "+", ((1, 3),))], "ATG")
    assert a.annotate(1, "A", "+CC")[0]["effect"] == "indel"
    assert a.annotate(1, "A", "-T")[0]["effect"] == "indel"
    with pytest.raises(ValueError, match="mismatch"):
        a.annotate(1, "C", "T")
    with pytest.raises(ValueError, match="outside"):
        a.annotate(9, "A", "T")
    with pytest.raises(ValueError, match="phase"):
        parse_gff("c\tx\tCDS\t1\t3\t.\t+\t1\tID=z")


def test_every_single_substitution_matches_an_independent_translation():
    # brute force: all 64 codons x 3 positions x 3 alts against a direct translation
    for codon, aa in CODON_TABLE.items():
        a = Annotator([Cds("c", "A", "", "+", ((1, 3),))], codon)
        for i in range(3):
            for alt in "ACGT":
                if alt == codon[i]:
                    continue
                alt_codon = codon[:i] + alt + codon[i + 1:]
                (r,) = a.annotate(i + 1, codon[i], alt)
                expect = "synonymous" if CODON_TABLE[alt_codon] == aa else (
                    "stop_gained" if CODON_TABLE[alt_codon] == "*" else "stop_lost" if aa == "*" else "missense")
                assert r["effect"] == expect and r["alt_codon"] == alt_codon


def test_codon_table_degeneracy_matches_the_genetic_code():
    from collections import Counter
    n = Counter(CODON_TABLE.values())
    assert n["*"] == 3 and n["M"] == 1 and n["W"] == 1
    assert n["L"] == n["S"] == n["R"] == 6
    assert n["A"] == n["G"] == n["P"] == n["T"] == n["V"] == 4
    assert n["I"] == 3 and all(n[x] == 2 for x in "CDEFHKNQY")


# ---------- the real SARS-CoV-2 annotation (NCBI MN908947.3) and well-known sites ----------

def _real():
    from pathlib import Path
    cds = parse_gff((Path(__file__).parent / "data" / "MN908947.3.gff3").read_text())
    return cds, Annotator(cds, "A" * 29903)


def test_real_gff_structure():
    cds, _ = _real()
    assert [c.gene for c in cds] == ["orf1ab", "S", "ORF3a", "E", "M", "ORF6", "ORF7a", "ORF8", "N", "ORF10"]
    orf1ab = cds[0]
    assert orf1ab.segments == ((266, 13468), (13468, 21555)) and orf1ab.length == 21291 and orf1ab.length % 3 == 0
    assert all(c.length % 3 == 0 for c in cds)


@pytest.mark.parametrize("pos, gene, codon, cpos", [
    (23403, "S", 614, 2),        # A23403G = S:D614G
    (14408, "orf1ab", 4715, 2),  # C14408T = P4715L (RdRp P323L), second segment after the frameshift
    (3037, "orf1ab", 924, 3),    # C3037T, synonymous F924F
    (25563, "ORF3a", 57, 3),     # G25563T = ORF3a Q57H
    (21563, "S", 1, 1),          # S start codon
    (28274, "N", 1, 1),
])
def test_real_gff_well_known_sites_map_to_the_published_codons(pos, gene, codon, cpos):
    _, a = _real()
    rows = a.annotate(pos, "A", "C")
    assert [(r["gene"], r["codon"], r["codon_position"]) for r in rows] == [(gene, codon, cpos)]


def test_real_gff_utr_intergenic_and_shared_join_base():
    _, a = _real()
    assert a.annotate(241, "A", "C")[0]["effect"] == "noncoding"        # 5' UTR
    assert a.annotate(26230, "A", "C")[0]["effect"] == "noncoding"      # between ORF3a and E
    assert len(a.annotate(13468, "A", "C")) == 2                        # base read twice at the frameshift


# ---------- tables and the command line ----------

HEADER = "REGION\tPOS\tREF\tALT\tREF_DP\tREF_RV\tREF_QUAL\tALT_DP\tALT_RV\tALT_QUAL\tALT_FREQ\tTOTAL_DP\tPVAL\tPASS\n"


def _tsv(*rows):
    return HEADER + "".join(f"r\t{p}\t{r}\t{a}\t0\t0\t0\t0\t0\t0\t{f}\t{d}\t0\tTRUE\n" for p, r, a, f, d in rows)


MINI_REF = "CCC" + "ATGGATTGG"   # gene A = positions 4-12 (plus strand), see GFF above


def test_annotate_table_and_gene_summary():
    from quasigraph import read_ivar_variants
    from quasigraph.annotation import annotate_table, gene_summary
    cds = parse_gff(GFF)[:1]
    ann = Annotator(cds, MINI_REF)
    v = read_ivar_variants(_tsv((7, "G", "T", 0.03, 500),      # missense, G>T
                                (9, "T", "C", 0.04, 500),      # synonymous, other class
                                (1, "C", "T", 0.05, 500),      # noncoding
                                (12, "G", "A", 0.9, 500),      # fixed difference (stop gained)
                                (10, "T", "A", 0.03, 50),      # too shallow
                                (5, "T", "+AA", 0.2, 500)))    # indel
    t = annotate_table(v, ann)
    assert sorted(t["position"]) == [1, 7, 9, 12]
    assert dict(zip(t["position"], t["kind"])) == {1: "minor", 7: "minor", 9: "minor", 12: "fixed_difference"}
    assert bool(t.set_index("position").loc[7, "GT_CA"]) and not t.set_index("position").loc[9, "GT_CA"]
    g = gene_summary(t, cds).set_index("gene")
    assert g.loc["A", "minor_variants"] == 2 and g.loc["A", "missense"] == 1 and g.loc["A", "synonymous"] == 1
    assert g.loc["A", "minor_variants_without_GT_CA"] == 1 and g.loc["A", "missense_without_GT_CA"] == 0
    assert g.loc["A", "stop_gained"] == 0          # the stop is a fixed difference, not counted
    assert g.loc["(noncoding)", "minor_variants"] == 1
    assert g.loc["A", "per_kb"] == pytest.approx(2 / 9 * 1000)


def test_cli_writes_tables_and_reports_reference_mismatch(tmp_path, capsys):
    import annotate_sites
    (tmp_path / "x.gff").write_text(GFF)
    (tmp_path / "ref.fa").write_text(">r\n" + MINI_REF + "\n")
    (tmp_path / "S1.variants.tsv").write_text(_tsv((7, "G", "T", 0.03, 500), (12, "G", "A", 0.9, 500)))
    out = str(tmp_path / "o")
    base = ["--gff", str(tmp_path / "x.gff"), "--ref", str(tmp_path / "ref.fa"), "--out", out]
    assert annotate_sites.main(["--variants", str(tmp_path / "S1.variants.tsv")] + base) == 0
    text = capsys.readouterr().out
    assert "D2Y" in (tmp_path / "o.variants.tsv").read_text() and "A:W3*" in text and (tmp_path / "o.genes.tsv").exists()
    (tmp_path / "S2.variants.tsv").write_text(_tsv((7, "C", "T", 0.03, 500)))     # wrong REF base
    assert annotate_sites.main(["--variants", str(tmp_path / "S2.variants.tsv")] + base) == 2
    assert "mismatch" in capsys.readouterr().err


# ---------- cases found by mutation testing ----------

def test_minus_strand_every_base_of_every_segment_is_used():
    # minus gene in two segments: genomic (7,12) then (1,6); coding order is 12..7 then 6..1
    c = Cds("m", "M", "", "-", ((1, 6), (7, 12)))
    assert c.coding_positions() == [12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1]
    # all 12 positions map to a codon; the last coding base (genomic 1) is codon 4 position 3
    a = Annotator([c], "ACGTACGTACGT")
    assert [(r["codon"], r["codon_position"]) for r in a.annotate(1, "A", "C")] == [(4, 3)]
    assert [(r["codon"], r["codon_position"]) for r in a.annotate(12, "T", "C")] == [(1, 1)]


def test_two_base_trailing_codon_is_incomplete():
    a = Annotator([Cds("c", "A", "", "+", ((1, 5),))], "ATGAA")
    assert a.annotate(5, "A", "C")[0]["effect"] == "incomplete_codon"
    assert a.annotate(3, "G", "T")[0]["effect"] == "missense"   # complete first codon ATG -> ATT (M->I)


def test_position_just_past_the_reference_is_rejected():
    a = Annotator([Cds("c", "A", "", "+", ((1, 3),))], "ATG")
    with pytest.raises(ValueError, match="outside"):
        a.annotate(4, "A", "C")
    with pytest.raises(ValueError, match="outside"):
        a.annotate(0, "A", "C")


def test_annotate_table_frequency_threshold_is_inclusive():
    from quasigraph import read_ivar_variants
    from quasigraph.annotation import annotate_table
    a = Annotator(parse_gff(GFF)[:1], MINI_REF)
    v = read_ivar_variants(_tsv((7, "G", "T", 0.02, 100), (9, "T", "C", 0.0199, 100)))
    assert list(annotate_table(v, a)["position"]) == [7]


# ---------- linked variants annotated as one haplotype ----------

def _plus(codons: str):
    return Annotator([Cds("c", "G", "", "+", ((1, len(codons)),))], codons)


def test_joint_effect_differs_from_single_effects_and_is_flagged_emergent():
    a = _plus("AAA")                                   # K
    (r,) = annotate_haplotype(a, [(1, "A", "G"), (2, "A", "T")])   # GAA=E alone, ATA=I alone, GTA=V together
    assert (r["alt_codon"], r["aa_change"], r["effect"]) == ("GTA", "K1V", "missense")
    assert r["single_effects"] == "K1E; K1I" and r["emergent"] and r["positions"] == "1,2"


def test_joint_effect_not_emergent_when_it_equals_one_single_effect():
    a = _plus("CAA")                                   # Q; pos1 C>G -> GAA (E); pos3 A>G -> CAG (Q, synonymous)
    (r,) = annotate_haplotype(a, [(1, "C", "G"), (3, "A", "G")])
    assert r["aa_change"] == "Q1E" and r["single_effects"] == "Q1E; Q1Q" and not r["emergent"]


def test_changes_in_different_codons_give_separate_rows_and_noncoding_rows():
    a = Annotator([Cds("c", "G", "", "+", ((4, 9),))], "CCC" + "GATGGT")      # D G in positions 4-9
    rows = annotate_haplotype(a, [(1, "C", "T"), (4, "G", "T"), (9, "T", "C")])
    assert [(r["positions"], r["effect"]) for r in rows] == [("1", "noncoding"), ("4", "missense"), ("9", "synonymous")]


def test_joint_minus_strand_codon():
    c = Cds("m", "M", "", "-", ((1, 3),))
    a = Annotator([c], "CAT")                          # reverse complement ATG (M)
    (r,) = annotate_haplotype(a, [(1, "C", "A"), (3, "T", "C")])
    # coding order is genomic 3,2,1 with complement: bases A,T,G. genomic 3 T>C -> coding A>G; genomic 1 C>A -> coding G>T
    assert (r["ref_codon"], r["alt_codon"]) == ("ATG", "GTT") and r["aa_change"] == "M1V"


def test_joint_input_checks():
    a = _plus("AAA")
    with pytest.raises(ValueError, match="same position"):
        annotate_haplotype(a, [(1, "A", "G"), (1, "A", "T")])
    with pytest.raises(ValueError, match="mismatch"):
        annotate_haplotype(a, [(1, "C", "G"), (2, "A", "T")])
    with pytest.raises(ValueError, match="single-base"):
        annotate_haplotype(a, [(1, "A", "+GG"), (2, "A", "T")])


def test_frameshift_shared_base_is_annotated_in_both_codons():
    a = Annotator([Cds("c", "O", "", "+", ((1, 6), (6, 11)))], "ATGGCTTGGAC")
    rows = annotate_haplotype(a, [(6, "T", "C")])
    assert sorted(r["codon"] for r in rows) == [2, 3]


def test_the_two_real_clusters_of_this_sample_reproduce_the_single_site_table():
    ref = ["A"] * 29903
    for pos, base in {**dict(zip(range(5741, 5744), "GAG")), **dict(zip(range(5744, 5747), "TAC")),
                      **dict(zip(range(26041, 26044), "ACT")), **dict(zip(range(26044, 26047), "CAA"))}.items():
        ref[pos - 1] = base
    from pathlib import Path
    ann = Annotator(parse_gff((Path(__file__).parent / "data" / "MN908947.3.gff3").read_text()), "".join(ref))
    orf3a = annotate_haplotype(ann, [(26042, "C", "T"), (26044, "C", "G"), (26046, "A", "G")])
    assert [(r["codon"], r["aa_change"]) for r in orf3a] == [(217, "T217I"), (218, "Q218E")]
    assert orf3a[1]["single_effects"] == "Q218E; Q218Q" and not orf3a[1]["emergent"]
    orf1ab = annotate_haplotype(ann, [(5743, "G", "C"), (5744, "T", "C")])
    assert [(r["codon"], r["aa_change"]) for r in orf1ab] == [(1826, "E1826D"), (1827, "Y1827H")]


def test_input_order_does_not_matter():
    a = Annotator([Cds("c", "G", "", "+", ((4, 9),))], "CCC" + "AAAGGT")
    fwd = annotate_haplotype(a, [(1, "C", "T"), (4, "A", "G"), (5, "A", "T"), (9, "T", "C")])
    rev = annotate_haplotype(a, [(9, "T", "C"), (5, "A", "T"), (4, "A", "G"), (1, "C", "T")])
    assert fwd == rev
    assert [r["positions"] for r in rev] == ["1", "4,5", "9"] and rev[1]["single_effects"] == "K1E; K1I"


def test_rows_are_ordered_by_genomic_position_even_on_the_minus_strand():
    a = Annotator([Cds("m", "M", "", "-", ((1, 6),))], "CATCAT")      # codon 1 = genomic 6..4, codon 2 = genomic 3..1
    rows = annotate_haplotype(a, [(6, "T", "C"), (1, "C", "G")])
    assert [(r["positions"], r["codon"]) for r in rows] == [("1", 2), ("6", 1)]
