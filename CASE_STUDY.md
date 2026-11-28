# Case study: one SARS-CoV-2 sample, and what the tool could and could not say

Date: 2026-11-28. Sample: SRR12447392 (Illumina paired-end, SARS-CoV-2, Washington DC, collected 2020-03-24), reference MN908947.3. Seven more runs from the same public study were used for comparison (see below). Everything here was produced with the scripts in this repository.

## Question

How much intra-host diversity does this sample carry, and are there minor haplotypes (variants that sit together on the same genome)?

## What the pipeline did

minimap2 → samtools → `ivar trim` with the ARTIC v3 primers → `ivar variants` (down to 1%) → site-level analysis (minor-allele frequency at depth ≥ 100). `--skip-cliquesnv` leaves out haplotype reconstruction.

## Results

**1. CliqueSNV could not run on this machine.** It stopped with `java.lang.OutOfMemoryError: Java heap space` on a 5 GB computer, with `-tf 0.01` and also with `-tf 0.05 -t 20 -cm fast`. It exits normally anyway and writes a placeholder haplotype with frequency 1.0; the error appears only in its JSON file. The earlier version of this tool read that placeholder as a real one-haplotype result. The pipeline now reads the JSON and refuses it. No haplotypes were reconstructed from real reads, so the haplotype-level views (FASTA mode, variation graph) have not been used on real data.

**2. Site-level diversity is low.** 17 variable sites at ≥ 2% (π ≈ 2.9 × 10⁻⁵) and 9 fixed differences from the reference. Of the 17, 10 are G>T or C>A changes (see 4); without them: 7 sites, π ≈ 1.5 × 10⁻⁵.

**3. Two clusters of minor variants appear linked.** Using the read pairs themselves (`linkage_check.py`, after primer trimming):

| | Sites | Fragments with the alt alleles | Behaviour |
|---|---|---|---|
| Cluster 1 | 5743 G>C, 5744 T>C (orf1ab E1826D, Y1827H) | 7 of 524 (1.3%) | always together; both strands; mates agree |
| Cluster 2 | 26042 C>T, 26044 C>G, 26046 A>G (ORF3a T217I, Q218E as a pair) | 8 of 311 (2.6%) | always together; both strands; two fragment start positions |

These are consistent with minor haplotypes. They are not proven: each rests on 7–8 fragments.

A read-level check (`read_check.py`, on the alignment before primer trimming) found no sign of mis-mapped reads for the ORF3a cluster: its 16 alt reads carry as few mismatches elsewhere in the read as 713 reference reads over the same sites (0.38 vs 0.47 per read) and have the same mapping quality (60.0 vs 58.5). The orf1ab cluster is within the thresholds (0.50 vs 0.45 mismatches per read; mapping quality 55.6 vs 59.1), but 2 of its 7 alt fragments also carry two adjacent changes at 5765–5766, which may be a larger linked variant or an artifact; this was not followed up. The check cannot exclude PCR or library artifacts, and the alt reads have not been compared with a sequence database (BLAST).

Annotated as one haplotype (`annotate_haplotype.py`): ORF3a T217I with Q218E, and orf1ab E1826D with Y1827H.

**4. A sequencing artifact dominates the low-frequency variants across samples.** Eight Illumina samples from the same study:

| Sample | Variable sites ≥ 2% | G>T or C>A | Left without them |
|---|---|---|---|
| SRR12447411 | 153 | 99% | 1 |
| SRR12447389 | 27 | 85% | 4 |
| SRR12447392 | 17 | 59% | 7 |
| SRR12447388 | 12 | 91% | 1 |
| SRR12447406 | 10 | 67% | 3 |
| SRR12447508 | 10 | 50% | 5 |
| SRR12447404, SRR12447451 | 2 each | 0% | 2 |

G>T and C>A are the usual signature of oxidative damage during library preparation, although real mutations of these types exist, so the tool reports every number both ways and removes nothing. All 44 stop-gained minor variants in the eight samples are of this class. SRR12447411 looks ten times more diverse than the rest only because of it.

**5. Six samples share one lineage; two do not.** `fixed_differences.py` lists the fixed differences (alternative allele ≥ 50%) of the eight samples: 33 distinct changes. Five are carried by exactly six samples: 241 (noncoding), 3037 (orf1ab F924F, synonymous), 14408 (orf1ab P4715L), 23403 (S D614G) and 25563 (ORF3a Q57H). SRR12447389 and SRR12447411 carry none of them and differ from every other sample by 12–18 changes (from each other by 13). Within the six, the closest pair (SRR12447404, SRR12447508) differs by one change. Another 25 changes occur in one sample only and one (orf1ab T265I) in four. A sample is recorded as lacking a change when iVar lists no such variant, which can also mean the position was covered poorly. So the eight samples are not one lineage, and differences in diversity between them are not evidence of biology.

**6. A variation graph from variant calls works on real data.** For SRR12447392 (`build_variant_graph.py`, no haplotype reconstruction): 26 variable positions (9 fixed differences, 17 minor alleles), 78 nodes, 104 edges, plus the `reference` and `consensus` paths and four read-supported paths. Each alternative allele carries its gene and amino-acid change, and each linked path its joint change. In the ORF3a window the path with all three alternative alleles (8 fragments, against 361 for the reference allele) reads T217I with Q218E; the 26046 A>G change is synonymous on its own and changes Q218E only through its codon neighbour. The orf1ab path (5743, 5744; 7 alternative fragments against 517 reference) reads E1826D with Y1827H. Node frequencies are per-site allele frequencies; a complete path through many bubbles has no frequency, because the combination of distant alleles is not observed.

## What was not done

- No haplotypes from real reads. The graph is built from variant calls and read pairs, so it shows only the groups of nearby sites that single read pairs span.
- The graph was built for one sample only; the lineage comparison is the only step done across all eight.
- The two clusters have not been checked against a sequence database, and the extra changes at 5765–5766 on two orf1ab alt fragments have not been investigated.
- No error model; the artifact reading rests on the substitution spectrum alone. Variants near iVar's 1% calling threshold are partly noise.
- One study, one platform, three weeks of human samples from two places: the set tests the tool and cannot show a biological pattern.
- Per-gene counts are counts of variants, not rates per available site, and not dN/dS.
- iVar frequencies count reads; the linkage check counts fragments. They differ slightly.

## What this suggests for the tool

An explicit error model for low-frequency calls (the artifact reading rests on the substitution spectrum alone); graphs for every sample and a way to compare them; and the same analysis on samples whose metadata can answer a biological question (a time series from one patient, or samples from different tissues or treatments), because this study's eight samples, two places and a few weeks cannot.
