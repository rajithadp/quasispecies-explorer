"""Check whether the reads behind a linked minor haplotype look like ordinary reads of the sample.

Example
  python read_check.py --alignment SRR12447392_work/trimmed.sam --variants SRR12447392.variants.tsv \\
      --ref MN908947.3.fasta --sites 26042,26044,26046 --sites 5743,5744 --out readcheck

For each group of sites: alt reads (alternative allele at every covered site) are compared with reference reads
over the same sites on mismatches elsewhere in the read, soft-clipping, mapping quality, indels and strand.
Writes <out>.alt_reads.tsv (one row per alt read) and <out>.alt_reads.fasta (the alt reads, to check with BLAST).
Fixed differences from the reference (variant frequency >= 50%) are not counted as mismatches.
"""
import argparse
import sys
from pathlib import Path

from quasigraph import read_ivar_variants
from quasigraph.linkage import alignment_lines, best_snv_per_position, parse_site_spec
from quasigraph.pipeline_utils import read_first_fasta_sequence
from quasigraph.readcheck import alt_fasta, compare_group, format_comparison, profile_group


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--alignment", required=True, help=".sam (or .bam via samtools)")
    ap.add_argument("--variants", required=True)
    ap.add_argument("--ref", required=True)
    ap.add_argument("--sites", action="append", required=True, help="comma-separated positions of one group; repeatable")
    ap.add_argument("--min-baseq", type=int, default=20)
    ap.add_argument("--min-mapq", type=int, default=0)
    ap.add_argument("--out", default="readcheck")
    args = ap.parse_args(argv)
    try:
        ref = read_first_fasta_sequence(Path(args.ref).read_text())
        variants = read_ivar_variants(Path(args.variants).read_text())
        known = best_snv_per_position(variants)
        fixed = set(variants.loc[(~variants["is_indel"]) & (variants["freq"] >= 0.5), "position"].astype(int))
        tsv = ["group\tread\tmate\tflag\tpos\tmapq\tcigar\tunexplained_mismatches\tmismatch_positions\tsoft_clipped\tstrand"]
        fasta = []
        for spec in args.sites:
            group = parse_site_spec(spec, known)
            label = "+".join(str(s.pos) for s in group)
            profiles = profile_group(alignment_lines(args.alignment), ref, group, fixed, args.min_baseq, args.min_mapq)
            print("\n" + format_comparison(f"{label} ({', '.join(f'{s.pos}{s.ref}>{s.alt}' for s in group)})",
                                           compare_group(profiles)))
            for p in profiles:
                if p.kind == "alt":
                    tsv.append("\t".join(map(str, [label, p.name, p.mate, p.flag, p.pos, p.mapq, p.cigar,
                                                   p.unexplained_mismatches, ",".join(map(str, p.mismatch_positions)) or "-",
                                                   p.soft_clipped, "rev" if p.reverse else "fwd"])))
            fasta.append(alt_fasta(profiles, label))
    except (ValueError, OSError, RuntimeError) as e:
        print(f"error: {e}", file=sys.stderr)
        return 2
    Path(f"{args.out}.alt_reads.tsv").write_text("\n".join(tsv) + "\n")
    Path(f"{args.out}.alt_reads.fasta").write_text("".join(fasta))
    print(f"\nwrote {args.out}.alt_reads.tsv and {args.out}.alt_reads.fasta")
    print("The FASTA can be pasted into NCBI BLAST (blastn, standard nucleotide collection) to see what the alt reads match.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
