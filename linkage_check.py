"""Check whether nearby minor variants sit on the same read fragments (physical linkage).

Examples
  # let the script group nearby polymorphic sites from the variants table
  python linkage_check.py --alignment real_haplotypes_work/trimmed.sam --variants real_haplotypes.variants.tsv --auto

  # test specific sites (alleles are taken from the variants table; override with POS:ALT)
  python linkage_check.py --alignment trimmed.sam --variants variants.tsv --sites 5743,5744 --sites 26042,26044,26046

The alignment can be .sam, or .bam through samtools. Results are printed and written to <out>.pairs.tsv,
<out>.patterns.tsv and <out>.alt_fragments.tsv (the reads behind every alt-carrying fragment).
"""
import argparse
import sys

from quasigraph import read_ivar_variants, site_table
from quasigraph.linkage import (
    alignment_lines,
    analyze_group,
    auto_groups,
    best_snv_per_position,
    collect_fragments,
    format_report,
    parse_site_spec,
    write_tables,
)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--alignment", required=True, help=".sam (or .bam via samtools) of the mapped, trimmed reads")
    ap.add_argument("--variants", required=True, help="iVar variants table (.tsv)")
    ap.add_argument("--sites", action="append", default=[], help="comma-separated positions, e.g. 5743,5744 or 5743:C,5744:C; repeatable")
    ap.add_argument("--auto", action="store_true", help="group nearby polymorphic sites automatically")
    ap.add_argument("--min-freq", type=float, default=0.02, help="--auto: minimum minor-allele frequency (default 0.02)")
    ap.add_argument("--min-depth", type=int, default=100, help="--auto: minimum read depth (default 100)")
    ap.add_argument("--window", type=int, default=200, help="--auto: max distance in bp between neighbouring sites in a group (default 200)")
    ap.add_argument("--min-baseq", type=int, default=20, help="ignore bases below this quality (default 20)")
    ap.add_argument("--min-mapq", type=int, default=0, help="ignore reads below this mapping quality (default 0)")
    ap.add_argument("--out", default="linkage", help="output prefix (default: linkage)")
    args = ap.parse_args(argv)

    variants = read_ivar_variants(open(args.variants).read())
    known = best_snv_per_position(variants)

    groups = []
    try:
        groups += [parse_site_spec(s, known) for s in args.sites]
    except ValueError as e:
        print(f"error: {e}", file=sys.stderr)
        return 2
    if args.auto:
        table = site_table(variants, min_depth=args.min_depth, snv_only=True)
        groups += auto_groups(table, known, args.min_freq, args.window)
    if not groups:
        print("Nothing to test: give --sites or --auto (and check that at least two polymorphic sites are close together).",
              file=sys.stderr)
        return 2
    print(f"Testing {len(groups)} group(s) of nearby sites ...", flush=True)

    positions = {s.pos for g in groups for s in g}
    meta = {}
    frags = collect_fragments(alignment_lines(args.alignment), positions, min_baseq=args.min_baseq, min_mapq=args.min_mapq, meta=meta)
    results = [analyze_group(frags, g, meta) for g in groups]
    print(format_report(results, len(frags)))
    for path in write_tables(args.out, results):
        print(f"wrote {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
