"""Build a variation graph from a reference and called variants, optionally with read-supported linkage paths.

Examples
  # graph of all variable sites (>= 2%, depth >= 100): reference + sample consensus paths
  python build_variant_graph.py --variants SRR12447392.variants.tsv --ref MN908947.3.fasta --out graph

  # add paths supported by the read pairs for groups of nearby sites you choose, or found automatically
  python build_variant_graph.py --variants X.variants.tsv --ref MN908947.3.fasta \\
      --alignment X_work/trimmed.sam --sites 26042,26044,26046 --sites 5743,5744 --out graph
  python build_variant_graph.py --variants X.variants.tsv --ref MN908947.3.fasta --alignment X_work/trimmed.sam --auto --out graph

Add --plot 26030-26060 to draw that window of the graph to <out>.png.
Writes <out>.gfa (GFA 1.0: segments with allele frequencies, links, and paths: reference, consensus, linked_*).
Frequencies on nodes are per-site allele frequencies. Paths other than reference, consensus and the read-supported
ones have NO frequency: without haplotypes the combination of alleles across distant sites is unknown.
"""
import argparse
import sys
from pathlib import Path

from quasigraph import read_ivar_variants, site_table
from quasigraph.linkage import (alignment_lines, analyze_group, auto_groups, best_snv_per_position, collect_fragments,
                                parse_site_spec)
from quasigraph.pipeline_utils import read_first_fasta_sequence
from quasigraph.variant_graph import build_variant_graph, to_gfa


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--variants", required=True)
    ap.add_argument("--ref", required=True)
    ap.add_argument("--alignment", help=".sam/.bam; needed for --sites / --auto (linkage evidence)")
    ap.add_argument("--sites", action="append", default=[], help="comma-separated positions of one linked group; repeatable")
    ap.add_argument("--auto", action="store_true", help="group nearby polymorphic sites automatically")
    ap.add_argument("--gff", help="GFF3 annotation: label alleles and linked paths with gene and amino-acid effect")
    ap.add_argument("--min-freq", type=float, default=0.02)
    ap.add_argument("--min-depth", type=int, default=100)
    ap.add_argument("--window", type=int, default=200, help="--auto: max distance between neighbouring sites of a group")
    ap.add_argument("--min-baseq", type=int, default=20)
    ap.add_argument("--plot", metavar="START-END", help="also draw this genome window (e.g. 26030-26060) to <out>.png")
    ap.add_argument("--out", default="graph")
    args = ap.parse_args(argv)
    try:
        fasta_text = Path(args.ref).read_text()
        ref = read_first_fasta_sequence(fasta_text)
        region = fasta_text.split("\n", 1)[0][1:].split()[0] if fasta_text.startswith(">") else "ref"
        variants = read_ivar_variants(Path(args.variants).read_text())
        known = best_snv_per_position(variants)
        groups = [parse_site_spec(s, known) for s in args.sites]
        if args.auto:
            groups += auto_groups(site_table(variants, min_depth=args.min_depth, snv_only=True), known, args.min_freq, args.window)
        evidence = []
        if groups:
            if not args.alignment:
                raise ValueError("--sites / --auto need --alignment (the reads that show which variants sit together)")
            frags = collect_fragments(alignment_lines(args.alignment), {s.pos for g in groups for s in g}, min_baseq=args.min_baseq)
            evidence = [(g, analyze_group(frags, g)["patterns"]) for g in groups]
        vg = build_variant_graph(ref, variants, args.min_freq, args.min_depth, region, evidence)
        if args.gff:
            from quasigraph.annotation import Annotator, parse_gff
            from quasigraph.graph_annotation import annotate_graph
            cds = parse_gff(Path(args.gff).read_text())
            if not cds:
                raise ValueError(f"no CDS features found in {args.gff}")
            annotate_graph(vg, Annotator(cds, ref))
    except (ValueError, OSError, RuntimeError) as e:
        print(f"error: {e}", file=sys.stderr)
        return 2
    window = None
    if args.plot:
        try:
            a, _, b = args.plot.partition("-")
            window = (int(a), int(b))
        except ValueError:
            print(f"error: --plot needs a window like 26030-26060, got {args.plot!r}", file=sys.stderr)
            return 2
    Path(f"{args.out}.gfa").write_text(to_gfa(vg))
    if window:
        from quasigraph.variant_graph_plot import plot_variant_graph
        try:
            plot_variant_graph(vg, *window).savefig(f"{args.out}.png", dpi=200)
        except ValueError as e:
            print(f"error: {e}", file=sys.stderr)
            return 2
        print(f"wrote {args.out}.png")
    st = vg.stats()
    print(f"Graph of {region}: {st['bubbles']} variable positions ({st['fixed_difference_bubbles']} fixed differences from the "
          f"reference, {st['minor_alleles']} minor alleles), {st['nodes']} nodes, {st['edges']} edges.")
    print("Skipped variants:", ", ".join(f"{k}={v}" for k, v in vg.skipped.items() if v) or "none")
    for e in vg.evidence:
        print(f"  path {e.name}: {e.fragments} read fragments" + (f"  -> {e.joint}" if e.joint else ""))
    for positions, why in vg.unplaced:
        print(f"  group {positions} NOT placed as a path: {why}")
    print(f"wrote {args.out}.gfa")
    print("Node frequencies are per-site allele frequencies; only the reference, consensus and linked_* paths are meaningful walks.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
