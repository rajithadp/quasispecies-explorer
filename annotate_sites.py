"""Label variants with the gene, codon and amino-acid change they cause.

Example
  python annotate_sites.py --variants samples/*.variants.tsv --gff MN908947.3.gff3 --ref MN908947.3.fasta --out annotated

For each variants table (one per sample) this prints a per-gene summary of the minor variants (default:
frequency >= 2%, depth >= 100) and writes <out>.variants.tsv (every variant with gene, codon, amino-acid
change and effect) and <out>.genes.tsv (per-gene counts). Fixed differences (frequency >= 50%) are annotated
too but not counted in the gene summary.
"""
import argparse
import sys
from pathlib import Path

import pandas as pd

from quasigraph import read_ivar_variants
from quasigraph.annotation import Annotator, annotate_table, gene_summary, parse_gff
from quasigraph.compare import run_id_from_path
from quasigraph.pipeline_utils import read_first_fasta_sequence


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--variants", nargs="+", required=True)
    ap.add_argument("--gff", required=True, help="GFF3 annotation (CDS features) of the reference")
    ap.add_argument("--ref", required=True, help="reference FASTA the reads were mapped to")
    ap.add_argument("--min-freq", type=float, default=0.02)
    ap.add_argument("--min-depth", type=int, default=100)
    ap.add_argument("--out", default="annotated")
    args = ap.parse_args(argv)

    try:
        cds = parse_gff(Path(args.gff).read_text())
        if not cds:
            raise ValueError(f"no CDS features found in {args.gff}")
        annotator = Annotator(cds, read_first_fasta_sequence(Path(args.ref).read_text()))
        tables, genes = [], []
        for path in args.variants:
            sample = run_id_from_path(path)
            t = annotate_table(read_ivar_variants(Path(path).read_text()), annotator, args.min_freq, args.min_depth)
            g = gene_summary(t, cds)
            t.insert(0, "sample", sample)
            g.insert(0, "sample", sample)
            tables.append(t)
            genes.append(g)
    except (ValueError, OSError) as e:
        print(f"error: {e}", file=sys.stderr)
        return 2

    all_t, all_g = pd.concat(tables), pd.concat(genes)
    for sample, g in all_g.groupby("sample", sort=False):
        t = all_t[all_t["sample"] == sample]
        print(f"\n=== {sample}: minor variants (frequency {args.min_freq}-0.5, depth >= {args.min_depth}) per gene ===")
        show = g[g["minor_variants"] > 0][["gene", "minor_variants", "missense", "synonymous", "stop_gained",
                                           "minor_variants_without_GT_CA", "missense_without_GT_CA"]]
        print(show.to_string(index=False) if len(show) else "none")
        fixed = t[t["kind"] == "fixed_difference"]
        names = [f"{r.gene}:{r.aa_change}" if r.aa_change else f"{r.position}{r.ref}>{r.alt}({r.effect})"
                 for r in fixed.itertuples()]
        print("fixed differences:", ", ".join(names) if names else "none")
    all_t.to_csv(f"{args.out}.variants.tsv", sep="\t", index=False)
    all_g.to_csv(f"{args.out}.genes.tsv", sep="\t", index=False)
    print(f"\nwrote {args.out}.variants.tsv and {args.out}.genes.tsv")
    print("Counts are of variants, not rates per available site; they are not dN/dS. Codon numbers for orf1ab "
          "count along the whole polyprotein, not nsp numbering.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
