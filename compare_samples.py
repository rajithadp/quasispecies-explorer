"""Compare diversity across samples from their iVar variants tables, with an artifact check.

Example
  python compare_samples.py --variants samples/*.variants.tsv --metadata study_runs.tsv --out compare

For every sample the table gives the number of variable sites and the nucleotide-diversity estimate
twice: with all minor variants, and without G>T / C>A changes (a recognised sequencing/library artifact
class; see quasigraph/compare.py). Metadata columns (any tab-separated table with a run_accession column)
are joined on the sample name. Writes <out>.tsv and <out>.png.
"""
import argparse
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

from quasigraph.compare import compare_samples


def make_figure(table: pd.DataFrame):
    t = table.sort_values("variable_sites_all", ascending=True)
    other = t["variable_sites_without_GT_CA"]
    ox = t["variable_sites_all"] - other
    fig, ax = plt.subplots(figsize=(7, 0.45 * len(t) + 1.5))
    ax.barh(t["sample"], other, label="other changes", color="#3b6ea5")
    ax.barh(t["sample"], ox, left=other, label="G>T / C>A (artifact-prone)", color="#d98a2b")
    ax.set_xlabel("variable sites (minor-allele frequency >= threshold)")
    ax.legend(loc="lower right", frameon=False)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    return fig


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--variants", nargs="+", required=True, help="iVar variants tables, one per sample")
    ap.add_argument("--metadata", help="tab-separated table with a run_accession column (optional)")
    ap.add_argument("--id-column", default="run_accession")
    ap.add_argument("--genome-length", type=int, default=29903)
    ap.add_argument("--min-freq", type=float, default=0.02)
    ap.add_argument("--min-depth", type=int, default=100)
    ap.add_argument("--out", default="compare")
    args = ap.parse_args(argv)

    meta = pd.read_csv(args.metadata, sep="\t", dtype=str) if args.metadata else None
    try:
        table = compare_samples(args.variants, args.genome_length, args.min_freq, args.min_depth, meta, args.id_column)
    except ValueError as e:
        print(f"error: {e}", file=sys.stderr)
        return 2
    keep = [c for c in ["sample", "collection_date", "country", "variable_sites_all", "variable_sites_without_GT_CA",
                        "oxidative_fraction", "pi_all", "pi_without_GT_CA", "fixed_differences"] if c in table.columns]
    show = table[keep].copy()
    for c in ("pi_all", "pi_without_GT_CA"):
        show[c] = show[c].map(lambda x: f"{x:.2e}")
    show["oxidative_fraction"] = show["oxidative_fraction"].map(lambda x: "n/a" if x != x else f"{x:.0%}")
    print(show.to_string(index=False))
    table.to_csv(f"{args.out}.tsv", sep="\t", index=False)
    make_figure(table).savefig(f"{args.out}.png", dpi=150)
    print(f"\nwrote {args.out}.tsv and {args.out}.png")
    print(f"Minor variants are SNVs with frequency {args.min_freq}-0.5 at depth >= {args.min_depth}. "
          "'oxidative_fraction' is the share that are G>T or C>A; a high share means diversity is largely artifact-prone.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
