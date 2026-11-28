"""Which fixed differences from the reference do the samples share (lineage-like grouping)?

Example
  python fixed_differences.py --variants samples/*.variants.tsv --ref MN908947.3.fasta --gff MN908947.3.gff3 --out fixed

Prints every fixed difference (alternative allele >= 50% at depth >= 100) with its gene and amino-acid change, the samples
that carry it, the groups of samples with identical sets, and the number of differing changes between each pair.
A sample "not carrying" a change may also simply not cover that position (see quasigraph/shared_fixed.py).
"""
import argparse
import sys
from pathlib import Path

from quasigraph.compare import run_id_from_path
from quasigraph.shared_fixed import difference_matrix, fixed_set, groups, shared_table
from quasigraph.pipeline_utils import read_first_fasta_sequence
from quasigraph.sites import read_ivar_variants


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--variants", nargs="+", required=True)
    ap.add_argument("--ref", help="reference FASTA (needed with --gff)")
    ap.add_argument("--gff", help="GFF3 to label each change with gene and amino-acid effect")
    ap.add_argument("--min-depth", type=int, default=100)
    ap.add_argument("--out", default="fixed")
    args = ap.parse_args(argv)
    try:
        ann = None
        if args.gff:
            if not args.ref:
                raise ValueError("--gff needs --ref")
            from quasigraph.annotation import Annotator, parse_gff
            cds = parse_gff(Path(args.gff).read_text())
            if not cds:
                raise ValueError(f"no CDS features found in {args.gff}")
            ann = Annotator(cds, read_first_fasta_sequence(Path(args.ref).read_text()))
        sets = {}
        for p in args.variants:
            name = run_id_from_path(p)
            if name in sets:
                raise ValueError(f"two variants tables give the sample name {name}")
            sets[name] = fixed_set(read_ivar_variants(Path(p).read_text()), args.min_depth)
        table = shared_table(sets, ann)
    except (ValueError, OSError) as e:
        print(f"error: {e}", file=sys.stderr)
        return 2
    n = len(sets)
    show = table.drop(columns=["carried_by", "not_in"])
    print(f"{len(table)} fixed differences in {n} samples (alternative allele >= 50%, depth >= {args.min_depth})\n")
    print(show.to_string(index=False))
    print("\nGroups of samples with identical sets:")
    for members, k in groups(sets):
        print(f"  {k} changes: {', '.join(members)}")
    print("\nChanges present in one sample but not the other:")
    print(difference_matrix(sets).to_string())
    table.to_csv(f"{args.out}.tsv", sep="\t", index=False)
    print(f"\nwrote {args.out}.tsv")
    print("'Not carrying' means iVar listed no such variant; it can also mean the position was poorly covered.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
