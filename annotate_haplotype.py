"""Annotate linked variants as ONE haplotype: what do the changes do together to each codon?

Example (the sites come from the linkage check; alleles are taken from the iVar table, or given as POS:ALT)
  python annotate_haplotype.py --variants SRR12447392.variants.tsv --gff MN908947.3.gff3 --ref MN908947.3.fasta \\
      --sites 26042,26044,26046 --sites 5743,5744

One table per group: the codons hit, the joint amino-acid change, what each change would do alone, and whether the
combination is emergent (gives an amino acid none of the single changes gives). Only valid if the variants really are
on the same genome: this script does not check that; use linkage_check.py first.
"""
import argparse
import sys
from pathlib import Path

from quasigraph import read_ivar_variants
from quasigraph.annotation import Annotator, annotate_haplotype, parse_gff
from quasigraph.linkage import best_snv_per_position, parse_site_spec
from quasigraph.pipeline_utils import read_first_fasta_sequence


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--variants", required=True, help="iVar variants table of the sample")
    ap.add_argument("--gff", required=True)
    ap.add_argument("--ref", required=True)
    ap.add_argument("--sites", action="append", required=True, help="comma-separated positions of one linked group; repeatable")
    args = ap.parse_args(argv)
    try:
        cds = parse_gff(Path(args.gff).read_text())
        ann = Annotator(cds, read_first_fasta_sequence(Path(args.ref).read_text()))
        known = best_snv_per_position(read_ivar_variants(Path(args.variants).read_text()))
        for spec in args.sites:
            group = parse_site_spec(spec, known)
            rows = annotate_haplotype(ann, [(s.pos, s.ref, s.alt) for s in group])
            print(f"\n=== haplotype: {', '.join(f'{s.pos}{s.ref}>{s.alt}' for s in group)} ===")
            for r in rows:
                where = f"{r['gene']} codon {r['codon']}" if r["gene"] else "noncoding"
                joint = f"{r['ref_codon']}>{r['alt_codon']} {r['aa_change']} ({r['effect']})" if r["aa_change"] else r["effect"]
                print(f"  positions {r['positions']:<12} {where:<18} {joint}")
                if r["positions"].count(",") or r["emergent"]:
                    print(f"      alone: {r['single_effects']}" + ("   <- EMERGENT: the combination gives an amino acid no single change gives" if r["emergent"] else ""))
    except (ValueError, OSError) as e:
        print(f"error: {e}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
