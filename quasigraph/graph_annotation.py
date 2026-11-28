"""Label the alleles and linked paths of a variant graph with their gene and amino-acid effect."""
from __future__ import annotations

from .annotation import Annotator, annotate_haplotype
from .variant_graph import VariantGraph


def _label(row: dict) -> str:
    if row["effect"] in ("noncoding", "indel"):
        return row["effect"]
    if row["aa_change"]:
        return f"{row['gene']}:{row['aa_change']}:{row['effect']}"
    return f"{row['gene']}:{row['effect']}"


def annotate_graph(vg: VariantGraph, annotator: Annotator) -> VariantGraph:
    """Fill VNode.annotation (each alternative allele alone) and EvidencePath.joint (the alternative alleles of a path
    together, with changes in one codon combined). Reference alleles are not annotated. A base read in two codon contexts
    (a frameshift join) gets both, separated by ' | '."""
    if annotator.ref != vg.reference:
        raise ValueError("the annotation reference differs from the graph's reference sequence")
    for n in vg.nodes.values():
        if n.kind != "allele" or n.matches_reference:
            continue
        n.annotation = " | ".join(_label(r) for r in annotator.annotate(n.start, vg.reference[n.start - 1], n.seq))
    for e in vg.evidence:
        subs = [(vg.nodes[i].start, vg.reference[vg.nodes[i].start - 1], vg.nodes[i].seq)
                for i in e.nodes if vg.nodes[i].kind == "allele" and not vg.nodes[i].matches_reference]
        if not subs:
            e.joint = ""
            continue
        rows = annotate_haplotype(annotator, subs)
        labels = [f"{r['gene']}:{r['aa_change']}" if r["aa_change"] else r["effect"] for r in rows]
        e.joint = ",".join(labels)
    return vg
