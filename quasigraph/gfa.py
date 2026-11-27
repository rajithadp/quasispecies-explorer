"""Export a VariationGraph as GFA 1.0 so other graph tools can read it.

Segments are the graph nodes (``FR`` tag = summed frequency of the haplotypes
through them), links join consecutive nodes, and each haplotype is a ``P`` path.
Alignment gaps ('-') are kept inside segment sequences as ordinary characters.
"""
from __future__ import annotations

from .graph import END, START, VariationGraph


def to_gfa(vg: VariationGraph, names=None) -> str:
    G = vg.graph
    if names is None:
        names = [f"hap{i + 1}" for i in range(len(vg.haplotypes))]
    if len(names) != len(vg.haplotypes):
        raise ValueError("names must match the number of distinct haplotypes in the graph")

    inner = [n for n in G.nodes if n not in (START, END)]
    inner.sort(key=lambda n: (n[0], -G.nodes[n]["weight"], n[1]))
    ids = {n: f"s{i + 1}" for i, n in enumerate(inner)}

    lines = ["H\tVN:Z:1.0"]
    for n in inner:
        lines.append(f"S\t{ids[n]}\t{G.nodes[n]['allele']}\tFR:f:{G.nodes[n]['weight']:.6f}")
    for u, v, d in G.edges(data=True):
        if u in (START, END) or v in (START, END):
            continue
        lines.append(f"L\t{ids[u]}\t+\t{ids[v]}\t+\t0M\tFR:f:{d['weight']:.6f}")
    for name, path, freq in zip(names, vg.paths, vg.freqs):
        segs = ",".join(f"{ids[n]}+" for n in path)
        lines.append(f"P\t{name}\t{segs}\t*\tFR:f:{float(freq):.6f}")
    return "\n".join(lines) + "\n"
