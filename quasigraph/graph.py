"""A simple sequence variation graph for a set of aligned haplotypes.

The alignment is cut into blocks. A block is either a stretch of columns where
all haplotypes agree (one node) or a maximal run of variable columns (a
"bubble" with one node per distinct allele). Each haplotype is a path from
START to END. Node and edge weights are the summed frequencies of the
haplotypes passing through them.
"""
from __future__ import annotations

from dataclasses import dataclass

import networkx as nx
import numpy as np

from .metrics import collapse_haplotypes, encode

START = "START"
END = "END"


@dataclass
class VariationGraph:
    graph: nx.DiGraph
    haplotypes: list
    freqs: np.ndarray
    blocks: list  # (start, end, is_variable), 0-based half-open columns
    paths: list  # per haplotype: list of node ids

    @property
    def n_bubbles(self) -> int:
        return sum(1 for _, _, var in self.blocks if var)

    def reconstruct(self, i: int) -> str:
        """Rebuild haplotype i by walking its path through the graph."""
        return "".join(self.graph.nodes[n]["allele"] for n in self.paths[i])

    def stats(self) -> dict:
        return {
            "nodes": self.graph.number_of_nodes() - 2,
            "edges": self.graph.number_of_edges(),
            "blocks": len(self.blocks),
            "bubbles": self.n_bubbles,
            "haplotypes": len(self.haplotypes),
        }


def build_variation_graph(seqs, freqs=None) -> VariationGraph:
    haps, f = collapse_haplotypes(seqs, freqs)
    A = encode(haps)
    L = A.shape[1]
    variable = (A != A[0]).any(axis=0)

    blocks = []
    start = 0
    for i in range(1, L + 1):
        if i == L or variable[i] != variable[start]:
            blocks.append((start, i, bool(variable[start])))
            start = i

    G = nx.DiGraph()
    G.add_node(START, block=-1, weight=1.0, allele="", variable=False)
    G.add_node(END, block=len(blocks), weight=1.0, allele="", variable=False)

    paths = []
    for seq, x in zip(haps, f):
        prev = START
        path = []
        for b, (st, en, var) in enumerate(blocks):
            allele = seq[st:en]
            node = (b, allele)
            if node not in G:
                G.add_node(node, block=b, start=st, end=en, allele=allele, variable=var, weight=0.0)
            G.nodes[node]["weight"] += float(x)
            if G.has_edge(prev, node):
                G[prev][node]["weight"] += float(x)
            else:
                G.add_edge(prev, node, weight=float(x))
            prev = node
            path.append(node)
        if G.has_edge(prev, END):
            G[prev][END]["weight"] += float(x)
        else:
            G.add_edge(prev, END, weight=float(x))
        paths.append(path)

    return VariationGraph(graph=G, haplotypes=haps, freqs=f, blocks=blocks, paths=paths)
