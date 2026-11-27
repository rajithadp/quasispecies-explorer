"""Draw a VariationGraph with matplotlib (blocks left to right, major allele on the middle line)."""
from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from .graph import END, START, VariationGraph  # noqa: E402


def _offsets(k: int):
    out = [0.0]
    step = 1
    while len(out) < k:
        out.append(float(step))
        if len(out) < k:
            out.append(float(-step))
        step += 1
    return out


def layout(vg: VariationGraph) -> dict:
    G = vg.graph
    nb = len(vg.blocks)
    pos = {START: (0.0, 0.0), END: (nb + 1.0, 0.0)}
    for b in range(nb):
        nodes = [n for n in G.nodes if n not in (START, END) and n[0] == b]
        nodes.sort(key=lambda n: -G.nodes[n]["weight"])
        for n, y in zip(nodes, _offsets(len(nodes))):
            pos[n] = (b + 1.0, y)
    return pos


def draw_graph(vg: VariationGraph, ax=None, max_label: int = 6):
    G = vg.graph
    pos = layout(vg)
    nb = len(vg.blocks)
    max_k = max(len([n for n in G.nodes if n not in (START, END) and n[0] == b]) for b in range(nb))

    if ax is None:
        width = min(18.0, max(7.0, 0.9 * (nb + 2)))  # capped so large graphs stay viewable
        fig, ax = plt.subplots(figsize=(width, max(3.0, 0.8 * max_k + 2.0)))
    else:
        fig = ax.figure

    for u, v, d in G.edges(data=True):
        (x1, y1), (x2, y2) = pos[u], pos[v]
        ax.plot([x1, x2], [y1, y2], color="#4c78a8", alpha=0.55,
                linewidth=0.8 + 7.0 * d["weight"], solid_capstyle="round", zorder=1)

    for n, (x, y) in pos.items():
        d = G.nodes[n]
        if n in (START, END):
            ax.scatter([x], [y], s=60, color="#444444", zorder=3)
            ax.text(x, y - 0.35, n.lower(), ha="center", va="top", fontsize=7, color="#444444")
            continue
        color = "#e45756" if d["variable"] else "#9aa5b1"
        ax.scatter([x], [y], s=150 + 900 * d["weight"], color=color, edgecolor="white", zorder=3)
        if d["variable"]:
            label = d["allele"] if len(d["allele"]) <= max_label else d["allele"][: max_label - 1] + "…"
        else:
            label = f"{d['end'] - d['start']} bp"
        ax.text(x, y, label, ha="center", va="center", fontsize=6, zorder=4)

    for b, (st, en, var) in enumerate(vg.blocks):
        if var:
            txt = f"{st + 1}" if en - st == 1 else f"{st + 1}–{en}"
            ax.text(b + 1.0, max_k / 2 + 0.9, txt, ha="center", fontsize=6, color="#e45756", rotation=90)

    ax.set_xlim(-0.5, nb + 1.5)
    ax.set_ylim(-max_k / 2 - 1.2, max_k / 2 + 2.6)
    ax.axis("off")
    ax.set_title("Variation graph (red = variable positions, width = haplotype frequency)", fontsize=9)
    fig.tight_layout()
    return fig


def plot_site_positions(sites, genome_length: int, min_freq: float):
    """Minor-allele frequency of every polymorphic position along the genome (log scale)."""
    poly = sites[~sites["fixed_difference"]]
    fig, ax = plt.subplots(figsize=(11, 3.2))
    above = poly[poly["minor_freq"] >= min_freq]
    below = poly[poly["minor_freq"] < min_freq]
    ax.vlines(poly["position"], 1e-4, poly["minor_freq"], color="#c8d1db", linewidth=0.8, zorder=1)
    ax.scatter(below["position"], below["minor_freq"], s=14, color="#9aa5b1", label=f"below {min_freq:.1%}", zorder=2)
    ax.scatter(above["position"], above["minor_freq"], s=26, color="#e45756", label=f"at or above {min_freq:.1%}", zorder=3)
    ax.axhline(min_freq, color="#444444", linestyle="--", linewidth=0.8)
    ax.set_yscale("log")
    ax.set_ylim(5e-3, 0.6)
    ax.set_xlim(0, genome_length)
    ax.set_xlabel("Genome position (bp)")
    ax.set_ylabel("Minor-allele frequency")
    ax.legend(loc="upper right", fontsize=8, frameon=False)
    fig.tight_layout()
    return fig


def plot_site_spectrum(spectrum):
    """Bar chart of how many polymorphic positions fall in each minor-allele-frequency bin."""
    fig, ax = plt.subplots(figsize=(5.5, 3.2))
    ax.bar(spectrum["minor_allele_frequency"], spectrum["sites"], color="#4c78a8")
    for i, n in enumerate(spectrum["sites"]):
        ax.text(i, n, str(n), ha="center", va="bottom", fontsize=8)
    ax.set_xlabel("Minor-allele frequency")
    ax.set_ylabel("Number of sites")
    fig.tight_layout()
    return fig
