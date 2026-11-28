"""Draw a window of a variant-call graph: backbone stretches, allele bubbles and read-supported linked paths.

Columns are spaced evenly (not by genome distance), so a window of a few variable sites stays readable.
Reference alleles sit on the middle line, alternative alleles above and below. Marker area grows with the allele's
frequency, and every allele is labelled with its base and frequency. Minor G>T / C>A alleles (artifact-prone) are
ringed in orange and labelled 'ox'. Linked paths are drawn as coloured lines through the alleles they cover (teal when every site is
reference, a warm colour when a path carries alternative alleles), with the fragment count in the legend.
"""
from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from .variant_graph import _OXIDATIVE, VariantGraph  # noqa: E402

MAJOR, MINOR, OX = "#1f4e79", "#9ecae1", "#d98a2b"
REF_PATH = "#2a9d8f"                                   # every path that is reference at all its sites
ALT_PATHS = ["#c1440e", "#7b2d8e", "#b8860b", "#a31621"]   # paths carrying alternative alleles, in order of appearance


def columns_in_region(vg: VariantGraph, start: int, end: int) -> list:
    """Indices of the graph columns that overlap the 1-based window [start, end]."""
    if start > end:
        raise ValueError("region start must not be after its end")
    out = []
    for i, col in enumerate(vg.columns):
        n = vg.nodes[col[0]]
        if n.end >= start and n.start <= end:
            out.append(i)
    return out


def layout_region(vg: VariantGraph, start: int, end: int) -> dict:
    """{node id: (x, y)} for the columns in the window; x = column rank, y = 0 for the reference allele / backbone."""
    cols = columns_in_region(vg, start, end)
    pos = {}
    for x, ci in enumerate(cols):
        col = vg.columns[ci]
        if vg.nodes[col[0]].kind == "backbone":
            pos[col[0]] = (float(x), 0.0)
            continue
        ref_node = col[0]                                  # the reference allele is always first in a bubble
        others = sorted(col[1:], key=lambda n: -vg.nodes[n].freq)
        pos[ref_node] = (float(x), 0.0)
        for k, n in enumerate(others):
            step = k // 2 + 1
            pos[n] = (float(x), float(step if k % 2 == 0 else -step))
    return pos


def _is_oxidative(vg: VariantGraph, node) -> bool:
    return (node.kind == "allele" and node.role == "minor" and not node.matches_reference
            and (vg.reference[node.start - 1], node.seq) in _OXIDATIVE)


def plot_variant_graph(vg: VariantGraph, start: int, end: int, ax=None, title: str | None = None):
    pos = layout_region(vg, start, end)
    if not pos:
        raise ValueError(f"no graph columns overlap the region {start}-{end}")
    if ax is None:
        width = max(6.0, 1.3 * (max(x for x, _ in pos.values()) + 2))
        height = max(3.0, 0.9 * (max(y for _, y in pos.values()) - min(y for _, y in pos.values()) + 3.5))
        fig, ax = plt.subplots(figsize=(width, height))
    else:
        fig = ax.figure
    for u, v in vg.edges:
        if u in pos and v in pos:
            (x0, y0), (x1, y1) = pos[u], pos[v]
            ax.plot([x0, x1], [y0, y1], color="#bbbbbb", lw=1, zorder=1)
    for nid, (x, y) in pos.items():
        n = vg.nodes[nid]
        if n.kind == "backbone":
            ax.text(x, y, f"{n.end - n.start + 1} bp\n{n.start}-{n.end}", ha="center", va="center", fontsize=7,
                    bbox=dict(boxstyle="round,pad=0.3", fc="#f2f2f2", ec="#999999"), zorder=3)
            continue
        ox = _is_oxidative(vg, n)
        ax.scatter([x], [y], s=120 + 900 * n.freq ** 0.5, c=MAJOR if n.role == "major" else MINOR,
                   edgecolors=OX if ox else "#333333", linewidths=2.5 if ox else 0.8, zorder=3)
        text = f"{n.start}{n.seq}  {100 * n.freq:.1f}%" + ("  ox" if ox else "")
        if n.annotation:
            first = n.annotation.split(" | ")[0].split(":")
            text += "\n" + (f"{first[0]} {first[1]}" if len(first) >= 3 else first[0])
        ax.text(x, y - 0.38, text, ha="center", va="top", fontsize=7, zorder=4)
    n_alt = 0
    for k, e in enumerate(vg.evidence):
        pts = [pos[n] for n in e.nodes if n in pos]
        if "A" in e.pattern:
            color = ALT_PATHS[n_alt % len(ALT_PATHS)]
            n_alt += 1
        else:
            color = REF_PATH
        if len(pts) < 2:
            continue
        off = 0.06 * (k % 4 - 1.5)
        ax.plot([p[0] for p in pts], [p[1] + off for p in pts], color=color, lw=3,
                alpha=0.85, zorder=2, label=f"{e.name.replace('linked_', '')}: {e.fragments} fragments" + (f" ({e.joint.replace(',', ', ')})" if e.joint else ""))
    ax.set_xlim(-0.7, max(x for x, _ in pos.values()) + 0.7)
    ax.set_ylim(min(y for _, y in pos.values()) - 1.1, max(y for _, y in pos.values()) + 0.9)
    ax.axis("off")
    ax.set_title(title or f"{vg.region}: positions {start}-{end}", fontsize=10)
    handles, labels = ax.get_legend_handles_labels()
    if handles:
        ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.02), ncol=2, fontsize=7, frameon=False,
                  title="read-supported paths (A = alternative allele, R = reference allele)", title_fontsize=7)
    fig.tight_layout()
    return fig
