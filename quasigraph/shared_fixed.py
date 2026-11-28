"""Which fixed differences from the reference do samples share?

A fixed difference is a single-base change whose alternative allele reaches frequency 0.5 or more at depth >=
min_depth (the same definition as sites.py). Samples that carry the same set belong to the same lineage-like group;
samples that differ in a few of them are close but not identical. This is a descriptive comparison, not a tree.

Limit: a sample is recorded as NOT carrying a change when iVar lists no such variant there. That is the reference
allele, but it can also be a position the sample did not cover well (iVar lists only variants, not coverage), so a
single missing change is weaker evidence than a change that is present.
"""
from __future__ import annotations

from itertools import combinations

import pandas as pd

from .annotation import Annotator
from .graph_annotation import _label


def fixed_set(variants: pd.DataFrame, min_depth: int = 100) -> set:
    """{(position, ref, alt)} for SNVs at frequency >= 0.5 and depth >= min_depth."""
    v = variants[(~variants["is_indel"]) & (variants["freq"] >= 0.5) & (variants["depth"] >= min_depth)]
    return {(int(r.position), r.ref, r.alt) for r in v.itertuples()}


def shared_table(sets: dict, annotator: Annotator | None = None) -> pd.DataFrame:
    """One row per fixed difference seen in any sample: how many samples carry it, which, and which do not."""
    samples = list(sets)
    allchanges = sorted(set().union(*sets.values())) if sets else []
    rows = []
    for pos, ref, alt in allchanges:
        carried = [s for s in samples if (pos, ref, alt) in sets[s]]
        label = ""
        if annotator is not None:
            label = " | ".join(_label(r) for r in annotator.annotate(pos, ref, alt))
        rows.append({"position": pos, "ref": ref, "alt": alt, "annotation": label, "n_samples": len(carried),
                     "carried_by": ",".join(carried), "not_in": ",".join(s for s in samples if s not in carried)})
    t = pd.DataFrame(rows, columns=["position", "ref", "alt", "annotation", "n_samples", "carried_by", "not_in"])
    return t.sort_values(["n_samples", "position"], ascending=[False, True], kind="stable").reset_index(drop=True)


def groups(sets: dict) -> list:
    """Samples with exactly the same set of fixed differences: [(samples, n_changes)], largest group first."""
    by: dict = {}
    for s, fs in sets.items():
        by.setdefault(frozenset(fs), []).append(s)
    out = [(tuple(m), len(k)) for k, m in by.items()]
    return sorted(out, key=lambda g: (-len(g[0]), g[0]))


def difference_matrix(sets: dict) -> pd.DataFrame:
    """Number of fixed differences present in one sample but not the other (symmetric difference)."""
    names = list(sets)
    m = pd.DataFrame(0, index=names, columns=names, dtype=int)
    for a, b in combinations(names, 2):
        m.loc[a, b] = m.loc[b, a] = len(sets[a] ^ sets[b])
    return m
