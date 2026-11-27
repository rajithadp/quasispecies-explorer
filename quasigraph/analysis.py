"""One call that runs every analysis step; used by the CLI and the Streamlit app."""
from __future__ import annotations

import pandas as pd

from .graph import build_variation_graph
from .metrics import site_entropy, summarize, top_variable_sites


def analyze(
    seqs,
    freqs=None,
    top_k: int = 10,
    rhl_threshold: float = 0.05,
    min_minor_freq: float = 0.0,
) -> dict:
    """Summary metrics, per-site entropy, top variable sites and the variation graph.

    ``min_minor_freq`` (default 0 = off) and ``rhl_threshold`` are explicit
    parameters, never hidden constants.
    """
    summary = summarize(seqs, freqs, rhl_threshold=rhl_threshold, min_minor_freq=min_minor_freq)
    ent = site_entropy(seqs, freqs, min_minor_freq=min_minor_freq)
    return {
        "summary": summary,
        "site_entropy": pd.DataFrame({"position": range(1, len(ent) + 1), "entropy_bits": ent}),
        "top_sites": top_variable_sites(seqs, freqs, k=top_k, min_minor_freq=min_minor_freq),
        "graph": build_variation_graph(seqs, freqs),
    }
