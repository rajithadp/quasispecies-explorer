"""Compare diversity across samples, with and without an artifact-prone substitution class.

G>T and C>A (the same change read on the other strand) are the signature of oxidative damage
(8-oxoguanine) introduced during library preparation, and low-frequency calls are often dominated by
them. They are also genuine SARS-CoV-2 mutations, so they are reported separately rather than removed:
every summary is given with all minor variants and without this class (a sensitivity analysis).

Definitions
  minor variant   SNV with frequency in [min_freq, 0.5) at depth >= min_depth
  fixed difference position where the alternative allele reaches 0.5 or more: never dropped
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

from .sites import read_ivar_variants, site_summary, site_table

OXIDATIVE_CLASSES = {("G", "T"), ("C", "A")}


def _is_class(variants: pd.DataFrame) -> pd.Series:
    return pd.Series([(r, a) in OXIDATIVE_CLASSES for r, a in zip(variants["ref"], variants["alt"])],
                     index=variants.index, dtype=bool)


def minor_variants(variants: pd.DataFrame, min_freq: float = 0.02, min_depth: int = 100) -> pd.DataFrame:
    """SNV rows with min_freq <= frequency < 0.5 and depth >= min_depth (both bounds inclusive/exclusive as stated)."""
    v = variants
    return v[(~v["is_indel"]) & (v["freq"] >= min_freq) & (v["freq"] < 0.5) & (v["depth"] >= min_depth)]


def substitution_counts(variants: pd.DataFrame, min_freq: float = 0.02, min_depth: int = 100) -> dict:
    """Counts of 'REF>ALT' among minor variants."""
    m = minor_variants(variants, min_freq, min_depth)
    out: dict = {}
    for r, a in zip(m["ref"], m["alt"]):
        out[f"{r}>{a}"] = out.get(f"{r}>{a}", 0) + 1
    return out


def drop_classes(variants: pd.DataFrame) -> pd.DataFrame:
    """Remove G>T / C>A rows that are minor (frequency < 0.5), whatever their depth.

    Fixed differences (frequency >= 0.5), other substitution classes and indels are kept.
    """
    drop = _is_class(variants) & (~variants["is_indel"]) & (variants["freq"] < 0.5)
    return variants[~drop]


def summarize_sample(variants: pd.DataFrame, genome_length: int = 29903, min_freq: float = 0.02,
                     min_depth: int = 100) -> dict:
    """Per-sample numbers with all minor variants and without G>T / C>A."""
    m = minor_variants(variants, min_freq, min_depth)
    n_ox = int(_is_class(m).sum())
    a = site_summary(site_table(variants, min_depth=min_depth), genome_length, min_freq)
    b = site_summary(site_table(drop_classes(variants), min_depth=min_depth), genome_length, min_freq)
    return {
        "minor_variants": int(len(m)),
        "oxidative_class": n_ox,
        "oxidative_fraction": n_ox / len(m) if len(m) else float("nan"),
        "variable_sites_all": a["variable_sites"],
        "variable_sites_without_GT_CA": b["variable_sites"],
        "pi_all": a["nucleotide_diversity_estimate"],
        "pi_without_GT_CA": b["nucleotide_diversity_estimate"],
        "fixed_differences": a["fixed_differences"],
        "median_depth": a["median_depth"],
    }


def run_id_from_path(path) -> str:
    """'samples/SRR1.variants.tsv' -> 'SRR1'."""
    return Path(path).name.split(".")[0]


def compare_samples(paths, genome_length: int = 29903, min_freq: float = 0.02, min_depth: int = 100,
                    metadata: pd.DataFrame | None = None, id_column: str = "run_accession") -> pd.DataFrame:
    """One row per variants table; metadata columns are joined on the sample name (left join, order kept)."""
    rows = []
    for p in paths:
        s = summarize_sample(read_ivar_variants(Path(p).read_text()), genome_length, min_freq, min_depth)
        rows.append({"sample": run_id_from_path(p), **s})
    table = pd.DataFrame(rows)
    if metadata is not None:
        if id_column not in metadata.columns:
            raise ValueError(f"metadata has no '{id_column}' column (columns: {', '.join(metadata.columns)})")
        meta = metadata.drop_duplicates(subset=id_column).rename(columns={id_column: "sample"})
        table = table.merge(meta, on="sample", how="left")
    return table
