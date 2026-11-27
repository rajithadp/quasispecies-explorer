"""Site-level quasispecies analysis from an iVar variants table.

When linked haplotypes are not available (short reads, or CliqueSNV found none),
the standard description of within-host diversity is per-site: how many positions
are polymorphic and at what frequency. This module works from ``ivar variants``
output (columns POS, REF, ALT, ALT_FREQ, TOTAL_DP, PASS, ...).

Definitions, per genome position:
  allele frequencies  reference allele = 1 - sum(ALT_FREQ), plus each ALT allele
  minor_freq          1 - (frequency of the most common allele)
  entropy_bits        Shannon entropy of the allele frequencies
  heterozygosity      1 - sum(p^2) = probability two random reads differ here
A position whose most common allele is not the reference (ALT_FREQ >= 0.5) is a
"fixed difference" from the reference (e.g. lineage-defining mutations), not diversity.

Limits: only PASS variants above iVar's own threshold are present, so very rare
variants and sequencing noise below that threshold are invisible, and variants
near the threshold are partly noise. Linkage between sites cannot be seen.
"""
from __future__ import annotations

import io

import numpy as np
import pandas as pd

REQUIRED = ["POS", "REF", "ALT", "ALT_FREQ", "TOTAL_DP", "PASS"]


def read_ivar_variants(text: str) -> pd.DataFrame:
    """PASS rows of an iVar variants table with tidy column names."""
    try:
        raw = pd.read_csv(io.StringIO(text), sep="\t", dtype=str)
    except Exception as e:  # pandas raises several types for malformed input
        raise ValueError(f"could not read the table as tab-separated text: {e}") from e
    missing = [c for c in REQUIRED if c not in raw.columns]
    if missing:
        raise ValueError(f"not an iVar variants table; missing column(s): {', '.join(missing)}")
    raw = raw[raw["PASS"].str.strip().str.upper() == "TRUE"]
    df = pd.DataFrame(
        {
            "position": raw["POS"].astype(int),
            "ref": raw["REF"],
            "alt": raw["ALT"],
            "freq": raw["ALT_FREQ"].astype(float),
            "depth": raw["TOTAL_DP"].astype(int),
        }
    ).reset_index(drop=True)
    df["is_indel"] = df["alt"].str.startswith(("+", "-"))
    return df


def site_table(variants: pd.DataFrame, min_depth: int = 0, snv_only: bool = True) -> pd.DataFrame:
    """One row per variant position with allele-frequency statistics (no frequency filter)."""
    v = variants[variants["depth"] >= min_depth]
    if snv_only:
        v = v[~v["is_indel"]]
    rows = []
    for pos, g in v.groupby("position", sort=True):
        alts = g["freq"].to_numpy(dtype=float)
        ref_freq = max(0.0, 1.0 - float(alts.sum()))
        p = np.append(alts, ref_freq)
        p = p / p.sum()
        nz = p[p > 0]
        rows.append(
            {
                "position": int(pos),
                "ref": g["ref"].iloc[0],
                "alt": ",".join(g["alt"]),
                "alt_freq": float(alts.sum()),
                "minor_freq": float(1.0 - p.max()),
                "entropy_bits": float(-(nz * np.log2(nz)).sum()),
                "heterozygosity": float(1.0 - (p**2).sum()),
                "depth": int(g["depth"].iloc[0]),
                "fixed_difference": bool(alts.sum() >= 0.5),
            }
        )
    return pd.DataFrame(
        rows,
        columns=["position", "ref", "alt", "alt_freq", "minor_freq", "entropy_bits",
                 "heterozygosity", "depth", "fixed_difference"],
    )


def frequency_spectrum(sites: pd.DataFrame) -> pd.DataFrame:
    """Counts of polymorphic positions by minor-allele frequency bin (fixed differences excluded)."""
    edges = [(0.0, 0.01), (0.01, 0.02), (0.02, 0.05), (0.05, 0.2), (0.2, 0.5001)]
    labels = ["<1%", "1-2%", "2-5%", "5-20%", "20-50%"]
    poly = sites[~sites["fixed_difference"]]
    counts = [int(((poly["minor_freq"] >= lo) & (poly["minor_freq"] < hi)).sum()) for lo, hi in edges]
    return pd.DataFrame({"minor_allele_frequency": labels, "sites": counts})


def site_summary(sites: pd.DataFrame, genome_length: int, min_freq: float = 0.02) -> dict:
    """Summary at an explicit minor-allele-frequency threshold.

    nucleotide_diversity_estimate = sum of heterozygosity over polymorphic sites
    / genome_length; positions without a called variant contribute zero.
    """
    if len(sites) and genome_length < int(sites["position"].max()):
        raise ValueError(f"genome length {genome_length} is shorter than the last variant position {int(sites['position'].max())}")
    poly = sites[(~sites["fixed_difference"]) & (sites["minor_freq"] >= min_freq)]
    return {
        "min_minor_freq": min_freq,
        "variable_sites": int(len(poly)),
        "fixed_differences": int(sites["fixed_difference"].sum()),
        "mean_minor_freq": float(poly["minor_freq"].mean()) if len(poly) else 0.0,
        "nucleotide_diversity_estimate": float(poly["heterozygosity"].sum() / genome_length),
        "median_depth": int(sites["depth"].median()) if len(sites) else 0,
        "genome_length": int(genome_length),
    }
