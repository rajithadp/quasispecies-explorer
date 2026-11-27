"""Diversity metrics for a viral quasispecies given as aligned haplotypes.

A quasispecies is described by a set of haplotypes (aligned sequences of equal
length) and their relative frequencies. Alignment gaps ('-') are treated as
ordinary characters, so an indel counts as a difference.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass

import numpy as np
import pandas as pd


def collapse_haplotypes(seqs, freqs=None):
    """Merge identical sequences and normalise frequencies.

    Returns (unique_seqs, freqs) sorted by decreasing frequency, freqs sum to 1.
    If ``freqs`` is None, every input sequence counts equally.
    """
    seqs = [s.upper() for s in seqs]
    if not seqs:
        raise ValueError("no sequences given")
    if len({len(s) for s in seqs}) != 1:
        raise ValueError("sequences must be aligned: all must have the same length")
    if len(seqs[0]) == 0:
        raise ValueError("sequences are empty")
    if freqs is None:
        f = np.full(len(seqs), 1.0 / len(seqs))
    else:
        f = np.asarray(freqs, dtype=float)
        if len(f) != len(seqs):
            raise ValueError("number of frequencies must match number of sequences")
        if (f < 0).any() or f.sum() <= 0:
            raise ValueError("frequencies must be non-negative with a positive sum")
        f = f / f.sum()
    merged: dict[str, float] = {}
    for s, x in zip(seqs, f):
        merged[s] = merged.get(s, 0.0) + float(x)
    items = sorted(merged.items(), key=lambda kv: (-kv[1], kv[0]))
    return [s for s, _ in items], np.array([x for _, x in items])


def encode(seqs) -> np.ndarray:
    """Haplotypes x sites matrix of character codes."""
    return np.array([np.frombuffer(s.encode("ascii"), dtype=np.uint8) for s in seqs])


def pairwise_distances(seqs) -> np.ndarray:
    """Hamming distance matrix between aligned sequences."""
    A = encode(seqs)
    D = np.zeros((len(A), len(A)), dtype=int)
    for i in range(len(A)):
        D[i] = (A != A[i]).sum(axis=1)
    return D


def site_entropy(seqs, freqs=None, min_minor_freq: float = 0.0) -> np.ndarray:
    """Shannon entropy (bits) of the allele distribution at each alignment column.

    ``min_minor_freq`` is an explicit, optional noise filter: a column whose
    combined non-major allele frequency is below it is reported as 0. The default
    (0.0) applies no filter.
    """
    haps, f = collapse_haplotypes(seqs, freqs)
    A = encode(haps)
    syms = np.unique(A)
    # allele frequency per column; normalise each column so that rounding error in the
    # haplotype frequencies cannot give an invariant column a tiny non-zero entropy
    P = np.array([(A == sym).T @ f for sym in syms])  # alleles x sites
    P = P / P.sum(axis=0)
    with np.errstate(divide="ignore", invalid="ignore"):
        terms = np.where(P > 0, P * np.log2(P), 0.0)
    H = -terms.sum(axis=0)
    H[H < 1e-12] = 0.0
    if min_minor_freq > 0:
        H[(1.0 - P.max(axis=0)) < min_minor_freq] = 0.0
    return H


def rare_haplotype_load(freqs, threshold: float = 0.05) -> dict:
    """Share of the population carried by rare haplotypes.

    Returns both readings of "rare haplotype load":
    - mass: total frequency of haplotypes below ``threshold``
    - count_fraction: fraction of distinct haplotypes below ``threshold``
    Frequencies are renormalised to sum to 1.
    """
    f = np.asarray(freqs, dtype=float)
    f = f / f.sum()
    rare = f < threshold
    return {
        "mass": float(f[rare].sum()),
        "count_fraction": float(rare.mean()),
        "n_rare": int(rare.sum()),
    }


def mutations_vs_consensus(seqs, freqs=None) -> list:
    """Differences between each input sequence and the frequency-weighted consensus."""
    seqs = [s.upper() for s in seqs]
    haps, f = collapse_haplotypes(seqs, freqs)
    A = encode(haps)
    consensus = np.empty(A.shape[1], dtype=np.uint8)
    for j in range(A.shape[1]):
        vals = np.unique(A[:, j])
        weights = [f[A[:, j] == v].sum() for v in vals]
        consensus[j] = vals[int(np.argmax(weights))]
    return [int((np.frombuffer(s.encode("ascii"), dtype=np.uint8) != consensus).sum()) for s in seqs]


@dataclass
class QuasispeciesSummary:
    n_input_sequences: int
    n_haplotypes: int
    alignment_length: int
    segregating_sites: int
    dominant_haplotype_freq: float
    haplotype_entropy_bits: float
    effective_haplotypes: float
    nucleotide_diversity: float  # expected pairwise differences per site
    mean_pairwise_distance: float  # expected pairwise differences per genome
    mean_site_entropy_bits: float
    rhl_mass: float  # total frequency of haplotypes below the RHL threshold
    rhl_count_fraction: float  # fraction of distinct haplotypes below the RHL threshold
    n_rare_haplotypes: int

    def as_dict(self) -> dict:
        return asdict(self)


def summarize(seqs, freqs=None, rhl_threshold: float = 0.05, min_minor_freq: float = 0.0) -> QuasispeciesSummary:
    """Size and diversity of a quasispecies.

    - haplotype_entropy_bits: Shannon entropy of the haplotype frequencies
    - effective_haplotypes: 2**entropy, the number of equally common haplotypes
      that would give the same entropy
    - nucleotide_diversity: frequency-weighted expected number of differences
      between two randomly drawn genomes, divided by alignment length
    """
    n_input = len(seqs)
    haps, f = collapse_haplotypes(seqs, freqs)
    A = encode(haps)
    L = A.shape[1]
    D = pairwise_distances(haps)
    expected_dist = float(f @ D @ f)
    p = f[f > 0]
    h_entropy = float(-(p * np.log2(p)).sum())
    ent = site_entropy(haps, f, min_minor_freq=min_minor_freq)
    rhl = rare_haplotype_load(f, rhl_threshold)
    return QuasispeciesSummary(
        n_input_sequences=n_input,
        n_haplotypes=len(haps),
        alignment_length=L,
        segregating_sites=int((A != A[0]).any(axis=0).sum()),
        dominant_haplotype_freq=float(f[0]),
        haplotype_entropy_bits=h_entropy,
        effective_haplotypes=float(2.0 ** h_entropy),
        nucleotide_diversity=expected_dist / L,
        mean_pairwise_distance=expected_dist,
        mean_site_entropy_bits=float(ent.mean()),
        rhl_mass=rhl["mass"],
        rhl_count_fraction=rhl["count_fraction"],
        n_rare_haplotypes=rhl["n_rare"],
    )


def top_variable_sites(seqs, freqs=None, k: int = 10, min_minor_freq: float = 0.0) -> pd.DataFrame:
    """The k most variable alignment columns (1-based positions)."""
    haps, f = collapse_haplotypes(seqs, freqs)
    A = encode(haps)
    ent = site_entropy(haps, f, min_minor_freq=min_minor_freq)
    rows = []
    for pos in np.argsort(-ent)[:k]:
        if ent[pos] == 0:
            break
        freq_by_allele: dict[str, float] = {}
        for ch, x in zip(A[:, pos], f):
            c = chr(ch)
            freq_by_allele[c] = freq_by_allele.get(c, 0.0) + float(x)
        major = max(freq_by_allele, key=freq_by_allele.get)
        rows.append(
            {
                "position": int(pos) + 1,
                "entropy_bits": float(ent[pos]),
                "n_alleles": len(freq_by_allele),
                "major_allele": major,
                "major_allele_freq": freq_by_allele[major],
            }
        )
    return pd.DataFrame(
        rows,
        columns=["position", "entropy_bits", "n_alleles", "major_allele", "major_allele_freq"],
    )
