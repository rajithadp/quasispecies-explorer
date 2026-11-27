"""Simulate a toy quasispecies: one master sequence plus mutants.

Mutations are concentrated on a few "hotspot" sites so that the hotspot
analysis has something to find. This is test data, not a biological model.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

BASES = "ACGT"


@dataclass
class Quasispecies:
    names: list
    seqs: list
    freqs: np.ndarray
    hotspots: np.ndarray  # 0-based positions


def simulate_quasispecies(
    length: int = 300,
    n_haplotypes: int = 12,
    n_hotspots: int = 8,
    hotspot_fraction: float = 0.7,
    max_mutations: int = 6,
    seed: int = 0,
) -> Quasispecies:
    if not 0 < n_hotspots < length:
        raise ValueError("n_hotspots must be between 1 and length - 1")
    if n_haplotypes < 1:
        raise ValueError("n_haplotypes must be at least 1")
    rng = np.random.default_rng(seed)
    master = "".join(rng.choice(list(BASES), size=length))
    hotspots = np.sort(rng.choice(length, size=n_hotspots, replace=False))

    p = np.full(length, (1 - hotspot_fraction) / (length - n_hotspots))
    p[hotspots] = hotspot_fraction / n_hotspots
    p = p / p.sum()

    seqs = [master]
    seen = {master}
    attempts = 0
    while len(seqs) < n_haplotypes:
        attempts += 1
        if attempts > 10000:
            raise RuntimeError("could not generate enough distinct haplotypes; relax parameters")
        k = int(rng.integers(1, max_mutations + 1))
        sites = rng.choice(length, size=min(k, length), replace=False, p=p)
        chars = list(master)
        for s in sites:
            chars[s] = rng.choice([b for b in BASES if b != chars[s]])
        cand = "".join(chars)
        if cand not in seen:
            seen.add(cand)
            seqs.append(cand)

    freqs = np.sort(rng.dirichlet(np.linspace(3.0, 0.5, n_haplotypes)))[::-1]
    names = ["master"] + [f"variant_{i}" for i in range(1, n_haplotypes)]
    return Quasispecies(names=names, seqs=seqs, freqs=freqs, hotspots=hotspots)
