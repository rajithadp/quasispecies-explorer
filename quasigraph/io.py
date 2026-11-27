"""Read and write aligned haplotypes as FASTA, keeping their frequencies.

Frequencies are read from the header in either of two forms:

* ``freq=0.25``        (written by this package)
* ``_fr_0.25``         (CliqueSNV's default ``>{id}_fr_{frequency}`` defline)

A header containing ``synthetic=true`` marks haplotypes that were generated
from allele frequencies instead of observed in reads. The dashboard uses this
to refuse haplotype-level statistics for such data.

If no header has a frequency, every record is weighted equally and the result
says so (``weights_inferred``). If only some headers have one, that is an error.
"""
from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np

_FREQ = re.compile(r"(?:freq=|_fr_)([0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?)")
_SYNTH = re.compile(r"synthetic=(?:true|1|yes)", re.IGNORECASE)


@dataclass
class HaplotypeSet:
    names: list
    seqs: list
    freqs: object  # np.ndarray, or None when the file carries no frequencies
    synthetic: bool
    weights_inferred: bool


def parse_haplotypes(text: str) -> HaplotypeSet:
    names, seqs, freqs, synth = [], [], [], False
    current = None
    for raw in text.splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            header = line[1:].strip()
            m = _FREQ.search(header)
            freqs.append(float(m.group(1)) if m else None)
            synth = synth or bool(_SYNTH.search(header))
            names.append(header.split()[0] if header else f"seq{len(names) + 1}")
            seqs.append([])
            current = seqs[-1]
        else:
            if current is None:
                raise ValueError("FASTA text must start with a '>' header line")
            current.append(line)
    seqs = ["".join(parts).upper() for parts in seqs]
    if not seqs:
        raise ValueError("no sequences found")
    n_with = sum(f is not None for f in freqs)
    if n_with == 0:
        return HaplotypeSet(names, seqs, None, synth, True)
    if n_with != len(freqs):
        raise ValueError("some headers have a frequency (freq=... or _fr_...) and others do not")
    return HaplotypeSet(names, seqs, np.array(freqs), synth, False)


def parse_fasta(text: str):
    """Return (names, seqs, freqs_or_None). Kept for backward compatibility."""
    hs = parse_haplotypes(text)
    return hs.names, hs.seqs, hs.freqs


def read_fasta(path):
    return parse_fasta(Path(path).read_text())


def to_fasta(names, seqs, freqs=None, synthetic: bool = False) -> str:
    lines = []
    for i, (name, seq) in enumerate(zip(names, seqs)):
        header = f">{name}"
        if freqs is not None:
            header += f" freq={float(freqs[i]):.6f}"
        if synthetic:
            header += " synthetic=true"
        lines.append(header)
        lines.append(seq)
    return "\n".join(lines) + "\n"


def write_fasta(path, names, seqs, freqs=None, synthetic: bool = False) -> None:
    Path(path).write_text(to_fasta(names, seqs, freqs, synthetic))
