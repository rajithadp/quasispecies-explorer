"""Pure-Python helpers for the FASTQ-to-haplotype pipeline (testable without the bioinformatics tools)."""
from __future__ import annotations

import csv
import io
import json
import random
from dataclasses import dataclass

from .io import parse_haplotypes


def convert_cliquesnv_fasta(text: str):
    """CliqueSNV FASTA (``>1_fr_0.58``) -> (names, seqs, freqs), most frequent first.

    Names become ``hap_1, hap_2, ...`` and the frequency from the header is kept.
    Raises ValueError if the file has no frequencies, because then it is not
    CliqueSNV output in the expected format.
    """
    hs = parse_haplotypes(text)
    if hs.freqs is None:
        raise ValueError("no '_fr_<frequency>' found in headers; expected CliqueSNV 'short' defline format")
    order = sorted(range(len(hs.seqs)), key=lambda i: -hs.freqs[i])
    seqs = [hs.seqs[i] for i in order]
    freqs = [float(hs.freqs[i]) for i in order]
    names = [f"hap_{i + 1}" for i in range(len(seqs))]
    return names, seqs, freqs


@dataclass
class Variant:
    pos: int  # 1-based
    ref: str
    alt: str
    freq: float

    @property
    def is_indel(self) -> bool:
        return self.alt.startswith("+") or self.alt.startswith("-") or len(self.ref) != 1 or len(self.alt) != 1


def parse_ivar_variants(tsv_text: str, min_freq: float = 0.02) -> list:
    """Variants from an ``ivar variants`` TSV: PASS == TRUE and ALT_FREQ >= min_freq."""
    out = []
    for row in csv.DictReader(io.StringIO(tsv_text), delimiter="\t"):
        if str(row["PASS"]).strip().upper() != "TRUE":
            continue
        freq = float(row["ALT_FREQ"])
        if freq < min_freq:
            continue
        out.append(Variant(int(row["POS"]), row["REF"], row["ALT"], freq))
    return out


def read_first_fasta_sequence(text: str) -> str:
    seq, started = [], False
    for line in text.splitlines():
        line = line.strip()
        if line.startswith(">"):
            if started:
                break
            started = True
        elif line and started:
            seq.append(line)
    if not seq:
        raise ValueError("no sequence found in reference FASTA")
    return "".join(seq).upper()


def build_synthetic_cloud(ref_seq: str, variants, n: int, seed: int = 0):
    """Generate ``n`` SYNTHETIC haplotypes by sampling each variant site independently.

    This ignores linkage between sites: haplotype counts and rare-haplotype load
    computed from this output depend on ``n`` and are NOT properties of the
    sample. Only the per-site allele frequencies are real. Indels are skipped.
    Returns (names, seqs, skipped_indels).
    """
    rng = random.Random(seed)
    snvs, skipped = {}, 0
    for v in variants:
        if v.is_indel:
            skipped += 1
            continue
        if not 1 <= v.pos <= len(ref_seq):
            raise ValueError(f"variant position {v.pos} is outside the reference (length {len(ref_seq)})")
        snvs.setdefault(v.pos - 1, []).append(v)

    seqs = []
    for _ in range(n):
        chars = list(ref_seq)
        for idx, alts in snvs.items():
            r = rng.random()
            cum = 0.0
            for v in sorted(alts, key=lambda x: -x.freq):
                cum += v.freq
                if r < cum:
                    chars[idx] = v.alt
                    break
        seqs.append("".join(chars))
    return [f"hap_{i + 1}" for i in range(n)], seqs, skipped


CLIQUESNV_ERRORS = {
    1: "too many cliques",
    2: "average clique degree is too high",
    3: "too many cliques of cliques",
    4: "timeout reached",
    5: "all haplotypes are below the frequency threshold (-tf)",
    6: "runtime error: an exception inside CliqueSNV (see the stack trace in its console output)",
}


def read_cliquesnv_json(text: str) -> dict:
    """Status of a CliqueSNV run from its JSON output.

    CliqueSNV exits normally even when it fails: it prints a stack trace, writes a
    placeholder haplotype (frequency 1.0) and records the failure only in this JSON
    (``error`` and ``errorCode``). Returns {"error", "error_code", "found", "failed"}.
    """
    j = json.loads(text)
    error = str(j.get("error", "none"))
    code = int(j.get("errorCode", 0) or 0)
    return {
        "error": error,
        "error_code": code,
        "found": int(j.get("foundHaplotypes", 0) or 0),
        "failed": error.lower() != "none" or code != 0,
    }


def describe_cliquesnv_error(info: dict) -> str:
    detail = CLIQUESNV_ERRORS.get(info["error_code"], "unknown error code")
    return f"CliqueSNV reported '{info['error']}' (code {info['error_code']}): {detail}"
