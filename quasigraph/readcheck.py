"""Do the reads that carry a linked minor haplotype otherwise look like normal reads of the sample?

A cluster of linked minor variants can be real, or reads from another source (a different sequence, a contaminant,
a paralogous region) that map badly onto the reference and so show several 'variants' at once. Such reads usually
carry extra mismatches elsewhere in the read and more soft-clipping than ordinary reads over the same sites.

For each group of sites this module splits the reads covering them into
  alt    every covered group site shows the alternative allele (at least one site covered)
  ref    every covered group site shows the reference allele
  other  anything else (mixed, or bases below the quality cutoff only)
and compares alt with ref reads: mismatches to the reference OUTSIDE the group sites and outside known fixed
differences ("unexplained mismatches"), soft-clipped bases, mapping quality, insertions/deletions and strand.
It also writes the alt reads as FASTA (in the reference's orientation, as stored in the SAM file) so they can be checked with BLAST.

What this can and cannot show: it can reveal reads that do not fit the reference. It cannot reveal a PCR or
library artifact that produces normal-looking reads, and with a handful of alt reads the comparison is descriptive.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from statistics import mean, median

from .linkage import SAM_EXCLUDE, bases_at_sites, parse_cigar, reference_end


@dataclass
class ReadProfile:
    name: str
    flag: int
    pos: int
    mapq: int
    cigar: str
    length: int
    states: dict                      # group position -> 'A' (alt) / 'R' (ref) / 'O' (other base)
    unexplained_mismatches: int = 0
    examined_bases: int = 0
    soft_clipped: int = 0
    insertions: int = 0
    deletions: int = 0
    seq: str = ""
    qual: str = ""
    mismatch_positions: list = field(default_factory=list)

    @property
    def reverse(self) -> bool:
        return bool(self.flag & 16)

    @property
    def mate(self) -> int:
        return 1 if self.flag & 64 else 2 if self.flag & 128 else 0

    @property
    def kind(self) -> str:
        v = set(self.states.values())
        if v == {"A"}:
            return "alt"
        if v == {"R"}:
            return "ref"
        return "other"


def profile_read(fields: list, ref: str, group: list, exclude: set, min_baseq: int = 20) -> ReadProfile | None:
    """Profile one SAM record against group sites; None if it does not cover any group site."""
    pos, cigar = int(fields[3]), fields[5]
    ops = parse_cigar(cigar)
    seq, qual = fields[9], fields[10]
    sites = sorted(s.pos for s in group)
    obs = bases_at_sites(pos, ops, seq, qual, sites, min_baseq)
    if not obs:
        return None
    by_pos = {s.pos: s for s in group}
    states = {}
    for p, b in obs.items():
        states[p] = "A" if b == by_pos[p].alt else "R" if b == by_pos[p].ref else "O"
    rp = ReadProfile(fields[0], int(fields[1]), pos, int(fields[4]), cigar, len(seq), states, seq=seq, qual=qual)
    rpos, qpos = pos, 0
    for length, op in ops:
        if op in "M=X":
            for i in range(length):
                p, q = rpos + i, qpos + i
                if q >= len(seq) or p > len(ref) or p in exclude:
                    continue
                if qual != "*" and ord(qual[q]) - 33 < min_baseq:
                    continue
                b = seq[q].upper()
                if b == "N":
                    continue
                rp.examined_bases += 1
                if b != ref[p - 1]:
                    rp.unexplained_mismatches += 1
                    rp.mismatch_positions.append(p)
            rpos += length
            qpos += length
        elif op == "S":
            rp.soft_clipped += length
            qpos += length
        elif op == "I":
            rp.insertions += 1
            qpos += length
        elif op in "DN":
            rp.deletions += op == "D"
            rpos += length
    return rp


def profile_group(lines, ref: str, group: list, exclude: set, min_baseq: int = 20, min_mapq: int = 0) -> list:
    """Profiles of all usable reads that cover at least one site of the group."""
    exclude = set(exclude) | {s.pos for s in group}
    lo, hi = min(s.pos for s in group), max(s.pos for s in group)
    out = []
    for line in lines:
        if not line or line[0] == "@":
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 11 or int(f[1]) & SAM_EXCLUDE or f[5] == "*" or int(f[4]) < min_mapq:
            continue
        pos = int(f[3])
        if reference_end(pos, parse_cigar(f[5])) <= lo or pos > hi:
            continue
        rp = profile_read(f, ref, group, exclude, min_baseq)
        if rp is not None:
            out.append(rp)
    return out


def _stats(reads: list) -> dict:
    if not reads:
        return {"n": 0}
    mm = [r.unexplained_mismatches for r in reads]
    clip = [r.soft_clipped for r in reads]
    return {
        "n": len(reads),
        "mismatches_mean": mean(mm), "mismatches_median": median(mm),
        "share_2plus_mismatches": sum(m >= 2 for m in mm) / len(reads),
        "soft_clipped_mean": mean(clip), "share_clipped_10plus": sum(c >= 10 for c in clip) / len(reads),
        "mapq_mean": mean(r.mapq for r in reads),
        "indels_per_read": sum(r.insertions + r.deletions for r in reads) / len(reads),
        "reverse_share": sum(r.reverse for r in reads) / len(reads),
        "length_mean": mean(r.length for r in reads),
    }


def compare_group(profiles: list) -> dict:
    """Stats for alt and ref reads and plain-language flags (heuristic, descriptive)."""
    alt = [p for p in profiles if p.kind == "alt"]
    ref = [p for p in profiles if p.kind == "ref"]
    a, r = _stats(alt), _stats(ref)
    flags = []
    if a["n"] and r["n"]:
        if a["mismatches_mean"] - r["mismatches_mean"] >= 1.0:
            flags.append("alt reads carry at least one more unexplained mismatch per read than reference reads: "
                         "possible mis-mapped reads from another sequence")
        if a["share_clipped_10plus"] - r["share_clipped_10plus"] >= 0.25:
            flags.append("alt reads are soft-clipped (10+ bases) much more often than reference reads: "
                         "parts of them do not align to the reference")
        if r["mapq_mean"] - a["mapq_mean"] >= 10:
            flags.append("alt reads have clearly lower mapping quality than reference reads")
    notes = []
    if r["n"] and r["share_clipped_10plus"] >= 0.9:
        notes.append("almost every reference read has 10+ soft-clipped bases (typical of amplicon reads after primer trimming), "
                     "so the clipping comparison cannot discriminate here; repeat on the untrimmed alignment (mapped.sorted.bam)")
    return {"alt": a, "ref": r, "n_other": len(profiles) - a["n"] - r["n"], "flags": flags, "notes": notes}


def alt_fasta(profiles: list, label: str) -> str:
    """Alt reads as FASTA. SAM stores every read in the reference's orientation already (reverse-strand reads are
    stored reverse-complemented), so the sequences are written as stored; the header says which strand was sequenced."""
    lines = []
    for p in profiles:
        if p.kind != "alt":
            continue
        lines.append(f">{label}|{p.name}/{p.mate}|pos{p.pos}|{'rev' if p.reverse else 'fwd'}\n{p.seq}")
    return "\n".join(lines) + ("\n" if lines else "")


def format_comparison(label: str, result: dict) -> str:
    a, r = result["alt"], result["ref"]
    out = [f"=== {label} ===", f"reads: alt={a['n']}  ref={r['n']}  other={result['n_other']}"]
    if not a["n"] or not r["n"]:
        out.append("Not enough reads of both kinds to compare.")
        return "\n".join(out)
    rows = [("unexplained mismatches per read (mean)", "mismatches_mean", "{:.2f}"),
            ("reads with 2+ unexplained mismatches", "share_2plus_mismatches", "{:.0%}"),
            ("soft-clipped bases per read (mean)", "soft_clipped_mean", "{:.1f}"),
            ("reads with 10+ soft-clipped bases", "share_clipped_10plus", "{:.0%}"),
            ("mapping quality (mean)", "mapq_mean", "{:.1f}"),
            ("indels per read", "indels_per_read", "{:.2f}"),
            ("reverse-strand reads", "reverse_share", "{:.0%}")]
    out.append(f"{'':<42}{'alt':>8}{'ref':>8}")
    for name, key, fmt in rows:
        out.append(f"{name:<42}{fmt.format(a[key]):>8}{fmt.format(r[key]):>8}")
    if result["flags"]:
        out += ["FLAGS:"] + [f"  - {f}" for f in result["flags"]]
    else:
        out.append("No difference in mismatches, clipping or mapping quality: nothing here suggests mis-mapped reads. "
                   "This does not exclude a PCR or library artifact.")
    out += [f"NOTE: {n}" for n in result.get("notes", [])]
    return "\n".join(out)
