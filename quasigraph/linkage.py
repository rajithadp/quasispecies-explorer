"""Are two variant alleles on the same read fragment? (physical linkage from reads)

Per-site frequencies cannot tell whether nearby variants sit on the same viral
genome. Reads can: if the alternative alleles of two sites are found together on
one read (or on the two mates of one read pair) far more often than chance, they
are linked, i.e. one minor haplotype carries both.

Input is a SAM file (or BAM via ``samtools view``) and the iVar variants table.
For every fragment (read name) we record which allele it shows at each site of
interest, merging the two mates and ignoring bases below a quality cutoff.

Limits: only sites within one read pair's span (about the insert size, a few
hundred bp) can be compared; co-occurrence on a fragment is evidence of linkage
but a shared mapping or PCR artifact can mimic it, and the verdicts below are
heuristics, not proof.
"""
from __future__ import annotations

import math
import re
import subprocess
from bisect import bisect_left
from collections import Counter
from dataclasses import dataclass
from itertools import combinations

_CIGAR = re.compile(r"(\d+)([MIDNSHP=X])")
SAM_EXCLUDE = 4 | 256 | 512 | 2048  # unmapped, secondary, QC fail, supplementary


@dataclass(frozen=True)
class Site:
    pos: int  # 1-based
    ref: str
    alt: str


def parse_cigar(cigar: str):
    return [(int(n), op) for n, op in _CIGAR.findall(cigar)]


def reference_end(pos: int, ops) -> int:
    """First reference position after the alignment (exclusive)."""
    return pos + sum(n for n, op in ops if op in "MDN=X")


def bases_at_sites(pos: int, ops, seq: str, qual: str, sites_sorted, min_baseq: int = 20) -> dict:
    """Base this read shows at each requested reference position (CIGAR-aware).

    Positions inside a deletion, an intron skip or an unaligned (soft-clipped)
    stretch are simply absent. Bases below ``min_baseq`` are dropped.
    """
    out = {}
    rpos, qpos = pos, 0
    for length, op in ops:
        if op in "M=X":
            i = bisect_left(sites_sorted, rpos)
            while i < len(sites_sorted) and sites_sorted[i] < rpos + length:
                s = sites_sorted[i]
                q = qpos + (s - rpos)
                if q < len(seq) and (qual == "*" or ord(qual[q]) - 33 >= min_baseq):
                    out[s] = seq[q].upper()
                i += 1
            rpos += length
            qpos += length
        elif op in "IS":
            qpos += length
        elif op in "DN":
            rpos += length
    return out


def collect_fragments(lines, site_positions, region=None, min_baseq: int = 20, min_mapq: int = 0, meta=None) -> dict:
    """{read name: {position: base or None}}; None marks mates that disagree (ambiguous).

    If a dict is passed as ``meta`` it is filled with, per fragment, the fragment extent
    (start, length) and the reads (flag, position, mapq, cigar) that showed a base at a tested site.
    """
    sites_sorted = sorted(set(site_positions))
    if not sites_sorted:
        return {}
    frags: dict = {}
    for line in lines:
        if not line or line[0] == "@":
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 11:
            continue
        if int(f[1]) & SAM_EXCLUDE:
            continue
        if region is not None and f[2] != region:
            continue
        if int(f[4]) < min_mapq or f[5] == "*":
            continue
        pos = int(f[3])
        ops = parse_cigar(f[5])
        if reference_end(pos, ops) <= sites_sorted[0] or pos > sites_sorted[-1]:
            continue
        obs = bases_at_sites(pos, ops, f[9], f[10], sites_sorted, min_baseq)
        if not obs:
            continue
        if meta is not None:
            pnext, tlen = int(f[7]), int(f[8])
            extent = (min(pos, pnext), abs(tlen)) if (tlen != 0 and pnext > 0) else (pos, reference_end(pos, ops) - pos)
            meta.setdefault(f[0], {"extent": extent, "reads": []})["reads"].append((int(f[1]), pos, int(f[4]), f[5], dict(obs)))
        d = frags.setdefault(f[0], {})
        for p, b in obs.items():
            if p not in d:
                d[p] = b
            elif d[p] is not None and d[p] != b:
                d[p] = None
    return frags


def alignment_lines(path: str):
    """Yield SAM text lines from a .sam file, or from a .bam/.cram through samtools."""
    if path.lower().endswith((".bam", ".cram")):
        try:
            proc = subprocess.Popen(["samtools", "view", "-F", "0xB04", path], stdout=subprocess.PIPE, text=True)
        except FileNotFoundError as e:
            raise RuntimeError("samtools was not found; give a .sam file or install samtools") from e
        yield from proc.stdout
        if proc.wait() != 0:
            raise RuntimeError(f"samtools view failed on {path}")
    else:
        with open(path) as fh:
            yield from fh


# ---------- statistics ----------

def fisher_exact_2x2(a: int, b: int, c: int, d: int) -> float:
    """Two-sided Fisher exact test p-value for [[a, b], [c, d]] (hypergeometric, via log-gamma)."""
    n = a + b + c + d
    r1, c1 = a + b, a + c
    if n == 0:
        return 1.0

    def logp(x):
        return (math.lgamma(r1 + 1) + math.lgamma(n - r1 + 1) + math.lgamma(c1 + 1) + math.lgamma(n - c1 + 1)
                - math.lgamma(n + 1) - math.lgamma(x + 1) - math.lgamma(r1 - x + 1)
                - math.lgamma(c1 - x + 1) - math.lgamma(n - r1 - c1 + x + 1))

    lo, hi = max(0, r1 + c1 - n), min(r1, c1)
    obs = logp(a)
    total = sum(math.exp(logp(x)) for x in range(lo, hi + 1) if logp(x) <= obs + 1e-7)
    return min(1.0, total)


def pair_verdict(aa: int, ar: int, ra: int, rr: int, p_value: float) -> str:
    n = aa + ar + ra + rr
    if n == 0:
        return "no fragment covers both sites (too far apart for read pairs)"
    a1, a2 = aa + ar, aa + ra
    if a1 + a2 == 0:
        return "no alternative allele observed on these fragments"
    expected = a1 * a2 / n
    p21 = aa / a1 if a1 else 0.0
    p12 = aa / a2 if a2 else 0.0
    if aa >= 3 and p21 >= 0.9 and p12 >= 0.9:
        return "LINKED: the alternative alleles travel together"
    if aa >= 3 and (p21 >= 0.9 or p12 >= 0.9):
        return "NESTED: one allele almost always occurs together with the other"
    if aa == 0:
        if expected >= 5:
            return f"MUTUALLY EXCLUSIVE: never together although ~{expected:.1f} co-occurrences expected by chance"
        return "no co-occurrence, but the alleles are too rare to tell whether they are linked"
    if p_value < 0.01 and aa > expected:
        return f"partly linked: co-occur more often than chance (expected {expected:.1f})"
    return "consistent with independent alleles (not linked)"


@dataclass
class PairResult:
    site1: Site
    site2: Site
    aa: int  # alt at both
    ar: int  # alt at site1 only
    ra: int  # alt at site2 only
    rr: int  # ref at both
    p_site2_given_1: float
    p_site1_given_2: float
    r_squared: float
    fisher_p: float
    verdict: str

    @property
    def n_both(self) -> int:
        return self.aa + self.ar + self.ra + self.rr


def _states(frag: dict, group: list) -> list:
    """A = alt, R = ref, O = another base, - = not observed (or ambiguous), per site of the group."""
    st = []
    for s in group:
        b = frag.get(s.pos)
        if b is None:
            st.append("-")
        elif b == s.alt:
            st.append("A")
        elif b == s.ref:
            st.append("R")
        else:
            st.append("O")
    return st


def fragment_summary(frags: dict, meta: dict, group: list) -> dict:
    """Facts about the reads behind alt-carrying versus reference-only fragments (artifact check)."""
    alt, ref = [], []
    for q, frag in frags.items():
        st = _states(frag, group)
        if "A" in st:
            alt.append(q)
        elif "R" in st and "O" not in st:
            ref.append(q)

    def stats(qs):
        reads = [r for q in qs for r in meta[q]["reads"]]
        return {
            "fragments": len(qs),
            "distinct_extents": len({meta[q]["extent"] for q in qs}),
            "starts": Counter(meta[q]["extent"][0] for q in qs),
            "both_mates_overlap": sum(len(meta[q]["reads"]) >= 2 for q in qs),
            "reads": len(reads),
            "reverse_reads": sum(bool(r[0] & 16) for r in reads),
            "duplicate_flagged": sum(bool(r[0] & 1024) for r in reads),
        }

    rows = [(q, "".join(_states(frags[q], group)), meta[q]["extent"], meta[q]["reads"]) for q in alt]
    return {"alt": stats(alt), "ref": stats(ref), "alt_rows": rows}


def analyze_group(frags: dict, group: list, meta=None) -> dict:
    """Per-site counts, haplotype patterns and pairwise linkage for one group of sites."""
    patterns: Counter = Counter()
    per_site = {s.pos: Counter() for s in group}
    states = []
    for frag in frags.values():
        st = _states(frag, group)
        if all(c == "-" for c in st):
            continue
        states.append(st)
        for s, c in zip(group, st):
            per_site[s.pos][c] += 1
        patterns["".join(st)] += 1

    pairs = []
    for i, j in combinations(range(len(group)), 2):
        c = Counter()
        for st in states:
            if st[i] in "AR" and st[j] in "AR":
                c[st[i] + st[j]] += 1
        aa, ar, ra, rr = c["AA"], c["AR"], c["RA"], c["RR"]
        a1, a2 = aa + ar, aa + ra
        n = aa + ar + ra + rr
        denom = a1 * (ra + rr) * a2 * (ar + rr)
        r2 = ((aa * rr - ar * ra) ** 2 / denom) if denom else float("nan")
        p = fisher_exact_2x2(aa, ar, ra, rr)
        pairs.append(PairResult(group[i], group[j], aa, ar, ra, rr,
                                aa / a1 if a1 else float("nan"), aa / a2 if a2 else float("nan"),
                                r2, p, pair_verdict(aa, ar, ra, rr, p)))
    res = {"group": group, "per_site": per_site, "patterns": patterns, "pairs": pairs}
    if meta is not None:
        res["fragment_summary"] = fragment_summary(frags, meta, group)
        disc, disc_rows = {}, []
        for q, frag in frags.items():
            for s in group:
                if s.pos in frag and frag[s.pos] is None:  # the two mates show different bases here
                    bases = [(r[0], r[4].get(s.pos, "?")) for r in meta[q]["reads"]]
                    d = disc.setdefault(s.pos, {"fragments": 0, "alt_on_one_mate": 0})
                    d["fragments"] += 1
                    d["alt_on_one_mate"] += any(b == s.alt for _, b in bases)
                    disc_rows.append((q, s.pos, bases))
        res["discordant"] = disc
        res["discordant_rows"] = disc_rows
    return res


# ---------- choosing sites from the variants table ----------

def best_snv_per_position(variants) -> dict:
    """{position: Site} using, at each position, the SNV allele with the highest frequency."""
    sites = {}
    best = {}
    for row in variants.itertuples():
        if row.is_indel or len(row.ref) != 1 or len(row.alt) != 1:
            continue
        if row.position not in best or row.freq > best[row.position]:
            best[row.position] = row.freq
            sites[row.position] = Site(int(row.position), row.ref.upper(), row.alt.upper())
    return sites


def parse_site_spec(spec: str, known: dict) -> list:
    """'5743,5744' or '5743:C,5744:C' -> [Site, ...]; the allele after ':' overrides the table's."""
    group = []
    for token in spec.split(","):
        token = token.strip()
        if not token:
            continue
        pos_s, _, alt = token.partition(":")
        pos = int(pos_s)
        if pos not in known:
            raise ValueError(f"position {pos} is not an SNV in the variants table")
        s = known[pos]
        group.append(Site(pos, s.ref, alt.upper() if alt else s.alt))
    if len(group) < 2:
        raise ValueError(f"a site group needs at least two positions: {spec!r}")
    return group


def auto_groups(sites_table, known: dict, min_freq: float, window: int) -> list:
    """Chains of polymorphic sites (minor frequency >= min_freq) whose neighbours are <= window bp apart."""
    keep = sites_table[(~sites_table["fixed_difference"]) & (sites_table["minor_freq"] >= min_freq)]
    positions = sorted(int(p) for p in keep["position"] if int(p) in known)
    groups, current = [], []
    for p in positions:
        if current and p - current[-1] > window:
            groups.append(current)
            current = []
        current.append(p)
    if current:
        groups.append(current)
    return [[known[p] for p in g] for g in groups if len(g) >= 2]


# ---------- reporting ----------

def format_report(results: list, n_fragments: int, top_patterns: int = 8) -> str:
    out = [f"{n_fragments} read fragments overlap at least one tested site.\n"]
    for k, res in enumerate(results, 1):
        g = res["group"]
        out.append(f"=== Group {k}: " + ", ".join(f"{s.pos} {s.ref}>{s.alt}" for s in g) + " ===")
        out.append("Per site (fragments showing alt / covering the site):")
        for s in g:
            c = res["per_site"][s.pos]
            cov = c["A"] + c["R"]
            frac = f"{c['A'] / cov:.2%}" if cov else "n/a"
            out.append(f"  {s.pos}: {c['A']} / {cov}  ({frac})")
        if "discordant" in res:
            if res["discordant"]:
                out.append("Fragments where the two mates show different bases (excluded from the counts above):")
                for pos, d in sorted(res["discordant"].items()):
                    out.append(f"  {pos}: {d['fragments']} fragment(s), {d['alt_on_one_mate']} with the alt allele on one mate only")
            else:
                out.append("Fragments where the two mates show different bases: none")
        out.append("Allele patterns (A = alt, R = ref, O = other, - = not covered), most common first:")
        for pat, n in res["patterns"].most_common(top_patterns):
            out.append(f"  {pat}  {n}")
        fs = res.get("fragment_summary")
        if fs:
            for label, key in (("Alt-carrying fragments", "alt"), ("Reference-only fragments", "ref")):
                x = fs[key]
                if x["fragments"]:
                    out.append(f"{label}: n={x['fragments']} | distinct (start,length)={x['distinct_extents']} | both mates over the sites: "
                               f"{x['both_mates_overlap']}/{x['fragments']} | reverse-strand reads: {x['reverse_reads']}/{x['reads']} | "
                               f"duplicate-flagged reads: {x['duplicate_flagged']}")
        if fs and fs["alt"]["fragments"]:
            out.append("Fragment start positions (a different start usually means a different amplicon):")
            out.append("  alt-carrying:   " + _fmt_starts(fs["alt"]["starts"], 8))
            out.append("  reference-only: " + _fmt_starts(fs["ref"]["starts"], 8) + " (most common first)")
        out.append("Pairs:")
        for p in res["pairs"]:
            out.append(f"  {p.site1.pos} & {p.site2.pos}: both alt={p.aa}, alt1 only={p.ar}, alt2 only={p.ra}, "
                       f"neither={p.rr}  | P(2|1)={_fmt(p.p_site2_given_1)}  P(1|2)={_fmt(p.p_site1_given_2)}  "
                       f"Fisher p={p.fisher_p:.2g}")
            out.append(f"      -> {p.verdict}")
        out.append("")
    out.append("Heuristic verdicts. Co-occurrence on a fragment suggests physical linkage but does not rule out a shared "
               "mapping or PCR artifact; check strand balance and amplicon boundaries before concluding.")
    return "\n".join(out)


def _fmt_starts(counter, k: int) -> str:
    if not counter:
        return "none"
    items = counter.most_common(k)
    text = ", ".join(f"{start} x{n}" for start, n in items)
    return text + (f", ... (+{len(counter) - k} more starts)" if len(counter) > k else "")


def _fmt(x: float) -> str:
    return "n/a" if x != x else f"{x:.2f}"


def write_tables(prefix: str, results: list) -> list:
    paths = [f"{prefix}.pairs.tsv", f"{prefix}.patterns.tsv"]
    with open(paths[0], "w") as fh:
        fh.write("group\tpos1\tpos2\tboth_alt\talt1_only\talt2_only\tneither\tp2_given_1\tp1_given_2\tr_squared\tfisher_p\tverdict\n")
        for k, res in enumerate(results, 1):
            for p in res["pairs"]:
                fh.write(f"{k}\t{p.site1.pos}\t{p.site2.pos}\t{p.aa}\t{p.ar}\t{p.ra}\t{p.rr}\t"
                         f"{_fmt(p.p_site2_given_1)}\t{_fmt(p.p_site1_given_2)}\t{_fmt(p.r_squared)}\t{p.fisher_p:.4g}\t{p.verdict}\n")
    if any("fragment_summary" in r for r in results):
        path = f"{prefix}.alt_fragments.tsv"
        paths.append(path)
        with open(path, "w") as fh:
            fh.write("group\tfragment\tpattern\textent_start\textent_length\treads_over_sites\treads(flag:pos:mapq:cigar)\n")
            for k, res in enumerate(results, 1):
                for q, pat, (start, length), reads in res.get("fragment_summary", {}).get("alt_rows", []):
                    fh.write(f"{k}\t{q}\t{pat}\t{start}\t{length}\t{len(reads)}\t"
                             + ";".join(f"{r[0]}:{r[1]}:{r[2]}:{r[3]}" for r in reads) + "\n")
    if any(r.get("discordant_rows") for r in results):
        path = f"{prefix}.discordant.tsv"
        paths.append(path)
        with open(path, "w") as fh:
            fh.write("group\tfragment\tsite\tbase_per_read(flag:base)\n")
            for k, res in enumerate(results, 1):
                for q, pos, bases in res.get("discordant_rows", []):
                    fh.write(f"{k}\t{q}\t{pos}\t" + ";".join(f"{fl}:{b}" for fl, b in bases) + "\n")
    with open(paths[1], "w") as fh:
        fh.write("group\tsites\tpattern\tfragments\n")
        for k, res in enumerate(results, 1):
            sites = ",".join(str(s.pos) for s in res["group"])
            for pat, n in res["patterns"].most_common():
                fh.write(f"{k}\t{sites}\t{pat}\t{n}\n")
    return paths
