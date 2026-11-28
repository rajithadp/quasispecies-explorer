"""A variation graph built directly from a reference and its called variants (no haplotype reconstruction needed).

Model
  * Backbone: the reference sequence, cut at every variable position.
  * Bubble: at each variable position one node per allele: the reference allele and each alternative allele
    called at >= min_freq. A node carries the allele's frequency (iVar ALT_FREQ for alternatives; one minus the
    sum of all alternative frequencies for the reference allele) and a role: 'major' (the most frequent allele)
    or 'minor'. A bubble whose major allele is an alternative is a fixed difference from the reference.
  * Edges join every node of one column to every node of the next column.
  * Paths: 'reference' (reference alleles everywhere), 'consensus' (major allele everywhere), and one path per
    read-supported allele pattern of a linked group of sites (see linkage.py), tagged with its fragment count.

What the weights mean. A frequency on a node is the marginal frequency of that allele at that site. The graph
does NOT say how often any complete path occurs: without haplotypes the combination of alleles across distant
sites is unknown, so paths other than reference, consensus and the read-supported ones carry no frequency.
Evidence paths cover only the sites of their group, and are written only when no other variable site lies between
the first and last site of the group (otherwise the walk would have to guess an allele there).

Only SNVs are represented; indels are counted and skipped.
"""
from __future__ import annotations

from dataclasses import dataclass, field


_OXIDATIVE = {("G", "T"), ("C", "A")}


@dataclass
class VNode:
    id: str
    kind: str                 # 'backbone' or 'allele'
    start: int                # 1-based, inclusive
    end: int
    seq: str
    freq: float | None = None
    role: str | None = None   # 'major' / 'minor' for alleles
    matches_reference: bool | None = None
    annotation: str = ""      # gene, amino-acid change and effect of this allele alone (see graph_annotation.py)


@dataclass
class EvidencePath:
    name: str
    sites: list               # positions
    pattern: str              # e.g. 'AAA' (A = alternative, R = reference)
    fragments: int
    nodes: list               # node ids of the walk, first to last site of the group
    joint: str = ""           # amino-acid changes of the alternative alleles of this path together


@dataclass
class VariantGraph:
    reference: str
    region: str
    nodes: dict = field(default_factory=dict)       # id -> VNode, in genome order
    edges: list = field(default_factory=list)       # (u, v)
    columns: list = field(default_factory=list)     # list of lists of node ids, in genome order
    paths: dict = field(default_factory=dict)       # name -> node ids (reference, consensus)
    evidence: list = field(default_factory=list)    # EvidencePath
    skipped: dict = field(default_factory=dict)
    unplaced: list = field(default_factory=list)    # (sites, reason)

    def sequence(self, node_ids) -> str:
        return "".join(self.nodes[n].seq for n in node_ids)

    def stats(self) -> dict:
        alleles = [n for n in self.nodes.values() if n.kind == "allele"]
        bubbles = [c for c in self.columns if c and self.nodes[c[0]].kind == "allele"]
        return {
            "nodes": len(self.nodes), "edges": len(self.edges), "bubbles": len(bubbles),
            "fixed_difference_bubbles": sum(
                1 for c in bubbles if any(self.nodes[n].role == "major" and not self.nodes[n].matches_reference for n in c)),
            "minor_alleles": sum(1 for n in alleles if n.role == "minor" and not n.matches_reference),
            "evidence_paths": len(self.evidence), "unplaced_groups": len(self.unplaced),
        }


def _alleles_at(position: int, ref_base: str, alts: dict) -> list:
    """[(base, freq, matches_reference)] with the reference allele first."""
    ref_freq = max(0.0, 1.0 - sum(alts.values()))
    return [(ref_base, ref_freq, True)] + [(b, f, False) for b, f in alts.items()]


def build_variant_graph(reference: str, variants, min_freq: float = 0.02, min_depth: int = 100,
                        region: str = "ref", evidence=None) -> VariantGraph:
    """Graph from the reference and an iVar table (as read by read_ivar_variants).

    evidence: optional list of (sites, patterns) with sites = [Site(pos, ref, alt), ...] and patterns = Counter of
    strings over 'A' (alt), 'R' (reference), 'O' (other) and '-' (not observed), as produced by analyze_group.
    """
    ref = reference.upper()
    vg = VariantGraph(reference=ref, region=region)

    # 1. alleles per variable position
    per_pos: dict[int, dict] = {}
    skipped = {"indels": 0, "below_min_depth": 0, "no_alt_at_min_freq": 0, "reference_mismatch": 0, "other_base": 0}
    all_alt: dict[int, float] = {}
    for row in variants.itertuples():
        if row.is_indel or len(row.ref) != 1 or len(row.alt) != 1:
            skipped["indels"] += 1
            continue
        if row.depth < min_depth:
            skipped["below_min_depth"] += 1
            continue
        pos = int(row.position)
        if pos < 1 or pos > len(ref) or ref[pos - 1] != row.ref.upper():
            skipped["reference_mismatch"] += 1
            continue
        if row.alt.upper() not in "ACGT":
            skipped["other_base"] += 1
            continue
        all_alt[pos] = all_alt.get(pos, 0.0) + float(row.freq)
        if row.freq >= min_freq:
            per_pos.setdefault(pos, {})[row.alt.upper()] = float(row.freq)
    skipped["no_alt_at_min_freq"] = sum(1 for p in all_alt if p not in per_pos)
    vg.skipped = skipped

    # 2. columns: backbone segments and bubbles
    prev: list = []
    cursor = 1

    def connect(col):
        for u in prev:
            for v in col:
                vg.edges.append((u, v))

    def add_backbone(a, b):
        nonlocal prev
        if b < a:
            return
        n = VNode(f"b{a}_{b}", "backbone", a, b, ref[a - 1:b])
        vg.nodes[n.id] = n
        vg.columns.append([n.id])
        connect([n.id])
        prev = [n.id]

    for pos in sorted(per_pos):
        add_backbone(cursor, pos - 1)
        alleles = _alleles_at(pos, ref[pos - 1], per_pos[pos])
        major = max(alleles, key=lambda a: (a[1], not a[2]))[0]     # a tie goes to the alternative (as in sites.py: >= 0.5 is fixed)
        col = []
        for base, f, is_ref in alleles:
            n = VNode(f"a{pos}{base}", "allele", pos, pos, base, f, "major" if base == major else "minor", is_ref)
            vg.nodes[n.id] = n
            col.append(n.id)
        vg.columns.append(col)
        connect(col)
        prev = col
        cursor = pos + 1
    add_backbone(cursor, len(ref))

    # 3. reference and consensus paths (walk one node per column)
    vg.paths["reference"] = [c[0] for c in vg.columns]          # the reference allele is always first in a bubble
    vg.paths["consensus"] = [c[0] if vg.nodes[c[0]].kind == "backbone"
                             else next(n for n in c if vg.nodes[n].role == "major") for c in vg.columns]

    # 4. read-supported evidence paths
    bubble_positions = sorted(per_pos)
    for sites, patterns in (evidence or []):
        positions = [s.pos for s in sites]
        lo, hi = min(positions), max(positions)
        missing = [s.pos for s in sites if s.pos not in per_pos or s.alt not in per_pos[s.pos]]
        if missing:
            vg.unplaced.append((positions, f"site(s) {missing} are not bubbles in this graph (frequency or depth below the thresholds)"))
            continue
        between = [p for p in bubble_positions if lo < p < hi and p not in positions]
        if between:
            vg.unplaced.append((positions, f"other variable site(s) {between} lie between the group's sites"))
            continue
        for pattern, count in sorted(patterns.items(), key=lambda kv: -kv[1]):
            if len(pattern) != len(sites) or any(c not in "AR" for c in pattern):
                continue
            chosen = {s.pos: (s.alt if c == "A" else ref[s.pos - 1]) for s, c in zip(sites, pattern)}
            walk = []
            for col in vg.columns:
                n0 = vg.nodes[col[0]]
                if n0.kind == "backbone":
                    if lo < n0.start and n0.end < hi:
                        walk.append(col[0])
                else:
                    if lo <= n0.start <= hi:
                        walk.append(f"a{n0.start}{chosen[n0.start]}")
            vg.evidence.append(EvidencePath(f"linked_{'_'.join(map(str, positions))}_{pattern}", positions,
                                            pattern, int(count), walk))
    return vg


def to_gfa(vg: VariantGraph) -> str:
    """GFA 1.0. Segment tags: SN region, SO 0-based offset, LN length; for alleles FR frequency, rl role,
    mr (1 if the allele equals the reference), ox (1 on a minor G>T / C>A allele, the artifact-prone class), an (gene:amino-acid change:effect of the allele alone,
    when annotated). Linked paths carry fc (fragments) and, when annotated, jn (joint amino-acid changes). Links into an allele carry FR = that allele's frequency.
    Paths: reference, consensus, and linked_* with fc = number of read fragments."""
    lines = ["H\tVN:Z:1.0"]
    for n in vg.nodes.values():
        tags = f"SN:Z:{vg.region}\tSO:i:{n.start - 1}\tLN:i:{len(n.seq)}"
        if n.kind == "allele":
            tags += f"\tFR:f:{n.freq:.6f}\trl:Z:{n.role}\tmr:i:{int(bool(n.matches_reference))}"
            if n.role == "minor" and not n.matches_reference and (vg.reference[n.start - 1], n.seq) in _OXIDATIVE:
                tags += "\tox:i:1"   # minor G>T / C>A allele: the artifact-prone class (see compare.py)
            if n.annotation:
                tags += f"\tan:Z:{n.annotation}"
        lines.append(f"S\t{n.id}\t{n.seq}\t{tags}")
    for u, v in vg.edges:
        t = vg.nodes[v]
        extra = f"\tFR:f:{t.freq:.6f}" if t.kind == "allele" else ""
        lines.append(f"L\t{u}\t+\t{v}\t+\t0M{extra}")
    for name, walk in vg.paths.items():
        lines.append(f"P\t{name}\t{','.join(n + '+' for n in walk)}\t*")
    for e in vg.evidence:
        lines.append(f"P\t{e.name}\t{','.join(n + '+' for n in e.nodes)}\t*\tfc:i:{e.fragments}"
                     + (f"\tjn:Z:{e.joint}" if e.joint else ""))
    return "\n".join(lines) + "\n"
