"""Annotate single-nucleotide variants with the gene, codon and amino-acid effect they cause.

Input: a GFF3 file (CDS features), the reference sequence, and variants (position, REF, ALT).
Handles CDS made of several segments (e.g. SARS-CoV-2 orf1ab, which is a join of two segments with a
-1 ribosomal frameshift: the base at the join is read twice and so can have two codon contexts) and genes
on the minus strand.

Effects: synonymous, missense, stop_gained, stop_lost, noncoding (no CDS covers the position), indel
(not a substitution), incomplete_codon, unknown (codon contains a non-ACGT base).
Each effect is the effect of this one substitution on its own codon, with the reference as background.
Codon numbers count within the CDS as annotated (for orf1ab: the whole polyprotein, not nsp numbering).
"""
from __future__ import annotations

from dataclasses import dataclass
from urllib.parse import unquote

_BASES = "TCAG"
_AAS = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON_TABLE = {a + b + c: _AAS[16 * i + 4 * j + k]
               for i, a in enumerate(_BASES) for j, b in enumerate(_BASES) for k, c in enumerate(_BASES)}
_COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}


@dataclass(frozen=True)
class Cds:
    cds_id: str
    gene: str
    product: str
    strand: str
    segments: tuple  # ((start, end), ...) 1-based inclusive, in file order

    def coding_positions(self) -> list[int]:
        """Genomic positions in the order the codons are read."""
        if self.strand == "-":
            out: list[int] = []
            for s, e in sorted(self.segments, reverse=True):
                out += range(e, s - 1, -1)
            return out
        out = []
        for s, e in self.segments:
            out += range(s, e + 1)
        return out

    @property
    def length(self) -> int:
        return sum(e - s + 1 for s, e in self.segments)


def _attrs(field: str) -> dict:
    out = {}
    for part in field.strip().split(";"):
        if "=" in part:
            k, v = part.split("=", 1)
            out[k] = unquote(v)
    return out


def parse_gff(text: str) -> list[Cds]:
    """CDS features of a GFF3 file; rows sharing an ID are the segments of one CDS."""
    order: list[str] = []
    seg: dict[str, list] = {}
    meta: dict[str, tuple] = {}
    for line in text.splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "CDS":
            continue
        strand, phase = f[6], f[7]
        if strand not in "+-":
            raise ValueError(f"CDS with unsupported strand '{strand}': {line[:80]}")
        if phase not in ("0", "."):
            raise ValueError(f"CDS with phase {phase} is not supported (only phase 0): {line[:80]}")
        a = _attrs(f[8])
        cid = a.get("ID") or f"cds_{len(order)}"
        if cid not in seg:
            order.append(cid)
            seg[cid] = []
            meta[cid] = (a.get("gene") or a.get("Name") or cid, a.get("product", ""), strand)
        seg[cid].append((int(f[3]), int(f[4])))
    return [Cds(c, *meta[c][:1], meta[c][1], meta[c][2], tuple(seg[c])) for c in order]


class Annotator:
    def __init__(self, cds_list: list[Cds], reference: str):
        self.cds_list = cds_list
        self.ref = reference.upper()
        self._by_pos: dict[int, list] = {}
        self._coding: dict[str, list[int]] = {}
        for c in cds_list:
            pos = c.coding_positions()
            self._coding[c.cds_id] = pos
            for idx, p in enumerate(pos):
                self._by_pos.setdefault(p, []).append((c, idx))

    def _base(self, pos: int, strand: str, override: str | None = None) -> str:
        b = override if override is not None else self.ref[pos - 1]
        return _COMP.get(b, "N") if strand == "-" else b

    def annotate(self, pos: int, ref: str, alt: str) -> list[dict]:
        """One dict per codon context of the substitution (usually exactly one)."""
        if pos < 1 or pos > len(self.ref):
            raise ValueError(f"position {pos} is outside the reference (length {len(self.ref)})")
        if self.ref[pos - 1] != ref.upper():
            raise ValueError(f"reference mismatch at {pos}: variants table says {ref}, reference FASTA has "
                             f"{self.ref[pos - 1]} (is this the same reference the reads were mapped to?)")
        if len(alt) != 1 or alt.upper() not in "ACGT" or alt.startswith(("+", "-")):
            return [self._row(None, None, "indel")]
        hits = self._by_pos.get(pos)
        if not hits:
            return [self._row(None, None, "noncoding")]
        rows = []
        for c, idx in hits:
            coding = self._coding[c.cds_id]
            start = idx - idx % 3
            codon_pos = coding[start:start + 3]
            if len(codon_pos) < 3:
                rows.append(self._row(c, idx, "incomplete_codon"))
                continue
            ref_codon = "".join(self._base(p, c.strand) for p in codon_pos)
            alt_codon = "".join(self._base(p, c.strand, alt.upper() if p == pos else None) for p in codon_pos)
            if "N" in ref_codon or "N" in alt_codon:
                rows.append(self._row(c, idx, "unknown"))
                continue
            ra, aa = CODON_TABLE[ref_codon], CODON_TABLE[alt_codon]
            effect = ("synonymous" if ra == aa else "stop_gained" if aa == "*" else
                      "stop_lost" if ra == "*" else "missense")
            rows.append(self._row(c, idx, effect, ref_codon, alt_codon, ra, aa))
        return rows

    @staticmethod
    def _row(c, idx, effect, ref_codon="", alt_codon="", ra="", aa="") -> dict:
        codon = (idx // 3 + 1) if c is not None and idx is not None else None
        return {
            "gene": c.gene if c else "", "product": c.product if c else "",
            "codon": codon, "codon_position": (idx % 3 + 1) if c is not None and idx is not None else None,
            "ref_codon": ref_codon, "alt_codon": alt_codon, "ref_aa": ra, "alt_aa": aa,
            "aa_change": f"{ra}{codon}{aa}" if ra and aa else "", "effect": effect,
        }


# ---------- tables ----------

def annotate_table(variants, annotator: "Annotator", min_freq: float = 0.02, min_depth: int = 100):
    """Annotate every SNV with frequency >= min_freq and depth >= min_depth (iVar rows, one per alt allele).

    'kind' is fixed_difference (frequency >= 0.5) or minor. 'GT_CA' flags the artifact-prone class
    (see compare.py). A base read in two codon contexts gives two rows.
    """
    import pandas as pd

    from .compare import _is_class

    v = variants[(~variants["is_indel"]) & (variants["freq"] >= min_freq) & (variants["depth"] >= min_depth)]
    flags = _is_class(v)
    rows = []
    for (_, r), gt_ca in zip(v.iterrows(), flags):
        for ann in annotator.annotate(int(r["position"]), r["ref"], r["alt"]):
            rows.append({"position": int(r["position"]), "ref": r["ref"], "alt": r["alt"], "freq": float(r["freq"]),
                         "depth": int(r["depth"]), "kind": "fixed_difference" if r["freq"] >= 0.5 else "minor",
                         "GT_CA": bool(gt_ca), **ann})
    cols = ["position", "ref", "alt", "freq", "depth", "kind", "GT_CA", "gene", "codon", "codon_position",
            "ref_codon", "alt_codon", "aa_change", "effect", "product"]
    return pd.DataFrame(rows, columns=cols)


def gene_summary(table, cds_list: list[Cds]):
    """Minor variants per gene (and noncoding), all and without G>T / C>A, with effect counts.

    Counts are of variants, not of independent sites at risk, so they are not dN/dS: compare genes
    only with this in mind (per_kb divides by the CDS length).
    """
    import pandas as pd

    minor = table[table["kind"] == "minor"]
    genes = [c.gene for c in cds_list]
    rows = []
    for g in genes + [""]:
        sub = minor[minor["gene"] == g]
        length = next((c.length for c in cds_list if c.gene == g), None)
        keep = sub[~sub["GT_CA"]]
        rows.append({
            "gene": g or "(noncoding)", "length_nt": length,
            "minor_variants": len(sub), "missense": int((sub["effect"] == "missense").sum()),
            "synonymous": int((sub["effect"] == "synonymous").sum()),
            "stop_gained": int((sub["effect"] == "stop_gained").sum()),
            "minor_variants_without_GT_CA": len(keep),
            "missense_without_GT_CA": int((keep["effect"] == "missense").sum()),
            "synonymous_without_GT_CA": int((keep["effect"] == "synonymous").sum()),
            "per_kb": (len(sub) / length * 1000) if length else float("nan"),
        })
    return pd.DataFrame(rows)


# ---------- linked variants annotated as one haplotype ----------

def annotate_haplotype(annotator: "Annotator", substitutions) -> list[dict]:
    """Effect of several substitutions that sit on the SAME genome, one row per affected codon.

    substitutions: iterable of (position, ref, alt) single-base changes. Changes that fall in one codon are applied
    together. 'single_effects' lists what each change would do alone; 'emergent' is True when the joint amino acid
    is not the amino acid any single change gives (the combination does something no member does). Changes in
    noncoding sequence get one row each (gene empty). A base read in two codon contexts (a frameshift join) is
    reported in both.
    """
    subs = [(int(p), r.upper(), a.upper()) for p, r, a in substitutions]
    if len({p for p, _, _ in subs}) != len(subs):
        raise ValueError("two substitutions at the same position cannot be on the same haplotype")
    for p, r, a in subs:
        if len(a) != 1 or a not in "ACGT":
            raise ValueError(f"only single-base substitutions can be combined, got {p} {r}>{a}")
        if p < 1 or p > len(annotator.ref):
            raise ValueError(f"position {p} is outside the reference (length {len(annotator.ref)})")
        if annotator.ref[p - 1] != r:
            raise ValueError(f"reference mismatch at {p}: given {r}, reference FASTA has {annotator.ref[p - 1]}")

    codons: dict = {}   # (cds_id, codon_start_index) -> (cds, [(pos, ref, alt, idx)])
    rows: list[dict] = []
    for p, r, a in subs:
        hits = annotator._by_pos.get(p)
        if not hits:
            rows.append({**Annotator._row(None, None, "noncoding"), "positions": str(p), "single_effects": "noncoding",
                         "emergent": False})
            continue
        for c, idx in hits:
            codons.setdefault((c.cds_id, idx - idx % 3), (c, []))[1].append((p, r, a, idx))
    for (_, start), (c, members) in codons.items():
        coding = annotator._coding[c.cds_id]
        codon_pos = coding[start:start + 3]
        positions = ",".join(str(p) for p, _, _, _ in sorted(members))
        if len(codon_pos) < 3:
            rows.append({**Annotator._row(c, start, "incomplete_codon"), "positions": positions,
                         "single_effects": "incomplete_codon", "emergent": False})
            continue
        change = {p: a for p, _, a, _ in members}
        ref_codon = "".join(annotator._base(p, c.strand) for p in codon_pos)
        alt_codon = "".join(annotator._base(p, c.strand, change.get(p)) for p in codon_pos)
        if "N" in ref_codon or "N" in alt_codon:
            rows.append({**Annotator._row(c, start, "unknown"), "positions": positions, "single_effects": "unknown",
                         "emergent": False})
            continue
        ra, aa = CODON_TABLE[ref_codon], CODON_TABLE[alt_codon]
        effect = ("synonymous" if ra == aa else "stop_gained" if aa == "*" else "stop_lost" if ra == "*" else "missense")
        singles = []
        for p, _, a, _ in sorted(members):
            sc = "".join(annotator._base(q, c.strand, a if q == p else None) for q in codon_pos)
            singles.append(f"{ra}{start // 3 + 1}{CODON_TABLE[sc]}")
        single_aas = {s[-1] for s in singles}
        rows.append({**Annotator._row(c, start, effect, ref_codon, alt_codon, ra, aa), "positions": positions,
                     "single_effects": "; ".join(singles), "emergent": len(members) > 1 and aa not in single_aas})
    rows.sort(key=lambda r: (int(r["positions"].split(",")[0])))
    return rows
