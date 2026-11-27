"""Command line demo: python -m quasigraph [--fasta aligned.fasta] [--out results]"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from .analysis import analyze
from .io import read_fasta, write_fasta
from .plot import draw_graph
from .simulate import simulate_quasispecies


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description="Quantify a quasispecies from aligned haplotypes.")
    p.add_argument("--fasta", help="aligned haplotypes (optional header token freq=0.25); simulated if omitted")
    p.add_argument("--out", default="results", help="output folder (default: results)")
    p.add_argument("--seed", type=int, default=0, help="seed for the simulation")
    p.add_argument("--length", type=int, default=300, help="simulated genome length")
    p.add_argument("--haplotypes", type=int, default=12, help="number of simulated haplotypes")
    args = p.parse_args(argv)

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    if args.fasta:
        names, seqs, freqs = read_fasta(args.fasta)
    else:
        sim = simulate_quasispecies(length=args.length, n_haplotypes=args.haplotypes, seed=args.seed)
        names, seqs, freqs = sim.names, sim.seqs, sim.freqs
        write_fasta(out / "haplotypes.fasta", names, seqs, freqs)
        print(f"Simulated hotspots (1-based): {[int(h) + 1 for h in sim.hotspots]}")

    res = analyze(seqs, freqs)
    (out / "summary.json").write_text(json.dumps(res["summary"].as_dict(), indent=2))
    res["site_entropy"].to_csv(out / "site_entropy.csv", index=False)
    res["top_sites"].to_csv(out / "top_sites.csv", index=False)
    fig = draw_graph(res["graph"])
    fig.savefig(out / "graph.png", dpi=150)

    print(json.dumps(res["summary"].as_dict(), indent=2))
    print("\nMost variable sites:")
    print(res["top_sites"].to_string(index=False))
    print(f"\nGraph: {res['graph'].stats()}")
    print(f"Files written to {out.resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
