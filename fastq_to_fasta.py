"""Map paired-end reads, call variants, and reconstruct haplotypes.

What changed compared with the first version:
  * no shell=True; every tool is called with an argument list
  * all intermediate files live in a per-run work directory (stale files can't be reused)
  * CliqueSNV haplotype frequencies are KEPT in the output headers (freq=...)
  * if CliqueSNV finds no linked haplotypes, no haplotypes are invented: the filtered
    iVar variants are written for site-level analysis, and the run is recorded in
    run_info.json. The old frequency-sampled "cloud" is available only with
    --simulate-cloud N and is clearly marked synthetic=true in the output headers.
  * primer trimming is optional (--primers FILE, or --artic-v3 for SARS-CoV-2)
  * the random seed is fixed and recorded
"""
import argparse
import json
import shlex
import shutil
import subprocess
import sys
import urllib.request
from pathlib import Path

from quasigraph.io import write_fasta
from quasigraph.pipeline_utils import (
    build_synthetic_cloud,
    convert_cliquesnv_fasta,
    describe_cliquesnv_error,
    parse_ivar_variants,
    read_cliquesnv_json,
    read_first_fasta_sequence,
)

ARTIC_V3_URL = (
    "https://raw.githubusercontent.com/artic-network/artic-ncov2019/master/"
    "primer_schemes/nCoV-2019/V3/nCoV-2019.bed"
)


def run(cmd):
    print("Running:", " ".join(shlex.quote(str(c)) for c in cmd), flush=True)
    subprocess.run([str(c) for c in cmd], check=True)


def run_pipe(cmd1, cmd2):
    """cmd1 | cmd2, without a shell. Raises if either command fails."""
    print("Running:", " ".join(shlex.quote(str(c)) for c in cmd1), "|",
          " ".join(shlex.quote(str(c)) for c in cmd2), flush=True)
    p1 = subprocess.Popen([str(c) for c in cmd1], stdout=subprocess.PIPE)
    p2 = subprocess.Popen([str(c) for c in cmd2], stdin=p1.stdout)
    p1.stdout.close()
    rc2 = p2.wait()
    rc1 = p1.wait()
    if rc1 or rc2:
        raise subprocess.CalledProcessError(rc1 or rc2, cmd1 if rc1 else cmd2)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Map raw FASTQ reads and reconstruct haplotypes.")
    ap.add_argument("--fastq1", required=True)
    ap.add_argument("--fastq2", required=True)
    ap.add_argument("--ref", required=True, help="reference FASTA")
    ap.add_argument("--out", default="real_haplotypes.fasta", help="output haplotype FASTA")
    ap.add_argument("--workdir", help="work directory (default: <out>_work); recreated on every run")
    ap.add_argument("--primers", help="primer BED file for 'ivar trim' (skip trimming if omitted)")
    ap.add_argument("--artic-v3", action="store_true", help="download the ARTIC v3 SARS-CoV-2 primer BED and trim")
    ap.add_argument("--min-variant-freq", type=float, default=0.02, help="minimum ALT_FREQ for site-level output (default 0.02)")
    ap.add_argument("--tf", type=float, default=0.001, help="CliqueSNV -tf: minimum haplotype frequency (default 0.001; lower is slower)")
    ap.add_argument("--skip-cliquesnv", action="store_true",
                    help="do not run CliqueSNV (site-level variants only); useful for many samples or low-memory machines")
    ap.add_argument("--simulate-cloud", type=int, metavar="N", default=0,
                    help="ONLY if CliqueSNV finds no haplotypes: write N SYNTHETIC haplotypes sampled from allele frequencies")
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args(argv)

    out = Path(args.out).resolve()
    work = Path(args.workdir).resolve() if args.workdir else out.with_name(out.stem + "_work")
    if work.exists():
        shutil.rmtree(work)
    work.mkdir(parents=True)
    ref = Path(args.ref).resolve()

    primers = Path(args.primers).resolve() if args.primers else None
    if args.artic_v3:
        primers = work / "primers.bed"
        print(f"Downloading {ARTIC_V3_URL}")
        urllib.request.urlretrieve(ARTIC_V3_URL, primers)

    # 1. map
    bam = work / "mapped.sorted.bam"
    run_pipe(["minimap2", "-ax", "sr", ref, args.fastq1, args.fastq2], ["samtools", "sort", "-o", bam, "-"])
    run(["samtools", "index", bam])

    # 2. trim primers (optional)
    final_bam = bam
    if primers is not None:
        run(["ivar", "trim", "-p", work / "trimmed", "-i", bam, "-b", primers])
        final_bam = work / "trimmed.sorted.bam"
        run(["samtools", "sort", "-o", final_bam, work / "trimmed.bam"])
        run(["samtools", "index", final_bam])
    else:
        print("WARNING: no primer file given; amplicon primers (if any) are NOT trimmed.", file=sys.stderr)

    # 3. site-level variants
    run_pipe(
        ["samtools", "mpileup", "-aa", "-A", "-d", "600000", "-B", "-Q", "0", "--reference", ref, final_bam],
        ["ivar", "variants", "-p", work / "variants", "-q", "20", "-t", "0.01", "-r", ref],
    )
    variants_tsv = work / "variants.tsv"
    variants = parse_ivar_variants(variants_tsv.read_text(), min_freq=args.min_variant_freq)
    print(f"{len(variants)} PASS variants with ALT_FREQ >= {args.min_variant_freq}")

    # 4. linked haplotypes with CliqueSNV
    sam = work / "trimmed.sam"
    run(["samtools", "view", "-h", "-o", sam, final_bam])
    cs_dir = work / "cliquesnv_out"
    if args.skip_cliquesnv:
        print("Skipping CliqueSNV (--skip-cliquesnv).")
    else:
        try:
            run(["cliquesnv", "-m", "snv-illumina", "-in", sam, "-outDir", cs_dir, "-tf", args.tf, "-fdf", "short", "-log"])
        except subprocess.CalledProcessError:
            print("CliqueSNV failed.", file=sys.stderr)

    info = {"tf": args.tf, "min_variant_freq": args.min_variant_freq, "n_site_variants": len(variants),
            "seed": args.seed, "primers_trimmed": primers is not None, "cliquesnv_skipped": args.skip_cliquesnv}
    # CliqueSNV exits normally even when it fails; the failure is recorded only in its JSON output.
    cs_failed = False
    json_files = sorted(cs_dir.glob("*.json")) if cs_dir.exists() else []
    if json_files:
        cs_status = read_cliquesnv_json(json_files[0].read_text())
        if cs_status["failed"]:
            cs_failed = True
            msg = describe_cliquesnv_error(cs_status)
            info["cliquesnv_error"] = msg
            print(f"WARNING: {msg}.\nIts haplotype output is a placeholder and will NOT be used. "
                  "Re-run CliqueSNV by hand on the .sam file in the work directory to see the full error.", file=sys.stderr)
    fasta_files = [] if cs_failed else (sorted(cs_dir.glob("*.fasta")) if cs_dir.exists() else [])
    names = seqs = freqs = None
    if fasta_files:
        names, seqs, freqs = convert_cliquesnv_fasta(fasta_files[0].read_text())

    if seqs and len(seqs) > 1:
        write_fasta(out, names, seqs, freqs)
        info.update(haplotype_source="cliquesnv", n_haplotypes=len(seqs))
        print(f"CliqueSNV found {len(seqs)} haplotypes; written to {out}")
    else:
        sites_out = out.with_suffix(".variants.tsv")
        shutil.copy(variants_tsv, sites_out)
        info.update(haplotype_source="none", n_haplotypes=len(seqs or []))
        print(("CliqueSNV was skipped.\n" if args.skip_cliquesnv else
               "CliqueSNV did not produce valid haplotypes (see the warning above).\n" if cs_failed
               else "CliqueSNV found no linked haplotypes above the threshold.\n") +
              f"Site-level variants are in {sites_out}. Haplotype-level metrics are not available for this sample.")
        if args.simulate_cloud <= 0:
            print("No haplotypes were invented (use --simulate-cloud N only if you want a clearly marked synthetic cloud).")
        else:
            names, syn_seqs, skipped = build_synthetic_cloud(
                read_first_fasta_sequence(ref.read_text()), variants, args.simulate_cloud, args.seed)
            write_fasta(out, names, syn_seqs, synthetic=True)
            info.update(haplotype_source="synthetic", n_haplotypes=len(syn_seqs), skipped_indels=skipped)
            print(f"WARNING: wrote {len(syn_seqs)} SYNTHETIC haplotypes to {out} (headers carry synthetic=true). "
                  "Haplotype counts and rare-haplotype load from this file are not properties of the sample.")

    info_path = out.with_suffix(".run_info.json")
    info_path.write_text(json.dumps(info, indent=2))
    print(f"Run details: {info_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
