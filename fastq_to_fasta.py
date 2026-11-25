import subprocess
import argparse
import os
import random
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def run_cmd(cmd):
    print(f"Running: {cmd}")
    subprocess.run(cmd, shell=True, check=True)

def main():
    parser = argparse.ArgumentParser(description="Map raw FASTQ reads and reconstruct haplotypes.")
    parser.add_argument("--fastq", required=True, help="Path to raw FASTQ file (single or R1)")
    parser.add_argument("--ref", required=True, help="Path to reference FASTA")
    parser.add_argument("--out", default="reconstructed_haplotypes.fasta", help="Output FASTA file")
    parser.add_argument("--num_haps", type=int, default=100, help="Number of haplotypes to simulate")
    args = parser.parse_args()

    # 1. Map reads with minimap2
    bam_file = "mapped.sorted.bam"
    run_cmd(f"minimap2 -ax sr {args.ref} {args.fastq} | samtools sort -o {bam_file} -")
    run_cmd(f"samtools index {bam_file}")

    # 2. Call variants with iVar
    # We use a very low threshold (1%) to capture the minor variants (quasispecies)
    run_cmd(f"samtools mpileup -aa -A -d 600000 -B -Q 0 --reference {args.ref} {bam_file} | ivar variants -p variants -q 20 -t 0.01 -r {args.ref}")

    # 3. Load reference and variants
    ref_record = next(SeqIO.parse(args.ref, "fasta"))
    ref_seq = list(str(ref_record.seq))
    
    # Read the iVar TSV
    if not os.path.exists("variants.tsv"):
        print("No variants found. Exiting.")
        return
        
    variants = pd.read_csv("variants.tsv", sep="\t")
    # Filter for variants that passed iVar's statistical tests
    variants = variants[variants['PASS'] == True]

    print(f"Found {len(variants)} valid variants. Reconstructing {args.num_haps} haplotypes...")

    # 4. Reconstruct haplotypes based on allele frequencies
    haplotypes = []
    for i in range(args.num_haps):
        hap_seq = ref_seq.copy()
        for _, var in variants.iterrows():
            pos = int(var['POS']) - 1 # Convert to 0-indexed
            alt_allele = var['ALT']
            freq = float(var['ALT_FREQ'])
            
            # Randomly apply the mutation based on its frequency
            if random.random() < freq:
                hap_seq[pos] = alt_allele
                
        haplotypes.append(SeqRecord(Seq("".join(hap_seq)), id=f"read_{i+1}", description=""))

    # 5. Write output
    SeqIO.write(haplotypes, args.out, "fasta")
    print(f"✅ Success! Written to {args.out}")

if __name__ == "__main__":
    main()
