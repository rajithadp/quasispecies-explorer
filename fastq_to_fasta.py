import subprocess
import argparse
import os
import random
import glob
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def run_cmd(cmd):
    print(f"Running: {cmd}")
    subprocess.run(cmd, shell=True, check=True)

def main():
    parser = argparse.ArgumentParser(description="Map raw FASTQ reads and reconstruct haplotypes.")
    parser.add_argument("--fastq1", required=True, help="Path to raw FASTQ R1")
    parser.add_argument("--fastq2", required=True, help="Path to raw FASTQ R2")
    parser.add_argument("--ref", required=True, help="Path to reference FASTA")
    parser.add_argument("--out", default="real_haplotypes.fasta", help="Output FASTA")
    parser.add_argument("--num_haps", type=int, default=100, help="Number of haplotypes to simulate if fallback is used")
    args = parser.parse_args()

    # 1. Download primers (ARTIC V3) if not present
    if not os.path.exists("primers.bed"):
        run_cmd("curl -L -o primers.bed https://raw.githubusercontent.com/artic-network/artic-ncov2019/master/primer_schemes/nCoV-2019/V3/nCoV-2019.bed")

    # 2. Map paired-end reads with minimap2
    bam_file = "mapped.sorted.bam"
    run_cmd(f"minimap2 -ax sr {args.ref} {args.fastq1} {args.fastq2} | samtools sort -o {bam_file} -")
    run_cmd(f"samtools index {bam_file}")

    # 3. Trim primers with iVar
    run_cmd(f"ivar trim -p trimmed -i {bam_file} -b primers.bed")
    run_cmd(f"samtools sort -o trimmed.sorted.bam trimmed.bam")
    run_cmd(f"samtools index trimmed.sorted.bam")

    # 4. Call variants with iVar (for fallback)
    print("Running iVar variants...")
    run_cmd(f"samtools mpileup -aa -A -d 600000 -B -Q 0 --reference {args.ref} trimmed.sorted.bam | ivar variants -p variants -q 20 -t 0.01 -r {args.ref}")

    # 5. Try CliqueSNV
    sam_file = "trimmed.sam"
    run_cmd(f"samtools view -h -o {sam_file} trimmed.sorted.bam")
    output_dir = "cliquesnv_out"
    
    print("Running CliqueSNV...")
    try:
        run_cmd(f"cliquesnv -m snv-illumina -in {sam_file} -outDir {output_dir} -tf 0.001")
    except subprocess.CalledProcessError:
        print("CliqueSNV failed. Proceeding to fallback.")

    # 6. Check CliqueSNV output
    fasta_files = glob.glob(os.path.join(output_dir, "*.fasta"))
    use_fallback = True

    if len(fasta_files) > 0:
        hap_file = fasta_files[0]
        records = list(SeqIO.parse(hap_file, "fasta"))
        
        # If CliqueSNV found more than 1 sequence, use it!
        if len(records) > 1:
            print(f"✅ CliqueSNV successfully found {len(records)} haplotypes.")
            for i, rec in enumerate(records):
                rec.id = f"hap_{i+1}"
                rec.description = ""
            SeqIO.write(records, args.out, "fasta")
            use_fallback = False

    # 7. Fallback to Frequency-Based Simulation
    if use_fallback:
        print("⚠️ CliqueSNV returned only 1 sequence. Falling back to frequency-based reconstruction...")
        
        ref_record = next(SeqIO.parse(args.ref, "fasta"))
        ref_seq = list(str(ref_record.seq))
        
        if os.path.exists("variants.tsv"):
            variants = pd.read_csv("variants.tsv", sep="\t")
            variants = variants[variants['PASS'] == True]
            
            # FILTER: Only keep variants above 2% to prevent haplotype explosion
            variants = variants[variants['ALT_FREQ'] >= 0.02]
            
            print(f"Reconstructing {args.num_haps} haplotypes based on {len(variants)} valid iVar variants (>2%)...")
            haplotypes = []
            
            for i in range(args.num_haps):
                hap_seq = ref_seq.copy()
                for _, var in variants.iterrows():
                    pos = int(var['POS']) - 1
                    alt = var['ALT']
                    freq = float(var['ALT_FREQ'])
                    
                    if random.random() < freq:
                        hap_seq[pos] = alt
                        
                haplotypes.append(SeqRecord(Seq("".join(hap_seq)), id=f"hap_{i+1}", description=""))
                
            SeqIO.write(haplotypes, args.out, "fasta")
            print(f"✅ Successfully created {args.num_haps} simulated haplotypes. Written to {args.out}")

if __name__ == "__main__":
    main()