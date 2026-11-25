import random

random.seed(42)

# 1. The "Master" reference sequence (100 bp)
ref_seq = "ATGCGTACGTAGCTAGCTAGCTAGCTAGCTAGCATCGATCGATCGTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAG"

haplotypes = []
for i in range(20):
    mut_seq = list(ref_seq)
    
    # Variant 1: Position 30 (0-indexed 29, A -> G) at 30% frequency (6/20)
    if i < 6: 
        mut_seq[29] = "G"
        
    # Variant 2: Position 60 (0-indexed 59, C -> T) at 10% frequency (2/20)
    if i >= 18: 
        mut_seq[59] = "T"
        
    # Add one random unique SNP (1/20 frequency)
    if i == 5:
        mut_seq[79] = "A"
        
    seq_str = "".join(mut_seq)
    haplotypes.append(f">haplotype_{i+1}\n{seq_str}\n")

# Write to FASTA
with open("simulated_haplotypes.fasta", "w") as f:
    f.writelines(haplotypes)

print("✅ Success! Created simulated_haplotypes.fasta with 20 sequences.")
