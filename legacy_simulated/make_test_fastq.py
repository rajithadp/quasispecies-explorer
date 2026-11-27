import random

ref = "ATGCGTACGTAGCTAGCTAGCTAGCTAGCTAGCATCGATCGATCGTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAG"
reads = []

def rev_comp(seq):
    comp = {'A':'T', 'T':'A', 'C':'G', 'G':'C', 'N':'N'}
    return "".join(comp.get(b, b) for b in reversed(seq))

# Generate 100 reads of 100bp
for i in range(100):
    start = random.randint(0, len(ref)-100)
    read = list(ref[start:start+100])
    
    # Inject variant at position 50 for 20% of reads
    if start <= 50 < start+100 and i < 20:
        read[50-start] = "T" 
        
    read_str = "".join(read)
    
    # Mix forward and reverse strands to pass iVar's strand bias filter
    if i % 2 == 0:
        read_str = rev_comp(read_str)
        
    reads.append(f"@read_{i}\n{read_str}\n+\n{'I'*100}\n")

with open("test_reads.fastq", "w") as f:
    f.writelines(reads)
print("✅ Created test_reads.fastq with 100bp mixed-strand reads")
