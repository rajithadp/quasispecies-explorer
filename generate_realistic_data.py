import pandas as pd
import random

random.seed(42)

# Base sequence (simulating a 200bp region of Spike)
base_seq = "ATGCGTACGTAGCTAGCTAGCTAGCTAGCTAGCATCGATCGATCGTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAG"
# Let's make it exactly 200bp for a realistic feel
base_seq = base_seq * 2 
base_seq = base_seq[:200]

haplotypes = []
metadata = []

# Create 50 haplotypes (representing 50 reads from a patient)
for i in range(50):
    seq = list(base_seq)
    
    # Hotspot 1: 40% frequency (20/50)
    if i < 20:
        seq[50] = "T"
        
    # Hotspot 2: 10% frequency (5/50)
    if 20 <= i < 25:
        seq[120] = "G"
        
    # Hotspot 3: Rare variant (1/50)
    if i == 49:
        seq[180] = "C"
        
    hap_id = f"hap_{i+1}"
    haplotypes.append(f">{hap_id}\n{''.join(seq)}\n")
    
    # Create metadata (e.g., sampling temperature and host)
    # We are simulating a scenario where higher temperature correlates with slightly more diversity
    temp = random.choice([15, 20, 25, 30]) if i % 2 == 0 else random.choice([10, 15, 20])
    host = random.choice(["Human", "Bat", "Human", "Human"])
    metadata.append({"HaplotypeID": hap_id, "Temperature_C": temp, "Host": host})

# Write FASTA
with open("realistic_haplotypes.fasta", "w") as f:
    f.writelines(haplotypes)

# Write Metadata CSV
df_meta = pd.DataFrame(metadata)
df_meta.to_csv("metadata.csv", index=False)

print("✅ Created realistic_haplotypes.fasta and metadata.csv")
