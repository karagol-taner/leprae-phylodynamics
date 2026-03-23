from Bio import AlignIO

print("--- Calculating Zoonotic SNP Distance (M. leprae clade) ---")

aln_file = "/content/drive/MyDrive/KCL/Leprosy/results_3/tree/alignment_fixed.fasta"
print("Loading core SNP alignment into memory...")
alignment = AlignIO.read(aln_file, "fasta")

squirrel_seq = None
human_seq = None

for record in alignment:
    # 1. We swap to SRR3672752 (Confirmed Brownsea Island M. leprae)
    if "SRR3672752" in record.id:
        squirrel_seq = str(record.seq).upper()
        print("✅ Found True M. leprae Squirrel (SRR3672752)")
    # 2. Medieval Winchester Human
    elif "SRR847036" in record.id:
        human_seq = str(record.seq).upper()
        print("✅ Found Medieval Winchester Human (SRR847036)")

if squirrel_seq and human_seq:
    mutations = 0
    valid_sites = 0

    for sq_base, hu_base in zip(squirrel_seq, human_seq):
        if sq_base in "ACGT" and hu_base in "ACGT":
            valid_sites += 1
            if sq_base != hu_base:
                mutations += 1

    print("\n" + "="*50)
    print("🔬 FINAL TRANSMISSION METRICS 🔬")
    print("="*50)
    print(f"Total core SNPs compared: {valid_sites:,}")
    print(f"Exact mutations between 11th Century Human and Modern Squirrel: {mutations}")
    print("="*50)
else:
    print("❌ Could not locate sequences.")
    
from Bio import AlignIO

print("--- Calculating Evolutionary Calibration Points ---")

aln_file = "/content/drive/MyDrive/KCL/Leprosy/results_3/tree/alignment_fixed.fasta"
alignment = AlignIO.read(aln_file, "fasta")

# Samples for comparison:
# 1. Modern Human (Global): SRR3330052
# 2. Medieval Human (Winchester): SRR847036
# 3. Modern Squirrel (UK): SRR3672752

seqs = {rec.id.split('/')[-1].replace('.final.bam', ''): str(rec.seq).upper() for rec in alignment}

def get_dist(s1, s2):
    m, v = 0, 0
    for a, b in zip(s1, s2):
        if a in "ACGT" and b in "ACGT":
            v += 1
            if a != b: m += 1
    return m

# A. The 1000-year Human Drift
dist_human_drift = get_dist(seqs["SRR3330052"], seqs["SRR847036"])
# B. Modern Human vs. Modern Squirrel
dist_modern_jump = get_dist(seqs["SRR3330052"], seqs["SRR3672752"])

print(f"\n1000-Year Human Drift (Modern vs Medieval): {dist_human_drift} SNPs")
print(f"Modern Zoonotic Similarity (Modern Human vs Squirrel): {dist_modern_jump} SNPs")
print(f"Historical Link (Medieval Human vs Squirrel): 405 SNPs (Previously calculated)")
