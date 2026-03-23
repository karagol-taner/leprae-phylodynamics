import gzip

VCF_FILE = "/content/drive/MyDrive/KCL/Leprosy/results_3/core_snps.vcf.gz"
FASTA_OUT = "/content/drive/MyDrive/KCL/Leprosy/results_3/tree/alignment_fixed.fasta"

print("--- Converting Core SNPs to FASTA (Resolving Reference Alleles) ---")

samples = []
seqs = {}

with gzip.open(VCF_FILE, 'rt') as f:
    for line in f:
        if line.startswith('##'):
            continue
        if line.startswith('#CHROM'):
            samples = line.strip().split('\t')[9:]
            for s in samples:
                seqs[s] = []
            print(f"Found {len(samples)} samples. Reconstructing matrix...")
            continue

        parts = line.strip().split('\t')
        ref = parts[3]
        alt = parts[4].split(',')[0]

        for i, s in enumerate(samples):
            gt_field = parts[9+i].split(':')[0]

            # THE FIX: Treat '.' (missing in a variants-only merge) as Reference
            if gt_field in ['0/0', '0', '.', './.']:
                seqs[s].append(ref)
            # 1/1 or 1 = Alternate mutation
            elif gt_field in ['1/1', '1']:
                seqs[s].append(alt)
            # True missing/fail
            else:
                seqs[s].append('N')

with open(FASTA_OUT, 'w') as out:
    for s in samples:
        out.write(f">{s}\n{''.join(seqs[s])}\n")

# Verify lengths
lengths = set([len(seqs[s]) for s in samples])
print(f"✅ Success! Matrix reconstructed. All sequences are exactly {list(lengths)[0]} base pairs long.")
