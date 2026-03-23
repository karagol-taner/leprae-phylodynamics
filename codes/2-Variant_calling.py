import os

OUT_DIR = "/content/drive/MyDrive/KCL/Leprosy/results_3"
REF = "/content/drive/MyDrive/KCL/Leprosy/analysis/reference/NC_002677.1.fasta"

# Ensure output directories exist
os.makedirs(f"{OUT_DIR}/bam", exist_ok=True)
os.makedirs(f"{OUT_DIR}/vcf", exist_ok=True)

def run_universal_mapping(input_dir):
    # Find all unique sample accessions in the directory
    all_files = [f for f in os.listdir(input_dir) if f.endswith('.fastq.gz')]
    samples = set([f.split('_')[0].split('.')[0] for f in all_files])

    for s in samples:
        print(f"\n--- Processing: {s} ---")
        r1 = os.path.join(input_dir, f"{s}_1.fastq.gz")
        r2 = os.path.join(input_dir, f"{s}_2.fastq.gz")
        single = os.path.join(input_dir, f"{s}.fastq.gz")

        bam_namesrt = f"{OUT_DIR}/bam/{s}.namesrt.bam"
        bam_fixmate = f"{OUT_DIR}/bam/{s}.fixmate.bam"
        bam_sort = f"{OUT_DIR}/bam/{s}.sorted.bam"
        bam_final = f"{OUT_DIR}/bam/{s}.final.bam"
        vcf_out = f"{OUT_DIR}/vcf/{s}.vcf.gz"

        # Determine Sequencing Mode & Align
        if os.path.exists(r1) and os.path.exists(r2):
            print("  Mode: Paired-end")
            !bwa mem -t 4 {REF} {r1} {r2} | samtools sort -n -o {bam_namesrt}
            !samtools fixmate -m {bam_namesrt} {bam_fixmate}
            !samtools sort -o {bam_sort} {bam_fixmate}
            !rm {bam_namesrt} {bam_fixmate}

        elif os.path.exists(single) or os.path.exists(r1):
            print("  Mode: Single-end (aDNA)")
            # Single-end files don't use fixmate
            target_fastq = single if os.path.exists(single) else r1
            !bwa mem -t 4 {REF} {target_fastq} | samtools sort -o {bam_sort}
        else:
            print(f"  ⚠️ Could not find valid FastQ files for {s}. Skipping.")
            continue

        print("  Marking Duplicates and Calling Variants...")
        # Mark duplicates
        !samtools markdup -r {bam_sort} {bam_final}
        !samtools index {bam_final}
        !rm {bam_sort}

        # Call variants (Haploid mode)
        !bcftools mpileup -Ou -f {REF} {bam_final} | \
         bcftools call -mv --ploidy 1 -Oz -o {vcf_out}
        !bcftools index {vcf_out}

# Execute the pipeline on the pristine data
run_universal_mapping("/content/drive/MyDrive/KCL/Leprosy/Human_SRA/Medieval")
run_universal_mapping("/content/drive/MyDrive/KCL/Leprosy/Human_SRA/Modern")
run_universal_mapping("/content/drive/MyDrive/KCL/Leprosy/ENA")

# Check missing

import os
import time

OUT_DIR = "/content/drive/MyDrive/KCL/Leprosy/results_3"
REF = "/content/drive/MyDrive/KCL/Leprosy/analysis/reference/NC_002677.1.fasta"

os.makedirs(f"{OUT_DIR}/bam", exist_ok=True)
os.makedirs(f"{OUT_DIR}/vcf", exist_ok=True)

def run_universal_mapping(input_dir):
    all_files = [f for f in os.listdir(input_dir) if f.endswith('.fastq.gz')]
    samples = set([f.split('_')[0].split('.')[0] for f in all_files])

    for s in samples:
        vcf_out = f"{OUT_DIR}/vcf/{s}.vcf.gz"
        bam_final = f"{OUT_DIR}/bam/{s}.final.bam"

        # 1. FULL SKIP: If VCF is done, skip entirely
        if os.path.exists(vcf_out):
            print(f"\n--- Skipping {s}: VCF already exists ---")
            continue

        print(f"\n--- Processing: {s} ---")
        start_time = time.time()

        # 2. PARTIAL SKIP: If BAM is done, skip alignment
        if not os.path.exists(bam_final):
            r1 = os.path.join(input_dir, f"{s}_1.fastq.gz")
            r2 = os.path.join(input_dir, f"{s}_2.fastq.gz")
            single = os.path.join(input_dir, f"{s}.fastq.gz")

            bam_namesrt = f"{OUT_DIR}/bam/{s}.namesrt.bam"
            bam_fixmate = f"{OUT_DIR}/bam/{s}.fixmate.bam"
            bam_sort = f"{OUT_DIR}/bam/{s}.sorted.bam"

            # Determine Sequencing Mode & Align
            if os.path.exists(r1) and os.path.exists(r2):
                print("  Mode: Paired-end")
                !bwa mem -t 4 {REF} {r1} {r2} | samtools sort -n -o {bam_namesrt}
                !samtools fixmate -m {bam_namesrt} {bam_fixmate}
                !samtools sort -o {bam_sort} {bam_fixmate}
                !rm -f {bam_namesrt} {bam_fixmate}

            elif os.path.exists(single) or os.path.exists(r1):
                print("  Mode: Single-end (aDNA)")
                target_fastq = single if os.path.exists(single) else r1
                !bwa mem -t 4 {REF} {target_fastq} | samtools sort -o {bam_sort}
            else:
                print(f"  ⚠️ Could not find valid FastQ files for {s}. Skipping.")
                continue

            print("  Marking Duplicates...")
            !samtools markdup -r {bam_sort} {bam_final}
            !samtools index {bam_final}
            !rm -f {bam_sort}
        else:
            print("  BAM already exists. Skipping alignment and jumping to variant calling...")

        # 3. VARIANT CALLING (with the broken pipe fix)
        bcf_temp = f"{OUT_DIR}/vcf/{s}.temp.bcf"

        print("  Generating Pileup...")
        !bcftools mpileup -Ou -f {REF} {bam_final} -o {bcf_temp}

        print("  Calling Variants (Haploid mode)...")
        !bcftools call -mv --ploidy 1 -Oz -o {vcf_out} {bcf_temp}
        !bcftools index -f {vcf_out}

        !rm -f {bcf_temp}

        elapsed = round((time.time() - start_time) / 60, 2)
        print(f"  ✅ Finished {s} in {elapsed} minutes.")

# Execute the pipeline
run_universal_mapping("/content/drive/MyDrive/KCL/Leprosy/Human_SRA/Medieval")
run_universal_mapping("/content/drive/MyDrive/KCL/Leprosy/Human_SRA/Modern")
run_universal_mapping("/content/drive/MyDrive/KCL/Leprosy/ENA")
