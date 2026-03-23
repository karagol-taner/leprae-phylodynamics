import os
import subprocess
import time

# ============================================================================
# STEP 1: INSTALL SRA TOOLKIT
# ============================================================================

print("\n[1/4] Installing SRA Toolkit...")

try:
    # Check if already installed
    result = subprocess.run(["prefetch", "--version"],
                          capture_output=True, check=False)
    if result.returncode == 0:
        print("  ✓ SRA Toolkit already installed")
    else:
        raise FileNotFoundError
except FileNotFoundError:
    print("  Installing SRA Toolkit...")

    # Download and install
    subprocess.run([
        "wget", "-q",
        "https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz"
    ], check=True)

    subprocess.run(["tar", "-xzf", "sratoolkit.current-ubuntu64.tar.gz"], check=True)

    # Move binaries to /usr/local/bin
    subprocess.run([
        "bash", "-c",
        "mv sratoolkit.*/bin/* /usr/local/bin/"
    ], check=True)

    # Verify installation
    result = subprocess.run(["prefetch", "--version"], capture_output=True)
    if result.returncode == 0:
        print("  ✓ SRA Toolkit installed successfully")
    else:
        print("  ✗ Installation failed")
        exit(1)

# ============================================================================
# STEP 2: INSTALL PYSRADB
# ============================================================================

print("\n[2/4] Installing pysradb...")
subprocess.run(["pip", "install", "-q", "pysradb"], check=True)
print("  ✓ pysradb installed")

# ============================================================================
# STEP 3: SETUP DIRECTORIES
# ============================================================================

print("\n[3/4] Setting up directories...")

HUMAN_MED_DIR = "/content/drive/MyDrive/KCL/Leprosy/Human_SRA/Medieval"
HUMAN_MOD_DIR = "/content/drive/MyDrive/KCL/Leprosy/Human_SRA/Modern"

os.makedirs(HUMAN_MED_DIR, exist_ok=True)
os.makedirs(HUMAN_MOD_DIR, exist_ok=True)

print(f"  ✓ Medieval: {HUMAN_MED_DIR}")
print(f"  ✓ Modern: {HUMAN_MOD_DIR}")

# ============================================================================
# STEP 4: GET ACCESSION LISTS
# ============================================================================

print("\n[4/4] Fetching accession lists...")

import pysradb
db = pysradb.SRAweb()

medieval_df = db.sra_metadata('PRJNA200950', detailed=True)
all_medieval = medieval_df['run_accession'].unique().tolist()

modern_df = db.sra_metadata('PRJNA317287', detailed=True)
all_modern = modern_df['run_accession'].unique().tolist()

print(f"  ✓ Medieval: {len(all_medieval)} total samples")
print(f"  ✓ Modern: {len(all_modern)} total samples")

# ============================================================================
# CHECK WHAT'S COMPLETED
# ============================================================================

def get_completed(output_dir):
    """Find completed downloads"""
    completed = set()
    if os.path.exists(output_dir):
        for f in os.listdir(output_dir):
            if f.endswith('.fastq.gz') or f.endswith('.fastq'):
                acc = f.split('_')[0].split('.')[0]
                completed.add(acc)
    return completed

medieval_done = get_completed(HUMAN_MED_DIR)
modern_done = get_completed(HUMAN_MOD_DIR)

print(f"\nCurrent status:")
print(f"  Medieval: {len(medieval_done)}/{len(all_medieval)} complete")
print(f"  Modern: {len(modern_done)}/{len(all_modern)} complete")

# Get remaining
medieval_remaining = [acc for acc in all_medieval if acc not in medieval_done]
modern_remaining = [acc for acc in all_modern if acc not in modern_done]

print(f"\nTo download:")
print(f"  Medieval: {len(medieval_remaining)} samples")
print(f"  Modern: {len(modern_remaining)} samples")

if not medieval_remaining and not modern_remaining:
    print("\n✅ ALL SAMPLES ALREADY DOWNLOADED!")
    exit(0)

# ============================================================================
# DOWNLOAD FUNCTION
# ============================================================================

def download_direct(accession, output_dir):
    """Download directly to Google Drive"""

    try:
        # 1. Prefetch
        result = subprocess.run(
            ["prefetch", accession, "-O", output_dir],
            capture_output=True,
            text=True,
            timeout=1200
        )
        if result.returncode != 0:
            return False, "prefetch_failed"

        # 2. Convert to FASTQ
        result = subprocess.run(
            ["fasterq-dump", accession,
             "-O", output_dir,
             "--split-files",
             "--threads", "2"],
            capture_output=True,
            text=True,
            timeout=3600  # 1 hour
        )
        if result.returncode != 0:
            return False, "conversion_failed"

        # 3. Compress
        compressed = False
        for suffix in ["_1.fastq", "_2.fastq", ".fastq"]:
            fastq = os.path.join(output_dir, f"{accession}{suffix}")
            if os.path.exists(fastq):
                subprocess.run(["gzip", fastq], timeout=900, check=False)
                compressed = True

        if not compressed:
            return False, "no_output"

        # 4. Cleanup
        sra_dir = os.path.join(output_dir, accession)
        if os.path.exists(sra_dir):
            subprocess.run(["rm", "-rf", sra_dir], check=False)

        return True, "success"

    except subprocess.TimeoutExpired:
        # Cleanup on timeout
        sra_dir = os.path.join(output_dir, accession)
        if os.path.exists(sra_dir):
            subprocess.run(["rm", "-rf", sra_dir], check=False)
        return False, "timeout"
    except Exception as e:
        return False, str(e)[:50]

# ============================================================================
# DOWNLOAD MEDIEVAL REMAINING
# ============================================================================

if medieval_remaining:
    print(f"\n{'='*70}")
    print(f"DOWNLOADING {len(medieval_remaining)} MEDIEVAL SAMPLES")
    print(f"{'='*70}\n")

    med_success = 0
    med_failed = 0

    for i, acc in enumerate(medieval_remaining, 1):
        print(f"[{i}/{len(medieval_remaining)}] {acc}...", end=" ", flush=True)
        success, reason = download_direct(acc, HUMAN_MED_DIR)

        if success:
            print("✓")
            med_success += 1
        else:
            print(f"✗ ({reason})")
            med_failed += 1

        time.sleep(2)

    print(f"\nMedieval: ✓ {med_success} | ✗ {med_failed}")

# ============================================================================
# DOWNLOAD MODERN REMAINING
# ============================================================================

if modern_remaining:
    print(f"\n{'='*70}")
    print(f"DOWNLOADING {len(modern_remaining)} MODERN SAMPLES")
    print(f"{'='*70}\n")

    mod_success = 0
    mod_failed = 0

    for i, acc in enumerate(modern_remaining, 1):
        print(f"[{i}/{len(modern_remaining)}] {acc}...", end=" ", flush=True)
        success, reason = download_direct(acc, HUMAN_MOD_DIR)

        if success:
            print("✓")
            mod_success += 1
        else:
            print(f"✗ ({reason})")
            mod_failed += 1

        time.sleep(2)

        # Progress update every 10
        if i % 10 == 0:
            print(f"\n  Progress: {i}/{len(modern_remaining)} | ✓ {mod_success} | ✗ {mod_failed}\n")

    print(f"\nModern: ✓ {mod_success} | ✗ {mod_failed}")

# ============================================================================
# FINAL SUMMARY
# ============================================================================

print(f"\n{'='*70}")
print("FINAL SUMMARY")
print(f"{'='*70}\n")

# Re-check completion
medieval_final = get_completed(HUMAN_MED_DIR)
modern_final = get_completed(HUMAN_MOD_DIR)

print(f"Medieval: {len(medieval_final)}/{len(all_medieval)} complete")
print(f"Modern: {len(modern_final)}/{len(all_modern)} complete")

# Storage
def get_size_gb(directory):
    total = 0
    if os.path.exists(directory):
        for root, dirs, files in os.walk(directory):
            for f in files:
                fp = os.path.join(root, f)
                if os.path.exists(fp):
                    total += os.path.getsize(fp)
    return total / (1024**3)

med_size = get_size_gb(HUMAN_MED_DIR)
mod_size = get_size_gb(HUMAN_MOD_DIR)

print(f"\nStorage:")
print(f"  Medieval: {med_size:.2f} GB")
print(f"  Modern: {mod_size:.2f} GB")
print(f"  Total: {med_size + mod_size:.2f} GB")

# Check if complete
if len(medieval_final) == len(all_medieval) and len(modern_final) == len(all_modern):
    print(f"\n🎉 SUCCESS! ALL {len(medieval_final) + len(modern_final)} SAMPLES DOWNLOADED!")
else:
    missing = (len(all_medieval) - len(medieval_final)) + (len(all_modern) - len(modern_final))
    print(f"\n⚠️  {missing} samples still missing - may need manual retry")

print("\n✅ Download complete!")
