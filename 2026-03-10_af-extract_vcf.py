#!/usr/bin/env python3
"""
vcf_to_af_annotated.py

Heavily-annotated converter: VCF (bgzipped) -> samples_short.txt + variants_af.tsv

What this script does (high level)
- Reads a bgzipped VCF (hard-coded path near the top)
- Uses bcftools to list sample names (bcftools query -l)
- Streams per-variant FORMAT fields (AD and DP for each sample) via bcftools query
- Computes per-sample AF = alt_count / (ref_count + alt_count) using AD when present
- Writes:
    - samples_short.txt  (one sample per line)
    - variants_af.tsv    (tab-delimited: CHROM POS ID REF ALT <sample1_AF> <sample2_AF> ...)

Intended usage / constraints
- No CLI — edit the top-of-script hard-coded VCF / paths to point to your files.
- Requires: bcftools on PATH and a VCF with FORMAT:AD populated per sample (best case).
- If AD missing or malformed for a sample/site, writes "NA" for that sample's AF at that variant.
- For multi-allelic sites, this script uses only the first ALT column in the AD parsing logic as written.
  (A more complete implementation would treat each ALT separately and output multiple rows / alt-specific AFs.)

Why this script exists
- Many downstream plotting and correlation scripts expect a simple variant x sample AF matrix (TSV)
  rather than parsing VCFs repeatedly. This script streamlines creation of that matrix.

Major limitations (be aware)
- Uses bcftools query formatting that produces one line per variant where the FORMAT fields
  are expanded inline. The format string is constructed by repeating a "[\t%AD\t%DP]" segment
  per sample — this is brittle if FORMAT ordering or presence varies among samples.
- For large VCFs (> millions of variants) the produced TSV can be very large; consider filtering VCF
  upstream (e.g., bcftools view -i 'TYPE="snp" && QUAL>30' ...) prior to using this script.
- Multi-allelic AD handling is not robust; see improvements below.
"""

# Standard library imports
import subprocess, csv, math, sys, os

# ======= HARDCODED PATHS (EDIT THESE ONLY) =======
# Change these to match your environment before running.
VCF = "/home/jgd/Documents/bioinformatics_working/output/merged_proc-enu_samples.shortnames.vcf.gz"
OUTDIR = "/home/jgd/Documents/bioinformatics_working/output"
SAMPLES_FILE = os.path.join(OUTDIR, "samples_short.txt")
AF_FILE = os.path.join(OUTDIR, "variants_af.tsv")
# ================================================

def die(msg, code=1):
    """Print an error and exit. Small helper to centralize failure behavior."""
    print("ERROR:", msg, file=sys.stderr)
    sys.exit(code)

# -----------------------------
# Preliminary environment checks
# -----------------------------
# Validate the VCF path and that bcftools is available.
if not os.path.exists(VCF):
    die(f"VCF not found: {VCF}")

# Use shutil.which to detect bcftools on PATH. Import here to avoid top-level import requirement.
import shutil
if not shutil.which("bcftools"):
    die("bcftools not found on PATH. Install or add bcftools to PATH.")

# Ensure output directory exists
os.makedirs(OUTDIR, exist_ok=True)

# -----------------------------
# 1) Extract sample names (one-per-line)
# -----------------------------
# We use "bcftools query -l" which lists sample IDs from the VCF header in the order they appear.
# This reproduces the sample order we will use when streaming FORMAT fields.
try:
    p = subprocess.run(["bcftools", "query", "-l", VCF], check=True, stdout=subprocess.PIPE, text=True)
except subprocess.CalledProcessError as e:
    die(f"bcftools query -l failed: {e}")

samples = [s.strip() for s in p.stdout.splitlines() if s.strip()]
if not samples:
    die("No samples found in VCF header (bcftools query -l returned empty).")

# Write the sample list to a simple text file for downstream consumers.
with open(SAMPLES_FILE, "w") as sfh:
    for s in samples:
        sfh.write(s + "\n")
print("Wrote sample list:", SAMPLES_FILE, file=sys.stderr)

# -----------------------------
# 2) Stream AD/DP fields and build AF matrix
# -----------------------------
# Plan:
# - Build a bcftools query format string that prints CHROM,POS,ID,REF,ALT, then for each sample prints AD and DP
# - Invoke bcftools query and read stdout line-by-line
# - Parse the per-sample AD fields (expected format: REF,ALT[,ALT2,...]) and compute AF = ALT/(REF+ALT)
# - Write TSV rows: CHROM POS ID REF ALT <AF_sample1> <AF_sample2> ...

# Notes about the format string:
# - We use a repeated segment '[\t%AD\t%DP]' * len(samples) to emit AD and DP for each sample.
#   The outer code expects these placeholders to expand into textual AD and DP values for each sample in order.
# - This approach relies on bcftools query producing AD/DP even when missing (produces '.' or empty) so that the per-sample
#   field positions align. However, FORMAT availability can vary; if AD absent for a sample, we will see '.' and produce NA.
#
# Possible alternative approaches (more robust):
# - Use cyvcf2 or pysam in Python to parse per-sample FORMAT arrays explicitly (recommended for complex VCFs or multi-allelic handling).
# - Use "bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT[\t%AD\t%DP]\n'" (note the format with brackets used by bcftools to iterate samples).
#   The script below uses an equivalent approach (string concatenation).
fmt = '%CHROM\t%POS\t%ID\t%REF\t%ALT' + '[\\t%AD\\t%DP]' * len(samples) + '\\n'
# Explanation of fmt:
# - %CHROM, %POS, %ID, %REF, %ALT produce the first five metadata columns.
# - The bracket expression [\t%AD\t%DP] instructs bcftools to output AD then DP for each sample in sample order.
#   Repeating it len(samples) times ensures we have as many AD/DP tokens as samples.
#
# Construct command
cmd = ["bcftools", "query", "-f", fmt, VCF]

# Launch the subprocess and open the output file for writing
proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True, bufsize=1)

# Prepare TSV writer
with open(AF_FILE, "w", newline='') as outfh:
    writer = csv.writer(outfh, delimiter='\t', lineterminator='\n')
    # write header: metadata columns + sample IDs
    writer.writerow(["CHROM","POS","ID","REF","ALT"] + samples)

    line_count = 0
    # We read bcftools stdout line-by-line to keep memory low (streaming)
    for line in proc.stdout:
        line_count += 1
        parts = line.rstrip("\n").split("\t")
        # sanity: a valid bcftools query line should have at least 5 metadata columns
        if len(parts) < 5:
            # skip malformed / header-like lines (defensive)
            continue
        chrom, pos, vid, ref, alt = parts[:5]
        rest = parts[5:]
        af_row = [chrom, pos, vid, ref, alt]

        # rest is expected to contain AD, DP, AD, DP, ... for each sample in sample order.
        # We'll parse AD token and ignore DP for AF computation (but DP is available for diagnostics).
        # We iterate by sample; 'i' advances by 2 for AD and DP tokens.
        i = 0
        for si in range(len(samples)):
            ad_field = rest[i] if i < len(rest) else '.'
            alt_count = None
            ref_count = None
            # Parse AD if present and not '.' sentinel
            if ad_field and ad_field != '.':
                # AD is often formatted as "ref,alt" or "ref,alt1,alt2,..." for multi-allelic sites.
                parts_ad = ad_field.split(',')
                try:
                    # We attempt to parse numeric counts robustly; '.' elements become None
                    if len(parts_ad) >= 2:
                        ref_count = float(parts_ad[0]) if parts_ad[0] != '.' else None
                        # For biallelic sites typical alt_count is parts_ad[1]. If multiple ALTs exist,
                        # this code uses the first ALT (parts_ad[1]) — see 'Limitations' below.
                        alt_count = float(parts_ad[1]) if parts_ad[1] != '.' else None
                except Exception:
                    # If parsing fails (non-numeric content), treat as missing
                    alt_count = None

            af_val = None
            # Compute AF when alt_count is available; denominator uses ref_count when available
            if alt_count is not None:
                denom = (ref_count if ref_count is not None else 0.0) + alt_count
                if denom > 0:
                    af_val = alt_count / denom
            # Format AF numeric with 6 decimal places or 'NA' if missing
            af_row.append("{:.6f}".format(af_val) if (af_val is not None and not math.isnan(af_val)) else "NA")
            i += 2  # advance past AD and DP tokens

        writer.writerow(af_row)

        # Progress report every 10k variants — helpful for long runs
        if (line_count % 10000) == 0:
            print(f"Processed {line_count} variants...", file=sys.stderr)

# Wait for bcftools process to exit and capture its return code
ret = proc.wait()
if ret != 0:
    # We continue on non-zero exit but warn the user — often bcftools returns 0 even on partial outputs.
    print("Warning: bcftools query exited with status", ret, file=sys.stderr)

print("Wrote AF matrix:", AF_FILE, file=sys.stderr)

# -----------------------------
# End of script — Additional notes & recommended improvements
# -----------------------------
#
# Improvements you should consider (non-exhaustive)
# 1) Robust multi-allelic handling:
#    - Current script reads AD and uses only AD[1] (first ALT). For sites with multiple ALTs you may want to:
#      * Expand multi-ALT variants into multiple rows (one per ALT) with per-ALT AF.
#      * Or compute combined AF across all ALTs (sum of all alt ADs / (ref + sum alts)) if that suits your downstream analysis.
#
# 2) Use a VCF parser (cyvcf2, pysam) rather than shelling out to bcftools for more robust handling:
#    - cyvcf2 allows direct access to variant.format("AD") as arrays and simplifies per-sample looping.
#    - Example advantage: simpler handling of missing FORMAT fields and multi-allelic AD arrays.
#
# 3) Guard against inconsistent FORMAT fields:
#    - Some VCFs may not include AD at all, or may include per-sample tags in a different order.
#    - Consider verifying the FORMAT column per-variant or using bcftools query’s per-sample bracket notation properly:
#         bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT[\t%AD]\n' file.vcf.gz
#      which automatically repeats %AD for each sample in the file.
#
# 4) Performance & filtering:
#    - For very large VCFs, filter upstream (e.g., restrict to SNPs, PASS, chromosome regions) to reduce output size.
#    - Consider compressing the AF output (gzip) or storing selected variants only.
#
# 5) Output formats:
#    - The resulting variants_af.tsv has string "NA" for missing AFs. Downstream scripts should treat "NA" as missing.
#    - You might prefer empty cells or numeric NaN; using pandas to read and write would allow consistent numeric NA handling.
#
# 6) Header and sample ordering:
#    - The sample order in the header is the same as bcftools query -l output. Make sure downstream tools expect this order.
#    - If you want shortnames or map sample IDs to other labels, produce a sample map file and use it to rename columns.
#
# 7) Error handling:
#    - This script prints progress and warns on bcftools non-zero exit.
#    - For production pipelines, add logging, try/except around the bcftools process, and optionally retry logic.
#
# 8) Unit tests:
#    - Create a small mock VCF (few variants, few samples) with explicit AD/DP to validate the script behavior and edge cases:
#       * biallelic site with AD present
#       * site with missing AD for some samples
#       * multi-allelic site with AD containing 3+ values
#       * indel or MNP (should be skipped if REF/ALT lengths != 1 depending on your use-case)
#
# 9) Use tabix/bcftools view to prefilter:
#    - Example:
#         bcftools view -f PASS -v snps input.vcf.gz -Oz -o filtered.vcf.gz
#         bcftools index filtered.vcf.gz
#      Then point this script at filtered.vcf.gz to reduce noise.
#
# If you want, I can:
# - Convert this into a more robust Python implementation using cyvcf2 (supports AD arrays, multi-allelic, indexed access).
# - Add CLI argument parsing (argparse) and options for alt-handling (first-alt vs expand-multiallelic vs sum-alt).
# - Implement streaming gzip output and optional per-chromosome partitioning for large VCFs.
#
# End of annotations.
