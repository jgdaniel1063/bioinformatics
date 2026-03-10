#!/usr/bin/env bash
#
# bcftools_mpileup_local_annotated.sh
#
# Local, non-cluster version of the bcftools mpileup + call per-sample script.
# All SLURM and cluster-specific lines have been removed. This runs sequentially
# over sample subdirectories under ROOT_DIR and writes per-sample compressed VCFs.
#
# Usage:
#   Edit the ROOT_DIR, REF_FASTA, and OUTPUT_DIR variables below, then run:
#     bash bcftools_mpileup_local_annotated.sh
#
# Important notes:
#  - This script runs on a single machine. For many samples consider parallelizing
#    externally (GNU parallel, xargs -P, or a job scheduler).
#  - The script expects one BAM per sample directory (it takes the first *.bam found).
#  - It ensures BAMs are indexed (.bai) and creates an index with samtools index if missing.
#  - Outputs are written as bgzipped VCF (.vcf.gz) and tabix-indexed.
#  - The script uses an intermediate BCF for efficient piping: mpileup -> bcftools call -> bcf -> vcf.gz
#
# Safety & reproducibility recommendations:
#  - Verify reference FASTA (.fa) and its index (.fai) match the contig names in BAMs.
#  - Capture tool versions (bcftools, samtools, tabix) for provenance (this script writes a versions.txt).
#  - Check BAM files are coordinate-sorted and consistent (samtools flagstat can help).
#
set -euo pipefail
umask 002

# -------------------------
# User-editable paths
# -------------------------
ROOT_DIR="/nfs/turbo/umms-jshavit/jgdaniel/06-03-2024_proc-enu_paper/11-19-2025_final_results/gatk_workflow/star_alignment"
REF_FASTA="/nfs/turbo/umms-jshavit/jgdaniel/reference/danRer11.fa"
OUTPUT_DIR="/nfs/turbo/umms-jshavit/jgdaniel/06-03-2024_proc-enu_paper/raw_output/bcftools"  # Top-level folder for bcftools results

# Tools (assume on PATH)
BCFTOOLS_BIN="${BCFTOOLS_BIN:-bcftools}"
SAMTOOLS_BIN="${SAMTOOLS_BIN:-samtools}"
TABIX_BIN="${TABIX_BIN:-tabix}"

# Optional tuning: mpileup/call parameters you may want to adjust
MPILEUP_OPTS="-f ${REF_FASTA} -Ou"    # -Ou to output uncompressed BCF for piping
CALL_OPTS="-mv -Ob"                  # -m multi-allelic caller, -v variants-only, -Ob output BCF

# Create output dir
mkdir -p "${OUTPUT_DIR}"

# -------------------------
# Sanity checks
# -------------------------
for tool in "$BCFTOOLS_BIN" "$SAMTOOLS_BIN" "$TABIX_BIN"; do
  if ! command -v "$tool" >/dev/null 2>&1; then
    echo "[ERROR] Required tool not found on PATH: $tool" >&2
    exit 1
  fi
done

# Record tool versions for provenance
VERSIONS_FILE="${OUTPUT_DIR}/tools_versions.txt"
{
  echo "Run: $(date -u +"%Y-%m-%d %H:%M:%SZ")"
  echo "bcftools: $($BCFTOOLS_BIN --version 2>&1 | head -n1)"
  echo "samtools: $($SAMTOOLS_BIN --version 2>&1 | head -n1)"
  echo "tabix: $($TABIX_BIN --version 2>&1 | head -n1 || true)"
} > "${VERSIONS_FILE}"

# Ensure reference FASTA index exists (.fai). Create with samtools faidx if missing.
if [[ ! -s "${REF_FASTA}.fai" ]]; then
  echo "[INFO] Reference FASTA index not found, creating: ${REF_FASTA}.fai"
  "${SAMTOOLS_BIN}" faidx "${REF_FASTA}"
fi

# -------------------------
# Discover sample directories (deterministic ordering)
# -------------------------
# Use sort to ensure stable ordering. For reproducible mapping of samples to jobs,
# a precomputed sample list (one-per-line) is preferred; here we auto-discover.
mapfile -t SAMPLE_DIRS < <(find "${ROOT_DIR}" -mindepth 1 -maxdepth 1 -type d -print | sort)

TOTAL_SAMPLES=${#SAMPLE_DIRS[@]}
echo "[INFO] Found ${TOTAL_SAMPLES} sample directories under ${ROOT_DIR}"

if [[ ${TOTAL_SAMPLES} -eq 0 ]]; then
  echo "[ERROR] No sample directories found under ${ROOT_DIR}. Exiting." >&2
  exit 1
fi

# -------------------------
# Loop over samples sequentially (adjust to parallelize if desired)
# -------------------------
for SAMPLE_DIR in "${SAMPLE_DIRS[@]}"; do
  SAMPLE_NAME=$(basename "${SAMPLE_DIR}")
  SAMPLE_OUTPUT_DIR="${OUTPUT_DIR}/${SAMPLE_NAME}"
  mkdir -p "${SAMPLE_OUTPUT_DIR}"

  # Find the first BAM file in the sample directory (adjust glob if you prefer a naming pattern)
  BAM_FILE=$(find "$SAMPLE_DIR" -maxdepth 2 -type f -name "*.bam" | sort | head -n 1 || true)

  if [[ -z "$BAM_FILE" ]]; then
    echo "[WARN] No BAM file found for sample ${SAMPLE_NAME} in ${SAMPLE_DIR}; skipping."
    continue
  fi

  echo "----------------------------------------"
  echo "[START] Sample: ${SAMPLE_NAME}  ($(date))"
  echo "  BAM: ${BAM_FILE}"
  echo "  Output dir: ${SAMPLE_OUTPUT_DIR}"

  # Ensure BAM index exists; create if missing
  if [[ ! -s "${BAM_FILE}.bai" && ! -s "${BAM_FILE%.bam}.bai" ]]; then
    echo "[INFO] BAM index missing for ${BAM_FILE}; creating with samtools index..."
    "${SAMTOOLS_BIN}" index "${BAM_FILE}"
  fi

  # Construct file paths for intermediate and final outputs
  TMP_BCF="${SAMPLE_OUTPUT_DIR}/${SAMPLE_NAME}_bcftools.bcf"
  OUTPUT_VCF_GZ="${SAMPLE_OUTPUT_DIR}/${SAMPLE_NAME}_bcftools.vcf.gz"

  # Remove old outputs if present (be careful — you may prefer to skip existing samples)
  rm -f "${TMP_BCF}" "${OUTPUT_VCF_GZ}" "${OUTPUT_VCF_GZ}.tbi"

  # Run mpileup -> call pipeline:
  #  - mpileup: produce uncompressed BCF (-Ou) for piping performance
  #  - call: call variants, output BCF (-Ob)
  #
  # Adjust MPILEUP_OPTS and CALL_OPTS above for mapping/base quality filters or depth caps:
  # e.g., MPILEUP_OPTS="-f ${REF_FASTA} -Ou -q 20 -Q 13 -d 2000"
  set -x
  ${BCFTOOLS_BIN} mpileup ${MPILEUP_OPTS} "${BAM_FILE}" \
    | ${BCFTOOLS_BIN} call ${CALL_OPTS} -o "${TMP_BCF}"
  set +x

  # Check the temporary BCF was created
  if [[ ! -s "${TMP_BCF}" ]]; then
    echo "[ERROR] BCF not produced for ${SAMPLE_NAME} (TMP_BCF missing). Skipping." >&2
    continue
  fi

  # Convert to compressed VCF (bgzip) and index with tabix (VCF spec)
  echo "[INFO] Converting BCF to gzipped VCF and indexing..."
  ${BCFTOOLS_BIN} view -Oz -o "${OUTPUT_VCF_GZ}" "${TMP_BCF}"
  "${TABIX_BIN}" -p vcf "${OUTPUT_VCF_GZ}"

  # Optionally remove the intermediate BCF to save space
  rm -f "${TMP_BCF}"

  echo "[DONE] Sample: ${SAMPLE_NAME}  ($(date))"
  echo "  Output VCF: ${OUTPUT_VCF_GZ}"
done

echo "========================================"
echo "[ALL DONE] Processed ${TOTAL_SAMPLES} samples."
echo "Outputs written under: ${OUTPUT_DIR}"
echo "Provenance file: ${VERSIONS_FILE}"
echo "========================================"
