#!/usr/bin/env bash
#
# run_irfinder_hardcoded_annotated.sh
#
# Heavily annotated, hard-coded script to run IRFinder on a set of BAMs and
# produce the per-sample IRFinder outputs expected by the downstream
# Δlogit(IR) pairwise meta-analysis pipeline.
#
# What this does (high level)
#  - Looks for BAM files in a hard-coded directory (one BAM per sample).
#  - Ensures each BAM is indexed (.bai); creates index with samtools if missing.
#  - Runs IRFinder in BAM mode for each sample using a pre-built IRFinder reference
#    directory (IRFINDER_REF_DIR). The IRFinder output for each sample is placed
#    under OUT_ROOT/samples/<sample_name>/ and will include IRFinder-IR-nondir.txt.
#  - Performs simple success/failure logging and produces a sample list file
#    pointing at the produced IRFinder outputs for use by downstream R scripts.
#
# IMPORTANT: this script is intentionally hard-coded. Edit the paths below to
# reflect your environment before running. The script assumes you have a working
# IRFinder reference already built (see earlier scripts that build the reference).
#
# Typical IRFinder usage (examples):
#   # Build a reference (done separately; not in this script)
#   IRFinder -m BuildRefProcess -r /path/to/IRFinder_ref genome.fa annotation.gtf
#
#   # Run IRFinder on one BAM (BAM input mode):
#   IRFinder -m BAM -r /path/to/IRFinder_ref -d outdir sample.bam
#
# Notes:
#  - IRFinder supports different invocation styles across versions. The form used
#    here ("-m BAM -r <ref> -d <outdir> <bam>") is common. If your IRFinder version
#    requires different flags, edit IRFINDER_CMD below accordingly.
#  - IRFinder requires reasonably recent R and some R packages internally; make
#    sure the IRFinder installation you call is functional in your environment.
#
set -euo pipefail
shopt -s nullglob

# ---------------------------
# HARD-CODED SETTINGS (EDIT THESE)
# ---------------------------
# Path to IRFinder executable (set to full path if not on PATH)
IRFINDER_BIN="/home/jgd/miniconda3/envs/irfinder_env/bin/IRFinder"

# Path to an existing IRFinder reference directory (produced by BuildRefProcess or BuildRefFromSTARRef)
IRFINDER_REF_DIR="/home/jgd/Documents/2026-04-01_proc-enu_bioinformatics_results/ref/irfinder"

# Directory containing input BAM files (one BAM per sample)
BAM_DIR="/home/jgd/Documents/2026-04-01_proc-enu_bioinformatics_results/alignments/bams"

# Output root (this will get a timestamped subdirectory for each run)
OUT_BASE="/home/jgd/PROJECT/irfinder_counts"   # parent directory
TS="$(date +%Y%m%d_%H%M%S)"
OUT_ROOT="${OUT_BASE}/run_${TS}"

# Number of parallel jobs for running IRFinder (if you have many cores)
# Set to 1 for sequential runs; set higher to run multiple samples in parallel (GNU parallel required).
JOBS=4

# Samtools binary (used to index BAMs if .bai missing)
SAMTOOLS_BIN="samtools"

# Where to write a simple manifest of produced sample outputs (used by downstream scripts)
MANIFEST="${OUT_ROOT}/samples_list.txt"

# Toggle: compress IRFinder outputs after run (keeps disk usage reasonable)
COMPRESS_OUTPUTS=true
# ---------------------------
# End of hard-coded settings
# ---------------------------

# ---------------------------
# Basic safety & environment checks
# ---------------------------
echo "[INFO] Starting IRFinder batch run"
echo "[INFO] OUT_ROOT = ${OUT_ROOT}"
mkdir -p "${OUT_ROOT}/samples"
mkdir -p "${OUT_ROOT}/logs"

# Check IRFinder executable
if [[ ! -x "${IRFINDER_BIN}" ]]; then
  echo "[ERROR] IRFinder binary not found or not executable at: ${IRFINDER_BIN}" >&2
  echo "[HINT] Edit IRFINDER_BIN at top of script to point to the IRFinder executable." >&2
  exit 1
fi

# Check IRFinder reference
if [[ ! -d "${IRFINDER_REF_DIR}" ]]; then
  echo "[ERROR] IRFinder reference directory not found: ${IRFINDER_REF_DIR}" >&2
  echo "[HINT] Build IRFinder reference first (see IRFinder BuildRefProcess / BuildRefFromSTARRef)." >&2
  exit 1
fi

# Check bam dir
if [[ ! -d "${BAM_DIR}" ]]; then
  echo "[ERROR] BAM directory not found: ${BAM_DIR}" >&2
  exit 1
fi

# Check samtools
if ! command -v "${SAMTOOLS_BIN}" >/dev/null 2>&1; then
  echo "[ERROR] samtools not found as '${SAMTOOLS_BIN}'. Please install samtools or change SAMTOOLS_BIN variable." >&2
  exit 1
fi

# Check parallel availability if JOBS>1
if [[ "${JOBS}" -gt 1 ]]; then
  if ! command -v parallel >/dev/null 2>&1; then
    echo "[WARN] GNU parallel not found; falling back to sequential execution (JOBS=1)" >&2
    JOBS=1
  fi
fi

echo "[INFO] IRFinder: ${IRFINDER_BIN}"
echo "[INFO] IRFinder ref: ${IRFINDER_REF_DIR}"
echo "[INFO] Using samtools: ${SAMTOOLS_BIN}"
echo "[INFO] JOBS: ${JOBS}"

# ---------------------------
# Build list of input BAMs
# ---------------------------
# We accept BAMs ending with .bam (case-sensitive). Adjust glob if you use .BAM etc.
mapfile -t BAM_FILES < <(find "${BAM_DIR}" -maxdepth 1 -type f -name '*.bam' -print | sort)

if [[ ${#BAM_FILES[@]} -eq 0 ]]; then
  echo "[ERROR] No BAM files found in ${BAM_DIR} (expected files ending in .bam)" >&2
  exit 1
fi

echo "[INFO] Found ${#BAM_FILES[@]} BAM files to process"
for b in "${BAM_FILES[@]}"; do
  echo "  - $(basename "${b}")"
done

# ---------------------------
# Helper: run IRFinder on one sample
# ---------------------------
run_one_sample() {
  local bam="$1"
  local sample
  sample="$(basename "${bam}" .bam)"   # sample name derived from BAM filename (strip .bam)
  local sample_out="${OUT_ROOT}/samples/${sample}"
  local logf="${OUT_ROOT}/logs/${sample}.irfinder.log"

  echo "[INFO] Processing sample: ${sample}" | tee -a "${logf}"

  # Ensure BAM is indexed; IRFinder can require index depending on version/mode.
  if [[ ! -s "${bam}.bai" && ! -s "${bam%.bam}.bai" ]]; then
    echo "[INFO] BAM index missing for ${bam}. Creating index with samtools index..." | tee -a "${logf}"
    "${SAMTOOLS_BIN}" index "${bam}" 2>> "${logf}" || { echo "[ERROR] samtools index failed for ${bam}" | tee -a "${logf}"; return 2; }
  fi

  # Remove sample output dir if present (we want a fresh run), then create it
  rm -rf "${sample_out}"
  mkdir -p "${sample_out}"

  # Compose IRFinder command.
  # Common IRFinder invocation (BAM mode):
  #   IRFinder -m BAM -r <IRFinder_ref_dir> -d <outdir> <input.bam>
  #
  # Note: IRFinder versions and flags may differ. If your version requires different flags,
  # edit the command below.
  cmd=( "${IRFINDER_BIN}" -m BAM -r "${IRFINDER_REF_DIR}" -d "${sample_out}" "${bam}" )

  echo "[INFO] Running IRFinder for ${sample} ..." | tee -a "${logf}"
  echo "[DEBUG] CMD: ${cmd[*]}" >> "${logf}"

  # Execute IRFinder and capture exit code
  if "${cmd[@]}" >> "${logf}" 2>&1; then
    echo "[INFO] IRFinder finished for ${sample}" | tee -a "${logf}"
  else
    local rc=$?
    echo "[ERROR] IRFinder failed for ${sample} with exit code ${rc}. See ${logf}" >&2
    return ${rc}
  fi

  # Validate expected output exists (IRFinder-IR-nondir.txt)
  local irfile="${sample_out}/IRFinder-IR-nondir.txt"
  if [[ ! -s "${irfile}" ]]; then
    echo "[ERROR] Expected IRFinder output not found for ${sample}: ${irfile}" | tee -a "${logf}"
    return 3
  fi

  # Optional: compress the sample output directory (to save space)
  if [[ "${COMPRESS_OUTPUTS}" = true ]]; then
    echo "[INFO] Compressing IRFinder output: ${sample_out}" | tee -a "${logf}"
    # Create a tar.gz of the sample folder (keeps original for safety; remove original if desired)
    tar -C "${OUT_ROOT}/samples" -czf "${OUT_ROOT}/samples/${sample}.irfinder.tar.gz" "${sample}" >> "${logf}" 2>&1 \
      || echo "[WARN] Compression failed for ${sample}" | tee -a "${logf}"
    # Optionally remove the uncompressed directory after tar succeeded:
    # rm -rf "${sample_out}"
    # But we leave the directory intact because downstream scripts expect the IRFinder-IR-nondir.txt under samples/<sample>
  fi

  return 0
}

# ---------------------------
# Run IRFinder for all samples (parallel or sequential)
# ---------------------------
if [[ "${JOBS}" -gt 1 ]]; then
  echo "[INFO] Running IRFinder in parallel with ${JOBS} jobs"
  # Export necessary variables and functions for GNU parallel
  export IRFINDER_BIN IRFINDER_REF_DIR OUT_ROOT SAMTOOLS_BIN COMPRESS_OUTPUTS
  export -f run_one_sample
  printf '%s\n' "${BAM_FILES[@]}" | parallel -j "${JOBS}" --halt now,fail=1 run_one_sample {}
  rc=$?
  if [[ $rc -ne 0 ]]; then
    echo "[ERROR] Parallel IRFinder processing failed (exit ${rc})" >&2
    exit $rc
  fi
else
  echo "[INFO] Running IRFinder sequentially"
  for bam in "${BAM_FILES[@]}"; do
    if ! run_one_sample "${bam}"; then
      echo "[ERROR] IRFinder failed for BAM: ${bam}" >&2
      exit 1
    fi
  done
fi

# ---------------------------
# Build a manifest of produced sample folders (for downstream scripts)
# ---------------------------
echo "[INFO] Building manifest of sample outputs to ${MANIFEST}"
: > "${MANIFEST}"
for sample_dir in "${OUT_ROOT}/samples"/*; do
  # skip tar.gz files (if compressing)
  if [[ -d "${sample_dir}" ]]; then
    sample="$(basename "${sample_dir}")"
    irfile="${sample_dir}/IRFinder-IR-nondir.txt"
    if [[ -s "${irfile}" ]]; then
      echo -e "${sample}\t${irfile}" >> "${MANIFEST}"
    else
      echo -e "${sample}\tMISSING" >> "${MANIFEST}"
    fi
  fi
done

echo "[INFO] Done. Outputs written under: ${OUT_ROOT}"
echo "[INFO] Sample manifest: ${MANIFEST}"
echo "[INFO] Sample directories (first 20):"
ls -1 "${OUT_ROOT}/samples" | head -n 20

# Final message with guidance for downstream usage
cat <<EOF

SUCCESS: IRFinder run complete (or at least attempted) for ${#BAM_FILES[@]} samples.

Downstream:
 - The Δlogit pairwise R script expects sample folders at:
     ${OUT_ROOT}/samples/<sample>/IRFinder-IR-nondir.txt

 - You can copy or point your pairwise script to:
     SAMPLES_DIR = ${OUT_ROOT}/samples

 - Inspect logs for any failed samples:
     ls -1 ${OUT_ROOT}/logs/*.irfinder.log
     tail -n 200 ${OUT_ROOT}/logs/<sample>.irfinder.log

Caveats / troubleshooting:
 - If IRFinder fails for a sample, common causes include:
    * BAM not coordinate-sorted (IRFinder needs properly sorted BAM)
    * Missing/mismatched FASTA/GTF in the IRFinder reference (contig name mismatch)
    * Insufficient read coverage in the sample for IRFinder heuristics
    * Version mismatch between IRFinder and the reference build method (try BuildRefProcess fallback)
 - If you need me to adapt this script to:
    * Run on a compute cluster (Slurm) or dispatch per-sample jobs,
    * Use different IRFinder invocation flags for your version,
    * Pre-filter BAMs (e.g., remove duplicates) before IRFinder,
   tell me and I will produce the modified script.

EOF

exit 0
