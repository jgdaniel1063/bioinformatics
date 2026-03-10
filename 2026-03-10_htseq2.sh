#!/usr/bin/env bash
#
# htseq2_run_hardcoded_annotated.sh
#
# Heavily annotated, hard-coded script to run HTSeq-count over a directory of BAM files,
# produce per-sample count files, and assemble a merged counts matrix suitable for downstream
# DE analysis (DESeq2/edgeR, etc.).
#
# Design goals & notes (read this before running)
# - Hard-coded paths and parameters live near the top of the script — edit them before running.
# - The script attempts to be safe and informative:
#     * Checks for required tools (htseq-count, samtools, bgzip/tabix optional)
#     * Verifies the annotation GTF is present
#     * Handles both single-end and paired-end BAMs. For paired-end, HTSeq-count requires
#       reads to be name-sorted (use -r name). For single-end you can use position-sorted BAMs.
#     * Supports GNU parallel for concurrent runs; falls back to xargs if missing.
# - The script writes a timestamped run directory (RUN_DIR) under OUT_BASE and creates:
#     RUN_DIR/counts/<sample>.counts.txt    - HTSeq raw output per sample
#     RUN_DIR/counts_matrix.tsv              - merged matrix (genes x samples)
#     RUN_DIR/logs/...                       - per-sample logs + run.log
#     RUN_DIR/failed_samples.txt             - list of samples that failed
#
# Performance & I/O
# - Sorting BAMs (especially name-sorting) can be very I/O and CPU heavy. If your aligner
#   produced name-sorted BAMs already, point BAM_DIR at those to save time.
# - For large cohorts, consider running per-sample work on a cluster (this script is for local).
#
# HTSeq-count invocation notes
# - Common recommended invocation for paired-end RNA-seq:
#     samtools sort -n in.bam -o in.namesorted.bam
#     htseq-count -f bam -r name -s <0/1/2> -t exon -i gene_id in.namesorted.bam > out.counts
#   - -r name: iterate by read name so paired reads get processed together
#   - -s (stranded): 0=unstranded, 1=stranded, 2=reverse
#   - -t (feature type): usually "exon"
#   - -i (attribute id): usually "gene_id"
#
# - For single-end or when using position-sorted BAMs use -r pos.
#
# Requirements
# - htseq-count (installed as part of HTSeq/pip/bioconda)
# - samtools (for sorting/indexing)
# - GNU parallel (optional, speeds up batch processing)
#
# -----------------------------------------------------------------------------
set -euo pipefail
shopt -s nullglob

# ==============================
# USER-EDITABLE HARD-CODED VALUES
# Edit these to match your environment BEFORE running the script.
# ==============================

# Directory containing input BAM files (one BAM per sample).
# The script will search for *.bam files (case-sensitive).
BAM_DIR="/path/to/bam_dir"

# Path to GTF annotation file (required). Prefer coordinate-sorted GTF that matches the FASTA used for alignment.
GTF_FILE="/path/to/annotation.gtf"

# Output base directory. A timestamped RUN_DIR will be created under this.
OUT_BASE="/path/to/htseq_runs"

# HTSeq-count and samtools commands (set full paths here if not on PATH)
HTSEQ_CMD="${HTSEQ_CMD:-htseq-count}"   # e.g., /home/user/miniconda3/bin/htseq-count
SAMTOOLS_CMD="${SAMTOOLS_CMD:-samtools}" # e.g., /home/user/miniconda3/bin/samtools

# Analysis parameters — set to appropriate values for your library prep
PAIRED=1               # 1 = paired-end library (use name-sorted BAM + -r name), 0 = single-end (use -r pos)
STRANDED=0             # 0 = unstranded, 1 = stranded, 2 = reverse-stranded
FEATURETYPE="exon"     # feature type in GTF to count (commonly 'exon')
ATTR="gene_id"         # attribute in GTF to use as gene identifier (commonly 'gene_id' or 'gene_name')
MODE="union"           # htseq -m mode: union (default), intersection-strict, intersection-nonempty

# Parallelism: number of concurrent htseq-count processes. If GNU parallel is available we'll use it.
JOBS=4

# Whether to keep intermediate name-sorted BAMs created by this script (set to 0 to delete them)
KEEP_NAME_SORTED_BAMS=0

# How to name per-sample output files
RUN_SUFFIX="htseq"     # used in output filenames: <sample>.<RUN_SUFFIX>.counts.txt

# ==============================
# END OF USER-EDITABLE VALUES
# ==============================

# ------------------------------
# Derived paths and run setup
# ------------------------------
TS="$(date +%Y%m%d_%H%M%S)"
RUN_DIR="${OUT_BASE}/run_${TS}"
COUNTS_DIR="${RUN_DIR}/counts"
LOG_DIR="${RUN_DIR}/logs"
TMP_DIR="${RUN_DIR}/tmp"

mkdir -p "${COUNTS_DIR}" "${LOG_DIR}" "${TMP_DIR}" || { echo "Failed to create run directories"; exit 1; }

RUN_LOG="${LOG_DIR}/run.log"
: > "${RUN_LOG}"   # create/clear

# Helper to log to run.log with timestamp
log() {
  local msg ts
  ts="$(date '+%F %T')"
  msg="[$ts] $*"
  echo "${msg}" | tee -a "${RUN_LOG}" >&2
}

log "START HTSeq-count batch run"
log "BAM_DIR=${BAM_DIR}"
log "GTF_FILE=${GTF_FILE}"
log "OUT_BASE=${OUT_BASE}"
log "RUN_DIR=${RUN_DIR}"
log "PAIRED=${PAIRED} STRANDED=${STRANDED} FEATURETYPE=${FEATURETYPE} ATTR=${ATTR} MODE=${MODE}"
log "JOBS=${JOBS} HTSEQ_CMD=${HTSEQ_CMD} SAMTOOLS_CMD=${SAMTOOLS_CMD}"

# ------------------------------
# Tool checks
# ------------------------------
command -v "${HTSEQ_CMD}" >/dev/null 2>&1 || { echo "[ERROR] htseq-count not found: ${HTSEQ_CMD}" | tee -a "${RUN_LOG}"; exit 1; }
command -v "${SAMTOOLS_CMD}" >/dev/null 2>&1 || { echo "[ERROR] samtools not found: ${SAMTOOLS_CMD}" | tee -a "${RUN_LOG}"; exit 1; }

# Optionally check GNU parallel
PARALLEL_CMD=""
if command -v parallel >/dev/null 2>&1; then
  PARALLEL_CMD="parallel"
  log "GNU parallel found; will use up to ${JOBS} concurrent jobs"
else
  log "GNU parallel not found; will fall back to xargs (limited concurrency)"
fi

# Validate GTF
if [[ ! -s "${GTF_FILE}" ]]; then
  log "❌ GTF file not found: ${GTF_FILE}"
  exit 1
fi

# ------------------------------
# Discover BAM files
# ------------------------------
mapfile -t BAM_FILES < <(find "${BAM_DIR}" -maxdepth 1 -type f -name "*.bam" | sort)

if [[ ${#BAM_FILES[@]} -eq 0 ]]; then
  log "❌ No BAM files found in ${BAM_DIR} (expected *.bam). Exiting."
  exit 1
fi

log "Found ${#BAM_FILES[@]} BAM files (showing up to 40):"
for ((i=0;i<${#BAM_FILES[@]} && i<40;i++)); do
  log "  ${BAM_FILES[i]}"
done
if [[ ${#BAM_FILES[@]} -gt 40 ]]; then
  log "  ... and $(( ${#BAM_FILES[@]} - 40)) more"
fi

# ------------------------------
# Prepare per-sample worker script (helper)
# ------------------------------
HELPER="${TMP_DIR}/htseq_worker.sh"
cat > "${HELPER}" <<'SH'
#!/usr/bin/env bash
set -euo pipefail
sample_bam="$1"
sample_name="$2"
counts_out="$3"
log_out="$4"
# User variables substituted by parent script: HTSEQ_CMD, SAMTOOLS_CMD, PAIRED, STRANDED, FEATURETYPE, ATTR, MODE, TMP_DIR, KEEP_NAME_SORTED_BAMS
# The parent will export those variables into the environment for this helper to use.

{
  echo "[$(date '+%F %T')] START sample=${sample_name} bam=${sample_bam}"
  # Ensure BAM exists and is readable
  if [[ ! -s "${sample_bam}" ]]; then
    echo "ERROR: BAM not found: ${sample_bam}"
    exit 2
  fi

  # For paired-end libraries HTSeq-count expects name-sorted BAM (if using -r name)
  # For single-end or when using -r pos, position-sorted BAM is fine.
  if [[ "${PAIRED}" -eq 1 ]]; then
    sort_mode="name"
    sort_flag="-n"
    # Create a name-sorted BAM in TMP if the BAM filename does not look name-sorted
    # Heuristic: if filename contains ".namesorted." assume it's already name-sorted
    if [[ "${sample_bam}" == *".namesorted."* || "${sample_bam}" == *".namesorted.bam" ]]; then
      sorted_bam="${sample_bam}"
      made_temp_sorted=0
    else
      sorted_bam="${TMP_DIR}/${sample_name}.namesorted.bam"
      echo "  Name-sorting BAM to ${sorted_bam} (this may take time) ..."
      ${SAMTOOLS_CMD} sort -n -o "${sorted_bam}" "${sample_bam}" || { echo "ERROR: samtools sort -n failed"; exit 3; }
      made_temp_sorted=1
    fi
    htseq_r_opt="-r name"
  else
    sort_mode="pos"
    # Assume coordinate-sorted BAM is OK; if not sorted you may want to sort with samtools sort (no -n)
    sorted_bam="${sample_bam}"
    made_temp_sorted=0
    htseq_r_opt="-r pos"
  fi

  # Build htseq-count command
  cmd="${HTSEQ_CMD} -f bam ${htseq_r_opt} -s ${STRANDED} -t ${FEATURETYPE} -i ${ATTR} -m ${MODE} \"${sorted_bam}\""
  echo "  Running: ${cmd}"
  # Execute and write to counts_out
  # Note: redirect stdout to counts_out, stderr to log_out
  eval ${cmd} > "${counts_out}" 2> "${log_out}.htseq.stderr.log" || rc=$? || true
  rc=${rc:-0}
  if [[ ${rc} -ne 0 ]]; then
    echo "ERROR: htseq-count failed with exit ${rc} (see ${log_out}.htseq.stderr.log)"
    # Clean up if temp sorted created
    if [[ ${made_temp_sorted} -eq 1 && "${KEEP_NAME_SORTED_BAMS}" -eq 0 ]]; then
      rm -f "${sorted_bam}"
    fi
    exit ${rc}
  fi

  # Optionally keep or remove temporary name-sorted BAM
  if [[ ${made_temp_sorted} -eq 1 ]]; then
    if [[ "${KEEP_NAME_SORTED_BAMS}" -eq 1 ]]; then
      echo "  Kept temporary namesorted BAM: ${sorted_bam}"
    else
      rm -f "${sorted_bam}"
    fi
  fi

  echo "[$(date '+%F %T')] DONE sample=${sample_name} rc=0"
} >> "${log_out}" 2>&1
SH

chmod +x "${HELPER}"

# Export env vars used by the helper script
export HTSEQ_CMD SAMTOOLS_CMD PAIRED STRANDED FEATURETYPE ATTR MODE TMP_DIR KEEP_NAME_SORTED_BAMS

# ------------------------------
# Build list of sample jobs
# ------------------------------
JOBLIST="${TMP_DIR}/htseq_jobs.tsv"
: > "${JOBLIST}"
FAILED_SAMPLES="${RUN_DIR}/failed_samples.txt"
: > "${FAILED_SAMPLES}"

for bam in "${BAM_FILES[@]}"; do
  sample=$(basename "${bam}")
  sample="${sample%.bam}"
  counts_out="${COUNTS_DIR}/${sample}.${RUN_SUFFIX}.counts.txt"
  log_out="${LOG_DIR}/${sample}.log"
  echo -e "${bam}\t${sample}\t${counts_out}\t${log_out}" >> "${JOBLIST}"
done

# ------------------------------
# Run per-sample jobs (parallel or sequential)
# ------------------------------
log "Launching per-sample HTSeq-count jobs (JOBS=${JOBS})"
if [[ -n "${PARALLEL_CMD}" ]]; then
  # Use GNU parallel -0 style not needed here; pass TSV columns
  "${PARALLEL_CMD}" -j "${JOBS}" --colsep '\t' bash "${HELPER}" {1} {2} {3} {4} :::: "${JOBLIST}" || true
else
  # fallback: sequential with while/read
  while IFS=$'\t' read -r bam sample counts_out log_out; do
    bash "${HELPER}" "${bam}" "${sample}" "${counts_out}" "${log_out}" || {
      echo "${sample}" >> "${FAILED_SAMPLES}"
      log "Sample ${sample} failed, recorded in ${FAILED_SAMPLES}"
    }
  done < "${JOBLIST}"
fi

# Note: GNU parallel + helper returns non-zero if any job failed; we captured that with || true and record failed samples separately.
# Collect failed samples by scanning logs for 'ERROR' or using FAILED_SAMPLES file
if [[ ! -s "${FAILED_SAMPLES}" ]]; then
  # attempt to discover failures by looking for non-zero exit logs
  while IFS=$'\t' read -r bam sample counts_out log_out; do
    if [[ ! -s "${counts_out}" ]]; then
      echo "${sample}" >> "${FAILED_SAMPLES}"
    fi
  done < "${JOBLIST}"
fi

NUM_FAILED=$(wc -l < "${FAILED_SAMPLES}" | tr -d ' ')
if [[ "${NUM_FAILED}" -gt 0 ]]; then
  log "WARNING: ${NUM_FAILED} sample(s) failed. See ${FAILED_SAMPLES} and per-sample logs under ${LOG_DIR}"
else
  log "All samples produced counts successfully (no failures detected)"
fi

# ------------------------------
# Merge per-sample count files into a single matrix
# HTSeq outputs have consistent gene ordering (provided same GTF); we assume that here.
# If gene orders differ, you should use a join-based merge keyed on gene IDs.
# ------------------------------
log "Merging per-sample count files into matrix"

# Build list of produced count files in the order of JOBLIST
mapfile -t COUNT_FILES < <(awk -F'\t' '{print $3}' "${JOBLIST}")

# Verify at least one count file exists
if [[ ${#COUNT_FILES[@]} -eq 0 ]]; then
  log "❌ No count files found to merge. Exiting."
  exit 1
fi

# Use the first count file to get the gene order
first="${COUNT_FILES[0]}"
if [[ ! -s "${first}" ]]; then
  log "❌ First count file missing: ${first}"
  exit 1
fi

# Extract gene names (1st column) from first file
cut -f1 "${first}" > "${TMP_DIR}/genes.list"

# For each count file, extract second column and paste together
paste_cmd=(paste)
paste_cmd+=("${TMP_DIR}/genes.list")
header="gene"
samples=()
for cf in "${COUNT_FILES[@]}"; do
  # Check file exists
  if [[ ! -s "${cf}" ]]; then
    log "WARNING: count file missing, filling with zeros: ${cf}"
    # create zero column matching genes.list
    awk '{print 0}' "${TMP_DIR}/genes.list" > "${TMP_DIR}/$(basename "${cf}").zeros"
    paste_cmd+=("${TMP_DIR}/$(basename "${cf}").zeros")
    samples+=("$(basename "${cf}")")
  else
    # htseq-count outputs format: <gene> <count>
    paste_cmd+=(<(cut -f2 "${cf}"))
    samples+=("$(basename "${cf}" | sed -e "s/\.${RUN_SUFFIX}\.counts\.txt$//" -e 's/\.counts\.txt$//'))
  fi
done

# Compose header line (gene + sample names)
header_line="gene"
for s in "${samples[@]}"; do
  header_line="${header_line}\t${s}"
done

# Use awk to join pasted columns (we used process substitution so call paste via an array)
# But since paste_cmd contains process substitutions which are expanded only in bash, use eval to run it.
# Build the paste command string:
paste_str="paste ${TMP_DIR}/genes.list"
for cf in "${COUNT_FILES[@]}"; do
  if [[ -s "${cf}" ]]; then
    paste_str+=" <(cut -f2 '${cf}')"
  else
    paste_str+=" ${TMP_DIR}/$(basename "${cf}").zeros"
  fi
done

# Execute paste via bash -c to allow process substitution
set +e
bash -c "${paste_str} > '${TMP_DIR}/counts_merged_body.tsv'"
rc_paste=$?
set -e
if [[ ${rc_paste} -ne 0 ]]; then
  log "ERROR: merging columns with paste failed (rc=${rc_paste}). Falling back to safe join-based merge."

  # Fallback robust join: iterative join on gene column
  cp "${first}" "${TMP_DIR}/merge_0.tsv"
  awk '{print $1 "\t" $2}' "${first}" > "${TMP_DIR}/merge_0.tsv"
  idx=1
  for cf in "${COUNT_FILES[@]:1}"; do
    awk '{print $1 "\t" $2}' "${cf}" > "${TMP_DIR}/col_${idx}.tsv"
    join -t $'\t' -a1 -a2 -e "0" -o auto -t $'\t' "${TMP_DIR}/merge_0.tsv" "${TMP_DIR}/col_${idx}.tsv" > "${TMP_DIR}/merge_1.tsv"
    mv "${TMP_DIR}/merge_1.tsv" "${TMP_DIR}/merge_0.tsv"
    idx=$((idx+1))
  done
  mv "${TMP_DIR}/merge_0.tsv" "${TMP_DIR}/counts_merged_body.tsv"
fi

# Write header + body to final counts matrix
echo -e "${header_line}" > "${RUN_DIR}/counts_matrix.tsv"
cat "${TMP_DIR}/counts_merged_body.tsv" >> "${RUN_DIR}/counts_matrix.tsv"
log "Counts matrix written to: ${RUN_DIR}/counts_matrix.tsv"

# Also write a CSV for convenience
awk 'BEGIN{FS=OFS="\t"} NR==1{print; next} {print}' "${RUN_DIR}/counts_matrix.tsv" > "${RUN_DIR}/counts_matrix.csv"

# ------------------------------
# Basic QC / summary
# ------------------------------
# Compute library sizes (sum of counts per sample)
awk -v OFS="\t" 'NR==1{for(i=2;i<=NF;i++) header[i]=$i; next} {for(i=2;i<=NF;i++) sums[i]+=($i+0)} END{print "sample","libSize"; for(i=2;i in header;i++) print header[i],sums[i]}' "${RUN_DIR}/counts_matrix.tsv" > "${RUN_DIR}/library_sizes.tsv"

log "Library sizes written to: ${RUN_DIR}/library_sizes.tsv"

# ------------------------------
# Done
# ------------------------------
log "HTSeq-count batch run COMPLETE"
log "RUN_DIR: ${RUN_DIR}"
log "Per-sample counts: ${COUNTS_DIR}"
log "Merged counts matrix: ${RUN_DIR}/counts_matrix.tsv"
log "Failed samples (if any): ${FAILED_SAMPLES}"

# Provide next-step suggestions
cat <<EOF > "${RUN_DIR}/README_next_steps.txt"
HTSeq2 run summary (generated at ${TS}):

- Per-sample HTSeq counts are in:
    ${COUNTS_DIR}/<sample>.${RUN_SUFFIX}.counts.txt

- Merged counts matrix:
    ${RUN_DIR}/counts_matrix.tsv
    (First column: gene id; subsequent columns: sample counts)

- Logs (per-sample and run-level):
    ${LOG_DIR}

Notes and recommendations:
- Inspect ${RUN_DIR}/library_sizes.tsv to verify library sizes and detect outliers.
- If you performed paired-end counting, ensure HTSeq was run with '-r name' on name-sorted BAMs.
- Before DE analysis:
    * Remove HTSeq internal summary rows (lines starting with '__') if you don't want them.
    * Use DESeq2/edgeR on the merged matrix; do NOT use HTSeq per-sample raw output where sample names don't match.
- If gene order mismatches appear, re-run this script and set COUNTS_MATRIX or adjust merging strategy (join-based) to be robust.

EOF

log "Wrote next-step notes to: ${RUN_DIR}/README_next_steps.txt"

exit 0
