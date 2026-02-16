#!/usr/bin/env bash
set -uo pipefail  # Removed -e to allow manual error handling for GATK

# --- Config (adjust as needed) ---
ROOT_DIR="/media/jgd/jd-1tb/bioinformatics/2024-04-01_proc-enu_paper/bam/star"
REF_FASTA="/media/jgd/jd-1tb/bioinformatics/reference/danRer11.fa"
REF_DICT="/media/jgd/jd-1tb/bioinformatics/reference/danRer11.dict"
DBSNP="/media/jgd/jd-1tb/bioinformatics/reference/danRer11_115_knownsites.raw.vcf"
KNOWN_INDELS="/media/jgd/jd-1tb/bioinformatics/reference/shavit.hq.indel.norm.cleaned.vcf.gz"
OUTPUT_DIR="/media/jgd/jd-1tb/bioinformatics/2024-04-01_proc-enu_paper/output/gatk"
THREADS=12  # Number of threads to use (adjust based on your CPU)
MEMORY=40   # Max memory in GB (adjust based on your RAM)
# -------------------------------

JOB_TITLE="${SLURM_JOB_NAME:-${JOB_NAME:-local}}"
TIMESTAMP_GLOBAL="$(date +%Y%m%d_%H%M%S)"

# Create output dir and global log file early for reliable logging
mkdir -p "$OUTPUT_DIR"
GLOBAL_LOG="${OUTPUT_DIR}/${JOB_TITLE}_${TIMESTAMP_GLOBAL}.log"

# Find all BAMs recursively within ROOT_DIR (any subdir depth)
mapfile -d '' BAM_LIST < <(find "$ROOT_DIR" -type f -name "*.bam" -print0)
TOTAL_SAMPLES=${#BAM_LIST[@]}

if [[ $TOTAL_SAMPLES -eq 0 ]]; then
  log "ERROR: No BAM files found under $ROOT_DIR"
  exit 1
fi

log() {
  local msg="$*"
  local ts
  ts="$(date '+%Y-%m-%d %H:%M:%S')"
  # Always append to global log file
  printf '%s [%s] %s\n' "$ts" "$JOB_TITLE" "$msg" >> "$GLOBAL_LOG"
  # Append to sample log file if set
  if [[ -n "${SAMPLE_LOG:-}" ]]; then
    printf '%s [%s] %s\n' "$ts" "$JOB_TITLE" "$msg" >> "$SAMPLE_LOG"
  fi
  # Always print to stdout
  printf '%s [%s] %s\n' "$ts" "$JOB_TITLE" "$msg"
}

# Trap for unhandled errors (now that -e is removed)
trap 'log "ERROR: Unhandled failure at line $LINENO"; exit 1' ERR

process_bam() {
  local BAM_FILE="$1"
  local BAM_PATH_REL
  BAM_PATH_REL=$(realpath --relative-to="$ROOT_DIR" "$BAM_FILE")
  local SAMPLE_BASENAME
  SAMPLE_BASENAME="$(basename "$BAM_FILE" .bam)"
  local SAMPLE_DIR
  SAMPLE_DIR="$(dirname "$BAM_FILE")"

  # Unique sample identifier for output folder and renaming
  local SAMPLE_NAME_REL="${BAM_PATH_REL//\//_}"
  SAMPLE_NAME_REL="${SAMPLE_NAME_REL%.bam}"

  # -- Check for existing output directory to enable autoresume --
  local SAMPLE_OUTPUT_DIR
  local SAMPLE_TIMESTAMP
  local existing_dirs
  mapfile -d '' existing_dirs < <(find "$OUTPUT_DIR" -maxdepth 1 -type d -name "${SAMPLE_NAME_REL}_[0-9]*" -print0 2>/dev/null | sort -z)
  if [[ ${#existing_dirs[@]} -gt 0 ]]; then
    # Use the latest existing directory (sorted by name, which includes timestamp)
    SAMPLE_OUTPUT_DIR="${existing_dirs[-1]}"
    # Extract timestamp from directory name (e.g., from "sample_20260214_120000" get "20260214_120000")
    SAMPLE_TIMESTAMP="${SAMPLE_OUTPUT_DIR##*_}"
    log "Resuming from existing output dir: $SAMPLE_OUTPUT_DIR"
  else
    # No existing dir; create new one
    SAMPLE_TIMESTAMP="$(date +%Y%m%d_%H%M%S)"
    SAMPLE_OUTPUT_DIR="${OUTPUT_DIR}/${SAMPLE_NAME_REL}_${SAMPLE_TIMESTAMP}"
    mkdir -p "$SAMPLE_OUTPUT_DIR"
    log "Starting new run; created output dir: $SAMPLE_OUTPUT_DIR"
  fi

  find "$SAMPLE_OUTPUT_DIR" -type f -size 0 -delete

  # Rename BAM if it does not already match its containing folder
  local EXPECTED_BAM="${SAMPLE_DIR}/$(basename "$SAMPLE_DIR").bam"
  if [[ "$BAM_FILE" != "$EXPECTED_BAM" ]]; then
    log "Renaming $BAM_FILE -> $EXPECTED_BAM"
    mv "$BAM_FILE" "$EXPECTED_BAM"
    BAM_FILE="$EXPECTED_BAM"
  fi

  # Set per-sample log file (append to existing if resuming)
  SAMPLE_LOG="${SAMPLE_OUTPUT_DIR}/${SAMPLE_NAME_REL}_${SAMPLE_TIMESTAMP}_${JOB_TITLE}.log"

  log "Starting processing of sample: $SAMPLE_NAME_REL"
  log "BAM: $BAM_FILE"
  log "Output dir: $SAMPLE_OUTPUT_DIR"
  log "REF_FASTA: $REF_FASTA"
  log "REF_DICT: $REF_DICT"

  # -- AddOrReplaceReadGroups --
  local RG_BAM="${SAMPLE_OUTPUT_DIR}/${SAMPLE_NAME_REL}_rg.bam"
  if [[ ! -s "$RG_BAM" ]] || ! samtools quickcheck "$RG_BAM"; then
    local TMP_RG_BAM="${RG_BAM}.tmp"
    log "Running AddOrReplaceReadGroups -> $(basename "$RG_BAM")"
    gatk --java-options "-Xmx${MEMORY}g" AddOrReplaceReadGroups \
      -I "$BAM_FILE" \
      -O "$TMP_RG_BAM" \
      -RGID 1 -RGLB library -RGPL ILLUMINA -RGPU unit -RGSM "$SAMPLE_NAME_REL" \
      --VALIDATION_STRINGENCY LENIENT &>> "$SAMPLE_LOG"
    if [[ $? -eq 0 && -s "$TMP_RG_BAM" ]] && samtools quickcheck "$TMP_RG_BAM"; then
      mv "$TMP_RG_BAM" "$RG_BAM"
      log "Indexing $(basename "$RG_BAM")"
      samtools index "$RG_BAM"  # Index the BAM
    else
      rm -f "$TMP_RG_BAM"
      log "AddOrReplaceReadGroups failed for $SAMPLE_NAME_REL"
      exit 1
    fi
  else
    log "Skipping AddOrReplaceReadGroups; found valid $RG_BAM"
    if [[ ! -f "${RG_BAM}.bai" ]]; then
      log "Indexing $(basename "$RG_BAM") (existing)"
      samtools index "$RG_BAM"  # Ensure existing BAM is indexed
    fi
  fi

  # -- MarkDuplicates --
  local MARKED_BAM="${SAMPLE_OUTPUT_DIR}/${SAMPLE_NAME_REL}_marked.bam"
  local METRICS="${SAMPLE_OUTPUT_DIR}/${SAMPLE_NAME_REL}_dup_metrics.txt"
  if [[ ! -s "$MARKED_BAM" ]] || ! samtools quickcheck "$MARKED_BAM"; then
    local TMP_MARKED_BAM="${MARKED_BAM}.tmp"
    log "Running MarkDuplicates -> $(basename "$MARKED_BAM")"
    gatk --java-options "-Xmx${MEMORY}g" MarkDuplicates \
      -I "$RG_BAM" \
      -O "$TMP_MARKED_BAM" \
      -M "$METRICS" \
      --VALIDATION_STRINGENCY LENIENT &>> "$SAMPLE_LOG"
    if [[ $? -eq 0 && -s "$TMP_MARKED_BAM" ]] && samtools quickcheck "$TMP_MARKED_BAM"; then
      mv "$TMP_MARKED_BAM" "$MARKED_BAM"
      log "Indexing $(basename "$MARKED_BAM")"
      samtools index "$MARKED_BAM"  # Index the BAM
    else
      rm -f "$TMP_MARKED_BAM"
      log "MarkDuplicates failed for $SAMPLE_NAME_REL"
      exit 1
    fi
  else
    log "Skipping MarkDuplicates; found valid $MARKED_BAM"
    if [[ ! -f "${MARKED_BAM}.bai" ]]; then
      log "Indexing $(basename "$MARKED_BAM") (existing)"
      samtools index "$MARKED_BAM"  # Ensure existing BAM is indexed
    fi
  fi

  # -- SplitNCigarReads --
  local SPLIT_BAM="${SAMPLE_OUTPUT_DIR}/${SAMPLE_NAME_REL}_split.bam"
  if [[ ! -s "$SPLIT_BAM" ]] || ! samtools quickcheck "$SPLIT_BAM"; then
    local TMP_SPLIT_BAM="${SPLIT_BAM}.tmp"
    log "Running SplitNCigarReads -> $(basename "$SPLIT_BAM")"
    gatk --java-options "-Xmx${MEMORY}g" SplitNCigarReads \
      -R "$REF_FASTA" \
      -I "$MARKED_BAM" \
      -O "$TMP_SPLIT_BAM" \
      --read-validation-stringency LENIENT &>> "$SAMPLE_LOG"
    if [[ $? -eq 0 && -s "$TMP_SPLIT_BAM" ]] && samtools quickcheck "$TMP_SPLIT_BAM"; then
      mv "$TMP_SPLIT_BAM" "$SPLIT_BAM"
      log "Indexing $(basename "$SPLIT_BAM")"
      samtools index "$SPLIT_BAM"  # Index the BAM
    else
      rm -f "$TMP_SPLIT_BAM"
      log "SplitNCigarReads failed for $SAMPLE_NAME_REL"
      exit 1
    fi
  else
    log "Skipping SplitNCigarReads; found valid $SPLIT_BAM"
    if [[ ! -f "${SPLIT_BAM}.bai" ]]; then
      log "Indexing $(basename "$SPLIT_BAM") (existing)"
      samtools index "$SPLIT_BAM"  # Ensure existing BAM is indexed
    fi
  fi

  # -- BQSR --
  local RECAL_TABLE="${SAMPLE_OUTPUT_DIR}/${SAMPLE_NAME_REL}_recal.table"
  local RECAL_BAM="${SAMPLE_OUTPUT_DIR}/${SAMPLE_NAME_REL}_recal.bam"
  local INPUT_BAM
  if [[ -f "$DBSNP" && -f "$KNOWN_INDELS" ]]; then
    if [[ ! -s "$RECAL_BAM" ]] || ! samtools quickcheck "$RECAL_BAM"; then
      local TMP_RECAL_BAM="${RECAL_BAM}.tmp"
      log "Running BaseRecalibrator -> $(basename "$RECAL_TABLE")"
      gatk --java-options "-Xmx${MEMORY}g" BaseRecalibrator \
        -R "$REF_FASTA" \
        -I "$SPLIT_BAM" \
        --known-sites "$DBSNP" \
        --known-sites "$KNOWN_INDELS" \
        -O "$RECAL_TABLE" \
        --read-validation-stringency LENIENT &>> "$SAMPLE_LOG"
      if [[ $? -ne 0 ]]; then
        log "BaseRecalibrator failed for $SAMPLE_NAME_REL"
        exit 1
      fi
      log "Running ApplyBQSR -> $(basename "$RECAL_BAM")"
      gatk --java-options "-Xmx${MEMORY}g" ApplyBQSR \
        -R "$REF_FASTA" \
        -I "$SPLIT_BAM" \
        --bqsr-recal-file "$RECAL_TABLE" \
        -O "$TMP_RECAL_BAM" \
        --read-validation-stringency LENIENT &>> "$SAMPLE_LOG"
      if [[ $? -eq 0 && -s "$TMP_RECAL_BAM" ]] && samtools quickcheck "$TMP_RECAL_BAM"; then
        mv "$TMP_RECAL_BAM" "$RECAL_BAM"
        log "Indexing $(basename "$RECAL_BAM")"
        samtools index "$RECAL_BAM"  # Index the BAM
      else
        rm -f "$TMP_RECAL_BAM"
        log "ApplyBQSR failed for $SAMPLE_NAME_REL"
        exit 1
      fi
    else
      log "Skipping BaseRecalibrator/ApplyBQSR; found valid $RECAL_BAM"
      if [[ ! -f "${RECAL_BAM}.bai" ]]; then
        log "Indexing $(basename "$RECAL_BAM") (existing)"
        samtools index "$RECAL_BAM"  # Ensure existing BAM is indexed
      fi
    fi
    INPUT_BAM="$RECAL_BAM"
  else
    INPUT_BAM="$SPLIT_BAM"
    log "Known sites not found; skipping BQSR"
  fi

  # -- HaplotypeCaller (GVCF) --
  local OUTPUT_GVCF="${SAMPLE_OUTPUT_DIR}/${SAMPLE_NAME_REL}_variants_${SAMPLE_TIMESTAMP}.g.vcf"
  local TMP_GVCF="${OUTPUT_GVCF}.tmp"
  if [[ -s "$OUTPUT_GVCF" ]] && grep -q "^#CHROM" "$OUTPUT_GVCF"; then
    log "GVCF for $SAMPLE_NAME_REL appears complete; skipping."
  else
    log "Running HaplotypeCaller -> $(basename "$OUTPUT_GVCF")"
    gatk --java-options "-Xmx${MEMORY}g" HaplotypeCaller \
      -R "$REF_FASTA" \
      -I "$INPUT_BAM" \
      -O "$TMP_GVCF" \
      --dont-use-soft-clipped-bases true \
      --standard-min-confidence-threshold-for-calling 20.0 \
      --emit-ref-confidence GVCF \
      --native-pair-hmm-threads "$THREADS" \
      --read-validation-stringency LENIENT &>> "$SAMPLE_LOG"
    if [[ $? -eq 0 && -s "$TMP_GVCF" ]] && grep -q "^#CHROM" "$TMP_GVCF"; then
      mv "$TMP_GVCF" "$OUTPUT_GVCF"
      log "Completed $SAMPLE_NAME_REL: Produced complete $OUTPUT_GVCF"
    else
      rm -f "$TMP_GVCF"
      log "HaplotypeCaller failed or produced incomplete GVCF for $SAMPLE_NAME_REL"
      exit 1
    fi
  fi

  unset SAMPLE_LOG
}

# Process all BAMs locally (no SLURM)
for bam in "${BAM_LIST[@]}"; do
  process_bam "$bam"
done

log "All done (timestamp: $TIMESTAMP_GLOBAL)"