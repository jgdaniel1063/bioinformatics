#!/usr/bin/env bash
# kallisto_quant.sh - batch quantify FASTQ files with kallisto
# - STAR-like features: timestamped runs, per-sample folders, master run log, CSV summary
# - Builds kallisto index automatically if missing (requires TRANSCRIPTOME_FASTA)
# - Detects paired-end files; optionally supports single-end (requires FRAGMENT_LENGTH & SD)
# - Per-sample logs (tee'd), verbose timestamped events
# - Output organization: OUTPUT_ROOT/RESULTS_NAME/kallisto_YYYYMMDD_HHMMSS/<sample>/
#
# Edit the "EDIT THESE" block below before running.
# Usage:
#   chmod +x kallisto_quant.sh
#   ./kallisto_quant.sh
#
# Requirements:
# - kallisto in PATH
# - GNU coreutils for mktemp/sort; bash 4+ recommended

set -euo pipefail

# ----------------- EDIT THESE -----------------
# Input / output
INPUT_DIR="/media/jgd/jd-1tb/bioinformatics/2024-04-01_proc-enu_paper/trimming"  # folder with FASTQ files (fastp output)
OUTPUT_ROOT="/home/jgd/Documents/bioinformatics_working/output"                    # parent for results
RESULTS_NAME="kallisto_results"                                                   # named folder under OUTPUT_ROOT that will contain timestamped runs

# Kallisto index / transcriptome
KALLISTO_INDEX="/media/jgd/jd-1tb/bioinformatics/reference/kallisto_index"  # path to index (if already built) OR where to write it
TRANSCRIPTOME_FASTA="/media/jgd/jd-1tb/bioinformatics/reference/mrna.fa.gz"   # REQUIRED if index missing (e.g. /path/to/transcripts.fa or .fa.gz)

# Single-end support (only used if you have SE FASTQs)
SINGLE_MODE=0            # set to 1 to allow running kallisto quant on single-end FASTQs
FRAGMENT_LENGTH=300      # required when SINGLE_MODE=1
FRAGMENT_SD=30           # required when SINGLE_MODE=1

# Kallisto / parallelism settings
JOBS=4                  # number of concurrent kallisto jobs
THREADS_PER_JOB=3       # threads per kallisto process (-t); JOBS*THREADS_PER_JOB <= CPU cores ideally
BOOTSTRAPS=100          # number of bootstrap samples (set 0 to skip bootstraps)
KMER=31                 # k-mer length to use if building index

FORCE=0                 # set to 1 to overwrite per-sample outputs within a run
# ---------------------------------------------

# helper: timestamp
ts() { date -u +%Y-%m-%dT%H:%M:%SZ; }

# quick checks
if ! command -v kallisto >/dev/null 2>&1; then
  echo "[$(ts)] ERROR: kallisto not found in PATH. Install kallisto and retry." >&2
  exit 1
fi

# create named results folder and timestamped run dir
RESULTS_DIR="${OUTPUT_ROOT%/}/${RESULTS_NAME%/}"
mkdir -p "$RESULTS_DIR"
TIMESTAMP="$(date +%Y%m%d_%H%M%S)"
RUN_DIR="${RESULTS_DIR%/}/kallisto_${TIMESTAMP}"
mkdir -p "$RUN_DIR"
echo "[$(ts)] Writing kallisto run outputs to: $RUN_DIR"

# write pointer file (exFAT-safe)
printf "%s\n" "$RUN_DIR" > "${RESULTS_DIR%/}/kallisto_latest.txt" 2>/dev/null || \
  echo "[$(ts)] Warning: could not write kallisto_latest.txt in $RESULTS_DIR (permission?). Continuing..."

INPUT_DIR=$(realpath "$INPUT_DIR")
RUN_DIR=$(realpath "$RUN_DIR")

MASTER_LOG="$RUN_DIR/kallisto_run.log"
SUMMARY_CSV="$RUN_DIR/run_summary.csv"

echo "KALLISTO RUN $(ts) - INPUT: $INPUT_DIR - RESULTS_DIR: $RESULTS_DIR - RUN: $RUN_DIR" | tee -a "$MASTER_LOG"
echo "sample,mode,input1,input2,output_dir,start_time,end_time,elapsed_seconds,exit_code,log_file" > "$SUMMARY_CSV"

# CPU safety: adjust JOBS if oversubscribing
NPROC=$(nproc --all || echo 1)
if (( JOBS * THREADS_PER_JOB > NPROC )); then
  max_jobs=$(( NPROC / THREADS_PER_JOB ))
  if (( max_jobs < 1 )); then
    max_jobs=1
    THREADS_PER_JOB=$(( NPROC / max_jobs ))
    (( THREADS_PER_JOB < 1 )) && THREADS_PER_JOB=1
  fi
  echo "[$(ts)] Warning: JOBS*THREADS_PER_JOB > available cores (${NPROC}). Adjusting JOBS -> ${max_jobs}, THREADS_PER_JOB -> ${THREADS_PER_JOB}." | tee -a "$MASTER_LOG"
  JOBS=${max_jobs}
fi

# Detect/prepare kallisto index
index_present() {
  local idx="$1"
  [[ -f "$idx" ]] && return 0
  return 1
}

maybe_uncompress_fa() {
  local ref="$1"
  if [[ "$ref" == *.gz ]]; then
    local tmpf
    tmpf="$(mktemp --tmpdir kallisto_ref_XXXXXX.fa)"
    echo "[$(ts)] Uncompressing transcriptome $ref -> $tmpf" | tee -a "$MASTER_LOG"
    gzip -dc "$ref" > "$tmpf"
    echo "$tmpf"
  else
    echo "$ref"
  fi
}

if ! index_present "$KALLISTO_INDEX"; then
  echo "[$(ts)] kallisto index not found at: $KALLISTO_INDEX" | tee -a "$MASTER_LOG"
  if [[ -z "$TRANSCRIPTOME_FASTA" ]]; then
    echo "[$(ts)] ERROR: TRANSCRIPTOME_FASTA not set. To build an index set TRANSCRIPTOME_FASTA to transcripts.fa (can be gzipped)." | tee -a "$MASTER_LOG"
    exit 1
  fi

  # uncompress if needed
  FASTA_TO_USE="$TRANSCRIPTOME_FASTA"
  CLEAN_TMP_REF=0
  if [[ "$TRANSCRIPTOME_FASTA" == *.gz ]]; then
    FASTA_TO_USE="$(maybe_uncompress_fa "$TRANSCRIPTOME_FASTA")"
    CLEAN_TMP_REF=1
  fi

  echo "[$(ts)] Building kallisto index at: $KALLISTO_INDEX (k=$KMER)" | tee -a "$MASTER_LOG"
  kallisto index -i "$KALLISTO_INDEX" -k "$KMER" "$FASTA_TO_USE" 2>&1 | tee -a "$MASTER_LOG"
  rc=${PIPESTATUS[0]:-0}
  if (( rc != 0 )); then
    echo "[$(ts)] ERROR: kallisto index build failed (exit $rc)." | tee -a "$MASTER_LOG"
    [[ "$CLEAN_TMP_REF" -eq 1 ]] && rm -f "$FASTA_TO_USE"
    exit 1
  fi
  echo "[$(ts)] kallisto index built: $KALLISTO_INDEX" | tee -a "$MASTER_LOG"
  [[ "$CLEAN_TMP_REF" -eq 1 ]] && rm -f "$FASTA_TO_USE"
fi

# collect fastq files (case-insensitive)
mapfile -d $'\0' -t ALL_FILES < <(find "$INPUT_DIR" -maxdepth 1 -type f -regextype posix-extended -iregex ".*\.(fastq|fq)(\.gz)?$" -print0 | sort -z)

# derive sample base (strip gz, fq/fastq and common pair tokens)
sample_base() {
  local f bn
  bn="$(basename "$f")"
  bn="${bn%.gz}"
  bn="${bn%.*}"
  bn="${bn%_R1_001}"; bn="${bn%_R2_001}"
  bn="${bn%_R1}"; bn="${bn%_R2}"
  bn="${bn%_1}"; bn="${bn%_2}"
  echo "$bn"
}

# find mate by replacing common tokens
find_mate() {
  local f dir bn candidate cand
  dir="$(dirname "$f")"
  bn="$(basename "$f")"
  local -a froms=("_R1_001" "_R1" "_1" "R1" "_read1" "_READ1" "_r1")
  local -a tos=( "_R2_001" "_R2" "_2" "R2" "_read2" "_READ2" "_r2" )
  for i in "${!froms[@]}"; do
    candidate="${bn/${froms[i]}/${tos[i]}}"
    if [[ "$candidate" != "$bn" ]]; then
      [[ -f "$dir/$candidate" ]] && { echo "$dir/$candidate"; return 0; }
      [[ -f "$dir/${candidate}.gz" ]] && { echo "$dir/${candidate}.gz"; return 0; }
    fi
  done

  cand="${bn//R1/R2}"
  if [[ "$cand" != "$bn" ]]; then
    [[ -f "$dir/$cand" ]] && { echo "$dir/$cand"; return 0; }
    [[ -f "$dir/${cand}.gz" ]] && { echo "$dir/${cand}.gz"; return 0; }
  fi

  cand="${bn//_1/_2}"
  if [[ "$cand" != "$bn" ]]; then
    [[ -f "$dir/$cand" ]] && { echo "$dir/$cand"; return 0; }
    [[ -f "$dir/${cand}.gz" ]] && { echo "$dir/${cand}.gz"; return 0; }
  fi

  return 1
}

# concurrency control
wait_for_slot() {
  while :; do
    local n
    n=$(jobs -pr | wc -l)
    if (( n < JOBS )); then
      break
    fi
    sleep 0.4
  done
}

# run kallisto for paired-end
process_pe() {
  local r1="$1" r2="$2" sample="$3"
  local sample_dir logf start_time start_ts end_time end_ts elapsed exit_code kallisto_cmd
  sample_dir="$RUN_DIR/$sample"
  mkdir -p "$sample_dir"
  logf="$sample_dir/${sample}.kallisto.log"

  if [[ -d "$sample_dir/abundance.tsv" || -f "$sample_dir/abundance.tsv" ]] && [[ "$FORCE" -ne 1 ]]; then
    echo "[$(ts)] Skipping paired sample $sample (output exists): $sample_dir" | tee -a "$MASTER_LOG"
    return
  fi

  echo "[$(ts)] START kallisto (PE) -> $sample   in: $r1 , $r2" | tee -a "$MASTER_LOG"
  start_time="$(ts)"; start_ts="$(date +%s)"

  # build kallisto command
  kallisto_cmd=(kallisto quant -i "$KALLISTO_INDEX" -o "$sample_dir" -t "$THREADS_PER_JOB" --bias)
  if (( BOOTSTRAPS > 0 )); then
    kallisto_cmd+=( -b "$BOOTSTRAPS" )
  fi
  kallisto_cmd+=( "$r1" "$r2" )

  # run and tee
  "${kallisto_cmd[@]}" 2>&1 | tee -a "$logf"
  exit_code=${PIPESTATUS[0]:-0}

  end_time="$(ts)"; end_ts="$(date +%s)"
  elapsed=$(( end_ts - start_ts ))

  printf '%s,%s,%s,%s,%s,%s,%s,%s,%s\n' \
    "$sample" "PE" "$r1" "$r2" "$sample_dir" "$start_time" "$end_time" "$elapsed" "$exit_code" \
    >> "$SUMMARY_CSV"

  echo "[$(ts)] DONE kallisto (PE) -> $sample   elapsed=${elapsed}s exit=${exit_code}" | tee -a "$MASTER_LOG"
}

# run kallisto for single-end (requires FRAGMENT_LENGTH & FRAGMENT_SD)
process_se() {
  local fq="$1" sample="$2"
  local sample_dir logf start_time start_ts end_time end_ts elapsed exit_code kallisto_cmd
  sample_dir="$RUN_DIR/$sample"
  mkdir -p "$sample_dir"
  logf="$sample_dir/${sample}.kallisto.log"

  if [[ "$SINGLE_MODE" -ne 1 ]]; then
    echo "[$(ts)] WARNING: single-end sample $sample detected but SINGLE_MODE=0. Skipping." | tee -a "$MASTER_LOG"
    return
  fi

  if [[ -z "$FRAGMENT_LENGTH" || -z "$FRAGMENT_SD" ]]; then
    echo "[$(ts)] ERROR: FRAGMENT_LENGTH/FRAGMENT_SD must be set for single-end kallisto quant. Skipping $sample." | tee -a "$MASTER_LOG"
    return
  fi

  if [[ -f "$sample_dir/abundance.tsv" ]] && [[ "$FORCE" -ne 1 ]]; then
    echo "[$(ts)] Skipping single sample $sample (output exists): $sample_dir" | tee -a "$MASTER_LOG"
    return
  fi

  echo "[$(ts)] START kallisto (SE) -> $sample   in: $fq" | tee -a "$MASTER_LOG"
  start_time="$(ts)"; start_ts="$(date +%s)"

  kallisto_cmd=(kallisto quant -i "$KALLISTO_INDEX" -o "$sample_dir" -t "$THREADS_PER_JOB" --single -l "$FRAGMENT_LENGTH" -s "$FRAGMENT_SD" --bias)
  if (( BOOTSTRAPS > 0 )); then
    kallisto_cmd+=( -b "$BOOTSTRAPS" )
  fi
  kallisto_cmd+=( "$fq" )

  "${kallisto_cmd[@]}" 2>&1 | tee -a "$logf"
  exit_code=${PIPESTATUS[0]:-0}

  end_time="$(ts)"; end_ts="$(date +%s)"
  elapsed=$(( end_ts - start_ts ))

  printf '%s,%s,%s,%s,%s,%s,%s,%s,%s\n' \
    "$sample" "SE" "$fq" "" "$sample_dir" "$start_time" "$end_time" "$elapsed" "$exit_code" \
    >> "$SUMMARY_CSV"

  echo "[$(ts)] DONE kallisto (SE) -> $sample   elapsed=${elapsed}s exit=${exit_code}" | tee -a "$MASTER_LOG"
}

declare -A seen

# iterate over FASTQ files and dispatch jobs
for f in "${ALL_FILES[@]}"; do
  [[ -n "${f// }" ]] || continue
  [[ -n "${seen[$f]:-}" ]] && continue

  bn="$(basename "$f")"

  if [[ "$bn" =~ (_R1_001|_R1|_1|R1|_read1|_READ1|_r1) ]]; then
    mate="$(find_mate "$f" || true)"
    sample="$(sample_base "$f")"
    seen["$f"]=1
    if [[ -n "$mate" && -f "$mate" ]]; then
      seen["$mate"]=1
      wait_for_slot
      process_pe "$f" "$mate" "$sample" &
    else
      wait_for_slot
      process_se "$f" "$sample" &
    fi
  else
    # if it's likely an R2 and R1 exists skip it; else treat as SE
    if [[ "$bn" =~ (_R2_001|_R2|_2|R2|_read2|_READ2|_r2) ]]; then
      maybe_r1="${f//R2/R1}"
      maybe_r1="${maybe_r1//_2/_1}"
      if [[ -f "$maybe_r1" || -f "${maybe_r1}.gz" ]]; then
        seen["$f"]=1
        continue
      fi
    fi
    sample="$(sample_base "$f")"
    seen["$f"]=1
    wait_for_slot
    process_se "$f" "$sample" &
  fi
done

# wait for remaining jobs
wait

echo "[$(ts)] All kallisto jobs finished. Results: $RUN_DIR" | tee -a "$MASTER_LOG"
echo "Summary CSV: $SUMMARY_CSV" | tee -a "$MASTER_LOG"