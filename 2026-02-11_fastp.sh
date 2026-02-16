#!/usr/bin/env bash
# 2026-02-11_fastp.sh - batch trim FASTQ files with fastp (hard-coded dirs)
# This variant DOES NOT attempt to create a POSIX symlink on OUTPUT_ROOT (exFAT-safe).
# Instead it writes a fastp_latest.txt file containing the latest run path.
#
# Edit INPUT_DIR and OUTPUT_ROOT below if needed, then:
#   chmod +x 2026-02-11_fastp.sh
#   ./2026-02-11_fastp.sh
#
# Defaults tuned for "medium" samples on a 12-core / 48 GB machine:
#   JOBS=4, FASTP_THREADS=3  (4 * 3 = 12)
#
# Requirements: fastp in PATH

set -euo pipefail

# ----------------- EDIT THESE -----------------
INPUT_DIR="/media/jgd/jd-1tb/bioinformatics/2024-04-01_proc-enu_paper/reads"
OUTPUT_ROOT="/media/jgd/jd-1tb/bioinformatics/2024-04-01_proc-enu_paper/output"
JOBS=4
FASTP_THREADS=3
FORCE=0           # set to 1 to overwrite existing outputs within this run
# ---------------------------------------------

if ! command -v fastp >/dev/null 2>&1; then
  echo "fastp not found in PATH. Install fastp (conda/mamba/apt) and retry." >&2
  exit 1
fi

# create timestamped output directory
TIMESTAMP="$(date +%Y%m%d_%H%M%S)"
OUTPUT_DIR="${OUTPUT_ROOT%/}/fastp_${TIMESTAMP}"
mkdir -p "$OUTPUT_DIR"
echo "[$(date -u +%Y-%m-%dT%H:%M:%SZ)] Writing run outputs to: $OUTPUT_DIR"

# write a fallback pointer file instead of a symlink (exFAT-friendly)
printf "%s\n" "$OUTPUT_DIR" > "${OUTPUT_ROOT%/}/fastp_latest.txt" || {
  echo "Warning: could not write fastp_latest.txt in $OUTPUT_ROOT (permission?). Continuing..."
}

INPUT_DIR=$(realpath "$INPUT_DIR")
OUTPUT_DIR=$(realpath "$OUTPUT_DIR")

MASTER_LOG="$OUTPUT_DIR/fastp_run.log"
SUMMARY_CSV="$OUTPUT_DIR/run_summary.csv"

# header for master log and CSV
echo "FASTP RUN $(date -u +%Y-%m-%dT%H:%M:%SZ) - INPUT: $INPUT_DIR - OUTPUT: $OUTPUT_DIR" | tee -a "$MASTER_LOG"
echo "sample,mode,input1,input2,output1,output2,start_time,end_time,elapsed_seconds,exit_code,log_file" > "$SUMMARY_CSV"

# Safety: don't oversubscribe CPU
NPROC=$(nproc --all || echo 1)
if (( JOBS * FASTP_THREADS > NPROC )); then
  max_jobs=$(( NPROC / FASTP_THREADS ))
  if (( max_jobs < 1 )); then
    max_jobs=1
    FASTP_THREADS=$(( NPROC / max_jobs ))
    (( FASTP_THREADS < 1 )) && FASTP_THREADS=1
  fi
  echo "[$(date -u +%Y-%m-%dT%H:%M:%SZ)] Warning: JOBS*FASTP_THREADS > available cores (${NPROC}). Adjusting JOBS -> ${max_jobs}, threads -> ${FASTP_THREADS}." | tee -a "$MASTER_LOG"
  JOBS=${max_jobs}
fi

# collect fastq files (case-insensitive extensions)
mapfile -d $'\0' -t ALL_FILES < <(find "$INPUT_DIR" -maxdepth 1 -type f -regextype posix-extended -iregex ".*\.(fastq|fq)(\.gz)?$" -print0 | sort -z)

# derive sample base (strip gz, fq/fastq and common pair tokens)
sample_base() {
  local f="$1"
  local bn
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
  local f="$1"
  local dir bn candidate cand
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

# process paired-end sample (verbose logging)
process_pe() {
  local r1="$1" r2="$2" sample="$3"
  local out1 out2 html json logf
  out1="$OUTPUT_DIR/${sample}_R1.fastq.gz"
  out2="$OUTPUT_DIR/${sample}_R2.fastq.gz"
  html="$OUTPUT_DIR/${sample}.fastp.html"
  json="$OUTPUT_DIR/${sample}.fastp.json"
  logf="$OUTPUT_DIR/${sample}.fastp.log"

  if [[ -f "$out1" && -f "$out2" && "$FORCE" -ne 1 ]]; then
    echo "[$(date -u +%Y-%m-%dT%H:%M:%SZ)] Skipping paired sample $sample (outputs exist)." | tee -a "$MASTER_LOG"
    return
  fi

  local start_time start_ts end_time end_ts elapsed exit_code
  start_time="$(date -u +%Y-%m-%dT%H:%M:%SZ)"; start_ts="$(date +%s)"
  echo "[$start_time] START fastp (PE) -> $sample   in: $r1 , $r2" | tee -a "$MASTER_LOG"
  fastp \
    -i "$r1" -I "$r2" \
    -o "$out1" -O "$out2" \
    -w "$FASTP_THREADS" \
    --detect_adapter_for_pe \
    --qualified_quality_phred 15 \
    --length_required 30 \
    -h "$html" -j "$json" \
    2>&1 | tee -a "$logf"
  exit_code=${PIPESTATUS[0]:-0}
  end_time="$(date -u +%Y-%m-%dT%H:%M:%SZ)"; end_ts="$(date +%s)"
  elapsed=$(( end_ts - start_ts ))

  printf '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n' \
    "$sample" "PE" "$r1" "$r2" "$out1" "$out2" "$start_time" "$end_time" "$elapsed" "$exit_code" \
    >> "$SUMMARY_CSV"

  echo "[$end_time] DONE  fastp (PE) -> $sample   elapsed=${elapsed}s exit=${exit_code}" | tee -a "$MASTER_LOG"
}

# process single-end sample (verbose logging)
process_se() {
  local fq="$1" sample="$2"
  local out html json logf
  out="$OUTPUT_DIR/${sample}.trim.fastq.gz"
  html="$OUTPUT_DIR/${sample}.fastp.html"
  json="$OUTPUT_DIR/${sample}.fastp.json"
  logf="$OUTPUT_DIR/${sample}.fastp.log"

  if [[ -f "$out" && "$FORCE" -ne 1 ]]; then
    echo "[$(date -u +%Y-%m-%dT%H:%M:%SZ)] Skipping single sample $sample (output exists)." | tee -a "$MASTER_LOG"
    return
  fi

  local start_time start_ts end_time end_ts elapsed exit_code
  start_time="$(date -u +%Y-%m-%dT%H:%M:%SZ)"; start_ts="$(date +%s)"
  echo "[$start_time] START fastp (SE) -> $sample   in: $fq" | tee -a "$MASTER_LOG"
  fastp \
    -i "$fq" -o "$out" \
    -w "$FASTP_THREADS" \
    --detect_adapter_for_pe \
    --qualified_quality_phred 15 \
    --length_required 30 \
    -h "$html" -j "$json" \
    2>&1 | tee -a "$logf"
  exit_code=${PIPESTATUS[0]:-0}
  end_time="$(date -u +%Y-%m-%dT%H:%M:%SZ)"; end_ts="$(date +%s)"
  elapsed=$(( end_ts - start_ts ))

  printf '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n' \
    "$sample" "SE" "$fq" "" "$out" "" "$start_time" "$end_time" "$elapsed" "$exit_code" \
    >> "$SUMMARY_CSV"

  echo "[$end_time] DONE  fastp (SE) -> $sample   elapsed=${elapsed}s exit=${exit_code}" | tee -a "$MASTER_LOG"
}

declare -A seen

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

# wait for remaining background jobs
wait

echo "[$(date -u +%Y-%m-%dT%H:%M:%SZ)] All fastp jobs finished. Results: $OUTPUT_DIR" | tee -a "$MASTER_LOG"
echo "Summary CSV: $SUMMARY_CSV" | tee -a "$MASTER_LOG"