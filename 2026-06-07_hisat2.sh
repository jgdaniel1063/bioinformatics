#!/usr/bin/env bash
# Sassy HISAT2 pipeline for Jeff
# Now with a SANITY METER that fills as samples finish.

set -euo pipefail

###############################################
# COLORS & UTILITIES
###############################################
RED="\033[0;31m"
GREEN="\033[0;32m"
YELLOW="\033[1;33m"
BLUE="\033[0;34m"
MAGENTA="\033[0;35m"
CYAN="\033[0;36m"
RESET="\033[0m"

timestamp() { date +"%Y-%m-%d %H:%M:%S"; }

say() {
    local color="$1"; shift
    echo -e "${color}[$(timestamp)] $*${RESET}"
}

###############################################
# USER SETTINGS
###############################################
INPUT_DIR="/home/jgd/Documents/bioinformatics_working/output/fastp_latest"
OUTPUT_ROOT="/home/jgd/Documents/bioinformatics_working/output"
HISAT2_INDEX="/path/to/hisat2/index/prefix"
THREADS=8
JOBS=3
FORCE=0

say "$CYAN" "Welcome back, Jeff. Initializing sassy HISAT2 pipeline."
say "$CYAN" "Input directory: $INPUT_DIR"
say "$CYAN" "HISAT2 index prefix: $HISAT2_INDEX"
say "$CYAN" "Max parallel samples: $JOBS, threads per HISAT2: $THREADS"

###############################################
# OUTPUT SETUP
###############################################
TS=$(date +%Y%m%d_%H%M%S)
OUTPUT_DIR="${OUTPUT_ROOT%/}/hisat2_${TS}"
mkdir -p "$OUTPUT_DIR"

say "$GREEN" "Output directory: $OUTPUT_DIR"
echo "$OUTPUT_DIR" > "${OUTPUT_ROOT%/}/hisat2_latest.txt"

SUMMARY_CSV="$OUTPUT_DIR/run_summary.csv"
echo "sample,input1,input2,bam,start,end,elapsed,exit" > "$SUMMARY_CSV"

###############################################
# SANITY METER SETUP
###############################################
TOTAL_SAMPLES=$(find "$INPUT_DIR" -type f -name "*_R1.fastq.gz" | wc -l)
COMPLETED=0

draw_sanity_meter() {
    local filled=$(( COMPLETED * 20 / TOTAL_SAMPLES ))
    local empty=$(( 20 - filled ))
    local bar="$(printf "%${filled}s" | tr ' ' '█')$(printf "%${empty}s" | tr ' ' '░')"
    say "$MAGENTA" "Sanity Meter: [$bar]  $COMPLETED / $TOTAL_SAMPLES samples aligned"
}

###############################################
# JOB SLOT MANAGER
###############################################
wait_for_slot() {
    while :; do
        local running
        running=$(jobs -pr | wc -l || echo 0)
        say "$MAGENTA" "Active jobs: $running / $JOBS. Deep breaths."
        (( running < JOBS )) && break
        sleep 60
    done
}

###############################################
# HISAT2 RUNNER (SASSY)
###############################################
run_hisat2() {
    local r1="$1" r2="$2" sample="$3"
    local bam logf

    bam="$OUTPUT_DIR/${sample}.sorted.bam"
    logf="$OUTPUT_DIR/${sample}.log"

    if [[ -f "$bam" && "$FORCE" -eq 0 ]]; then
        say "$YELLOW" "Skipping $sample — BAM already exists. Past Jeff was on fire."
        return
    fi

    say "$BLUE" "Launching HISAT2 for sample: $sample"
    say "$BLUE" "  R1: $r1"
    say "$BLUE" "  R2: $r2"
    say "$BLUE" "  Output BAM: $bam"

    local start end elapsed exit_code
    start=$(date +%s)

    hisat2 \
        -x "$HISAT2_INDEX" \
        -1 "$r1" -2 "$r2" \
        --threads "$THREADS" \
        --dta \
        2>&1 | tee "$logf" | samtools view -bS - \
        | samtools sort -@ "$THREADS" -o "$bam"

    exit_code=${PIPESTATUS[0]}
    end=$(date +%s)
    elapsed=$(( end - start ))

    if [[ "$exit_code" -eq 0 ]]; then
        say "$GREEN" "Sample $sample aligned and sorted in ${elapsed}s. Gorgeous work."
    else
        say "$RED" "Sample $sample failed with exit code $exit_code. Drama."
    fi

    say "$BLUE" "Indexing BAM for sample: $sample"
    samtools index "$bam"

    COMPLETED=$(( COMPLETED + 1 ))
    draw_sanity_meter

    printf "%s,%s,%s,%s,%s,%s,%s,%s\n" \
        "$sample" "$r1" "$r2" "$bam" \
        "$(date -u -d @$start +%FT%TZ)" \
        "$(date -u -d @$end +%FT%TZ)" \
        "$elapsed" "$exit_code" >> "$SUMMARY_CSV"
}

###############################################
# MAIN LOOP
###############################################
say "$CYAN" "Scanning $INPUT_DIR for trimmed FASTQs."

while IFS= read -r r1; do
    [[ -z "$r1" ]] && continue

    bn=$(basename "$r1")
    sample="${bn%_R1.fastq.gz}"
    r2="${r1/_R1.fastq.gz/_R2.fastq.gz}"

    if [[ ! -f "$r2" ]]; then
        say "$RED" "Missing R2 for sample $sample — skipping."
        continue
    fi

    say "$YELLOW" "Detected sample: $sample. Time to align."
    wait_for_slot
    run_hisat2 "$r1" "$r2" "$sample" &

done < <(find "$INPUT_DIR" -type f -name "*_R1.fastq.gz" | sort)

say "$CYAN" "All HISAT2 jobs launched. Waiting for the alignment gods."
wait

say "$GREEN" "All HISAT2 jobs finished. BAMs sorted, indexed, and ready for variant calling."
say "$GREEN" "Output directory: $OUTPUT_DIR"
say "$GREEN" "Summary CSV: $SUMMARY_CSV"
say "$CYAN" "Your sanity meter is full. Go celebrate, king."
