#!/usr/bin/env bash
#
# run_fastqc_multiqc_annotated.sh
#
# Heavily annotated wrapper to run FastQC on a directory of sequencing files
# (FASTQ / FASTQ.GZ; single- or paired-end), record per-file success/failure,
# and aggregate results with MultiQC. The script runs FastQC in parallel using
# GNU parallel or xargs as a fallback and organizes outputs into a timestamped
# run directory for reproducibility.
#
# Goals / Design:
#  - Discover input FASTQ files (optionally recursively)
#  - Run FastQC on each file in parallel and capture per-file status
#  - Aggregate outputs into a structured run directory
#  - Run MultiQC over the FastQC results to produce a single HTML report
#  - Provide robust logging and a simple status TSV for downstream checks
#
# Important notes:
#  - FastQC analyzes individual read files. For paired-end data you should run
#    FastQC on both mates; downstream MultiQC can show side-by-side results.
#  - This script DOES NOT perform trimming or adapter removal. Run trimming
#    before FastQC if you want to evaluate post-trimming quality.
#  - FastQC and MultiQC must be installed and on PATH (conda/pip recommended).
#
# Example usage:
#   ./run_fastqc_multiqc_annotated.sh
#   ./run_fastqc_multiqc_annotated.sh -i /data/reads -o /results -t 8 --recursive
#
# Output layout:
#  <out_parent>/<timestamp>/
#    run.log                 combined stdout/stderr
#    fastqc/                 per-file FastQC outputs (zips and html)
#    multiqc_report/         MultiQC report (multiqc_report.html)
#    fastqc_status.tsv       tab-separated: file<TAB>OK|FAIL
#    failed_fastqc.txt       list of files that failed FastQC
#    run_helper_script.sh    helper script used to run FastQC in parallel (saved for debugging)
#
# Recommendations:
#  - Use a recent FastQC and MultiQC version.
#  - For reproducibility, pin package versions in a conda env and record env details.
#  - For very large datasets, consider batching files per-run or submitting to a cluster.
#
set -euo pipefail

# ---------------------------
# Basic defaults & usage
# ---------------------------
progname=$(basename "$0")

usage() {
  cat <<USAGE
Usage: $progname [options]

Options:
  -i, --input DIR       Input directory containing FASTQ/FQ files
                        (default: /home/jgd/Documents/2026-04-01_proc-enu_bioinformatics_results/reads)
  -o, --out-parent DIR  Parent output directory; a timestamped subdir will be created inside
                        (default: /home/jgd/Documents/bioinformatics_working/output)
  -t, --threads N       Number of parallel jobs (default: 12)
  -r, --recursive       Search input directory recursively (default: false)
  -k, --keep            Keep per-sample FastQC HTML files (default: keep)
  -h, --help            Show this help message and exit

Example:
  $progname -t 8 --recursive

Notes:
 - This script only runs FastQC and MultiQC; it does not perform trimming/quality filtering.
 - For paired-end data you will see two FastQC entries per sample (R1 and R2).
USAGE
}

# ---------------------------
# Defaults (edit if needed)
# ---------------------------
INPUT_DIR="/home/jgd/Documents/2026-04-01_proc-enu_bioinformatics_results/reads"
OUT_PARENT="/home/jgd/Documents/bioinformatics_working/output"
THREADS=12
RECURSIVE=false
KEEP=true   # keep per-sample FastQC HTMLs; set to false to remove them after MultiQC

# ---------------------------
# Parse CLI args (simple)
# ---------------------------
while [[ $# -gt 0 ]]; do
  case "$1" in
    -i|--input) INPUT_DIR="$2"; shift 2 ;;
    -o|--out-parent) OUT_PARENT="$2"; shift 2 ;;
    -t|--threads) THREADS="$2"; shift 2 ;;
    -r|--recursive) RECURSIVE=true; shift ;;
    -k|--keep) KEEP=true; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown option: $1"; usage; exit 1 ;;
  esac
done

# ---------------------------
# Minimal environment & tools checks
# ---------------------------
command -v fastqc >/dev/null 2>&1 || { echo "ERROR: fastqc not found in PATH. Install fastqc and retry." >&2; exit 1; }
command -v multiqc >/dev/null 2>&1 || { echo "ERROR: multiqc not found in PATH. Install multiqc and retry." >&2; exit 1; }
# GNU parallel is optional but recommended; fallback to xargs is provided
if command -v parallel >/dev/null 2>&1; then
  PARALLEL_BIN="$(command -v parallel)"
  PARALLEL_AVAILABLE=true
else
  PARALLEL_AVAILABLE=false
fi

# ---------------------------
# Create run directories
# ---------------------------
timestamp=$(date +"%Y%m%d_%H%M%S")
OUT_DIR="${OUT_PARENT%/}/$timestamp"

# Avoid collision (extremely unlikely)
if [ -d "$OUT_DIR" ]; then
  suffix=$(date +%s%N)
  OUT_DIR="${OUT_DIR}_$suffix"
fi

mkdir -p "$OUT_DIR"
FASTQC_OUT="$OUT_DIR/fastqc"
MULTIQC_OUT="$OUT_DIR/multiqc_report"
mkdir -p "$FASTQC_OUT" "$MULTIQC_OUT"

# Log and status files
LOG="$OUT_DIR/run.log"
STATUS_TSV="$OUT_DIR/fastqc_status.tsv"
FAILED_LIST="$OUT_DIR/failed_fastqc.txt"
HELPER_SCRIPT="$OUT_DIR/run_helper_script.sh"

# Initialize log & status
: > "$LOG"
printf "file\tstatus\n" > "$STATUS_TSV"
: > "$FAILED_LIST"

# Logging helper to prepend timestamps to messages (writes to both console and log)
log() {
  printf '[%s] %s\n' "$(date '+%F %T')" "$*" | tee -a "$LOG"
}

log "Starting run"
log "Input directory: $INPUT_DIR"
log "Output directory: $OUT_DIR"
log "Threads: $THREADS"
log "Recursive: $RECURSIVE"
log "Keep per-sample HTML: $KEEP"

# ---------------------------
# Discover FASTQ files
# ---------------------------
# We accept files with extensions: .fastq .fastq.gz .fq .fq.gz (case-insensitive).
# Use find -print0 and mapfile to robustly handle filenames with whitespace/newlines.
if [ "$RECURSIVE" = true ]; then
  FIND_MAXDEPTH=""
else
  FIND_MAXDEPTH="-maxdepth 1"
fi

# If the input directory doesn't exist, stop early with a helpful message.
if [ ! -d "$INPUT_DIR" ]; then
  log "ERROR: Input directory does not exist: $INPUT_DIR"
  exit 1
fi

# Collect files into an array (null-delimited)
mapfile -d '' FILES < <(
  find "$INPUT_DIR" $FIND_MAXDEPTH -type f \( -iname '*.fastq' -o -iname '*.fastq.gz' -o -iname '*.fq' -o -iname '*.fq.gz' \) -print0 2>/dev/null \
  | sort -z
)

if [ "${#FILES[@]}" -eq 0 ]; then
  log "No FASTQ/FQ files found in '$INPUT_DIR'. Exiting."
  exit 1
fi

log "Found ${#FILES[@]} FASTQ/FQ file(s) (showing up to 10):"
for ((i=0;i<${#FILES[@]} && i<10; i++)); do
  log "  ${FILES[i]}"
done
if [ "${#FILES[@]}" -gt 10 ]; then
  log "  ... and $(( ${#FILES[@]} - 10 )) more files"
fi

# ---------------------------
# Prepare helper script to run FastQC on a single file
# ---------------------------
# Rationale:
#  - Running FastQC inside a small helper script makes it easy to handle per-file
#    logging and status updates atomically.
#  - Using tee/append to the STATUS_TSV may risk interleaved writes; we keep operations
#    simple and rely on the helper being run one per process by parallel/xargs.
cat > "$HELPER_SCRIPT" <<'SH'
#!/usr/bin/env bash
set -euo pipefail
# Helper invoked as: helper.sh <file>
f="$1"
# FASTQC_OUT, LOG, STATUS_TSV, FAILED_LIST are exported into the environment
# Run FastQC with single thread for this invocation; the caller controls parallelism.
if fastqc --threads 1 -o "$FASTQC_OUT" "$f" 2>>"$LOG"; then
  printf '%s\tOK\n' "$f" >> "$STATUS_TSV"
else
  printf '%s\tFAIL\n' "$f" >> "$STATUS_TSV"
  printf '%s\n' "$f" >> "$FAILED_LIST"
fi
SH
chmod +x "$HELPER_SCRIPT"

# Export variables used by the helper (GNU parallel or xargs will inherit env)
export FASTQC_OUT LOG STATUS_TSV FAILED_LIST

# ---------------------------
# Run FastQC in parallel
# ---------------------------
log "Running FastQC in parallel (this may take a while)..."

# Two modes:
#  - GNU parallel: handles null-delimited input cleanly and is robust for many files
#  - xargs fallback: widely available; we use -0 / -n1 / -P for parallel execution
if $PARALLEL_AVAILABLE; then
  log "Using GNU parallel: $PARALLEL_BIN"
  # print files as null-delimited and pipe into parallel -0
  printf '%s\0' "${FILES[@]}" | parallel -0 -j "$THREADS" --no-notice "$HELPER_SCRIPT" {}
else
  log "GNU parallel not found; falling back to xargs (may be slightly less robust)"
  printf '%s\0' "${FILES[@]}" | xargs -0 -n1 -P "$THREADS" -I{} bash "$HELPER_SCRIPT" {}
fi

log "FastQC step finished. Statuses written to: $STATUS_TSV"

# ---------------------------
# Summarize failures
# ---------------------------
fail_count=0
if [ -f "$FAILED_LIST" ]; then
  # wc -l can return spaces; trim them
  fail_count=$(wc -l < "$FAILED_LIST" | tr -d '[:space:]' || echo 0)
fi

if [ "$fail_count" -gt 0 ]; then
  log "WARNING: $fail_count file(s) failed FastQC. See: $FAILED_LIST"
  head -n 50 "$FAILED_LIST" | sed 's/^/  /' >> "$LOG"
else
  log "All FastQC jobs reported OK in helper script (note: this indicates FastQC invocation returned 0)."
fi

# ---------------------------
# Aggregate with MultiQC
# ---------------------------
log "Running MultiQC to aggregate FastQC outputs..."
# MultiQC scans the FASTQC_OUT directory and creates a report. We log stderr to LOG.
if multiqc "$FASTQC_OUT" -n multiqc_report -o "$MULTIQC_OUT" 2>>"$LOG"; then
  log "MultiQC report created: ${MULTIQC_OUT}/multiqc_report.html"
else
  log "MultiQC returned non-zero exit (see $LOG). MultiQC may still have produced partial results."
fi

# ---------------------------
# Optionally remove per-sample HTML files
# ---------------------------
if [ "$KEEP" != true ]; then
  log "Removing per-sample FastQC HTML files (keeping zip files and MultiQC output)..."
  find "$FASTQC_OUT" -type f -iname '*.html' -print0 | xargs -0 -r rm -f --
fi

# ---------------------------
# Post-run notes and advice
# ---------------------------
log "Run completed. Outputs are in: $OUT_DIR"
if [ "$fail_count" -gt 0 ]; then
  log "Failed files list: $FAILED_LIST"
fi

cat <<'END_ADVICE' | tee -a "$LOG"

ADVICE & TROUBLESHOOTING
------------------------
1) FastQC failures:
   - If FastQC returns a non-zero exit code for a file, the helper records it as FAIL.
   - Common causes:
       * corrupt or truncated FASTQ (incomplete gzip file)
       * permissions issues reading files
       * extremely large single reads or malformed FASTA/FASTQ headers
   - Inspect the individual FastQC stderr lines in the run.log for the specific error.

2) Paired-end data:
   - This script treats files independently. You will see two results per sample (R1 and R2).
   - For per-sample QC summaries you may want to group paired files before downstream steps.

3) Trimming & adapter removal:
   - Run adapter trimming (e.g., TrimGalore/Cutadapt) before FastQC if you need to assess post-trim quality.
   - Rerun FastQC on trimmed files to compare pre/post trimming.

4) Performance:
   - FastQC does not scale linearly with threads across multiple files; running many single-threaded jobs
     in parallel (via GNU parallel) is usually the most efficient approach.
   - I/O can be the bottleneck on networked filesystems; consider copying files to local scratch when possible.

5) Reproducibility:
   - Record tool versions used:
        fastqc --version
        multiqc --version
     Consider writing them into ${OUT_DIR}/tools_versions.txt for provenance.

6) Extending the script:
   - Add optional pre-filtering (e.g., only run FastQC on R1 files, or group by sample name).
   - Produce a summary CSV with FastQC module pass/fail counts by parsing the FastQC zip contents.
   - Support compressed stdout-only logs for long-running runs.

7) Example pipeline step:
   - Typical fastq QC and trimming flow:
        1) run this FastQC to inspect raw reads
        2) run trimming (TrimGalore/Cutadapt)
        3) run this FastQC again on trimmed files and compare with MultiQC

END_ADVICE

log "Done."
