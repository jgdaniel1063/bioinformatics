#!/bin/bash
#
# dexseq_count_bam_folder_annotated.sh
#
# Heavily annotated wrapper to run DEXSeq's dexseq_count.py across a folder
# of BAM files, producing one per-sample count file suitable for DEXSeq/R.
#
# Key responsibilities:
#  - validate environment and inputs
#  - locate BAM files (recursive)
#  - call dexseq_count.py on each BAM to produce exon-level counts
#  - report success / failure for each sample
#
# Notes:
#  - dexseq_count.py is a Python script shipped with the DEXSeq R package.
#    Many installations expose it via a conda environment (activate first),
#    or it can be called using the full path inside the R package installation:
#      $R_HOME/library/DEXSeq/python_scripts/dexseq_count.py
#
#  - dexseq_count.py assumptions:
#     * Input GTF must be the "flattened" GTF produced by dexseq_prepare_annotation.py.
#     * Input BAM should be coordinate-sorted and (for best results) indexed (.bai).
#     * The "-s" argument controls strandness: "yes", "no", or "reverse".
#     * For paired-end libraries additional flags may be needed (see dexseq_count.py -h).
#
#  - This annotated version improves robustness (quoting, safer file iteration)
#    and explains caveats / recommended checks at each step.
#

# -------------------- HARDCODED VALUES (edit here if needed) --------------------
# Important: change these to match your environment before running.
BAM_FOLDER="/home/jgd/Documents/bioinformatics_working/2025-11-06_enu-proc/align/star-ensembl"  # folder to search for *.bam
FLATTENED_GTF="/home/jgd/Documents/bioinformatics_working/ref/danio_rerio.grcz11.115.dexseq.gtf"  # flattened GTF produced by dexseq_prepare_annotation.py
OUTPUT_DIR="/home/jgd/Documents/bioinformatics_working/output/dexseq_counts"  # where to write *_dexseq_counts.txt files
DEXSEQ_COUNT_SCRIPT="/home/jgd/miniconda3/lib/R/library/DEXSeq/python_scripts/dexseq_count.py"  # explicit path to dexseq_count.py (change if different)
STRANDED="no"  # "yes" / "no" / "reverse" — set according to library prep
# -------------------------------------------------------------------------------

# -------------------- RUNTIME SAFETY & ENVIRONMENT CHECKS -----------------------
# These checks fail fast with clear messages if requirements are not met.
# They help avoid long-running jobs that later error with obscure messages.

# Prefer python3 where available. Many systems have `python` symlinked to python2.
PYTHON=$(command -v python3 || command -v python || true)
if [ -z "$PYTHON" ]; then
  echo "ERROR: Python (python3 or python) not found in PATH. Activate your conda env or install Python."
  exit 1
fi

# Verify the dexseq_count.py exists and is executable (some installs are not executable; we still can call with interpreter)
if [ ! -f "$DEXSEQ_COUNT_SCRIPT" ]; then
  echo "ERROR: dexseq_count.py not found at '$DEXSEQ_COUNT_SCRIPT'."
  echo "  -> Tip: if you installed DEXSeq in an R environment, look for 'python_scripts/dexseq_count.py' under your R library path."
  echo "  -> Example (inside R): system.file('python_scripts/dexseq_count.py', package='DEXSeq')"
  exit 1
fi

# Check flattened GTF exists — dexseq_count.py requires the flattened annotation
if [ ! -f "$FLATTENED_GTF" ]; then
  echo "ERROR: Flattened GTF '$FLATTENED_GTF' not found. Run dexseq_prepare_annotation.py first."
  exit 1
fi

# Create output directory (mkdir -p is idempotent and safe)
mkdir -p "$OUTPUT_DIR" || { echo "ERROR: Could not create output directory '$OUTPUT_DIR'"; exit 1; }

# Ensure BAM folder exists
if [ ! -d "$BAM_FOLDER" ]; then
  echo "ERROR: BAM folder '$BAM_FOLDER' does not exist or is not a directory."
  exit 1
fi

# -------------------------------------------------------------------------------

# Use find with -print0 and a while read -d '' loop to safely handle filenames with spaces/newlines.
# This avoids word-splitting issues with for-loops over plain command substitution.
echo "Searching for BAM files (recursively) in: $BAM_FOLDER"
found_any=0
# We'll populate an array BAM_FILES for robust iteration later (Bash arrays handle whitespace safely).
declare -a BAM_FILES=()
while IFS= read -r -d '' f; do
  BAM_FILES+=("$f")
  found_any=1
done < <(find "$BAM_FOLDER" -type f -name "*.bam" -print0 | sort -z)

if [ "$found_any" -eq 0 ]; then
  echo "ERROR: No BAM files found in '$BAM_FOLDER'. Check the folder and pattern (*.bam)."
  exit 1
fi

echo "Found ${#BAM_FILES[@]} BAM files."
echo "Using flattened GTF: $FLATTENED_GTF"
echo "Output directory: $OUTPUT_DIR"
echo "Stranded: $STRANDED"
echo "dexseq_count.py: $DEXSEQ_COUNT_SCRIPT"
echo "Python interpreter: $PYTHON"

# -------------------- OPTIONAL PRE-CHECKS (recommended) ------------------------
# These checks are not mandatory but strongly recommended to detect common problems.

# 1) Ensure BAMs appear coordinate-sorted. dexseq_count.py expects coordinate-sorted BAMs.
#    We check the @HD SO: coordinate header tag if present; this is a light-weight heuristic.
echo "Performing lightweight BAM header checks (sorted? indexed?) ..."
bad_sort=0
for BAM in "${BAM_FILES[@]}"; do
  # Extract the "SO" value from SAM header if present. Use samtools view -H if available.
  if command -v samtools >/dev/null 2>&1; then
    so=$(samtools view -H "$BAM" 2>/dev/null | awk -F'\t' '/\@HD/ { for(i=1;i<=NF;i++) if($i ~ /^SO:/) {split($i,a,":"); print a[2]; exit}}')
    if [ -n "$so" ] && [ "$so" != "coordinate" ]; then
      echo "  WARNING: $BAM header SO=$so (expected 'coordinate'). Consider running 'samtools sort' to coordinate-sort the BAM."
      bad_sort=1
    fi
    # Check for index (.bai) adjacent to BAM; dexseq_count.py does not strictly require a .bai, but many tools expect it.
    if [ ! -f "${BAM}.bai" ] && [ ! -f "${BAM%.bam}.bai" ]; then
      echo "  NOTE: index not found for $BAM (no .bai). For some downstream workflows indexing improves random-access performance."
    fi
  else
    # samtools not found — skip header check but warn
    echo "  NOTE: samtools not found; skipping BAM header checks. Install samtools to enable header/index checks."
    break
  fi
done
if [ "$bad_sort" -eq 1 ]; then
  echo "  Recommended: ensure BAMs are coordinate-sorted. Continue at your own risk."
fi

# 2) Suggest activating conda env where DEXSeq is installed (uncomment if you use Conda)
# echo "If DEXSeq installed inside a conda env, consider: source activate <env-name>"

# -------------------------------------------------------------------------------

# -------------------- ITERATE & COUNT -----------------------------------------
# We iterate over BAM_FILES and call dexseq_count.py. Each call produces one output file
# named <BAM_BASENAME>_dexseq_counts.txt in OUTPUT_DIR.
#
# We use "$PYTHON" to explicitly run with python3 if available.
# We pass "-f bam -r pos -s STRANDED" which are common flags:
#   -f bam : input format is BAM
#   -r pos : read order mode — "pos" recommended for coordinate-sorted BAMs
#   -s     : strandedness ("yes", "no", "reverse")
#
# Caveats & tips:
#  - For paired-end libraries, dexseq_count.py accepts additional options to handle paired reads.
#    See 'dexseq_count.py -h' for details (e.g., --paired).
#  - If dexseq_count.py is not executable, calling it as: python /path/to/dexseq_count.py ... works.
# -------------------------------------------------------------------------------

# Example: parallelization note
# If you want to run multiple counts in parallel (e.g., on a multi-core machine),
# consider running this loop through GNU parallel or a background-job pattern.
# For safety, this script runs sequentially to preserve log ordering.

for BAM in "${BAM_FILES[@]}"; do
  # Use basename and remove .bam suffix robustly
  BAM_BASE=$(basename "$BAM")
  # Remove final .bam extension (case-insensitive)
  BAM_STEM="${BAM_BASE%.[bB][aA][mM]}"  # safe removal of .bam or .BAM
  OUTPUT_FILE="${OUTPUT_DIR}/${BAM_STEM}_dexseq_counts.txt"

  echo "----------------------------------------"
  echo "Counting reads for sample: $BAM_STEM"
  echo "  BAM: $BAM"
  echo "  OUT: $OUTPUT_FILE"

  # Run dexseq_count.py and capture stdout/stderr to per-sample log to help debugging
  SAMPLE_LOG="${OUTPUT_DIR}/${BAM_STEM}_dexseq_count.log"
  # Example invocation; adjust flags for your library type:
  "$PYTHON" "$DEXSEQ_COUNT_SCRIPT" -f bam -r pos -s "$STRANDED" "$FLATTENED_GTF" "$BAM" "$OUTPUT_FILE" \
    > "$SAMPLE_LOG" 2>&1

  RET=$?
  if [ $RET -eq 0 ]; then
    echo "  ✓ Success: counts saved to $OUTPUT_FILE"
    echo "    log: $SAMPLE_LOG"
  else
    echo "  ✗ ERROR: dexseq_count.py failed for $BAM (exit code $RET)."
    echo "    Inspect $SAMPLE_LOG for details (first lines):"
    head -n 40 "$SAMPLE_LOG" || true
    # Decide policy: continue with other samples or exit early. We continue but flag an overall error.
    overall_error=1
  fi
done

# -------------------- SUMMARY --------------------------------------------------
if [ -n "${overall_error:-}" ]; then
  echo
  echo "One or more samples failed. Inspect individual *_dexseq_count.log files in $OUTPUT_DIR"
  exit 2
else
  echo
  echo "DEXSeq counting complete. Check $OUTPUT_DIR for count files and per-sample logs."
  exit 0
fi

# -------------------- POST-RUN CHECKS & RECOMMENDATIONS -----------------------
# - Inspect one or two produced count files: they should be two-column text with exonID and integer counts.
# - Ensure the flattened GTF used here is the same assembly/annotation used when aligning reads.
# - Before running DEXSeq in R, you may want to remove trailing summary rows starting with '_' (some dexseq_count variants include them).
#   The R script in the earlier message sanitizes these lines.
# - To speed up large datasets: run counting in parallel using GNU parallel:
#     printf "%s\0" "${BAM_FILES[@]}" | parallel -0 -j 8 ./dexseq_count_bam_folder_annotated.sh --single-task {}
#   Or convert the loop to background jobs and throttle concurrency with a job counter.
# - Alternative tools: for exon-level summarization, you can use featureCounts (subread) with a bins GTF or custom scripts.
#
# -------------------- TROUBLESHOOTING ----------------------------------------
# Common error messages and tips:
#  - "ValueError: ..." from dexseq_count.py: often due to malformed GTF or unexpected BAM format.
#  - "No reads overlapping features": check chromosome naming (chr1 vs 1) mismatches between BAM and GTF.
#  - "Permission denied": ensure dexseq_count.py is readable and python interpreter can execute it.
#  - If dexseq_count.py prints lots of stderr about SAM flags, confirm library type: single vs paired, strandedness.
#
# -------------------- OPTIONAL IMPROVEMENTS ----------------------------------
# - Add CLI argument parsing (getopts) to allow specifying BAM_FOLDER, GTF, OUTPUT_DIR etc. at runtime.
# - Add a --dry-run to only list found BAMs without executing counting.
# - Add automatic detection of strandedness by sampling a BAM via RSeQC or by heuristics (non-trivial).
# - Add resource-aware parallelization (e.g., xargs -P or GNU parallel) and per-job resource limits.
# - Validate output: confirm each OUTPUT_FILE exists and contains expected header/format before moving on.
#
# End of annotated wrapper script.
