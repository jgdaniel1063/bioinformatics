#!/bin/bash
#
# flatten_gtf_for_dexseq_annotated.sh
#
# Heavily-annotated wrapper to run dexseq_prepare_annotation.py and produce the
# "flattened" GTF required by the DEXSeq counting pipeline.
#
# Purpose
# -------
# DEXSeq requires a special "flattened" GTF that collapses overlapping exons
# into disjoint exon bins so reads can be counted at the exon-bin level.
# This script:
#  - validates environment and inputs
#  - handles gzipped GTF input transparently
#  - invokes dexseq_prepare_annotation.py (from DEXSeq)
#  - performs basic post-checks on output
#
# Important notes & caveats
# ------------------------
# - dexseq_prepare_annotation.py is distributed with the DEXSeq R package
#   in the 'python_scripts' subdirectory of the installed package, or may be
#   available on PATH if you installed it into a conda env.
# - Historically some DEXSeq python scripts required Python 2; in modern
#   Bioconductor/DEXSeq installs they should work with Python 3. Confirm by
#   running dexseq_prepare_annotation.py --help.
# - The input GTF should match the assembly used to align reads, and should
#   contain standard attributes like gene_id. If GTF uses nonstandard fields,
#   dexseq_prepare_annotation.py may fail.
# - It's best if input GTF is coordinate-sorted by chromosome and start.
#
# Recommended workflow:
#  1. Prepare / obtain a validated annotation GTF for your assembly (e.g. GRCz11).
#  2. Run this script to make the flattened GTF.
#  3. Use dexseq_count.py on BAMs with the flattened GTF to produce counts for R/DEXSeq.

set -euo pipefail
IFS=$'\n\t'

# -------------------- CONFIG (edit if needed) --------------------
# Paths below are hard-coded for convenience. You can replace these with
# command-line arguments or env vars if you prefer.
GTF_FILE="/home/jgd/Documents/bioinformatics_working/ref/danio_rerio.grcz11.115.sorted.gtf.gz"   # Input GTF (can be gzipped)
FLATTENED_GTF="/home/jgd/Documents/bioinformatics_working/ref/danio_rerio.grcz11.115.dexseq.gtf" # Output flattened GTF
DEXSEQ_SCRIPT="/home/jgd/miniconda3/bin/dexseq_prepare_annotation.py"                           # Path to dexseq_prepare_annotation.py
# -----------------------------------------------------------------

# -------------------- SAFETY: helper functions --------------------
log() { printf '[%s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*"; }

err() {
  echo "ERROR: $*" >&2
  exit 1
}

usage() {
  cat <<EOF
Usage: $(basename "$0") [--gtf /path/to/annotations.gtf[.gz]] [--out /path/to/flattened.gtf] [--dexseq /path/to/dexseq_prepare_annotation.py]
Options:
  --gtf       Input GTF (can be gzipped; default: $GTF_FILE)
  --out       Output flattened GTF (default: $FLATTENED_GTF)
  --dexseq    Path to dexseq_prepare_annotation.py (default: $DEXSEQ_SCRIPT)
  -h, --help  Show this help and exit
EOF
}

# Simple CLI parsing so the file is editable but also usable on the command line.
while [ $# -gt 0 ]; do
  case "$1" in
    --gtf) GTF_FILE="$2"; shift 2 ;;
    --out) FLATTENED_GTF="$2"; shift 2 ;;
    --dexseq) DEXSEQ_SCRIPT="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) err "Unknown option: $1";;
  esac
done

# Expand to absolute paths (best-effort)
GTF_FILE="$(realpath -s "$GTF_FILE")"
FLATTENED_GTF="$(realpath -s "$FLATTENED_GTF")"
DEXSEQ_SCRIPT="$(realpath -s "$DEXSEQ_SCRIPT" 2>/dev/null || true)"

log "Input GTF : $GTF_FILE"
log "Flattened : $FLATTENED_GTF"
log "DEXSeq script (user-provided): ${DEXSEQ_SCRIPT:-<not-set>}"

# -------------------- ENVIRONMENT CHECKS --------------------
# 1) Check input exists
if [ ! -f "$GTF_FILE" ]; then
  err "Input GTF file '$GTF_FILE' not found."
fi

# 2) Locate dexseq_prepare_annotation.py
#    Accept either an explicit file path or a script available on PATH.
if [ ! -x "$DEXSEQ_SCRIPT" ]; then
  if command -v dexseq_prepare_annotation.py >/dev/null 2>&1; then
    DEXSEQ_SCRIPT="$(command -v dexseq_prepare_annotation.py)"
    log "Using dexseq_prepare_annotation.py from PATH: $DEXSEQ_SCRIPT"
  else
    # If the provided path exists but is not executable, we'll attempt to call it with the Python interpreter.
    if [ -f "$DEXSEQ_SCRIPT" ]; then
      log "dexseq_prepare_annotation.py found at $DEXSEQ_SCRIPT (not marked executable) — will call with Python interpreter."
    else
      # Print actionable advice on how to find the script in an R installation
      log "dexseq_prepare_annotation.py not found at provided path and not on PATH."
      log "If you installed DEXSeq in R, try running this inside R to find it:"
      log "  system.file('python_scripts/dexseq_prepare_annotation.py', package='DEXSeq')"
      err "dexseq_prepare_annotation.py not found. Install DEXSeq or set --dexseq to the script path."
    fi
  fi
fi

# 3) Prefer python3 if available (many modern systems)
PYTHON=$(command -v python3 || command -v python || true)
if [ -z "$PYTHON" ]; then
  err "Python not found in PATH. Activate environment with Python or install Python."
fi
log "Using Python interpreter: $PYTHON"

# -------------------- PREPARE TEMP DIR --------------------
# Use a temporary directory to hold intermediate uncompressed GTF if needed.
TMPDIR=$(mktemp -d -t dexseq_flatten_XXXXXXXX) || err "Could not create temp dir"
cleanup() { rm -rf "$TMPDIR"; }
trap cleanup EXIT

# -------------------- HANDLE GZIPPED INPUT --------------------
# dexseq_prepare_annotation.py can typically accept gzipped input, but some
# installations or versions may not. We support both: if input is gzipped we
# will attempt to call the script directly. If that fails, we uncompress to a temp file.
is_gzipped=0
if [[ "$GTF_FILE" == *.gz ]]; then
  is_gzipped=1
fi

# Try directly invoking the script with the given GTF. If dexseq_prepare fails due
# to compressed input not supported, we'll try with an uncompressed temporary copy.
call_dexseq() {
  # arguments: input_gtf output_gtf
  local in_gtf="$1"
  local out_gtf="$2"
  # We call via explicit python to avoid relying on executable bit or shebang issues.
  log "Running dexseq_prepare_annotation.py with input: $in_gtf"
  "$PYTHON" "$DEXSEQ_SCRIPT" "$in_gtf" "$out_gtf"
}

# If output file already exists, offer to backup or exit.
if [ -f "$FLATTENED_GTF" ]; then
  log "Output flattened GTF already exists at $FLATTENED_GTF"
  backup="${FLATTENED_GTF}.$(date +%Y%m%d_%H%M%S).bak"
  log "Backing up existing file to: $backup"
  mv -v "$FLATTENED_GTF" "$backup"
fi

# First attempt: call the script as-is (lets the script handle gzipped input if supported)
if call_dexseq "$GTF_FILE" "$FLATTENED_GTF"; then
  log "dexseq_prepare_annotation.py completed successfully (first attempt)."
else
  # If the call failed and input was gzipped, try uncompressing to a temp file and retry
  if [ "$is_gzipped" -eq 1 ]; then
    log "First attempt failed. Trying an uncompressed copy of the GTF (some dexseq versions don't accept .gz)."
    TMP_GTF="$TMPDIR/$(basename "${GTF_FILE%.gz}")"
    log "Uncompressing to: $TMP_GTF"
    if command -v zcat >/dev/null 2>&1; then
      zcat "$GTF_FILE" > "$TMP_GTF"
    else
      gzip -dc "$GTF_FILE" > "$TMP_GTF"
    fi
    if call_dexseq "$TMP_GTF" "$FLATTENED_GTF"; then
      log "dexseq_prepare_annotation.py completed successfully (after uncompressing)."
    else
      err "dexseq_prepare_annotation.py failed even after uncompressing input. See stderr above for details."
    fi
  else
    err "dexseq_prepare_annotation.py failed. If GTF is gzipped, try uncompressing and rerunning."
  fi
fi

# -------------------- POST-RUN VALIDATION --------------------
# Basic checks on resulting flattened GTF to ensure it's plausible:
if [ ! -s "$FLATTENED_GTF" ]; then
  err "Flattened GTF was produced but is empty: $FLATTENED_GTF"
fi

# Check for expected 'exon' or 'exon_id' patterns in the flattened GTF:
# DEXSeq flattened GTF commonly contains lines with 'exon' features and attributes like 'exon_id' or 'group_id'.
grep -m 5 -E -i 'exon|exon_id|group_id' "$FLATTENED_GTF" >/dev/null 2>&1 || {
  log "Warning: Could not find expected exon-related lines in flattened GTF. Inspect the file:"
  log "  head -n 20 $FLATTENED_GTF"
}

# Inform user of success and location
log "Flattening complete. Output: $FLATTENED_GTF"

# -------------------- RECOMMENDATIONS & TROUBLESHOOTING --------------------
# Common problems and fixes:
#
# 1) Script not found:
#    - If dexseq_prepare_annotation.py isn't found, locate it inside your R library:
#      In R: system.file('python_scripts/dexseq_prepare_annotation.py', package='DEXSeq')
#    - If using conda, activate the environment with DEXSeq installed (e.g. conda activate <env>).
#
# 2) Permission / Python version issues:
#    - If you see syntax errors, ensure you're running an appropriate Python version.
#    - Historically some DEXSeq scripts used Python 2 constructs. Most modern versions support Python 3.
#
# 3) GTF formatting errors:
#    - Check that GTF contains expected attributes like gene_id. Example:
#        head -n 20 annotations.gtf | awk '{print $3; print $9; exit}'
#    - If the file is not sorted by chromosome/position, you can sort with:
#        zcat file.gtf.gz | sort -k1,1 -k4,4n | gzip -c > file.sorted.gtf.gz
#
# 4) Chromosome naming mismatch:
#    - If your BAM files use "chr1" but your GTF uses "1", you'll see no overlaps later.
#    - Normalize names early (e.g., with awk) if necessary.
#
# 5) Large files / runtime:
#    - Flattening large whole-genome GTFs can take several minutes to hours depending on CPU and memory.
#    - Run on a machine with sufficient RAM; monitor progress (script may print status to stdout).
#
# -------------------- OPTIONAL IMPROVEMENTS (suggestions) --------------------
# - Add an interactive --yes flag or --force to control automatic backup / overwrite behavior.
# - Provide a --dry-run to print commands without executing.
# - Add logging to a file in addition to stdout for reproducibility:
#       LOGFILE="/path/to/flatten_gtf.log"; exec > >(tee -a "$LOGFILE") 2>&1
# - Add validation against a gene list (ensure flattened bins exist for genes of interest).
# - Containerize (Docker/Singularity) to ensure consistent DEXSeq/Python versions.
#
# Done.
