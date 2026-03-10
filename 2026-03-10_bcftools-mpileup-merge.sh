#!/usr/bin/env bash
#
# merge_vcfs_local_annotated.sh
#
# Purpose (local / non-cluster)
# -----------------------------
# - Recursively find VCF/VCF.GZ/BCF files under INPUT_DIR
# - Ensure each file is bgzipped and indexed (tabix for .vcf.gz; bcftools index for .bcf)
# - Merge all valid files with bcftools merge into a single compressed, indexed VCF
#
# Differences from the original:
# - All SLURM / module / cluster-specific directives removed
# - Added defensive checks, provenance reporting, and careful handling of filenames
# - Annotated heavily with guidance and alternatives
#
# Usage:
#   Edit INPUT_DIR, OUTPUT_DIR, OUTPUT_FILE variables below, then run:
#     bash merge_vcfs_local_annotated.sh
#
# Notes / recommendations:
# - Prefer using pre-computed, bgzipped+tabix-indexed VCFs. Recompressing large files is I/O-heavy
#   and can take a long time; ensure you have disk space and patience.
# - For reproducibility, record tool versions (this script writes a versions.txt file).
# - If merging many samples (hundreds), bcftools merge may need lots of memory; consider merging in batches.
#
set -euo pipefail
umask 002

# -------------------------
# User-editable configuration
# -------------------------
INPUT_DIR="/nfs/turbo/umms-jshavit/jgdaniel/06-03-2024_proc-enu_paper/11-19-2025_final_results/bcftools-mpileup_workflow/snp_calling"
OUTPUT_DIR="/nfs/turbo/umms-jshavit/jgdaniel/06-03-2024_proc-enu_paper/raw_output"
OUTPUT_FILE="${OUTPUT_DIR}/merged_proc-enu_samples.vcf.gz"

# Tools: assume available on PATH; override here if you want to point to specific installs
BCFTOOLS="${BCFTOOLS:-bcftools}"
TABIX="${TABIX:-tabix}"
BGZIP="${BGZIP:-bgzip}"     # bgzip is usually part of htslib
SAMTOOLS="${SAMTOOLS:-samtools}"

# Temporary working area for recompression (may hold large temporary files)
TMPDIR="${TMPDIR:-$(mktemp -d "${OUTPUT_DIR}/merge_tmp.XXXX")}"

# Minimal sanity thresholds
MIN_FILES_REQUIRED=1

# -------------------------
# Sanity checks: tools
# -------------------------
for tool in "$BCFTOOLS" "$TABIX" "$BGZIP" "$SAMTOOLS"; do
  if ! command -v "$tool" >/dev/null 2>&1; then
    echo "[ERROR] Required tool not found on PATH or by override: $tool" >&2
    echo "  Install it or set the variable at the top of this script to the full path." >&2
    exit 1
  fi
done

# Record tool versions for provenance
mkdir -p "$OUTPUT_DIR"
VERSIONS="${OUTPUT_DIR}/merge_tools_versions.txt"
{
  echo "Timestamp: $(date -u +"%Y-%m-%d %H:%M:%SZ")"
  "$BCFTOOLS" --version 2>&1 | head -n 1
  "$TABIX" --version 2>&1 || true
  "$BGZIP" --version 2>&1 || true
  "$SAMTOOLS" --version 2>&1 | head -n 1
} > "$VERSIONS"

echo "[INFO] Tool versions written to: $VERSIONS"

# -------------------------
# Find candidate input files (robust to spaces/newlines)
# -------------------------
mapfile -d '' INPUT_FILES < <(find "$INPUT_DIR" -type f \( -iname "*.vcf" -o -iname "*.vcf.gz" -o -iname "*.bcf" \) -print0 | sort -z)

NUM_FOUND=${#INPUT_FILES[@]}
if [[ $NUM_FOUND -lt $MIN_FILES_REQUIRED ]]; then
  echo "[ERROR] No VCF/VCF.GZ/BCF files found under: $INPUT_DIR" >&2
  exit 1
fi

echo "[INFO] Found ${NUM_FOUND} candidate input files (first 25 shown):"
for ((i=0;i<NUM_FOUND && i<25;i++)); do
  printf "  %s\n" "${INPUT_FILES[i]}"
done
if [[ $NUM_FOUND -gt 25 ]]; then
  echo "  ... and $((NUM_FOUND-25)) more"
fi

# -------------------------
# Helper: safe recompress a .vcf.gz if tabix fails
# -------------------------
# Recompressing a potentially non-BGZF gzip into proper bgzip can be I/O intensive.
# This function creates a temporary decompressed copy and bgzips it in TMPDIR.
recompress_to_bgzip() {
  local src="$1"
  local dst="$2"   # full path to final bgzipped file
  local tmp_unz="${TMPDIR}/$(basename "${src%.*}").tmp.vcf"
  echo "[INFO] Recompressing (gunzip -> bgzip) '$src' -> '$dst' (tmp: $tmp_unz)"
  set -o pipefail
  if gunzip -c "$src" > "$tmp_unz"; then
    if "$BGZIP" -c "$tmp_unz" > "${dst}"; then
      rm -f "$tmp_unz"
      return 0
    else
      rm -f "$tmp_unz"
      return 2
    fi
  else
    rm -f "$tmp_unz"
    return 3
  fi
  set +o pipefail
}

# -------------------------
# Normalize / validate input files: ensure bgzip + index present
# Build a list of VALID_FILES to pass to bcftools merge
# -------------------------
VALID_FILES=()
echo "[INFO] Validating and preparing input files (compress/index if needed)..."
for f in "${INPUT_FILES[@]}"; do
  skip=false
  original="$f"

  # For BCF, ensure index (.csi) exists (bcftools index creates .csi)
  if [[ "${f,,}" == *.bcf ]]; then
    if [[ ! -f "${f}.csi" && ! -f "${f}.csi" ]]; then
      echo "[INFO] Indexing BCF: $f"
      if ! "$BCFTOOLS" index "$f"; then
        echo "[WARN] Failed to index BCF: $f; skipping" >&2
        skip=true
      fi
    fi
    if [[ "$skip" = false ]]; then
      VALID_FILES+=("$f")
    fi

  # For compressed VCF (.vcf.gz)
  elif [[ "${f,,}" == *.vcf.gz ]]; then
    # Check for a tabix index (.tbi or .csi)
    if [[ -f "${f}.tbi" || -f "${f}.csi" ]]; then
      # good
      VALID_FILES+=("$f")
    else
      # Try to index directly; if tabix fails, try to recompress from the underlying gzip stream
      echo "[INFO] Index missing for $f; attempting tabix -p vcf ..."
      if "$TABIX" -p vcf "$f" 2>/dev/null; then
        VALID_FILES+=("$f")
      else
        echo "[WARN] tabix failed for $f; attempting to recompress to proper bgzip..."
        dst="${TMPDIR}/$(basename "$f")"
        if recompress_to_bgzip "$f" "$dst"; then
          if "$TABIX" -p vcf "$dst"; then
            # Move recompressed file into OUTPUT_DIR/input_recompressed/ to not overwrite originals
            mkdir -p "${OUTPUT_DIR}/recompressed_inputs"
            mv "$dst" "${OUTPUT_DIR}/recompressed_inputs/$(basename "$dst")"
            dst2="${OUTPUT_DIR}/recompressed_inputs/$(basename "$dst")"
            VALID_FILES+=("$dst2")
            echo "[INFO] Recompressed and indexed -> $dst2"
          else
            echo "[WARN] tabix failed on recompressed file $dst; skipping original $f" >&2
            rm -f "$dst"
          fi
        else
          echo "[WARN] Failed to recompress $f; skipping" >&2
        fi
      fi
    fi

  # For uncompressed VCF (.vcf)
  elif [[ "${f,,}" == *.vcf ]]; then
    # Compress with bgzip and index; write compressed file next to original with .vcf.gz
    compressed="${f}.gz"
    echo "[INFO] Compressing uncompressed VCF: $f -> $compressed (bgzip)"
    if "$BGZIP" -c "$f" > "$compressed"; then
      if "$TABIX" -p vcf "$compressed"; then
        VALID_FILES+=("$compressed")
      else
        echo "[WARN] tabix failed on newly compressed file $compressed; skipping" >&2
        rm -f "$compressed"
      fi
    else
      echo "[WARN] bgzip failed on $f; skipping" >&2
    fi

  else
    echo "[WARN] Unknown file extension, skipping: $f"
  fi
done

NUM_VALID=${#VALID_FILES[@]}
echo "[INFO] ${NUM_VALID} valid files ready for merge."

if [[ $NUM_VALID -lt 1 ]]; then
  echo "[ERROR] No valid input files to merge. Exiting." >&2
  # Clean up TMPDIR
  rm -rf "$TMPDIR"
  exit 1
fi

# -------------------------
# Merge files with bcftools merge
# -------------------------
echo "[INFO] Merging files into: $OUTPUT_FILE"
# Create parent dir
mkdir -p "$(dirname "$OUTPUT_FILE")"

# It is often preferable to run bcftools merge with a limited number of files at once
# for memory reasons; for simplicity we call merge on all files here.
# Output compressed VCF (-Oz) for convenience.
set -x
"$BCFTOOLS" merge -Oz -o "$OUTPUT_FILE" "${VALID_FILES[@]}"
merge_rc=$?
set +x

if [[ $merge_rc -ne 0 ]]; then
  echo "[ERROR] bcftools merge failed with exit code ${merge_rc}" >&2
  rm -rf "$TMPDIR"
  exit $merge_rc
fi

# Index the output
echo "[INFO] Indexing merged VCF: $OUTPUT_FILE"
"$TABIX" -p vcf "$OUTPUT_FILE"

echo "[OK] Merge completed successfully."
echo "Output: $OUTPUT_FILE"
echo "Provenance: $VERSIONS"
echo "Temporary working dir (removed): $TMPDIR"

# Cleanup
rm -rf "$TMPDIR"

# End of script
