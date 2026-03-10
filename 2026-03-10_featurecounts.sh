#!/usr/bin/env bash
#
# featurecounts_recursive_annotated.sh
#
# Heavily annotated wrapper around featureCounts to:
#  - discover BAMs recursively under a directory
#  - run comprehensive counts (gene-level, exon-level, junctions)
#  - optionally write BAMs annotated with read assignment tags
#  - capture summary files for downstream QC
#
# This annotated version:
#  - explains every major step and rationales for the flags used
#  - points out common pitfalls and alternatives
#  - shows safer ways to iterate over file lists (handles spaces in paths)
#  - suggests follow-up analyses and improvements
#
# Note: This script is intended for interactive use / single-node runs.
# For large-scale batch processing, convert to a per-sample job submitted to a scheduler
# or parallelize with GNU parallel.
#
# Requirements:
#  - featureCounts (Subread package)
#  - samtools (for BAM indexing and BAM outputs)
#  - bash 4+ recommended for arrays
#
# Usage:
#  - Edit the HARD-CODED VALUES section below to point to your BAM folder and GTF.
#  - Run: bash featurecounts_recursive_annotated.sh
#
# -----------------------------------------------------------------------------

set -euo pipefail
shopt -s nullglob

log() { printf '[%s] %s\n' "$(date +%F\ %T)" "$*"; }

# ---------------------------------------------------------------------------
# USER-EDITABLE VALUES (change these for your environment)
# ---------------------------------------------------------------------------
# Directory to search (recursively) for BAM files.
BAM_FOLDER="/home/jgd/Documents/bioinformatics_working/2025-11-06_enu-proc/align/star-ensembl"

# Path to GTF (can be gzipped). Prefer a coordinate-sorted GTF.
GTF_FILE="/home/jgd/Documents/bioinformatics_working/ref/danio_rerio.grcz11.115.sorted.gtf.gz"

# Prefix used for output files. Script will append descriptive suffixes.
OUTPUT_PREFIX="featurecounts"

# Number of threads to give featureCounts. Choose according to CPU core availability.
THREADS=12

# Strandedness: 0 = unstranded, 1 = stranded (reads from the same strand), 2 = reverse stranded.
# IMPORTANT: Set this to the correct value for your library prep; using the wrong value will mis-assign reads.
STRANDED=0

# If your reads are single-end set PAIRED=0, otherwise 1 for paired-end.
PAIRED=1
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# SANITY / TOOL CHECKS
# ---------------------------------------------------------------------------
# Fail early with clear message if tools aren't present.
if ! command -v featureCounts >/dev/null 2>&1; then
  echo "ERROR: featureCounts not found. Install Subread (conda install -c bioconda subread)" >&2
  exit 1
fi
if ! command -v samtools >/dev/null 2>&1; then
  echo "ERROR: samtools not found. Install samtools (conda install -c bioconda samtools)" >&2
  exit 1
fi

log "featureCounts: $(command -v featureCounts)"
log "samtools: $(command -v samtools)"

# Verify GTF exists
if [[ ! -f "$GTF_FILE" ]]; then
  echo "ERROR: GTF file '$GTF_FILE' not found." >&2
  exit 1
fi

# ---------------------------------------------------------------------------
# DISCOVER BAM FILES (robust to spaces/newlines)
# ---------------------------------------------------------------------------
# Original script used: find ... | tr '\n' ' ' which breaks with spaces.
# Here we build a bash array using NUL separators which handles arbitrary filenames.
mapfile -d '' BAM_FILES < <(find "$BAM_FOLDER" -type f -name "*.bam" -print0)

if [[ ${#BAM_FILES[@]} -eq 0 ]]; then
  echo "ERROR: No BAM files found in '$BAM_FOLDER'." >&2
  exit 1
fi

# Log how many BAMs found; don't print full list if it's long.
log "Found ${#BAM_FILES[@]} BAM file(s) under: $BAM_FOLDER"
# Example: print first few
for ((i=0;i<${#BAM_FILES[@]} && i<5;i++)); do
  log "  BAM[$i]=${BAM_FILES[i]}"
done
if [[ ${#BAM_FILES[@]} -gt 5 ]]; then
  log "  ... (and $(( ${#BAM_FILES[@]} - 5 )) more)"
fi

# ---------------------------------------------------------------------------
# COMMON OPTIONS: explanation & composition
# ---------------------------------------------------------------------------
# We create a set of featureCounts options and explain why each is included.
#
# Key flags used here:
#  - -p : paired-end mode. If your data are single-end, do NOT set -p.
#  - -B : count only fragments where both ends are properly mapped (paired-end only).
#  - -C : exclude chimeric fragments (that map to different chromosomes / chimeric).
#  - -M : count multi-mapping reads (reads with multiple alignments). Without -M
#         multi-mappers are ignored by default; with -M they're counted but flagged.
#  - -O : allow reads to be assigned to overlapping features (counts may then double-count).
#  - --fraction : when multi-mappers are counted (-M), distribute fractional counts.
#  - -s : strandedness. 0 = unstranded, 1 = stranded, 2 = reverse stranded.
#  - -t : feature type in the GTF (e.g., exon). featureCounts will group exons to genes
#  - -g : attribute to use as gene id (e.g., gene_id, gene_name). Use the attribute present in your GTF.
#  - -T : threads (for parallel counting)
#
# Advice:
#  - Choose -s carefully. Many RNA-seq preps are unstranded (0) or stranded (1 or 2).
#  - Use -p only for paired-end libraries.
#  - If you want conservative gene counts (one-count-per-fragment), use -p -B and avoid -O if genes overlap.
#  - If transcripts / exons overlap and you want to avoid double-counting, omit -O.
# ---------------------------------------------------------------------------

# Build the options string based on our variables
FC_OPTS=()
if [[ "$PAIRED" -eq 1 ]]; then
  FC_OPTS+=("-p" "-B" "-C")
fi
# Count multi-mapping reads and use fractional assignment
FC_OPTS+=("-M" "--fraction")
# Allow reads overlapping multiple features (useful for exon-level when features overlap)
FC_OPTS+=("-O")
# Strandedness
FC_OPTS+=("-s" "$STRANDED")
# Feature type and gene attribute
FC_OPTS+=("-t" "exon" "-g" "gene_id")
# Threads
FC_OPTS+=("-T" "$THREADS")
# Annotation file
FC_OPTS+=("-a" "$GTF_FILE")

log "featureCounts options: ${FC_OPTS[*]}"

# ---------------------------------------------------------------------------
# OUTPUT DIRECTORY & PREFIX
# ---------------------------------------------------------------------------
# By default OUTPUT_PREFIX is just a string. We create a directory to hold outputs
# and produce files with descriptive suffixes.
OUTPUT_DIR="$(pwd)/${OUTPUT_PREFIX}_results"
mkdir -p "$OUTPUT_DIR"
log "Outputs will be written to: $OUTPUT_DIR"

# We will write intermediate and final filenames here.
GENE_COUNTS="${OUTPUT_DIR}/${OUTPUT_PREFIX}_gene_counts.txt"
FEATURE_COUNTS="${OUTPUT_DIR}/${OUTPUT_PREFIX}_feature_counts.txt"
JUNCTION_COUNTS="${OUTPUT_DIR}/${OUTPUT_PREFIX}_junction_counts.txt"
ASSIGNED_BAM_DIR="${OUTPUT_DIR}/assigned_bams"
mkdir -p "$ASSIGNED_BAM_DIR"

# ---------------------------------------------------------------------------
# 1) GENE-LEVEL COUNTS
# ---------------------------------------------------------------------------
log "Running gene-level counts..."

# featureCounts accepts multiple BAMs as separate arguments. We pass the array safely.
# The -o option writes the main table. featureCounts also generates a *.summary file next to it.
featureCounts "${FC_OPTS[@]}" -o "$GENE_COUNTS" "${BAM_FILES[@]}"
fc_rc=$?
if [[ $fc_rc -ne 0 ]]; then
  echo "ERROR: Gene-level featureCounts failed (exit $fc_rc)" >&2
  exit $fc_rc
fi
log "Gene-level counts written to: $GENE_COUNTS"

# ---------------------------------------------------------------------------
# 2) FEATURE-LEVEL (EXON) COUNTS
# ---------------------------------------------------------------------------
log "Running feature-level (exon) counts..."

# -f tells featureCounts to produce counts at the feature level (exon).
featureCounts "${FC_OPTS[@]}" -f -o "$FEATURE_COUNTS" "${BAM_FILES[@]}"
fc_rc=$?
if [[ $fc_rc -ne 0 ]]; then
  echo "ERROR: Feature-level featureCounts failed (exit $fc_rc)" >&2
  exit $fc_rc
fi
log "Feature-level counts written to: $FEATURE_COUNTS"

# ---------------------------------------------------------------------------
# 3) JUNCTION COUNTS
# ---------------------------------------------------------------------------
# featureCounts supports "-J" to count reads spanning splice junctions when run with appropriate options.
# Note: the exact behavior depends on aligner flags and whether junctions are reported as features.
log "Running junction counts..."
featureCounts "${FC_OPTS[@]}" -J -o "$JUNCTION_COUNTS" "${BAM_FILES[@]}"
fc_rc=$?
if [[ $fc_rc -ne 0 ]]; then
  echo "ERROR: Junction featureCounts failed (exit $fc_rc)" >&2
  exit $fc_rc
fi
log "Junction counts written to: $JUNCTION_COUNTS"

# ---------------------------------------------------------------------------
# 4) PER-BAM ASSIGNED BAMs (BAMs WITH READ ASSIGNMENT TAGS)
# ---------------------------------------------------------------------------
# The -R BAM option makes featureCounts write a BAM with additional tags describing
# assignment (e.g., gene assignment information). This can be very useful for
# inspecting which reads were assigned and why.
#
# Important:
#  - featureCounts -R BAM writes output in the current directory by default; we direct outputs into ASSIGNED_BAM_DIR.
#  - This step can be slow (it reprocesses each BAM) and will write one BAM per input BAM.
#  - Some versions of featureCounts may write a SAM; we keep .bam extension and use samtools to index after.
# ---------------------------------------------------------------------------
log "Generating BAMs with read assignment tags (this may take a while)..."
for bam in "${BAM_FILES[@]}"; do
  # basename without .bam (handles filenames with multiple dots)
  bambase="$(basename "$bam")"
  bamstem="${bambase%.bam}"
  outbam="${ASSIGNED_BAM_DIR}/${bamstem}_assigned.bam"

  log "  Assigning reads for: $bam -> $outbam"
  # featureCounts -R BAM -o output -T threads input.bam
  # Note: featureCounts -o expects an output filename for the table; when used with -R BAM and a single input,
  # the assigned BAM is written to a file named like <input>.featureCounts.bam (behavior depends on version).
  # To avoid relying on that naming, we run featureCounts with a per-sample temporary working directory.
  tmpdir="$(mktemp -d "${ASSIGNED_BAM_DIR}/fc_tmp_${bamstem}.XXXX")"
  pushd "$tmpdir" >/dev/null
  # run featureCounts for this single BAM; -o name we set to tmp table
  featureCounts "${FC_OPTS[@]}" -R BAM -o "${bamstem}_assigned_counts.txt" "$bam"
  rc=$?
  if [[ $rc -ne 0 ]]; then
    echo "WARNING: featureCounts assignment failed for $bam (exit $rc). See $tmpdir for details." >&2
    popd >/dev/null
    rm -rf "$tmpdir"
    continue
  fi

  # featureCounts typically writes a BAM in the working dir with name like <input>.featureCounts.bam or similar.
  # Search for a BAM in the tmpdir and move it to outbam
  assigned_candidates=( "$tmpdir"/*.bam )
  if [[ ${#assigned_candidates[@]} -eq 0 ]]; then
    echo "WARNING: No assigned BAM found in $tmpdir for $bam. Inspect featureCounts output files." >&2
    popd >/dev/null
    rm -rf "$tmpdir"
    continue
  fi
  # Move the first candidate to the desired destination
  mv -v "${assigned_candidates[0]}" "$outbam"
  # Index the output BAM for convenience
  samtools index -@ "$THREADS" "$outbam"
  popd >/dev/null
  # Clean up temporary dir
  rm -rf "$tmpdir"
  log "  Wrote assigned BAM: $outbam"
done

# ---------------------------------------------------------------------------
# 5) PARSE/AGGREGATE SUMMARIES
# ---------------------------------------------------------------------------
# featureCounts writes a *.summary file beside each output table. We attempt to find and collect them.
log "Collecting featureCounts summary files..."

# locate summary files (search in OUTPUT_DIR)
shopt -s globstar
summary_files=()
while IFS= read -r -d '' f; do summary_files+=( "$f" ); done < <(find "$OUTPUT_DIR" -maxdepth 3 -type f -name "*.summary" -print0 2>/dev/null)

if [[ ${#summary_files[@]} -gt 0 ]]; then
  log "Found ${#summary_files[@]} summary file(s)."
  # Concatenate or parse per need. Example: print filenames and a short extract
  for sf in "${summary_files[@]}"; do
    log "  Summary: $sf"
    # show the header and the assignment line (if present)
    echo "----- $sf -----"
    head -n 20 "$sf"
    echo
  done
else
  log "No summary files found under $OUTPUT_DIR"
fi

# ---------------------------------------------------------------------------
# RECOMMENDED NEXT STEPS / NOTES
# ---------------------------------------------------------------------------
cat <<'EOF'

NOTES & RECOMMENDATIONS:
 - Strandness:
     Ensure STRANDED is set correctly for your library prep: 0 (unstranded), 1 (stranded), 2 (reverse).
     Using an incorrect value can dramatically change assignment numbers.

 - Paired-end:
     If your data are single-end, set PAIRED=0 and remove -p -B -C from FC_OPTS accordingly.
     The script currently uses PAIRED=1 by default.

 - Overlapping features (-O):
     If you set -O, featureCounts can assign the same read to multiple features (counts may not be exclusive).
     If you require exclusive counting per gene, omit -O, but then reads overlapping multiple features may be unassigned.

 - Multi-mappers (-M, --fraction):
     Counting multi-mapping reads can increase sensitivity but complicates downstream normalization.
     --fraction apportions fractional counts among mapping loci. If you prefer to ignore multimappers, remove -M.

 - Performance:
     For many BAMs, consider running featureCounts per-sample in parallel (GNU parallel) and then merging tables.
     featureCounts itself is multi-threaded for large input lists (-T), but per-sample parallelism can be more efficient
     when you have many samples and many CPUs.

 - Output formats & downstream:
     The gene-level output table can be read into R (edgeR/DESeq2) after minor reshaping.
     Consider producing a sample table (sample IDs and file paths) to pair with counts for reproducible analyses.

 - Quality control:
     Inspect the *.summary files for assignment fractions (assigned, ambiguous, unassigned).
     High unassigned fraction often indicates mismatched annotation (e.g., wrong GTF), wrong strandedness, or poor mapping.

 - Reproducibility:
     Save the versions of tools used (featureCounts, samtools) and the exact GTF used:
         featureCounts -v
         samtools --version
     Consider saving checksums for the GTF and a README in $OUTPUT_DIR documenting options.

 - Alternative tools:
     If you need transcript-level abundance (not gene counts), consider Salmon/Kallisto + tximport workflows.
     For junction-centric analysis, specialized tools (RegTools, MAJIQ) may provide richer junction metrics.

EOF

log "FeatureCounts run complete. Outputs in: $OUTPUT_DIR"
log "Gene counts: $GENE_COUNTS"
log "Feature counts: $FEATURE_COUNTS"
log "Junction counts: $JUNCTION_COUNTS"
log "Assigned BAMs (if generated): $ASSIGNED_BAM_DIR"

# End of script
