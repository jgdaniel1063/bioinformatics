#!/usr/bin/env bash
#
# build_fresh_star_and_irfinder_no_conda_annotated.sh
#
# Heavily annotated script to build a fresh STAR genome index and an IRFinder
# reference, without activating any conda environment (assumes binaries are on PATH
# or that you set full paths to STAR/SAMTOOLS/IRFinder).
#
# Purpose
# -------
# - Create an uncompressed FASTA + GTF working copy in a temporary directory
# - Create a samtools FASTA index (.fai)
# - Build a STAR genome index (genomeGenerate) suitable for mapping RNA-seq reads
# - Use the STAR index as input to IRFinder's reference-building helper, falling
#   back to IRFinder's direct BuildRefProcess if the STAR-based step fails
#
# Why this exists (context)
# ------------------------
# Some downstream tools (mapping, intron retention detection) require properly
# prepared genome and annotation indices. STAR is a common aligner; IRFinder
# provides a tool to generate a reference tailored to intron retention detection.
#
# Key assumptions
# ---------------
# - You have enough disk space for an uncompressed genome FASTA and indexes.
# - STAR, samtools, and IRFinder binaries are installed and callable (or set as
#   explicit absolute paths in the variables below).
# - The GTF is coordinate-sorted (recommended) and matches the genome FASTA contig
#   names (chr1 vs 1 mismatches are a common source of errors).
#
# Reproducibility
# ---------------
# Record exact versions of the tools used (samtools, STAR, IRFinder) and checksums
# of the input FASTA/GTF. Consider running inside a pinned conda environment or
# container image for reproducibility.
#
# Safety
# ------
# - This script deletes (rm -rf) the target index directories before building them.
#   Do not point STAR_INDEX_DIR or IRFINDER_REF_DIR at important existing data unless
#   you intend to replace it.
# - Temporary working files are placed under TMPDIR and removed upon success.
#
# Invocation
# ----------
# Paste and run:
#   bash build_fresh_star_and_irfinder_no_conda_annotated.sh
#
# Edit the top-of-file variables (GENOME_FASTA_GZ, ANNOTATION_GTF_GZ, STAR_INDEX_DIR, IRFINDER_REF_DIR)
# before running if the placeholders are not correct.
#
# -----------------------------------------------------------------------------

set -euo pipefail               # strict: exit on error, undefined var is error, pipeline returns failure
# It's useful to see which command failed. The trap at the end could be added for even more context.

# ----------------------- CONFIGURE / EDIT IF NEEDED -----------------------
# Hard-coded placeholders (edit these to point to your files before running)
GENOME_FASTA_GZ="/home/jgd/Documents/2026-04-01_proc-enu_bioinformatics_results/ref/danio_rerio.grcz11.dna.primary_assembly.fa.gz"
ANNOTATION_GTF_GZ="/home/jgd/Documents/2026-04-01_proc-enu_bioinformatics_results/ref/danio_rerio.grcz11.115.sorted.gtf.gz"

# Where STAR should write its index directory (will be removed & re-created)
STAR_INDEX_DIR="/home/jgd/Documents/2026-04-01_proc-enu_bioinformatics_results/ref/star"

# Where IRFinder reference will be created (will be removed & re-created)
IRFINDER_REF_DIR="/home/jgd/Documents/2026-04-01_proc-enu_bioinformatics_results/ref/irfinder"

# Build parameters - tune for your reads and machine
READ_LENGTH=150            # approximate maximum read length in your data; STAR sjdbOverhang = READ_LENGTH - 1
THREADS=12                 # threads used by STAR and some other steps

# Binaries - if blank, the script tries to find them on PATH; you can set full paths
IRFINDER_BIN="/home/jgd/miniconda3/envs/irfinder_env/bin/IRFinder"   # set to IRFinder executable if not on PATH
STAR_BIN=""    # leave blank to auto-detect STAR or STAR-avx2 in PATH
SAMTOOLS_BIN="samtools"     # path or command name for samtools

# Temporary directory and logging locations
TMPDIR_BASE="/home/jgd/Documents/bioinformatics_working/tmp"   # base scratch dir - create if needed
TMPDIR="$(mktemp -d -p "${TMPDIR_BASE:-/tmp}" irbuild_XXXX)"
LOG="${TMPDIR}/build_all_indices.log"
# -------------------------------------------------------------------------

# Print basic provenance
echo "Log: ${LOG}"
echo "Temp dir: ${TMPDIR}"
echo

# -------------------------
# Detect STAR executable
# -------------------------
# STAR has variants: STAR, STAR-avx2 (some systems install avx2-optimized binary).
# If the user didn't hard-code STAR_BIN, try to find a suitable binary on PATH.
if [ -z "${STAR_BIN}" ]; then
  if command -v STAR >/dev/null 2>&1; then
    STAR_BIN="$(command -v STAR)"
  elif command -v STAR-avx2 >/dev/null 2>&1; then
    STAR_BIN="$(command -v STAR-avx2)"
  else
    echo "ERROR: STAR not found on PATH. Please install STAR or set STAR_BIN to full path." | tee "${LOG}"
    exit 1
  fi
fi

# -------------------------
# Verify samtools
# -------------------------
if ! command -v "${SAMTOOLS_BIN}" >/dev/null 2>&1; then
  echo "ERROR: samtools not found on PATH. Please install samtools or set SAMTOOLS_BIN to full path." | tee -a "${LOG}"
  exit 1
fi

# -------------------------
# Verify IRFinder binary exists and is executable
# -------------------------
if [ ! -x "${IRFINDER_BIN}" ]; then
  echo "ERROR: IRFinder not executable at ${IRFINDER_BIN}. Please ensure the path is correct and the file is executable." | tee -a "${LOG}"
  exit 1
fi

# -------------------------
# Report configuration
# -------------------------
echo "Using STAR: ${STAR_BIN}" | tee "${LOG}"
echo "Using samtools: $(command -v ${SAMTOOLS_BIN})" | tee -a "${LOG}"
echo "Using IRFinder: ${IRFINDER_BIN}" | tee -a "${LOG}"
echo "READ_LENGTH=${READ_LENGTH}, THREADS=${THREADS}" | tee -a "${LOG}"
echo "----" | tee -a "${LOG}"

# -------------------------
# Prepare uncompressed copies in TMPDIR
# -------------------------
# STAR & IRFinder prefer plain FASTA/GTF files. We copy/uncompress to TMPDIR to avoid modifying originals.
GENOME_FA="${TMPDIR}/genome.fa"
ANNOTATION_GTF="${TMPDIR}/annotation.gtf"

echo "Uncompressing FASTA and GTF to temporary files (this may take some time)..." | tee -a "${LOG}"
if [[ "${GENOME_FASTA_GZ}" == *.gz ]]; then
  # Using gunzip -c ensures we don't change the original file. Redirect stderr to log.
  gunzip -c "${GENOME_FASTA_GZ}" > "${GENOME_FA}" 2>> "${LOG}"
else
  cp -a "${GENOME_FASTA_GZ}" "${GENOME_FA}" 2>> "${LOG}"
fi

if [[ "${ANNOTATION_GTF_GZ}" == *.gz ]]; then
  gunzip -c "${ANNOTATION_GTF_GZ}" > "${ANNOTATION_GTF}" 2>> "${LOG}"
else
  cp -a "${ANNOTATION_GTF_GZ}" "${ANNOTATION_GTF}" 2>> "${LOG}"
fi

echo "Prepared files:" | tee -a "${LOG}"
ls -lh "${GENOME_FA}" "${ANNOTATION_GTF}" | tee -a "${LOG}"
echo "----" | tee -a "${LOG}"

# -------------------------
# Create FASTA index (.fai)
# -------------------------
echo "Creating samtools faidx..." | tee -a "${LOG}"
"${SAMTOOLS_BIN}" faidx "${GENOME_FA}" 2>> "${LOG}"
echo "faidx: ${GENOME_FA}.fai" | tee -a "${LOG}"
echo "----" | tee -a "${LOG}"

# -------------------------
# STAR genomeGenerate
# -------------------------
# Remove any existing index dir then run STAR genomeGenerate into it.
# Warning: STAR index can be large; ensure enough disk space (many GB).
echo "Removing existing STAR index dir (if any) and creating fresh index at ${STAR_INDEX_DIR}..." | tee -a "${LOG}"
rm -rf "${STAR_INDEX_DIR}"
mkdir -p "${STAR_INDEX_DIR}"

# Important STAR parameters to consider:
#  --sjdbOverhang: typically READ_LENGTH - 1 (captures maximum splice junction length for your reads)
#  --genomeSAindexNbases: default tuned for large genomes; can be lowered for small genomes (e.g., bacteria)
#  --runThreadN: number of threads
#
# For zebrafish (moderate genome size) defaults are usually fine. For very long read lengths
# adjust sjdbOverhang accordingly. Too large an sjdbOverhang increases memory usage.
echo "Running STAR genomeGenerate (sjdbOverhang=$((READ_LENGTH-1))) ..." | tee -a "${LOG}"
"${STAR_BIN}" --runThreadN "${THREADS}" \
  --runMode genomeGenerate \
  --genomeDir "${STAR_INDEX_DIR}" \
  --genomeFastaFiles "${GENOME_FA}" \
  --sjdbGTFfile "${ANNOTATION_GTF}" \
  --sjdbOverhang $((READ_LENGTH-1)) >> "${LOG}" 2>&1

echo "STAR genomeGenerate finished. Index at: ${STAR_INDEX_DIR}" | tee -a "${LOG}"
echo "----" | tee -a "${LOG}"

# -------------------------
# Build IRFinder reference
# -------------------------
# Preferred approach: IRFinder can import information from an existing STAR index
# using BuildRefFromSTARRef (faster if STAR index is present).
# If that fails (IRFinder version mismatch, different expectations), fallback to
# BuildRefProcess from raw FASTA+GTF (slower but usually works).
echo "Removing any existing IRFinder ref dir and building IRFinder from STAR reference at ${IRFINDER_REF_DIR}..." | tee -a "${LOG}"
rm -rf "${IRFINDER_REF_DIR}"

echo "Running IRFinder BuildRefFromSTARRef (this may take a while)..." | tee -a "${LOG}"
set +e
"${IRFINDER_BIN}" -m BuildRefFromSTARRef -r "${IRFINDER_REF_DIR}" -x "${STAR_INDEX_DIR}" >> "${LOG}" 2>&1
RC=$?
set -e

if [ ${RC} -eq 0 ]; then
  echo "IRFinder BuildRefFromSTARRef succeeded. Reference created at: ${IRFINDER_REF_DIR}" | tee -a "${LOG}"
  echo "Build log saved at: ${LOG}" | tee -a "${LOG}"
  # Cleanup temporary files
  rm -rf "${TMPDIR}"
  exit 0
fi

# If we reach here, the STAR-based IRFinder build failed.
echo "IRFinder BuildRefFromSTARRef FAILED (exit ${RC}). See log ${LOG}." | tee -a "${LOG}"
echo "Attempting fallback: BuildRefProcess from local FASTA+GTF (positional args)..." | tee -a "${LOG}"

# Fallback: BuildRefProcess. The exact IRFinder subcommand names/options vary by IRFinder versions,
# so this simpler, positional form is often more robust across versions.
set +e
"${IRFINDER_BIN}" -m BuildRefProcess -r "${IRFINDER_REF_DIR}" "${GENOME_FA}" "${ANNOTATION_GTF}" >> "${LOG}" 2>&1
RC2=$?
set -e

if [ ${RC2} -eq 0 ]; then
  echo "IRFinder BuildRefProcess succeeded. Reference created at: ${IRFINDER_REF_DIR}" | tee -a "${LOG}"
  rm -rf "${TMPDIR}"
  exit 0
fi

# Both IRFinder attempts failed — keep temp for debugging
echo "IRFinder BuildRefProcess also FAILED (exit ${RC2}). Keeping temporary files for debugging: ${TMPDIR}" | tee -a "${LOG"
echo "Last 200 lines of log (${LOG}):" | tee -a "${LOG}"
tail -n 200 "${LOG}" | tee -a "${LOG}"
echo
echo "Please inspect the log and the temporary files under ${TMPDIR}. If you want help diagnosing,"
echo " paste the tail of the log into the chat or open an issue with the IRFinder/STAR versions used."
exit 1

# -------------------------
# End of script
# -------------------------

# -------------------------
# Detailed troubleshooting & NOTES (non-executable)
# -------------------------
: <<'DOC_NOTES'

If something goes wrong, here's a checklist and explanation for common failure modes:

1) STAR not found or incompatible binary:
   - STAR avx2/avx512 binaries are optimized for specific CPU instruction sets. Running
     the wrong binary on an unsupported CPU may fail. Use the generic STAR if unsure.
   - Check version: STAR --version
   - If you need a specific STAR build, install via Bioconda: mamba install -c bioconda star

2) GTF/FASTA contig name mismatches:
   - A frequent source of errors is that the GTF uses 'chr1' while the FASTA uses '1' (or vice versa).
     STAR will still build an index, but IRFinder or other tools may later fail to fetch contigs.
   - Inspect the first few headers:
       zcat file.fa.gz | head -n 5
       zcat file.gtf.gz | grep -m 20 '^' | head -n 10
   - Normalize names beforehand (awk/perl) or create a mapping.

3) Memory / Disk usage:
   - STAR index building can require tens of GB depending on genome size and options.
   - Ensure TMPDIR is on a disk with sufficient space. You can change TMPDIR_BASE to a large local disk.

4) IRFinder failures:
   - IRFinder's BuildRefFromSTARRef expects a particular layout/version of STAR index; mismatches may occur.
   - Consult IRFinder documentation for version-specific instructions.
   - For many IRFinder versions the BuildRefProcess fallback will work when given raw FASTA and GTF.

5) Permission issues:
   - If building indexes as a non-root user, ensure you have write permission to STAR_INDEX_DIR and IRFINDER_REF_DIR.
   - The script removes existing directories with 'rm -rf' — be cautious if pointing at shared resources.

6) Reproducibility:
   - Record:
       samtools --version
       STAR --version
       IRFinder --version
     and the md5sum of your input FASTA/GTF:
       md5sum genome.fa.gz annotation.gtf.gz
   - Consider saving these into a small manifest inside the index directories.

7) Try a dry run / interactive debugging:
   - Run STAR genomeGenerate with fewer threads to test, e.g., --runThreadN 2
   - Inspect ${LOG} for STAR stderr output (it is verbose and contains resource usage hints).

8) If you want me to:
   - Convert this to an arg-parsing script with safer path checks and optional resume behavior.
   - Add an option to keep or remove temporary FASTA/GTF copies on success/failure.
   - Implement name-normalization between FASTA and GTF (add 'chr' or remove it).
   - Add a small "sanity check" BAM mapping after build (map a few reads to confirm index validity).
DOC_NOTES
