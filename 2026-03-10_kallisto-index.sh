#!/usr/bin/env bash
#
# build_kallisto_index_annotated.sh
#
# Heavily annotated script to build a kallisto index from a transcriptome FASTA.
#
# Summary:
#  - Input: transcriptome FASTA (transcripts, not the genome)
#  - Output: kallisto index file (.idx)
#  - The script performs basic checks, logs the run, and records provenance information.
#
# Why you need a transcriptome index:
#  - Kallisto performs pseudo-alignment of reads to transcripts. It requires an index
#    built over transcript sequences (FASTA). Do not use a genome FASTA here.
#
# Quick usage:
#  - Edit the TRANSCRIPTS_FA and KALLISTO_INDEX variables below then run:
#      bash build_kallisto_index_annotated.sh
#
# Reproducibility tips:
#  - Record the kallisto version and input file checksums (the script does this).
#  - Consider running inside a pinned conda environment or container so the same binary is used later.
#
# Notes on performance & parameters:
#  - The kallisto default k-mer length (k) is typically 31; this works for common read lengths (>=50bp).
#    For very short reads, you may want to choose a smaller k (kallisto index supports a --kmer-length / -k option).
#  - Kallisto index building is relatively fast for typical transcriptomes (minutes), but large transcriptomes
#    or network-mounted filesystems may be slower.
#
# Safety:
#  - Ensure the FASTA contains transcript sequences (mRNA) and not genomic sequences; otherwise quantification
#    results will be uninterpretable.
#
# -----------------------------------------------------------------------------

set -euo pipefail

########################################
# USER-EDITABLE CONFIGURATION
########################################

# Path to transcriptome FASTA (can be gzipped). MUST be transcript sequences (e.g., Ensembl/RefSeq transcripts).
TRANSCRIPTS_FA="/home/jgd/Documents/bioinformatics_working/reference/mrna.fa.gz"

# Destination index file path. Pick a meaningful name including the assembly and transcript source/version.
KALLISTO_INDEX="/home/jgd/Documents/bioinformatics_working/reference/danRer11.kallisto.idx"

# Directory for logs and small provenance files. The script creates it if necessary.
OUT_DIR="/home/jgd/Documents/bioinformatics_working/output/kallisto"

# Optional: change this to pass extra flags to kallisto index (e.g., "-k 23" to set kmer length)
# See `kallisto index --help` for available options on your installed version.
KALLISTO_EXTRA_OPTS="${KALLISTO_EXTRA_OPTS:-}"

########################################
# NO EDITS BELOW THIS LINE (unless you know what you're changing)
########################################

mkdir -p "${OUT_DIR}"

RUN_ID="$(date +%Y%m%d_%H%M%S)"
LOG="${OUT_DIR}/kallisto_index_${RUN_ID}.log"

# Send stdout/stderr to the log and to the terminal (tee)
exec > >(tee -a "${LOG}") 2>&1

echo "=============================================="
echo " K A L L I S T O   I N D E X   B U I L D"
echo "=============================================="
echo "[INFO] Start: $(date -u +"%Y-%m-%d %H:%M:%SZ")     (UTC)"
echo "[INFO] Transcript FASTA: ${TRANSCRIPTS_FA}"
echo "[INFO] Output index:     ${KALLISTO_INDEX}"
echo "[INFO] Log file:         ${LOG}"
echo

########################################
# SANITY CHECKS
########################################

echo "[INFO] Checking kallisto availability..."
if ! command -v kallisto >/dev/null 2>&1; then
    echo "[ERROR] kallisto not found in PATH. Install it first (e.g. conda/pip) and ensure it is on PATH." >&2
    echo "[HINT] conda example: mamba install -c bioconda kallisto" >&2
    exit 1
fi

# Capture kallisto version for provenance
KALLISTO_BIN="$(command -v kallisto)"
KALLISTO_VER="$("${KALLISTO_BIN}" version 2>&1 | head -n1 || true)"
echo "[INFO] kallisto binary: ${KALLISTO_BIN}"
echo "[INFO] kallisto version: ${KALLISTO_VER}"

echo "[INFO] Verifying input FASTA exists and is readable..."
if [[ ! -f "${TRANSCRIPTS_FA}" ]]; then
    echo "[ERROR] Transcript FASTA not found: ${TRANSCRIPTS_FA}" >&2
    exit 1
fi
if [[ ! -r "${TRANSCRIPTS_FA}" ]]; then
    echo "[ERROR] Transcript FASTA exists but is not readable: ${TRANSCRIPTS_FA}" >&2
    exit 1
fi

########################################
# PROVENANCE: checksum & size
########################################

echo "[INFO] Recording input FASTA checksum and size (may be slow for large files)..."
# If the FASTA is gzipped, md5sum will process the compressed bytes; that's still useful as a fingerprint.
# If you want checksum of the uncompressed content, compute md5sum after gunzip -c, but that is I/O heavy.
if command -v md5sum >/dev/null 2>&1; then
    md5sum "${TRANSCRIPTS_FA}" | awk '{print "[PROV] md5:", $1, $2}'
else
    echo "[PROV] md5sum not available; skipping checksum."
fi
# Report file size
stat --printf='[PROV] size: %s bytes %n\n' "${TRANSCRIPTS_FA}" || ls -lh "${TRANSCRIPTS_FA}"

########################################
# PREP: check destination
########################################

# If index already exists, warn and back it up (or abort). By default we overwrite.
if [[ -f "${KALLISTO_INDEX}" ]]; then
    echo "[WARN] Target index already exists: ${KALLISTO_INDEX}"
    BACKUP="${KALLISTO_INDEX}.bak.${RUN_ID}"
    echo "[INFO] Backing up existing index to: ${BACKUP}"
    mv -v "${KALLISTO_INDEX}" "${BACKUP}"
fi

########################################
# BUILD THE INDEX
########################################

echo "[INFO] Building kallisto index... (this may take a few minutes)"
echo "[INFO] Extra options: ${KALLISTO_EXTRA_OPTS}"

# Many versions of kallisto accept gzipped FASTA input directly; if your installation does not accept gz,
# you can uncompress to a temporary file before calling kallisto. Here we pass the path as-is.
# If you need to set a custom k-mer length, add e.g. "-k 23" to KALLISTO_EXTRA_OPTS above.
kallisto index -i "${KALLISTO_INDEX}" ${KALLISTO_EXTRA_OPTS} "${TRANSCRIPTS_FA}"

RET=$?
if [[ $RET -ne 0 ]]; then
    echo "[ERROR] kallisto index failed with exit code ${RET}" >&2
    exit $RET
fi

echo
echo "[INFO] Finished kallisto index build successfully."
echo "[INFO] Index created at: ${KALLISTO_INDEX}"

########################################
# POST-RUN CHECKS & NOTES
########################################

# Verify index file exists and non-empty
if [[ -s "${KALLISTO_INDEX}" ]]; then
    echo "[CHECK] Index file exists and is non-empty."
    if command -v md5sum >/dev/null 2>&1; then
        md5sum "${KALLISTO_INDEX}" | awk '{print "[PROV] index md5:", $1, $2}'
    fi
else
    echo "[ERROR] Index file was not created or is empty: ${KALLISTO_INDEX}" >&2
    exit 1
fi

# Print recommended next steps for users
cat <<EOF

[INFO] Recommendations / next steps:
 - Use this index with kallisto quant:
     kallisto quant -i ${KALLISTO_INDEX} -o output_dir -b 100 reads_1.fastq.gz reads_2.fastq.gz

 - If your reads are short (<50bp), consider rebuilding the index with a smaller k-mer length:
     kallisto index -i ${KALLISTO_INDEX} -k 23 ${TRANSCRIPTS_FA}

 - For reproducibility, record:
     * kallisto version: ${KALLISTO_VER}
     * transcript FASTA checksum and this log file: ${LOG}

 - If you encounter unexpected results during quantification, check that:
     * the FASTA contains transcript sequences (header lines usually contain transcript IDs)
     * transcript IDs match the IDs used by downstream annotations (e.g., tx2gene mapping)

EOF

echo "[INFO] End: $(date -u +"%Y-%m-%d %H:%M:%SZ")"
echo "=============================================="
