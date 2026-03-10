#!/usr/bin/env bash
#
# make_vcf_for_eqtl.sh
#
# Produce a normalized, filtered, bgzipped & indexed VCF suitable for downstream
# eQTL workflows, and convert it to PLINK/dosage files. Safe, self-contained,
# with sensible defaults. Edit variables or pass args to override.
#
# Usage:
#   bash make_vcf_for_eqtl.sh -i input.vcf.gz -r /path/to/ref.fa -o /path/to/outdir
#
# Example (using your files):
#   bash make_vcf_for_eqtl.sh \
#     -i /home/jgd/Documents/.../gatk-ensembl_filtered.20260224.nodash.vep.vcf.gz \
#     -r /home/jgd/Documents/.../danio_rerio.grcz11.dna.primary_assembly.fa.uncompressed.fa \
#     -o /home/jgd/Documents/bioinformatics_working/output
#
set -euo pipefail

# -----------------------
# Default params (edit or override via CLI)
# -----------------------
INPUT_VCF="/media/jgd/jd-linux1tb/2026-04-01_proc-enu_bioinformatics_results/snp_calling/gatk-ensembl/vep_annotated/gatk-ensembl_filtered.20260224.nodash.vep.vcf.gz"
REF_FA="/home/jgd/Documents/2026-04-01_proc-enu_bioinformatics_results/ref/danio_rerio.grcz11.dna.primary_assembly.fa.uncompressed.fa"
OUT_DIR="/home/jgd/Documents/bioinformatics_working/output"
MIN_QUAL=30         # filtering threshold
PLINK_GENO=0.05     # plink --geno threshold (variant missingness)
PLINK_MAF=0.01      # plink --maf threshold
PLINK_MIND=0.1      # plink --mind threshold (sample missingness)
KEEP_INTERMEDIATE=0 # 0 = remove large intermediate plain VCFs, 1 = keep
TMP_PREFIX="/tmp/make_vcf_$$"

# -----------------------
# Helpers / parse CLI
# -----------------------
usage() {
  cat <<EOF
Usage: $0 -i INPUT_VCF -r REF_FA -o OUT_DIR [options]

Required:
  -i INPUT_VCF   bgzipped VCF (or .vcf.gz that will be converted to bgzf)
  -r REF_FA      reference FASTA (unmasked primary assembly recommended)
  -o OUT_DIR     output directory

Options:
  -q MIN_QUAL    minimum QUAL to keep (default ${MIN_QUAL})
  --geno         plink --geno threshold for variant missingness (default ${PLINK_GENO})
  --maf          plink --maf threshold (default ${PLINK_MAF})
  --mind         plink --mind threshold (default ${PLINK_MIND})
  --keep-temp    keep intermediate plain VCFs (default remove)
  -h             show this help
EOF
  exit 1
}

# Minimal CLI parsing (supports long options for plink thresholds)
while [[ $# -gt 0 ]]; do
  case "$1" in
    -i) INPUT_VCF="$2"; shift 2;;
    -r) REF_FA="$2"; shift 2;;
    -o) OUT_DIR="$2"; shift 2;;
    -q) MIN_QUAL="$2"; shift 2;;
    --geno) PLINK_GENO="$2"; shift 2;;
    --maf) PLINK_MAF="$2"; shift 2;;
    --mind) PLINK_MIND="$2"; shift 2;;
    --keep-temp) KEEP_INTERMEDIATE=1; shift;;
    -h|--help) usage;;
    *) echo "Unknown arg: $1"; usage;;
  esac
done

if [[ -z "${INPUT_VCF}" || -z "${REF_FA}" || -z "${OUT_DIR}" ]]; then
  echo "ERROR: required args missing." >&2
  usage
fi

mkdir -p "${OUT_DIR}"

# -----------------------
# Tool checks
# -----------------------
for t in bcftools samtools bgzip tabix plink; do
  if ! command -v "${t}" >/dev/null 2>&1; then
    echo "ERROR: required tool not found in PATH: ${t}" >&2
    echo "Please install bcftools, samtools, bgzip, tabix, plink." >&2
    exit 2
  fi
done

# -----------------------
# Paths for outputs
# -----------------------
TMP_BGZ="${TMP_PREFIX}.input.bgz.vcf.gz"
NORMALIZED="${OUT_DIR}/normalized.vcf.gz"
FILTERED_GZ="${OUT_DIR}/filtered.vcf.gz"
FILTERED_PLAIN="${OUT_DIR}/filtered.vcf"
TAGGED_GZ="${OUT_DIR}/filtered.tags.vcf.gz"
PLINK_PREFIX="${OUT_DIR}/plink_filtered"

echo "INPUT_VCF: ${INPUT_VCF}"
echo "REF_FA:    ${REF_FA}"
echo "OUT_DIR:   ${OUT_DIR}"
echo "MIN_QUAL:  ${MIN_QUAL}"
echo "PLINK params: --geno ${PLINK_GENO} --maf ${PLINK_MAF} --mind ${PLINK_MIND}"
echo

# -----------------------
# Ensure reference FASTA index exists
# -----------------------
if [[ ! -f "${REF_FA}.fai" ]]; then
  echo "Indexing reference FASTA: samtools faidx ${REF_FA}"
  samtools faidx "${REF_FA}"
fi

# -----------------------
# Ensure input VCF is proper BGZF + index (make fresh bgzf copy if necessary)
# -----------------------
echo "Ensuring input VCF is BGZF/indexed..."
if tabix -l "${INPUT_VCF}" >/dev/null 2>&1; then
  echo "Input VCF already indexed and BGZF: using ${INPUT_VCF}"
  TMP_BGZ="${INPUT_VCF}"
else
  echo "Creating BGZF copy -> ${TMP_BGZ}"
  gunzip -c "${INPUT_VCF}" | bgzip -c > "${TMP_BGZ}"
  tabix -p vcf "${TMP_BGZ}"
  echo "Created and indexed TMP BGZF."
fi

# -----------------------
# Normalize: split multiallelic alleles and left-align using the provided reference
# -----------------------
echo "Normalizing (split multiallelics and left-align) using ${REF_FA} ..."
bcftools norm -m -both -f "${REF_FA}" "${TMP_BGZ}" -Oz -o "${NORMALIZED}"
tabix -p vcf "${NORMALIZED}"
echo "Normalized VCF: ${NORMALIZED}"

# -----------------------
# Filter: keep PASS & QUAL>=MIN_QUAL. Write plain VCF then compress with bgzip to avoid BGZF issues.
# -----------------------
echo "Filtering PASS & QUAL>=${MIN_QUAL} -> writing plain VCF ${FILTERED_PLAIN} ..."
bcftools view -f PASS -i "QUAL>=${MIN_QUAL}" "${NORMALIZED}" -Ov -o "${FILTERED_PLAIN}"

echo "Compressing filtered VCF -> ${FILTERED_GZ} and indexing ..."
bgzip -c "${FILTERED_PLAIN}" > "${FILTERED_GZ}"
tabix -p vcf "${FILTERED_GZ}"
echo "Filtered VCF (bgz): ${FILTERED_GZ}"

# -----------------------
# Annotate MAF (fill-tags) and produce tagged bgz
# -----------------------
echo "Annotating MAF -> plain tagged VCF then bgzip -> ${TAGGED_GZ} ..."
bcftools +fill-tags "${FILTERED_GZ}" -Ov -o "${OUT_DIR}/filtered.tags.vcf" -- -t MAF
bgzip -c "${OUT_DIR}/filtered.tags.vcf" > "${TAGGED_GZ}"
tabix -p vcf "${TAGGED_GZ}"
rm -f "${OUT_DIR}/filtered.tags.vcf"
echo "Tagged VCF: ${TAGGED_GZ}"

# -----------------------
# Convert to PLINK files and run PLINK QC and export additive dosages
# -----------------------
echo "Converting ${TAGGED_GZ} -> PLINK bed/bim/fam ..."
# plink (v1.9+) accepts --vcf; use --vcf and --make-bed
plink --vcf "${TAGGED_GZ}" --make-bed --out "${PLINK_PREFIX}.raw" || { echo "plink conversion failed"; exit 3; }

echo "Running PLINK QC: --geno ${PLINK_GENO} --maf ${PLINK_MAF} --mind ${PLINK_MIND} ..."
plink --bfile "${PLINK_PREFIX}.raw" --geno "${PLINK_GENO}" --maf "${PLINK_MAF}" --mind "${PLINK_MIND}" --make-bed --out "${PLINK_PREFIX}.qc" || { echo "plink QC failed"; exit 4; }

echo "Exporting additive dosages (numeric genotypes) -> ${PLINK_PREFIX}.qc.dosage.raw ..."
# --recode A produces .raw file with additive allele counts (0/1/2)
plink --bfile "${PLINK_PREFIX}.qc" --recode A --out "${PLINK_PREFIX}.qc.dosage" || { echo "plink recode A failed"; exit 5; }

echo "PLINK outputs (prefix ${PLINK_PREFIX}.*):"
ls -lh "${PLINK_PREFIX}".raw.* "${PLINK_PREFIX}".qc.* "${PLINK_PREFIX}".qc.dosage.* 2>/dev/null || true

# -----------------------
# Clean up (optionally)
# -----------------------
if [[ "${KEEP_INTERMEDIATE}" -eq 0 ]]; then
  echo "Removing intermediate plain VCF(s) and tmp files..."
  rm -f "${FILTERED_PLAIN}" "${TMP_PREFIX}"* || true
fi

echo
echo "DONE."
echo "Files you can use for eQTL tools:"
echo "  Normalized VCF:   ${NORMALIZED}"
echo "  Filtered VCF:     ${FILTERED_GZ}"
echo "  Tagged VCF (MAF): ${TAGGED_GZ}"
echo "  PLINK (QC'd):     ${PLINK_PREFIX}.qc.{bed,bim,fam}"
echo "  Dosage/raw table: ${PLINK_PREFIX}.qc.dosage.raw"
echo
echo "Notes / next steps:"
echo " - Use the PLINK dosage/raw file (samples in columns) as genotype input for Matrix eQTL or convert to the format your chosen eQTL tool requires."
echo " - For cis-eQTL, restrict variants by distance from each gene (e.g., +/-1Mb) prior to running associations."
echo " - If you need the genotype matrix in a different format, tell me which tool and I will provide the exact conversion command."
