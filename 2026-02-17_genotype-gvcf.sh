#!/usr/bin/env bash
set -euo pipefail

# combine_and_genotype_recursive.sh (configured for your directories)
# - Recursively finds all per-sample gVCFs (*.g.vcf or *.g.vcf.gz) under GVCF_ROOT
# - bgzips & tabix-indexes any plain .g.vcf files (or indexes missing .g.vcf.gz)
# - Runs GATK CombineGVCFs to produce a combined .g.vcf.gz
# - Runs GATK GenotypeGVCFs on the combined gvcf to produce a joint VCF
# - Creates a timestamped run folder and writes a log
#
# This copy uses the directories you provided.

# ---------------- CONFIG ----------------
REF_FASTA="/media/jgd/jd-1tb/bioinformatics/reference/danRer11.fa"
GVCF_ROOT="/media/jgd/jd-1tb/bioinformatics/2024-04-01_proc-enu_paper/output/gatk"
OUTBASE="/media/jgd/jd-1tb/bioinformatics/2024-04-01_proc-enu_paper/output"
GATK="gatk"    # adjust if you need a full path
BGZIP="bgzip"
TABIX="tabix"
# Optional: set to 1 to overwrite existing outputs
FORCE_OVERWRITE=0
# -----------------------------------------

if [[ ! -f "$REF_FASTA" ]]; then
  echo "ERROR: REF_FASTA not found: $REF_FASTA" >&2
  exit 1
fi
if ! command -v "$GATK" >/dev/null 2>&1; then
  echo "ERROR: GATK not found at '$GATK' (or not in PATH)" >&2
  exit 1
fi
if ! command -v "$BGZIP" >/dev/null 2>&1 || ! command -v "$TABIX" >/dev/null 2>&1; then
  echo "ERROR: bgzip/tabix (htslib) required but not found" >&2
  exit 1
fi

TIMESTAMP="$(date +%Y%m%d_%H%M%S)"
OUTDIR="${OUTBASE%/}/joint_${TIMESTAMP}"
mkdir -p "$OUTDIR"
LOGFILE="${OUTDIR}/combine_genotype_${TIMESTAMP}.log"

# safe mktemp usage
TMPDIR="$(mktemp -d "${OUTDIR}/tmp.XXXXXX")"
trap 'rc=$?; rm -rf "$TMPDIR"; exit $rc' EXIT

echo "Run timestamp: $TIMESTAMP" | tee -a "$LOGFILE"
echo "GVCF root: $GVCF_ROOT" | tee -a "$LOGFILE"
echo "Output dir: $OUTDIR" | tee -a "$LOGFILE"
echo "Reference FASTA: $REF_FASTA" | tee -a "$LOGFILE"
echo "Using GATK: $(command -v "$GATK")" | tee -a "$LOGFILE"
echo "--------------------------------------------------" | tee -a "$LOGFILE"

# 1) Find gVCF files (case-insensitive for extension)
echo "[1/5] Scanning for g.vcf / g.vcf.gz files..." | tee -a "$LOGFILE"
mapfile -t RAW_GVCFS < <(find "$GVCF_ROOT" -type f \( -iname "*.g.vcf" -o -iname "*.g.vcf.gz" \) -print 2>/dev/null | sort)

N=${#RAW_GVCFS[@]}
if [[ $N -eq 0 ]]; then
  echo "ERROR: No g.vcf or g.vcf.gz files found under $GVCF_ROOT" | tee -a "$LOGFILE" >&2
  exit 1
fi
echo "Found $N files" | tee -a "$LOGFILE"

# 2) Ensure bgzipped + tabix-indexed; produce list of .g.vcf.gz paths
GZ_LIST="${TMPDIR}/gz_list.txt"
: > "$GZ_LIST"

echo "[2/5] Preparing (bgzip + tabix) gVCFs as needed..." | tee -a "$LOGFILE"
for f in "${RAW_GVCFS[@]}"; do
  # Normalize path
  f="$(readlink -f "$f")"
  if [[ "$f" == *.gz ]]; then
    gz="$f"
    if [[ ! -f "${gz}.tbi" ]]; then
      echo "Indexing (tabix) $gz" | tee -a "$LOGFILE"
      "$TABIX" -p vcf "$gz"
    else
      echo "Found gzipped+indexed: $gz" | tee -a "$LOGFILE"
    fi
  else
    gz="${f}.gz"
    if [[ -f "$gz" ]]; then
      echo "Found pre-existing $gz for $f" | tee -a "$LOGFILE"
    else
      echo "Bgzipping $f -> $gz" | tee -a "$LOGFILE"
      "$BGZIP" -c "$f" > "$gz"
    fi
    if [[ ! -f "${gz}.tbi" ]]; then
      echo "Indexing (tabix) $gz" | tee -a "$LOGFILE"
      "$TABIX" -p vcf "$gz"
    fi
  fi
  echo "$gz" >> "$GZ_LIST"
done

TOTAL=$(wc -l < "$GZ_LIST")
echo "Prepared $TOTAL gzipped gVCFs (paths saved to $GZ_LIST)" | tee -a "$LOGFILE"

# 3) CombineGVCFs
COMBINED_GVCF="${OUTDIR}/combined.g.vcf.gz"
if [[ -s "$COMBINED_GVCF" && $FORCE_OVERWRITE -eq 0 ]]; then
  echo "[3/5] Combined GVCF already exists: $COMBINED_GVCF (skipping combine)" | tee -a "$LOGFILE"
else
  echo "[3/5] Running GATK CombineGVCFs -> $COMBINED_GVCF" | tee -a "$LOGFILE"
  # Build args safely
  args=()
  while IFS= read -r g; do
    args+=(--variant "$g")
  done < "$GZ_LIST"

  ( set -x; "$GATK" CombineGVCFs -R "$REF_FASTA" "${args[@]}" -O "$COMBINED_GVCF" ) >>"$LOGFILE" 2>&1
  rc=$?
  if [[ $rc -ne 0 || ! -s "$COMBINED_GVCF" ]]; then
    echo "ERROR: CombineGVCFs failed (rc=$rc). See $LOGFILE" | tee -a "$LOGFILE" >&2
    exit 1
  fi
  # Index combined gVCF
  if [[ ! -f "${COMBINED_GVCF}.tbi" ]]; then
    echo "Indexing combined gVCF with tabix" | tee -a "$LOGFILE"
    "$TABIX" -p vcf "$COMBINED_GVCF"
  fi
fi

# 4) GenotypeGVCFs
JOINT_VCF="${OUTDIR}/joint_genotyped.vcf.gz"
if [[ -s "$JOINT_VCF" && $FORCE_OVERWRITE -eq 0 ]]; then
  echo "[4/5] Joint genotyped VCF already exists: $JOINT_VCF (skipping genotype)" | tee -a "$LOGFILE"
else
  echo "[4/5] Running GATK GenotypeGVCFs -> $JOINT_VCF" | tee -a "$LOGFILE"
  ( set -x; "$GATK" GenotypeGVCFs -R "$REF_FASTA" -V "$COMBINED_GVCF" -O "$JOINT_VCF" ) >>"$LOGFILE" 2>&1
  rc=$?
  if [[ $rc -ne 0 || ! -s "$JOINT_VCF" ]]; then
    echo "ERROR: GenotypeGVCFs failed (rc=$rc). See $LOGFILE" | tee -a "$LOGFILE" >&2
    exit 1
  fi
  # Index joint VCF
  if [[ ! -f "${JOINT_VCF}.tbi" ]]; then
    echo "Indexing joint VCF with tabix" | tee -a "$LOGFILE"
    "$TABIX" -p vcf "$JOINT_VCF"
  fi
fi

# 5) Finish
echo "--------------------------------------------------" | tee -a "$LOGFILE"
echo "SUCCESS: Combined GVCF: $COMBINED_GVCF" | tee -a "$LOGFILE"
echo "SUCCESS: Joint VCF:        $JOINT_VCF" | tee -a "$LOGFILE"
echo "Logfile: $LOGFILE" | tee -a "$LOGFILE"
echo "Temporary dir (removed on exit): $TMPDIR" | tee -a "$LOGFILE"