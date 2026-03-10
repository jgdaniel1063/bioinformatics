#!/usr/bin/env bash
#
# haplotype_from_bam_gene_FULL_annotated.sh
#
# Heavily annotated version of the haplotype extraction pipeline (Linux, no Slurm).
#
# Purpose
# -------
# Given:
#  - one RNA-seq BAM (aligned to danRer11 / GRCz11)
#  - a cohort VCF (bgzipped + indexed) that contains genotype calls for the sample
#  - a reference FASTA (bgzipped/indexed)
#  - a gene annotation GTF (.gz + .tbi)
#
# This script:
#  1) Finds the genomic interval for a gene by searching the GTF for gene_name == GENE_QUERY
#     (simple grep-first heuristic; see notes below for robust alternatives).
#  2) Subsets the cohort VCF to that interval and the specified sample.
#  3) Uses whatshap to do read-backed phasing of variants in that interval (phases GT fields).
#  4) Optionally haplotags the input BAM with haplotype labels using whatshap haplotag.
#  5) Optionally produces haplotype-specific consensus FASTA sequences for the region
#     using bcftools consensus (compatible with bcftools 1.19).
#
# Output layout (under OUTPUT_ROOT/RUN_TAG):
#   inputs/        : saved interval metadata
#   vcf/           : region-subset VCF + phased VCF
#   bam/           : optionally haplotagged BAM
#   haplotypes/    : optionally per-haplotype consensus FASTA
#   logs/          : logs from whatshap and other steps
#   tmp/           : temporary files
#
# Design goals
# ------------
# - Self-contained single-script runnable on a typical bioinformatics Linux host
# - Defensive checks for required programs and files
# - Conservative error handling and informative messages
# - Practical defaults suitable for interactive use; convert to a CLI tool for production use
#
# Safety / scaling notes
# ----------------------
# - This script fetches data for a single gene region; it is NOT optimized for genome-wide runs.
# - For many genes or many samples, parallelize per-gene or per-chromosome jobs (e.g., GNU parallel or Slurm arrays).
# - For reproducibility, pin versions of samtools/bcftools/whatshap and record them in the run directory.
#
# ------------------------------------------------------------------------------

set -euo pipefail          # strict error handling
shopt -s nullglob          # make globs that don't match expand to nothing
# Trap to print informative error including failing command and line
trap 'ec=$?; echo "❌ ERROR (exit=$ec) at line $LINENO: $BASH_COMMAND" >&2; exit $ec' ERR

# Simple logger function used throughout; writes to stderr with timestamp
log(){ echo "[$(date '+%F %T')] $*" >&2; }

# Ensure per-user installs are visible (useful when whatshap is installed via `pip --user`)
export PATH="$HOME/.local/bin:$PATH"
hash -r

# ============================================================
# USER EDITS — change only these variables for a run
# ============================================================

# Inputs: BAM (single-sample), VCF (cohort), sample ID in the VCF, reference FASTA, GTF
BAM="/home/jgd/Documents/bioinformatics_working/star_aligned/cohort1&2/bam/no_thrombosis/ai08nothrom1.bam"
VCF="/home/jgd/Documents/bioinformatics_working/snp_calling/gatk/cohort3/gvcf/combinegvcf_genotypegvcf/cohort.raw.vcf.gz"
VCF_SAMPLE="ai08nothrom1"
REF_FA="/home/jgd/Documents/bioinformatics_working/reference/danRer11.fa"
GTF="/home/jgd/Documents/bioinformatics_working/reference/danio_rerio.grcz11.110.gtf.gz"

# Gene to extract: case-insensitive match against gene_name attribute in GTF
GENE_QUERY="mmp9"          # use gene_name string (e.g., 'mmp9')
PAD_BP=2000                # padding added to both sides of matched interval
MAX_INTERVAL_BP=5000000    # safety threshold to prevent huge intervals

# whatshap options: array expanded when calling whatshap
WHATSHAP_PHASE_OPTS=(--ignore-read-groups --indels)

# toggles: whether to haplotag BAM and/or create consensus FASTA for haplotypes
DO_HAPLOTAG_BAM=1
DO_CONSENSUS_FASTA=1

# Where to write all outputs for this run (each run gets a timestamped subdirectory)
OUTPUT_ROOT="/home/jgd/Documents/bioinformatics_working/output"

# Optional Python virtualenv activation (leave blank if not used)
VENV="${VENV:-}"   # e.g. /home/user/envs/haplo_env

# Threads for samtools indexing and similar tasks
THREADS="${THREADS:-4}"

# ============================================================
# SETUP: create run-specific directory and record inputs
# ============================================================
RUN_TAG="haplotype_from_bam_${GENE_QUERY}_$(date '+%Y%m%d_%H%M%S')"
RUN_DIR="${OUTPUT_ROOT}/${RUN_TAG}"

log "RUN_DIR=${RUN_DIR}"
log "BAM=${BAM}"
log "VCF=${VCF}"
log "VCF_SAMPLE=${VCF_SAMPLE}"
log "REF_FA=${REF_FA}"
log "GTF=${GTF}"
log "GENE_QUERY=${GENE_QUERY}  PAD_BP=${PAD_BP}"

# create directory structure; be permissive (parents)
mkdir -p "${OUTPUT_ROOT}"
mkdir -p "${RUN_DIR}"/{inputs,vcf,bam,haplotypes,logs,tmp}
export TMPDIR="${TMPDIR:-${RUN_DIR}/tmp}"
mkdir -p "$TMPDIR"

# ============================================================
# OPTIONAL: ACTIVATE VENV (if provided)
# ============================================================
if [[ -n "${VENV}" ]]; then
  if [[ -f "${VENV}/bin/activate" ]]; then
    # shellcheck disable=SC1090
    source "${VENV}/bin/activate"
    log "Activated venv: ${VENV}"
  else
    echo "❌ VENV set but activate not found: ${VENV}/bin/activate" >&2
    exit 1
  fi
fi

# ============================================================
# TOOL CHECKS — fail early with helpful message if missing
# ============================================================
need(){ command -v "$1" >/dev/null 2>&1 || { echo "❌ Missing required tool: $1" >&2; exit 1; }; }
need samtools
need bcftools
need bgzip
need tabix
need whatshap

log "samtools:  $(command -v samtools)"
log "bcftools:  $(command -v bcftools)"
log "whatshap:  $(command -v whatshap)"

# Note: whatshap must support 'phase' and 'haplotag' subcommands; prefer version >= 1.0

# ============================================================
# VALIDATION — verify files and indexes exist and are non-empty
# ============================================================
[[ -s "${BAM}" ]] || { echo "❌ Missing BAM: ${BAM}" >&2; exit 1; }
[[ -s "${VCF}" ]] || { echo "❌ Missing VCF: ${VCF}" >&2; exit 1; }
[[ -s "${REF_FA}" ]] || { echo "❌ Missing REF_FA: ${REF_FA}" >&2; exit 1; }
[[ -s "${GTF}" ]] || { echo "❌ Missing GTF: ${GTF}" >&2; exit 1; }

# Check BAM index (either .bam.bai or .bai beside the file)
[[ -s "${BAM}.bai" || -s "${BAM%.bam}.bai" ]] || { echo "❌ Missing BAM index (.bai) for: ${BAM}" >&2; exit 1; }
# Check VCF index (.tbi or .csi)
[[ -s "${VCF}.tbi" || -s "${VCF}.csi" ]] || { echo "❌ Missing VCF index (.tbi/.csi) for: ${VCF}" >&2; exit 1; }
# Check FASTA index (.fai)
[[ -s "${REF_FA}.fai" ]] || { echo "❌ Missing FASTA index (.fai) for: ${REF_FA}" >&2; exit 1; }

# Confirm sample exists in VCF
if ! bcftools query -l "${VCF}" | grep -Fxq "${VCF_SAMPLE}"; then
  echo "❌ VCF sample not found: ${VCF_SAMPLE}" >&2
  log "Available samples (first 50):"
  bcftools query -l "${VCF}" | head -n 50 >&2
  exit 1
fi

# ============================================================
# FIND GENE INTERVAL (simple grep-first on gene_name in GTF)
# ============================================================
# We use a simple and fast heuristic: grep for gene_name "GENE_QUERY" (case-insensitive),
# then use awk to compute min(start) and max(end) across matched features.
#
# Caveats:
# - This relies on GTF attribute key being exactly 'gene_name "X"' (common for Ensembl GTFs).
# - If the gene appears on multiple contigs (unexpected), the script aborts.
# - For robust behavior in pipelines consider using a proper GTF parser (e.g., gffutils, pybedtools)
#   or a gene-to-coordinate lookup table (e.g., precomputed gene intervals).
GTF_STREAM="gunzip -c"
[[ "${GTF}" != *.gz ]] && GTF_STREAM="cat"

log "Searching GTF for gene_name \"${GENE_QUERY}\" (case-insensitive; grep-first) ..."

# The awk accumulates min and max coordinates; it also detects if matches occur on multiple contigs.
read -r CHROM START END MATCH_COUNT <<<"$(
  ${GTF_STREAM} "${GTF}" \
  | grep -i -F "gene_name \"${GENE_QUERY}\"" \
  | awk -F'\t' '
      BEGIN{min=0;max=0;cnt=0;chrom="";multi=0}
      NF<5{next}
      {
        cnt++
        if (min==0 || $4<min) min=$4
        if (max==0 || $5>max) max=$5
        if (chrom=="" ) chrom=$1
        else if ($1!=chrom) { multi=1 }
      }
      END{
        if (cnt==0) { print "NA 0 0 0"; exit 0 }
        if (multi==1) { print "MULTI", min, max, cnt; exit 0 }
        print chrom, min, max, cnt
      }'
)"

# If no matches, abort with helpful diagnostics
if [[ "${MATCH_COUNT}" -eq 0 || "${CHROM}" == "NA" ]]; then
  echo "❌ No matching records found in GTF for: ${GENE_QUERY}" >&2
  log "Sanity checks (first lines matching):"
  log "  ${GTF_STREAM} ${GTF} | grep -i -F 'gene_name \"${GENE_QUERY}\"' | head"
  log "  ${GTF_STREAM} ${GTF} | grep -i '${GENE_QUERY}' | head"
  exit 1
fi

# If matches span multiple contigs, avoid guessing — abort for manual inspection
if [[ "${CHROM}" == "MULTI" ]]; then
  echo "❌ Matched gene_name \"${GENE_QUERY}\" on multiple contigs in GTF (unexpected). Aborting." >&2
  log "Show first matches to inspect:"
  ${GTF_STREAM} "${GTF}" | grep -i -F "gene_name \"${GENE_QUERY}\"" | head -n 20 >&2
  exit 1
fi

# Pad the interval and sanity-check its size
START_PAD=$(( START - PAD_BP )); [[ "${START_PAD}" -lt 1 ]] && START_PAD=1
END_PAD=$(( END + PAD_BP ))
INTERVAL_BP=$(( END_PAD - START_PAD + 1 ))

if [[ "${INTERVAL_BP}" -gt "${MAX_INTERVAL_BP}" ]]; then
  echo "❌ Interval too large (${INTERVAL_BP} bp). Likely wrong match. Aborting." >&2
  echo "   Matched ${MATCH_COUNT} features on ${CHROM}:${START}-${END}" >&2
  exit 1
fi

REGION="${CHROM}:${START_PAD}-${END_PAD}"
log "GTF match count: ${MATCH_COUNT}"
log "Resolved gene interval: ${CHROM}:${START}-${END}"
log "Using padded region:    ${REGION}"

# Save interval metadata for provenance
printf "%s\t%s\t%s\t%s\t%s\n" "${GENE_QUERY}" "${CHROM}" "${START}" "${END}" "${REGION}" \
  > "${RUN_DIR}/inputs/gene_interval.tsv"

# ============================================================
# CONTIG STYLE FIXUPS: handle 'chr' prefix differences
# ============================================================
# Many pipelines mix contig name styles: "1" vs "chr1". If the contig name from the GTF
# is not present in the BAM and/or VCF, try the alternate form by adding/removing "chr".
bam_has_contig(){
  local contig="$1"
  # samtools idxstats prints contigs present; exit code 0 if found, 1 otherwise
  samtools idxstats "${BAM}" | awk -v c="${contig}" '$1==c {found=1} END{exit(found?0:1)}'
}
vcf_has_contig(){
  local contig="$1"
  # Parse VCF header lines like: ##contig=<ID=chr1,length=...>
  bcftools view -h "${VCF}" | awk -v c="${contig}" '
    $0 ~ /^##contig=<ID=/ {
      id=$0; sub(/.*ID=/,"",id); sub(/[,>].*/,"",id);
      if (id==c) {found=1}
    }
    END{exit(found?0:1)}'
}

CONTIG="${CHROM}"
if ! bam_has_contig "${CONTIG}" || ! vcf_has_contig "${CONTIG}"; then
  alt="${CONTIG}"
  if [[ "${CONTIG}" == chr* ]]; then
    alt="${CONTIG#chr}"
  else
    alt="chr${CONTIG}"
  fi

  log "Contig ${CONTIG} not found in BAM and/or VCF; trying alternate: ${alt}"
  if bam_has_contig "${alt}" && vcf_has_contig "${alt}"; then
    CONTIG="${alt}"
    REGION="${CONTIG}:${START_PAD}-${END_PAD}"
    log "Using alternate contig. New region: ${REGION}"
  else
    echo "❌ Region contig not found in BAM/VCF: ${CHROM} (also tried ${alt})" >&2
    log "BAM contigs sample:"
    samtools idxstats "${BAM}" | head -n 10 >&2
    log "VCF contigs sample:"
    bcftools view -h "${VCF}" | grep '^##contig=' | head -n 10 >&2
    exit 1
  fi
fi

# Save the region used
echo "${REGION}" > "${RUN_DIR}/inputs/region.txt"

# ============================================================
# SUBSET VCF (region + sample) — create per-run region/sample VCF
# ============================================================
SUB_VCF="${RUN_DIR}/vcf/${VCF_SAMPLE}.${GENE_QUERY}.region.vcf.gz"
log "Subsetting VCF -> ${SUB_VCF}"

# bcftools view -r region -s sample -Oz compress and output bgzipped VCF
bcftools view -r "${REGION}" -s "${VCF_SAMPLE}" -Oz -o "${SUB_VCF}" "${VCF}"
# Index the new VCF for random access
tabix -f -p vcf "${SUB_VCF}"

# Count records to know if there are variants to phase
NREC="$(bcftools view -H "${SUB_VCF}" | wc -l | tr -d ' ')"
log "Subset VCF records in region: ${NREC}"
if [[ "${NREC}" -eq 0 ]]; then
  echo "⚠️  No variants in region for ${VCF_SAMPLE}. Nothing to phase." >&2
  echo "   Region: ${REGION}" >&2
  exit 0
fi

# ============================================================
# PHASE WITH WHATSHAP — read-backed phasing using BAM
# ============================================================
PHASED_VCF="${RUN_DIR}/vcf/${VCF_SAMPLE}.${GENE_QUERY}.phased.vcf.gz"
PHASED_UNZ="${RUN_DIR}/vcf/${VCF_SAMPLE}.${GENE_QUERY}.phased.vcf"
PHASE_LOG="${RUN_DIR}/logs/whatshap_phase.log"

log "Phasing with whatshap -> ${PHASED_VCF}"
log "whatshap options: ${WHATSHAP_PHASE_OPTS[*]}"

# whatshap phase usage:
#  whatshap phase [options] -o out.vcf input.vcf sample.bam
# We redirect stdout/stderr to a log
whatshap phase \
  "${WHATSHAP_PHASE_OPTS[@]}" \
  --reference "${REF_FA}" \
  --output "${PHASED_UNZ}" \
  "${SUB_VCF}" \
  "${BAM}" \
  >"${PHASE_LOG}" 2>&1

# Compress and index the phased VCF
bgzip -f -c "${PHASED_UNZ}" > "${PHASED_VCF}"
rm -f "${PHASED_UNZ}"
tabix -f -p vcf "${PHASED_VCF}"

# Quick check: count phased GTs (those using '|' instead of '/')
PHASED_CT="$(bcftools view -H "${PHASED_VCF}" | awk -F'\t' '
  {
    fmt=$9; sample=$10;
    n=split(fmt,keys,":"); m=split(sample,vals,":");
    gt=""; for(i=1;i<=n;i++){ if(keys[i]=="GT"){ gt=vals[i] } }
    if (gt ~ /\|/) phased++
    total++
  }
  END{print phased+0, total+0}'
)"
log "Phased GT count / total in region: ${PHASED_CT}"
# Explanation: the awk prints "phased total" where phased is number of records with phased GT.

# ============================================================
# OPTIONAL: HAPLOTAG BAM using whatshap haplotag (adds HP tag)
# ============================================================
if [[ "${DO_HAPLOTAG_BAM}" -eq 1 ]]; then
  HT_BAM="${RUN_DIR}/bam/${VCF_SAMPLE}.${GENE_QUERY}.haplotag.bam"
  HT_LOG="${RUN_DIR}/logs/whatshap_haplotag.log"

  log "Haplotagging BAM -> ${HT_BAM}"
  # whatshap haplotag usage:
  #  whatshap haplotag --output out.bam phased.vcf input.bam
  whatshap haplotag \
    --reference "${REF_FA}" \
    --output "${HT_BAM}" \
    --regions "${REGION}" \
    "${PHASED_VCF}" \
    "${BAM}" \
    >"${HT_LOG}" 2>&1

  # Index the haplotagged BAM; newer samtools versions accept -@ for threads
  samtools index -@ "${THREADS}" "${HT_BAM}"
fi

# ============================================================
# OPTIONAL: GENERATE HAPLOTYPE-SPECIFIC CONSENSUS FASTA
# Notes about bcftools consensus compatibility:
#  - bcftools consensus (old/new) command-line syntax changed across versions.
#  - Here we use the approach compatible with bcftools 1.19 as documented:
#      samtools faidx REF region | bcftools consensus --fasta-ref REF -s SAMPLE -H hap_index PHASED_VCF > hap.fa
#  - The region FASTA is provided on stdin; bcftools consensus uses the FASTA header to know the coordinates.
# ============================================================
if [[ "${DO_CONSENSUS_FASTA}" -eq 1 ]]; then
  HAP1_FA="${RUN_DIR}/haplotypes/${VCF_SAMPLE}.${GENE_QUERY}.hap1.fa"
  HAP2_FA="${RUN_DIR}/haplotypes/${VCF_SAMPLE}.${GENE_QUERY}.hap2.fa"
  REF_REG_FA="${RUN_DIR}/haplotypes/${VCF_SAMPLE}.${GENE_QUERY}.ref.fa"

  log "Writing reference sequence for region -> ${REF_REG_FA}"
  # Extract the region sequence; fasta header will be like >chr:start-end
  samtools faidx "${REF_FA}" "${REGION}" > "${REF_REG_FA}"

  log "Building consensus haplotypes from phased VCF (bcftools-1.19 compatible)"
  # Haplotype 1 (-H 1) and haplotype 2 (-H 2)
  samtools faidx "${REF_FA}" "${REGION}" \
    | bcftools consensus --fasta-ref "${REF_FA}" -s "${VCF_SAMPLE}" -H 1 "${PHASED_VCF}" \
    > "${HAP1_FA}"

  samtools faidx "${REF_FA}" "${REGION}" \
    | bcftools consensus --fasta-ref "${REF_FA}" -s "${VCF_SAMPLE}" -H 2 "${PHASED_VCF}" \
    > "${HAP2_FA}"
  # Note:
  # - If a sample is homozygous at many sites or many sites are unphased, you may see identical hap1/hap2.
  # - Consider soft-masking or annotating positions with low genotype quality before making consensus sequences.
fi

# ============================================================
# FINISH: print summary and relevant artifact paths
# ============================================================
log "DONE"
log "Region: ${REGION}"
log "Subset VCF:  ${SUB_VCF}"
log "Phased VCF:  ${PHASED_VCF}"
[[ "${DO_HAPLOTAG_BAM}" -eq 1 ]] && log "Haplotag BAM: ${RUN_DIR}/bam/${VCF_SAMPLE}.${GENE_QUERY}.haplotag.bam"
[[ "${DO_CONSENSUS_FASTA}" -eq 1 ]] && log "Hap1 FASTA:   ${RUN_DIR}/haplotypes/${VCF_SAMPLE}.${GENE_QUERY}.hap1.fa"
[[ "${DO_CONSENSUS_FASTA}" -eq 1 ]] && log "Hap2 FASTA:   ${RUN_DIR}/haplotypes/${VCF_SAMPLE}.${GENE_QUERY}.hap2.fa"

# ==============================================================================
# FOLLOW-UP NOTES, TROUBLESHOOTING, AND SUGGESTED ENHANCEMENTS
# ==============================================================================

# 1) Gene interval determination:
#    - This script uses a simple grep + awk on the GTF to locate gene_name matches.
#      This is fast and works for standard Ensembl-style GTFs, but fails if:
#        * your GTF uses different attribute keys (e.g., gene_symbol instead of gene_name)
#        * the gene has multiple loci (paralogs) resulting in matches on multiple contigs
#      For robust behavior consider:
#        * Querying a precomputed gene BED file (one line per gene) that you create once.
#        * Using Python with HTSeq/gffutils/pybedtools to parse GTF and lookup by attribute.

# 2) Contig name mismatches:
#    - Many alignment/VCF/annotation pipelines differ in contig naming (chr-prefixed or not).
#    - We try a simple alternate ('chr' + or without 'chr') guess. If you repeatedly have mismatches,
#      normalize contig names across your files up-front with a consistent convention.

# 3) Phasing considerations:
#    - whatshap uses reads (BAM) to phase heterozygous variants. It is read-backed and works well when:
#        * the reads span multiple heterozygous sites (longer reads, paired-end fragments)
#        * genotype calls are available in the VCF (not only gVCF records)
#    - Whatshap may produce phased blocks separated by regions without read links; inspect the .vcf and log.
#    - For RNA-seq data, read coverage may be uneven and phasing may be partial — expect many unphased sites.

# 4) Haplotagging:
#    - whatshap haplotag writes an HP tag (haplotype) per read into the BAM. Downstream tools can use HP/PS tags.
#    - Haplotagging is reversible (you can re-generate the original BAM) but be careful to use copies, not original files.

# 5) Consensus FASTA caveats:
#    - bcftools consensus will substitute alleles from the VCF into the reference. If VCF contains multi-allelic sites
#      or complex indels, consensus behavior depends on bcftools version. Always inspect the resulting FASTA.
#    - Consider masking low-confidence positions (e.g., genotype quality or read depth) prior to consensus.

# 6) Multi-allelic and structural variants:
#    - This pipeline focuses on small variants (SNPs and short indels). Structural variants and complex events
#      are not handled; consensus generation for those requires additional logic.

# 7) Parallelization and scaling:
#    - To phase many genes/samples, split the gene list into tasks and run this script in parallel (e.g., GNU parallel,
#      Slurm array jobs). Keep in mind temporary file isolation (each run_dir should be unique).
#    - When running many instances, consider using a shared fast disk for temporary files or a per-node tmpdir.

# 8) Reproducibility:
#    - Save versions of tools used:
#        samtools --version
#        bcftools --version
#        whatshap --version
#      Consider writing them to RUN_DIR/logs/tools_versions.txt
#    - Consider containerizing (Docker/Singularity) to ensure identical tool versions.

# 9) Improvements you may ask me to implement:
#    - Convert this script into a parameterized CLI with getopt/argparse-like checking.
#    - Expand multi-allelic records to per-ALT outputs, producing one consensus per allele.
#    - Add per-variant filters before phasing (min DP, min GQ) for higher-confidence haplotypes.
#    - Produce a small HTML report summarizing reads supporting each haplotype and phased block lengths.
#
# End of script.
