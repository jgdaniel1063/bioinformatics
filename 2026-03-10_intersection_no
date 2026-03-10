#!/usr/bin/env bash
#
# enu_family_HET_present_absent_GLOBAL_BG_SUBTRACTION_STRICT_local_annotated.sh
#
# Heavily annotated local version of the family-based VCF intersection and
# background-subtraction script — all cluster/module/SLURM-specific lines removed.
#
# Differences from the cluster version:
#  - Removed module load / SLURM directives.
#  - Assumes required tools (bcftools, tabix, bgzip, awk, sort, etc.) are on PATH.
#  - Performs the same per-family operations (present/absent/specific/noBG/unique)
#    and optional VEP annotation, writing everything into a timestamped RUN_DIR.
#
# Usage:
#   - Edit the configuration block below (COHORT_VCF, OUT_BASE, thresholds, VEP paths).
#   - Ensure bcftools/tabix/bgzip and other utilities are installed and available on PATH.
#   - Run on a workstation or compute node:
#       bash enu_family_HET_present_absent_GLOBAL_BG_SUBTRACTION_STRICT_local_annotated.sh
#
# Notes:
#  - This script can be run interactively; for many families you may want to run the
#    per-family loop in parallel (GNU parallel or separate jobs).
#  - Verify you have enough disk space for intermediate family VCFs and VEP outputs.
#
set -euo pipefail
shopt -s nullglob
trap 'ec=$?; echo "❌ ERROR (exit=$ec) at line $LINENO: $BASH_COMMAND" >&2; exit $ec' ERR

# Ensure consistent locale (helps sorting/awk behavior)
export LC_ALL=C
export LANG=C

# -----------------------------------------------------------------------------
# Simple logger helper
# -----------------------------------------------------------------------------
log(){
  local ts msg
  ts="$(date '+%F %T')"
  msg="[$ts] $*"
  if [[ -n "${LOGFILE:-}" ]]; then
    echo "${msg}" | tee -a "${LOGFILE}" >&2
  else
    echo "${msg}" >&2
  fi
}

# =============================================================================
# CONFIGURATION / KNOBS (edit these)
# =============================================================================

# Inputs / outputs (hard-coded)
COHORT_VCF="/nfs/turbo/umms-jshavit/jgdaniel/06-03-2024_proc-enu_paper/11-19-2025_final_results/gatk_workflow/snp_calling/gatk-ensembl_filtered.20260224.nodash.vcf.gz"
OUT_BASE="/nfs/turbo/umms-jshavit/jgdaniel/06-03-2024_proc-enu_paper/raw_output"
BGZIP_THREADS="8"    # threads to pass to bcftools view -Oz or other bgzip-enabled steps

# Regexes to detect 'throm' and 'nothrom' tokens in sample identifiers.
THROM_RE='(thrombosis|throm)'
NOTHROM_RE='(nothrombosis|nothrom)'

# Site-level filters
SITE_FILTER_ENABLE=1
SITE_REQUIRE_FILTER_PASS=0
SITE_TYPES_ENABLE=1
SITE_TYPES="snps"
SITE_BIALLELIC_ONLY=1
SITE_MIN_QUAL=80

# Strict-present thresholds (heterozygous evidence)
MIN_DP_PRESENT=75
MIN_GQ_PRESENT=30
MIN_ALT_AD_PRESENT=20
MIN_AB_HET=0.45
MAX_AB_HET=0.55

# Strict-absent thresholds (homozygous-reference evidence)
MIN_DP_ABSENT=75
MAX_ALT_ADSENT=0
MAX_ALT_AD_ABSENT=0
MAX_AB_ABSENT=0.01

# Global background settings
BG_ENABLE=1
BG_MIN_DP=10
BG_MIN_SAMPLES_ALL=3

# Global-unique filter (optional)
UNIQUE_ENABLE=1
UNIQUE_DP=20

# Emit per-bucket VCFs?
EMIT_BUCKET_VCFS=1

# VEP (optional)
VEP_ENABLE=1
VEP_THREADS=6
VEP_CACHE_DIR="${HOME}/vep_cache_grcz11"
VEP_EXEC="/home/jgdaniel/miniconda3/envs/vep_env/share/ensembl-vep-115.2-1/vep"
VEP_SPECIES="danio_rerio"
VEP_ASSEMBLY="GRCz11"
VEP_EXTRA_ARGS="--symbol --canonical --no_stats"

# =============================================================================
# Tool checks (fail early if tools missing)
# =============================================================================
BCT="bcftools"
TABIX="tabix"
BGZIP="bgzip"

need(){ command -v "$1" >/dev/null 2>&1 || { echo "❌ Missing tool: $1" >&2; exit 1; }; }

# Ensure required utilities are available
need "${BCT}"; need "${TABIX}"; need "${BGZIP}"
need awk; need sort; need uniq; need wc; need head; need paste; need mktemp; need mkdir; need cut; need tee; need grep

# -----------------------------------------------------------------------------
# VEP setup — prefer explicit VEP_EXEC if provided; else try to activate vep_env (if conda present)
# -----------------------------------------------------------------------------
if [[ "${VEP_ENABLE}" -eq 1 ]]; then
  VEP_CMD=""
  if [[ -n "${VEP_EXEC:-}" && -x "${VEP_EXEC}" ]]; then
    VEP_CMD="${VEP_EXEC}"
    log "VEP: using explicit VEP_EXEC=${VEP_EXEC}"
  else
    # Try to activate conda environment if available (not cluster-specific)
    if [[ -f "${HOME}/miniconda3/etc/profile.d/conda.sh" ]]; then
      # shellcheck disable=SC1090
      source "${HOME}/miniconda3/etc/profile.d/conda.sh"
      if conda activate vep_env >/dev/null 2>&1 && command -v vep >/dev/null 2>&1; then
        VEP_CMD="vep"
        log "VEP: using vep from activated vep_env"
      fi
    fi
    if [[ -z "${VEP_CMD}" ]]; then
      echo "❌ vep not found (neither VEP_EXEC nor vep in vep_env). Set VEP_EXEC or install VEP." >&2
      exit 1
    fi
  fi
  if [[ -n "${VEP_CACHE_DIR}" && ! -d "${VEP_CACHE_DIR}" ]]; then
    echo "❌ VEP cache directory not found: ${VEP_CACHE_DIR}" >&2
    exit 1
  fi
fi

# -----------------------------------------------------------------------------
# Basic input validation
# -----------------------------------------------------------------------------
[[ -s "${COHORT_VCF}" ]] || { echo "❌ Missing COHORT_VCF: ${COHORT_VCF}" >&2; exit 1; }
[[ -s "${COHORT_VCF}.tbi" || -s "${COHORT_VCF}.csi" ]] || { echo "❌ Missing VCF index (.tbi/.csi): ${COHORT_VCF}" >&2; exit 1; }

# Ensure FORMAT/AD exists in header
HDR_TMP="$(mktemp)"
"${BCT}" view -h "${COHORT_VCF}" > "${HDR_TMP}"
if ! grep -q '^##FORMAT=<ID=AD' "${HDR_TMP}"; then
  echo "❌ This script requires FORMAT/AD in the VCF (allelic depths). Aborting." >&2
  grep -n '^##FORMAT=' "${HDR_TMP}" | head -n 80 >&2 || true
  rm -f "${HDR_TMP}"
  exit 1
fi
rm -f "${HDR_TMP}"

# -----------------------------------------------------------------------------
# Prepare run directory
# -----------------------------------------------------------------------------
if [[ -z "${OUT_BASE}" ]]; then
  echo "❌ OUT_BASE is empty. Please set OUT_BASE to a valid directory path." >&2
  exit 1
fi
if [[ "${OUT_BASE}" == "/" ]]; then
  echo "❌ OUT_BASE is '/' — refusing to operate at root level." >&2
  exit 1
fi
mkdir -p "${OUT_BASE}" || { echo "❌ Failed to create OUT_BASE: ${OUT_BASE}" >&2; exit 1; }

TS="$(date '+%Y%m%d_%H%M%S')"
RUN_DIR="${OUT_BASE}/run_${TS}"
mkdir -p "${RUN_DIR}"/{logs,tmp,per_family,summary,inputs,global}
LOGFILE="${RUN_DIR}/logs/run.log"
: > "${LOGFILE}"

log "RUN_DIR=${RUN_DIR}"
log "COHORT_VCF=${COHORT_VCF}"

# =============================================================================
# Helper functions: set-based operations on site lists (CHROM POS TSVs)
# =============================================================================
normalize_sites(){
  local in_tsv="$1" out_tsv="$2"
  if [[ ! -s "${in_tsv}" ]]; then : > "${out_tsv}"; return 0; fi
  local tmp; tmp="$(mktemp)"
  sort -u -k1,1 -k2,2n "${in_tsv}" > "${tmp}"
  mv -f "${tmp}" "${out_tsv}"
}

intersect_sites(){
  local a_tsv="$1" b_tsv="$2" out_tsv="$3"
  if [[ ! -s "${a_tsv}" || ! -s "${b_tsv}" ]]; then : > "${out_tsv}"; return 0; fi
  awk 'BEGIN{FS=OFS="\t"} NR==FNR{seen[$1 FS $2]=1; next} {k=$1 FS $2; if(seen[k]) print $1,$2}' \
    "${a_tsv}" "${b_tsv}" > "${out_tsv}"
  normalize_sites "${out_tsv}" "${out_tsv}"
}

subtract_sites(){
  local a_tsv="$1" b_tsv="$2" out_tsv="$3"
  if [[ ! -s "${a_tsv}" ]]; then : > "${out_tsv}"; return 0; fi
  if [[ ! -s "${b_tsv}" ]]; then cp -f "${a_tsv}" "${out_tsv}"; return 0; fi
  awk 'BEGIN{FS=OFS="\t"} NR==FNR{ban[$1 FS $2]=1; next} {k=$1 FS $2; if(!ban[k]) print $1,$2}' \
    "${b_tsv}" "${a_tsv}" > "${out_tsv}"
  normalize_sites "${out_tsv}" "${out_tsv}"
}

sites_to_regions(){
  local sites_tsv="$1" regions_tsv="$2"
  if [[ -s "${sites_tsv}" ]]; then
    awk 'BEGIN{FS=OFS="\t"} {
      chrom=$1; pos=$2;
      start=pos-1;
      if (start<0) start=0;
      print chrom, start, pos
    }' "${sites_tsv}" > "${regions_tsv}"
  else
    : > "${regions_tsv}"
  fi
}

# =============================================================================
# make_present_absent_lists_het: strict HET/HOMREF logic implemented via bcftools query + AWK
# =============================================================================
make_present_absent_lists_het(){
  local in_vcf="$1"
  local samples_csv="$2"
  local out_present_tsv="$3"
  local out_absent_tsv="$4"

  : > "${out_present_tsv}"
  : > "${out_absent_tsv}"

  if [[ -z "${samples_csv}" ]]; then
    return 0
  fi

  "${BCT}" query \
    -s "${samples_csv}" \
    -f '%CHROM\t%POS[\t%GT\t%DP\t%GQ\t%AD]\n' \
    "${in_vcf}" \
  | awk -v min_dp_p="${MIN_DP_PRESENT}" \
        -v min_gq_p="${MIN_GQ_PRESENT}" \
        -v min_alt_p="${MIN_ALT_AD_PRESENT}" \
        -v min_ab_p="${MIN_AB_HET}" \
        -v max_ab_p="${MAX_AB_HET}" \
        -v min_dp_a="${MIN_DP_ABSENT}" \
        -v max_alt_a="${MAX_ALT_ADSENT:-${MAX_ALT_AD_ABSENT:-0}}" \
        -v max_ab_a="${MAX_AB_ABSENT}" '
      function is_het(gt){ return (gt=="0/1"||gt=="1/0"||gt=="0|1"||gt=="1|0") }
      function is_homref(gt){ return (gt=="0/0"||gt=="0|0") }
      function parse_ad(ad, a, n){
        REFAD=-1; ALTAD=-1
        if (ad=="" || ad==".") return 0
        n=split(ad,a,","); if(n<2) return 0
        if (a[1]==""||a[1]=="."||a[2]==""||a[2]==".") return 0
        REFAD=a[1]+0; ALTAD=a[2]+0
        return 1
      }
      {
        chrom=$1; pos=$2
        any_present=0
        all_absent=1
        for (i=3; i<=NF; i+=4) {
          gt=$(i); dp=$(i+1); gq=$(i+2); ad=$(i+3)
          if (dp==""||dp==".") dp=0
          if (gq==""||gq==".") gq=0
          dp=dp+0; gq=gq+0
          okad=parse_ad(ad)
          if (!okad) { all_absent=0; continue }
          tot=REFAD+ALTAD
          if (tot<=0) { all_absent=0; continue }
          ab=ALTAD/tot
          if (is_het(gt) && dp>=min_dp_p && (min_gq_p<=0 || gq>=min_gq_p) && ALTAD>=min_alt_p && ab>=min_ab_p && ab<=max_ab_p) {
            any_present=1
          }
          if (dp < min_dp_a) { all_absent=0; continue }
          if (!(is_homref(gt) && ALTAD<=max_alt_a && ab<=max_ab_a)) { all_absent=0 }
        }
        if (any_present==1) print chrom "\t" pos >> "'"${out_present_tsv}"'"
        if (all_absent==1)  print chrom "\t" pos >> "'"${out_absent_tsv}"'"
      }'
  normalize_sites "${out_present_tsv}" "${out_present_tsv}"
  normalize_sites "${out_absent_tsv}" "${out_absent_tsv}"
}

# =============================================================================
# build_global_background: GT-based prevalence across cohort
# =============================================================================
build_global_background(){
  local global_vcf="$1"
  local all_samples_csv="$2"
  local out_bg_sites="$3"

  : > "${out_bg_sites}"

  "${BCT}" query \
    -s "${all_samples_csv}" \
    -f '%CHROM\t%POS[\t%GT\t%DP]\n' \
    "${global_vcf}" \
  | awk -v mindp="${BG_MIN_DP}" -v k="${BG_MIN_SAMPLES_ALL}" '
      function is_nonref(gt){
        return (gt=="0/1"||gt=="1/0"||gt=="0|1"||gt=="1|0"||gt=="1/1"||gt=="1|1")
      }
      BEGIN{FS=OFS="\t"}
      {
        chrom=$1; pos=$2
        c=0
        for(i=3;i<=NF;i+=2){
          gt=$(i); dp=$(i+1)
          if(dp==""||dp==".") dp=0
          dp+=0
          if(dp < mindp) continue
          if(is_nonref(gt)) c++
        }
        if(c>=k) print chrom,pos
      }' > "${out_bg_sites}"

  normalize_sites "${out_bg_sites}" "${out_bg_sites}"
}

# =============================================================================
# filter_sites_unique_global_gt: ensure site not nonref in OTHER samples
# =============================================================================
filter_sites_unique_global_gt(){
  local global_vcf="$1"
  local all_samples_csv="$2"
  local focal_samples_csv="$3"
  local in_sites_tsv="$4"
  local out_sites_tsv="$5"
  local tmp_prefix="$6"

  : > "${out_sites_tsv}"
  [[ -s "${in_sites_tsv}" ]] || return 0

  local regions="${tmp_prefix}.regions.tsv"
  sites_to_regions "${in_sites_tsv}" "${regions}"
  [[ -s "${regions}" ]] || return 0

  "${BCT}" view -R "${regions}" "${global_vcf}" -Ou \
    | "${BCT}" query -s "${all_samples_csv}" -f '%CHROM\t%POS[\t%GT\t%DP]\n' - \
    | awk -v mindp="${UNIQUE_DP}" -v samples="${all_samples_csv}" -v focal="${focal_samples_csv}" '
        function is_nonref(gt){
          return (gt=="0/1"||gt=="1/0"||gt=="0|1"||gt=="1|0"||gt=="1/1"||gt=="1|1")
        }
        BEGIN{
          FS=OFS="\t"
          n=split(samples,S,",")
          m=split(focal,F,",")
          for(i=1;i<=m;i++) if(F[i]!="") focal_set[F[i]]=1
        }
        {
          chrom=$1; pos=$2
          other_has=0
          for(i=3;i<=NF;i+=2){
            idx=int((i-3)/2)+1
            sname=S[idx]
            if(sname in focal_set) continue
            gt=$(i); dp=$(i+1)
            if(dp==""||dp==".") dp=0
            dp+=0
            if(dp < mindp) continue
            if(is_nonref(gt)){ other_has=1; break }
          }
          if(other_has==0) print chrom,pos
        }' > "${out_sites_tsv}"

  normalize_sites "${out_sites_tsv}" "${out_sites_tsv}"
}

# =============================================================================
# emit_bucket_vcf: create per-bucket VCFs from site lists
# =============================================================================
emit_bucket_vcf(){
  local fam_vcf_gz="$1"
  local sites_tsv="$2"
  local out_vcf_gz="$3"
  local label="$4"
  local tmpdir="$5"

  local regions="${tmpdir}/${label}.regions.tsv"
  sites_to_regions "${sites_tsv}" "${regions}"

  if [[ ! -s "${regions}" ]]; then
    "${BCT}" view -h "${fam_vcf_gz}" | "${BCT}" view -Oz -o "${out_vcf_gz}" -
  else
    "${BCT}" view -R "${regions}" -Oz --threads "${BGZIP_THREADS}" -o "${out_vcf_gz}" "${fam_vcf_gz}"
  fi
  "${TABIX}" -f -p vcf "${out_vcf_gz}"
}

# =============================================================================
# run_vep_on_vcf: annotate VCF with VEP
# =============================================================================
run_vep_on_vcf(){
  local in_vcf="$1" out_vep_vcf_gz="$2" tmpdir="$3"

  mkdir -p "${tmpdir}"
  local tmp_out
  tmp_out="$(mktemp -p "${tmpdir}" --suffix=".vep.vcf")" || tmp_out="${tmpdir}/$(basename "${out_vep_vcf_gz%.gz}").tmp.vep.vcf"
  local tmp_log="${tmpdir}/$(basename "${out_vep_vcf_gz%.vcf.gz}").vep.log"

  log "    VEP: annotating ${in_vcf} -> ${out_vep_vcf_gz} (threads=${VEP_THREADS})"
  log "    VEP: cache=${VEP_CACHE_DIR:-(none)} extra_args='${VEP_EXTRA_ARGS}'"

  local vep_args=(--input_file "${in_vcf}" --output_file "${tmp_out}" --vcf --offline --cache --fork "${VEP_THREADS}" --force_overwrite)
  if [[ -n "${VEP_CACHE_DIR}" ]]; then
    vep_args+=( --dir_cache "${VEP_CACHE_DIR}" )
  fi
  if [[ -n "${VEP_SPECIES}" ]]; then
    vep_args+=( --species "${VEP_SPECIES}" )
  fi
  if [[ -n "${VEP_ASSEMBLY}" ]]; then
    vep_args+=( --assembly "${VEP_ASSEMBLY}" )
  fi
  if [[ -n "${VEP_EXTRA_ARGS}" ]]; then
    read -r -a extra <<< "${VEP_EXTRA_ARGS}"
    vep_args+=( "${extra[@]}" )
  fi

  if [[ -z "${VEP_CMD:-}" ]]; then
    log "    VEP: VEP_CMD not defined; skipping annotation"
    return 1
  fi
  if "${VEP_CMD}" "${vep_args[@]}" 2> "${tmp_log}"; then
    "${BGZIP}" -c "${tmp_out}" > "${out_vep_vcf_gz}"
    "${TABIX}" -f -p vcf "${out_vep_vcf_gz}"
    rm -f "${tmp_out}"
    log "    VEP: completed ${out_vep_vcf_gz} (log ${tmp_log})"
    return 0
  else
    local rc=$?
    log "    VEP: FAILED on ${in_vcf} (rc=${rc}). See ${tmp_log}"
    [[ -f "${tmp_out}" ]] && rm -f "${tmp_out}"
    return ${rc}
  fi
}

# =============================================================================
# Detect sample list, families, and controls
# =============================================================================
SAMPLES_TSV="${RUN_DIR}/inputs/cohort_samples.tsv"
"${BCT}" query -l "${COHORT_VCF}" | sort > "${SAMPLES_TSV}"
NSAMP="$(wc -l < "${SAMPLES_TSV}" | tr -d ' ')"
log "Cohort samples: ${NSAMP}"

WT_SAMPLES_FILE="${RUN_DIR}/inputs/wt_samples.txt"
PC_SAMPLES_FILE="${RUN_DIR}/inputs/pc_samples.txt"

grep -Ei '(^|[^A-Za-z0-9])wt([^A-Za-z0-9]|$)|(^wt[0-9]+)' "${SAMPLES_TSV}" | sort -u > "${WT_SAMPLES_FILE}" 2>/dev/null || :
if [[ ! -s "${WT_SAMPLES_FILE}" ]]; then
  grep -Ei 'wt' "${SAMPLES_TSV}" | sort -u > "${WT_SAMPLES_FILE}" 2>/dev/null || :
fi

grep -Ei '(^|[^A-Za-z0-9])pc([^A-Za-z0-9]|$)|(^pc[0-9]+)' "${SAMPLES_TSV}" | sort -u > "${PC_SAMPLES_FILE}" 2>/dev/null || :
if [[ ! -s "${PC_SAMPLES_FILE}" ]]; then
  grep -Ei 'pc' "${SAMPLES_TSV}" | sort -u > "${PC_SAMPLES_FILE}" 2>/dev/null || :
fi

WT_N="$(wc -l < "${WT_SAMPLES_FILE}" 2>/dev/null || echo 0 | tr -d ' ')" || WT_N=0
PC_N="$(wc -l < "${PC_SAMPLES_FILE}" 2>/dev/null || echo 0 | tr -d ' ')" || PC_N=0
log "Detected global controls: WT=${WT_N}, PC=${PC_N}"

# Build families list from sample names using token removal
FAMILIES_TSV="${RUN_DIR}/inputs/families.tsv"
awk -v IGNORECASE=1 -v throm_re="${THROM_RE}" -v nothrom_re="${NOTHROM_RE}" '
  {
    s=$0
    if (match(s, nothrom_re)) { fam=substr(s, 1, RSTART-1); if (fam!="") print fam; next }
    if (match(s, throm_re))   { fam=substr(s, 1, RSTART-1); if (fam!="") print fam; next }
  }
' "${SAMPLES_TSV}" | sort -u > "${FAMILIES_TSV}"

# Add wt/pc families if controls detected
if [[ -s "${WT_SAMPLES_FILE}" ]]; then echo "wt" >> "${FAMILIES_TSV}"; fi
if [[ -s "${PC_SAMPLES_FILE}" ]]; then echo "pc" >> "${FAMILIES_TSV}"; fi
sort -u "${FAMILIES_TSV}" -o "${FAMILIES_TSV}"

NFAM="$(wc -l < "${FAMILIES_TSV}" | tr -d ' ')"
log "Detected families: ${NFAM}"
[[ "${NFAM}" -gt 0 ]] || { log "❌ No families detected. Check regex / sample names."; exit 1; }

ALL_SAMPLES_CSV="$(paste -sd, "${SAMPLES_TSV}")"

# =============================================================================
# Build global filtered VCF once
# =============================================================================
GLOBAL_DIR="${RUN_DIR}/global"
mkdir -p "${GLOBAL_DIR}"/{vcf,lists,tmp}

GLOBAL_VCF="${GLOBAL_DIR}/vcf/global.filtered.vcf.gz"
SITE_ARGS=()
if [[ "${SITE_FILTER_ENABLE}" -eq 1 ]]; then
  [[ "${SITE_TYPES_ENABLE}" -eq 1 ]] && SITE_ARGS+=( -v "${SITE_TYPES}" )
  SITE_ARGS+=( -i "QUAL>=${SITE_MIN_QUAL}" )
  [[ "${SITE_REQUIRE_FILTER_PASS}" -eq 1 ]] && SITE_ARGS+=( -f "PASS" )
  [[ "${SITE_BIALLELIC_ONLY}" -eq 1 ]] && SITE_ARGS+=( -m2 -M2 )
fi

log "GLOBAL: building global filtered VCF once (SITE_ARGS: ${SITE_ARGS[*]:-(none)})"
"${BCT}" view "${SITE_ARGS[@]}" -s "${ALL_SAMPLES_CSV}" -Oz --threads "${BGZIP_THREADS}" -o "${GLOBAL_VCF}" "${COHORT_VCF}"
"${TABIX}" -f -p vcf "${GLOBAL_VCF}"
log "GLOBAL: global VCF records=$("${BCT}" view -H "${GLOBAL_VCF}" | wc -l | tr -d ' ')"

BG_SITES="${GLOBAL_DIR}/lists/background_sites.tsv"
if [[ "${BG_ENABLE}" -eq 1 ]]; then
  log "GLOBAL: building background sites = NONREF-by-GT in >=${BG_MIN_SAMPLES_ALL} samples (DP>=${BG_MIN_DP})"
  build_global_background "${GLOBAL_VCF}" "${ALL_SAMPLES_CSV}" "${BG_SITES}"
  log "GLOBAL: background_sites=$(wc -l < "${BG_SITES}" | tr -d ' ')"
else
  : > "${BG_SITES}"
  log "GLOBAL: BG_ENABLE=0 (no background subtraction)"
fi

# =============================================================================
# Prepare summary TSV header
# =============================================================================
SUMMARY_TSV="${RUN_DIR}/summary/family_counts.tsv"
echo -e "family\tthrom_samples\tnothrom_samples\tfam_vcf_records\tthrom_present\tthrom_absent\tnothrom_present\tnothrom_absent\tthrom_specific\tthrom_noBG\tthrom_unique\tnothrom_specific\tnothrom_noBG\tnothrom_unique" > "${SUMMARY_TSV}"

# =============================================================================
# Main per-family loop
# =============================================================================
while read -r FAM; do
  [[ -n "${FAM}" ]] || continue
  log "=== FAMILY ${FAM}"

  FAM_DIR="${RUN_DIR}/per_family/${FAM}"
  mkdir -p "${FAM_DIR}"/{inputs,lists,vcf,logs,tmp}
  mkdir -p "${FAM_DIR}/vcf/buckets"

  THROM_SAMPLES_FILE="${FAM_DIR}/inputs/throm_samples.txt"
  NOTHROM_SAMPLES_FILE="${FAM_DIR}/inputs/nothrom_samples.txt"

  if [[ "${FAM}" == "wt" ]]; then
    if [[ -s "${WT_SAMPLES_FILE}" ]]; then cp -f "${WT_SAMPLES_FILE}" "${NOTHROM_SAMPLES_FILE}"; else : > "${NOTHROM_SAMPLES_FILE}"; fi
    : > "${THROM_SAMPLES_FILE}"
  elif [[ "${FAM}" == "pc" ]]; then
    if [[ -s "${PC_SAMPLES_FILE}" ]]; then cp -f "${PC_SAMPLES_FILE}" "${THROM_SAMPLES_FILE}"; else : > "${THROM_SAMPLES_FILE}"; fi
    : > "${NOTHROM_SAMPLES_FILE}"
  else
    awk -v IGNORECASE=1 -v fam="${FAM}" -v nothrom_re="${NOTHROM_RE}" \
      '$0 ~ ("^" fam) && $0 ~ nothrom_re {print}' "${SAMPLES_TSV}" | sort -u > "${NOTHROM_SAMPLES_FILE}"
    awk -v IGNORECASE=1 -v fam="${FAM}" -v throm_re="${THROM_RE}" -v nothrom_re="${NOTHROM_RE}" \
      '$0 ~ ("^" fam) && $0 ~ throm_re && $0 !~ nothrom_re {print}' "${SAMPLES_TSV}" | sort -u > "${THROM_SAMPLES_FILE}"
  fi

  if [[ "${FAM}" != "wt" && -s "${WT_SAMPLES_FILE}" ]]; then cat "${WT_SAMPLES_FILE}" >> "${NOTHROM_SAMPLES_FILE}"; fi
  if [[ "${FAM}" != "pc" && -s "${PC_SAMPLES_FILE}" ]]; then cat "${PC_SAMPLES_FILE}" >> "${THROM_SAMPLES_FILE}"; fi

  if [[ -s "${NOTHROM_SAMPLES_FILE}" ]]; then sort -u "${NOTHROM_SAMPLES_FILE}" -o "${NOTHROM_SAMPLES_FILE}"; fi
  if [[ -s "${THROM_SAMPLES_FILE}" ]]; then sort -u "${THROM_SAMPLES_FILE}" -o "${THROM_SAMPLES_FILE}"; fi

  N_THROM="$(wc -l < "${THROM_SAMPLES_FILE}" 2>/dev/null || echo 0 | tr -d ' ')" || N_THROM=0
  N_NOTHROM="$(wc -l < "${NOTHROM_SAMPLES_FILE}" 2>/dev/null || echo 0 | tr -d ' ')" || N_NOTHROM=0
  log "  THROM samples:   ${N_THROM}"
  log "  NOTHROM samples: ${N_NOTHROM}"

  THROM_CSV="$(paste -sd, "${THROM_SAMPLES_FILE}" 2>/dev/null || true)"
  NOTHROM_CSV="$(paste -sd, "${NOTHROM_SAMPLES_FILE}" 2>/dev/null || true)"
  if [[ -n "${THROM_CSV}" && -n "${NOTHROM_CSV}" ]]; then BOTH_CSV="${THROM_CSV},${NOTHROM_CSV}"; elif [[ -n "${THROM_CSV}" ]]; then BOTH_CSV="${THROM_CSV}"; else BOTH_CSV="${NOTHROM_CSV}"; fi

  if [[ "${N_THROM}" -eq 0 && "${N_NOTHROM}" -eq 0 ]]; then
    log "  ⚠️  Skipping ${FAM} (no samples in throm nor nothrom/control arms)"
    echo -e "${FAM}\t${N_THROM}\t${N_NOTHROM}\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0" >> "${SUMMARY_TSV}"
    continue
  fi

  FAM_VCF="${FAM_DIR}/vcf/${FAM}.family.vcf.gz"
  log "  Building family VCF (SITE_ARGS: ${SITE_ARGS[*]:-(none)})"
  "${BCT}" view "${SITE_ARGS[@]}" -s "${BOTH_CSV}" -Oz --threads "${BGZIP_THREADS}" -o "${FAM_VCF}" "${COHORT_VCF}"
  "${TABIX}" -f -p vcf "${FAM_VCF}"

  FAM_NREC="$("${BCT}" view -H "${FAM_VCF}" | wc -l | tr -d ' ')"
  log "  Family VCF records: ${FAM_NREC}"
  [[ "${FAM_NREC}" -gt 0 ]] || { log "❌ Family VCF has ZERO records. Loosen SITE_* knobs first."; exit 1; }

  THROM_PRESENT="${FAM_DIR}/lists/throm.present.tsv"
  THROM_ABSENT="${FAM_DIR}/lists/throm.absent.tsv"
  NOTHROM_PRESENT="${FAM_DIR}/lists/nothrom.present.tsv"
  NOTHROM_ABSENT="${FAM_DIR}/lists/nothrom.absent.tsv"

  log "  Building THROM present/absent (STRICT HET logic)..."
  if [[ -n "${THROM_CSV}" ]]; then
    make_present_absent_lists_het "${FAM_VCF}" "${THROM_CSV}" "${THROM_PRESENT}" "${THROM_ABSENT}"
  else
    : > "${THROM_PRESENT}"; : > "${THROM_ABSENT}"
  fi

  log "  Building NOTHROM present/absent (STRICT HET logic)..."
  if [[ -n "${NOTHROM_CSV}" ]]; then
    make_present_absent_lists_het "${FAM_VCF}" "${NOTHROM_CSV}" "${NOTHROM_PRESENT}" "${NOTHROM_ABSENT}"
  else
    : > "${NOTHROM_PRESENT}"; : > "${NOTHROM_ABSENT}"
  fi

  THROM_SPEC="${FAM_DIR}/lists/throm_specific.tsv"
  NOTHROM_SPEC="${FAM_DIR}/lists/nothrom_specific.tsv"
  intersect_sites "${THROM_PRESENT}"   "${NOTHROM_ABSENT}" "${THROM_SPEC}"
  intersect_sites "${NOTHROM_PRESENT}" "${THROM_ABSENT}"   "${NOTHROM_SPEC}"

  THROM_NOBG="${FAM_DIR}/lists/throm_specific_noBG.tsv"
  NOTHROM_NOBG="${FAM_DIR}/lists/nothrom_specific_noBG.tsv"
  subtract_sites "${THROM_SPEC}"   "${BG_SITES}" "${THROM_NOBG}"
  subtract_sites "${NOTHROM_SPEC}" "${BG_SITES}" "${NOTHROM_NOBG}"

  THROM_UNI="${FAM_DIR}/lists/throm_unique.tsv"
  NOTHROM_UNI="${FAM_DIR}/lists/nothrom_unique.tsv"
  if [[ "${UNIQUE_ENABLE}" -eq 1 ]]; then
    log "  Global-unique filtering (GT-based) against OTHER samples (DP>=${UNIQUE_DP})..."
    filter_sites_unique_global_gt "${GLOBAL_VCF}" "${ALL_SAMPLES_CSV}" "${BOTH_CSV}" "${THROM_NOBG}" "${THROM_UNI}" "${FAM_DIR}/tmp/throm_uni"
    filter_sites_unique_global_gt "${GLOBAL_VCF}" "${ALL_SAMPLES_CSV}" "${BOTH_CSV}" "${NOTHROM_NOBG}" "${NOTHROM_UNI}" "${FAM_DIR}/tmp/nothrom_uni"
  else
    cp -f "${THROM_NOBG}" "${THROM_UNI}"
    cp -f "${NOTHROM_NOBG}" "${NOTHROM_UNI}"
  fi

  C_TP="$(wc -l < "${THROM_PRESENT}" | tr -d ' ')"
  C_TA="$(wc -l < "${THROM_ABSENT}" | tr -d ' ')"
  C_NP="$(wc -l < "${NOTHROM_PRESENT}" | tr -d ' ')"
  C_NA="$(wc -l < "${NOTHROM_ABSENT}" | tr -d ' ')"
  C_TS="$(wc -l < "${THROM_SPEC}" | tr -d ' ')"
  C_TB="$(wc -l < "${THROM_NOBG}" | tr -d ' ')"
  C_TU="$(wc -l < "${THROM_UNI}" | tr -d ' ')"
  C_NS="$(wc -l < "${NOTHROM_SPEC}" | tr -d ' ')"
  C_NB="$(wc -l < "${NOTHROM_NOBG}" | tr -d ' ')"
  C_NU="$(wc -l < "${NOTHROM_UNI}" | tr -d ' ')"

  log "  Evidence counts:"
  log "    THROM: present=${C_TP} absent=${C_TA} specific=${C_TS} noBG=${C_TB} unique=${C_TU}"
  log "    NOTHROM: present=${C_NP} absent=${C_NA} specific=${C_NS} noBG=${C_NB} unique=${C_NU}"

  echo -e "${FAM}\t${N_THROM}\t${N_NOTHROM}\t${FAM_NREC}\t${C_TP}\t${C_TA}\t${C_NP}\t${C_NA}\t${C_TS}\t${C_TB}\t${C_TU}\t${C_NS}\t${C_NB}\t${C_NU}" >> "${SUMMARY_TSV}"

  if [[ "${EMIT_BUCKET_VCFS}" -eq 1 ]]; then
    log "  Emitting per-bucket VCFs (no snpEff)..."

    declare -a BUCKETS=(
      "throm_present|${THROM_PRESENT}"
      "throm_absent|${THROM_ABSENT}"
      "nothrom_present|${NOTHROM_PRESENT}"
      "nothrom_absent|${NOTHROM_ABSENT}"
      "throm_specific|${THROM_SPEC}"
      "nothrom_specific|${NOTHROM_SPEC}"
      "throm_specific_noBG|${THROM_NOBG}"
      "nothrom_specific_noBG|${NOTHROM_NOBG}"
      "throm_unique|${THROM_UNI}"
      "nothrom_unique|${NOTHROM_UNI}"
    )

    for entry in "${BUCKETS[@]}"; do
      label="${entry%%|*}"
      sites="${entry#*|}"

      out_vcf="${FAM_DIR}/vcf/buckets/${FAM}.${label}.vcf.gz"
      emit_bucket_vcf "${FAM_VCF}" "${sites}" "${out_vcf}" "${label}" "${FAM_DIR}/tmp"

      nrec="$("${BCT}" view -H "${out_vcf}" | wc -l | tr -d ' ')"
      log "    bucket=${label} records=${nrec}"

      if [[ "${VEP_ENABLE}" -eq 1 ]]; then
        if [[ "${nrec}" -gt 0 ]]; then
          out_vcf_vep="${out_vcf%.vcf.gz}.vep.vcf.gz"
          if run_vep_on_vcf "${out_vcf}" "${out_vcf_vep}" "${FAM_DIR}/tmp"; then
            log "    bucket=${label} VEP annotated: ${out_vcf_vep}"
          else
            log "    bucket=${label} VEP annotation FAILED (see per-bucket log in ${FAM_DIR}/tmp)"
          fi
        else
          log "    bucket=${label} empty — skipping VEP"
        fi
      fi
    done
  fi

done < "${FAMILIES_TSV}"

log "DONE"
log "Summary TSV: ${SUMMARY_TSV}"
log "Background sites: ${BG_SITES}"
log "Run dir: ${RUN_DIR}"
log "Log file: ${LOGFILE}"

# End of script

: <<'END_NOTES'
Troubleshooting checklist and suggestions:
 - Verify COHORT_VCF contains FORMAT/AD (allelic depths). If not, re-run variant caller with AD output.
 - If family detection fails, set THROM_RE / NOTHROM_RE to patterns that match your sample names, or supply families manually.
 - For multi-allelic sites the AWK AD parsing may not capture the intended ALT index; consider preprocessing to split multi-allelic sites.
 - To speed up runs, parallelize per-family iterations (GNU parallel) or run per-family scripts on a cluster.
 - To change strictness, adjust MIN_DP_*/MIN_GQ_*/MIN_ALT_AD_* and AB thresholds.
 - If you want, I can: (a) convert AWK parsing into Python for clearer multi-allelic support, (b) add an option to output one-row-per-allele, or (c) add a dry-run mode that prints family sample lists and counts only.

END_NOTES
