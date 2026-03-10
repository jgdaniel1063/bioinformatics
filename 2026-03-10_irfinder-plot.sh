#!/usr/bin/env bash
#
# 2026-02-03_irfinder_pseudorep_logit_pairwise_diffIR_FULL_annotated.sh
#
# Extremely heavily annotated version of the original pipeline that:
#  - performs per-family paired Δlogit(IR) tests using IRFinder outputs,
#  - estimates dispersion by leave-one-out from other families (robust MAD),
#  - performs Stouffer meta-analysis across families,
#  - writes per-family and meta results and some diagnostic plots.
#
# This annotation explains:
#  - the intent of each block of code,
#  - expected file layouts and naming conventions,
#  - failure modes and how to debug them,
#  - suggestions for improvements and reproducibility,
#  - and where to edit paths / tokens to adapt to your environment.
#
# Run (example):
#   bash 2026-02-03_irfinder_pseudorep_logit_pairwise_diffIR_FULL_annotated.sh
#
# IMPORTANT NOTES / ASSUMPTIONS
# ---------------------------
# - This script is written for interactive Linux usage (no Slurm). It expects Rscript and R packages
#   data.table, ggplot2 to be available in the environment where Rscript runs.
# - IRFinder outputs must be structured as:
#       ${IRFINDER_ROOT}/samples/<sample_name>/IRFinder-IR-nondir.txt
#   The R code attempts to auto-detect the IR ratio column and the intron identifier column.
# - Family/sample naming: the script tries to detect throm vs nothrom sample folders by regexes
#   applied to the folder names under SAMPLES_DIR. If detection fails, edit THROM_TOKEN_REGEX /
#   NOTHROM_TOKEN_REGEX or provide samples manually (script would need modest adaptation).
# - The pipeline uses a pseudoreplication approach where each family contributes a paired delta:
#       Δlogit(IR) = logit(IR_throm) - logit(IR_nothrom)
#   and variance (sigma) for each intron is estimated from other families' deltas using a robust MAD.
#
# OUTPUT
# ------
# Under RUN_DIR (timestamped) you will find:
#   pairs.tsv                          : family -> throm_sample, nothrom_sample
#   per_family/diffIR_<family>.tsv.gz  : per-family test table
#   meta/meta_stouffer.tsv.gz          : combined meta-analysis table (Stouffer)
#   plots/volcano_<family>.png         : per-family volcano plots
#   plots/volcano_META.png             : meta volcano plot
#   run logs                           : Rscript stdout/log capture
#   LATEST -> RUN_DIR                  : symlink in OUT_BASE pointing to latest run
#
# If you want help interpreting results or errors, paste the Rscript stdout log and the tail of the Rscript-written files.
#
set -euo pipefail
shopt -s nullglob

# --------------------------------------------------------------------
# HARD-CODED PATHS (EDIT THESE BEFORE RUNNING)
# --------------------------------------------------------------------
# IRFINDER_ROOT must contain a 'samples' subdirectory with per-sample IRFinder outputs:
#   ${IRFINDER_ROOT}/samples/<sample_name>/IRFinder-IR-nondir.txt
IRFINDER_ROOT="/home/jgd/PROJECT/irfinder_counts/run_20260113_120037_39787519"
SAMPLES_DIR="${IRFINDER_ROOT}/samples"

# Base directory where timestamped runs will be created
OUT_BASE="/home/jgd/PROJECT/irfinder_diff_pseudorep_logit"
TS="$(date +%Y%m%d_%H%M%S)"
RUN_DIR="${OUT_BASE}/run_${TS}"

# Optionally run R from an activated conda env: set the env name here (leave empty to use system R)
CONDA_ENV_R=""   # e.g. "r_env" or "irfinder_env"

# Tokens used to detect throm / nothrom sample folders from their names
# These are regular expressions used with grep -Ei (case-insensitive)
THROM_TOKEN_REGEX='throm|thrombosis'
NOTHROM_TOKEN_REGEX='nothrom|no[_-]?throm|no[_-]?thrombosis|nonthrom|not[_-]?thrombosis'

# Families (pseudoreplicates) listed explicitly. The script will look for sample folders starting with these tokens.
FAMILIES=(ai08 au01 bo01 bq01 bq06 bqm2)
# --------------------------------------------------------------------

# Create output directories and a LATEST symlink for convenience. The directories under RUN_DIR:
#   logs/, per_family/, meta/, plots/, tmp/
mkdir -p "${RUN_DIR}"/{logs,per_family,meta,plots,tmp}
mkdir -p "${OUT_BASE}"
ln -sfn "${RUN_DIR}" "${OUT_BASE}/LATEST"

# --------------------------------------------------------------------
# Conda activation helper (optional)
# --------------------------------------------------------------------
# If you want to run R inside a specific conda environment, set CONDA_ENV_R above.
# This function finds the user's conda installation and sources the activation script,
# then activates the named environment. If conda is missing or activate fails, the script aborts.
conda_activate_if_set() {
  if [[ -n "${CONDA_ENV_R}" ]]; then
    if [[ -f "${HOME}/miniconda3/etc/profile.d/conda.sh" ]]; then
      # shellcheck disable=SC1091
      source "${HOME}/miniconda3/etc/profile.d/conda.sh"
    elif [[ -f "${HOME}/anaconda3/etc/profile.d/conda.sh" ]]; then
      # shellcheck disable=SC1091
      source "${HOME}/anaconda3/etc/profile.d/conda.sh"
    else
      echo "ERROR: CONDA_ENV_R is set but conda.sh not found under ~/miniconda3 or ~/anaconda3" >&2
      exit 2
    fi
    conda activate "${CONDA_ENV_R}" >/dev/null 2>&1 || { echo "ERROR: failed to conda activate ${CONDA_ENV_R}" >&2; exit 2; }
  fi
}

# --------------------------------------------------------------------
# Build pairs.tsv automatically by looking for throm/nothrom sample folders
# --------------------------------------------------------------------
# Format of pairs.tsv:
#   family <TAB> throm_sample <TAB> nothrom_sample
#
# This auto-detection is convenient but fragile if your folder names do not follow the expected patterns.
# If auto-detection fails, adjust FAMILIES / THROM_TOKEN_REGEX / NOTHROM_TOKEN_REGEX, or create pairs.tsv manually.
PAIRS_TSV="${RUN_DIR}/pairs.tsv"
{
  echo -e "family\tthrom_sample\tnothrom_sample"
  for fam in "${FAMILIES[@]}"; do
    # We search folder names that start with the family token then contain a throm or nothrom token
    throm="$(ls -1 "${SAMPLES_DIR}" 2>/dev/null | grep -Ei "^${fam}.*(${THROM_TOKEN_REGEX})" | head -n 1 || true)"
    nothrom="$(ls -1 "${SAMPLES_DIR}" 2>/dev/null | grep -Ei "^${fam}.*(${NOTHROM_TOKEN_REGEX})" | head -n 1 || true)"
    if [[ -z "${throm}" || -z "${nothrom}" ]]; then
      echo "ERROR: could not auto-detect throm/nothrom sample folders for family '${fam}' in ${SAMPLES_DIR}" >&2
      echo "  throm='${throm}' nothrom='${nothrom}'" >&2
      echo "  Inspect ${SAMPLES_DIR} and adjust regex or FAMILIES." >&2
      exit 2
    fi
    echo -e "${fam}\t${throm}\t${nothrom}"
  done
} > "${PAIRS_TSV}"

# Annotated note:
# - pairs.tsv will be used by the R script below to decide which sample pairs to compare.
# - If your sample names include spaces or unusual characters, the grep may fail; sanitize folder names first.

# --------------------------------------------------------------------
# Create the R script that performs the core differential IR analysis
# We write a fully self-contained R script into RUN_DIR, then run it.
# --------------------------------------------------------------------
R_SCRIPT="${RUN_DIR}/diffIR_pseudorep_logit_pairwise.R"
cat > "${R_SCRIPT}" <<'RSCRIPT'
# ---------------------------
# R: diffIR_pseudorep_logit_pairwise.R
# ---------------------------
# The embedded R script is heavily commented inline; it:
#  - reads pairs.tsv and IRFinder per-sample IR files,
#  - auto-detects IR column and intron identifier,
#  - computes per-family Δlogit(IR),
#  - estimates per-intron sigma by leave-one-out robust MAD,
#  - computes Z-scores, p-values and BH-adjusted q-values per family,
#  - performs Stouffer meta-analysis across families,
#  - writes per-family and meta results and generates simple volcano plots.
#
suppressPackageStartupMessages({
  library(data.table)   # fast table IO and operations
  library(ggplot2)      # plotting
})

# Command-line args: SAMPLES_DIR, PAIRS_TSV, RUN_DIR
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) stop("Usage: Rscript diffIR_pseudorep_logit_pairwise.R <SAMPLES_DIR> <PAIRS_TSV> <RUN_DIR>")
SAMPLES_DIR <- normalizePath(args[[1]], mustWork=TRUE)
PAIRS_TSV   <- normalizePath(args[[2]], mustWork=TRUE)
RUN_DIR     <- normalizePath(args[[3]], mustWork=TRUE)

# Prepare output dirs under RUN_DIR (the shell script already created them, but ensure they exist)
dir.create(file.path(RUN_DIR, "per_family"), showWarnings=FALSE, recursive=TRUE)
dir.create(file.path(RUN_DIR, "meta"),       showWarnings=FALSE, recursive=TRUE)
dir.create(file.path(RUN_DIR, "plots"),      showWarnings=FALSE, recursive=TRUE)

pairs <- fread(PAIRS_TSV)
stopifnot(all(c("family","throm_sample","nothrom_sample") %in% names(pairs)))

# ----------------------------
# Helper functions (documented)
# ----------------------------
# logit: transform probability to log-odds
logit <- function(p) log(p/(1-p))

# clamp01: clamp probabilities into (eps, 1-eps) interval to avoid Inf in logit
clamp01 <- function(x, eps=1e-6) pmin(pmax(x, eps), 1-eps)

# robust_sigma: robust estimate of spread for a numeric vector
#   uses MAD scaled to approximate sigma (1.4826*mad),
#   falls back to sd or floor if MAD is not available.
robust_sigma <- function(x, floor=0.10) {
  x <- x[is.finite(x)]
  if (length(x) < 2) return(floor)
  s <- 1.4826 * median(abs(x - median(x)))
  if (!is.finite(s) || s <= 0) s <- sd(x)
  if (!is.finite(s) || s <= 0) s <- floor
  max(s, floor)
}

# pick_ir_col: heuristics to find the IR column name in IRFinder's output (many possible names)
pick_ir_col <- function(dt) {
  cn <- names(dt)
  # common names include IRratio, IRRatio, IR_ratio, irratio, etc.
  cand <- cn[grepl("^IRratio$|^IRRatio$|^IR_ratio$|^irratio$", cn, ignore.case=TRUE)]
  if (length(cand) > 0) return(cand[[1]])
  cand <- cn[grepl("ir.*ratio|ratio.*ir", cn, ignore.case=TRUE)]
  if (length(cand) > 0) return(cand[[1]])
  stop("Could not find IR ratio column in IRFinder file.")
}

# make_intron_id: create a reproducible intron identifier from columns (prefer explicit ID columns)
make_intron_id <- function(dt) {
  cn <- names(dt)
  if ("Intron" %in% cn) return(as.character(dt[["Intron"]]))
  if ("intron" %in% cn) return(as.character(dt[["intron"]]))
  if ("Name" %in% cn)   return(as.character(dt[["Name"]]))
  if ("name" %in% cn)   return(as.character(dt[["name"]]))

  # otherwise attempt Chr:Start-End[:Strand]
  chr_col <- cn[grepl("^Chr$|^CHROM$|^Chrom$|^chrom$|^chr$", cn)]
  st_col  <- cn[grepl("^Start$|^POS$|^pos$|^start$", cn)]
  en_col  <- cn[grepl("^End$|^end$", cn)]
  str_col <- cn[grepl("^Strand$|^strand$", cn)]

  if (length(chr_col)>0 && length(st_col)>0 && length(en_col)>0) {
    chr <- as.character(dt[[chr_col[[1]]]])
    st  <- as.integer(dt[[st_col[[1]]]])
    en  <- as.integer(dt[[en_col[[1]]]])
    if (length(str_col)>0) {
      str <- as.character(dt[[str_col[[1]]]])
      return(paste0(chr, ":", st, "-", en, ":", str))
    } else {
      return(paste0(chr, ":", st, "-", en))
    }
  }
  stop("Could not build intron_id (no Intron/Name or coordinate columns found).")
}

# read_irfinder_ir: read the IRFinder-IR-nondir.txt file for one sample and return a 2-column data.table
#   with columns intron_id and IR (raw numeric ratio)
read_irfinder_ir <- function(sample_name) {
  f <- file.path(SAMPLES_DIR, sample_name, "IRFinder-IR-nondir.txt")
  if (!file.exists(f)) stop(paste0("Missing file: ", f))
  dt <- fread(f)
  intron_id <- make_intron_id(dt)
  ir_col <- pick_ir_col(dt)
  ir <- as.numeric(dt[[ir_col]])
  data.table(intron_id=intron_id, IR=ir)
}

# ----------------------------
# Load all samples involved in pairs
# ----------------------------
all_samples <- unique(c(pairs$throm_sample, pairs$nothrom_sample))
ir_list <- setNames(vector("list", length(all_samples)), all_samples)

for (s in all_samples) {
  ir_list[[s]] <- read_irfinder_ir(s)
}

# Merge all per-sample IRs into a wide table keyed on intron_id
wide <- Reduce(function(a,b){
  merge(a, b, by="intron_id", all=TRUE)
}, lapply(names(ir_list), function(s){
  dt <- copy(ir_list[[s]])
  setnames(dt, "IR", s)
  dt
}))

# Clamp IR values into (eps, 1-eps) and compute logit per sample
eps <- 1e-6
for (s in all_samples) {
  wide[[s]] <- clamp01(as.numeric(wide[[s]]), eps=eps)
  wide[[paste0("L_", s)]] <- logit(wide[[s]])
}

# ----------------------------
# Compute per-family deltas: DELTA_<family> = L_throm - L_nothrom
# ----------------------------
delta_cols <- character(nrow(pairs))
for (i in seq_len(nrow(pairs))) {
  fam <- pairs$family[[i]]
  tS  <- pairs$throm_sample[[i]]
  nS  <- pairs$nothrom_sample[[i]]
  col <- paste0("DELTA_", fam)
  wide[[col]] <- wide[[paste0("L_", tS)]] - wide[[paste0("L_", nS)]]
  delta_cols[[i]] <- col
}

# ----------------------------
# Per-family tests: for each family compute sigma from other families (leave-one-out),
# then compute z, p, q for that family's delta.
# ----------------------------
for (i in seq_len(nrow(pairs))) {
  fam <- pairs$family[[i]]
  tS  <- pairs$throm_sample[[i]]
  nS  <- pairs$nothrom_sample[[i]]
  dcol <- paste0("DELTA_", fam)
  other_dcols <- setdiff(delta_cols, dcol)

  # sigma: for every intron compute robust sigma across other deltas (leave-one-out)
  sig <- wide[, {
    dvals <- unlist(.SD, use.names=FALSE)
    .(sigma = robust_sigma(dvals, floor=0.10))
  }, by=intron_id, .SDcols=other_dcols]

  tmp <- merge(
    wide[, .(
      intron_id,
      IR_throm       = get(tS),
      IR_nothrom     = get(nS),
      logit_throm    = get(paste0("L_", tS)),
      logit_nothrom  = get(paste0("L_", nS)),
      delta          = get(dcol)
    )],
    sig, by="intron_id", all.x=TRUE
  )

  # Z-score and p/q values under normal approx (two-sided)
  tmp[, z := delta / sigma]
  tmp[, p := 2*pnorm(-abs(z))]
  tmp[, q := p.adjust(p, method="BH")]

  out_path <- file.path(RUN_DIR, "per_family", paste0("diffIR_", fam, ".tsv.gz"))
  fwrite(tmp[order(p)], out_path, sep="\t")

  # Simple volcano plot: Δlogit vs -log10(p)
  pvol <- ggplot(tmp, aes(x=delta, y=-log10(p))) +
    geom_point(size=0.6) +
    labs(title=paste0("Δlogit(IR) volcano: ", fam),
         x="Δlogit(IR) (throm - nothrom)",
         y="-log10(p)") +
    theme_bw(base_size=12)

  ggsave(filename=file.path(RUN_DIR, "plots", paste0("volcano_", fam, ".png")),
         plot=pvol, width=7, height=5, dpi=200)
}

# ----------------------------
# Meta-analysis (Stouffer) across families
# ----------------------------
# Compute a global robust sigma across all delta columns, then compute per-family z and combine via Stouffer.
sig_global <- wide[, {
  dvals <- unlist(.SD, use.names=FALSE)
  .(sigma_global = robust_sigma(dvals, floor=0.10))
}, by=intron_id, .SDcols=delta_cols]

meta <- merge(wide[, .(intron_id)], sig_global, by="intron_id", all.x=TRUE)

for (i in seq_len(nrow(pairs))) {
  fam <- pairs$family[[i]]
  dcol <- paste0("DELTA_", fam)
  meta[[paste0("z_", fam)]] <- wide[[dcol]] / meta[["sigma_global"]]
  meta[[paste0("delta_", fam)]] <- wide[[dcol]]
}

zcols <- paste0("z_", pairs$family)
dcols <- paste0("delta_", pairs$family)

# k = number of finite z-values available for the intron (some families missing data cause NAs)
meta[, k := rowSums(is.finite(.SD)), .SDcols=zcols]
meta[, z_sum := rowSums(.SD, na.rm=TRUE), .SDcols=zcols]
meta[, z_meta := z_sum / sqrt(pmax(k, 1))]
meta[, p_meta := 2*pnorm(-abs(z_meta))]
meta[, q_meta := p.adjust(p_meta, method="BH")]

# sign_consistency counts how many families agree on sign
meta[, sign_consistency := {
  zs <- unlist(.SD, use.names=FALSE)
  zs <- zs[is.finite(zs)]
  if (length(zs) == 0) NA_integer_ else max(sum(zs > 0), sum(zs < 0))
}, by=intron_id, .SDcols=zcols]

meta[, delta_mean := rowMeans(.SD, na.rm=TRUE), .SDcols=dcols]

meta_out <- meta[order(p_meta)]
meta_path <- file.path(RUN_DIR, "meta", "meta_stouffer.tsv.gz")
fwrite(meta_out, meta_path, sep="\t")

# Meta volcano
pmeta <- ggplot(meta_out, aes(x=delta_mean, y=-log10(p_meta))) +
  geom_point(size=0.6) +
  labs(title="Meta (Stouffer) Δlogit(IR) volcano",
       x="Mean Δlogit(IR) across families",
       y="-log10(meta p)") +
  theme_bw(base_size=12)

ggsave(filename=file.path(RUN_DIR, "plots", "volcano_META.png"),
       plot=pmeta, width=7, height=5, dpi=200)

# Save top hits (e.g., top 200)
topN <- 200
top_hits <- meta_out[1:min(topN, nrow(meta_out)),
                     .(intron_id, delta_mean, z_meta, p_meta, q_meta, sign_consistency, k)]
fwrite(top_hits, file.path(RUN_DIR, "meta", "top_hits.tsv"), sep="\t")

cat("OK\n")
cat("Per-family: ", file.path(RUN_DIR, "per_family"), "\n", sep="")
cat("Meta:       ", meta_path, "\n", sep="")
cat("Plots:      ", file.path(RUN_DIR, "plots"), "\n", sep="")
RSCRIPT

# --------------------------------------------------------------------
# Execute the R script
# --------------------------------------------------------------------
# Optionally activate a conda R environment before running Rscript
conda_activate_if_set

# Capture Rscript version used and then run the analysis script, teeing stdout to a log
Rscript --version >"${RUN_DIR}/logs/Rscript.version.txt" 2>&1 || true

# Run the R script. The R script writes its own files under RUN_DIR.
Rscript "${R_SCRIPT}" "${SAMPLES_DIR}" "${PAIRS_TSV}" "${RUN_DIR}" \
  | tee "${RUN_DIR}/logs/Rscript.stdout.log"

echo
echo "[DONE] RUN_DIR=${RUN_DIR}"
echo "[DONE] LATEST -> ${OUT_BASE}/LATEST"
echo "[DONE] pairs.tsv=${PAIRS_TSV}"
echo "[DONE] See per-family results under: ${RUN_DIR}/per_family"
echo "[DONE] See meta results under: ${RUN_DIR}/meta"
echo "[DONE] Plots under: ${RUN_DIR}/plots"
echo

# --------------------------------------------------------------------
# VERY DETAILED ANNOTATION, TROUBLESHOOTING & SUGGESTIONS
# --------------------------------------------------------------------
: <<'END_NOTES'
This section contains extended guidance for debugging, validating results,
and improving the pipeline. It does not execute as part of the script.

1) Input file expectations (IRFinder)
   - For each sample, the script expects:
       ${SAMPLES_DIR}/${sample_name}/IRFinder-IR-nondir.txt
     This file is created by IRFinder's 'IRFinder -i' output (the IR-nondir set).
   - The R code tries to auto-detect the IR ratio column (e.g., IRratio). If your file
     uses a different header, inspect the first few lines:
       head -n 20 ${SAMPLES_DIR}/${sample_name}/IRFinder-IR-nondir.txt

2) Introns identification
   - The script constructs intron identifiers using an explicit 'Intron' or 'Name' column
     if present; otherwise it constructs "Chr:Start-End[:Strand]".
   - Ensure coordinate columns exist in the IRFinder output; if not, modify make_intron_id()
     in the R script to match your file column names.

3) Family/sample detection (pairs.tsv)
   - The shell code uses grep regex on folder names to auto-detect throm / nothrom samples.
   - If your sample folder names do not follow the assumed convention, either:
     * update THROM_TOKEN_REGEX / NOTHROM_TOKEN_REGEX, or
     * create a manual pairs.tsv and skip auto-build (you can replace PAIRS_TSV with your file).

4) Statistical model & assumptions
   - Δlogit(IR) computed per intron per family is treated as approximately normal with variance sigma^2.
   - Sigma for each intron is estimated robustly from the other families' deltas (leave-one-out),
     using a MAD-based estimator scaled to approximate standard deviation (1.4826*mad).
   - If family sample sizes are extremely small or a family has missing values for many introns,
     sigma estimates may be noisy. The code floors sigma at 0.10 to avoid excessive z-values.
   - The pipeline uses two-sided z-tests and BH correction per family, then Stouffer to combine z-scores.
   - This is an approximate approach; for full statistical rigor consider mixed-effect models
     or hierarchical Bayesian models to share information across families.

5) Common failure modes
   - Missing files: the R script will error if an expected IRFinder output is not present.
     Fix: ensure SAMPLES_DIR is correct and that IRFinder outputs were generated.
   - IR column missing: pick_ir_col() may fail if IRFinder output format differs. Inspect column names
     and adapt the regex or edit the function to pick the correct column.
   - Intron id not constructed: make_intron_id() tries multiple common names. If your output uses
     other column names, add them to the detection logic.
   - All deltas NA or zero: may happen if IR values are exactly 0/1 or if clamp01 fails; check IR distributions.

6) Validation & sanity checks you should run
   - View NAs and counts:
       zcat ${RUN_DIR}/per_family/diffIR_*.tsv.gz | head
     Inspect the header and a few lines.
   - For top hits, check read-level support using IRFinder original metrics or by examining BAMs.
   - Plot simple scatter plots of IR_throm vs IR_nothrom for top introns to ensure sensible effects.

7) Improvements & extensions
   - Replace leave-one-out sigma with pooled empirical Bayes variance estimates across introns,
     or use a variance-stabilizing transform tailored to IR data.
   - Add covariates (batch, sample-level QC metrics) if families have technical confounders.
   - Expand per-family plots: label top significant introns, use density plots of delta distribution.
   - Add FDR control across meta-analysis stage (we compute q_meta here via BH).
   - Add sample/contig whitelist or filter introns by minimum coverage / junction support before testing.

8) Reproducibility
   - Note the R version and package versions in logs (Rscript --version is saved).
   - Save checksums of the IRFinder outputs used for the run if you want cryptographic provenance.
   - Consider wrapping this workflow in a SnakeMake or Nextflow pipeline for robust re-running.

9) Help & debugging
   - If Rscript fails: inspect ${RUN_DIR}/logs/Rscript.stdout.log for the R-side error trace.
   - If you need me to adapt the R code (for example to accept a custom IR column name,
     or to change the variance estimator), tell me which change you want and I can produce a patch.

END_NOTES
