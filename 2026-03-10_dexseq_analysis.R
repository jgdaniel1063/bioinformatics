#!/usr/bin/env Rscript
#
# 2026-03-04_dexseq-analysis_annotated.R
#
# DEXSeq pipeline with inline sample definitions.
# This heavily annotated version explains rationale, edge-cases, troubleshooting,
# recommended extensions, and interpretation notes for every major step.
#
# Usage:
#   Rscript 2026-03-04_dexseq-analysis_annotated.R 2>&1 | tee ~/dexseq_run_annotated.log


# ======================== CONFIGURATION ========================
# Paths and thresholds are centralized here for easy editing.
# Verify these paths are correct before running.
FLATTENED_GTF  <- "/home/jgd/Documents/bioinformatics_working/ref/danio_rerio.grcz11.115.dexseq.gtf"
OUTPUT_DIR     <- "/home/jgd/Documents/bioinformatics_working/output/dexseq_results"
RESULTS_RDS    <- file.path(OUTPUT_DIR, "dexseq_results.rds")
RESULTS_TSV    <- file.path(OUTPUT_DIR, "dexseq_results.tsv")

# Statistical threshold for reporting significance (adjust if you prefer stricter or more lenient)
FDR_THRESHOLD  <- 0.05

# Directory containing per-sample DEXSeq count files (created by dexseq_count.py)
COUNTS_DIR     <- "/home/jgd/Documents/bioinformatics_working/2025-11-06_enu-proc/counts/dexseq/dexseq_counts"

# -----------------------
# Notes on reproducibility:
# - Record package versions (sessionInfo()) and include in the run log.
# - Consider pinning Bioconductor/DEXSeq versions in a manifest or Dockerfile.
# -----------------------

# ======================== SAMPLE LIST ==========================
# Each row: c(sample_name, condition, filename_in_COUNTS_DIR)
# To include/exclude samples simply comment/uncomment relevant rows.
#
# Important: sample_name must be unique; filenames must exist in COUNTS_DIR.
# The 'condition' is the experimental factor used for the exon-usage interaction
# (e.g., "thrombosis" vs "no_thrombosis"). You can expand to more complex
# designs (see notes later).
SAMPLES <- list(
  # --- ai08 family ---
 #  c("ai08nothrom1", "no_thrombosis",  "ai08nothrom1_dexseq_counts.txt"),
  # c("ai08nothrom2", "no_thrombosis",  "ai08nothrom2_dexseq_counts.txt"),
  # c("ai08nothrom3", "no_thrombosis",  "ai08nothrom3_dexseq_counts.txt"),
  # c("ai08throm1",   "thrombosis",     "ai08throm1_dexseq_counts.txt"),
  # c("ai08throm2",   "thrombosis",     "ai08throm2_dexseq_counts.txt"),
  # c("ai08throm3",   "thrombosis",     "ai08throm3_dexseq_counts.txt")
  # --- au01 family ---
  # c("au01nothrom1", "no_thrombosis",  "au01nothrom1_dexseq_counts.txt"),
  # c("au01throm1",   "thrombosis",     "au01throm1_dexseq_counts.txt"),
  # --- bo01 family ---
  # c("bo01nothrom1", "no_thrombosis",  "bo01nothrom1_dexseq_counts.txt"),
  # c("bo01throm1",   "thrombosis",     "bo01throm1_dexseq_counts.txt"),
  # --- bq01 family ---
  # c("bq01nothrom1", "no_thrombosis",  "bq01nothrom1_dexseq_counts.txt"),
  # c("bq01nothrom2", "no_thrombosis",  "bq01nothrom2_dexseq_counts.txt"),
  # c("bq01nothrom3", "no_thrombosis",  "bq01nothrom3_dexseq_counts.txt"),
  # c("bq01throm1",   "thrombosis",     "bq01throm1_dexseq_counts.txt"),
  # c("bq01throm2",   "thrombosis",     "bq01throm2_dexseq_counts.txt"),
  # c("bq01throm3",   "thrombosis",     "bq01throm3_dexseq_counts.txt")
  # --- bq06 family ---
  c("bq06nothrom1", "no_thrombosis",  "bq06nothrom1_dexseq_counts.txt"),
  c("bq06nothrom2", "no_thrombosis",  "bq06nothrom2_dexseq_counts.txt"),
  c("bq06nothrom3", "no_thrombosis",  "bq06nothrom3_dexseq_counts.txt"),
  c("bq06throm1",   "thrombosis",     "bq06throm1_dexseq_counts.txt"),
  c("bq06throm2",   "thrombosis",     "bq06throm2_dexseq_counts.txt"),
  c("bq06throm3",   "thrombosis",     "bq06throm3_dexseq_counts.txt")
  # --- bqm2 family ---
  # c("bqm2nothrom1", "no_thrombosis",  "bqm2nothrom1_dexseq_counts.txt"),
  # c("bqm2throm1",   "thrombosis",     "bqm2throm1_dexseq_counts.txt")
  # --- positive controls ---
  # c("pc04",         "thrombosis",     "pc04_dexseq_counts.txt"),
  # c("pc05",         "thrombosis",     "pc05_dexseq_counts.txt"),
  # c("pc06",         "thrombosis",     "pc06_dexseq_counts.txt"),
  # --- wild-type controls ---
  # c("wt01",         "no_thrombosis",  "wt01_dexseq_counts.txt"),
  # c("wt02",         "no_thrombosis",  "wt02_dexseq_counts.txt"),
  # c("wt03",         "no_thrombosis",  "wt03_dexseq_counts.txt")
)
# ===============================================================

# Suppress package startup messages to keep the log readable
suppressPackageStartupMessages({
  library(DEXSeq)     # core package for exon-level differential usage
  library(ggplot2)    # plotting
  library(pheatmap)   # heatmaps of exon counts
  library(dplyr)      # data manipulation (used for volcano)
})

# Create output directory (recursive) if necessary
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ----- helpers -----
# Small logging wrapper to prefix timestamps for each message — useful in long runs.
logmsg <- function(...) {
  cat(sprintf("[%s] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), ..., "\n")
}

# sanitize_count_file
# Purpose:
#  - Some count files produced by dexseq_count.py include trailing summary rows
#    starting with '_' (underscore). DEXSeq's import functions do not expect
#    those and may fail. This helper reads chunked blocks and writes only the
#    non-underscore lines to a temporary sanitized file used for import.
#
# Implementation notes:
#  - Uses readLines in chunks (n = 10000) to avoid reading very large files into memory.
#  - The function is conservative and will preserve any line that does not start with
#    optional whitespace followed by underscore: '^\\s*_'
sanitize_count_file <- function(infile, outfile) {
  con_in  <- file(infile, "r")
  con_out <- file(outfile, "w")
  on.exit({ try(close(con_in), silent = TRUE); try(close(con_out), silent = TRUE) })
  repeat {
    lines <- readLines(con_in, n = 10000, warn = FALSE)
    if (length(lines) == 0) break
    keep <- lines[!grepl("^\\s*_", lines)]
    if (length(keep) > 0) writeLines(keep, con_out)
  }
  invisible(TRUE)
}

# -----------------------
# Design notes:
# DEXSeq models the count y_ijk for exon i, sample j, gene k using generalized
# linear models with a sample effect, exon effect and an interaction term
# condition:exon that tests whether exon usage changes across condition.
#
# In this script we set the full model:
#   ~ sample + exon + condition:exon
# and reduced:
#   ~ sample + exon
# and testForDEU(dxd, reducedModel = ~ sample + exon) performs the likelihood
# ratio test for the interaction term.
# -----------------------

# ----- main -----
tryCatch({

  # ---------- 1. Parse sample list ----------
  # Basic validations and building of sampleData used by DEXSeq
  if (length(SAMPLES) == 0) stop("SAMPLES list is empty. Uncomment at least two samples.")

  # Extract columns from the SAMPLES list-of-character vectors
  sample_names <- vapply(SAMPLES, `[`, character(1), 1)
  conditions   <- vapply(SAMPLES, `[`, character(1), 2)
  filenames    <- vapply(SAMPLES, `[`, character(1), 3)
  count_files  <- file.path(COUNTS_DIR, filenames)

  # Check for duplicate sample names (DEXSeq requires unique row.names)
  if (any(duplicated(sample_names))) {
    stop("Duplicate sample names: ", paste(sample_names[duplicated(sample_names)], collapse = ", "))
  }

  # Check the count files exist before doing heavy work — fail early with a clear message
  missing_files <- count_files[!file.exists(count_files)]
  if (length(missing_files) > 0) {
    stop("Count files not found:\n", paste("  ", missing_files, collapse = "\n"))
  }

  # sampleData: DEXSeq expects rownames = sample names and a dataframe of covariates
  sampleData <- data.frame(
    row.names = sample_names,
    condition = factor(conditions)
  )

  logmsg("Samples included:", length(sample_names))
  logmsg(" ", paste(sample_names, collapse = ", "))
  logmsg("Condition summary:")
  print(table(sampleData$condition))

  # Save the sample table used for the run — good provenance practice
  sample_table_copy <- file.path(OUTPUT_DIR, "sample_table_used.tsv")
  write.table(
    data.frame(sample = sample_names, condition = conditions, count_file = count_files),
    file = sample_table_copy, sep = "\t", quote = FALSE, row.names = FALSE
  )
  logmsg("Wrote sample table:", sample_table_copy)

  # ---------- 2. Sanitize count files ----------
  logmsg("Sanitizing count files (removing '_' summary rows)...")
  san_tmp <- tempfile("dexseq_count_sanitized_")
  dir.create(san_tmp)
  san_files <- character(length(count_files))
  for (i in seq_along(count_files)) {
    outfile <- file.path(san_tmp, basename(count_files[i]))
    sanitize_count_file(count_files[i], outfile)
    san_files[i] <- outfile
  }
  logmsg("Preview first sanitized file:")
  cat(paste0(readLines(san_files[1], n = 6), collapse = "\n"), "\n")

  # ---------- 3. Build DEXSeqDataSet ----------
  # NOTE: Different DEXSeq versions or R versions may have slightly different
  # function signatures for DEXSeqDataSetFromHTSeq. We call it inside try() and
  # fallback to a named-argument invocation if necessary. If all fail, we print
  # the function formals to aid debugging.
  logmsg("Building DEXSeq dataset...")
  dxd <- NULL
  try({
    # Common simple call: positional args (older examples use this style)
    dxd <- DEXSeqDataSetFromHTSeq(
      san_files, sampleData,
      design = ~ sample + exon + condition:exon,
      flattenedfile = FLATTENED_GTF
    )
  }, silent = TRUE)

  if (is.null(dxd)) {
    try({
      # More explicit named-argument style (some versions require 'countFiles' and 'sampleData' names)
      dxd <- DEXSeqDataSetFromHTSeq(
        countFiles    = san_files,
        sampleData    = sampleData,
        design        = ~ sample + exon + condition:exon,
        flattenedfile = FLATTENED_GTF
      )
    }, silent = TRUE)
  }

  if (is.null(dxd)) {
    logmsg("Failed to construct DEXSeq dataset. Signature for debugging:")
    print(formals(DEXSeqDataSetFromHTSeq))
    stop("Could not construct DEXSeqDataSetFromHTSeq; adjust call for your DEXSeq version.")
  }

  # ---------- 4. DEXSeq workflow ----------
  # Key steps:
  #  - estimateSizeFactors: normalizes for library size differences (by size factors)
  #  - estimateDispersions: models dispersion across exons (critical step)
  #  - testForDEU: performs LRT for exon-by-condition interaction
  #  - estimateExonFoldChanges: optional, computes fold-changes per exon relative to condition
  #
  logmsg("Estimating size factors...")
  dxd <- estimateSizeFactors(dxd)

  logmsg("Estimating dispersions...")
  # Some DEXSeq versions accept a formula argument to control dispersion estimation;
  # try both styles and fall back if needed.
  try({ dxd <- estimateDispersions(dxd, formula = ~ sample + exon + condition:exon) }, silent = TRUE)
  if (!"dispersionFunction" %in% names(mcols(dxd))) {
    # If the previous call didn't set dispersions as expected, attempt the simpler call
    try({ dxd <- estimateDispersions(dxd) }, silent = TRUE)
  }

  logmsg("Testing for differential exon usage (DEU)...")
  # Reduced model excludes the condition:exon interaction so the LRT tests that interaction.
  dxd <- testForDEU(dxd, reducedModel = ~ sample + exon)

  logmsg("Estimating exon fold changes...")
  # fitExpToVar = "condition" aligns fold-change estimates to the condition factor
  dxd <- estimateExonFoldChanges(dxd, fitExpToVar = "condition")

  # ---------- 5. Extract & save results ----------
  # DEXSeqResults contains per-exon test statistics including pvalue and padj
  res    <- DEXSeqResults(dxd)
  res_df <- as.data.frame(res)
  logmsg("Tested exon bins:", nrow(res_df))

  # Save the results object as RDS for reproducibility and future re-loading
  saveRDS(res, file = RESULTS_RDS)
  logmsg("Saved RDS ->", RESULTS_RDS)

  # Some DEXSeq result columns are list-columns (complex objects). For TSV export,
  # we flatten any list-columns by concatenating elements with ';' to avoid writing R lists.
  list_cols <- names(res_df)[sapply(res_df, is.list)]
  if (length(list_cols) > 0) {
    logmsg("Flattening list-columns:", paste(list_cols, collapse = ", "))
    for (col in list_cols) {
      res_df[[col]] <- vapply(res_df[[col]], function(x) {
        if (is.null(x)) return(NA_character_)
        paste(as.character(unlist(x)), collapse = ";")
      }, character(1))
    }
  }

  # Ensure an exon identifier column exists for TSV consumers (some DEXSeq versions include rownames)
  if (!"exon_id" %in% names(res_df)) {
    res_df <- cbind(exon_id = rownames(res_df), res_df)
  }

  # Write a tab-separated, human-readable TSV of results
  write.table(res_df, file = RESULTS_TSV, sep = "\t", quote = FALSE, row.names = FALSE)
  logmsg("Wrote TSV ->", RESULTS_TSV)

  # ---------- 6. Plots & diagnostics ----------
  # Identify significant exon bins based on adjusted p-value threshold (FDR)
  sig_exons <- res_df[!is.na(res_df$padj) & res_df$padj < FDR_THRESHOLD, ]
  logmsg("Significant exon bins (padj <", FDR_THRESHOLD, "):", nrow(sig_exons))

  # MA plot: useful global visual of log fold-change vs mean expression
  try({
    png(file.path(OUTPUT_DIR, "ma_plot.png"), width = 1600, height = 1200, res = 200)
    # plotMA expects a DEXSeqResults object; `res` is that object stored earlier.
    plotMA(res, ylim = c(-5, 5), main = "MA Plot")
    dev.off()
    logmsg("MA plot written.")
  }, silent = TRUE)

  # Volcano plot: custom ggplot visualising log2 fold and -log10(padj)
  # Note: DEXSeq's exact fold-change column name can vary by version. We attempt to find a matching column.
  fc_col <- intersect(c("log2fold", "log2fold_treatment_control"), colnames(res_df))[1]
  if (!is.na(fc_col)) {
    try({
      png(file.path(OUTPUT_DIR, "volcano_plot.png"), width = 1600, height = 1200, res = 200)
      volcano_data <- res_df %>% mutate(significant = !is.na(padj) & padj < FDR_THRESHOLD)
      print(ggplot(volcano_data, aes(x = .data[[fc_col]], y = -log10(padj))) +
              geom_point(aes(color = significant), alpha = 0.6) +
              scale_color_manual(values = c("grey", "red")) +
              labs(title = "Volcano Plot", x = "Log2 Fold Change", y = "-Log10 Adjusted P-value") +
              theme_minimal())
      dev.off()
      logmsg("Volcano plot written.")
    }, silent = TRUE)
  } else {
    logmsg("Warning: could not find a log2 fold-change column for volcano plot. Columns present:", paste(colnames(res_df), collapse = ", "))
  }

  # Per-gene exon usage plots (top 5 significant genes) — plotDEXSeq uses the DEXSeqDataSet
  if (nrow(sig_exons) > 0) {
    top_genes <- head(unique(sig_exons$groupID), 5)
    for (g in top_genes) {
      try({
        png(file.path(OUTPUT_DIR, paste0("exon_usage_", g, ".png")),
            width = 1600, height = 1200, res = 200)
        # plotDEXSeq visualises counts per exon across conditions; it requires the DEXSeq dataset
        plotDEXSeq(dxd, g, legend = TRUE, main = paste("Exon Usage:", g))
        dev.off()
      }, silent = TRUE)
    }
    logmsg("Exon usage plots generated.")
  } else {
    logmsg("No significant genes found — skipping exon usage plots.")
  }

  # Heatmap of significant exon bins — helpful to visually inspect patterns across samples
  try({
    if (nrow(sig_exons) > 0) {
      # counts(dxd) returns the raw count matrix (rows = exon bins, cols = samples)
      mat <- counts(dxd)[rownames(sig_exons), , drop = FALSE]
      annot_col <- data.frame(condition = sampleData$condition)
      rownames(annot_col) <- rownames(sampleData)
      pheatmap(mat, annotation_col = annot_col,
               filename = file.path(OUTPUT_DIR, "exon_heatmap.png"),
               show_rownames = FALSE)
      logmsg("Exon heatmap written.")
    }
  }, silent = TRUE)

  # Dispersion plot — inspect model fit and dispersion estimates
  try({
    png(file.path(OUTPUT_DIR, "dispersion_plot.png"), width = 1600, height = 1200, res = 200)
    plotDispEsts(dxd)
    dev.off()
    logmsg("Dispersion plot written.")
  }, silent = TRUE)

  # P-value histogram — helps detect p-value distribution anomalies (e.g., conservative behaviour)
  try({
    png(file.path(OUTPUT_DIR, "pvalue_histogram.png"), width = 1600, height = 1200, res = 200)
    hist(res_df$pvalue, breaks = 50, main = "P-value Distribution", xlab = "P-value")
    dev.off()
    logmsg("P-value histogram written.")
  }, silent = TRUE)

  # Final log message with location of results
  logmsg("DEXSeq analysis complete. Results in:", OUTPUT_DIR)

  # Save session information for reproducibility — package versions, R version, etc.
  try({
    session_info_file <- file.path(OUTPUT_DIR, "sessionInfo.txt")
    writeLines(capture.output(sessionInfo()), con = session_info_file)
    logmsg("Saved sessionInfo ->", session_info_file)
  }, silent = TRUE)

}, error = function(e) {
  # Centralized error handler: prints traceback and the error message and exits with non-zero status.
  cat("ERROR during DEXSeq run:\n")
  traceback()
  cat("Error message:", conditionMessage(e), "\n")
  quit(status = 1)
})

# ========================= EXPLANATORY NOTES & RECOMMENDATIONS =========================
#
# 1) Input count file format:
#    - Produced by: dexseq_count.py script from DEXSeq pipeline (Python).
#    - Expected to be a 2-column file: exon_id <tab> counts
#    - Some files include a trailing summary like "_total" lines; sanitize_count_file removes those.
#
# 2) Experimental design & formula:
#    - Current (simple) design: ~ sample + exon + condition:exon (reduced: ~ sample + exon)
#    - This effectively tests whether the proportion of counts assigned to an exon
#      (relative to its gene) changes with 'condition'.
#    - For more complex experiments, include additional covariates in sampleData:
#        sampleData$batch <- factor(batch_vector)
#      and modify the design accordingly. Example:
#        design = ~ sample + exon + batch:exon + condition:exon
#      but be careful: adding many interaction terms consumes degrees of freedom.
#
# 3) Replicates, balance, and power:
#    - DEXSeq requires replicate samples per condition to estimate dispersion reliably.
#    - Unbalanced designs (very few replicates in one condition) reduce power and can lead
#      to unstable dispersion estimates. Aim for >=3 biological replicates per condition.
#
# 4) Filtering:
#    - This script does not perform pre-filtering on low-count exons. Low-count bins can
#      inflate dispersion estimates and produce noisy tests. Consider pre-filtering:
#         keep <- rowSums(counts(dxd) >= 10) >= min_samples
#         dxd <- dxd[keep, ]
#      Adjust thresholds to your experiment and desired sensitivity.
#
# 5) Multiple testing:
#    - Results include pvalue and padj (Benjamini-Hochberg adjusted p-values).
#    - Use padj < FDR_THRESHOLD to select significant exon bins. You may tune FDR_THRESHOLD.
#
# 6) Interpretation:
#    - A significant exon bin means the relative usage of that exon (vs other exons of the gene)
#      differs between conditions. This is not necessarily differential gene expression.
#    - Follow-up validation: inspect per-exon plots (plotDEXSeq) and consider RT-PCR for confirmation.
#
# 7) Performance & parallelization:
#    - estimateDispersions can be slow for large datasets. DEXSeq uses internal optimization; to gain speed:
#       * Use BiocParallel to parallelize (register MulticoreParam / SnowParam) if supported.
#       * Run on a machine with sufficient memory; counts(dxd) can be large if many exons.
#
# 8) Versioning & reproducibility:
#    - Save sessionInfo() (done above). Consider full containerization (Docker/Singularity)
#      with Bioconductor/DEXSeq pinned to specific versions.
#
# 9) Robustness:
#    - If DEXSeqDataSetFromHTSeq fails, inspect its formals (printed by this script) and adapt call.
#    - If count file names or columns don't match expectations, open a sanitized file and inspect the header.
#
# 10) Extensions you might want:
#    - Add an argparse-like interface (e.g., optparse or argparse) so you can pass paths and thresholds
#      from the command line instead of editing the script.
#    - Add pre-filtering of low-count exons (recommended) with a tunable threshold.
#    - Export per-gene aggregate results (combine exon-level p-values with Simes' method or use transcript models).
#    - Integrate with tximport/DEXSeq pipeline if you prefer transcript-level inputs.
#    - Produce HTML report (e.g., rmarkdown) that includes tables, plots, and sessionInfo.
#
# If you'd like, I can:
#  - Add a CLI wrapper using optparse so the script accepts arguments (paths, thresholds).
#  - Implement pre-filtering of low-count exons with recommended defaults.
#  - Add BiocParallel configuration to speed dispersion estimation.
#  - Generate a minimal unit-test or example fixture (small synthetic count files) to validate the pipeline.
#
# Tell me which extension you prefer and I'll implement it.
#
# End of annotated script.
