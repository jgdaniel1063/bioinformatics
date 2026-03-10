#!/usr/bin/env Rscript
#
# compare_two_groups_hardcoded.R
#
# Hard-coded script to compare two groups of count files and rank features by log2FC
# and absolute expression difference. Supports:
#  - "statistical" testing with edgeR (estimate dispersion from data when possible)
#  - "manual dispersion" (supply BCV) when no replicates exist
#  - "no dispersion" mode: compute normalized fold-changes and absolute differences with no p-values
#
# Edit the HARD-CODED section below (paths, parameters) and run:
#   Rscript compare_two_groups_hardcoded.R
#
# Outputs (in OUTPUT_DIR):
#  - counts_matrix.csv           : merged raw counts matrix (genes x samples)
#  - logCPM_matrix.csv           : normalized logCPM matrix
#  - results_full.csv            : full ranked table (log2FC, meanA, meanB, absDiff, pval, FDR if available)
#  - results_ranked_log2fc.csv   : ranked by log2FC
#  - results_ranked_absdiff.csv  : ranked by absolute mean difference
#  - topN_by_log2fc.csv, topN_by_absdiff.csv
#  - PCA, MA, Volcano plots (PDF/PNG when applicable)
#  - session_info.txt
#
suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(pheatmap)
})

# =============================================================================
# HARD-CODED CONFIGURATION - EDIT THIS SECTION
# =============================================================================

# Input: either supply two directories with per-sample count files (HTSeq-style:
# two columns: gene \t count) OR provide a counts matrix file (genes x samples).
# When both are provided, counts matrix is used.
COUNTS_MATRIX <- ""                     # e.g. "/home/user/counts_matrix.tsv" (leave empty to use group dirs)
GROUPA_DIR    <- "/path/to/groupA_counts" # e.g. "/data/groupA_counts"
GROUPB_DIR    <- "/path/to/groupB_counts" # e.g. "/data/groupB_counts"
# Optionally only include files matching these substrings in the directories
GROUPA_PREFIX <- ""   # e.g. "ai08" (leave empty to take all files)
GROUPB_PREFIX <- ""   # e.g. "ai08"

# Output
OUTPUT_DIR <- "./compare_two_groups_output"

# Filtering & ranking parameters
MIN_COUNT <- 10        # minimum raw counts to consider a gene "expressed"
MIN_SAMPLES <- 2       # minimum number of samples with >= MIN_COUNT for gene to be kept
MIN_LOG2FC <- 0        # minimum absolute log2 fold-change for filtered outputs
PADJ_CUTOFF <- 0.05
TOP_N <- 200

# Statistical mode:
# Choose one of: "edgeR_estimate" (estimate from data if replicates exist),
#               "edgeR_manual"   (use manual BCV - see MANUAL_BCV),
#               "no_dispersion"  (no p-values; compute normalized fold changes only)
MODE <- "edgeR_estimate"

# If using manual mode, specify BCV (biological coefficient of variation). dispersion = BCV^2
MANUAL_BCV <- 0.2   # example: 0.2 -> dispersion = 0.04

# Plotting & misc
VST_PSEUDO <- 1     # pseudo-count for log2 transforms
SEED <- 20240301

# =============================================================================
# END HARD-CODED CONFIG
# =============================================================================

set.seed(SEED)
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

log_msg <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste0(...)))

# Helper: read per-sample 2-column count files (HTSeq-like). Returns named list of named integer vectors.
read_counts_dir <- function(dirpath, prefix = "") {
  if (!dir.exists(dirpath)) return(list())
  files <- list.files(dirpath, full.names = TRUE)
  if (nzchar(prefix)) files <- files[grepl(prefix, basename(files))]
  files <- files[file.info(files)$isdir == FALSE]
  out <- list()
  for (f in files) {
    # support gz (read.table handles gz by filename on many systems)
    df <- tryCatch(read.table(f, header = FALSE, sep = "\t", stringsAsFactors = FALSE, quote = "", comment.char = ""), error = function(e) NULL)
    if (is.null(df) || ncol(df) < 2) {
      warning("Skipping non-2-col file: ", f); next
    }
    df <- df[!grepl("^__", df[[1]]), , drop = FALSE]   # remove HTSeq summary rows
    sample <- tools::file_path_sans_ext(basename(f))
    out[[sample]] <- as.integer(round(as.numeric(df[[2]])))
    names(out[[sample]]) <- as.character(df[[1]])
  }
  return(out)
}

# Build counts matrix
if (nzchar(COUNTS_MATRIX)) {
  if (!file.exists(COUNTS_MATRIX)) stop("Counts matrix not found: ", COUNTS_MATRIX)
  log_msg("Reading counts matrix:", COUNTS_MATRIX)
  counts_df <- as.data.frame(read.table(COUNTS_MATRIX, header = TRUE, row.names = 1, check.names = FALSE, sep = "\t", stringsAsFactors = FALSE))
  counts_df[] <- lapply(counts_df, function(x) as.integer(round(as.numeric(x))))
} else {
  log_msg("Reading per-sample counts from directories:")
  log_msg("  Group A:", GROUPA_DIR, " prefix:", GROUPA_PREFIX)
  log_msg("  Group B:", GROUPB_DIR, " prefix:", GROUPB_PREFIX)
  a_list <- read_counts_dir(GROUPA_DIR, GROUPA_PREFIX)
  b_list <- read_counts_dir(GROUPB_DIR, GROUPB_PREFIX)
  if (length(a_list) + length(b_list) == 0) stop("No count files found in both group directories.")
  # union of genes
  all_genes <- sort(unique(unlist(lapply(c(a_list, b_list), names))))
  all_samples <- c(names(a_list), names(b_list))
  counts_mat <- matrix(0L, nrow = length(all_genes), ncol = length(all_samples), dimnames = list(all_genes, all_samples))
  for (s in names(a_list)) counts_mat[names(a_list[[s]]), s] <- a_list[[s]]
  for (s in names(b_list)) counts_mat[names(b_list[[s]]), s] <- b_list[[s]]
  counts_df <- as.data.frame(counts_mat)
}

# Create sample grouping vector
groupA_samples <- intersect(colnames(counts_df), names(read_counts_dir(GROUPA_DIR, GROUPA_PREFIX)))
groupB_samples <- intersect(colnames(counts_df), names(read_counts_dir(GROUPB_DIR, GROUPB_PREFIX)))
# If counts matrix was supplied, user must rely on directory filenames to derive groups; if none found, assume first half/second half
if (nzchar(COUNTS_MATRIX) && (length(groupA_samples) == 0 || length(groupB_samples) == 0)) {
  # fallback: split columns in half (not ideal; user should edit configuration)
  ncol_tot <- ncol(counts_df)
  groupA_samples <- colnames(counts_df)[1:floor(ncol_tot/2)]
  groupB_samples <- colnames(counts_df)[(floor(ncol_tot/2)+1):ncol_tot]
  log_msg("Warning: groups not auto-detected from directories; falling back to split-by-half. Edit HARD-CODED section if needed.")
}

if (length(groupA_samples) < 1 || length(groupB_samples) < 1) {
  stop("Unable to determine group samples. Ensure directories and prefixes or counts matrix are set correctly.")
}

log_msg("Group A samples (", length(groupA_samples), "): ", paste(groupA_samples, collapse = ", "))
log_msg("Group B samples (", length(groupB_samples), "): ", paste(groupB_samples, collapse = ", "))

# Reorder columns to groupA then groupB
sample_order <- c(groupA_samples, groupB_samples)
counts_df <- counts_df[, sample_order, drop = FALSE]

# Save raw counts
write.csv(as.data.frame(counts_df), file = file.path(OUTPUT_DIR, "counts_matrix.csv"), row.names = TRUE)

# Filtering
keep_by_count <- rowSums(counts_df >= MIN_COUNT) >= MIN_SAMPLES
log_msg("Filtering genes: keeping ", sum(keep_by_count), " / ", nrow(counts_df), " with >= ", MIN_COUNT, " counts in >=", MIN_SAMPLES, " samples")
counts_df_filt <- counts_df[keep_by_count, , drop = FALSE]

# Build DGEList and normalize
dge <- DGEList(counts = as.matrix(counts_df_filt))
dge <- calcNormFactors(dge)

# compute logCPM for mean expression metrics
logCPM_mat <- cpm(dge, log=TRUE, prior.count = 1)
meanA_logCPM <- rowMeans(logCPM_mat[, groupA_samples, drop = FALSE])
meanB_logCPM <- rowMeans(logCPM_mat[, groupB_samples, drop = FALSE])
absDiff_logCPM <- meanA_logCPM - meanB_logCPM
absDiffAbs_logCPM <- abs(absDiff_logCPM)

# Decide statistical approach
results_tbl <- NULL
if (MODE == "no_dispersion") {
  # No p-values; just compute normalized fold changes and absolute differences
  log2FC <- meanA_logCPM - meanB_logCPM  # on logCPM scale approximates log2FC
  results_tbl <- tibble(
    feature = rownames(logCPM_mat),
    log2FC = log2FC,
    meanA_logCPM = meanA_logCPM,
    meanB_logCPM = meanB_logCPM,
    absDiff = absDiff_logCPM,
    absDiffAbs = absDiffAbs_logCPM,
    pvalue = NA_real_,
    FDR = NA_real_
  )
  log_msg("Mode=no_dispersion: no statistical testing performed; only fold-changes/means computed.")
} else {
  # Use edgeR GLM framework. Determine if we can estimate dispersion from data.
  group_factor <- factor(c(rep("A", length(groupA_samples)), rep("B", length(groupB_samples))))
  design <- model.matrix(~ group_factor)
  colnames(design) <- c("Intercept", "GroupB_vs_A")
  dge <- DGEList(counts = as.matrix(counts_df_filt), group = group_factor)
  dge <- calcNormFactors(dge)

  can_estimate <- (any(table(group_factor) >= 2) && MODE == "edgeR_estimate")
  if (can_estimate) {
    log_msg("Estimating dispersions from data (replicates detected).")
    dge <- estimateDisp(dge, design, robust = TRUE)
    fit <- glmFit(dge, design)
    lrt <- glmLRT(fit, coef = "GroupB_vs_A")
    tt <- topTags(lrt, n = nrow(dge$counts), sort.by = "none")$table
    # note: edgeR reports logFC = log2(B/A)
    results_tbl <- as_tibble(tt, rownames = "feature") %>%
      mutate(
        meanA_logCPM = meanA_logCPM[feature],
        meanB_logCPM = meanB_logCPM[feature],
        absDiff = (meanA_logCPM - meanB_logCPM),
        absDiffAbs = abs(absDiff)
      ) %>% rename(log2FC = logFC, pvalue = PValue, FDR = FDR)
    dispersion_info <- list(type = "estimated", common_dispersion = dge$common.dispersion)
    log_msg("Dispersion estimated. Common dispersion: ", signif(dge$common.dispersion, 4))
  } else if (!is.na(MANUAL_BCV) && MANUAL_BCV > 0 && MODE == "edgeR_manual") {
    # Use manual BCV -> dispersion = BCV^2
    disp_val <- MANUAL_BCV^2
    log_msg("No replicates or manual mode chosen. Using manual BCV=", MANUAL_BCV, " -> dispersion=", signif(disp_val,4))
    fit <- glmFit(dge, design, dispersion = disp_val)
    lrt <- glmLRT(fit, coef = "GroupB_vs_A")
    tt <- topTags(lrt, n = nrow(dge$counts), sort.by = "none")$table
    results_tbl <- as_tibble(tt, rownames = "feature") %>%
      mutate(
        meanA_logCPM = meanA_logCPM[feature],
        meanB_logCPM = meanB_logCPM[feature],
        absDiff = (meanA_logCPM - meanB_logCPM),
        absDiffAbs = abs(absDiff)
      ) %>% rename(log2FC = logFC, pvalue = PValue, FDR = FDR)
    dispersion_info <- list(type = "manual", BCV = MANUAL_BCV, dispersion = disp_val)
  } else {
    # Can't estimate and manual BCV not provided or mode mismatch -> fallback: compute fold-changes only
    log_msg("Cannot estimate dispersion (no replicates) and manual BCV not provided/or not in manual mode.")
    log_msg("Falling back to no_dispersion behavior: will compute fold-changes only.")
    log2FC <- meanA_logCPM - meanB_logCPM
    results_tbl <- tibble(
      feature = rownames(logCPM_mat),
      log2FC = log2FC,
      meanA_logCPM = meanA_logCPM,
      meanB_logCPM = meanB_logCPM,
      absDiff = absDiff_logCPM,
      absDiffAbs = absDiffAbs_logCPM,
      pvalue = NA_real_,
      FDR = NA_real_
    )
    dispersion_info <- list(type = "none")
  }
}

# Common post-processing: rank and filter
results_tbl <- results_tbl %>% mutate(rank_by_log2fc = rank(-abs(log2FC), ties.method = "first"),
                                      rank_by_absdiff = rank(-abs(absDiff), ties.method = "first"))

# Save logCPM matrix
write.csv(as.data.frame(logCPM_mat), file = file.path(OUTPUT_DIR, "logCPM_matrix.csv"), row.names = TRUE)

# Save full results
write.csv(results_tbl, file = file.path(OUTPUT_DIR, "results_full.csv"), row.names = FALSE)

# Ranked outputs
rank_log2fc <- results_tbl %>% arrange(desc(abs(log2FC)))
rank_absdiff <- results_tbl %>% arrange(desc(abs(absDiff)))

write.csv(rank_log2fc, file = file.path(OUTPUT_DIR, "results_ranked_log2fc.csv"), row.names = FALSE)
write.csv(rank_absdiff, file = file.path(OUTPUT_DIR, "results_ranked_absdiff.csv"), row.names = FALSE)

# Apply user filters (padj, min log2fc) for filtered lists
filtered <- results_tbl
if (!all(is.na(results_tbl$FDR))) {
  filtered <- filtered %>% filter(is.na(FDR) | FDR <= PADJ_CUTOFF)
}
if (!is.null(MIN_LOG2FC) && MIN_LOG2FC > 0) {
  filtered <- filtered %>% filter(abs(log2FC) >= MIN_LOG2FC)
}
write.csv(filtered, file = file.path(OUTPUT_DIR, "results_filtered.csv"), row.names = FALSE)

# Top N lists
write.csv(head(rank_log2fc, TOP_N), file = file.path(OUTPUT_DIR, paste0("top", TOP_N, "_by_log2fc.csv")), row.names = FALSE)
write.csv(head(rank_absdiff, TOP_N), file = file.path(OUTPUT_DIR, paste0("top", TOP_N, "_by_absdiff.csv")), row.names = FALSE)

# Plots
# PCA on logCPM
pca <- prcomp(t(logCPM_mat))
pc_df <- data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2], sample = colnames(logCPM_mat),
                    group = factor(c(rep("A", length(groupA_samples)), rep("B", length(groupB_samples)))))
p_pca <- ggplot(pc_df, aes(PC1, PC2, color = group, label = sample)) +
  geom_point(size = 3) + geom_text(hjust = 0.5, vjust = -1, size = 2.5) +
  labs(title = "PCA (logCPM)") + theme_bw()
ggsave(filename = file.path(OUTPUT_DIR, "PCA_logCPM.pdf"), plot = p_pca, width = 7, height = 6)
ggsave(filename = file.path(OUTPUT_DIR, "PCA_logCPM.png"), plot = p_pca, width = 7, height = 6, dpi = 300)

# MA plot (if stats available)
if (!all(is.na(results_tbl$pvalue))) {
  ma_df <- results_tbl %>% mutate(meanExpr = (meanA_logCPM + meanB_logCPM)/2)
  p_ma <- ggplot(ma_df, aes(meanExpr, log2FC, color = (FDR <= PADJ_CUTOFF))) +
    geom_point(alpha = 0.5) +
    labs(title = "MA plot", x = "Mean logCPM", y = "log2FC (A - B)") +
    theme_bw()
  ggsave(file.path(OUTPUT_DIR, "MA_plot.pdf"), p_ma, width = 7, height = 5)
  ggsave(file.path(OUTPUT_DIR, "MA_plot.png"), p_ma, width = 7, height = 5, dpi = 300)
}

# Volcano (if p-values available)
if (!all(is.na(results_tbl$pvalue))) {
  volcano_df <- results_tbl %>% mutate(negLog10FDR = -log10((FDR + 1e-300)))
  p_vol <- ggplot(volcano_df, aes(log2FC, negLog10FDR, color = (FDR <= PADJ_CUTOFF))) +
    geom_point(alpha = 0.6) +
    labs(title = "Volcano (by FDR)", x = "log2FC (A - B)", y = "-log10 FDR") +
    theme_bw()
  ggsave(file.path(OUTPUT_DIR, "volcano.pdf"), p_vol, width = 8, height = 6)
  ggsave(file.path(OUTPUT_DIR, "volcano.png"), p_vol, width = 8, height = 6, dpi = 300)
}

# Save session info and summary
writeLines(capture.output(sessionInfo()), con = file.path(OUTPUT_DIR, "session_info.txt"))
summary_txt <- c(
  paste0("MODE=", MODE),
  paste0("GROUPA samples: ", paste(groupA_samples, collapse = ",")),
  paste0("GROUPB samples: ", paste(groupB_samples, collapse = ",")),
  paste0("Kept features: ", nrow(counts_df_filt)),
  paste0("MIN_COUNT=", MIN_COUNT, " MIN_SAMPLES=", MIN_SAMPLES),
  paste0("PADJ_CUTOFF=", PADJ_CUTOFF, " MIN_LOG2FC=", MIN_LOG2FC)
)
if (exists("dispersion_info")) summary_txt <- c(summary_txt, paste0("DISPERSION_INFO: ", paste(capture.output(str(dispersion_info)), collapse=" ")))
writeLines(summary_txt, con = file.path(OUTPUT_DIR, "run_summary.txt"))

log_msg("Done. Results and plots in: ", normalizePath(OUTPUT_DIR))
