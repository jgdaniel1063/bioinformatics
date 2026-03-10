#!/usr/bin/env Rscript
#
# edgeR_counts_analysis_hardcoded.R
#
# Hard-coded, heavily commented edgeR analysis script for count data produced
# by featureCounts or HTSeq. Edit the hard-coded paths and sample/condition
# mapping below, then run:
#   Rscript edgeR_counts_analysis_hardcoded.R
#
# Features:
#  - Reads featureCounts matrix or HTSeq per-sample files
#  - Builds DGEList, normalizes, filters low-count genes
#  - Optionally estimates dispersion when replicates exist
#  - If no replicates are available, allows using a manual BCV (biological coefficient of variation)
#    to set dispersion and continue testing (recommended only when no replicates)
#  - Produces results table, significant subset, normalized/log-CPM matrices, PCA, heatmap, MA, volcano plots
#  - Saves sessionInfo and intermediate objects for reproducibility
#
# WARNING about dispersion:
#  - Reliable dispersion estimation requires replicates. If you have no replicates
#    per group, edgeR cannot estimate dispersion from the data; you must supply a
#    plausible BCV (biological coefficient of variation). Typical BCV values:
#      * technical replicates: ~0.01
#      * model organisms / well-controlled: ~0.1
#      * human clinical samples: ~0.4
#    Choose carefully; results depend strongly on this choice.
#
suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
  library(dplyr)
  library(tibble)
})

# ============================
# USER: HARD-CODED CONFIGURATION
# ============================
COUNTS_MODE <- "featurecounts"   # "featurecounts" or "htseq"

FEATURECOUNTS_FILE <- "/path/to/featurecounts_matrix.txt"   # featureCounts -o output
HTSEQ_DIR <- "/path/to/htseq_counts_dir"                    # directory with one HTSeq file per sample

OUTPUT_DIR <- "/path/to/edgeR_output"                       # where results will be written

# Named vector mapping sample names (column names in counts) to conditions
SAMPLE_CONDITIONS <- c(
  "ai08nothrom1" = "no_thrombosis",
  "bq01nothrom1" = "no_thrombosis",
  "bq06nothrom1" = "no_thrombosis",
  "ai08throm1"   = "thrombosis",
  "bq01throm1"   = "thrombosis",
  "bq06throm1"   = "thrombosis"
)

# Dispersion options
AUTO_ESTIMATE_DISP <- TRUE     # if TRUE and replicates exist, estimate dispersion from data
DISPERSION_MANUAL_BCV <- 0.20  # BCV to use when no replicates exist (set to NA to force estimate only)
# Note: dispersion = BCV^2 when passing to edgeR glmFit as numeric scalar

# Filtering & DE thresholds
MIN_COUNTS <- 10               # keep genes with >= MIN_COUNTS in at least MIN_SAMPLES_KEEP samples
MIN_SAMPLES_KEEP <- 2
LOGFC_CUTOFF <- 1.0
FDR_CUTOFF <- 0.05

# Plotting & analysis knobs
TOP_N_HEATMAP <- 50
SEED <- 2024

# ============================
# SETUP
# ============================
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
set.seed(SEED)

logmsg <- function(...) message(sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste0(...)))
logmsg("edgeR analysis starting. OUTPUT_DIR=", OUTPUT_DIR)

# ============================
# READ COUNTS
# ============================
if (COUNTS_MODE == "featurecounts") {
  if (!file.exists(FEATURECOUNTS_FILE)) stop("FEATURECOUNTS_FILE not found: ", FEATURECOUNTS_FILE)
  logmsg("Reading featureCounts table: ", FEATURECOUNTS_FILE)
  fc <- read.delim(FEATURECOUNTS_FILE, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
  # detect GeneID-like column
  gid_col <- which(tolower(colnames(fc)) %in% c("geneid","gene_id","gene"))
  if (length(gid_col) == 0) gid_col <- 1
  gene_ids <- as.character(fc[[gid_col]])
  # Heuristic: columns that are numeric are counts
  numeric_cols <- sapply(fc, is.numeric)
  count_cols <- colnames(fc)[numeric_cols]
  if (length(count_cols) >= 1) {
    counts_df <- fc[, count_cols, drop = FALSE]
  } else {
    # fallback: columns after first two
    counts_df <- fc[, (gid_col+2):ncol(fc), drop = FALSE]
  }
  rownames(counts_df) <- gene_ids
} else if (COUNTS_MODE == "htseq") {
  if (!dir.exists(HTSEQ_DIR)) stop("HTSEQ_DIR not found: ", HTSEQ_DIR)
  logmsg("Reading HTSeq count files from: ", HTSEQ_DIR)
  files <- list.files(HTSEQ_DIR, pattern = "\\.txt$|\\.counts$|\\.tsv$|\\.ct$", full.names = TRUE)
  if (length(files) == 0) stop("No HTSeq count files found in HTSEQ_DIR: ", HTSEQ_DIR)
  count_list <- list()
  for (f in files) {
    sname <- tools::file_path_sans_ext(basename(f))
    df <- read.delim(f, header = FALSE, stringsAsFactors = FALSE, col.names = c("gene","count"))
    df <- df[!grepl("^__", df$gene), ]   # remove HTSeq summary lines
    count_list[[sname]] <- df$count
    names(count_list[[sname]]) <- df$gene
  }
  all_genes <- sort(unique(unlist(lapply(count_list, names))))
  counts_mat <- matrix(0L, nrow = length(all_genes), ncol = length(count_list),
                       dimnames = list(all_genes, names(count_list)))
  for (nm in names(count_list)) {
    counts_mat[names(count_list[[nm]]), nm] <- as.integer(count_list[[nm]])
  }
  counts_df <- as.data.frame(counts_mat)
} else {
  stop("COUNTS_MODE must be 'featurecounts' or 'htseq'")
}

# Ensure column order matches SAMPLE_CONDITIONS
count_samples <- colnames(counts_df)
logmsg("Detected samples in counts: ", paste(count_samples, collapse = ", "))

if (!all(names(SAMPLE_CONDITIONS) %in% count_samples)) {
  unmatched <- setdiff(names(SAMPLE_CONDITIONS), count_samples)
  if (length(unmatched) > 0) {
    logmsg("Some SAMPLE_CONDITIONS keys not found: ", paste(unmatched, collapse = ", "))
    # Try to map by substring (best-effort)
    mapping <- character(length(SAMPLE_CONDITIONS))
    names(mapping) <- names(SAMPLE_CONDITIONS)
    for (s in names(SAMPLE_CONDITIONS)) {
      hits <- grep(s, count_samples, ignore.case = TRUE, value = TRUE)
      if (length(hits) == 1) mapping[s] <- hits
      else mapping[s] <- ""
    }
    if (all(nzchar(mapping))) {
      counts_df <- counts_df[, mapping]
      colnames(counts_df) <- names(mapping)
      logmsg("Remapped count columns to SAMPLE_CONDITIONS names using substring heuristics.")
    } else {
      stop("Sample name mismatch. Make SAMPLE_CONDITIONS keys match column names in counts.")
    }
  } else {
    counts_df <- counts_df[, names(SAMPLE_CONDITIONS), drop = FALSE]
  }
} else {
  counts_df <- counts_df[, names(SAMPLE_CONDITIONS), drop = FALSE]
}

# ============================
# Build DGEList and filter
# ============================
group <- factor(as.vector(SAMPLE_CONDITIONS[colnames(counts_df)]))
logmsg("Group levels: ", paste(levels(group), collapse = ", "))

dge <- DGEList(counts = as.matrix(counts_df), group = group)
dge <- calcNormFactors(dge)

# Filter low counts
keep <- rowSums(cpm(dge) * 1 >= (MIN_COUNTS / 1)) >= MIN_SAMPLES_KEEP  # uses CPM heuristic; alternative: raw counts
# Simpler: keep rows with counts >= MIN_COUNTS in at least MIN_SAMPLES_KEEP samples
keep2 <- rowSums(dge$counts >= MIN_COUNTS) >= MIN_SAMPLES_KEEP
keep_final <- keep | keep2
logmsg("Filtering: kept ", sum(keep_final), " / ", nrow(dge$counts), " genes (MIN_COUNTS=", MIN_COUNTS, " in >=", MIN_SAMPLES_KEEP, " samples)")

dge <- dge[keep_final,, keep.lib.sizes=FALSE]
dge <- calcNormFactors(dge)

# ============================
# Design matrix
# ============================
design <- model.matrix(~ group)
colnames(design) <- gsub("group", "", colnames(design))
logmsg("Design matrix columns: ", paste(colnames(design), collapse = ", "))

# Determine replication structure
rep_table <- table(group)
logmsg("Per-group sample counts: ", paste(names(rep_table), rep_table, sep=":", collapse=", "))
has_replicates <- any(rep_table >= 2)

# ============================
# Dispersion estimation or manual BCV
# ============================
dispersion_used <- NA_real_
fit <- NULL
lrt <- NULL

if (AUTO_ESTIMATE_DISP && has_replicates) {
  logmsg("Estimating dispersion from data (some groups have replicates).")
  # estimate common/trended/tagwise dispersion
  dge <- estimateDisp(dge, design, robust=TRUE)
  dispersion_used <- "estimated"
  # Fit GLM and perform LRT on coefficient 2 (group effect)
  fit <- glmFit(dge, design)
  # Pick coef index for group effect (usually 2)
  coef_index <- 2
  lrt <- glmLRT(fit, coef = coef_index)
} else {
  if (!is.na(DISPERSION_MANUAL_BCV) && DISPERSION_MANUAL_BCV > 0) {
    # No replicates to estimate dispersion, use manual BCV
    disp_value <- DISPERSION_MANUAL_BCV^2
    logmsg("No replicates for reliable dispersion estimation; using manual BCV=", DISPERSION_MANUAL_BCV, " -> dispersion=", signif(disp_value,3))
    dispersion_used <- paste0("manual_BCV=", DISPERSION_MANUAL_BCV)
    # Fit GLM using provided dispersion (scalar)
    fit <- glmFit(dge, design, dispersion = disp_value)
    coef_index <- 2
    lrt <- glmLRT(fit, coef = coef_index)
  } else {
    # Cannot proceed
    stop("No replicates in groups and no DISPERSION_MANUAL_BCV provided. Set DISPERSION_MANUAL_BCV to a plausible BCV (e.g., 0.2).")
  }
}

# ============================
# Extract results
# ============================
top <- topTags(lrt, n = Inf)$table
top <- as.data.frame(top) %>% rownames_to_column(var = "feature")
# top contains logFC (log2 fold change), logCPM, LR, PValue, FDR

# Add mean expression (logCPM can be used)
# Compute logCPM matrix for all genes (after filtering)
logcpm_mat <- cpm(dge, log=TRUE, prior.count=2)
mean_logcpm <- rowMeans(logcpm_mat[ top$feature , , drop = FALSE], na.rm = TRUE)
top$mean_logCPM <- mean_logcpm

# Save tables
write.csv(top, file = file.path(OUTPUT_DIR, "edgeR_results_all.csv"), row.names = FALSE)
sig <- top %>% filter(!is.na(FDR) & FDR < FDR_CUTOFF & abs(logFC) >= LOGFC_CUTOFF)
write.csv(sig, file = file.path(OUTPUT_DIR, "edgeR_results_significant.csv"), row.names = FALSE)
logmsg("Significant genes (FDR<", FDR_CUTOFF, " & |logFC|>=", LOGFC_CUTOFF, "): ", nrow(sig))

# Save DGEList and fit objects for later inspection
saveRDS(dge, file = file.path(OUTPUT_DIR, "dge.rds"))
saveRDS(fit, file = file.path(OUTPUT_DIR, "glmfit.rds"))
saveRDS(lrt, file = file.path(OUTPUT_DIR, "lrt.rds"))

# ============================
# Normalized counts and PCA
# ============================
# Use logCPM for downstream plots
write.csv(as.data.frame(logcpm_mat), file = file.path(OUTPUT_DIR, "logCPM_matrix.csv"), row.names = TRUE)

# PCA using logCPM (transpose so samples are rows)
pca_res <- prcomp(t(logcpm_mat))
percentVar <- round(100 * (pca_res$sdev^2 / sum(pca_res$sdev^2)), 1)
pca_df <- data.frame(PC1 = pca_res$x[,1], PC2 = pca_res$x[,2], sample = rownames(pca_res$x), condition = group)
p_pca <- ggplot(pca_df, aes(PC1, PC2, color = condition, label = sample)) +
  geom_point(size = 3) + geom_text_repel(size = 3, max.overlaps = 20) +
  labs(title = "PCA (logCPM)", x = paste0("PC1: ", percentVar[1], "%"), y = paste0("PC2: ", percentVar[2], "%")) +
  theme_bw(base_size = 12)
ggsave(file.path(OUTPUT_DIR, "PCA_logCPM.pdf"), p_pca, width = 7, height = 6)
ggsave(file.path(OUTPUT_DIR, "PCA_logCPM.png"), p_pca, width = 7, height = 6, dpi = 300)

# ============================
# MA plot
# ============================
ma_df <- top
ma_df$negLog10FDR <- -log10(ma_df$FDR + 1e-300)
p_ma <- ggplot(ma_df, aes(mean_logCPM, logFC, color = FDR < FDR_CUTOFF)) +
  geom_point(alpha = 0.5, size = 1) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "grey60")) +
  labs(title = "MA plot (edgeR)", x = "Mean logCPM", y = "log2 fold change") +
  theme_bw(base_size = 12)
ggsave(file.path(OUTPUT_DIR, "MA_plot.pdf"), p_ma, width = 7, height = 5)
ggsave(file.path(OUTPUT_DIR, "MA_plot.png"), p_ma, width = 7, height = 5, dpi = 300)

# ============================
# Volcano plot (use FDR)
# ============================
volcano_df <- top
volcano_df$negLog10FDR <- -log10(volcano_df$FDR + 1e-300)
volcano_df$significant <- volcano_df$FDR < FDR_CUTOFF & abs(volcano_df$logFC) >= LOGFC_CUTOFF
p_vol <- ggplot(volcano_df, aes(logFC, negLog10FDR, color = significant)) +
  geom_point(alpha = 0.6, size = 1.2) +
  scale_color_manual(values = c("TRUE" = "#d73027", "FALSE" = "grey60")) +
  geom_vline(xintercept = c(-LOGFC_CUTOFF, LOGFC_CUTOFF), linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = -log10(FDR_CUTOFF), linetype = "dashed", color = "grey40") +
  labs(title = "Volcano plot (edgeR)", x = "log2 Fold Change", y = "-log10 FDR") +
  theme_bw(base_size = 12)
ggsave(file.path(OUTPUT_DIR, "volcano_plot.pdf"), p_vol, width = 8, height = 6)
ggsave(file.path(OUTPUT_DIR, "volcano_plot.png"), p_vol, width = 8, height = 6, dpi = 300)

# ============================
# Heatmap of top variable genes (by SD on logCPM)
# ============================
rv <- apply(logcpm_mat, 1, sd, na.rm = TRUE)
top_genes <- names(sort(rv, decreasing = TRUE))[1:min(TOP_N_HEATMAP, length(rv))]
heat_mat <- logcpm_mat[top_genes, , drop = FALSE]
heat_mat <- t(scale(t(heat_mat), center = TRUE, scale = FALSE))  # center by gene
annotation_col <- data.frame(condition = group)
rownames(annotation_col) <- colnames(heat_mat)
ann_colors <- list(condition = setNames(brewer.pal(min(8, length(unique(group))), "Set1"), unique(group)))
pheatmap(heat_mat, annotation_col = annotation_col, annotation_colors = ann_colors,
         show_rownames = TRUE, show_colnames = TRUE,
         filename = file.path(OUTPUT_DIR, "heatmap_top_var.pdf"), width = 8, height = 10)

# ============================
# Save session info and provenance
# ============================
writeLines(capture.output(sessionInfo()), con = file.path(OUTPUT_DIR, "session_info.txt"))
write.csv(as.data.frame(top), file = file.path(OUTPUT_DIR, "edgeR_results_all_with_meanlogCPM.csv"), row.names = FALSE)

logmsg("edgeR analysis complete. Outputs written to: ", OUTPUT_DIR)
logmsg("Dispersion used: ", dispersion_used)
