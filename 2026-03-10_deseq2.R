#!/usr/bin/env Rscript
#
# deseq2_counts_analysis_hardcoded.R
#
# Hard-coded local DESeq2 analysis script for count data produced by featureCounts or HTSeq.
# Edit the hard-coded paths and sample/condition mapping below, then run:
#   Rscript deseq2_counts_analysis_hardcoded.R
#
# Outputs (written to OUTPUT_DIR):
#  - dds.rds            : DESeqDataSet object after DESeq()
#  - vst.rds            : variance-stabilized transformed matrix (for plotting)
#  - deseq2_results.csv : full results table (all features)
#  - deseq2_sig.csv     : significant results (padj < PADJ)
#  - rlog_counts.csv    : VST-transformed counts (matrix)
#  - pca_plot.pdf/png
#  - heatmap_top_var.pdf
#  - ma_plot.pdf/png
#  - volcano_plot.pdf/png
#  - session_info.txt
#
# Notes:
# - This script is intentionally simple and designed for interactive runs on a laptop/desktop.
# - For production runs, convert to an argument-driven script and/or a workflow manager.
#
suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
  library(dplyr)
  library(tibble)
})

# ------------------------------
# USER: HARD-CODED CONFIGURATION
# Edit these values to match your environment
# ------------------------------
# Mode: "featurecounts" or "htseq"
COUNTS_MODE <- "featurecounts"

# If featurecounts: a single table (output of featureCounts -o)
FEATURECOUNTS_FILE <- "/path/to/featurecounts_matrix.txt"   # example: /home/user/counts/featurecounts.txt

# If htseq: directory containing one HTSeq counts file per sample (two-column tab text)
HTSEQ_DIR <- "/path/to/htseq_counts_dir"                    # example: /home/user/counts/htseq/

# Output directory for this run
OUTPUT_DIR <- "/path/to/deseq2_output"

# Sample to condition mapping (must exactly match sample names / column headers)
# Format: named vector where names are sample IDs and values are condition labels.
# For example: c("sampleA"="no_throm", "sampleB"="throm", ...)
SAMPLE_CONDITIONS <- c(
  "ai08nothrom1" = "no_thrombosis",
  "bq01nothrom1" = "no_thrombosis",
  "bq06nothrom1" = "no_thrombosis",
  "ai08throm1"   = "thrombosis",
  "bq01throm1"   = "thrombosis",
  "bq06throm1"   = "thrombosis"
)

# DE thresholds & plotting knobs
PADJ <- 0.05
LFC_MIN <- 1.0           # optional fold-change threshold for reporting/significant subsetting
TOP_N_HEATMAP <- 50      # top variable genes for heatmap
MIN_COUNTS_FILTER <- 10  # minimum counts threshold to keep a gene before DESeq

# Misc
SEED <- 42
THREADS <- 4

# ------------------------------
# Setup
# ------------------------------
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
set.seed(SEED)

logmsg <- function(...) message(sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste0(...)))

logmsg("OUTPUT_DIR = ", OUTPUT_DIR)
logmsg("COUNTS_MODE = ", COUNTS_MODE)

# ------------------------------
# Read counts depending on mode
# ------------------------------
if (COUNTS_MODE == "featurecounts") {
  if (!file.exists(FEATURECOUNTS_FILE)) stop("FEATURECOUNTS_FILE not found: ", FEATURECOUNTS_FILE)
  logmsg("Reading featureCounts table: ", FEATURECOUNTS_FILE)
  # featureCounts default: first column "Geneid", possibly "Length" second, then sample columns
  fc <- read.delim(FEATURECOUNTS_FILE, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
  # detect gene id column
  gid_col <- which(tolower(colnames(fc)) %in% c("geneid","gene_id","gene"))
  if (length(gid_col) == 0) gid_col <- 1
  gene_ids <- as.character(fc[[gid_col]])
  # find first sample column: skip typical annotation columns
  samp_cols <- setdiff(colnames(fc), colnames(fc)[1:which(colnames(fc)==colnames(fc)[gid_col])+1])
  # Better heuristic: keep columns that are numeric (counts)
  numeric_cols <- sapply(fc, is.numeric)
  count_cols <- colnames(fc)[numeric_cols]
  if (length(count_cols) >= 1) {
    counts_df <- fc[, count_cols, drop = FALSE]
  } else {
    # fallback: assume columns after first two
    counts_df <- fc[, (gid_col+2):ncol(fc), drop = FALSE]
  }
  rownames(counts_df) <- gene_ids
} else if (COUNTS_MODE == "htseq") {
  if (!dir.exists(HTSEQ_DIR)) stop("HTSEQ_DIR not found: ", HTSEQ_DIR)
  logmsg("Reading HTSeq count files from: ", HTSEQ_DIR)
  files <- list.files(HTSEQ_DIR, pattern = "\\.txt$|\\.counts$|\\.tsv$|\\.ct$", full.names = TRUE)
  if (length(files) == 0) stop("No count files found in HTSEQ_DIR: ", HTSEQ_DIR)
  # read each file: first column gene, second column count. Skip summary rows (like __no_feature)
  count_list <- list()
  sample_names <- character()
  for (f in files) {
    sname <- tools::file_path_sans_ext(basename(f))
    sample_names <- c(sample_names, sname)
    df <- read.delim(f, header = FALSE, stringsAsFactors = FALSE, col.names = c("gene","count"))
    # Keep only rows where gene doesn't start with "__"
    df <- df[!grepl("^__", df$gene), ]
    count_list[[sname]] <- df$count
    if (!exists("gene_order")) gene_order <- df$gene
    else {
      if (!all(gene_order == df$gene)) {
        warning("Gene order differs between HTSeq files; merging by gene name.")
      }
    }
    names(count_list[[sname]]) <- df$gene
  }
  # build merged matrix by union of genes
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

# ------------------------------
# Ensure sample names match SAMPLE_CONDITIONS keys
# ------------------------------
count_samples <- colnames(counts_df)
logmsg("Detected samples in counts: ", paste(count_samples, collapse = ", "))

if (!all(names(SAMPLE_CONDITIONS) %in% count_samples)) {
  # If SAMPLE_CONDITIONS names are not exact, try to match by substring
  unmatched <- setdiff(names(SAMPLE_CONDITIONS), count_samples)
  if (length(unmatched) > 0) {
    logmsg("Some SAMPLE_CONDITIONS keys not found in count columns: ", paste(unmatched, collapse = ", "))
    # attempt to match sample keys to column names by prefix or substring
    newmap <- character(length(SAMPLE_CONDITIONS))
    names(newmap) <- names(SAMPLE_CONDITIONS)
    for (s in names(SAMPLE_CONDITIONS)) {
      hits <- grep(s, count_samples, ignore.case = TRUE, value = TRUE)
      if (length(hits) == 1) {
        newmap[s] <- hits
      } else {
        newmap[s] <- ""
      }
    }
    # if mapping successful (no ""), remap SAMPLE_CONDITIONS
    if (all(nzchar(newmap))) {
      counts_df <- counts_df[, newmap]
      colnames(counts_df) <- names(newmap)  # restore desired names
      logmsg("Remapped count columns to SAMPLE_CONDITIONS names using substring heuristics.")
    } else {
      stop("Sample name mismatch. Ensure SAMPLE_CONDITIONS names exactly match count column names.")
    }
  }
} else {
  # reorder columns to match SAMPLE_CONDITIONS order
  counts_df <- counts_df[, names(SAMPLE_CONDITIONS), drop = FALSE]
}

# ------------------------------
# Build colData (sample table)
# ------------------------------
coldata <- data.frame(sample = colnames(counts_df),
                      condition = factor(as.vector(SAMPLE_CONDITIONS[colnames(counts_df)]),
                                         levels = unique(SAMPLE_CONDITIONS)),
                      row.names = colnames(counts_df),
                      stringsAsFactors = FALSE)
logmsg("Sample table:")
print(coldata)

# ------------------------------
# Prefilter low count genes (simple)
# ------------------------------
filter_keep <- rowSums(counts_df >= MIN_COUNTS_FILTER) >= 2
logmsg("Filtering genes: keeping ", sum(filter_keep), " features (out of ", nrow(counts_df), ") with >= ", MIN_COUNTS_FILTER, " counts in >=2 samples")
counts_df_filt <- counts_df[filter_keep, , drop = FALSE]

# ------------------------------
# Create DESeqDataSet and run DESeq
# ------------------------------
dds <- DESeqDataSetFromMatrix(countData = round(as.matrix(counts_df_filt)),
                              colData = coldata,
                              design = ~ condition)

# Re-level condition so reference is first level (if a specific reference desired, set it above)
# Example: relevel(dds$condition, ref = "no_thrombosis")
dds$condition <- relevel(dds$condition, ref = levels(dds$condition)[1])

# Run DESeq (Wald test by default)
logmsg("Running DESeq() ... (this may take a while)")
dds <- DESeq(dds, parallel = FALSE)

# Save dds
saveRDS(dds, file = file.path(OUTPUT_DIR, "dds.rds"))

# ------------------------------
# Transformations for plotting and PCA
# ------------------------------
logmsg("Computing variance stabilizing transformation (vst)")
vst <- vst(dds, blind = FALSE)
saveRDS(vst, file = file.path(OUTPUT_DIR, "vst.rds"))

vst_mat <- assay(vst)
write.csv(vst_mat, file = file.path(OUTPUT_DIR, "vst_counts_matrix.csv"), row.names = TRUE)

# PCA
pca_res <- prcomp(t(vst_mat))
percentVar <- round(100 * pca_res$sdev^2 / sum(pca_res$sdev^2), 1)
pca_df <- data.frame(PC1 = pca_res$x[,1], PC2 = pca_res$x[,2], sample = rownames(pca_res$x), condition = coldata$condition)

p_pca <- ggplot(pca_df, aes(PC1, PC2, color = condition, label = sample)) +
  geom_point(size = 3) +
  geom_text(hjust = 0.5, vjust = -1.2, size = 3) +
  labs(title = "PCA (VST)", x = paste0("PC1: ", percentVar[1], "%"), y = paste0("PC2: ", percentVar[2], "%")) +
  theme_bw(base_size = 12)

ggsave(filename = file.path(OUTPUT_DIR, "pca_plot.pdf"), plot = p_pca, width = 7, height = 6)
ggsave(filename = file.path(OUTPUT_DIR, "pca_plot.png"), plot = p_pca, width = 7, height = 6, dpi = 300)

# ------------------------------
# Differential results
# ------------------------------
# Compare second level vs first level (e.g., thrombosis vs no_thrombosis)
coef_name <- "condition"
res <- results(dds, alpha = PADJ)
res <- lfcShrink(dds, coef = resultsNames(dds)[2], type = if (requireNamespace("apeglm", quietly=TRUE)) "apeglm" else "normal")
res_df <- as.data.frame(res) %>% rownames_to_column(var = "feature")
res_df <- as_tibble(res_df)

# Add mean expression from vst
mean_vst <- rowMeans(vst_mat[res_df$feature, , drop = FALSE])
res_df$mean_vst <- mean_vst

# Save full results
write.csv(res_df, file = file.path(OUTPUT_DIR, "deseq2_results.csv"), row.names = FALSE)

# Subset significant
sig_df <- res_df %>% filter(!is.na(padj) & padj < PADJ & abs(log2FoldChange) >= LFC_MIN) %>% arrange(padj)
write.csv(sig_df, file = file.path(OUTPUT_DIR, "deseq2_significant.csv"), row.names = FALSE)
logmsg("Significant genes (padj<", PADJ, " & |LFC|>=", LFC_MIN, "): ", nrow(sig_df))

# Save normalized counts (VST) for convenience
write.csv(as.data.frame(vst_mat), file = file.path(OUTPUT_DIR, "vst_transformed_counts.csv"), row.names = TRUE)

# ------------------------------
# MA plot (using results)
# ------------------------------
ma_df <- res_df
ma_df$minusLog10Padj <- -log10(ma_df$padj + 1e-300)
p_ma <- ggplot(ma_df, aes(mean_vst, log2FoldChange, color = padj < PADJ)) +
  geom_point(alpha = 0.5, size = 1) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "grey60")) +
  labs(title = "MA plot (VST mean vs log2FC)", x = "Mean VST", y = "log2 fold change") +
  theme_bw(base_size = 12)

ggsave(file.path(OUTPUT_DIR, "ma_plot.pdf"), p_ma, width = 7, height = 5)
ggsave(file.path(OUTPUT_DIR, "ma_plot.png"), p_ma, width = 7, height = 5, dpi = 300)

# ------------------------------
# Volcano plot
# ------------------------------
volcano_df <- res_df
volcano_df$negLog10Padj <- -log10(volcano_df$padj + 1e-300)
volcano_df$significant <- volcano_df$padj < PADJ & abs(volcano_df$log2FoldChange) >= LFC_MIN

p_volcano <- ggplot(volcano_df, aes(log2FoldChange, negLog10Padj, color = significant)) +
  geom_point(alpha = 0.6, size = 1.2) +
  scale_color_manual(values = c("TRUE" = "#d73027", "FALSE" = "grey60")) +
  geom_vline(xintercept = c(-LFC_MIN, LFC_MIN), linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = -log10(PADJ), linetype = "dashed", color = "grey40") +
  labs(title = "Volcano plot", x = "log2 Fold Change", y = "-log10 adjusted p-value") +
  theme_bw(base_size = 12)

ggsave(file.path(OUTPUT_DIR, "volcano_plot.pdf"), p_volcano, width = 8, height = 6)
ggsave(file.path(OUTPUT_DIR, "volcano_plot.png"), p_volcano, width = 8, height = 6, dpi = 300)

# ------------------------------
# Heatmap of top variable genes
# ------------------------------
rv <- rowVars(vst_mat)
names(rv) <- rownames(vst_mat)
top_genes <- names(sort(rv, decreasing = TRUE))[1:min(TOP_N_HEATMAP, length(rv))]
heat_mat <- vst_mat[top_genes, , drop = FALSE]
# center rows
heat_mat <- t(scale(t(heat_mat), center = TRUE, scale = FALSE))
annotation_col <- data.frame(condition = coldata$condition)
rownames(annotation_col) <- colnames(heat_mat)

ann_colors <- list(condition = setNames(brewer.pal(min(8, length(unique(coldata$condition))), "Set1"),
                                       unique(coldata$condition)))

pheatmap(heat_mat, annotation_col = annotation_col, annotation_colors = ann_colors,
         show_rownames = TRUE, show_colnames = TRUE,
         filename = file.path(OUTPUT_DIR, "heatmap_top_var.pdf"), width = 8, height = 10)

# ------------------------------
# Save session info and objects
# ------------------------------
writeLines(capture.output(sessionInfo()), con = file.path(OUTPUT_DIR, "session_info.txt"))
saveRDS(res_df, file = file.path(OUTPUT_DIR, "deseq2_results_df.rds"))

logmsg("DESeq2 analysis complete. Outputs in: ", OUTPUT_DIR)
