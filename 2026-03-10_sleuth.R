#!/usr/bin/env Rscript
#
# sleuth_thrombosis_annotated.R
#
# Heavily annotated Sleuth differential expression analysis script:
# Thrombosis vs No Thrombosis using kallisto outputs.
#
# This script is a "patched" version of a previously provided analysis and
# includes extensive inline comments explaining each step, assumptions,
# common failure modes and suggestions for robustness, reproducibility,
# and follow-up analyses.
#
# Summary of workflow:
#  1. Validate sample directories and locate abundance.tsv / abundance.h5 files.
#  2. Optionally load tx2gene mapping and normalize transcript IDs (remove version suffixes).
#  3. Build a TPM matrix (transcript or gene level) for exploratory PCA/heatmaps and as mean expression.
#  4. Prepare sleuth sample table (s2c) and run sleuth_prep (with optional target mapping/aggregation).
#  5. Fit model and perform Wald test for the condition coefficient.
#  6. Extract sleuth results, add mean expression, perform multiple-testing adjustment if needed.
#  7. Save results and create diagnostic plots (PCA, volcano, MA, heatmaps).
#
# Key design choices and rationale:
#  - We normalize target IDs (strip version suffixes like ENST00000.1 -> ENST00000) to improve mapping
#    consistency between kallisto outputs and any tx2gene mapping you provide.
#  - We aggregate to gene-level TPMs if a tx2gene mapping is supplied. Aggregation via TPM-sum is
#    appropriate for exploratory plots and ranking; for differential testing one may prefer aggregated counts.
#  - We create a log2(TPM+1) transformed matrix (vst-like) for PCA and clustering; this is not a replacement
#    for a formal variance-stabilizing transform (e.g., DESeq2's vst) but is suitable for visuals.
#  - We add an optional low-expression filter to reduce noisy features and avoid LOESS warnings inside sleuth.
#
# Requirements:
#  - R (>= 4.0 recommended)
#  - R packages: sleuth, dplyr, tidyr, ggplot2, ggrepel, pheatmap, RColorBrewer
#  - Optional: rhdf5 (to read abundance.h5 if you don't produce abundance.tsv)
#
# Reproducibility:
#  - Save sessionInfo() output (this script writes session_info.txt).
#  - Pin package versions via renv or use a conda environment.
#
# Troubleshooting quick tips (read before running):
#  - If sleuth_prep fails to find abundance files, check sample directory names and that
#    kallisto outputs (abundance.tsv / abundance.h5) exist.
#  - If tx2gene mapping yields zero matches with kallisto target IDs, inspect ID formats:
#    Ensembl often appends version suffixes (.1, .2). This script removes those suffixes.
#  - If sleuth_reports many LOESS warnings, filter low-expression features or increase min_bootstrap (if applicable).

# -------------------------
# Load libraries (quietly)
# -------------------------
suppressPackageStartupMessages({
  # Core DE tool: sleuth (built for kallisto output)
  library(sleuth)
  # Tidyverse-style data handling
  library(dplyr)
  library(tidyr)
  # Plotting
  library(ggplot2)
  library(ggrepel)
  # Heatmaps
  library(pheatmap)
  library(RColorBrewer)
  # rhdf5: optional; required only if reading abundance.h5 files directly
  if (!requireNamespace("rhdf5", quietly = TRUE)) {
    # The script will still run if abundance.tsv exists; if you need HDF5 reading, install rhdf5
    # Bioconductor install: BiocManager::install("rhdf5")
    message("[NOTE] rhdf5 not available; abundance.h5 reading will fail if needed.")
  } else {
    library(rhdf5)
  }
})

# -------------------------
# Configuration (edit as needed)
# -------------------------
KALLISTO_DIR <- "/home/jgd/Documents/2026-04-01_proc-enu_bioinformatics_results/counts/kallisto"
OUTPUT_DIR   <- "/home/jgd/Documents/bioinformatics_working/output"

# Optional target -> gene mapping (two-column TSV). If NULL, transcript-level analysis only.
TX2GENE_FILE <- "/home/jgd/Documents/2026-04-01_proc-enu_bioinformatics_results/counts/kallisto/tx2gene.tsv"

# Reference level in the condition factor (used in design)
REFERENCE_LEVEL <- "no_thrombosis"

# Significance thresholds (tweakable)
PADJ_CUTOFF  <- 0.05
LFC_CUTOFF   <- 1.0
TOP_N_LABELS  <- 10

# Filtering low-expression features to avoid LOESS extrapolation warnings in sleuth
ENABLE_LOW_EXPR_FILTER <- TRUE
LOW_EXPR_TPM <- 0.1
MIN_SAMPLES_FOR_TPM <- 2

# -------------------------
# Sample assignments (must match sample directory names under KALLISTO_DIR)
# -------------------------
# This object maps sample directory names (kallisto output folder names) to condition labels.
# Keep these exact names in KALLISTO_DIR.
SAMPLE_CONDITIONS <- c(
  "bq06nothrom1" = "no_thrombosis",
  "bq01nothrom1" = "no_thrombosis",
  "ai08nothrom1" = "no_thrombosis",
  "bq06throm1"   = "thrombosis",
  "bq01throm1"   = "thrombosis",
  "ai08throm1"   = "thrombosis"
)

# -------------------------
# Setup output directory
# -------------------------
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

message("=== Sleuth Thrombosis Analysis (patched & annotated) ===")
message("Kallisto directory: ", KALLISTO_DIR)
message("Output directory:   ", OUTPUT_DIR)

# =============================================================================
# STEP 1: Validate samples and locate abundance.tsv / abundance.h5 files
# =============================================================================
message("\n[1/7] Validating samples and locating abundance.tsv/h5 files...")

clean_names <- names(SAMPLE_CONDITIONS)
missing_from_dir <- character(0)
sample_paths <- character(length(clean_names))
names(sample_paths) <- clean_names

# For each expected sample name, check if a directory exists and contains abundance.tsv or abundance.h5
for (s in clean_names) {
  tsv_path <- file.path(KALLISTO_DIR, s, "abundance.tsv")
  h5_path  <- file.path(KALLISTO_DIR, s, "abundance.h5")
  if (file.exists(tsv_path)) {
    sample_paths[s] <- tsv_path
  } else if (file.exists(h5_path)) {
    # point to directory when passing to sleuth_prep (it can read .h5 when given directory)
    sample_paths[s] <- file.path(KALLISTO_DIR, s)
  } else {
    missing_from_dir <- c(missing_from_dir, s)
  }
}

# If any samples are missing, stop with an informative message
if (length(missing_from_dir) > 0) {
  stop(
    "The following samples listed in SAMPLE_CONDITIONS were not found under KALLISTO_DIR\n(looking for abundance.tsv or abundance.h5):\n",
    paste(" ", missing_from_dir, collapse = "\n"),
    "\nEnsure KALLISTO_DIR is correct and each sample dir contains abundance.tsv or abundance.h5."
  )
}

message("  Samples to analyze (", length(clean_names), "):")
for (s in clean_names) message("    - ", s, "  [", SAMPLE_CONDITIONS[s], "]")

# =============================================================================
# STEP 2: Optional tx2gene mapping and ID normalization
# =============================================================================
message("\n[2/7] Loading optional tx2gene mapping and normalizing target IDs...")

tx2gene <- NULL
if (!is.null(TX2GENE_FILE)) {
  if (!file.exists(TX2GENE_FILE)) {
    stop("TX2GENE_FILE specified but not found: ", TX2GENE_FILE)
  }
  # Expect header or no-header with two columns. We attempt to be flexible.
  tx2gene <- read.delim(TX2GENE_FILE, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  if (ncol(tx2gene) < 2) stop("TX2GENE_FILE must have at least two columns: target_id, gene_id")
  tx2gene <- tx2gene[, 1:2]
  colnames(tx2gene) <- c("target_id", "gene_id")

  # Remove version suffixes (e.g., ENSDART00000000001.1 -> ENSDART00000000001) to match kallisto outputs
  tx2gene$target_id <- sub("\\.\\d+$", "", tx2gene$target_id)

  # Warn about ambiguous mappings: same target_id mapping to multiple gene_ids
  mm <- tx2gene %>% group_by(target_id) %>% summarise(n = n()) %>% filter(n > 1)
  if (nrow(mm) > 0) {
    warning("Some target_id values map to multiple gene_id values. Example (first 10): ",
            paste(head(paste0(mm$target_id, " (", mm$n, ")"), 10), collapse = ", "))
  }

  # Remove duplicate rows (keep first) and duplicate targets (keep first)
  tx2gene <- tx2gene[!duplicated(tx2gene), ]
  tx2gene <- tx2gene[!duplicated(tx2gene$target_id), ]
  message("  Loaded tx2gene mapping with ", nrow(tx2gene), " unique rows (targets normalized).")
} else {
  message("  No tx2gene mapping provided; performing transcript-level analysis.")
}

# =============================================================================
# STEP 3: Build TPM matrix for PCA/heatmaps and mean expression calculation
# =============================================================================
message("\n[3/7] Building TPM matrix for PCA and downstream plotting...")

# Helper to read tpm from an abundance.tsv or from abundance.h5 (best-effort)
read_tpm <- function(path_or_dir) {
  # path_or_dir can be a directory (kallisto output dir) or a direct path to abundance.tsv
  if (file.info(path_or_dir)$isdir) {
    tsv <- file.path(path_or_dir, "abundance.tsv")
    h5  <- file.path(path_or_dir, "abundance.h5")
    if (file.exists(tsv)) {
      path <- tsv
      from_h5 <- FALSE
    } else if (file.exists(h5)) {
      path <- h5
      from_h5 <- TRUE
    } else {
      stop("Expected abundance.tsv or abundance.h5 in directory: ", path_or_dir)
    }
  } else {
    path <- path_or_dir
    from_h5 <- grepl("\\.h5$", path, ignore.case = TRUE)
  }

  if (!from_h5) {
    # Read the plain text TSV produced by kallisto
    df <- read.delim(path, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
    # Kallisto may use different column names; be tolerant
    if (!("target_id" %in% colnames(df))) {
      if ("target" %in% colnames(df)) df$target_id <- df$target
      else if ("transcript_id" %in% colnames(df)) df$target_id <- df$transcript_id
      else stop("abundance.tsv missing 'target_id' column in: ", path)
    }
    if (!("tpm" %in% colnames(df))) stop("abundance.tsv missing 'tpm' column in: ", path)
    df <- df[, c("target_id", "tpm")]
  } else {
    # Attempt to read HDF5. This code does a best-effort autodetection of ids and tpm datasets.
    if (!requireNamespace("rhdf5", quietly = TRUE)) {
      stop("abundance.tsv not found and rhdf5 package not available to read abundance.h5. Install rhdf5 or generate TSVs.")
    }
    h5_path <- path
    # List HDF5 contents to find ids/tpm datasets (kallisto HDF5 structure can vary by version)
    contents <- rhdf5::h5ls(h5_path)
    # Heuristics to find ids and tpm datasets
    ids_row  <- contents %>% filter(name %in% c("ids", "targets", "target_ids", "target_id"))
    tpm_row  <- contents %>% filter(name %in% c("tpm", "TPM", "est_tpm"))
    ids_path <- NULL; tpm_path <- NULL
    if (nrow(ids_row) > 0) {
      ids_path <- paste0(ids_row$group[1], "/", ids_row$name[1])
    }
    if (nrow(tpm_row) > 0) {
      tpm_path <- paste0(tpm_row$group[1], "/", tpm_row$name[1])
    }
    if (is.null(ids_path) || is.null(tpm_path)) {
      stop("Could not auto-detect ids/tpm datasets in HDF5 file: ", h5_path,
           "\nHDF5 contents:\n", paste(capture.output(print(contents)), collapse = "\n"),
           "\nConsider generating abundance.tsv files with kallisto quant (default behavior).")
    }
    ids <- rhdf5::h5read(h5_path, ids_path)
    tpm <- rhdf5::h5read(h5_path, tpm_path)
    df <- data.frame(target_id = as.character(ids), tpm = as.numeric(tpm), stringsAsFactors = FALSE)
  }

  # Normalize target IDs by removing trailing version suffixes such as ".1", ".2"
  df$target_id <- sub("\\.\\d+$", "", df$target_id)
  return(df)
}

# Read TPM for each sample and assemble a wide matrix
tpm_list <- lapply(sample_paths, read_tpm)
names(tpm_list) <- clean_names

# Convert list-of-data.frames to a single wide data frame by full joining on target_id
tpm_wide <- Reduce(function(a, b) full_join(a, b, by = "target_id"),
                   lapply(names(tpm_list), function(nm) {
                     df <- tpm_list[[nm]]
                     colnames(df) <- c("target_id", nm)
                     df
                   }))

# Replace NAs (missing targets) with 0 TPM
tpm_wide[is.na(tpm_wide)] <- 0

# If tx2gene mapping provided, aggregate TPMs per gene by summing transcript TPMs
if (!is.null(tx2gene)) {
  message("  Aggregating transcripts to genes via tx2gene...")
  # Left join to bring gene_id alongside TPM columns
  tpm_wide <- left_join(tx2gene, tpm_wide, by = "target_id")
  # Drop rows lacking gene_id (e.g., targets that didn't match mapping)
  tpm_wide <- tpm_wide[!is.na(tpm_wide$gene_id), ]
  # Sum TPMs per gene (summing across transcripts)
  tpm_gene <- tpm_wide %>%
    select(-target_id) %>%
    group_by(gene_id) %>%
    summarise(across(everything(), sum))
  # Build matrix: rows = genes, cols = samples
  rownames_mat <- tpm_gene$gene_id
  tpm_mat <- as.matrix(tpm_gene[, -1])
  rownames(tpm_mat) <- rownames_mat
} else {
  # Transcript-level TPM matrix
  message("  Running transcript-level analysis (no tx2gene provided).")
  rownames_mat <- tpm_wide$target_id
  tpm_mat <- as.matrix(tpm_wide[, -1])
  rownames(tpm_mat) <- rownames_mat
}

# Ensure column names are the sample names in correct order
colnames(tpm_mat) <- clean_names

message("  TPM matrix dimensions: ", nrow(tpm_mat), " features x ", ncol(tpm_mat), " samples")

# Diagnostic checks: ensure mapping between TPM features and tx2gene if provided
if (!is.null(tx2gene)) {
  tpm_features <- rownames(tpm_mat)
  tx_targets   <- unique(tx2gene$target_id)
  matched <- sum(tpm_features %in% tx_targets)
  message("  TPM features: ", length(tpm_features))
  message("  tx2gene target_id entries: ", length(tx_targets))
  message("  Matched targets (TPM ∩ tx2gene): ", matched, " / ", length(tpm_features))
  if (matched == 0) {
    stop("No matches between kallisto target_ids and tx2gene target_id after normalization. Re-check ID formats.")
  }
}

# Optional low-expression filter: keep features with TPM >= LOW_EXPR_TPM in at least MIN_SAMPLES_FOR_TPM samples
if (ENABLE_LOW_EXPR_FILTER) {
  keep_fe <- rowSums(tpm_mat >= LOW_EXPR_TPM) >= MIN_SAMPLES_FOR_TPM
  message("  Low-expression filter: keeping ", sum(keep_fe), " features (removed ", sum(!keep_fe), ")")
  tpm_mat <- tpm_mat[keep_fe, , drop = FALSE]
}

# Log-transform for visualization (simple variance-stabilization-like transform)
vst_mat <- log2(tpm_mat + 1)

# =============================================================================
# STEP 4: Prepare sleuth sample table (s2c) and run sleuth_prep
# =============================================================================
message("\n[4/7] Preparing sleuth sample table and running sleuth_prep...")

# Build sample-to-condition table (s2c) expected by sleuth
s2c <- data.frame(sample = clean_names,
                  condition = factor(SAMPLE_CONDITIONS[clean_names],
                                     levels = c(REFERENCE_LEVEL, "thrombosis")),
                  path = as.character(sample_paths[clean_names]),
                  stringsAsFactors = FALSE)
rownames(s2c) <- s2c$sample

# sleuth expects a path to the kallisto output directory for each sample (directory containing abundance.h5)
# If we pointed to abundance.tsv earlier, derive the directory with dirname()
s2c$sleuth_path <- ifelse(file.info(s2c$path)$isdir, s2c$path, dirname(s2c$path))

# Prepare target_mapping only if tx2gene is present
target_mapping <- NULL
if (!is.null(tx2gene)) {
  target_mapping <- tx2gene
}

# Call sleuth_prep. We pass a simple ~ condition formula and (optionally) an aggregation column 'gene_id'
# so sleuth aggregates transcript-level estimates to gene-level. aggregation_column expects a column name
# in the provided target_mapping (in our tx2gene we named it 'gene_id').
so <- sleuth_prep(s2c %>% select(sample, condition, path = sleuth_path),
                  ~ condition,
                  target_mapping = target_mapping,
                  aggregation_column = if (!is.null(target_mapping)) "gene_id" else NULL,
                  extra_bootstrap_summary = FALSE,
                  verbose = TRUE)

# After sleuth_prep run, inspect target_mapping in 'so' to ensure mapping is correct.
message("\n[5/7] Inspecting so$target_mapping ...")
if (!is.null(so$target_mapping)) {
  tm <- so$target_mapping
  # normalize IDs there too for clarity
  if ("target_id" %in% colnames(tm)) {
    tm$target_id <- sub("\\.\\d+$", "", tm$target_id)
  }
  message("so$target_mapping: rows=", nrow(tm), ", columns=", ncol(tm))
  print(head(tm, 12))
  present_in_tpm <- sum(tm$target_id %in% rownames(vst_mat))
  message("Mapping targets present in TPM matrix: ", present_in_tpm, " / ", nrow(tm))
  unmatched_in_mapping <- setdiff(rownames(vst_mat), tm$target_id)
  if (length(unmatched_in_mapping) > 0) {
    message("Example TPM features not present in so$target_mapping (first 10): ",
            paste(head(unmatched_in_mapping, 10), collapse = ", "))
  } else {
    message("All TPM features are represented in so$target_mapping (based on names).")
  }
} else {
  message("so$target_mapping is NULL (transcript-level analysis or no mapping provided).")
}

# =============================================================================
# STEP 5: Fit model and test using sleuth (use s2c sample table for model matrix)
# =============================================================================
message("\n[6/7] Fitting sleuth model and running Wald test...")

# Fit model: specify full model formula name "full" (sleuth stores model under that name)
so <- sleuth_fit(so, ~ condition, "full")

# Use the s2c data.frame to construct model.matrix and determine which coefficient corresponds to condition
mm <- model.matrix(~ condition, data = s2c)
coef_candidates <- setdiff(colnames(mm), "(Intercept)")
if (length(coef_candidates) == 0) stop("Unable to determine coefficient name for condition. Check sample table.")
coef_name <- coef_candidates[1]  # typically "conditionthrombosis" if reference is no_thrombosis
message("  Using coefficient for test: ", coef_name)

# Run Wald test for the coefficient of interest
so <- sleuth_wt(so, which_beta = coef_name)

# =============================================================================
# STEP 6: Extract results, add mean expression, adjust p-values if needed
# =============================================================================
message("\n[7/7] Extracting results and computing summaries...")

# Extract results using sleuth_results; 'wt' corresponds to Wald test (sleuth_wt)
res_tbl <- sleuth_results(so, coef_name, "wt")

# Normalize column naming for feature identifiers. sleuth can name features differently
if ("target_id" %in% colnames(res_tbl)) {
  res_tbl <- res_tbl %>% rename(feature = target_id)
} else if ("feature" %in% colnames(res_tbl)) {
  # already ok
} else if ("target" %in% colnames(res_tbl)) {
  res_tbl <- res_tbl %>% rename(feature = target)
} else if ("gene_id" %in% colnames(res_tbl)) {
  res_tbl <- res_tbl %>% rename(feature = gene_id)
} else {
  # as a last resort, rename first column to 'feature'
  colnames(res_tbl)[1] <- "feature"
}

# sleuth returns 'b' for effect size (approx log fold change) and 'pval' and/or 'qval'
if (!("b" %in% colnames(res_tbl))) stop("sleuth_results missing 'b' column (effect size).")
if (!("qval" %in% colnames(res_tbl))) {
  warning("sleuth_results missing 'qval' column. Falling back to BH correction on 'pval'.")
  if (!("pval" %in% colnames(res_tbl))) stop("sleuth_results missing both 'pval' and 'qval'. Cannot continue.")
  res_tbl$qval <- p.adjust(res_tbl$pval, method = "BH")
}

# Add mean expression information from vst_mat if the feature identifiers overlap
res_tbl$mean_log2TPM <- NA_real_
present_features <- intersect(res_tbl$feature, rownames(vst_mat))
if (length(present_features) > 0) {
  # compute rowMeans only for those features present
  res_tbl$mean_log2TPM[res_tbl$feature %in% present_features] <-
    rowMeans(vst_mat[res_tbl$feature[res_tbl$feature %in% present_features], , drop = FALSE])
}

message("  Sleuth results rows: ", nrow(res_tbl))

# -------------------------
# Save results and create plots
# -------------------------
message("\nSaving results and generating plots in: ", OUTPUT_DIR)

out_df <- res_tbl %>% rename(log2FoldChange = b, padj = qval) %>% arrange(padj)
write.csv(out_df, file.path(OUTPUT_DIR, "sleuth_all_results.csv"), row.names = FALSE)
message("  Saved: sleuth_all_results.csv  (", nrow(out_df), " features)")

# Significant results based on thresholds set near top
sig_df <- out_df %>% filter(!is.na(padj) & padj < PADJ_CUTOFF & abs(log2FoldChange) >= LFC_CUTOFF) %>% arrange(padj)
write.csv(sig_df, file.path(OUTPUT_DIR, "sleuth_significant_results.csv"), row.names = FALSE)
message("  Saved: sleuth_significant_results.csv  (", nrow(sig_df), " features)")

# Split up and down regulated
up_df   <- sig_df %>% filter(log2FoldChange > 0)
down_df <- sig_df %>% filter(log2FoldChange < 0)
write.csv(up_df,   file.path(OUTPUT_DIR, "sleuth_upregulated_in_thrombosis.csv"), row.names = FALSE)
write.csv(down_df, file.path(OUTPUT_DIR, "sleuth_downregulated_in_thrombosis.csv"), row.names = FALSE)
message("  Saved: sleuth_upregulated_in_thrombosis.csv   (", nrow(up_df), " features)")
message("  Saved: sleuth_downregulated_in_thrombosis.csv (", nrow(down_df), " features)")

# -------------------------
# PCA plot (log2 TPM matrix)
# -------------------------
# We use the vst_mat (log2(TPM+1)) for PCA - transpose so samples are rows
pca_res <- prcomp(t(vst_mat), center = TRUE, scale. = FALSE)
percent_var <- round(100 * (pca_res$sdev^2 / sum(pca_res$sdev^2)), 1)
pca_df <- data.frame(PC1 = pca_res$x[,1], PC2 = pca_res$x[,2], sample = rownames(pca_res$x),
                     condition = as.character(s2c$condition[rownames(pca_res$x)]), stringsAsFactors = FALSE)

p_pca <- ggplot(pca_df, aes(PC1, PC2, color = condition, label = sample)) +
  geom_point(size = 4, alpha = 0.85) +
  geom_text_repel(size = 3, max.overlaps = 20) +
  scale_color_manual(values = c("no_thrombosis" = "#4575b4", "thrombosis" = "#d73027")) +
  labs(title = "PCA - Thrombosis vs. No Thrombosis",
       x = paste0("PC1: ", percent_var[1], "% variance"),
       y = paste0("PC2: ", percent_var[2], "% variance"),
       color = "Condition") +
  theme_bw(base_size = 13)

ggsave(file.path(OUTPUT_DIR, "pca_plot.pdf"), p_pca, width = 7, height = 5)
ggsave(file.path(OUTPUT_DIR, "pca_plot.png"), p_pca, width = 7, height = 5, dpi = 300)
message("  Saved: pca_plot.pdf / .png")

# -------------------------
# Volcano plot
# -------------------------
volcano_df <- out_df
volcano_df <- volcano_df[!is.na(volcano_df$padj), ]
volcano_df$significance <- "Not significant"
volcano_df$significance[volcano_df$padj < PADJ_CUTOFF & volcano_df$log2FoldChange >=  LFC_CUTOFF] <- "Up in Thrombosis"
volcano_df$significance[volcano_df$padj < PADJ_CUTOFF & volcano_df$log2FoldChange <= -LFC_CUTOFF] <- "Down in Thrombosis"
top_features <- head(volcano_df[order(volcano_df$padj), "feature"], TOP_N_LABELS)
volcano_df$label <- ifelse(volcano_df$feature %in% top_features, volcano_df$feature, NA)

p_volcano <- ggplot(volcano_df, aes(log2FoldChange, -log10(padj), color = significance, label = label)) +
  geom_point(alpha = 0.5, size = 1.5) +
  geom_text_repel(size = 2.8, na.rm = TRUE, max.overlaps = 30, segment.size = 0.3) +
  geom_vline(xintercept = c(-LFC_CUTOFF, LFC_CUTOFF), linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = -log10(PADJ_CUTOFF), linetype = "dashed", color = "grey40") +
  scale_color_manual(values = c("Not significant" = "grey60", "Up in Thrombosis" = "#d73027", "Down in Thrombosis" = "#4575b4")) +
  labs(title = "Volcano Plot - Thrombosis vs. No Thrombosis",
       x = "log2 Fold Change (Thrombosis / No Thrombosis)",
       y = "-log10 Adjusted P-value",
       color = NULL) +
  theme_bw(base_size = 13)

ggsave(file.path(OUTPUT_DIR, "volcano_plot.pdf"), p_volcano, width = 8, height = 6)
ggsave(file.path(OUTPUT_DIR, "volcano_plot.png"), p_volcano, width = 8, height = 6, dpi = 300)
message("  Saved: volcano_plot.pdf / .png")

# -------------------------
# MA plot (Mean vs log2FC)
# -------------------------
pdf(file.path(OUTPUT_DIR, "ma_plot.pdf"), width = 7, height = 5)
plot(volcano_df$mean_log2TPM, volcano_df$log2FoldChange,
     xlab = "Mean log2(TPM + 1)", ylab = "log2 Fold Change",
     main = "MA Plot - Thrombosis vs. No Thrombosis",
     pch = 20, col = ifelse(volcano_df$padj < PADJ_CUTOFF, "#d73027", "grey60"))
abline(h = c(-LFC_CUTOFF, LFC_CUTOFF), lty = 2, col = "grey40")
dev.off()
message("  Saved: ma_plot.pdf")

# -------------------------
# Heatmap of top features (if any)
# -------------------------
if (nrow(sig_df) > 0) {
  top50 <- head(sig_df$feature, 50)
  top50_present <- intersect(top50, rownames(vst_mat))
  if (length(top50_present) == 0) {
    message("  No top features found in TPM matrix — skipping heatmap.")
  } else {
    mat <- vst_mat[top50_present, , drop = FALSE]
    # center rows for heatmap (z-score not used here; centering is common)
    mat <- mat - rowMeans(mat)
    annotation_col <- data.frame(Condition = s2c$condition)
    rownames(annotation_col) <- s2c$sample
    ann_colors <- list(Condition = c(no_thrombosis = "#4575b4", thrombosis = "#d73027"))
    pdf(file.path(OUTPUT_DIR, "heatmap_top50.pdf"), width = 10, height = 12)
    pheatmap(mat, annotation_col = annotation_col, annotation_colors = ann_colors,
             show_rownames = TRUE, show_colnames = TRUE, cluster_rows = TRUE, cluster_cols = TRUE,
             scale = "none", color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
             main = "Top DE Features - Thrombosis vs. No Thrombosis", fontsize_row = 7)
    dev.off()
    message("  Saved: heatmap_top50.pdf")
  }
} else {
  message("  Skipping heatmap — no significant features found.")
}

# -------------------------
# Sample distance heatmap
# -------------------------
sample_dists <- dist(t(vst_mat))
dist_mat <- as.matrix(sample_dists)
pdf(file.path(OUTPUT_DIR, "sample_distance_heatmap.pdf"), width = 8, height = 7)
pheatmap(dist_mat, clustering_distance_rows = sample_dists, clustering_distance_cols = sample_dists,
         col = colorRampPalette(rev(brewer.pal(9, "Blues")))(255), main = "Sample Distance Heatmap")
dev.off()
message("  Saved: sample_distance_heatmap.pdf")

# -------------------------
# Session info for reproducibility
# -------------------------
message("\nWriting session info...")
sink(file.path(OUTPUT_DIR, "session_info.txt"))
sessionInfo()
sink()

message("\n=== Analysis Complete ===")
message("Results in: ", OUTPUT_DIR)
message("  Significant features (padj < ", PADJ_CUTOFF, ", |LFC| >= ", LFC_CUTOFF, "): ", nrow(sig_df))
message("  Up in thrombosis:   ", nrow(up_df))
message("  Down in thrombosis: ", nrow(down_df))

# -------------------------
# Final notes and suggested improvements
# -------------------------
# - If you plan to do formal gene-level DE testing, consider aggregating counts (not TPMs) and running
#   DESeq2 or edgeR with appropriate normalization. Sleuth works with kallisto's estimated counts,
#   but many users prefer transcript-level methods followed by tximport summarization to counts.
#
# - For many samples or complex designs:
#     * Add covariates (batch, library_prep) to the model: ~ batch + condition
#     * Check sample distances and PCA for batch/confounders before testing.
#
# - For very low-expression features, LOESS fits in sleuth can issue warnings; we filter them here,
#   but another approach is to increase bootstrap/robustness parameters or use different variance modeling.
#
# - Consider saving the 'so' sleuth object (saveRDS(so, "so.rds")) for later inspection or reproducibility.
#
# - If you have abundance.h5 files only and rhdf5 cannot read them, re-run kallisto quant with --plaintext
#   to generate abundance.tsv files that are easier to parse.
#
# - If you want per-gene transcripts breakdown for significant genes, extract transcript-level results
#   from the sleuth object prior to aggregation (so$obs or sleuth_table(so, '...')).
#
# End of script.
