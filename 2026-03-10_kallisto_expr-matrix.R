#!/usr/bin/env Rscript
#
# tximport_kallisto_expr_annotated.R
#
# Purpose
# -------
# Import transcript-level quantifications produced by Kallisto (abundance.tsv or abundance.h5),
# collapse to gene-level counts/TPM using tximport and tx-to-gene mapping, normalize using edgeR's
# TMM, compute log-CPM (log-counts per million) and write a gene x sample expression matrix (TSV).
#
# This heavily annotated version documents assumptions, failure modes, parameter choices,
# and suggestions for improvements/reproducibility. Edit the top-of-file paths to match
# your environment or convert them to CLI options (recommended).
#
# Outputs
# -------
# - out_expr: gene x sample expression matrix (log-CPM) written as TSV
# - out_sample_table: table mapping sample_id -> kallisto path (useful provenance)
#
# Requirements
# ------------
# - R (>= 4.0 recommended)
# - Bioconductor packages: tximport, edgeR
# - CRAN package: readr (used indirectly here; base write.table is used)
#
# Notes about tximport/countsFromAbundance:
# - countsFromAbundance = "lengthScaledTPM" attempts to correct estimated counts by transcript length
#   and the library size (via TPM) to produce gene-level counts that are more comparable across
#   samples than raw estimated counts. This is appropriate when you want gene-level counts for
#   differential expression methods that expect counts (like edgeR/DESeq2), but it is not identical
#   to true fragment counts.
# - Alternatives:
#     * "no": use the raw estimated counts from kallisto (may be biased by length)
#     * "scaledTPM": scale TPM to library size but not length scale (older behavior)
# ----------------------------
# USER-EDITABLE PATHS (top)
# ----------------------------
# Change these to match your filesystem. Consider converting to command-line
# arguments using 'argparse' or 'optparse' for pipeline flexibility.
kallisto_dir <- "/home/jgd/Documents/bioinformatics_working/counts/kallisto_quant_ak-cwb"
tx2gene_file <- "/home/jgd/Documents/bioinformatics_working/ref/tx2gene.tsv"

# NOTE: Original script wrote outputs to /Users/jeffdaniel/...; below paths
# likely should be changed to a project output directory. Update as needed.
out_expr <- "/Users/jeffdaniel/Documents/bioinformatics_working/output/expr_matrix.tsv"
out_sample_table <- "/Users/jeffdaniel/Documents/bioinformatics_working/output/sample_table.tsv"

# ----------------------------
# Package installation (best-effort)
# ----------------------------
# The script attempts to install missing dependencies automatically. In production
# pipelines, prefer using a pinned environment (conda, renv, Docker) rather than
# installing packages at runtime — installation may prompt / be slow or fail on CI.
required_bioc <- c("tximport","edgeR")
required_cran <- c("readr")

# BiocManager is used to install Bioconductor packages if they're missing.
if (!requireNamespace("BiocManager", quietly=TRUE)) {
  install.packages("BiocManager", repos="https://cloud.r-project.org")
}

for (p in required_bioc) {
  # requireNamespace is non-intrusive and doesn't attach the package; we install only if not available.
  if (!suppressWarnings(requireNamespace(p, quietly=TRUE))) {
    # In non-interactive contexts, BiocManager::install(..., ask=FALSE) is appropriate.
    BiocManager::install(p, ask=FALSE, update=FALSE)
  }
}
for (p in required_cran) {
  if (!suppressWarnings(requireNamespace(p, quietly=TRUE))) {
    install.packages(p, repos="https://cloud.r-project.org")
  }
}

# ----------------------------
# Load libraries
# ----------------------------
# Using library() attaches packages — tximport supplies the import function; edgeR for TMM and cpm.
library(tximport)
library(readr)   # readr isn't used below heavily, but is often useful; keep for later extension.
library(edgeR)

# ----------------------------
# Discover kallisto output files
# ----------------------------
# This assumes a typical output layout: each sample has a directory containing 'abundance.tsv'
# (text) or 'abundance.h5' (binary). We collect both and use tximport's auto-detection.
files <- list.files(path=kallisto_dir, pattern="abundance.tsv$|abundance.h5$", recursive=TRUE, full.names=TRUE)

# Safety checks
if (length(files) == 0) {
  stop("No kallisto 'abundance.tsv' or 'abundance.h5' files found under: ", kallisto_dir,
       "\nCheck kallisto_dir and the file organization.")
}

# Heuristic to derive sample IDs:
# We take the basename of the parent directory for each abundance file. This is common
# when kallisto results are in <sample>/abundance.tsv. If your layout differs, adjust.
sample_ids <- basename(dirname(files))
names(files) <- sample_ids

# Important note:
# - If your directory contains both abundance.h5 and abundance.tsv for the same sample,
#   the 'files' vector will contain two paths. tximport's input expects a single file per sample
#   (for type="kallisto" it can accept abundance.h5). To avoid duplicates, you may want to
#   deduplicate by sample name, preferring abundance.h5 where present. Example improvement:
#
#   prefer_h5 <- grepl("abundance.h5$", files)
#   files_df <- data.frame(path = files, sample = sample_ids, is_h5 = prefer_h5, stringsAsFactors = FALSE)
#   files_df <- files_df[order(files_df$sample, -files_df$is_h5), ]
#   files_unique <- tapply(files_df$path, files_df$sample, function(x) x[1])
#   files <- unname(files_unique); names(files) <- names(files_unique)
#
# For brevity we proceed with the simple basename(dirname(...)) approach which often suffices.

# ----------------------------
# Read tx2gene mapping
# ----------------------------
# tx2gene must be a two-column file mapping transcript IDs (TXNAME) to gene IDs (GENEID).
# Typical format created from GTF with a short awk/perl script. tximport expects a data.frame
# with two columns: transcript identifier and gene identifier.
if (!file.exists(tx2gene_file)) {
  stop("tx2gene file not found: ", tx2gene_file, "\nCreate tx2gene as two-column TSV: TXNAME<TAB>GENEID")
}
tx2g <- read.delim(tx2gene_file, header=FALSE, stringsAsFactors=FALSE)
if (ncol(tx2g) < 2) stop("tx2gene file must have at least two columns: tx and gene")
colnames(tx2g) <- c("TXNAME","GENEID")

# NOTE: If your transcript IDs have version suffixes (ENST00000.1), tximport has `ignoreTxVersion=TRUE`
# to remove dot-number suffixes when matching txnames. Using ignoreTxVersion=TRUE below is helpful
# if tx names in tx2gene and kallisto outputs differ by version suffix only.

# ----------------------------
# tximport: import and summarize to genes
# ----------------------------
# countsFromAbundance: choose method carefully. Options:
#  - "no": use kallisto estimated counts (raw)
#  - "scaledTPM": scale TPMs to library size to obtain approximate counts
#  - "lengthScaledTPM": recommended in many workflows — produces counts that adjust for
#     both transcript length and library size, suitable for gene-level differential analysis.
#
# In this script we use "lengthScaledTPM". If you plan to run downstream tests that expect integer
# counts (edgeR, DESeq2), these are accepted; both DESeq2 and edgeR authors describe tximport
# workflows in their vignettes.
txi <- tximport(files, type="kallisto",
                tx2gene=tx2g,
                countsFromAbundance="lengthScaledTPM",
                ignoreTxVersion=TRUE)

# txi is a list with elements: abundance (TPM-like), counts (estimated counts per gene), length (effective length)
# Check components for sanity:
if (!("counts" %in% names(txi))) stop("tximport failed to produce 'counts' element. Inspect txi object.")

# ----------------------------
# edgeR normalization and log-CPM
# ----------------------------
# Create a DGEList from txi counts. Note that tximport's counts are not integer reads but scaled
# estimates — it's still common to treat them as counts for normalization and differential analysis.
dge <- DGEList(counts=txi$counts)

# calcNormFactors performs TMM normalization (Trimmed Mean of M-values), which adjusts for
# compositional differences between libraries. This is standard prior to computing logCPM.
dge <- calcNormFactors(dge, method="TMM")

# cpm(..., log=TRUE, prior.count=1) computes log2-CPM with a small prior count to avoid -Inf.
# The prior is added on counts scale before normalization; choice of prior.count can affect low-count behavior.
logcpm <- cpm(dge, log=TRUE, prior.count=1)

# Convert to data.frame and add gene_id column (rownames)
expr_df <- as.data.frame(logcpm)
expr_df <- cbind(gene_id=rownames(expr_df), expr_df)

# ----------------------------
# Write outputs (TSV)
# ----------------------------
# It's a good idea to write also an RDS or RData snapshot of the txi/dge objects for reproducibility,
# but here we write a TSV matrix for downstream tools (Python/R).
# Note: write.table defaults are used; for very large matrices consider fwrite (data.table) for speed.
out_dir <- dirname(out_expr)
if (!dir.exists(out_dir)) dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)

write.table(expr_df, file=out_expr, sep="\t", row.names=FALSE, quote=FALSE)
# Also write a sample table with mapping to kallisto directories for provenance
sample_tbl <- data.frame(sample_id=sample_ids, kallisto_path=unname(files), stringsAsFactors=FALSE)
sample_tbl_dir <- dirname(out_sample_table)
if (!dir.exists(sample_tbl_dir)) dir.create(sample_tbl_dir, recursive=TRUE, showWarnings=FALSE)
write.table(sample_tbl, file=out_sample_table, sep="\t", row.names=FALSE, quote=FALSE)

cat("Wrote expression matrix to:", out_expr, "\n")
cat("Wrote sample table to:", out_sample_table, "\n")

# ----------------------------
# Post-run recommendations & troubleshooting
# ----------------------------
# 1) Confirm sample ordering:
#    - The columns of out_expr (after gene_id) will be in the order of 'files' discovered
#      using list.files + basename(dirname(files)). If you rely on a particular sample order,
#      build 'files' as a named vector where names(files) are sample IDs and entries are paths,
#      and pass that to tximport. Example:
#        samples <- c("s1"="/path/to/s1/abundance.h5", "s2"="/path/to/s2/abundance.h5")
#        txi <- tximport(samples, type="kallisto", tx2gene=...)
#
# 2) Verify txnames mapping:
#    - If tximport fails to match transcript IDs between the kallisto output and tx2gene,
#      set ignoreTxVersion=TRUE (done above) OR adapt tx2gene to remove version suffixes.
#
# 3) Check countsFromAbundance choice:
#    - Different options lead to different count estimates. If you want to use transcript-level
#      summarize without length correction (e.g., for some multi-mapped cases), use countsFromAbundance="no".
#
# 4) Reproducibility:
#    - Save R session info: writeLines(capture.output(sessionInfo()), "sessionInfo.txt")
#    - Consider storing txi and dge as RDS: saveRDS(txi, file="txi.rds"); saveRDS(dge, file="dge.rds")
#
# 5) Integration with downstream DE analysis:
#    - For gene-level differential expression, use edgeR/limma-voom or DESeq2 on either:
#       a) txi$counts with calcNormFactors (as here), or
#       b) use tximport with countsFromAbundance="lengthScaledTPM" and treat as counts.
#    - See tximport and DESeq2/edgeR vignettes for recommended workflows.
#
# 6) Large datasets:
#    - For many samples, reading all abundance files and creating txi might be slow; ensure sufficient memory.
#    - Consider pre-filtering transcripts/genes with low counts before normalization if necessary.
#
# 7) Error handling:
#    - If tximport throws errors that abundance.h5 can't be read, ensure the HDF5 bindings are installed
#      and that the files are not corrupted. Kallisto's abundance.h5 requires rhdf5 (not used directly here
#      because tximport handles HDF5).
#
# 8) Next steps:
#    - If you'd like me to:
#       * Convert this script to accept CLI args (path arguments, countsFromAbundance choice, output dir),
#       * Add RDS snapshots and sessionInfo saving,
#       * Deduplicate abundance.h5 vs abundance.tsv per sample (prefer HDF5),
#       * Add gene-level filtering and compact expression outputs (TPM, raw counts),
#      tell me which and I will implement it.
#
# End of annotated script.
