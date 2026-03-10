#!/usr/bin/env python3
"""
AF Matrix Plotting and Variant-Gene Correlation Helper — Annotated

What this script does
---------------------
- Loads an allele-fraction (AF) matrix TSV where rows are variants and
  columns after the first five metadata columns (CHROM, POS, ID, REF, ALT)
  are per-sample allele fractions (AF = alt / (ref + alt) or similar).
- Produces per-sample histograms, a pooled AF density, PCA of samples based on
  AF patterns, a heatmap of the most variable variants, and optionally:
  - computes correlations between top-variable variants and genes in an
    expression matrix and writes correlation TSV + a handful of scatterplots.

Important constraints and assumptions
------------------------------------
- AF_FILE must be a tab-separated file with at least these five leading columns:
    CHROM  POS  ID  REF  ALT  <sample1>  <sample2> ...
  If your AF file uses different column names or ordering, adjust the code
  that slices [:5] vs [5:].
- Per-sample AF columns must be numeric or NA / '.' / 'NA' strings convertible to NaN.
- EXPR_MATRIX (optional): expects genes x samples table with gene IDs as row
  index (first column) and sample names as headers. Sample names must match AF
  sample column names exactly. The script attempts to reorder expression
  columns if they are a superset of AF samples.
- The script treats NaN as missing and imputes per-variant mean for PCA only
  (simple and acceptable for exploratory PCA; replace with more careful
  imputation if needed).

Performance and scaling notes
-----------------------------
- For whole-genome AF matrices with tens or hundreds of thousands of variants,
  the heatmap and the variant-gene correlation loop can be slow and memory-heavy.
  Use small top_variant_count for correlation (default 100) and top_n for heatmap.
- PCA uses sklearn.PCA if available (recommended) and otherwise falls back to SVD.
- Expression correlation currently computes pairwise correlations in Python loops
  (variant × gene) on the top-variant subset; for large gene sets use optimized
  vectorized methods or sparse approaches.

Safety & reproducibility
------------------------
- This script is intentionally small and dependency-light: pandas, numpy,
  seaborn, matplotlib, sklearn (optional) are used.
- For reproducibility, record package versions (e.g., pip freeze) or run inside
  a conda environment / container.

How to adapt
------------
- Convert constants AF_FILE / EXPR_MATRIX / OUTPUT_DIR to argparse args for more
  flexible usage in pipelines.
- Add logging to capture progress for long runs (e.g., logging module with INFO/DEBUG).
- Add a --top-variants / --top-genes CLI option to control correlation scope.
"""

# Standard libraries
import os
import sys
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# -----------------------------------------------------------------------------
# EDIT THESE to change source files and output target
# -----------------------------------------------------------------------------
AF_FILE = "/home/jgd/Documents/bioinformatics_working/counts/kallisto_quant_ak-cwb/variants_af.tsv"
# Set EXPR_MATRIX to a real gene expression matrix (genes x samples) if you want
# expression correlation scatter plots. Otherwise leave as placeholder or empty.
EXPR_MATRIX = "/home/jgd/Documents/bioinformatics_working/counts/kallisto_quant_ak-cwb/expr_matrix.tsv"

OUTPUT_DIR = "/home/jgd/Documents/bioinformatics_working/output"
PLOT_DIR = OUTPUT_DIR
os.makedirs(PLOT_DIR, exist_ok=True)

# Output TSV for top correlations (written if expression-based correlations are computed)
CORR_TSV = os.path.join(OUTPUT_DIR, "variant_gene_correlations.tsv")
# -----------------------------------------------------------------------------

# Set global plotting style for consistent figures
sns.set(style="whitegrid")
plt.rcParams.update({"figure.dpi": 150, "font.size": 10})

# -------------------------
# Helper: load AF matrix
# -------------------------
def load_af(af_path):
    """
    Load AF file and split metadata vs numeric AF matrix.

    Expected format:
      CHROM  POS  ID  REF  ALT  <sample1> <sample2> ...

    Returns:
      meta: pandas.DataFrame with first 5 columns (unchanged types)
      af: pandas.DataFrame of floats (variants x samples), same row order as meta
      samples: list of sample column names (order preserved)

    Notes:
     - We coerce AF columns to numeric using pandas.to_numeric(errors='coerce')
       so '.' / 'NA' become NaN.
     - This keeps the code robust to common AF file encodings.
    """
    df = pd.read_csv(af_path, sep="\t", na_values=["NA", "."], dtype=str)
    if df.shape[1] < 6:
        raise ValueError("AF file must have at least 6 columns (CHROM,POS,ID,REF,ALT,<samples...>)")
    # Keep first five columns as metadata (strings preserved)
    meta = df.iloc[:, :5].copy().reset_index(drop=True)
    # Rest of columns are per-sample AF; coerce to numeric
    af = df.iloc[:, 5:].apply(pd.to_numeric, errors="coerce")
    af.index = meta.index  # ensure identical index for row-aligned slicing
    samples = list(af.columns)
    return meta, af, samples

# -------------------------
# Per-sample histogram grid
# -------------------------
def plot_sample_histograms(af, samples, out_png, max_cols=4):
    """
    Small multiples: one histogram + KDE for each sample.

    Arguments:
      af: DataFrame variants x samples (floats with NaN)
      samples: ordered list of sample column names
      out_png: path to save PNG
      max_cols: number of columns per page in the small-multiple grid

    Behavior:
      - For samples with no numeric AF values, place a placeholder.
      - Limits x-axis to [0,1] because AF is a fraction.
      - Uses seaborn.histplot with kde for quick visual inspection.
    """
    n = len(samples)
    cols = min(max_cols, n) or 1
    rows = int(math.ceil(n / cols))
    fig, axes = plt.subplots(rows, cols, figsize=(4*cols, 3*rows), constrained_layout=True)
    # axes could be 2D or 1D depending on rows/cols; flatten to 1D array for iteration
    axes = np.array(axes).reshape(-1)
    i = -1
    for i, s in enumerate(samples):
        ax = axes[i]
        data = af[s].dropna()
        if len(data) == 0:
            ax.text(0.5, 0.5, "no data", ha='center', va='center')
            ax.set_title(s)
            ax.set_xlim(0,1)
            continue
        sns.histplot(data, bins=40, kde=True, stat="density", ax=ax, color="C0")
        ax.set_title(s)
        ax.set_xlim(0,1)
    # hide leftover (empty) axes if any
    for j in range(i+1, len(axes)):
        axes[j].axis('off')
    fig.suptitle("Per-sample AF histograms (AF = alt/(ref+alt))", fontsize=14)
    fig.savefig(out_png, dpi=150)
    plt.close(fig)

# -------------------------
# AF density across all samples
# -------------------------
def plot_af_density_all(af, out_png):
    """
    Flatten AF values across all samples and plot a single aggregated AF density.
    Useful to inspect allele frequency spectrum (AF distribution across cohort).
    """
    all_af = af.values.flatten()
    all_af = all_af[~np.isnan(all_af)]
    plt.figure(figsize=(6,4))
    if len(all_af) == 0:
        plt.text(0.5,0.5,"No AF values",ha='center',va='center')
    else:
        sns.histplot(all_af, bins=100, kde=True, stat="density", color="C1")
        plt.xlim(0,1)
        plt.xlabel("Alternate allele fraction (AF)")
        plt.ylabel("Density")
        plt.title("Allele frequency spectrum (all samples)")
    plt.tight_layout()
    plt.savefig(out_png, dpi=150)
    plt.close()

# -------------------------
# PCA on samples using AF
# -------------------------
def compute_pca(af, n_components=2):
    """
    Compute PCA on the AF matrix.

    Implementation details:
      - Input af is variants x samples. PCA should operate on samples, so we take
        the transpose to shape: samples x variants.
      - Missing AF values (NaN) are imputed with the variant (column) mean before PCA.
        This is a quick imputation suitable for exploratory analysis.
      - Uses sklearn.PCA if available; otherwise uses numpy.linalg.svd fallback.

    Returns:
      pcs: numpy array shape (n_samples, n_components)
      var: explained variance ratio (length n_components)
    """
    # Impute missing values by per-variant mean (axis=0 when samples are columns -> variants columns)
    A = af.fillna(af.mean(axis=0)).values.T  # shape: samples x variants
    try:
        from sklearn.decomposition import PCA
        pca = PCA(n_components=n_components)
        pcs = pca.fit_transform(A)
        var = pca.explained_variance_ratio_
    except Exception:
        # fallback SVD: center columns, get left singular vectors scaled by singular values
        A_centered = A - A.mean(axis=0)
        U, S, Vt = np.linalg.svd(A_centered, full_matrices=False)
        # approximate PCs: left singular vectors times singular values
        pcs = U[:, :n_components] * S[:n_components]
        var = (S**2) / np.sum(S**2)
        var = var[:n_components]
    return pcs, var

def plot_pca(af, samples, out_png):
    """
    Scatterplot of PC1 vs PC2 with sample labels.

    Tips:
     - For many samples label overlaps may clutter the plot; consider saving interactive plots
       or plotting labels selectively (top PCs extremes).
    """
    pcs, var = compute_pca(af, n_components=2)
    plt.figure(figsize=(6,5))
    x = pcs[:,0]; y = pcs[:,1]
    sns.scatterplot(x=x, y=y)
    for i, s in enumerate(samples):
        plt.text(x[i], y[i], s, fontsize=8, alpha=0.8)
    plt.xlabel(f"PC1 ({var[0]*100:.1f}%)" if len(var)>0 else "PC1")
    plt.ylabel(f"PC2 ({var[1]*100:.1f}%)" if len(var)>1 else "PC2")
    plt.title("PCA of samples (based on AF matrix)")
    plt.tight_layout()
    plt.savefig(out_png, dpi=150)
    plt.close()

# -------------------------
# Heatmap of top-variable variants
# -------------------------
def plot_top_variants_heatmap(af, samples, out_png, top_n=200):
    """
    Select top_n variants by variance across samples and plot a heatmap.

    Visualization choices:
      - z-score each variant (row-wise) to highlight relative differences across samples.
      - For variants with zero variance, std will be zero -> NaN; replaced with 0 for display.
      - For large top_n, the plot becomes tall; figsize is scaled heuristically.
    """
    var_series = af.var(axis=1, skipna=True)
    top_idx = var_series.sort_values(ascending=False).head(top_n).index
    sub = af.loc[top_idx]
    # z-score by variant (row-wise): (x - mean)/std
    sub_z = sub.subtract(sub.mean(axis=1), axis=0).divide(sub.std(axis=1).replace(0, np.nan), axis=0)
    plt.figure(figsize=(max(8, len(samples)*0.2), max(6, top_n*0.03)))
    # fill NaN with zero for visualization
    sns.heatmap(sub_z.fillna(0), cmap="vlag", center=0, xticklabels=True, yticklabels=False)
    plt.xlabel("Samples")
    plt.title(f"Heatmap of top {top_n} most variable variants (z-scored by variant)")
    plt.tight_layout()
    plt.savefig(out_png, dpi=150)
    plt.close()

# -------------------------
# Expression-based variant -> gene correlation and scatterplots
# -------------------------
def optionally_plot_expr_scatter(af_meta, af, samples, expr_path, out_dir, max_plots=6, top_variant_count=100):
    """
    For a small set of top-variable variants, compute Pearson correlations with genes
    in the expression matrix and produce:
      - variant_gene_correlations.tsv (sorted by absolute correlation)
      - up to max_plots scatter plots of AF vs expression with regression line.

    Implementation details & caveats:
      - expr_path must be a genes x samples TSV with geneID in first column (index_col=0).
      - Samples in expression matrix must match AF samples. The function attempts to
        reorder expression columns to AF sample order when expr contains a superset.
      - Correlation is computed only on pairwise non-missing sample subsets.
      - This is a nested loop (variants × genes) limited to top_variant_count for speed.
      - For larger searches, consider using matrix algebra or optimized correlation libraries.
    """
    if not os.path.exists(expr_path):
        print("Expression matrix not found at", expr_path, "- skipping expression-based plots", file=sys.stderr)
        return

    # Load expression; genes as row index
    expr = pd.read_csv(expr_path, sep="\t", index_col=0)
    expr_samples = list(expr.columns)
    # Reorder expression columns to match AF samples if possible
    if samples != expr_samples:
        if set(samples) <= set(expr_samples):
            expr = expr[samples]
            expr_samples = list(expr.columns)
            print("Reordered expression matrix to AF sample order.", file=sys.stderr)
        else:
            print("Sample names in expression matrix do not match AF samples; skipping expr plots", file=sys.stderr)
            return

    # Convert to numpy arrays for speed: expr_vals: genes x samples, af_vals: variants x samples
    expr_vals = expr.values
    expr_index = expr.index.tolist()
    af_vals = af.values
    # Compute per-variant variance and choose top indices to test
    var_var = np.nanvar(af_vals, axis=1)
    if top_variant_count <= 0:
        print("top_variant_count must be > 0; skipping expr plots", file=sys.stderr)
        return
    top_variant_idx = np.argsort(var_var)[-top_variant_count:]
    correlations = []

    # For each selected variant compute correlation with every gene (vectorized alternatives possible)
    for vi in top_variant_idx:
        v = af_vals[vi, :]
        # skip variants with almost all missing values (needs at least a few points to correlate)
        if np.isnan(v).sum() > (len(v) - 3):
            continue
        for gi in range(expr_vals.shape[0]):
            g = expr_vals[gi, :]
            mask = (~np.isnan(v)) & (~np.isnan(g))
            if mask.sum() < 4:
                continue
            r = np.corrcoef(v[mask], g[mask])[0,1]
            if np.isnan(r):
                continue
            correlations.append((abs(r), float(r), int(vi), int(gi)))

    if not correlations:
        print("No valid correlations found; skipping expr scatterplots", file=sys.stderr)
        return

    # Sort by absolute correlation descending
    correlations.sort(reverse=True)
    # Write correlation TSV for inspection (cap size for safety)
    max_write = min(len(correlations), 100000)
    out_tsv = os.path.join(out_dir, "variant_gene_correlations.tsv")
    with open(out_tsv, 'w') as oh:
        oh.write("abs_r\tpearson_r\tvariant_index\tgene_index\tvariant_label\tgene\n")
        for idx, (abs_r, r, vi, gi) in enumerate(correlations[:max_write]):
            try:
                variant_row = af_meta.iloc[vi]
                var_label = f"{variant_row['CHROM']}:{variant_row['POS']}:{variant_row['ID']}"
            except Exception:
                var_label = f"var_idx_{vi}"
            gene = expr.index[gi]
            oh.write(f"{abs_r:.6f}\t{r:.6f}\t{vi}\t{gi}\t{var_label}\t{gene}\n")
    print(f"Wrote correlation TSV to {out_tsv}", file=sys.stderr)

    # Make scatterplots for top correlated pairs (up to max_plots)
    made = 0
    for abs_r, r, vi, gi in correlations[:max_plots]:
        try:
            variant_row = af_meta.iloc[vi]
            var_label = f"{variant_row['CHROM']}:{variant_row['POS']}:{variant_row['ID']}"
            v_af = pd.Series(af_vals[vi, :], index=samples)
        except Exception:
            var_label = f"var_idx_{vi}"
            v_af = pd.Series(af_vals[vi, :], index=samples)
        gene = expr.index[gi]
        g_expr = expr.loc[gene]
        # align and drop nan pairs
        df = pd.concat([v_af, g_expr], axis=1, join='inner').dropna()
        if df.shape[0] < 4:
            continue
        plt.figure(figsize=(5,4))
        sns.regplot(x=df.iloc[:,0], y=df.iloc[:,1], scatter_kws={'s':20}, line_kws={'color':'red'})
        plt.xlabel("AF (alt/(ref+alt))")
        plt.ylabel("Expression")
        plt.title(f"{var_label} vs {gene}\nr={r:.3f}")
        out_png = os.path.join(out_dir, f"top_variant_gene_scatter_{made+1}.png")
        plt.tight_layout()
        plt.savefig(out_png, dpi=150)
        plt.close()
        made += 1
    print(f"Wrote {made} expression scatterplots to {out_dir}", file=sys.stderr)

# -------------------------
# Main orchestration
# -------------------------
def main():
    # Validate AF_FILE exists
    if not os.path.exists(AF_FILE):
        print("AF file not found:", AF_FILE, file=sys.stderr)
        sys.exit(1)
    print("Loading AF matrix from", AF_FILE, file=sys.stderr)
    meta, af, samples = load_af(AF_FILE)
    print(f"AF matrix: {af.shape[0]} variants x {af.shape[1]} samples", file=sys.stderr)

    # Basic exploratory plots written into OUTPUT_DIR
    plot_sample_histograms(af, samples, os.path.join(PLOT_DIR, "af_sample_histograms.png"))
    plot_af_density_all(af, os.path.join(PLOT_DIR, "af_density_all_samples.png"))
    plot_pca(af, samples, os.path.join(PLOT_DIR, "af_pca.png"))
    plot_top_variants_heatmap(af, samples, os.path.join(PLOT_DIR, "af_top_variants_heatmap.png"), top_n=200)

    # Optional expression-based correlation & scatterplots
    if EXPR_MATRIX and EXPR_MATRIX != "path/to/expr_matrix.tsv":
        # Keep the top-variant subset small (100) to keep runtime reasonable
        optionally_plot_expr_scatter(meta, af, samples, EXPR_MATRIX, PLOT_DIR, max_plots=6, top_variant_count=100)
    else:
        print("EXPR_MATRIX not set (left as placeholder) — skipping expression-based plots", file=sys.stderr)

    print("All plots and optional correlation outputs written to", OUTPUT_DIR, file=sys.stderr)

if __name__ == "__main__":
    main()
