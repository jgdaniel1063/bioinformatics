#!/usr/bin/env python3
"""
results_plots_from_tsv_full_annotated.py

Heavily annotated version of results_plots_from_tsv_full.py.

Purpose
-------
Read a tab-separated summary table (results_summary.tsv produced by analysis scripts)
and produce per-sample and aggregate diagnostic plots. Per-sample output includes a
set of PNG pages and a combined multi-page PDF. Aggregate plots include mutation
spectrum, genotype stacks, heatmap of metrics, MA / volcano / p-value diagnostics.

This annotated file explains:
 - expected input shapes and column names
 - fallback behaviors and heuristics used to detect LFC/p-value columns
 - plotting choices and parameterization
 - common failure modes and how to troubleshoot them
 - suggestions for further improvements (CLI, caching, styling, reproducibility)

Notes / Assumptions
-------------------
 - The input TSV is expected to be a "wide" table where each row is a sample and
   columns contain scalar metrics (e.g., mean_mapq), lists/dicts JSON-encoded
   (e.g., mut_six, mut_context96, vaf_list), and file paths (bam, counts_path).
 - The script tolerates many column name variants by heuristics (detect_lfc_col,
   detect_pval_col, detect_basemean_col). If those heuristics fail, aggregate
   diagnostics that rely on those columns are skipped silently.
 - Plotting requires seaborn / pandas / numpy / matplotlib. If missing, the
   script will raise on import (explicit error). Per-sample plotting is robust to
   missing specific metrics.
 - Filenames and output locations are hard-coded near the top. Convert to CLI
   with argparse if you want flexibility.
 - For best reproducibility, pin package versions or run inside a container.

Outputs
-------
Per sample (saved into SAMPLES_DIR):
 - <sample>_metrics.png         : horizontal bar chart of numeric metrics
 - <sample>_vaf.png             : VAF histogram / KDE or a single mean bar
 - <sample>_mutational_spectrum.png : 6-class mutation bar chart (derived from mut_six or mut_context96)
 - <sample>_genotype.png        : stacked genotype counts
 - <sample>_summary.pdf         : multi-page PDF combining the four PNG pages (one per page)

Aggregate (saved into OUTPUT_DIR):
 - mutational_spectrum_6class.png
 - genotype_counts_stacked.png
 - metrics_heatmap_zscore.png
 - ma_plot.png (if LFC + baseMean detected)
 - volcano_plot.png (if LFC + pvalue/padj detected)
 - pvalue_histogram.png, qq_plot.png

Also copies per-sample PNGs into RESULTS_PARENT/<sample>/plots/ if toggled, and writes a manifest.

Large-scale notes
-----------------
 - For many samples (>200), per-sample PDF generation can be slow; consider parallelizing.
 - Heatmap and bar widths scale with number of samples; plots attempt to compute a readable size,
   but manual tuning may be needed for very large cohort sizes.
 - The script tries to be forgiving with input column types (strings, JSON-encoded lists/dicts).
   It uses parse_possible_json to convert string representations into Python objects.
"""

from pathlib import Path
import json
import ast
import re
import sys
import math
import shutil
import traceback
import textwrap
import argparse

# Force non-interactive matplotlib backend suitable for headless servers / batch runs
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

# ---------------- HARD-CODED PATHS / TOGGLES ----------------
# Edit these to point to your input TSV and desired output locations.
# Consider converting these into argparse arguments for scripted runs.
INPUT_TSV = Path("/home/jgd/Documents/bioinformatics_working/results/enu_results/results_summary.tsv")
OUTPUT_DIR = Path("/home/jgd/Documents/bioinformatics_working/output/summary_results/plots")
SAMPLES_DIR = OUTPUT_DIR / "samples"
RESULTS_PARENT = Path("/home/jgd/Documents/bioinformatics_working/results/enu_results/results_summary")
RESULTS_SUBDIRS = ["plots", "qc", "logs", "raw"]

# Visual parameters
DPI = 300
GROUP_REGEX = r"\d+$"          # used by sample_to_group: strip trailing digits to group samples into families
MAX_SAMPLE_METRICS = 200       # if a sample has many numeric metrics, only top N displayed by magnitude
HEATMAP_TOP_N = 50             # pick top-N variable metrics for the heatmap

# Toggles: set False to skip computationally / visually heavy outputs
INCLUDE_AGGREGATE_DIAGNOSTICS = True   # include MA/volcano/pvalue hist/qq and heatmap
COPY_PAGES_TO_RESULTS_PARENT = True    # copy per-sample PNGs to RESULTS_PARENT/<sample>/plots/
# -----------------------------------------------------------

# imports & checks for required plotting/data libraries
try:
    import pandas as pd
    import numpy as np
    import seaborn as sns
except Exception as e:
    # Fail fast: these libs are required for the script's plotting features
    print("ERROR: Required packages missing: pandas, numpy, seaborn, matplotlib", file=sys.stderr)
    raise

# Global seaborn/matplotlib styling tuned for readable figures
sns.set(style="whitegrid")
plt.rcParams.update({
    "font.size": 12,
    "axes.titlesize": 14,
    "axes.labelsize": 12,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "legend.fontsize": 10,
    "figure.titlesize": 15
})

# ---------------- Utilities ----------------
def parse_possible_json(s):
    """
    Robustly parse a field that may be:
      - already a dict/list
      - a JSON string (double quotes)
      - a single-quoted Python string representation of a list/dict (e.g. "{'A':1}")
      - an ast literal (e.g., Python repr)
      - a plain string (returned as-is)
    Returns: parsed object or the original string fallback
    Usage: pass strings from TSV columns like 'mut_six', 'vaf_list', 'mut_context96'
    """
    if s is None:
        return None
    if isinstance(s, (dict, list)):
        return s
    s = str(s).strip()
    if s == "":
        return None
    # Try standard JSON
    try:
        return json.loads(s)
    except Exception:
        pass
    # If it looks like JSON but with single quotes, try replacing quotes
    if (s.startswith("{") or s.startswith("[")) and ("'" in s):
        try:
            return json.loads(s.replace("'", '"'))
        except Exception:
            pass
    # Try Python literal evaluation
    try:
        return ast.literal_eval(s)
    except Exception:
        pass
    # Give up — return raw string
    return s

def sanitize_fname(s):
    """
    Create a filesystem-safe filename fragment from a sample name.
    - Keeps alphanumeric, dot, dash, underscore.
    - Truncates to 220 characters for safety.
    """
    return re.sub(r"[^A-Za-z0-9._-]+", "_", str(s))[:220]

def sample_to_group(sample):
    """
    Heuristic grouping function for family/group-level aggregation.
    - Based on GROUP_REGEX: strips trailing digits.
    - Example: 'bq01nothrom1' -> 'bq01nothrom'
    Useful to compute family-level summary stats.
    """
    return re.sub(GROUP_REGEX, "", str(sample))

def numeric_cols_from_df(df, sample_col):
    """
    Identify numeric columns in the dataframe by coercing each column to numeric
    and counting non-NA values. Excludes known non-numeric columns in 'exclude'.
    Returns list of column names that contain numeric values.
    """
    exclude = {sample_col, "bam", "raw", "vaf_list", "mut_six", "mut_context96", "counts_path"}
    nums = []
    for c in df.columns:
        if c in exclude:
            continue
        coerced = pd.to_numeric(df[c], errors="coerce")
        # require at least one non-NA numeric to consider this a numeric column
        if coerced.notna().sum() > 0:
            nums.append(c)
    return nums

def bar_with_palette(ax, keys, vals, palette="Set2", rotation=0):
    """
    Helper to draw a categorical bar chart with seaborn palette.
    - Returns the Axes for further annotation.
    - keys: list of labels
    - vals: list of numeric values (same length)
    """
    colors = sns.color_palette(palette, n_colors=len(keys))
    idx = list(range(len(keys)))
    ax.bar(idx, vals, color=colors)
    ax.set_xticks(idx)
    ax.set_xticklabels(keys, rotation=rotation, ha="center")
    return ax

# Helpers to detect common LFC / p-value / baseMean columns despite variation in naming
def detect_lfc_col(df):
    """
    Heuristically find a log2 fold-change column in df.
    Checks a list of common names first, then falls back to presence of 'lfc' or 'log2' in column names.
    Returns column name or None if not found.
    """
    candidates = ["log2FoldChange","log2fold","log2fold_treatment_control","log2Fold"]
    for c in candidates:
        if c in df.columns:
            return c
    # fallback: look for 'lfc' or 'log2' in names
    for c in df.columns:
        if "lfc" in c.lower() or "log2" in c.lower():
            return c
    return None

def detect_pval_col(df):
    """
    Find a p-value-like column in df; prefers 'padj' (adjusted) then 'pvalue' etc.
    Returns column name or None.
    """
    # prefer padj then pvalue
    if "padj" in df.columns:
        return "padj"
    if "adjusted_pvalue" in df.columns:
        return "adjusted_pvalue"
    if "pvalue" in df.columns:
        return "pvalue"
    # fallback: any column starting with 'p' and numeric
    for c in df.columns:
        if c.lower().startswith("p") and pd.to_numeric(df[c], errors="coerce").notna().sum()>0:
            return c
    return None

def detect_basemean_col(df):
    """
    Find a baseMean / mean column to use on the A-axis for an MA-plot.
    Returns column name or None.
    """
    candidates = ["baseMean","mean"]
    for c in candidates:
        if c in df.columns:
            return c
    # fallback: any col with 'mean' or 'avg' and numeric
    for c in df.columns:
        if "mean" in c.lower() or "avg" in c.lower():
            if pd.to_numeric(df[c], errors="coerce").notna().sum()>0:
                return c
    return None

# ---------------- Aggregate diagnostic plots ----------------
def aggregate_mutational_spectrum(records, outdir: Path):
    """
    Aggregate per-sample 6-class mutation spectra (mut_six) and produce a barplot.
    - Records: list of dicts (rows), each may contain 'mut_six' as dict or JSON string.
    - If mut_six is missing but mut_context96 exists, consider deriving 6-class from context.
    Returns path to saved image or None if no data.
    """
    outdir.mkdir(parents=True, exist_ok=True)
    six_all = {}
    for rec in records:
        six = rec.get("mut_six")
        if isinstance(six, str):
            six = parse_possible_json(six) or {}
        if isinstance(six, dict):
            for k,v in six.items():
                try:
                    six_all[k] = six_all.get(k,0) + int(v)
                except Exception:
                    continue
    if not six_all:
        return None
    keys = ["C>A","C>G","C>T","T>A","T>C","T>G"]
    vals = [six_all.get(k,0) for k in keys]
    fig, ax = plt.subplots(figsize=(8,5), constrained_layout=True)
    bar_with_palette(ax, keys, vals, palette="Set2", rotation=0)
    ax.set_ylabel("Count")
    ax.set_title("Aggregate 6-class mutational spectrum")
    out = outdir / "mutational_spectrum_6class.png"
    fig.savefig(out, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    return out

def aggregate_genotype_stacked(records, outdir: Path):
    """
    Produce a stacked bar plot showing genotype counts (hom_ref / het / hom_alt / missing) per sample.
    - Useful to spot samples with many missing genotypes or unusual genotype distributions.
    - Colors are chosen to be distinguishable; for colorblind-friendly choices adjust palette.
    """
    outdir.mkdir(parents=True, exist_ok=True)
    rows = []
    for rec in records:
        sample = rec.get("sample")
        if sample is None:
            continue
        rows.append({
            "sample": sample,
            "hom_ref": int(rec.get("gt_hom_ref") or 0),
            "het": int(rec.get("gt_het") or 0),
            "hom_alt": int(rec.get("gt_hom_alt") or 0),
            "missing": int(rec.get("gt_missing") or 0)
        })
    if not rows:
        return None
    gtdf = pd.DataFrame(rows).set_index("sample")
    fig, ax = plt.subplots(figsize=(max(10, 0.3*len(gtdf)), 6), constrained_layout=True)
    ind = range(len(gtdf))
    bottom = np.zeros(len(gtdf))
    colors = ["#4C72B0","#55A868","#C44E52","#8172B2"]
    for i,col in enumerate(["hom_ref","het","hom_alt","missing"]):
        vals = gtdf[col].values
        ax.bar(ind, vals, bottom=bottom, label=col, color=colors[i])
        bottom = bottom + vals
    ax.set_xticks(ind)
    ax.set_xticklabels(gtdf.index.tolist(), rotation=90)
    ax.set_ylabel("Count")
    ax.set_title("Genotype counts per sample (stacked)")
    ax.legend()
    out = outdir / "genotype_counts_stacked.png"
    fig.savefig(out, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    return out

def aggregate_heatmap(df, sample_col, outdir: Path, top_n=HEATMAP_TOP_N):
    """
    Build a heatmap of z-score normalized top-N variable numeric metrics.
    - Identifies numeric columns via numeric_cols_from_df, filters columns with any non-NA values.
    - Computes variance across samples and picks the top-N most variable metrics.
    - Z-scores (mean=0, sd=1) computed column-wise to allow comparison across metrics.
    """
    outdir.mkdir(parents=True, exist_ok=True)
    metrics = numeric_cols_from_df(df, sample_col)
    if not metrics:
        return None
    mat = df.set_index(sample_col)[metrics].apply(pd.to_numeric, errors="coerce")
    mat = mat.loc[:, mat.notna().any(axis=0)]
    mat = mat[~mat.index.duplicated(keep="first")]  # drop duplicate sample names if any
    if mat.shape[1] == 0 or mat.shape[0] == 0:
        return None
    variances = mat.var(axis=0, skipna=True).sort_values(ascending=False)
    keep = variances.index[:min(top_n, len(variances))]
    mat2 = mat[keep]
    # Z-score normalization: (x - mean) / std ; use ddof=0 for population std
    matz = (mat2 - mat2.mean()) / mat2.std(ddof=0)
    fig, ax = plt.subplots(figsize=(max(8, 0.35 * matz.shape[0]), max(6, 0.25 * matz.shape[1])), constrained_layout=True)
    sns.heatmap(matz, cmap="vlag", center=0, xticklabels=True, yticklabels=True, ax=ax)
    ax.set_xlabel("metric")
    ax.set_ylabel("sample")
    ax.set_title("Z-score normalized metrics (heatmap)")
    out = outdir / "metrics_heatmap_zscore.png"
    fig.savefig(out, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    return out

def aggregate_ma_plot(df, outdir: Path):
    """
    Create an MA plot if a log fold-change column and a baseMean-like column exist.
    - Uses detect_lfc_col and detect_basemean_col heuristics to find suitable columns.
    - A-axis plotted as log10(baseMean + 1) to compress dynamic range.
    """
    outdir.mkdir(parents=True, exist_ok=True)
    lfc = detect_lfc_col(df)
    base = detect_basemean_col(df)
    if not lfc:
        return None
    # compute A (mean/expression) value
    if base and base in df.columns:
        df2 = df.copy()
        df2["_A"] = pd.to_numeric(df2[base], errors="coerce")
    else:
        # attempt to derive a mean from numeric columns - fallback and may be noisy
        df2 = df.copy()
        numeric_cols = [c for c in df.columns if pd.to_numeric(df[c], errors="coerce").notna().sum()>0]
        if numeric_cols:
            df2["_A"] = df2[numeric_cols].apply(pd.to_numeric, errors="coerce").mean(axis=1)
        else:
            return None
    df2["_M"] = pd.to_numeric(df2[lfc], errors="coerce")
    df2 = df2.dropna(subset=["_A","_M"])
    if df2.shape[0] == 0:
        return None
    fig, ax = plt.subplots(figsize=(8,6), constrained_layout=True)
    ax.scatter(np.log10(df2["_A"].values + 1.0), df2["_M"].values, s=10, alpha=0.6)
    ax.axhline(0, color="grey", linestyle="--")
    ax.set_xlabel(f"log10({base or 'mean'} + 1)")
    ax.set_ylabel(lfc)
    ax.set_title("MA plot")
    out = outdir / "ma_plot.png"
    fig.savefig(out, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    return out

def aggregate_volcano(df, outdir: Path):
    """
    Create a basic volcano plot if LFC and p-value columns exist.
    - Uses detect_lfc_col and detect_pval_col heuristics.
    - Negative log10 of p-values used; p-values clipped to avoid -log10(0).
    """
    outdir.mkdir(parents=True, exist_ok=True)
    lfc = detect_lfc_col(df)
    pcol = detect_pval_col(df)
    if not lfc or not pcol:
        return None
    df2 = df.copy()
    df2["_LFC"] = pd.to_numeric(df2[lfc], errors="coerce")
    df2["_P"] = pd.to_numeric(df2[pcol], errors="coerce")
    df2 = df2.dropna(subset=["_LFC","_P"])
    if df2.shape[0] == 0:
        return None
    fig, ax = plt.subplots(figsize=(8,6), constrained_layout=True)
    ax.scatter(df2["_LFC"], -np.log10(df2["_P"].clip(lower=1e-300)), s=8, alpha=0.6)
    ax.set_xlabel(lfc)
    ax.set_ylabel(f"-log10({pcol})")
    ax.set_title("Volcano plot")
    out = outdir / "volcano_plot.png"
    fig.savefig(out, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    return out

def aggregate_pvalue_hist_qq(df, outdir: Path):
    """
    Produce a p-value histogram and QQ plot. Returns tuple(paths) for histogram and QQ.
    - Uses detect_pval_col to choose column.
    - QQ plot implemented via -log10 sorted p-values vs sorted uniform expectations.
    """
    outdir.mkdir(parents=True, exist_ok=True)
    pcol = detect_pval_col(df)
    if not pcol:
        return (None, None)
    df2 = df.copy()
    df2["_P"] = pd.to_numeric(df2[pcol], errors="coerce")
    df2 = df2.dropna(subset=["_P"])
    if df2.shape[0] == 0:
        return (None, None)
    # histogram
    fig1, ax1 = plt.subplots(figsize=(8,5), constrained_layout=True)
    ax1.hist(df2["_P"], bins=50, color="C0")
    ax1.set_xlabel("P-value")
    ax1.set_ylabel("Frequency")
    ax1.set_title("P-value distribution")
    out_hist = outdir / "pvalue_histogram.png"
    fig1.savefig(out_hist, dpi=DPI, bbox_inches="tight")
    plt.close(fig1)
    # QQ plot: expected vs observed -log10(p)
    fig2, ax2 = plt.subplots(figsize=(6,6), constrained_layout=True)
    # Expected uniform points: use sorted quantiles
    n = df2.shape[0]
    exp = -np.log10(np.linspace(1/(n+1), n/(n+1), n))
    obs = -np.log10(np.sort(df2["_P"].clip(lower=1e-300)))
    mn = min(len(exp), len(obs))
    ax2.scatter(exp[:mn], obs[:mn], s=6)
    lim = max(exp.max(), obs.max())
    ax2.plot([0,lim],[0,lim], color="red", linestyle="--")
    ax2.set_xlabel("Expected -log10(p)")
    ax2.set_ylabel("Observed -log10(p)")
    ax2.set_title("QQ plot")
    out_qq = outdir / "qq_plot.png"
    fig2.savefig(out_qq, dpi=DPI, bbox_inches="tight")
    plt.close(fig2)
    return (out_hist, out_qq)

# ---------------- Per-sample single-page plot builders ----------------
def build_metrics_fig(sample, rec, numeric_cols, max_items=MAX_SAMPLE_METRICS):
    """
    Build a horizontal bar-chart figure summarizing numeric metrics for a single sample.
    - numeric_cols: list of numeric column names in the dataframe
    - rec: dict representing row from the TSV (values may be strings and must be coerced)
    - Items trimmed to top-N by magnitude to keep figure legible
    """
    metrics = {}
    for c in numeric_cols:
        if c in rec and rec.get(c) not in (None, ""):
            try:
                fv = float(rec.get(c))
            except Exception:
                continue
            if not math.isnan(fv):
                metrics[c] = fv
    # If too many metrics, keep the most extreme by absolute value (focus on what differs)
    if len(metrics) > max_items:
        metrics = dict(sorted(metrics.items(), key=lambda kv: abs(kv[1] if kv[1] is not None else 0), reverse=True)[:max_items])
    items = sorted(metrics.items(), key=lambda kv: kv[1] if kv[1] is not None else 0)
    names = [k for k,_ in items]
    vals = [v for _,v in items]
    wrapped = [textwrap.fill(n, width=40) for n in names]
    height = max(6, 0.15 * max(1, len(names)))
    fig, ax = plt.subplots(figsize=(11, height), constrained_layout=True)
    if items:
        y = list(range(len(names)))
        ax.barh(y, vals, color="C3", height=0.6)
        ax.set_yticks(y)
        ax.set_yticklabels(wrapped, fontsize=10)
        ax.invert_yaxis()
        ax.set_xlabel("value")
        ax.set_title(f"{sample}: numeric metrics", fontsize=14)
    else:
        # Fallback page for samples with no numeric metrics
        ax.text(0.5,0.5,"No numeric metrics available", ha="center", va="center")
        ax.axis("off")
    fig.subplots_adjust(left=0.22, right=0.98, top=0.92, bottom=0.05)
    return fig

def build_vaf_fig(sample, rec):
    """
    Build VAF histogram (KDE) page for a sample, or a single mean bar if only mean provided.
    - rec['vaf_list'] may be list-like or JSON-string. parse_possible_json handles both.
    """
    vlist = rec.get("vaf_list") or []
    if isinstance(vlist, str):
        vlist = parse_possible_json(vlist) or []
    fig, ax = plt.subplots(figsize=(11,6), constrained_layout=True)
    if vlist:
        try:
            vnums = [float(x) for x in vlist]
            sns.histplot(vnums, bins=40, kde=True, stat="density", ax=ax, color="C0")
            ax.set_xlim(0,1)
            ax.set_title(f"{sample}: VAF histogram + KDE")
            ax.set_xlabel("VAF")
            ax.set_ylabel("Density")
        except Exception:
            ax.text(0.5,0.5,"Bad VAF list", ha="center", va="center")
            ax.axis("off")
    elif rec.get("vaf_mean") not in (None,""):
        try:
            val = float(rec.get("vaf_mean"))
            ax.bar([0],[val], width=0.4, color="C0")
            ax.set_xticks([])
            ax.set_ylim(0,1)
            ax.set_title(f"{sample}: VAF mean: {val:.3f}")
        except Exception:
            ax.text(0.5,0.5,"No VAF data", ha="center", va="center")
            ax.axis("off")
    else:
        ax.text(0.5,0.5,"No VAF data", ha="center", va="center")
        ax.axis("off")
    fig.subplots_adjust(left=0.08, right=0.98, top=0.92, bottom=0.08)
    return fig

def build_mutation_fig(sample, rec):
    """
    Build 6-class mutational spectrum figure for a single sample:
    - preferrs rec['mut_six'] dict but will attempt to derive it from rec['mut_context96'] if necessary.
    - Shows a simple bar chart with the canonical six substitutions (pyrimidine-centered).
    """
    six = rec.get("mut_six")
    if isinstance(six, str):
        six = parse_possible_json(six) or {}
    # If six empty, try to derive from context dictionary (mut_context96)
    if not (isinstance(six, dict) and any(int(v or 0) for v in six.values())):
        ctx = rec.get("mut_context96")
        if isinstance(ctx, str):
            ctx = parse_possible_json(ctx) or {}
        derived = {}
        if isinstance(ctx, dict) and ctx:
            for k,v in ctx.items():
                try:
                    cnt = int(v)
                except Exception:
                    try: cnt = int(float(v))
                    except: cnt = 0
                # extract substitution from context key like A[C>T]G
                m = re.search(r"\[([ACGT]>[ACGT])\]", str(k))
                if m:
                    sub = m.group(1)
                else:
                    m2 = re.search(r"([ACGT]>[ACGT])", str(k))
                    sub = m2.group(1) if m2 else None
                if sub:
                    derived[sub] = derived.get(sub,0) + cnt
            six = derived
    fig, ax = plt.subplots(figsize=(11,6), constrained_layout=True)
    if six:
        keys = ["C>A","C>G","C>T","T>A","T>C","T>G"]
        vals = [int(six.get(k,0) or 0) for k in keys]
        bar_with_palette(ax, keys, vals, palette="Set2", rotation=0)
        ax.set_ylabel("Count")
        ax.set_title(f"{sample}: 6-class mutational spectrum")
        for lbl in ax.get_xticklabels():
            lbl.set_fontsize(11)
        fig.subplots_adjust(left=0.08, right=0.98, top=0.90, bottom=0.12)
    else:
        ax.text(0.5,0.5,"No mutational spectrum data", ha="center", va="center")
        ax.axis("off")
    return fig

def build_genotype_fig(sample, rec):
    """
    Build a single-row stacked genotype bar for a sample using gt_* columns.
    - If all counts are zero, shows a placeholder message.
    """
    gt_keys = [("gt_hom_ref","hom_ref"), ("gt_het","het"), ("gt_hom_alt","hom_alt"), ("gt_missing","missing")]
    vals = []
    labels = []
    for k,label in gt_keys:
        v = rec.get(k, None)
        try:
            v = int(v) if v not in (None,"") else 0
        except Exception:
            v = 0
        vals.append(v)
        labels.append(label)
    fig, ax = plt.subplots(figsize=(11,4), constrained_layout=True)
    if any(v != 0 for v in vals):
        left = 0
        colors = ["#4C72B0","#55A868","#C44E52","#8172B2"]
        for i,val in enumerate(vals):
            ax.barh([0], [val], left=left, color=colors[i], height=0.5, label=labels[i])
            left += val
        ax.set_yticks([])
        ax.set_xlabel("Count")
        ax.set_title(f"{sample}: Genotype counts (stacked)")
        ax.legend(loc='center left', bbox_to_anchor=(1.02, 0.5))
        fig.subplots_adjust(left=0.06, right=0.78, top=0.88, bottom=0.15)
    else:
        ax.text(0.5,0.5,"No genotype counts", ha="center", va="center")
        ax.axis("off")
    return fig

# ---------------- Copy helper ----------------
def copy_sample_pngs(samples_dir: Path, results_parent: Path):
    """
    Copy per-sample PNGs from samples_dir into results_parent/<sample>/plots/
    - Tries to infer sample name from file suffixes like _metrics.png etc.
    - Avoids overwriting by adding numeric suffixes if needed.
    - Writes and returns list of (src, dest) tuples for manifest writing.
    """
    moved = []
    if not samples_dir.exists():
        return moved
    for p in sorted(samples_dir.glob("*_*.png")):
        if not p.is_file():
            continue
        name = p.name
        # infer sample name from standard suffixes; fall back to splitting on last '_'
        for suf in ["_metrics.png","_vaf.png","_mutational_spectrum.png","_genotype.png","_summary.png"]:
            if name.endswith(suf):
                sample = name[:-len(suf)]
                break
        else:
            sample = name.rsplit("_",1)[0]
        dest_dir = results_parent / sample / "plots"
        dest_dir.mkdir(parents=True, exist_ok=True)
        dest = dest_dir / p.name
        if dest.exists():
            i = 2
            while True:
                newname = f"{sample}_{i}_{p.name}"
                newdest = dest_dir / newname
                if not newdest.exists():
                    dest = newdest
                    break
                i += 1
        try:
            shutil.copy2(str(p), str(dest))
            moved.append((p, dest))
            print(f"COPIED: {p} -> {dest}")
        except Exception:
            traceback.print_exc()
    return moved

# ---------------- Main ----------------
def main():
    """
    Main orchestration:
    - read TSV
    - parse JSON-like columns
    - compute aggregate diagnostics (optional)
    - create per-sample PNGs and PDFs
    - optionally copy pages into result folders and write manifest
    """
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    SAMPLES_DIR.mkdir(parents=True, exist_ok=True)

    if not INPUT_TSV.exists():
        print("Input TSV not found:", INPUT_TSV, file=sys.stderr)
        sys.exit(1)

    # Read TSV into DataFrame — keep as strings then parse complex columns below
    df = pd.read_csv(INPUT_TSV, sep="\t", dtype=str).fillna("")
    records = []
    # Convert each row to a dict and parse mut_six / mut_context96 / vaf_list into Python objects
    for _, row in df.iterrows():
        rec = row.to_dict()
        for col in ("mut_six","mut_context96","vaf_list"):
            if col in rec and rec[col] not in (None,""):
                rec[col] = parse_possible_json(rec[col])
        records.append(rec)

    # Ensure there's a 'sample' column; if not, fall back to using stem of BAM path if available
    if "sample" not in df.columns:
        if "bam" in df.columns:
            df["sample"] = df["bam"].apply(lambda x: Path(x).stem if x else "")
        else:
            df["sample"] = [f"sample_{i}" for i in range(len(df))]
    df["sample"] = df["sample"].astype(str)

    # Identify numeric metric columns for per-sample metrics page
    nums = numeric_cols_from_df(df, "sample")
    print("Detected numeric metrics:", nums)

    # Aggregate diagnostics (if requested) — each function returns path to created plot or None
    created = []
    if INCLUDE_AGGREGATE_DIAGNOSTICS:
        try:
            m1 = aggregate_mutational_spectrum(records, OUTPUT_DIR)
            if m1: created.append(m1)
        except Exception:
            traceback.print_exc()
        try:
            m2 = aggregate_genotype_stacked(records, OUTPUT_DIR)
            if m2: created.append(m2)
        except Exception:
            traceback.print_exc()
        try:
            m3 = aggregate_heatmap(df, "sample", OUTPUT_DIR, top_n=HEATMAP_TOP_N)
            if m3: created.append(m3)
        except Exception:
            traceback.print_exc()
        try:
            ma = aggregate_ma_plot(df, OUTPUT_DIR)
            if ma: created.append(ma)
        except Exception:
            traceback.print_exc()
        try:
            vol = aggregate_volcano(df, OUTPUT_DIR)
            if vol: created.append(vol)
        except Exception:
            traceback.print_exc()
        try:
            ph, qq = aggregate_pvalue_hist_qq(df, OUTPUT_DIR)
            if ph: created.append(ph)
            if qq: created.append(qq)
        except Exception:
            traceback.print_exc()

    # Per-sample pages & PDFs - iterate records which already contain parsed complex fields
    per_sample_pdfs = []
    for rec in records:
        sample = str(rec.get("sample") or "")
        if not sample:
            continue
        safe = sanitize_fname(sample)
        try:
            # build individual page figures
            fig_metrics = build_metrics_fig(sample, rec, nums)
            fig_vaf = build_vaf_fig(sample, rec)
            fig_mut = build_mutation_fig(sample, rec)
            fig_gt = build_genotype_fig(sample, rec)

            # paths for PNGs and PDF
            metrics_png = SAMPLES_DIR / f"{safe}_metrics.png"
            vaf_png = SAMPLES_DIR / f"{safe}_vaf.png"
            mut_png = SAMPLES_DIR / f"{safe}_mutational_spectrum.png"
            gt_png = SAMPLES_DIR / f"{safe}_genotype.png"
            pdf_path = SAMPLES_DIR / f"{safe}_summary.pdf"

            # Save individual PNGs
            fig_metrics.savefig(metrics_png, dpi=DPI, bbox_inches="tight"); plt.close(fig_metrics)
            fig_vaf.savefig(vaf_png, dpi=DPI, bbox_inches="tight"); plt.close(fig_vaf)
            fig_mut.savefig(mut_png, dpi=DPI, bbox_inches="tight"); plt.close(fig_mut)
            fig_gt.savefig(gt_png, dpi=DPI, bbox_inches="tight"); plt.close(fig_gt)

            # Create a multipage PDF by embedding the PNGs into PDF pages.
            # This method preserves layout and ensures the same visual pages are in PDF.
            with PdfPages(pdf_path) as pdf:
                for img in [metrics_png, vaf_png, mut_png, gt_png]:
                    fig_img = plt.figure(figsize=(11,8))
                    ax = fig_img.add_subplot(111)
                    im = plt.imread(str(img))
                    ax.imshow(im)
                    ax.axis("off")
                    pdf.savefig(fig_img, bbox_inches='tight')
                    plt.close(fig_img)

            per_sample_pdfs.append(pdf_path)
            print("Created PDF:", pdf_path)
        except Exception:
            # Traceback printed but script continues to next sample to avoid stopping a full run
            traceback.print_exc()

    # Copy per-sample PNGs into results parent directories for downstream per-sample inspection
    if COPY_PAGES_TO_RESULTS_PARENT:
        moved = copy_sample_pngs(SAMPLES_DIR, RESULTS_PARENT)
        if moved:
            manifest = RESULTS_PARENT / "moved_summary_manifest.txt"
            with manifest.open("w") as fh:
                for s,d in moved:
                    fh.write(f"{s}\t{d}\n")
            print("Wrote manifest:", manifest)

    # Final summary printed to stdout
    print("Done. Example PDFs:", per_sample_pdfs[:10])
    print("Aggregate plots created (if any):")
    for p in created[:40]:
        print(" ", p)

# Entry point
if __name__ == "__main__":
    main()
