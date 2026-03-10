#!/usr/bin/env python3
"""
SNP Density Ideogram — Zebrafish (Danio rerio, GRCz11) — Annotated
=================================================================
This is your original SNP density ideogram generator with heavy, line-by-line
and block-level annotations added. The goal is to make the script easier to
understand, maintain, validate, and extend.

What this script does (summary)
- Parses a VCF (gzipped or plain text) and collects heterozygous variant
  positions per sample (GT-based).
- Computes per-chromosome, fixed-size-window densities (SNP count per window).
- Merges user-defined sample pools (union of het positions).
- Produces:
    * Per-pool "density" ideograms (monochrome -> green progression)
    * Per-pair subtraction ideograms (diverging blue-white-red map; throm − nothrom)
    * Optional individual sample ideograms
- Saves figures to OUTPUT_BASE / timestamped folder.

High-level limitations & assumptions
- Only GT-based heterozygous calls are considered (FORMAT GT must be present).
- AD/DP-depth or allele-frequency filters are not applied — we trust caller to
  have filtered variants upstream. If you need depth-based filtering, add it
  to parse_vcf_per_sample().
- The script treats each heterozygote equally (no weighting by quality or VAF).
- Multi-allelic sites are treated only via GT field (the code does not split
  multi-allelic AD arrays into per-ALT calls).
- Chromosome naming expects 'chr1'..'chr25' or numeric '1'..'25' — the code will
  prepend 'chr' if missing. Non-standard contigs are ignored.
- The windowing is simple fixed-width bins. For other analyses (e.g., variable
  gene density), consider gene-based windows or GC/coverage normalization.
"""
# Standard library
import os
import gzip
from collections import defaultdict
from datetime import datetime

# Third-party — numerical & plotting
import numpy as np
import matplotlib
matplotlib.use("Agg")  # non-interactive backend suitable for HPC or CI
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import Normalize, LinearSegmentedColormap, TwoSlopeNorm

# ===========================================================================
#  SETTINGS  — edit these
# ===========================================================================
# Centralize tunable parameters. Prefer editing here or move to a YAML/JSON
# config for reproducible pipelines.

# Path to VCF file (gzipped VCF allowed)
VCF_FILE         = "/home/jgd/Documents/bioinformatics_working/2025-11-06_enu-proc/vcf/gatk/gatk-ensembl/vep_annotated/gatk-ensembl_filtered.20260224.nodash.vep.vcf.gz"

# WINDOW_SIZE: integer number of basepairs to collapse into a single density bin.
# Choose based on genome size & desired resolution. 1_000_000 (1 Mb) is coarse,
# useful for ideograms. For finer maps use 100_000 (100 kb) or 10_000 (10 kb)
# but be aware of memory/plotting consequences.
WINDOW_SIZE      = 1_000_000   # bp per bin (1 Mb)

# PERCENTILE_CAP: percentile at which the density colourmap saturates per pool.
# Lower values make the map saturate earlier; this avoids a few hotspots dominating
# the colour scale. Typical values: 75–99.
PERCENTILE_CAP   = 75

# Main title for figures (modify to include experiment/timepoint etc.)
TITLE            = "The number of SNPs within 1Mb window size"

# Output directory base — a timestamped subdirectory will be created here.
OUTPUT_BASE      = "/home/jgd/Documents/bioinformatics_working/output/snp_ideogram"

# Plot DPI — higher for publication quality (300+)
DPI              = 300

# Whether to create plots for each individual sample
PLOT_INDIVIDUALS = False

# ===========================================================================
#  SAMPLE POOLS — map pool name -> list of sample IDs as in VCF header
# ===========================================================================
# Pools can be used to aggregate signals (union of het positions). Pools with
# single members are also supported (treated like an individual).
SAMPLE_POOLS = {
    # ai08 ── throm / nothrom
    "ai08nothrom_pooled": [
        "ai08nothrom1",
        "ai08nothrom2",
        "ai08nothrom3",
    ],
    "ai08throm_pooled": [
        "ai08throm1",
        "ai08throm2",
        "ai08throm3",
    ],

    # bq01 ── throm / nothrom
    "bq01nothrom_pooled": [
        "bq01nothrom1",
        "bq01nothrom2",
        "bq01nothrom3",
    ],
    "bq01throm_pooled": [
        "bq01throm1",
        "bq01throm2",
        "bq01throm3",
    ],

    # bq06 ── throm / nothrom
    "bq06nothrom_pooled": [
        "bq06nothrom1",
        "bq06nothrom2",
        "bq06nothrom3",
    ],
    "bq06throm_pooled": [
        "bq06throm1",
        "bq06throm2",
        "bq06throm3",
    ],

    # pc ── no throm/nothrom split
    "pc_pooled": [
        "pc04",
        "pc05",
        "pc06",
    ],

    # wt ── no throm/nothrom split
    "wt_pooled": [
        "wt01",
        "wt02",
        "wt03",
    ],

    # Individuals — pools with single sample each side
    "au01nothrom": ["au01nothrom1"],
    "au01throm":   ["au01throm1"],
    "bo01nothrom": ["bo01nothrom1"],
    "bo01throm":   ["bo01throm1"],
    "bqm2nothrom": ["bqm2nothrom1"],
    "bqm2throm":   ["bqm2throm1"],
}

# ===========================================================================
#  SUBTRACTION PAIRS  — (throm_pool_name, nothrom_pool_name, output_label)
#  Each pair produces a diverging ideogram: throm − nothrom per window.
#  Blue  = nothrom has more het SNPs in that window
#  Red   = throm has more het SNPs in that window
# ===========================================================================
SUBTRACTION_PAIRS = [
    ("ai08throm_pooled",  "ai08nothrom_pooled",  "ai08_throm_minus_nothrom"),
    ("bq01throm_pooled",  "bq01nothrom_pooled",  "bq01_throm_minus_nothrom"),
    ("bq06throm_pooled",  "bq06nothrom_pooled",  "bq06_throm_minus_nothrom"),
    ("au01throm",         "au01nothrom",          "au01_throm_minus_nothrom"),
    ("bo01throm",         "bo01nothrom",          "bo01_throm_minus_nothrom"),
    ("bqm2throm",         "bqm2nothrom",          "bqm2_throm_minus_nothrom"),
]

# ===========================================================================
#  Zebrafish GRCz11 chromosome sizes (bp)
#  Important: sizes must match the assembly used to call variants.
# ===========================================================================
CHROM_SIZES = {
    "chr1":  58_044_006,
    "chr2":  59_545_122,
    "chr3":  64_804_009,
    "chr4":  62_045_047,
    "chr5":  62_070_934,
    "chr6":  60_406_253,
    "chr7":  75_753_327,
    "chr8":  56_349_073,
    "chr9":  56_724_980,
    "chr10": 45_117_831,
    "chr11": 46_279_481,
    "chr12": 49_438_694,
    "chr13": 52_002_682,
    "chr14": 46_786_950,
    "chr15": 47_855_872,
    "chr16": 53_323_654,
    "chr17": 54_320_357,
    "chr18": 51_520_288,
    "chr19": 49_120_674,
    "chr20": 56_398_776,
    "chr21": 45_024_088,
    "chr22": 47_768_291,
    "chr23": 45_365_134,
    "chr24": 43_214_605,
    "chr25": 44_073_028,
}

# CHROM_ORDER defines plotting order; here it's chr1..chr25
CHROM_ORDER = [f"chr{i}" for i in range(1, 26)]

# Colour maps — these are specified as hex colours, then converted to LinearSegmentedColormap
# Density colourmap: white → dark green → yellow → orange → red (for hotspots)
DENSITY_CMAP_COLORS = [
    "#F0F0F0",
    "#1A5C1A",
    "#2E8B2E",
    "#5AAF2A",
    "#9DC93B",
    "#E8E020",
    "#F0A500",
    "#F07800",
    "#E84800",
    "#E02000",
    "#CC0000",
]

# Diverging colourmap: blue → white → red. Blue indicates nothrom > throm (negative diffs),
# red indicates throm > nothrom (positive diffs).
DIVERGING_CMAP_COLORS = [
    "#2166AC",  # strong blue  (nothrom >> throm)
    "#4393C3",
    "#92C5DE",
    "#D1E5F0",
    "#F7F7F7",  # white        (no difference)
    "#FDDBC7",
    "#F4A582",
    "#D6604D",
    "#B2182B",  # strong red   (throm >> nothrom)
]

# Heterozygous genotype strings we accept in the GT field.
# We include both phased (|) and unphased (/) representations in either order.
HET_GENOTYPES = {"0/1", "1/0", "0|1", "1|0"}

# ===========================================================================
#  Helpers — I/O, parsing, pooling, density calculation, colormap builders
# ===========================================================================

def open_vcf(path):
    """
    Lightweight VCF opener that handles gzipped (.gz) VCFs (text mode).
    Note: this reads line-by-line and DOES NOT use pysam — for very large files
    pysam.VariantFile is faster and uses indexed access. This function is simple
    and portable.
    """
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def is_het(genotype_field):
    """
    Simple GT check to detect heterozygous calls.
    genotype_field is a single-sample genotype string from VCF columns[9+].
    FORMAT expected to contain GT as first field (GT:...).
    Returns True for 0/1, 1/0, 0|1, 1|0.
    Caveats:
      - Phased versus unphased separators handled.
      - Does not support multi-allelic alternatives like "0/2" or "1/2".
      - If GT uses '.' missing alleles (e.g., ./.), result is False.
    """
    gt = genotype_field.split(":")[0]
    return gt in HET_GENOTYPES


def parse_vcf_per_sample(path):
    """
    Parse the VCF and collect heterozygous variant positions per sample.

    Returns:
      sample_names: list of sample IDs (order preserved from VCF)
      data: dict mapping sample -> dict(chrom -> sorted list of positions)

    Implementation notes:
      - We only keep position (POS) for heterozygous genotypes. This is memory
        efficient compared to storing complete records, but still could be large
        for whole-genome VCFs with many hets.
      - We skip variants if FORMAT does not start with GT (robust but strict).
      - We normalize chromosome names by prepending 'chr' if missing (e.g., '1' -> 'chr1').
      - Non-chromosome contigs or chromosomes not in CHROM_SIZES are ignored.
    """
    sample_names = []
    data = {}

    with open_vcf(path) as fh:
        for line in fh:
            # Header line with sample names: '#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT sample1 sample2 ...'
            if line.startswith("#CHROM"):
                cols = line.rstrip("\n").split("\t")
                sample_names = cols[9:]
                # Using defaultdict(list) for positions by chromosome; will be converted to dict later
                data = {s: defaultdict(list) for s in sample_names}
                print(f"  Samples ({len(sample_names)}): {sample_names}")
                continue

            # Skip other header lines
            if line.startswith("#"):
                continue

            parts = line.rstrip("\n").split("\t")
            # If a non-standard VCF with <9 columns, skip defensively
            if len(parts) < 10:
                continue

            chrom   = parts[0]
            pos_str = parts[1]
            fmt     = parts[8]

            # Normalize chromosome naming: some VCFs use '1'..'25', others 'chr1'..'chr25'
            if not chrom.startswith("chr"):
                chrom = "chr" + chrom
            # Skip contigs not in CHROM_SIZES to avoid plotting unexpected scaffolds
            if chrom not in CHROM_SIZES:
                continue
            # Ensure FORMAT starts with GT — if GT is not first, this check is conservative.
            # A more flexible approach would parse FORMAT and find the GT index.
            if not fmt.startswith("GT"):
                continue

            # Parse POS as int; skip malformed positions
            try:
                pos = int(pos_str)
            except ValueError:
                continue

            # For each sample column, check if sample has het GT and record position
            for sample, gf in zip(sample_names, parts[9:]):
                if is_het(gf):
                    data[sample][chrom].append(pos)

    # Post-process: sort positions per chromosome and convert to plain dicts
    for sample in data:
        for chrom in data[sample]:
            data[sample][chrom].sort()
        data[sample] = dict(data[sample])  # convert defaultdict to dict

    return sample_names, data


def pool_samples(data, pool_members):
    """
    Merge (union) heterozygous positions across a list of sample names.
    Returns a dict chrom -> sorted unique positions.

    Important details:
      - Pooling uses a set to deduplicate positions that may appear in multiple members.
      - If a pool member is missing from 'data' (e.g., sample not in VCF),
        the function warns and skips that member.
      - This pooling is positional union — if you prefer intersections (shared hets),
        change 'merged[chrom].update(positions)' to an intersection logic.
    """
    merged = defaultdict(set)
    for sample in pool_members:
        if sample not in data:
            print(f"  ⚠ Pool member '{sample}' not found in VCF — skipping.")
            continue
        for chrom, positions in data[sample].items():
            merged[chrom].update(positions)
    # Return sorted lists for reproducible ordering
    return {chrom: sorted(pos_set) for chrom, pos_set in merged.items()}


def compute_density(positions, chrom_size, window):
    """
    Compute SNP counts per window for a single chromosome.

    Args:
      positions: sorted list of integer variant positions (1-based)
      chrom_size: length of chromosome in bp
      window: window size in bp

    Returns:
      numpy array of length n_bins with counts per bin (dtype int32)

    Implementation notes:
      - Ceil division for number of bins: n_bins = ceil(chrom_size / window).
      - Positions are binned with integer division (pos // window), which maps 1..window to bin 0..0
        Note: this treats positions as 1-based but uses floor division; it's consistent for fixed bins.
      - np.add.at is used to efficiently accumulate counts into counts array.
    """
    n_bins = max(1, -(-chrom_size // window))  # ceil division: avoid zero bins
    counts = np.zeros(n_bins, dtype=np.int32)
    if positions:
        pos_arr = np.array(positions, dtype=np.int64)
        # bin indices: convert positions to 0-based bins; clip to valid range
        bin_idx = np.clip((pos_arr // window).astype(np.int32), 0, n_bins - 1)
        np.add.at(counts, bin_idx, 1)
    return counts


def build_density_cmap():
    """Build a LinearSegmentedColormap from DENSITY_CMAP_COLORS for density plots."""
    return LinearSegmentedColormap.from_list("snp_density", DENSITY_CMAP_COLORS)


def build_diverging_cmap():
    """Build a LinearSegmentedColormap from DIVERGING_CMAP_COLORS for subtraction plots."""
    return LinearSegmentedColormap.from_list("snp_diverging", DIVERGING_CMAP_COLORS)


def percentile_norm(density_dict, percentile=95):
    """
    Build a Normalize instance (vmin=0, vmax=percentile of non-zero bin counts)
    so that colour saturation is robust to outliers.

    Returns:
      Normalize instance, vmin (0.0), vmax (float)
    If there are no non-zero bins, return a default Normalize(0,1).
    """
    all_nonzero = []
    for counts in density_dict.values():
        nz = counts[counts > 0]
        all_nonzero.extend(nz.tolist())
    if not all_nonzero:
        # Avoid division by zero and give a safe default scale
        return Normalize(vmin=0, vmax=1), 0.0, 1.0
    vmax = float(np.percentile(all_nonzero, percentile))
    vmax = max(vmax, 1.0)
    return Normalize(vmin=0, vmax=vmax, clip=True), 0.0, vmax


def diverging_norm(diff_dict, percentile=99):
    """
    Build a TwoSlopeNorm centered on 0, with symmetric vmin/vmax
    determined by a chosen percentile of absolute differences.

    Returns:
      TwoSlopeNorm instance, vmin (neg float), vmax (pos float)
    """
    all_vals = []
    for arr in diff_dict.values():
        all_vals.extend(arr.tolist())
    if not all_vals:
        return TwoSlopeNorm(vmin=-1, vcenter=0, vmax=1), -1, 1

    all_vals = np.array(all_vals)
    abs_max  = float(np.percentile(np.abs(all_vals), percentile))
    abs_max  = max(abs_max, 1.0)
    norm     = TwoSlopeNorm(vmin=-abs_max, vcenter=0, vmax=abs_max)
    return norm, -abs_max, abs_max


def print_diagnostics(label, density):
    """
    Print summary statistics for a density map (dictionary chrom -> counts array).
    Useful for quick QC: totals, distribution percentiles, colour saturation threshold.
    """
    total_snps  = sum(c.sum() for c in density.values())
    all_nonzero = []
    for counts in density.values():
        all_nonzero.extend(counts[counts > 0].tolist())
    if not all_nonzero:
        print(f"  {label}: no het SNPs found")
        return
    arr = np.array(all_nonzero)
    print(f"\n── {label} ────────────────────────────────────────────────────")
    print(f"  Total het SNPs : {total_snps:,}")
    print(f"  Window stats   : min={arr.min()}  median={np.median(arr):.1f}  "
          f"mean={arr.mean():.1f}  95th%={np.percentile(arr,95):.1f}  "
          f"99th%={np.percentile(arr,99):.1f}  max={arr.max()}")
    print(f"  Colour saturates at {PERCENTILE_CAP}th pct = "
          f"{np.percentile(arr, PERCENTILE_CAP):.1f} SNPs/window")
    print(f"────────────────────────────────────────────────────────────────")


# ===========================================================================
#  Density ideogram plotting
#  - This function draws chromosomes vertically with windows coloured by density.
# ===========================================================================
def plot_ideogram(label, density, output_folder, pool_members=None):
    """
    label: str used in title and file name
    density: dict chrom -> numpy array of counts per window
    output_folder: where to save PNG
    pool_members: optional list of samples that compose this pool (used in subtitle)
    """
    # Ordered chromosomes to include only those with data > 0
    ordered = [c for c in CHROM_ORDER
               if c in density and density[c].sum() > 0]
    if not ordered:
        print(f"  ⚠ No data for '{label}', skipping.")
        return

    # Build color scale and percentile-based normalization
    cmap             = build_density_cmap()
    norm, vmin, vmax = percentile_norm(density, PERCENTILE_CAP)

    # Layout geometry — tuned by eye for legibility
    max_chrom_mb  = max(CHROM_SIZES[c] for c in ordered) / 1e6
    n_chroms      = len(ordered)
    fig_w         = 13
    bar_h         = 0.36
    row_gap       = 0.14
    top_margin    = 1.7 if pool_members and len(pool_members) > 1 else 1.5
    bottom_margin = 0.5
    right_margin  = 2.8

    fig_h = top_margin + n_chroms * (bar_h + row_gap) + bottom_margin
    fig, ax = plt.subplots(figsize=(fig_w, fig_h), dpi=DPI)
    ax.set_axis_off()
    label_x = -2.5  # x-position for chromosome labels

    def chrom_y(idx):
        """Convert chromosome index to y coordinate for plotting."""
        return fig_h - top_margin - (idx + 1) * (bar_h + row_gap) + row_gap

    # Scale axis on top for Mb labels
    scale_y = fig_h - top_margin + 0.08
    ax.plot([0, max_chrom_mb], [scale_y, scale_y],
            color="black", linewidth=1.0, clip_on=False)
    for t in np.arange(0, max_chrom_mb + 1, 10):
        ax.plot([t, t], [scale_y, scale_y + 0.06],
                color="black", linewidth=0.8, clip_on=False)
        ax.text(t, scale_y + 0.14, f"{int(t)}Mb",
                ha="center", va="bottom", fontsize=8)

    # Title & subtitle
    subtitle  = f"pooled: {', '.join(pool_members)}" if pool_members and len(pool_members) > 1 else ""
    title_str = f"{TITLE}\n{label}"
    if subtitle:
        title_str += f"\n{subtitle}"
    ax.text(max_chrom_mb / 2, fig_h - 0.15, title_str,
            ha="center", va="top", fontsize=11,
            fontweight="bold", linespacing=1.5)

    # Draw chromosome bars and per-window coloured rectangles
    for idx, chrom in enumerate(ordered):
        y_bot   = chrom_y(idx)
        size_mb = CHROM_SIZES[chrom] / 1e6
        counts  = density[chrom]

        # Light rounded background box for chromosome
        ax.add_patch(mpatches.FancyBboxPatch(
            (0, y_bot), size_mb, bar_h,
            boxstyle="round,pad=0.0", linewidth=0.4,
            edgecolor="#AAAAAA", facecolor="#F0F0F0", zorder=1))

        # For each non-zero bin, draw a rectangle coloured according to norm(count)
        for bin_i, count in enumerate(counts):
            if count == 0:
                continue
            x0 = bin_i * WINDOW_SIZE / 1e6
            x1 = min((bin_i + 1) * WINDOW_SIZE / 1e6, size_mb)
            ax.add_patch(mpatches.Rectangle(
                (x0, y_bot), x1 - x0, bar_h,
                linewidth=0, facecolor=cmap(norm(count)), zorder=2))

        # Outlined border for chromosome bar
        ax.add_patch(mpatches.FancyBboxPatch(
            (0, y_bot), size_mb, bar_h,
            boxstyle="round,pad=0.0", linewidth=0.5,
            edgecolor="#555555", facecolor="none", zorder=3))

        # Chromosome label at left
        ax.text(label_x, y_bot + bar_h / 2, chrom,
                ha="right", va="center", fontsize=9)

    ax.set_xlim(label_x - 1, max_chrom_mb + right_margin)
    ax.set_ylim(0, fig_h)

    # Colourbar drawn manually as a vertical gradient of small rectangles so we can
    # position it aligned with the chromosome block (visual consistency).
    legend_x  = max_chrom_mb + 0.5
    y_bot_leg = chrom_y(n_chroms - 1)
    y_top_leg = chrom_y(0) + bar_h
    cb_h      = y_top_leg - y_bot_leg
    cb_w      = 0.4

    for i in range(100):
        frac = i / 100
        ax.add_patch(mpatches.Rectangle(
            (legend_x, y_bot_leg + frac * cb_h), cb_w, cb_h / 100,
            linewidth=0, facecolor=cmap(frac), zorder=4))

    ax.add_patch(mpatches.Rectangle(
        (legend_x, y_bot_leg), cb_w, cb_h,
        linewidth=0.5, edgecolor="#555555", facecolor="none", zorder=5))

    # Label colourbar ticks: 0, middle, >=vmax
    for frac, lbl in [(0.0, "0"), (0.5, f"{vmax/2:.0f}"), (1.0, f"≥{vmax:.0f}")]:
        ax.text(legend_x + cb_w + 0.2, y_bot_leg + frac * cb_h,
                lbl, ha="left", va="center", fontsize=8)

    ax.text(legend_x + cb_w / 2, y_bot_leg - 0.25,
            "SNPs/Mb", ha="center", va="top", fontsize=7)

    # Sanitize filename and save
    safe_name = label.replace("/", "_").replace(" ", "_")
    out_path  = os.path.join(output_folder, f"{safe_name}_snp_density.png")
    plt.tight_layout(pad=0)
    plt.savefig(out_path, dpi=DPI, bbox_inches="tight",
                facecolor="white", edgecolor="none")
    plt.close(fig)
    print(f"  ✓ Saved: {out_path}")


# ===========================================================================
#  Subtraction ideogram (diverging blue/red colourmap)
#  - diff_density: dict chrom -> array of (throm_count - nothrom_count)
# ===========================================================================
def plot_subtraction(label, diff_density, output_folder, throm_name, nothrom_name):
    """
    Produce a diverging ideogram showing where throm samples have more het SNPs
    (red) versus nothrom (blue). White represents zero difference.

    Notes:
      - diff_density should be signed integer arrays (positive => throm > nothrom).
      - TwoSlopeNorm centers the colormap on zero and scales symmetrically by percentile.
      - Very large outliers can still influence the percentile-based caps; inspect diagnostics.
    """
    ordered = [c for c in CHROM_ORDER if c in diff_density]
    if not ordered:
        print(f"  ⚠ No data for subtraction '{label}', skipping.")
        return

    cmap             = build_diverging_cmap()
    norm, vmin, vmax = diverging_norm(diff_density, percentile=99)

    max_chrom_mb  = max(CHROM_SIZES[c] for c in ordered) / 1e6
    n_chroms      = len(ordered)
    fig_w         = 13
    bar_h         = 0.36
    row_gap       = 0.14
    top_margin    = 1.9
    bottom_margin = 0.5
    right_margin  = 3.2

    fig_h = top_margin + n_chroms * (bar_h + row_gap) + bottom_margin
    fig, ax = plt.subplots(figsize=(fig_w, fig_h), dpi=DPI)
    ax.set_axis_off()
    label_x = -2.5

    def chrom_y(idx):
        return fig_h - top_margin - (idx + 1) * (bar_h + row_gap) + row_gap

    # scale axis (Mb ticks)
    scale_y = fig_h - top_margin + 0.08
    ax.plot([0, max_chrom_mb], [scale_y, scale_y],
            color="black", linewidth=1.0, clip_on=False)
    for t in np.arange(0, max_chrom_mb + 1, 10):
        ax.plot([t, t], [scale_y, scale_y + 0.06],
                color="black", linewidth=0.8, clip_on=False)
        ax.text(t, scale_y + 0.14, f"{int(t)}Mb",
                ha="center", va="bottom", fontsize=8)

    # Title: make clear interpretation of colours
    ax.text(max_chrom_mb / 2, fig_h - 0.15,
            f"SNP Density Difference (throm − nothrom)\n{label}\n"
            f"red = more in throm  |  blue = more in nothrom",
            ha="center", va="top", fontsize=11,
            fontweight="bold", linespacing=1.5)

    # Draw each chromosome row with signed rectangles coloured by norm(diff)
    for idx, chrom in enumerate(ordered):
        y_bot   = chrom_y(idx)
        size_mb = CHROM_SIZES[chrom] / 1e6
        diffs   = diff_density[chrom]

        # Base rectangle background
        ax.add_patch(mpatches.FancyBboxPatch(
            (0, y_bot), size_mb, bar_h,
            boxstyle="round,pad=0.0", linewidth=0.4,
            edgecolor="#AAAAAA", facecolor="#F7F7F7", zorder=1))

        for bin_i, diff in enumerate(diffs):
            if diff == 0:
                continue
            x0 = bin_i * WINDOW_SIZE / 1e6
            x1 = min((bin_i + 1) * WINDOW_SIZE / 1e6, size_mb)
            ax.add_patch(mpatches.Rectangle(
                (x0, y_bot), x1 - x0, bar_h,
                linewidth=0, facecolor=cmap(norm(diff)), zorder=2))

        # Outline
        ax.add_patch(mpatches.FancyBboxPatch(
            (0, y_bot), size_mb, bar_h,
            boxstyle="round,pad=0.0", linewidth=0.5,
            edgecolor="#555555", facecolor="none", zorder=3))

        # Chromosome label
        ax.text(label_x, y_bot + bar_h / 2, chrom,
                ha="right", va="center", fontsize=9)

    ax.set_xlim(label_x - 1, max_chrom_mb + right_margin)
    ax.set_ylim(0, fig_h)

    # Diverging colourbar drawn similarly to density bar
    legend_x  = max_chrom_mb + 0.5
    y_bot_leg = chrom_y(n_chroms - 1)
    y_top_leg = chrom_y(0) + bar_h
    cb_h      = y_top_leg - y_bot_leg
    cb_w      = 0.4

    for i in range(100):
        frac = i / 100
        ax.add_patch(mpatches.Rectangle(
            (legend_x, y_bot_leg + frac * cb_h), cb_w, cb_h / 100,
            linewidth=0, facecolor=cmap(frac), zorder=4))

    ax.add_patch(mpatches.Rectangle(
        (legend_x, y_bot_leg), cb_w, cb_h,
        linewidth=0.5, edgecolor="#555555", facecolor="none", zorder=5))

    # Label colourbar ticks: ≤vmin, 0, ≥vmax (vmin negative)
    for frac, lbl in [(0.0, f"≤{vmin:.0f}"), (0.5, "0"), (1.0, f"≥+{vmax:.0f}")]:
        ax.text(legend_x + cb_w + 0.2, y_bot_leg + frac * cb_h,
                lbl, ha="left", va="center", fontsize=8)

    ax.text(legend_x + cb_w / 2, y_bot_leg - 0.35,
            "SNPs/Mb\n(throm−nothrom)", ha="center", va="top", fontsize=7)

    safe_name = label.replace("/", "_").replace(" ", "_")
    out_path  = os.path.join(output_folder, f"{safe_name}_subtraction.png")
    plt.tight_layout(pad=0)
    plt.savefig(out_path, dpi=DPI, bbox_inches="tight",
                facecolor="white", edgecolor="none")
    plt.close(fig)
    print(f"  ✓ Saved: {out_path}")


# ===========================================================================
#  Main orchestration
#  - parse VCF, compute densities for individuals and pools, build subtraction maps,
#    print diagnostics, and create plots.
# ===========================================================================
def main():
    # Create timestamped output folder to avoid overwriting previous runs.
    timestamp     = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    output_folder = os.path.join(OUTPUT_BASE, timestamp)
    os.makedirs(output_folder, exist_ok=True)
    print(f"Output folder: {output_folder}")

    # Parse the VCF for per-sample heterozygous positions
    print(f"\nReading VCF: {VCF_FILE} ...")
    sample_names, data = parse_vcf_per_sample(VCF_FILE)

    # Quick per-sample counts to let the user assess sample-level call rates
    print(f"\nPer-sample het SNP counts:")
    for s in sample_names:
        total = sum(len(v) for v in data[s].values())
        print(f"  {s}: {total:,}")

    # ── compute density per individual sample ──────────────────────────────
    # all_densities[sample][chrom] = numpy array of counts per bin
    all_densities = {}
    for sample in sample_names:
        all_densities[sample] = {
            c: compute_density(data[sample].get(c, []), CHROM_SIZES[c], WINDOW_SIZE)
            for c in CHROM_ORDER
        }

    # ── compute density for each pool ─────────────────────────────────────
    pool_densities = {}
    for pool_name, members in SAMPLE_POOLS.items():
        print(f"\nPooling '{pool_name}': {members}")
        pooled_positions = pool_samples(data, members)
        pool_densities[pool_name] = {
            c: compute_density(pooled_positions.get(c, []), CHROM_SIZES[c], WINDOW_SIZE)
            for c in CHROM_ORDER
        }
        total = sum(len(v) for v in pooled_positions.values())
        print(f"  Pooled het SNP sites: {total:,}")

    # ── compute subtraction maps for user-specified pairs ──────────────────
    subtraction_maps = {}
    for throm_name, nothrom_name, sub_label in SUBTRACTION_PAIRS:
        if throm_name not in pool_densities:
            print(f"  ⚠ '{throm_name}' not in pool_densities — skipping subtraction.")
            continue
        if nothrom_name not in pool_densities:
            print(f"  ⚠ '{nothrom_name}' not in pool_densities — skipping subtraction.")
            continue
        diff = {}
        for chrom in CHROM_ORDER:
            # For safety, default to zero-length arrays sized by chromosome/window.
            t_counts = pool_densities[throm_name].get(chrom,
                np.zeros(max(1, -(-CHROM_SIZES[chrom] // WINDOW_SIZE)), dtype=np.int32))
            n_counts = pool_densities[nothrom_name].get(chrom,
                np.zeros(max(1, -(-CHROM_SIZES[chrom] // WINDOW_SIZE)), dtype=np.int32))
            # Align lengths (should be identical but guard anyway)
            min_len = min(len(t_counts), len(n_counts))
            diff[chrom] = (t_counts[:min_len].astype(np.int32)
                           - n_counts[:min_len].astype(np.int32))
        subtraction_maps[sub_label] = (diff, throm_name, nothrom_name)
        print(f"  Subtraction map computed: {sub_label}")

    # ── diagnostics ───────────────────────────────────────────────────────
    if PLOT_INDIVIDUALS:
        for sample in sample_names:
            print_diagnostics(sample, all_densities[sample])
    for pool_name in pool_densities:
        print_diagnostics(f"{pool_name} [POOL]", pool_densities[pool_name])

    # ── plot individuals ───────────────────────────────────────────────────
    if PLOT_INDIVIDUALS:
        print(f"\nPlotting {len(sample_names)} individual sample(s)...")
        for sample in sample_names:
            plot_ideogram(sample, all_densities[sample], output_folder)
    else:
        print(f"\nSkipping individual sample plots (PLOT_INDIVIDUALS = False)")

    # ── plot pools ─────────────────────────────────────────────────────────
    if SAMPLE_POOLS:
        print(f"\nPlotting {len(SAMPLE_POOLS)} pool ideogram(s)...")
        for pool_name, members in SAMPLE_POOLS.items():
            plot_ideogram(pool_name, pool_densities[pool_name],
                          output_folder, pool_members=members)

    # ── plot subtraction maps ──────────────────────────────────────────────
    if subtraction_maps:
        print(f"\nPlotting {len(subtraction_maps)} subtraction ideogram(s)...")
        for sub_label, (diff, throm_name, nothrom_name) in subtraction_maps.items():
            plot_subtraction(sub_label, diff, output_folder, throm_name, nothrom_name)

    print(f"\n✓ All done. Output in: {output_folder}")


# Module guard — allows importing functions for unit tests / interactive use
if __name__ == "__main__":
    main()

# ===========================================================================
# POST-SCRIPT: Suggestions, caveats, and possible extensions (not executed)
# ===========================================================================
# - Performance:
#   * For very large VCFs, using pysam.VariantFile and streaming only records
#     with GT present is faster and more memory efficient than pure text parsing.
#   * Consider limiting parse_vcf_per_sample to only variants on CHROM_ORDER
#     using region queries if your VCF is indexed (.tbi) and you have pysam.
#   * If window size is small, numpy arrays may become large; profile memory usage.
#
# - Filters & robustness:
#   * Add optional allele depth (AD) or DP filtering by parsing FORMAT and finding
#     the GT index dynamically (split FORMAT by ':' and find GT position).
#   * Support multi-allelic genotypes explicitly (e.g., detect "0/2", "1/2") if needed.
#   * Allow choosing union vs. intersection when pooling samples.
#
# - Plotting improvements:
#   * Use colorblind-friendly palettes or allow user to specify palettes.
#   * Add optional chromosome ordering by physical length or custom karyotype.
#   * Add options to output vector format (PDF/SVG) for downstream editing.
#   * For small-window, high-resolution plots, consider aggregating only chromosomes of interest.
#
# - CLI & config:
#   * Add argparse to expose VCF_FILE, WINDOW_SIZE, PERCENTILE_CAP, OUTPUT_BASE, and toggles.
#   * Move sample/pair configuration to YAML or JSON to separate code from data.
#
# - Tests:
#   * Create small synthetic VCF fixtures (few variants across chr1..chr3) to test:
#       - parse_vcf_per_sample handles missing GT/FORMAT correctly
#       - pool_samples merges sets correctly
#       - compute_density bins positions as expected across bin boundaries
#   * Use pytest to assert outputs (e.g., expected arrays).
#
# - Downstream:
#   * Save numeric density arrays to .npz or CSV so plotting can be redone without re-parsing the VCF.
#   * Add optional normalization by callable sites or coverage if you have per-window coverage maps.
#
# If you'd like, I can:
# - Convert this script to accept command-line arguments (argparse) and a YAML config file.
# - Replace the VCF text parser with pysam.VariantFile for speed and more robust FORMAT parsing.
# - Add AD/DP-based genotype filtering or VAF thresholds in parse_vcf_per_sample.
# - Produce unit tests (pytest) and small VCF fixtures to validate behavior.
# Tell me which of those you prefer and I'll implement it next.
