#!/usr/bin/env python3
"""
upset_genes_annotated.py — heavily annotated

Purpose
-------
Load multiple gene identifier lists (one gene ID per line), produce:
 - human-readable overlap summary printed to stdout
 - optional CSV files listing pairwise and all-set overlaps + membership matrix
 - an UpSet plot saved to disk to visualize intersections among the lists

This annotated version explains every block of logic, assumptions, failure modes,
and provides suggestions for robustness, reproducibility, and extensions.

Usage
-----
Edit the LIST_DIR, SETS, and OUTPUT_DIR variables below (or convert them into
command-line arguments) and run:

    python3 upset_genes_annotated.py

Dependencies
------------
- Python packages: pandas, matplotlib, upsetplot
  Install via pip: pip install pandas matplotlib upsetplot
  Or conda: mamba install -c conda-forge pandas matplotlib upsetplot

Design notes / assumptions
--------------------------
- Each input list is a plain text file, one identifier per line. Lines starting
  with '#' are treated as comments and ignored. Whitespace around IDs is stripped.
- Identifiers are compared as plain strings (case-sensitive). If you need
  case-insensitive matching or ID normalization (e.g., Ensembl version removal),
  preprocess the files accordingly or modify load_gene_list().
- For very large lists (hundreds of thousands of IDs) memory use can grow when
  building the full membership matrix — see recommendations near the end.
- The UpSet plot library (upsetplot) accepts Python sets and will compute
  intersections. The script uses from_contents(sets) to create the data structure.

Outputs
-------
- UpSet plot PNG at OUTPUT_FILE
- If SAVE_TABLES is True:
  - directory OUTPUT_DIR/overlap_tables containing:
    - pairwise overlap CSVs (one per pair of input lists)
    - overlap__ALL_SETS.csv (genes common to all input sets) if any
    - all_genes_membership.csv (boolean membership matrix with one row per unique gene)
"""

from pathlib import Path
import matplotlib
matplotlib.use("Agg")   # Use non-interactive backend for headless servers
import matplotlib.pyplot as plt
import pandas as pd

# Try importing upsetplot with a helpful error message if missing
try:
    from upsetplot import from_contents, UpSet
except Exception as e:
    raise RuntimeError(
        "The 'upsetplot' package is required but not installed. "
        "Install it with 'pip install upsetplot' or 'mamba install -c conda-forge upsetplot'. "
        f"Original error: {e}"
    )

# ---------------------------------------------------------------------
# CONFIGURATION — edit these paths and labels to match your data
# ---------------------------------------------------------------------
# LIST_DIR: directory containing plain text lists (one gene per line).
LIST_DIR   = Path("/home/jgd/Documents/bioinformatics_working/2025-11-06_enu-proc/lists")

# OUTPUT_DIR: where the UpSet plot and overlap tables (optional) will be written.
OUTPUT_DIR = Path("/home/jgd/Documents/bioinformatics_working/output")

# SETS: mapping label -> filepath for each list to include in the UpSet analysis.
# Use descriptive labels. Avoid extremely long labels as they may crowd the plot.
SETS = {
    "List 1 Label":  LIST_DIR / "list1.txt",
    "List 2 Label":  LIST_DIR / "list2.txt",
    "List 3 Label":  LIST_DIR / "list3.txt",
    "List 4 Label":  LIST_DIR / "list4.txt",
    "List 5 Label":  LIST_DIR / "list5.txt",
    "List 6 Label":  LIST_DIR / "list6.txt",
}

# Title written above the UpSet plot
TITLE       = "Gene Overlap Analysis"

# Where to save the UpSet figure (PNG). The script will create the parent dir if needed.
OUTPUT_FILE = OUTPUT_DIR / "upset_plot.png"

# If True, the script will generate CSV tables for pairwise overlaps and the membership matrix.
SAVE_TABLES = True
# ---------------------------------------------------------------------

# ---------------------------------------------------------------------
# Helper: robustly load one gene list file into a Python set
# ---------------------------------------------------------------------
def load_gene_list(filepath: Path) -> set:
    """
    Read a file and return a set of identifiers.

    - Strips leading/trailing whitespace from lines.
    - Ignores blank lines and lines beginning with '#' (simple comment support).
    - Returns a Python set for efficient set operations.

    Edge cases and suggestions:
    - If your IDs contain whitespace or complex separators, consider splitting
      lines differently or using a TSV with a dedicated column.
    - If you need normalized IDs (e.g. uppercase, remove version suffixes),
      perform normalization here (e.g., id.split('.')[0]).
    """
    genes = set()
    if not filepath.exists():
        # Fail fast with actionable error message
        raise FileNotFoundError(f"Input list not found: {filepath}")
    with filepath.open("r", encoding="utf-8") as fh:
        for line in fh:
            s = line.strip()
            # Skip empty and comment lines
            if not s or s.startswith("#"):
                continue
            # Optionally normalize ID (uncomment and adapt if needed):
            # s = s.split('.')[0]   # remove Ensembl version suffix like ENSG00000.1
            genes.add(s)
    return genes

# ---------------------------------------------------------------------
# Pretty-print a short summary of each input set and common intersections
# ---------------------------------------------------------------------
def print_summary(sets: dict):
    """
    Print counts per set and a preview of genes common to ALL sets.
    The function uses set intersection and union operations.
    """
    print("\n" + "=" * 60)
    print("OVERLAP SUMMARY")
    print("=" * 60)
    for label, genes in sets.items():
        print(f"  {label:30s}: {len(genes):>6,} genes")

    # Intersection across all provided sets. If there is only one set, this is that set.
    all_common = set.intersection(*sets.values()) if len(sets) > 0 else set()
    print(f"\n  Common to ALL sets: {len(all_common):,} genes")
    if all_common:
        preview = sorted(all_common)[:20]
        print("  " + ", ".join(preview) + (" ..." if len(all_common) > 20 else ""))
    print("=" * 60 + "\n")


# ---------------------------------------------------------------------
# Save CSV tables that list overlaps and the full membership matrix
# ---------------------------------------------------------------------
def save_tables(sets: dict):
    """
    Create an 'overlap_tables' directory under OUTPUT_DIR with:
      - pairwise overlap CSVs named overlap__<A>__<B>.csv (label spaces replaced by underscores)
      - overlap__ALL_SETS.csv for genes in the intersection of all sets
      - all_genes_membership.csv: boolean membership matrix (rows = genes, cols = set labels)

    Notes:
    - Filenames are derived from set labels; if labels contain special chars,
      they are replaced with underscores for safety.
    - For very large gene universes the membership CSV can be large; consider
      writing compressed CSV (gzip) instead if storage is a concern.
    """
    from itertools import combinations
    tables_dir = OUTPUT_DIR / "overlap_tables"
    tables_dir.mkdir(parents=True, exist_ok=True)

    # Pairwise overlaps
    for a, b in combinations(sets.keys(), 2):
        overlap = sorted(sets[a] & sets[b])
        if overlap:
            safe_a = a.replace(" ", "_")
            safe_b = b.replace(" ", "_")
            fname = tables_dir / f"overlap__{safe_a}__{safe_b}.csv"
            pd.DataFrame({"gene": overlap}).to_csv(fname, index=False)
            print(f"  Saved: {fname}")

    # All-set common genes (intersection across all sets)
    if sets:
        all_common = sorted(set.intersection(*sets.values()))
        if all_common:
            fname = tables_dir / "overlap__ALL_SETS.csv"
            pd.DataFrame({"gene": all_common}).to_csv(fname, index=False)
            print(f"  Saved: {fname}")

    # Full membership matrix (one row per gene, boolean columns for each set)
    all_genes = set.union(*sets.values()) if sets else set()
    membership = pd.DataFrame(index=sorted(all_genes))
    for label, genes in sets.items():
        # membership.index.isin(genes) returns boolean Series aligned to the membership index
        membership[label] = membership.index.isin(genes)
    fname = tables_dir / "all_genes_membership.csv"
    membership.to_csv(fname)
    print(f"  Saved: {fname}")


# ---------------------------------------------------------------------
# Build and save UpSet plot
# ---------------------------------------------------------------------
def plot_upset(sets: dict):
    """
    Build and save an UpSet plot using upsetplot.from_contents and UpSet.
    Important plotting options:
      - subset_size="count": display subset sizes by count
      - show_counts=True: overlay raw counts on bars
      - sort_by="cardinality": sort intersection sets by their size
      - min_subset_size=1: exclude empty intersections

    Customization suggestions:
      - To focus on the top-k intersections set `subset_size='sum'` and supply
        a custom `sort_by` or `sort_by='degree'`.
      - To produce vector output, change OUTPUT_FILE suffix to .pdf or .svg and
        ensure Matplotlib supports that backend (Agg supports these).
    """
    # Create upsetplot data structure from a dict of label -> set
    data = from_contents(sets)

    # Create figure and configure UpSet object
    fig = plt.figure(figsize=(14, 6))
    upset = UpSet(
        data,
        subset_size="count",
        show_counts=True,
        sort_by="cardinality",
        min_subset_size=1,
    )
    # The UpSet.plot method expects a Matplotlib figure object
    upset.plot(fig)

    # Title placement: place above plot; adjust y to avoid clipping
    fig.suptitle(TITLE, fontsize=13, y=1.02)

    # Ensure output dir exists and save PNG
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    # If you prefer PDF output for vector graphics, set OUTPUT_FILE with .pdf extension
    fig.savefig(OUTPUT_FILE, bbox_inches="tight", dpi=150)
    print(f"\n  UpSet plot saved: {OUTPUT_FILE}")
    plt.close(fig)


# ---------------------------------------------------------------------
# Main orchestration
# ---------------------------------------------------------------------
def main():
    # Validate inputs and load sets; any missing file raises FileNotFoundError
    print(f"Loading {len(SETS)} gene lists from {LIST_DIR} ...")
    sets = {}
    for label, path in SETS.items():
        try:
            genes = load_gene_list(Path(path))
        except FileNotFoundError as e:
            # Fail fast with a helpful message listing which file is missing
            raise RuntimeError(f"Failed to load list for label '{label}': {e}")
        sets[label] = genes

    # Print a human-readable summary showing per-list counts + common genes
    print_summary(sets)

    # Optionally save CSV overlap tables and the membership matrix
    if SAVE_TABLES:
        print("Saving overlap tables to disk ...")
        save_tables(sets)

    # Build and save the UpSet plot
    print("Generating UpSet plot ...")
    plot_upset(sets)

    print("\nDone.")

# ---------------------------------------------------------------------
# If executed as a script, run main()
# ---------------------------------------------------------------------
if __name__ == "__main__":
    main()

# ---------------------------------------------------------------------
# POST-SCRIPT — Practical tips, troubleshooting, and extensions
# ---------------------------------------------------------------------
#
# Troubleshooting:
# - "ModuleNotFoundError: No module named 'upsetplot'": install upsetplot via pip/conda.
# - Empty plot or no intersections shown: check that the SETS mapping points to files
#   that actually contain identifiers and that identifiers are normalized consistently
#   across files (case, whitespace, version suffixes).
# - Overcrowded labels in figure: shorten labels or use fewer sets per plot.
#
# Suggested extensions:
# 1) CLI / configuration:
#    - Replace hard-coded variables with argparse to pass input directory, set list,
#      output directory, and toggles on the command line.
#
# 2) ID normalization:
#    - Add optional normalization function applied to each ID (e.g., t.upper(), t.split('.')[0]).
#
# 3) Subsetting:
#    - Allow user to limit the plot to the top-k most abundant intersections, or
#      to exclude extremely small intersections from visualization.
#
# 4) Interactive exploration:
#    - Save the membership matrix; build an interactive dashboard (Dash, Streamlit)
#      that allows filtering by sets and inspecting the gene lists in each intersection.
#
# 5) Large-scale optimization:
#    - For many sets (>20) UpSet plots become less informative. Consider pairwise
#      heatmaps of Jaccard index or hierarchical clustering of sets instead:
#        J(A,B) = |A ∩ B| / |A ∪ B|
#
# 6) Output formats:
#    - For publication-quality figures, save as SVG or PDF (vector) and tweak style:
#        fig.savefig("upset_plot.pdf", bbox_inches="tight")
#
# 7) Unit tests:
#    - Create small test files (3-4 genes each) and assert expected intersections
#      and membership matrix rows in saved CSVs.
#
# If you'd like, I can:
# - Convert this into a CLI tool (argparse) with advanced options (normalization, output formats),
# - Add unit tests and a small set of example input lists,
# - Provide a Jupyter notebook illustrating interactive exploration of the membership matrix.
#
# End of annotated file.
