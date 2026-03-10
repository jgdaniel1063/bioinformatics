#!/usr/bin/env python3
"""
ENU Mutagenesis RNA-seq VCF Analyzer — Modifier/Suppressor Screen (Annotated)
==========================================================================
This file is your original analyzer with excessive, line-by-line annotations.
Annotations include:
 - Purpose and rationale for each block and threshold
 - Input/output expectations and shapes
 - Potential failure modes and how to detect/fix them
 - Notes on performance, scaling, and reproducibility
 - Suggestions for additional QC steps and unit tests

Keep in mind:
 - Comments and docstrings do not change runtime behavior (they are ignored by Python).
 - Where behaviour may be ambiguous (missing fields in VCF, nonstandard FORMAT tags),
   the code attempts to handle common cases but you should validate your VCF beforehand
   (e.g., using bcftools, GATK ValidateVariants, or a small script to check FORMAT tags).

Author: Annotated by ChatGPT for jgdaniel1063
Date: 2026-03-10 (provided in chat)
"""

# ===================================================================
# STANDARD LIBRARY IMPORTS
# ===================================================================
# These are standard modules used for file operations, system exit,
# container types, and type hints. Keep them minimal and explicit.
import os
import sys
from collections import Counter
from dataclasses import dataclass, field
from typing import Optional, List, Dict

# ===================================================================
# THIRD-PARTY DEPENDENCIES
# ===================================================================
# These must be installed in your Python environment. The script earlier
# lists pip requirements: pysam, pandas, matplotlib, numpy
#
# pysam: used to iterate through VCF records and access per-sample FORMAT fields.
# pandas/numpy: for tabular data summarization and numeric operations.
# matplotlib: to create figures for QC / exploratory plots.
#
# Important: pysam expects a properly indexed VCF (.tbi) when using some random-access
# features (not used here), but opening a gzipped VCF requires pysam to be compiled with
# bgzip / htslib support (standard in most installs).
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")  # non-interactive backend for headless servers / batch runs
import matplotlib.pyplot as plt
import pysam  # interface to htslib — use to read VCF/BCF files efficiently


# ===================================================================
# CONFIG — EDIT THESE
# ===================================================================
# All paths/parameters are centralized so you can tune the analysis quickly.
# Note the original script used absolute paths — consider switching to a small
# YAML/JSON config file for reproducible pipelines or integrate with snakemake.

INPUT_VCF = "/home/jgd/Documents/bioinformatics_working/2025-11-06_enu-proc/vcf/gatk/gatk-ensembl/vep_annotated/gatk-ensembl_filtered.20260224.nodash.vep.vcf.gz"
OUTPUT_DIR = "/home/jgd/Documents/bioinformatics_working/output/enu_suppressor_analysis_output"
ALL_VARIANTS_TSV = os.path.join(OUTPUT_DIR, "all_variants.tsv")
CANDIDATES_TSV = os.path.join(OUTPUT_DIR, "suppressor_candidates.tsv")
PLOT_PATH = os.path.join(OUTPUT_DIR, "enu_suppressor_analysis.png")

# -------------------------------------------------------------------
# SAMPLE SHEET
# -------------------------------------------------------------------
# Specify exact sample names that appear in the VCF header. Names must match
# the VCF sample IDs verbatim — otherwise the script will exit with an error.
#
# If you have many samples, prefer to load sample lists from a file to avoid
# typos and to enable sample reuse across scripts / pipelines.
NO_THROMBOSIS_SAMPLES = [
    "bq01nothrom1",
    "bq01nothrom2",
    "bq01nothrom3",
]

THROMBOSIS_SAMPLES = [
    "bq01throm1",
    "bq01throm2",
    "bq01throm3",
]

# -------------------------------------------------------------------
# GROUP CONCORDANCE
# -------------------------------------------------------------------
# Defines how strict you are about replicates agreeing on genotype pattern.
# Example:
#  - 1.0 requires all replicates show the expected genotype
#  - 0.67 for 3 replicates demands at least 2 agree (majority)
MIN_GROUP_CONCORDANCE = 0.67

# -------------------------------------------------------------------
# THRESHOLDS
# -------------------------------------------------------------------
# MIN_DEPTH rationale:
#   Set to a high value for pooled samples or many embryos (100x per embryo * number of embryos)
#   Original comment suggested MIN_DEPTH of 4000 for n=40 embryos (100x each), but this script
#   uses 400 as default — keep consistent with your experimental design.
#
# NOTE: In your supplied script, the header says 4000 but constant is 400 — double-check.
MIN_DEPTH = 400                # adjust according to your pool size: 100x per embryo * n embryos
MIN_CANDIDATE_SCORE = 50.0
REQUIRE_ENU_TIER_MAX = 4        # allow tiers 1-4 for general scoring; candidate filter narrows further

# GATK RNA-seq hard filter parameters (tightened values)
QD_MIN = 5.0
FS_MAX = 20.0
MQ_MIN = 50.0
MQRANKSUM_MIN = -5.0
READPOSRANKSUM_MIN = -4.0
SOR_MAX = 2.0

# ASE (allele-specific expression) skew threshold:
# Flag heterozygous sites in the "no-thrombosis" group whose VAF is far from 0.5
ASE_SKEW_THRESHOLD = 0.20

# ReadPosRankSum threshold for splice junction artifacts:
SPLICE_RPRS_CUTOFF = -3.0


# ===================================================================
# CONSTANTS — ENU TIER DEFINITIONS and SIGNATURES
# ===================================================================
# These sets define which base substitutions are most consistent with
# ENU mutagenesis on a collapsed-strand basis (A>G / T>C is the classical ENU signature).
ENU_TIER1 = {("A", "G"), ("T", "C")}
ENU_TIER2 = {("A", "T"), ("T", "A")}
ENU_TIER3 = {("A", "C"), ("T", "G"), ("G", "A"), ("C", "T")}
ENU_TIER4 = {("C", "A"), ("G", "T"), ("C", "G"), ("G", "C")}
ENU_ALL_TIERS = ENU_TIER1 | ENU_TIER2 | ENU_TIER3
ENU_ALL_TRANSITIONS = {("G", "A"), ("C", "T"), ("A", "G"), ("T", "C")}

# Strand-collapsed categories for plotting/summaries
STRAND_COLLAPSED = [
    ("A>G / T>C", {("A", "G"), ("T", "C")}),
    ("C>T / G>A", {("C", "T"), ("G", "A")}),
    ("A>T / T>A", {("A", "T"), ("T", "A")}),
    ("A>C / T>G", {("A", "C"), ("T", "G")}),
    ("C>A / G>T", {("C", "A"), ("G", "T")}),
    ("C>G / G>C", {("C", "G"), ("G", "C")}),
]

# RNA editing signatures: these transitions often occur due to deamination editing
RNA_EDITING_SIGNATURES = {
    ("A", "G"),  # A-to-I (+strand)
    ("T", "C"),  # A-to-I (-strand)
    ("C", "T"),  # C-to-U (+strand)
    ("G", "A"),  # C-to-U (-strand)
}


# ===================================================================
# DATA CLASS — VariantRecord
# ===================================================================
# Encapsulates a single variant's aggregated properties across groups and annotations.
# Using dataclasses improves readability and makes it straightforward to convert to dicts.
@dataclass
class VariantRecord:
    chrom: str
    pos: int
    ref: str
    alt: str
    qual: float
    filter_status: str

    # No-thrombosis (expected heterozygous: "suppressor carriers")
    no_thrombosis_genotype: str = "./."
    no_thrombosis_depth: int = 0
    no_thrombosis_ref_depth: int = 0
    no_thrombosis_alt_depth: int = 0
    no_thrombosis_vaf: float = 0.0
    no_thrombosis_het_count: int = 0
    no_thrombosis_total: int = 0

    # Thrombosis group (expected homozygous reference: lacks suppressor allele)
    thrombosis_genotype: str = "./."
    thrombosis_depth: int = 0
    thrombosis_ref_depth: int = 0
    thrombosis_alt_depth: int = 0
    thrombosis_vaf: float = 0.0
    thrombosis_homref_count: int = 0
    thrombosis_total: int = 0

    # Overall suppressor concordance metric (geometric mean between groups)
    suppressor_concordance: float = 0.0

    # Shared / VCF-level annotations (from INFO)
    qd: Optional[float] = None
    fs: Optional[float] = None
    mq: Optional[float] = None
    mq_rank_sum: Optional[float] = None
    read_pos_rank_sum: Optional[float] = None
    sor: Optional[float] = None
    is_snp: bool = False
    enu_tier: int = 0
    is_enu_transition: bool = False
    is_rna_editing_candidate: bool = False
    is_splice_junction_artifact: bool = False
    ase_skew: bool = False
    low_coverage_no_thrombosis: bool = False
    low_coverage_thrombosis: bool = False
    passes_gatk_hard_filters: bool = True
    meets_suppressor_genotype: bool = False
    priority_score: float = 0.0
    flags: list = field(default_factory=list)


# ===================================================================
# HELPERS — small utility functions
# ===================================================================

def safe_float(value):
    """
    Convert a VCF INFO value to float when possible.
    Returns None if the conversion fails, preserving downstream checks.
    Why: pysam returns tuples for some INFO fields or None when absent.
    """
    if value is None:
        return None
    try:
        return float(value)
    except (ValueError, TypeError):
        return None


def get_info(record, key):
    """
    Safe getter for INFO fields from a pysam.VariantRecord.
    If the INFO field is a tuple/list, return the first element.
    Return None if key missing.
    """
    try:
        val = record.info[key]
        if isinstance(val, tuple):
            return val[0]
        return val
    except (KeyError, TypeError):
        return None


def get_enu_tier(ref, alt):
    """
    Map a single-base substitution to an ENU tier (1..4) or 0 if not a recognized tier.
    NOTE: This is strand-aware only in the sense of explicit ref->alt,
    but many analyses collapse strands (A>G with T>C).
    """
    pair = (ref, alt)
    if pair in ENU_TIER1:
        return 1
    if pair in ENU_TIER2:
        return 2
    if pair in ENU_TIER3:
        return 3
    if pair in ENU_TIER4:
        return 4
    return 0


def genotype_string(gt_tuple):
    """
    Convert a GT tuple like (0,1) to string '0/1' or return './.' for missing.
    pysam uses None entries when genotype missing; preserve that formatting.
    """
    if gt_tuple is None:
        return "./."
    return "/".join(str(a) if a is not None else "." for a in gt_tuple)


def is_het(gt_str):
    """Return True if genotype string corresponds to a heterozygote."""
    return gt_str in ("0/1", "1/0")


def is_hom_ref(gt_str):
    """Return True if genotype string is homozygous reference."""
    return gt_str == "0/0"


# ===================================================================
# SAMPLE GROUP RESOLUTION — Validate samples are present in the VCF
# ===================================================================
def resolve_sample_groups(vcf_in):
    """
    Confirm the sample names specified above exist in the VCF header.
    Print status lines for each expected sample and abort with helpful diagnostics
    if there are missing or overlapping names.

    Arguments:
      vcf_in: pysam.VariantFile instance (already opened)

    Returns:
      tuple (no_thrombosis_samples_list, thrombosis_samples_list)
    """
    vcf_samples = set(vcf_in.header.samples)

    # --- Print the comparison being made (user-facing) ---
    print()
    print(f"  {'GROUP COMPARISON':=^60}")
    print(f"  No-thrombosis (expected HET):")
    for s in NO_THROMBOSIS_SAMPLES:
        status = "✓" if s in vcf_samples else "✗ NOT FOUND"
        print(f"    {s:<40s}  {status}")
    print(f"  Thrombosis (expected HOM-REF):")
    for s in THROMBOSIS_SAMPLES:
        status = "✓" if s in vcf_samples else "✗ NOT FOUND"
        print(f"    {s:<40s}  {status}")
    print(f"  Min group concordance:  {MIN_GROUP_CONCORDANCE:.0%} "
          f"({MIN_GROUP_CONCORDANCE:.0%} of replicates must agree)")
    print(f"  {'':=^60}")
    print()

    # --- Validate required conditions and produce helpful error messages ---
    missing_nt = [s for s in NO_THROMBOSIS_SAMPLES if s not in vcf_samples]
    missing_t = [s for s in THROMBOSIS_SAMPLES if s not in vcf_samples]
    overlap = set(NO_THROMBOSIS_SAMPLES) & set(THROMBOSIS_SAMPLES)

    errors = []
    if not NO_THROMBOSIS_SAMPLES:
        errors.append("NO_THROMBOSIS_SAMPLES is empty.")
    if not THROMBOSIS_SAMPLES:
        errors.append("THROMBOSIS_SAMPLES is empty.")
    if missing_nt:
        errors.append(f"No-thrombosis samples not found in VCF: {missing_nt}")
    if missing_t:
        errors.append(f"Thrombosis samples not found in VCF: {missing_t}")
    if overlap:
        errors.append(f"Samples appear in BOTH groups: {overlap}")

    if errors:
        # Print all errors and exit; better to fail early than to produce misleading results.
        print("  ERROR(S) in sample group configuration:")
        for e in errors:
            print(f"    • {e}")
        print(f"\n  Samples available in VCF: {sorted(vcf_samples)}")
        sys.exit(1)

    # --- Warn about unassigned samples in the VCF (not fatal) ---
    assigned = set(NO_THROMBOSIS_SAMPLES) | set(THROMBOSIS_SAMPLES)
    unassigned = vcf_samples - assigned
    if unassigned:
        print(f"  WARNING: {len(unassigned)} sample(s) in VCF not assigned "
              f"to any group (will be ignored): {sorted(unassigned)}")

    return list(NO_THROMBOSIS_SAMPLES), list(THROMBOSIS_SAMPLES)


# ===================================================================
# GATK HARD FILTERS — apply conservative thresholds
# ===================================================================
def passes_hard_filters(v):
    """
    Apply per-variant GATK-like thresholds (info is stored on VariantRecord v).
    This is a boolean sanity filter — it returns False if any metric indicates
    poor-quality evidence.

    Notes:
      - QD: quality-by-depth — low QD implies low confidence given depth
      - FS: FisherStrand — high FS implies strand bias
      - MQ: mapping quality
      - MQRankSum, ReadPosRankSum: rank-sum tests for read mapping / position biases
      - SOR: symmetric odds ratio for strand bias (another strand-bias metric)

    Careful: Not all VCFs contain all INFO fields. If an INFO field is None,
    we skip that check (i.e., treat missing as 'unknown', not an automatic fail).
    """
    if v.qd is not None and v.qd < QD_MIN:
        return False
    if v.fs is not None and v.fs > FS_MAX:
        return False
    if v.mq is not None and v.mq < MQ_MIN:
        return False
    if v.mq_rank_sum is not None and v.mq_rank_sum < MQRANKSUM_MIN:
        return False
    if v.read_pos_rank_sum is not None and v.read_pos_rank_sum < READPOSRANKSUM_MIN:
        return False
    if v.sor is not None and v.sor > SOR_MAX:
        return False
    return True


# ===================================================================
# PRIORITY SCORING — SUPPRESSOR LOGIC
# ===================================================================
def compute_priority(v):
    """
    Compute a 0..100 priority score for a VariantRecord v.

    The scoring is asymmetric and penalizes missing the suppressor genotype
    or poor quality features harshly (HARD KILL penalties).

    Rationale for weights:
     - Having the correct genotype pattern is essential (+30); not having it returns 0.
     - ENU tiering gives a strong boost because this project is an ENU screen.
     - Passing GATK filters is important (boost), but failing is only -10 (not total kill),
       because some good variants might miss an INFO metric due to caller differences.
     - Low coverage or failing genotype expectations is a HARD KILL (-100) keeping
       candidate list conservative.

    Notes:
      - This scoring is heuristic. Consider exposing weights as config options
        or computing them adaptively (e.g., via logistic regression trained on validated variants).
    """
    s = 0.0

    # --- CORE: suppressor genotype pattern must be met; otherwise the variant is irrelevant --
    if v.meets_suppressor_genotype:
        s += 30.0
    else:
        # HARD KILL: not meeting genotype pattern renders variant non-candidate
        return 0.0

    # Reward ENU-type substitutions for SNPs
    if v.is_snp:
        if v.enu_tier == 1:
            s += 25.0
        elif v.enu_tier == 2:
            s += 15.0
        elif v.enu_tier == 3:
            s += 10.0

    # GATK hard filter passing is weighted positively
    if v.passes_gatk_hard_filters:
        s += 15.0

    # Het VAF around 0.5 in no-thrombosis group is supportive of a real heterozygote
    if is_het(v.no_thrombosis_genotype) and 0.3 <= v.no_thrombosis_vaf <= 0.7:
        s += 10.0

    # Good depth in both groups
    if v.no_thrombosis_depth >= MIN_DEPTH and v.thrombosis_depth >= MIN_DEPTH:
        s += 10.0

    # PASS in VCF FILTER column gives a small boost (less reliable than INFO filters)
    if v.filter_status == "PASS":
        s += 5.0

    # High concordance across replicates gets a bonus
    if v.suppressor_concordance >= 0.9:
        s += 5.0

    # --- Penalties ---
    if v.is_rna_editing_candidate:
        s -= 25.0

    if v.is_splice_junction_artifact:
        s -= 20.0

    if v.ase_skew:
        s -= 15.0

    # HARD KILL: insufficient coverage in either group
    if v.low_coverage_no_thrombosis or v.low_coverage_thrombosis:
        s -= 100.0

    if not v.passes_gatk_hard_filters:
        s -= 10.0

    # Bound output to [0, 100] so it's simple to interpret
    return max(0.0, min(100.0, s))


# ===================================================================
# VCF PARSING — Extract and aggregate per-variant and per-group metrics
# ===================================================================

def extract_sample_fields(sample):
    """
    Extract DP, AD, GT, and compute VAF for a single pysam sample record.

    Parameters:
      sample: pysam.VariantRecordSample — behaves like a dict keyed by FORMAT tags

    Returns:
      gt_str: genotype string (e.g., '0/1', './.' if missing)
      dp: integer total depth (0 if missing)
      ref_depth: int (AD[0] or 0)
      alt_depth: int (AD[1] or 0)
      vaf: float alt_depth / (ref_depth + alt_depth) or 0.0 if total==0

    Caveats:
      - Some callers use separate fields (e.g., 'RO'/'AO' for samtools). This function
        assumes AD and DP are present. If your VCF uses other tags, add code to handle them.
      - AD may be None or missing; we convert gracefully to zeros.
    """
    # DP is sometimes missing; sample.get('DP', 0) is safest
    dp = sample.get("DP", 0) or 0

    # AD field is often a tuple of allele depths (ref, alt1, alt2, ...).
    # We only handle the first alt (bi-allelic assumption), which is correct for many screens.
    ad = sample.get("AD", (0, 0))
    if ad is None:
        ad = (0, 0)
    ref_depth = ad[0] if len(ad) > 0 else 0
    alt_depth = ad[1] if len(ad) > 1 else 0
    total = ref_depth + alt_depth

    # Avoid division by zero
    vaf = alt_depth / total if total > 0 else 0.0

    gt = sample.get("GT", (None, None))
    gt_str = genotype_string(gt) if gt else "./."
    return gt_str, dp, ref_depth, alt_depth, vaf


def aggregate_group(rec, sample_names):
    """
    Aggregate genotype & depth information across all samples in a group.

    Inputs:
      rec: pysam.VariantRecord for the variant
      sample_names: list of sample names to query (must exist in rec.samples)

    Returns:
      tuple: (consensus_gt, total_dp, total_ref, total_alt, mean_vaf,
              het_count, homref_count, n_samples)

    Notes:
      - total_dp is the sum of DP across the samples; this may double-count reads
        if the samples were pooled or derived from same library; the script assumes independent replicates.
      - mean_vaf is computed as sum_alt / (sum_ref + sum_alt) to reflect group-level allele balance.
      - consensus_gt is the most common genotype string across samples (majority call).
    """
    total_dp = 0
    total_ref = 0
    total_alt = 0
    het_count = 0
    homref_count = 0
    gt_strings = []

    for sname in sample_names:
        # Danger: rec.samples[sname] will raise KeyError if sname absent.
        # We checked sample existence earlier in resolve_sample_groups.
        sample = rec.samples[sname]
        gt_str, dp, rd, ad, vaf = extract_sample_fields(sample)
        total_dp += dp
        total_ref += rd
        total_alt += ad
        gt_strings.append(gt_str)
        if is_het(gt_str):
            het_count += 1
        if is_hom_ref(gt_str):
            homref_count += 1

    n = len(sample_names)
    total_allele_depth = total_ref + total_alt
    mean_vaf = total_alt / total_allele_depth if total_allele_depth > 0 else 0.0

    # Consensus genotype: the single most common genotype string; ties break arbitrarily.
    if gt_strings:
        gt_counter = Counter(gt_strings)
        consensus_gt = gt_counter.most_common(1)[0][0]
    else:
        consensus_gt = "./."

    return (consensus_gt, total_dp, total_ref, total_alt, mean_vaf,
            het_count, homref_count, n)


def parse_vcf():
    """
    Main VCF parsing loop.
    Iterates over every variant record and builds a VariantRecord instance per alt allele.

    Key operations per variant:
      - Identify SNPs vs. indels
      - Aggregate group-level genotype info for both groups
      - Apply genotype-based suppressor check (het in no-thromb, hom-ref in thromb)
      - Populate VariantRecord with INFO metrics and flags
      - Compute priority score

    Returns:
      list of VariantRecord objects
    """
    records = []
    # Open the VCF using pysam; this handles gzipped files natively
    vcf_in = pysam.VariantFile(INPUT_VCF)
    no_thromb_samples, thromb_samples = resolve_sample_groups(vcf_in)

    print(f"  No-thrombosis samples ({len(no_thromb_samples)}): {no_thromb_samples}")
    print(f"  Thrombosis samples ({len(thromb_samples)}): {thromb_samples}")

    # Iterate record-by-record; if the VCF is huge, consider limiting to PASS or high-Q variants first
    for rec in vcf_in:
        # For multiallelic sites, evaluate each ALT separately (this script treats each ALT as a separate candidate)
        for alt in rec.alts or []:
            # is_snp is a simple check for single-base ref and alt — multi-base or indels will be False
            is_snp = len(rec.ref) == 1 and len(alt) == 1

            # Aggregate metrics for both groups
            (nt_gt, nt_dp, nt_ref, nt_alt, nt_vaf,
             nt_het, nt_homref, nt_n) = aggregate_group(rec, no_thromb_samples)

            (t_gt, t_dp, t_ref, t_alt, t_vaf,
             t_het, t_homref, t_n) = aggregate_group(rec, thromb_samples)

            # --- Suppressor genotype check: we want het in no-thrombosis, hom-ref in thrombosis ---
            # Compute simple rates (#het / n) and (#hom-ref / n)
            nt_het_rate = nt_het / nt_n if nt_n > 0 else 0.0
            t_homref_rate = t_homref / t_n if t_n > 0 else 0.0

            # Concordance metric: geometric mean of the two rates; arbitrary but penalizes imbalance.
            # E.g., 1.0 if both are 1.0, 0.5 if one is 1.0 and the other 0.25 (sqrt(1 * 0.25)).
            concordance = (nt_het_rate * t_homref_rate) ** 0.5

            # Flag whether the variant meets the configured minimum concordance per group.
            meets_suppressor = (
                nt_het_rate >= MIN_GROUP_CONCORDANCE and
                t_homref_rate >= MIN_GROUP_CONCORDANCE
            )

            # Determine filter status field (VCF FILTER column). pysam stores rec.filter as an OrderedDict-like
            filters = rec.filter.keys()
            filter_str = ";".join(filters) if filters else "."

            # Instantiate VariantRecord — store fields we care about
            v = VariantRecord(
                chrom=rec.chrom,
                pos=rec.pos,
                ref=rec.ref,
                alt=alt,
                qual=rec.qual or 0.0,
                filter_status=filter_str,
                # no-thrombosis group
                no_thrombosis_genotype=nt_gt,
                no_thrombosis_depth=nt_dp,
                no_thrombosis_ref_depth=nt_ref,
                no_thrombosis_alt_depth=nt_alt,
                no_thrombosis_vaf=nt_vaf,
                no_thrombosis_het_count=nt_het,
                no_thrombosis_total=nt_n,
                # thrombosis group
                thrombosis_genotype=t_gt,
                thrombosis_depth=t_dp,
                thrombosis_ref_depth=t_ref,
                thrombosis_alt_depth=t_alt,
                thrombosis_vaf=t_vaf,
                thrombosis_homref_count=t_homref,
                thrombosis_total=t_n,
                # concordance
                suppressor_concordance=concordance,
                meets_suppressor_genotype=meets_suppressor,
                # info-level annotations
                qd=safe_float(get_info(rec, "QD")),
                fs=safe_float(get_info(rec, "FS")),
                mq=safe_float(get_info(rec, "MQ")),
                mq_rank_sum=safe_float(get_info(rec, "MQRankSum")),
                read_pos_rank_sum=safe_float(get_info(rec, "ReadPosRankSum")),
                sor=safe_float(get_info(rec, "SOR")),
                is_snp=is_snp,
            )

            # --- FLAGS: add human-friendly annotations stored in v.flags ---
            # 1. Low coverage per group
            # Use the aggregated DP — for pooled datasets this is appropriate; for per-embryo DP,
            # ensure your DP aggregation strategy is consistent with library prep.
            if nt_dp < MIN_DEPTH:
                v.low_coverage_no_thrombosis = True
                v.flags.append("LOW_COV_NO_THROMBOSIS")
            if t_dp < MIN_DEPTH:
                v.low_coverage_thrombosis = True
                v.flags.append("LOW_COV_THROMBOSIS")

            # 2. ENU tier — only makes sense for SNPs
            if is_snp:
                v.enu_tier = get_enu_tier(rec.ref, alt)
                if v.enu_tier == 1:
                    v.flags.append("ENU_TIER1_PRIMARY")
                elif v.enu_tier == 2:
                    v.flags.append("ENU_TIER2_STRONG")
                elif v.enu_tier == 3:
                    v.flags.append("ENU_TIER3_MODERATE")
                elif v.enu_tier == 4:
                    v.flags.append("ENU_TIER4_UNLIKELY")

                # is_enu_transition indicates a simple transition (Ti), used in Ti/Tv ratio calculations
                if (rec.ref, alt) in ENU_ALL_TRANSITIONS:
                    v.is_enu_transition = True

            # 3. RNA editing candidate: many A>G / C>T calls can be RNA editing rather than DNA mutation
            if is_snp and (rec.ref, alt) in RNA_EDITING_SIGNATURES:
                v.is_rna_editing_candidate = True
                v.flags.append("RNA_EDITING_CANDIDATE")

            # 4. ASE skew detection: heterozygotes in the no-thrombosis group with VAF far from 0.5
            if is_het(nt_gt):
                if nt_vaf < (0.5 - ASE_SKEW_THRESHOLD) or nt_vaf > (0.5 + ASE_SKEW_THRESHOLD):
                    v.ase_skew = True
                    v.flags.append("ASE_SKEW_NO_THROMBOSIS")

            # 5. Splice junction artifact detection: extreme ReadPosRankSum may indicate reads clipped at junctions
            if v.read_pos_rank_sum is not None and v.read_pos_rank_sum < SPLICE_RPRS_CUTOFF:
                v.is_splice_junction_artifact = True
                v.flags.append("SPLICE_JUNCTION_ARTIFACT")

            # 6. GATK-style hard filters (quality checks)
            v.passes_gatk_hard_filters = passes_hard_filters(v)
            if not v.passes_gatk_hard_filters:
                v.flags.append("FAILS_GATK_HARD_FILTER")

            # 7. Annotate whether the suppressor genotype pattern is met for this variant
            if meets_suppressor:
                v.flags.append("SUPPRESSOR_GENOTYPE")
            else:
                v.flags.append("NOT_SUPPRESSOR_GENOTYPE")

            # 8. Compute priority score given all flags & metrics
            v.priority_score = compute_priority(v)

            # Append to results
            records.append(v)

    vcf_in.close()
    return records


# ===================================================================
# BUILD DATAFRAME — Convert VariantRecord list to pandas DataFrame for summaries
# ===================================================================
def build_dataframe(records):
    """
    Convert a list of VariantRecord objects to a pandas DataFrame with
    flattened columns. This DataFrame is easier to save, filter, and plot.
    """
    rows = []
    for v in records:
        rows.append({
            "chrom": v.chrom,
            "pos": v.pos,
            "ref": v.ref,
            "alt": v.alt,
            "qual": v.qual,
            "filter": v.filter_status,
            # No-thrombosis group (suppressor carriers)
            "no_thromb_gt": v.no_thrombosis_genotype,
            "no_thromb_depth": v.no_thrombosis_depth,
            "no_thromb_ref_depth": v.no_thrombosis_ref_depth,
            "no_thromb_alt_depth": v.no_thrombosis_alt_depth,
            "no_thromb_vaf": round(v.no_thrombosis_vaf, 4),
            "no_thromb_het_count": v.no_thrombosis_het_count,
            "no_thromb_total": v.no_thrombosis_total,
            # Thrombosis group (affected)
            "thromb_gt": v.thrombosis_genotype,
            "thromb_depth": v.thrombosis_depth,
            "thromb_ref_depth": v.thrombosis_ref_depth,
            "thromb_alt_depth": v.thrombosis_alt_depth,
            "thromb_vaf": round(v.thrombosis_vaf, 4),
            "thromb_homref_count": v.thrombosis_homref_count,
            "thromb_total": v.thrombosis_total,
            # Suppressor concordance and flags
            "suppressor_concordance": round(v.suppressor_concordance, 4),
            "meets_suppressor_gt": v.meets_suppressor_genotype,
            # Annotations
            "QD": v.qd,
            "FS": v.fs,
            "MQ": v.mq,
            "MQRankSum": v.mq_rank_sum,
            "ReadPosRankSum": v.read_pos_rank_sum,
            "SOR": v.sor,
            "is_snp": v.is_snp,
            "enu_tier": v.enu_tier,
            "enu_transition": v.is_enu_transition,
            "rna_editing_candidate": v.is_rna_editing_candidate,
            "splice_artifact": v.is_splice_junction_artifact,
            "ase_skew": v.ase_skew,
            "low_cov_no_thromb": v.low_coverage_no_thrombosis,
            "low_cov_thromb": v.low_coverage_thrombosis,
            "passes_gatk_filters": v.passes_gatk_hard_filters,
            "priority_score": round(v.priority_score, 2),
            "flags": ";".join(v.flags) if v.flags else "NONE",
        })
    return pd.DataFrame(rows)


# ===================================================================
# CANDIDATE FILTER — SUPPRESSOR LOGIC
# ===================================================================
def get_candidates(df):
    """
    Apply hard filters to extract the suppressor candidate set.

    Conditions:
      1. Meets genotype pattern (het in no-thromb, hom-ref in thromb)
      2. Passes GATK hard filters
      3. Has sufficient depth in BOTH groups (no_thromb_depth >= MIN_DEPTH and thromb_depth >= MIN_DEPTH)
      4. Is a SNP
      5. Is ENU Tier 1..REQUIRE_ENU_TIER_MAX (inclusive)
      6. Has priority_score >= MIN_CANDIDATE_SCORE

    Returns:
      Filtered DataFrame sorted by priority_score descending.
    """
    mask = (
        (df["meets_suppressor_gt"]) &
        (df["passes_gatk_filters"]) &
        (df["no_thromb_depth"] >= MIN_DEPTH) &
        (df["thromb_depth"] >= MIN_DEPTH) &
        (df["priority_score"] >= MIN_CANDIDATE_SCORE) &
        (df["is_snp"]) &
        (df["enu_tier"] >= 1) &
        (df["enu_tier"] <= REQUIRE_ENU_TIER_MAX)
    )
    return df[mask].sort_values("priority_score", ascending=False)


# ===================================================================
# REPORT — Human-readable summary printed to stdout
# ===================================================================
def print_summary(df):
    """
    Print an extensive summary to the console that:
      - Documents the active filters/thresholds
      - Breaks down counts (SNPs, INDELs, filter pass/fail)
      - Shows suppressor-genotype distribution and VAF statistics
      - Displays ENU tier breakdown and mutation spectrum
      - Presents a filtering funnel and top candidates

    Notes:
      - Many computed percentages divide by 'total' which may be 0 if no variants found.
      - For reproducible pipelines consider writing this summary to a log file as well.
    """
    total = len(df)
    snps = df[df["is_snp"]]
    snp_total = len(snps)
    suppressor_gt = df[df["meets_suppressor_gt"]]

    print("=" * 72)
    print("  ENU MODIFIER / SUPPRESSOR SCREEN — RNA-seq VCF ANALYSIS")
    print("=" * 72)
    print()
    print("  SUPPRESSOR LOGIC:")
    print("    No-thrombosis group → HETEROZYGOUS  (carries suppressor)")
    print("    Thrombosis group    → HOM REFERENCE  (lacks suppressor)")
    print()

    # Active filters summary
    print(f"{'ACTIVE FILTERS':=^72}")
    print(f"  Min depth (per group):         {MIN_DEPTH:>8,}")
    print(f"  Min candidate score:           {MIN_CANDIDATE_SCORE:>8.1f}")
    print(f"  Min group concordance:         {MIN_GROUP_CONCORDANCE:>8.2f}")
    print(f"  Max ENU tier for candidates:   {REQUIRE_ENU_TIER_MAX:>8d}")
    print(f"  GATK: QD >= {QD_MIN}, FS <= {FS_MAX}, MQ >= {MQ_MIN}")
    print(f"  GATK: MQRankSum >= {MQRANKSUM_MIN}, ReadPosRankSum >= {READPOSRANKSUM_MIN}")
    print(f"  GATK: SOR <= {SOR_MAX}")
    print(f"  ASE skew threshold:            {ASE_SKEW_THRESHOLD:>8.2f}")
    print(f"  Splice artifact cutoff:        {SPLICE_RPRS_CUTOFF:>8.1f}")

    # Overall variant counts
    print(f"\n{'OVERALL VARIANT COUNTS':=^72}")
    print(f"  Total variants:                {total:>8,}")
    print(f"  SNPs:                          {len(snps):>8,}")
    print(f"  Indels:                        {len(df[~df['is_snp']]):>8,}")
    print(f"  Pass GATK hard filters:        {df['passes_gatk_filters'].sum():>8,}")
    print(f"  Fail GATK hard filters:        {(~df['passes_gatk_filters']).sum():>8,}")

    # Suppressor genotype distribution
    print(f"\n{'SUPPRESSOR GENOTYPE PATTERN':=^72}")
    if total > 0:
        print(f"  Meet suppressor genotype:      {len(suppressor_gt):>8,}  "
              f"({len(suppressor_gt) / total * 100:.1f}%)")
    else:
        print("  No variants to summarize.")

    print(f"  (Het in no-thromb AND hom-ref in thromb)")
    print()
    print(f"  --- No-thrombosis group genotypes ---")
    for gt, count in df["no_thromb_gt"].value_counts().items():
        pct = count / total * 100 if total > 0 else 0
        print(f"    {gt:<20s}  {count:>8,}  ({pct:5.1f}%)")
    print(f"\n  --- Thrombosis group genotypes ---")
    for gt, count in df["thromb_gt"].value_counts().items():
        pct = count / total * 100 if total > 0 else 0
        print(f"    {gt:<20s}  {count:>8,}  ({pct:5.1f}%)")

    # Coverage per group — basic stats
    print(f"\n{'COVERAGE BIAS (PER GROUP)':=^72}")
    for label, col_dp, col_low in [
        ("No-thrombosis", "no_thromb_depth", "low_cov_no_thromb"),
        ("Thrombosis", "thromb_depth", "low_cov_thromb"),
    ]:
        low_cov = df[col_low].sum()
        print(f"\n  {label} group:")
        if total > 0:
            print(f"    Depth < {MIN_DEPTH}:               {low_cov:>8,}  ({low_cov / total * 100:.1f}%)")
        else:
            print(f"    Depth < {MIN_DEPTH}:               {low_cov:>8,}")
        print(f"    Depth >= {MIN_DEPTH}:              {(~df[col_low]).sum():>8,}")
        # median/mean are nan-safe in pandas; display with fallback
        median = df[col_dp].median() if len(df) > 0 else 0
        mean = df[col_dp].mean() if len(df) > 0 else 0
        print(f"    Median depth:                {median:>8.0f}")
        print(f"    Mean depth:                  {mean:>8.1f}")

    # VAF comparison for suppressor genotype variants
    print(f"\n{'VAF COMPARISON BETWEEN GROUPS':=^72}")
    supp = suppressor_gt
    if len(supp) > 0:
        print(f"  Suppressor-genotype variants:  {len(supp):>8,}")
        print(f"\n  No-thrombosis group (expected ~0.5 het):")
        print(f"    Median VAF:                  {supp['no_thromb_vaf'].median():>8.4f}")
        print(f"    Mean VAF:                    {supp['no_thromb_vaf'].mean():>8.4f}")
        ase_ct = supp["ase_skew"].sum()
        print(f"    ASE skew sites:              {ase_ct:>8,}  "
              f"({ase_ct / len(supp) * 100:.1f}%)")
        print(f"\n  Thrombosis group (expected ~0.0 hom-ref):")
        print(f"    Median VAF:                  {supp['thromb_vaf'].median():>8.4f}")
        print(f"    Mean VAF:                    {supp['thromb_vaf'].mean():>8.4f}")
    else:
        print("  No variants meet the suppressor genotype pattern.")

    # ENU tiered signature for suppressor-genotype SNPs
    print(f"\n{'ENU MUTATION SIGNATURE (SUPPRESSOR GENOTYPE ONLY)':=^72}")
    supp_snps = supp[supp["is_snp"]] if len(supp) > 0 else pd.DataFrame()
    supp_snp_total = len(supp_snps)
    if supp_snp_total > 0:
        t1 = len(supp_snps[supp_snps["enu_tier"] == 1])
        t2 = len(supp_snps[supp_snps["enu_tier"] == 2])
        t3 = len(supp_snps[supp_snps["enu_tier"] == 3])
        t4 = len(supp_snps[supp_snps["enu_tier"] == 4])
        print(f"  Tier 1 — PRIMARY (A>G/T>C):    {t1:>8,}  ({t1 / supp_snp_total * 100:.1f}%)")
        print(f"  Tier 2 — STRONG  (A>T/T>A):    {t2:>8,}  ({t2 / supp_snp_total * 100:.1f}%)")
        print(f"  Tier 3 — MODERATE:             {t3:>8,}  ({t3 / supp_snp_total * 100:.1f}%)")
        print(f"  Tier 4 — UNLIKELY:             {t4:>8,}  ({t4 / supp_snp_total * 100:.1f}%)")
        enu_consistent = t1 + t2 + t3
        print(f"  ENU-consistent (Tier 1-3):     {enu_consistent:>8,}  ({enu_consistent / supp_snp_total * 100:.1f}%)")
        enu_trans = supp_snps["enu_transition"].sum()
        tv = supp_snp_total - enu_trans
        ti_tv = enu_trans / tv if tv > 0 else float("inf")
        print(f"  Ti/Tv ratio:                   {ti_tv:>8.2f}")
    else:
        print("  No suppressor-genotype SNPs found.")

    # Strand-collapsed mutation spectrum (all SNPs for context)
    if snp_total > 0:
        print(f"\n{'STRAND-COLLAPSED MUTATION SPECTRUM (ALL SNPs)':=^72}")
        spectrum_raw = Counter()
        for _, row in snps.iterrows():
            spectrum_raw[(row["ref"], row["alt"])] += 1

        for label, pairs in STRAND_COLLAPSED:
            count = sum(spectrum_raw.get(p, 0) for p in pairs)
            pct = count / snp_total * 100
            if pairs == ENU_TIER1:
                tag = " <-- ENU Tier 1 (PRIMARY)"
            elif pairs == ENU_TIER2:
                tag = " <-- ENU Tier 2 (STRONG)"
            elif label in ("C>T / G>A", "A>C / T>G"):
                tag = " <-- ENU Tier 3 (partial)"
            else:
                tag = ""
            print(f"  {label:<12s}  {count:>10,}  ({pct:5.1f}%){tag}")

        # Raw (per-strand) mutation spectrum for reference
        print(f"\n{'RAW MUTATION SPECTRUM (per-strand)':=^72}")
        for pair in sorted(spectrum_raw.keys()):
            count = spectrum_raw[pair]
            pct = count / snp_total * 100
            tier = get_enu_tier(pair[0], pair[1])
            tag = f" [Tier {tier}]" if tier > 0 else ""
            print(f"  {pair[0]}>{pair[1]}:  {count:>10,}  ({pct:5.1f}%){tag}")

    # RNA editing summary
    print(f"\n{'RNA EDITING CANDIDATE FLAGGING':=^72}")
    rna_edit = df["rna_editing_candidate"].sum()
    if total > 0:
        print(f"  RNA editing candidates:        {rna_edit:>8,}  ({rna_edit / total * 100:.1f}%)")
    else:
        print("  No variants to evaluate for RNA editing.")
    print(f"  WARNING: A>G / T>C overlaps ENU Tier 1!")
    print(f"  WARNING: G>A / C>T overlaps ENU Tier 3!")
    print(f"  Cross-reference with REDIportal / RADAR to resolve.")

    # Splice artifact summary
    print(f"\n{'SPLICE JUNCTION ARTIFACTS':=^72}")
    splice = df["splice_artifact"].sum()
    if total > 0:
        print(f"  Splice junction artifacts:     {splice:>8,}  ({splice / total * 100:.1f}%)")
    else:
        print("  No variants to evaluate for splice artifacts.")

    # Prioritization funnel — show how many variants survive each step
    print(f"\n{'SUPPRESSOR CANDIDATE PRIORITIZATION':=^72}")
    candidates = get_candidates(df)
    n_cand = len(candidates)

    supp_pass = df[df["meets_suppressor_gt"]]
    dp_pass = supp_pass[
        (supp_pass["no_thromb_depth"] >= MIN_DEPTH) &
        (supp_pass["thromb_depth"] >= MIN_DEPTH)
    ]
    gatk_pass = dp_pass[dp_pass["passes_gatk_filters"]]
    snp_pass = gatk_pass[gatk_pass["is_snp"]]
    tier_pass = snp_pass[
        (snp_pass["enu_tier"] >= 1) &
        (snp_pass["enu_tier"] <= REQUIRE_ENU_TIER_MAX)
    ]
    score_pass = tier_pass[tier_pass["priority_score"] >= MIN_CANDIDATE_SCORE]

    print(f"\n  FILTERING FUNNEL:")
    print(f"  Total variants:                {total:>8,}")
    print(f"  → suppressor genotype:         {len(supp_pass):>8,}")
    print(f"  → depth >= {MIN_DEPTH} (both):       {len(dp_pass):>8,}")
    print(f"  → pass GATK hard filters:      {len(gatk_pass):>8,}")
    print(f"  → SNPs only:                   {len(snp_pass):>8,}")
    print(f"  → ENU Tier 1-{REQUIRE_ENU_TIER_MAX} only:           {len(tier_pass):>8,}")
    print(f"  → score >= {MIN_CANDIDATE_SCORE}:              {len(score_pass):>8,}")
    print(f"  ═══════════════════════════════════════════")
    print(f"  FINAL SUPPRESSOR CANDIDATES:   {n_cand:>8,}")

    if n_cand > 0:
        high = candidates[candidates["priority_score"] >= 80]
        med = candidates[(candidates["priority_score"] >= 60) & (candidates["priority_score"] < 80)]
        low = candidates[candidates["priority_score"] < 60]
        print(f"\n  High priority (>= 80):         {len(high):>8,}")
        print(f"  Medium priority (60-79):       {len(med):>8,}")
        print(f"  Low priority (< 60):           {len(low):>8,}")
        print(f"  Mean priority score:           {candidates['priority_score'].mean():>8.1f}")

        print(f"\n  Candidates by ENU tier:")
        for tier in [1, 2]:
            t = candidates[candidates["enu_tier"] == tier]
            label = {1: "Tier 1 (PRIMARY)", 2: "Tier 2 (STRONG)"}[tier]
            print(f"    {label:<25s}  {len(t):>8,}")

    # Top candidates table
    print(f"\n{'TOP 20 SUPPRESSOR CANDIDATES':=^72}")
    print(f"  (het in no-thromb, hom-ref in thromb, ENU Tier 1-{REQUIRE_ENU_TIER_MAX},")
    print(f"   depth >= {MIN_DEPTH} both groups, score >= {MIN_CANDIDATE_SCORE})")
    top = candidates.head(20)
    cols = [
        "chrom", "pos", "ref", "alt",
        "no_thromb_gt", "no_thromb_vaf", "no_thromb_depth",
        "thromb_gt", "thromb_vaf", "thromb_depth",
        "enu_tier", "priority_score", "flags",
    ]
    if len(top) > 0:
        print(top[cols].to_string(index=False))
    else:
        print("  No variants meet all suppressor filter criteria.")

    print("\n" + "=" * 72)


# ===================================================================
# PLOTS — visual summaries saved as PNG (non-interactive)
# ===================================================================
def generate_plots(df):
    """
    Create a multi-panel figure summarizing:
     - VAF distributions for suppressor-genotype sites
     - Depth distributions for both groups
     - Strand-collapsed mutation spectrum colored by ENU tier
     - Priority score distribution among candidates
     - Filtering funnel bar chart
     - VAF scatter (no-thrombosis vs thrombosis) for suppressor-genotype sites

    The function uses matplotlib and saves the figure to PLOT_PATH.
    If df is empty, the function will generate placeholder panels.
    """
    snps = df[df["is_snp"]]
    suppressor = df[df["meets_suppressor_gt"]]
    candidates = get_candidates(df)

    fig, axes = plt.subplots(2, 3, figsize=(18, 10))
    fig.suptitle(
        f"ENU Modifier/Suppressor Screen — RNA-seq VCF Analysis\n"
        f"(Het in no-thromb, hom-ref in thromb, Tier 1-{REQUIRE_ENU_TIER_MAX}, "
        f"depth ≥ {MIN_DEPTH}, score ≥ {MIN_CANDIDATE_SCORE})",
        fontsize=13, fontweight="bold"
    )

    # Panel 1: VAF distribution for suppressor-genotype sites
    ax = axes[0, 0]
    supp_dp = suppressor[
        (suppressor["no_thromb_depth"] >= MIN_DEPTH) &
        (suppressor["thromb_depth"] >= MIN_DEPTH)
    ]
    if len(supp_dp) > 0:
        ax.hist(supp_dp["no_thromb_vaf"], bins=50, color="steelblue",
                edgecolor="black", alpha=0.7, label="No-thrombosis (het)")
        ax.hist(supp_dp["thromb_vaf"], bins=50, color="salmon",
                edgecolor="black", alpha=0.5, label="Thrombosis (hom-ref)")
        ax.axvline(0.5, color="red", linestyle="--", alpha=0.7, label="Expected het VAF")
        ax.axvline(0.0, color="blue", linestyle="--", alpha=0.4, label="Expected hom-ref VAF")
        # Shade ASE-skew regions
        ax.axvspan(0, 0.5 - ASE_SKEW_THRESHOLD, alpha=0.05, color="orange")
        ax.axvspan(0.5 + ASE_SKEW_THRESHOLD, 1.0, alpha=0.05, color="orange")
    ax.set_xlabel("Variant Allele Frequency")
    ax.set_ylabel("Count")
    ax.set_title("VAF: Suppressor Genotype Sites")
    ax.legend(fontsize=7)

    # Panel 2: Depth distribution for both groups (capped to avoid extreme outliers)
    ax = axes[0, 1]
    if len(df) > 0:
        cap_nt = df["no_thromb_depth"].quantile(0.99)
        cap_t = df["thromb_depth"].quantile(0.99)
        cap = max(cap_nt, cap_t)
        ax.hist(df["no_thromb_depth"].clip(upper=cap), bins=50, color="steelblue",
                edgecolor="black", alpha=0.6, label="No-thrombosis")
        ax.hist(df["thromb_depth"].clip(upper=cap), bins=50, color="salmon",
                edgecolor="black", alpha=0.5, label="Thrombosis")
        ax.axvline(MIN_DEPTH, color="red", linestyle="--", label=f"Min depth ({MIN_DEPTH})")
    ax.set_xlabel("Read Depth")
    ax.set_ylabel("Count")
    ax.set_title("Coverage Distribution (Both Groups)")
    ax.legend(fontsize=7)

    # Panel 3: Strand-collapsed mutation spectrum with tier coloring
    ax = axes[0, 2]
    if len(snps) > 0:
        spectrum_raw = Counter()
        for _, row in snps.iterrows():
            spectrum_raw[(row["ref"], row["alt"])] += 1
        labels = []
        counts = []
        colors = []
        tier_colors = {
            "Tier 1": "#e74c3c", "Tier 2": "#e67e22",
            "Tier 3": "#f1c40f", "Tier 4": "#3498db",
        }
        for label, pairs in STRAND_COLLAPSED:
            count = sum(spectrum_raw.get(p, 0) for p in pairs)
            labels.append(label)
            counts.append(count)
            if pairs == ENU_TIER1:
                colors.append(tier_colors["Tier 1"])
            elif pairs == ENU_TIER2:
                colors.append(tier_colors["Tier 2"])
            elif label in ("C>T / G>A", "A>C / T>G"):
                colors.append(tier_colors["Tier 3"])
            else:
                colors.append(tier_colors["Tier 4"])
        ax.bar(labels, counts, color=colors, edgecolor="black")
        from matplotlib.patches import Patch
        legend_items = [
            Patch(facecolor=tier_colors["Tier 1"], label="Tier 1 (PRIMARY)"),
            Patch(facecolor=tier_colors["Tier 2"], label="Tier 2 (STRONG)"),
            Patch(facecolor=tier_colors["Tier 3"], label="Tier 3 (MODERATE)"),
            Patch(facecolor=tier_colors["Tier 4"], label="Tier 4 (UNLIKELY)"),
        ]
        ax.legend(handles=legend_items, fontsize=7)
    ax.set_xlabel("Substitution Type (strand-collapsed)")
    ax.set_ylabel("Count")
    ax.set_title("ENU Tiered Mutation Spectrum")
    ax.tick_params(axis="x", rotation=45)

    # Panel 4: Priority score distribution (candidates)
    ax = axes[1, 0]
    if len(candidates) > 0:
        ax.hist(candidates["priority_score"], bins=30, color="darkorange",
                edgecolor="black", alpha=0.8)
    ax.axvline(MIN_CANDIDATE_SCORE, color="red", linestyle="--",
               label=f"Min score ({MIN_CANDIDATE_SCORE})")
    ax.set_xlabel("Priority Score")
    ax.set_ylabel("Count")
    ax.set_title("Priority Score (Suppressor Candidates)")
    ax.legend(fontsize=8)

    # Panel 5: Filtering funnel (log scale)
    ax = axes[1, 1]
    supp_pass = len(df[df["meets_suppressor_gt"]])
    dp_pass = len(df[
        (df["meets_suppressor_gt"]) &
        (df["no_thromb_depth"] >= MIN_DEPTH) &
        (df["thromb_depth"] >= MIN_DEPTH)
    ])
    gatk_pass = len(df[
        (df["meets_suppressor_gt"]) &
        (df["no_thromb_depth"] >= MIN_DEPTH) &
        (df["thromb_depth"] >= MIN_DEPTH) &
        (df["passes_gatk_filters"])
    ])
    snp_pass = len(df[
        (df["meets_suppressor_gt"]) &
        (df["no_thromb_depth"] >= MIN_DEPTH) &
        (df["thromb_depth"] >= MIN_DEPTH) &
        (df["passes_gatk_filters"]) &
        (df["is_snp"])
    ])
    tier_pass = len(df[
        (df["meets_suppressor_gt"]) &
        (df["no_thromb_depth"] >= MIN_DEPTH) &
        (df["thromb_depth"] >= MIN_DEPTH) &
        (df["passes_gatk_filters"]) &
        (df["is_snp"]) &
        (df["enu_tier"] >= 1) &
        (df["enu_tier"] <= REQUIRE_ENU_TIER_MAX)
    ])
    final = len(candidates)

    funnel_labels = [
        f"Total\n({len(df):,})",
        f"Suppressor\nGT\n({supp_pass:,})",
        f"Depth ≥\n{MIN_DEPTH}\n({dp_pass:,})",
        f"GATK\npass\n({gatk_pass:,})",
        f"SNPs\nonly\n({snp_pass:,})",
        f"Tier\n1-{REQUIRE_ENU_TIER_MAX}\n({tier_pass:,})",
        f"Score ≥\n{MIN_CANDIDATE_SCORE}\n({final:,})",
    ]
    funnel_vals = [len(df), supp_pass, dp_pass, gatk_pass, snp_pass, tier_pass, final]
    bar_colors = ["#95a5a6", "#1abc9c", "#3498db", "#2ecc71", "#e67e22", "#e74c3c", "#8e44ad"]
    ax.bar(range(len(funnel_vals)), funnel_vals, color=bar_colors, edgecolor="black")
    ax.set_xticks(range(len(funnel_labels)))
    ax.set_xticklabels(funnel_labels, fontsize=6)
    ax.set_ylabel("Variant Count")
    ax.set_title("Suppressor Filtering Funnel")
    ax.set_yscale("log")  # log scaling helps visualize wide ranges

    # Panel 6: VAF scatter (no-thrombosis vs thrombosis) for suppressor-gt variants
    ax = axes[1, 2]
    if len(supp_dp) > 0:
        tier_cmap = {1: "#e74c3c", 2: "#e67e22", 3: "#f1c40f", 0: "#95a5a6"}
        for tier in [0, 3, 2, 1]:  # draw tier 1 last so it appears on top
            subset = supp_dp[supp_dp["enu_tier"] == tier]
            if len(subset) > 0:
                label = {
                    1: "Tier 1 (PRIMARY)", 2: "Tier 2 (STRONG)",
                    3: "Tier 3 (MODERATE)", 0: "Other/Indel"
                }.get(tier, "Other")
                ax.scatter(
                    subset["thromb_vaf"],
                    subset["no_thromb_vaf"],
                    c=tier_cmap.get(tier, "#95a5a6"),
                    alpha=0.5, s=15, edgecolors="black", linewidths=0.3,
                    label=label,
                )
        ax.axhline(0.5, color="steelblue", linestyle="--", alpha=0.4, label="Het VAF = 0.5")
        ax.axvline(0.0, color="salmon", linestyle="--", alpha=0.4, label="Hom-ref VAF = 0")
        ax.set_xlabel("Thrombosis VAF (expected ~0)")
        ax.set_ylabel("No-Thrombosis VAF (expected ~0.5)")
        ax.set_title("Suppressor GT: VAF Group Comparison")
        ax.legend(fontsize=6)
    else:
        ax.text(0.5, 0.5, "No suppressor-GT variants\nwith sufficient depth",
                ha="center", va="center", transform=ax.transAxes)
        ax.set_title("Suppressor GT: VAF Group Comparison")

    plt.tight_layout(rect=[0, 0, 1, 0.90])
    plt.savefig(PLOT_PATH, dpi=150)
    plt.close()
    print(f"\n  Plots saved to: {PLOT_PATH}")


# ===================================================================
# MAIN ENTRYPOINT
# ===================================================================
def main():
    """
    Orchestrate the analysis:
      - Create output directory
      - Parse VCF and build DataFrame
      - Print summary
      - Save TSVs (all variants + candidates)
      - Generate plots
    """
    # Ensure output directory exists; use exist_ok so repeated runs don't crash
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    print("=" * 72)
    print("  ENU MODIFIER / SUPPRESSOR SCREEN ANALYZER")
    print("=" * 72)
    print(f"\n  Suppressor logic:")
    print(f"    No-thrombosis → HET (0/1)     = carries suppressor allele")
    print(f"    Thrombosis    → HOM-REF (0/0)  = lacks suppressor allele")
    print(f"\n  Input VCF:           {INPUT_VCF}")
    print(f"  Output dir:          {OUTPUT_DIR}")
    print(f"  Min depth filter:    {MIN_DEPTH} (100x * 40 embryos — adjust based on your experiment)")
    print(f"  Min candidate score: {MIN_CANDIDATE_SCORE}")
    print(f"  ENU tier max:        {REQUIRE_ENU_TIER_MAX}")
    print(f"\n  No-thrombosis group: {len(NO_THROMBOSIS_SAMPLES)} sample(s)")
    for s in NO_THROMBOSIS_SAMPLES:
        print(f"    → {s}")
    print(f"  Thrombosis group:    {len(THROMBOSIS_SAMPLES)} sample(s)")
    for s in THROMBOSIS_SAMPLES:
        print(f"    → {s}")
    print(f"  Group concordance:   {MIN_GROUP_CONCORDANCE:.0%}")

    # Parse VCF and compute per-variant metrics
    records = parse_vcf()

    if not records:
        print("  No variants found in VCF. Exiting.")
        sys.exit(0)

    print(f"  Parsed {len(records):,} variant records")

    # Convert to DataFrame for downstream filtering and plotting
    df = build_dataframe(records)

    # Print a comprehensive summary to stdout
    print_summary(df)

    # Save full table as TSV for downstream queries / record-keeping
    # Use tab separator — TSV is comfortable for large numeric columns and avoids CSV quoting issues with flags.
    df.to_csv(ALL_VARIANTS_TSV, sep="\t", index=False)
    print(f"\n  Full results saved to: {ALL_VARIANTS_TSV}")

    # Save suppressor candidates as TSV (this is the focused deliverable)
    candidates = get_candidates(df)
    candidates.to_csv(CANDIDATES_TSV, sep="\t", index=False)
    print(f"  {len(candidates):,} suppressor candidates saved to: {CANDIDATES_TSV}")

    # Generate QC + summary plots
    generate_plots(df)

    print("\n  Suppressor analysis complete!")


# Typical Python module guard so we can also import functions for unit testing
if __name__ == "__main__":
    main()

# ===================================================================
# POST-SCRIPT NOTES (not executed) — Additional suggestions
# ===================================================================
# - Reproducibility:
#   * Consider writing a small YAML configuration file for the sample groups and thresholds.
#   * Record the exact software versions used (pysam, pandas, numpy, matplotlib) to ensure
#     reproducibility. E.g. print their versions to the log.
#
# - Validation:
#   * Add a small unit-test script that constructs a tiny in-memory VCF (or uses a fixture)
#     and confirms parse_vcf(), aggregate_group(), and compute_priority() behave as intended.
#   * Test edge cases: missing AD, missing DP, multi-allelic sites, and GT fields like "0|1" (phased).
#
# - Performance:
#   * For very large VCFs (>10M variants), reading linearly is slow. Consider filtering upstream:
#     - Restrict to PASS variants and/or SNPs only using bcftools view -v snps -f PASS
#     - Use region-limited queries if you have candidate chromosomes
#   * Consider streaming and writing candidate rows directly rather than storing all VariantRecord objects
#     in memory if RAM is a constraint.
#
# - Improvements / Extensions:
#   * Add handling for multi-allelic AD arrays: compute per-ALT depths and treat each alt separately.
#   * Add optional cross-referencing with REDIportal or a local RNA-editing database to filter false positives.
#   * Use a simple logistic regression or rank-aggregation to combine features into a score that learns
#     from validated hits if you have labeled training data.
#   * Support VCFs that use RO/AO instead of AD (samtools legacy) by checking multiple tags in extract_sample_fields().
#
# - Logging:
#   * Replace print() with Python's logging module to get timestamps and adjustable verbosity.
#
# - Containerization:
#   * For portability, consider packaging this analysis in a Docker container with pinned package versions.
#
# End of annotations.
