#!/usr/bin/env python3
"""
results_summary_hardcoded_noargs_fix_patched_annotated.py

Heavily annotated version of the patched results summarizer.

Purpose
-------
Produce per-sample summary metrics combining BAM QC (flagstat / samtools/pysam),
variant-derived metrics (VAF distributions, genotype counts), and a mutational
spectrum (6-class and optional trinucleotide 96-context) derived from a VCF.

This annotated file keeps the original script logic and adds:
 - Extensive explanations for each function and block
 - Rationale for thresholds and fallbacks (AD -> GT)
 - Notes on common failure modes and troubleshooting
 - Suggestions for improvements (CLI, logging, unit tests, performance)
 - Guidance on reproducibility (versions, containerization)

Requirements
------------
- Python 3.8+
- cyvcf2 (for VCF parsing)
- pysam (optional, faster BAM parsing & reference access for contexts)
- samtools (on PATH) — used as fallback to compute some metrics
- featureCounts (optional) for gene-level counting
- matplotlib / seaborn / pandas / numpy (optional) for plotting

This script is intentionally "no-args" (hard-coded paths) to match the user's
existing workflow; annotate and then adapt to a CLI-based version if desired.
"""
# -------------------------
# Standard library imports
# -------------------------
from pathlib import Path
import shutil
import subprocess
import sys
import os
import json
import math
import re
import gzip
from statistics import mean, median, pstdev
from concurrent.futures import ProcessPoolExecutor, as_completed

# -------------------------
# Hard-coded project settings
# -------------------------
# These are intentionally embedded (mirrors user's environment). For production
# or containerized use, convert these to CLI args or a config file (YAML/JSON).
BAMROOT = Path("/home/jgd/Documents/bioinformatics_working/2025-11-06_enu-proc/align/star-ensembl")
VCF_PATH = Path("/home/jgd/Documents/bioinformatics_working/2025-11-06_enu-proc/vcf/gatk/gatk-ensembl/vep_annotated/gatk-ensembl_filtered.20260224.nodash.vep.vcf.gz")
OUT_DIR = Path("/home/jgd/Documents/bioinformatics_working/output")
SAMPLE_SHEET = Path("/home/jgd/Documents/bioinformatics_working/2025-11-06_enu-proc_bam/lists/2026-02-25_vcf-bam_samples.csv")
MAP_OUT = Path("/home/jgd/Documents/bioinformatics_working/2025-11-06_enu-proc_bam/lists/renaming_map_from_bams.tsv")
PLOTS_DIR = OUT_DIR / "plots"
REF_FA = Path("/home/jgd/Documents/bioinformatics_working/ref/danio_rerio.grcz11.bgz.fa.gz")  # recommended bgzipped and indexable
GTF = Path("/home/jgd/Documents/bioinformatics_working/ref/danio_rerio.grcz11.115.gtf.gz")

# SETTINGS: central place for toggles and thresholds. Edit as needed.
SETTINGS = {
    "compute_bam_metrics": True,
    "compute_mapq": True,
    "compute_depth": True,
    "compute_gene_counts": False,
    "gene_counts_stranded": 0,
    "compute_featurecounts": True,
    "make_plots": True,
    "jobs": 6,                 # number of parallel worker processes
    "threads": 12,             # threads to give featureCounts / similar
    "max_reads": None,         # limit reads scanned in BAM metrics (None = full scan)
    "max_variants": None,      # limit variants scanned in VCF loops (None = full)
    "max_reads_warn_threshold": 500000,
    "max_variants_warn_threshold": 100000,
    "min_depth": 10,           # minimum AD depth used to compute VAF
    "min_vaf": 0.0,            # minimum VAF to consider an alt allele present
    "mapq_threshold": 20
}

# -------------------------
# Optional heavy dependencies
# -------------------------
# We try importing cyvcf2 first (fast VCF parsing); pysam is optional but enables
# reference sequence access (trinucleotide contexts) and fast BAM iteration.
try:
    from cyvcf2 import VCF
except Exception:
    print("ERROR: cyvcf2 is required. Install with e.g. mamba install -c bioconda cyvcf2", file=sys.stderr)
    raise

try:
    import pysam
except Exception:
    pysam = None  # script will fall back to samtools for BAM introspection where possible

# plotting libs (optional)
try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import seaborn as sns
    import pandas as pd
    import numpy as np
except Exception:
    plt = None
    sns = None
    pd = None
    np = None

# -----------------------------------------------------------------------------
# Small helpers: shell command runner and robust numeric parsers
# -----------------------------------------------------------------------------
def run_cmd(cmd, capture_stdout=True):
    """
    Run a subprocess command and collect stdout/stderr.
    Returns (returncode, stdout_text, stderr_text).
    Note: cmd should be a list for safety.
    """
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE if capture_stdout else None,
                         stderr=subprocess.PIPE, text=True)
    out, err = p.communicate()
    return p.returncode, (out or ""), (err or "")

def int_from(s):
    """Parse an integer from string robustly (commas/whitespace allowed)."""
    if s is None:
        return None
    try:
        return int(re.sub(r"[,\s]", "", str(s)))
    except Exception:
        try:
            return int(float(s))
        except Exception:
            return None

def float_from(s):
    """Parse a float robustly from strings that may contain units or punctuation."""
    if s is None:
        return None
    try:
        return float(str(s))
    except Exception:
        try:
            return float(re.sub(r"[^\d\.\-eE]", "", str(s)))
        except Exception:
            return None

# -----------------------------------------------------------------------------
# featureCounts wrapper
# -----------------------------------------------------------------------------
def compute_gene_counts_featurecounts(bam_path, gtf_path, out_dir, sample_name, threads=1, stranded=0):
    """
    Run featureCounts to compute gene-level counts (optional).
    Returns (path_to_output, info_dict) or (None, {'skipped_reason': ...})
    - Verifies featureCounts on PATH
    - Writes output to out_dir/<sample>_featureCounts.txt
    - Parses the output to report the number of genes with counts > 0
    """
    outdir = Path(out_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    if not shutil.which("featureCounts"):
        return None, {"skipped_reason": "featureCounts not found on PATH"}
    safe_sample = re.sub(r"[^A-Za-z0-9._-]+", "_", str(sample_name))[:200]
    out_file = outdir / f"{safe_sample}_featureCounts.txt"
    # Build command
    cmd = ["featureCounts", "-T", str(threads), "-a", str(gtf_path), "-o", str(out_file), "-t", "exon", "-g", "gene_id"]
    try:
        cmd.extend(["-s", str(int(stranded))])
    except Exception:
        cmd.extend(["-s", "0"])
    cmd.append(str(bam_path))
    try:
        rc, out, err = run_cmd(cmd, capture_stdout=True)
    except Exception as e:
        return None, {"skipped_reason": f"featureCounts invocation failed: {e}"}
    if rc != 0:
        return None, {"skipped_reason": f"featureCounts failed (rc={rc}): {err.strip()}"}
    # Parse output to get how many genes had non-zero counts (basic QC)
    genes_with_counts = 0
    try:
        with open(out_file, "r", encoding="utf-8", errors="replace") as fh:
            header_seen = False
            for ln in fh:
                if ln.startswith("#") or ln.strip() == "":
                    continue
                if not header_seen and (ln.startswith("Geneid") or ln.lower().startswith("geneid")):
                    header_seen = True
                    continue
                if not header_seen:
                    continue
                parts = ln.rstrip("\n").split("\t")
                try:
                    cnt = int(parts[-1])
                except Exception:
                    cnt = 0
                if cnt > 0:
                    genes_with_counts += 1
    except Exception as e:
        return out_file, {"skipped_reason": f"failed to parse featureCounts output: {e}"}
    return out_file, {"genes_with_counts": genes_with_counts}

# -----------------------------------------------------------------------------
# Samtools flagstat wrapper: parse output lines into a small metrics dict
# -----------------------------------------------------------------------------
def run_flagstat(bam_path):
    """
    Run `samtools flagstat` and parse a few common metrics.
    - Returns dict containing parsed integers and the raw output under 'raw'
    - If samtools missing, raises RuntimeError
    Notes:
    - Parsing relies on typical flagstat English wording; localized samtools may differ.
    """
    if not shutil.which("samtools"):
        raise RuntimeError("samtools not found on PATH.")
    try:
        out = subprocess.check_output(["samtools", "flagstat", str(bam_path)], stderr=subprocess.STDOUT, text=True)
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"samtools flagstat failed: {e.output}")
    metrics = {"raw": out, "total": None, "mapped": None, "mapped_pct": None,
               "primary": None, "secondary": None, "duplicates": None,
               "paired": None, "properly_paired": None, "properly_paired_pct": None,
               "read1": None, "read2": None}
    for ln in out.splitlines():
        ln = ln.strip()
        # The following heuristics extract numbers from typical lines like:
        # "12345 + 0 in total (QC-passed reads + QC-failed reads)"
        if "in total" in ln and metrics["total"] is None:
            metrics["total"] = int_from(ln.split()[0])
        if "mapped (" in ln and metrics["mapped"] is None:
            metrics["mapped"] = int_from(ln.split()[0])
            m = re.search(r"mapped \(([\d\.]+)%", ln)
            if m:
                metrics["mapped_pct"] = float_from(m.group(1))
        if ln.endswith("primary") and metrics["primary"] is None:
            metrics["primary"] = int_from(ln.split()[0])
        if ln.endswith("secondary") and metrics["secondary"] is None:
            metrics["secondary"] = int_from(ln.split()[0])
        if ln.endswith("duplicates") and metrics["duplicates"] is None:
            metrics["duplicates"] = int_from(ln.split()[0])
        if "paired in sequencing" in ln and metrics["paired"] is None:
            metrics["paired"] = int_from(ln.split()[0])
        if "properly paired" in ln and metrics["properly_paired"] is None:
            metrics["properly_paired"] = int_from(ln.split()[0])
            m = re.search(r"properly paired \(([\d\.]+)%", ln)
            if m:
                metrics["properly_paired_pct"] = float_from(m.group(1))
        if "read1" in ln and metrics["read1"] is None:
            metrics["read1"] = int_from(ln.split()[0])
        if "read2" in ln and metrics["read2"] is None:
            metrics["read2"] = int_from(ln.split()[0])
    return metrics

# -----------------------------------------------------------------------------
# Mutational spectrum helpers
# -----------------------------------------------------------------------------
def canonicalize_substitution_and_context(ref_base, alt_base, tri):
    """
    Convert the substitution and trinucleotide context to pyrimidine-centered form:
    - For C/T reference bases, keep as-is (C>X or T>X).
    - For purine refs (A/G), complement the triplet and map substitution to complementary
      C/T space. This produces canonical classes used in many mutational signatures.
    Returns (substitution_str, context_str) where:
      substitution_str e.g. "C>T"
      context_str e.g. "A[C>T]G"
    """
    comp = {"A":"T","C":"G","G":"C","T":"A","N":"N"}
    ref_base = (ref_base or "N").upper()
    alt_base = (alt_base or "N").upper()
    tri = (tri or "").upper()
    if len(tri) != 3:
        tri = tri.ljust(3, "N")[:3]
    # If reference is a pyrimidine, we keep orientation
    if ref_base in ("C","T"):
        sub = f"{ref_base}>{alt_base}"
        context = f"{tri[0]}[{ref_base}>{alt_base}]{tri[2]}"
        return sub, context
    else:
        # reverse complement the tri-nucleotide and substitution
        tri_comp = "".join(comp.get(b,"N") for b in reversed(tri))
        ref_c = comp.get(ref_base,"N")
        alt_c = comp.get(alt_base,"N")
        sub = f"{ref_c}>{alt_c}"
        context = f"{tri_comp[0]}[{ref_c}>{alt_c}]{tri_comp[2]}"
        return sub, context

def compute_mutational_spectrum(vcf_path, sample_name, ref_fasta_path=None, max_variants=None, min_depth=10, min_vaf=0.0):
    """
    Compute per-sample mutational spectrum.
    - ALWAYS returns six-class spectrum (C>A,C>G,C>T,T>A,T>C,T>G) as dict 'six'.
    - Prefers allele depth (AD) values to decide whether the sample supports an ALT.
      If AD is missing, falls back to GENOTYPE (GT) presence.
    - If a reference fasta (and pysam available) is provided, computes trinucleotide
      contexts and returns a dictionary 'context96' mapping e.g. 'A[C>T]G' -> count.
    - Parameters:
        vcf_path: path to VCF (cyvcf2-compatible)
        sample_name: sample ID in VCF.samples
        max_variants: optional limit to number of variants scanned
        min_depth: minimum read depth (AD ref+alt) to consider site
        min_vaf: minimum variant allele frequency (alt/depth) to consider alt present
    - Returns:
        dict with keys: skipped (bool), reason (if skipped), six (dict), context96 (dict),
                       variants_examined (int)
    """
    v = VCF(str(vcf_path))
    if sample_name not in v.samples:
        v.close()
        return {"skipped": True, "reason": f"sample {sample_name} not in VCF"}
    idx = v.samples.index(sample_name)

    # Initialize 6-class counters
    six = {"C>A":0,"C>G":0,"C>T":0,"T>A":0,"T>C":0,"T>G":0}
    context96 = {}

    # Prepare reference if requested and pysam available
    ref = None
    use_context = False
    if ref_fasta_path and pysam is not None:
        try:
            ref = pysam.FastaFile(str(ref_fasta_path))
            use_context = True
        except Exception:
            use_context = False  # gracefully fall back; context96 will be empty

    examined = 0
    # Iterate variants; cyvcf2 yields Variant objects
    for variant in v:
        # We only consider single-nucleotide REF and single-nucleotide ALT alleles here.
        if variant.REF is None or len(variant.REF) != 1:
            continue
        alts = variant.ALT
        if not alts:
            continue

        # Attempt to read AD FORMAT once per variant for all samples
        try:
            ad_field = variant.format("AD")
        except Exception:
            ad_field = None

        for alt_index, alt in enumerate(alts):
            # skip non-SNV alt alleles
            if alt is None or len(alt) != 1:
                continue

            sample_has_alt = False

            # Prefer AD if present
            if ad_field is not None:
                try:
                    sample_ad = ad_field[idx]
                except Exception:
                    sample_ad = None
                if sample_ad is not None and len(sample_ad) >= 2:
                    # AD often contains ref followed by per-alt counts. For multi-allelics,
                    # alt counts can be AD[1], AD[2], ...
                    try:
                        ref_count = int(sample_ad[0])
                    except Exception:
                        ref_count = 0
                    try:
                        alt_count = int(sample_ad[1 + alt_index]) if (1 + alt_index) < len(sample_ad) else 0
                    except Exception:
                        alt_count = 0
                    depth = ref_count + alt_count
                    # require minimal evidence (depth and alt count/vaf)
                    if depth >= min_depth and alt_count > 0:
                        vaf = (alt_count / depth) if depth > 0 else 0.0
                        if vaf >= min_vaf:
                            sample_has_alt = True

            # If AD not available or insufficient, fallback to GT presence
            if not sample_has_alt:
                try:
                    gt = variant.genotypes[idx]
                except Exception:
                    gt = None
                if gt:
                    # cyvcf2 genotype list: [allele1, allele2, phased_flag?]
                    # handle haploid and diploid data robustly
                    a1 = gt[0] if len(gt) > 0 else None
                    a2 = gt[1] if len(gt) > 1 else None
                    alleles = [a for a in (a1, a2) if a is not None]
                    # skip missing allele encodings (-1)
                    if alleles and not any(a == -1 for a in alleles):
                        # In VCF, alt alleles are indexed 1..n; check if any allele equals this alt_index+1
                        if any(a == (alt_index + 1) for a in alleles):
                            sample_has_alt = True

            if not sample_has_alt:
                # sample does not carry this alt (by AD or GT fallback)
                continue

            # At this point the sample supports the alt (count it)
            # Canonicalize into pyrimidine-centered substitution and optionally gather context
            sub = None
            ctx = None
            if use_context:
                try:
                    pos0 = variant.start  # cyvcf2 uses 0-based start
                    chrom = variant.CHROM
                    # fetch tri-nucleotide centered on variant position
                    tri = ref.fetch(chrom, max(0, pos0 - 1), pos0 + 2)
                    if len(tri) == 3:
                        sub, ctx = canonicalize_substitution_and_context(variant.REF, alt, tri)
                        context96[ctx] = context96.get(ctx, 0) + 1
                except Exception:
                    # if anything fails (missing contig, fetch error), we silently fall back
                    sub = None

            # If we couldn't compute context-based canonicalization, do a simple mapping
            if sub is None:
                rb = (variant.REF or "N").upper()
                ab = (alt or "N").upper()
                if rb in ("C", "T"):
                    sub = f"{rb}>{ab}"
                else:
                    comp = {"A":"T","C":"G","G":"C","T":"A","N":"N"}
                    sub = f"{comp.get(rb,'N')}>{comp.get(ab,'N')}"

            if sub in six:
                six[sub] += 1

        examined += 1
        if max_variants and examined >= max_variants:
            break

    v.close()
    return {"skipped": False, "six": six, "context96": context96, "variants_examined": examined}

# -----------------------------------------------------------------------------
# Compute VAFs per sample (using AD when present)
# -----------------------------------------------------------------------------
def compute_sample_vafs(vcf_path, sample_name, min_depth=10, max_variants=None):
    """
    Walk through variants and collect per-site VAFs for a given sample.
    - Uses AD format if present; gracefully handles missing/malformed AD.
    - Returns (vaf_list (floats), counters dict) where counters include:
      total_seen, ad_missing, depth_filtered, vaf_count, ad_nonzero
    Notes:
    - The function aggregates multiple alt alleles' counts into alt (sum of alt AD fields)
      which is sensible for per-site VAF summaries in the presence of multiple ALTs.
    """
    v = VCF(str(vcf_path))
    if sample_name not in v.samples:
        v.close()
        raise KeyError(f"Sample '{sample_name}' not present in VCF samples")
    idx = v.samples.index(sample_name)
    vafs = []
    counters = {"total_seen": 0, "ad_missing": 0, "depth_filtered": 0, "vaf_count": 0, "ad_nonzero": 0}
    for variant in v:
        counters["total_seen"] += 1
        try:
            ad = variant.format("AD")
        except Exception:
            ad = None
        if ad is None:
            counters["ad_missing"] += 1
            continue
        try:
            sample_ad = ad[idx]
        except Exception:
            counters["ad_missing"] += 1
            continue
        if sample_ad is None or len(sample_ad) < 2:
            counters["ad_missing"] += 1
            continue
        # AD: first element = ref, following = per-alt counts
        try:
            ref = int(sample_ad[0])
            alt = sum(int(x) for x in sample_ad[1:])
        except Exception:
            counters["ad_missing"] += 1
            continue
        # count sites where alt coverage is non-zero (useful diagnostics)
        if alt > 0:
            counters["ad_nonzero"] += 1
        depth = ref + alt
        if depth < min_depth:
            counters["depth_filtered"] += 1
            continue
        if alt <= 0:
            continue
        vaf = alt / depth
        vafs.append(vaf)
        counters["vaf_count"] += 1
        if max_variants and counters["vaf_count"] >= max_variants:
            break
    v.close()
    return vafs, counters

def summarize_vafs(vafs):
    """
    Compute some basic summary statistics for a list of VAF floats.
    Returns dict with count, mean, median, min/max, stddev, quartiles, proportion >= 0.5
    """
    if not vafs:
        return {"count": 0, "mean": None, "median": None, "min": None, "max": None,
                "stddev": None, "q25": None, "q75": None, "prop_ge_0.5": None}
    s = sorted(vafs)
    cnt = len(s)
    q25 = s[int(0.25 * (cnt - 1))]
    q75 = s[int(0.75 * (cnt - 1))]
    return {"count": cnt, "mean": mean(s), "median": median(s), "min": s[0], "max": s[-1],
            "stddev": pstdev(s), "q25": q25, "q75": q75, "prop_ge_0.5": sum(1 for x in s if x >= 0.5) / cnt}

# -----------------------------------------------------------------------------
# Genotype summary using GT field (and robust handling of haploid/partial data)
# -----------------------------------------------------------------------------
def compute_genotype_summary(vcf_path, sample_name, max_variants=None):
    """
    Produce counts for genotype states across VCF sites:
    - hom_ref, het, hom_alt, missing, other, total
    Also computes nonref_count (any site where GT indicates at least one alt allele).
    Uses cyvcf2 Variant.genotypes for fast access.
    """
    v = VCF(str(vcf_path))
    if sample_name not in v.samples:
        v.close()
        raise KeyError(f"Sample '{sample_name}' not present in VCF samples")
    idx = v.samples.index(sample_name)
    counts = {"hom_ref": 0, "het": 0, "hom_alt": 0, "missing": 0, "other": 0, "total": 0, "nonref_count": 0}
    for variant in v:
        if max_variants and counts["total"] >= max_variants:
            break
        counts["total"] += 1
        try:
            gt = variant.genotypes[idx]
        except Exception:
            counts["missing"] += 1
            continue
        if not gt:
            counts["missing"] += 1
            continue
        # GT is typically a small list: [a1, a2, phased_flag?]
        a1 = gt[0] if len(gt) > 0 else None
        a2 = gt[1] if len(gt) > 1 else None
        alleles = [a for a in (a1, a2) if a is not None]
        if not alleles:
            counts["missing"] += 1
            continue
        # cyvcf2 encodes missing allele as -1
        if any(a == -1 for a in alleles):
            counts["missing"] += 1
            continue
        # nonref_count increments if any allele > 0 (ALT present)
        if any(a > 0 for a in alleles):
            counts["nonref_count"] += 1
        # classify genotype
        if all(a == 0 for a in alleles):
            counts["hom_ref"] += 1
        elif len(alleles) >= 2 and alleles[0] == alleles[1] and alleles[0] > 0:
            counts["hom_alt"] += 1
        elif (len(alleles) >= 2 and ((alleles[0] == 0 and alleles[1] > 0) or (alleles[1] == 0 and alleles[0] > 0))) or (len(alleles) >= 2 and alleles[0] != alleles[1]):
            # captures simple heterozygotes and multi-allelic differences
            counts["het"] += 1
        elif any(a > 0 for a in alleles):
            # fallback classification: treat as het if any alt present but not matching above
            counts["het"] += 1
        else:
            counts["other"] += 1
    v.close()
    return counts

# -----------------------------------------------------------------------------
# Coverage-based metrics and BAM scanning: using pysam if available; else samtools fallback.
# compute_bam_metrics returns reads_examined in both code paths so callers can report it.
# -----------------------------------------------------------------------------
def compute_coverage_metrics(bam_path, depth_thresholds=(10,20)):
    """
    Compute per-base coverage stats using samtools depth -aa.
    Returns mean depth and percentage of bases >= each threshold.
    Note: This can be slow for whole-genome BAMs; consider region-limited calls or using
    depth summary tools (mosdepth).
    """
    if not shutil.which("samtools"):
        raise RuntimeError("samtools not found on PATH.")
    cmd = ["samtools", "depth", "-aa", str(bam_path)]
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    total_bases = 0
    total_depth = 0
    thresh_counts = {t: 0 for t in depth_thresholds}
    for line in p.stdout:
        parts = line.split()
        if len(parts) < 3:
            continue
        try:
            d = int(parts[2])
        except Exception:
            continue
        total_bases += 1
        total_depth += d
        for t in depth_thresholds:
            if d >= t:
                thresh_counts[t] += 1
    p.stdout.close()
    p.wait()
    if total_bases == 0:
        return {"mean_depth": None, **{f"pct_ge_{t}": None for t in depth_thresholds}}
    mean_depth = total_depth / total_bases
    result = {"mean_depth": mean_depth}
    for t in depth_thresholds:
        result[f"pct_ge_{t}"] = (thresh_counts[t] / total_bases) * 100.0
    return result

def compute_bam_metrics(bam_path, max_reads=None, mapq_threshold=20):
    """
    Compute a few summary metrics scanning the BAM:
    - Insert size distribution (if paired)
    - Fraction reads below MAPQ threshold
    - Fraction of soft-clipped reads
    - Total reads examined
    Implementation:
      - Prefer pysam (faster, more robust). If pysam unavailable, fall back to `samtools view`.
      - Guarantees 'reads_examined' in returned dict for downstream consistency.
    """
    isizes = []
    low_mapq = 0
    total_mapq = 0
    softclip = 0
    total_reads = 0

    def _record_isize(v):
        try:
            v = abs(int(v))
            return v if v > 0 else None
        except Exception:
            return None

    # First try pysam if available
    if pysam is not None:
        try:
            with pysam.AlignmentFile(str(bam_path), "rb") as fh:
                # fetch(until_eof=True) reads all reads; be careful with huge BAMs
                for read in fh.fetch(until_eof=True):
                    if read.is_secondary or read.is_supplementary or read.is_unmapped:
                        continue
                    total_reads += 1
                    mq = read.mapping_quality if read.mapping_quality is not None else 0
                    if mq < mapq_threshold:
                        low_mapq += 1
                    total_mapq += 1
                    cig = read.cigartuples
                    if cig and any(op == 4 for op, l in cig):
                        softclip += 1
                    if read.is_paired and read.is_proper_pair:
                        tlen = abs(read.template_length)
                        if tlen and tlen > 0:
                            isizes.append(tlen)
                    if max_reads and total_reads >= max_reads:
                        break
        except Exception as e:
            # Log an error and fall back to samtools path below
            sys.stderr.write(f"pysam scan failed on {bam_path}: {e}\n")

    # If pysam is not available or didn't yield reads, try samtools view
    if pysam is None or (pysam is not None and total_reads == 0):
        if not shutil.which("samtools"):
            return {"insert_mean": None, "insert_median": None, "insert_std": None, "insert_count": 0,
                    "low_mapq_frac": None, "low_mapq_count": 0, "total_mapq_count": 0,
                    "pct_softclipped": None, "softclip_count": 0, "reads_examined": 0}
        # We'll exclude secondary alignments (0x100) and supplementary (0x800) using -F
        cmd = ["samtools", "view", "-F", "256", str(bam_path)]
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        for line in p.stdout:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 6:
                continue
            flag = int(parts[1])
            # filter out secondary (0x100), supplementary (0x800), unmapped (0x4)
            if flag & 0x100 or flag & 0x800 or flag & 0x4:
                continue
            total_reads += 1
            try:
                mq = int(parts[4])
            except Exception:
                mq = 0
            if mq < mapq_threshold:
                low_mapq += 1
            total_mapq += 1
            cigar = parts[5]
            if "S" in cigar:
                softclip += 1
            if len(parts) >= 9:
                tlen = _record_isize(parts[8])
                if tlen:
                    isizes.append(tlen)
            if max_reads and total_reads >= max_reads:
                break
        p.stdout.close()
        p.wait()

    # Compute insert size stats if any were observed
    insert_count = len(isizes)
    insert_mean = insert_median = insert_std = None
    if insert_count:
        insert_mean = float(sum(isizes)) / insert_count
        insert_median = float(median(isizes))
        try:
            insert_std = float(pstdev(isizes))
        except Exception:
            insert_std = None

    low_mapq_frac = (low_mapq / total_mapq) if total_mapq and total_mapq > 0 else None
    pct_softclipped = (softclip / total_reads) if total_reads and total_reads > 0 else None
    # Return robustly-typed dict; include reads_examined for callers
    return {"insert_mean": insert_mean, "insert_median": insert_median, "insert_std": insert_std, "insert_count": insert_count,
            "low_mapq_frac": low_mapq_frac, "low_mapq_count": low_mapq, "total_mapq_count": total_mapq,
            "pct_softclipped": pct_softclipped, "softclip_count": softclip,
            "reads_examined": total_reads}

# -----------------------------------------------------------------------------
# Plotting helpers (robust to missing plotting libs)
# -----------------------------------------------------------------------------
def generate_sample_mutation_plots(records, plots_dir: Path, dpi=150):
    """
    Create per-sample 6-class barplots and an aggregate plot.
    - Accepts a list of record dicts produced by process_sample_task.
    - If plotting libraries are missing (matplotlib/seaborn/pandas), skip gracefully.
    """
    if plt is None or sns is None:
        print("Plotting libs missing; skipping per-sample mutation plots.", file=sys.stderr)
        return []
    KEYS = ["C>A","C>G","C>T","T>A","T>C","T>G"]
    out = []
    ms_dir = Path(plots_dir) / "mut_spectra"
    ms_dir.mkdir(parents=True, exist_ok=True)

    muts = []
    # Extract per-sample six-class dictionaries (be permissive with JSON/string encodings)
    for rec in records:
        sample = str(rec.get("sample") or "")
        if not sample:
            continue
        six = rec.get("mut_six")
        if isinstance(six, str):
            try:
                import ast
                six = ast.literal_eval(six)
            except Exception:
                try:
                    import json as _json
                    six = _json.loads(six.replace("'", '"'))
                except Exception:
                    six = {}
        # If six-class is empty, attempt to derive from context96 if present
        if not isinstance(six, dict) or not any(int(six.get(k, 0) or 0) for k in six):
            ctx = rec.get("mut_context96")
            if isinstance(ctx, str):
                try:
                    import ast
                    ctx = ast.literal_eval(ctx)
                except Exception:
                    try:
                        import json as _json
                        ctx = _json.loads(ctx.replace("'", '"'))
                    except Exception:
                        ctx = {}
            if isinstance(ctx, dict) and ctx:
                derived = {k: 0 for k in KEYS}
                for k, v in ctx.items():
                    try:
                        cnt = int(v)
                    except Exception:
                        try:
                            cnt = int(float(v))
                        except Exception:
                            cnt = 0
                    # extract substitution string from context key (e.g., A[C>T]G)
                    m = re.search(r"\[([ACGT]>[ACGT])\]", str(k))
                    if not m:
                        m2 = re.search(r"([ACGT]>[ACGT])", str(k))
                        sub = m2.group(1) if m2 else None
                    else:
                        sub = m.group(1)
                    if sub in derived:
                        derived[sub] += cnt
                six = derived
            else:
                six = {k: 0 for k in KEYS}
        clean = {k: int(six.get(k, 0) or 0) for k in KEYS}
        muts.append((sample, clean))

        vals = [clean[k] for k in KEYS]
        fig, ax = plt.subplots(figsize=(6,4))
        try:
            sns.barplot(x=KEYS, y=vals, palette="Set2", ax=ax)
        except Exception:
            ax.bar(KEYS, vals, color="C0")
        ax.set_ylabel("Count")
        ax.set_title(f"{sample}: 6-class mutational spectrum")
        plt.tight_layout()
        fname = re.sub(r"[^A-Za-z0-9._-]+", "_", sample)[:220]
        p = ms_dir / f"{fname}_mutational_spectrum_6class.png"
        fig.savefig(p, dpi=dpi)
        plt.close(fig)
        out.append(p)

    # Aggregate plot
    if not muts:
        return out

    agg = {k: 0 for k in KEYS}
    for _, d in muts:
        for k in KEYS:
            agg[k] += int(d.get(k, 0) or 0)
    fig, ax = plt.subplots(figsize=(6,4))
    vals = [agg[k] for k in KEYS]
    try:
        sns.barplot(x=KEYS, y=vals, palette="Set2", ax=ax)
    except Exception:
        ax.bar(KEYS, vals, color="C1")
    ax.set_ylabel("Count (aggregate)")
    ax.set_title("Aggregate 6-class mutational spectrum")
    plt.tight_layout()
    p_agg = Path(plots_dir) / "mutational_spectrum_6class_aggregate.png"
    fig.savefig(p_agg, dpi=dpi)
    plt.close(fig)
    out.append(p_agg)

    # Heatmap of per-sample 6-class counts (if pandas available)
    try:
        import pandas as _pd
        dfm = _pd.DataFrame({s: d for s, d in muts}).T
        dfm = dfm.reindex(columns=KEYS).fillna(0).astype(int)
        fig, ax = plt.subplots(figsize=(max(6, 0.25 * dfm.shape[0]), 4))
        sns.heatmap(dfm, cmap="vlag", center=0, ax=ax, cbar_kws={"label":"count"})
        ax.set_xlabel("mutation")
        ax.set_ylabel("sample")
        ax.set_title("Per-sample 6-class mutation counts (heatmap)")
        plt.tight_layout()
        p_heat = Path(plots_dir) / "mutational_spectrum_6class_heatmap.png"
        fig.savefig(p_heat, dpi=dpi)
        plt.close(fig)
        out.append(p_heat)
    except Exception as e:
        print("Skipping mutation heatmap: " + str(e), file=sys.stderr)

    return out

def make_plots(records, plots_dir):
    """
    High-level plotting orchestration:
    - Ti/Tv bar plot
    - insert_median plot
    - low_mapq_frac plot
    - VAF violin per-sample
    - aggregate 6-class mutational spectrum
    """
    if plt is None or pd is None or sns is None:
        print("Plotting libraries missing; skipping plots.", file=sys.stderr)
        return
    plots_dir = Path(plots_dir)
    plots_dir.mkdir(parents=True, exist_ok=True)
    df = pd.DataFrame(records)
    for c in df.columns:
        # attempt to coerce numeric columns to numeric type for plotting convenience
        if c not in ("sample","bam","mut_six","mut_context96","vaf_list","counts_path"):
            try:
                df[c] = pd.to_numeric(df[c], errors="coerce")
            except Exception:
                pass
    sns.set(style="whitegrid")

    # Ti/Tv plot
    if "ti_tv_ratio" in df.columns:
        plot_df = df[["sample", "ti_tv_ratio"]].copy()
        plot_df["ti_tv_ratio"] = pd.to_numeric(plot_df["ti_tv_ratio"], errors="coerce")
        plot_df = plot_df.dropna(subset=["ti_tv_ratio"])
        if not plot_df.empty:
            plot_df = plot_df.sort_values("ti_tv_ratio", ascending=False)
            maxv = float(plot_df["ti_tv_ratio"].max())
            ylim_top = max(4.0, maxv * 1.1)
            plt.figure(figsize=(10,6))
            sns.barplot(x="sample", y="ti_tv_ratio", data=plot_df, palette="crest")
            plt.xticks(rotation=90)
            plt.ylabel("Ti/Tv")
            plt.ylim(0, ylim_top)
            plt.title("Transition/Transversion ratio per sample")
            plt.tight_layout()
            plt.savefig(plots_dir / "titv_per_sample.png")
            plt.close()

    # Insert size plot
    if "insert_median" in df.columns:
        plt.figure(figsize=(10,6))
        sns.barplot(x="sample", y="insert_median", data=df.sort_values("insert_median", ascending=False).fillna(0), palette="rocket")
        plt.xticks(rotation=90)
        plt.ylabel("Median insert size")
        plt.title("Median insert size per sample")
        plt.tight_layout()
        plt.savefig(plots_dir / "insert_median_per_sample.png")
        plt.close()

    # Low MAPQ fraction plot
    if "low_mapq_frac" in df.columns:
        plt.figure(figsize=(10,6))
        sns.barplot(x="sample", y="low_mapq_frac", data=df.sort_values("low_mapq_frac", ascending=False).fillna(0), palette="mako")
        plt.xticks(rotation=90)
        plt.ylabel("Fraction reads MAPQ < threshold")
        plt.title("Low MAPQ fraction per sample")
        plt.tight_layout()
        plt.savefig(plots_dir / "low_mapq_frac_per_sample.png")
        plt.close()

    # VAF violin per sample
    if any("vaf_list" in rec and rec["vaf_list"] for rec in records):
        rows = []
        for rec in records:
            vlist = rec.get("vaf_list")
            if vlist:
                for v in vlist:
                    try:
                        rows.append((rec["sample"], float(v)))
                    except Exception:
                        continue
        if rows:
            vdf = pd.DataFrame(rows, columns=["sample","vaf"])
            plt.figure(figsize=(max(8, 0.25 * vdf["sample"].nunique()), 6))
            sns.violinplot(x="sample", y="vaf", data=vdf, inner="box", cut=0)
            plt.xticks(rotation=90)
            plt.title("VAF distribution per sample")
            plt.tight_layout()
            plt.savefig(plots_dir / "vaf_violin_per_sample.png")
            plt.close()

    # Aggregate 6-class plot
    six_all = {}
    for rec in records:
        six = rec.get("mut_six")
        if isinstance(six, dict):
            for k,v in six.items():
                six_all[k] = six_all.get(k,0) + int(v)
    if six_all:
        keys = ["C>A","C>G","C>T","T>A","T>C","T>G"]
        vals = [six_all.get(k,0) for k in keys]
        plt.figure(figsize=(6,4))
        sns.barplot(x=keys, y=vals, palette="Set2")
        plt.ylabel("Count")
        plt.title("Aggregate 6-class mutational spectrum")
        plt.tight_layout()
        plt.savefig(plots_dir / "mutational_spectrum_6class.png")
        plt.close()

    try:
        sample_mut_plots = generate_sample_mutation_plots(records, plots_dir)
        if sample_mut_plots:
            print(f"Generated {len(sample_mut_plots)} mutation plot files under {plots_dir}")
    except Exception as e:
        print("Failed per-sample mutation plots:", e, file=sys.stderr)

    print(f"Plots written to {plots_dir}")

# -----------------------------------------------------------------------------
# Sample-sheet generation helpers (match VCF sample names to BAM filenames)
# -----------------------------------------------------------------------------
def strip_ext(name: str) -> str:
    """Trim common BAM/CRAM filename extensions."""
    for ext in (".sorted.bam", ".split.bam", ".bam", ".cram", ".bam.gz", ".cram.gz"):
        if name.endswith(ext):
            return name[:-len(ext)]
    return name

def trimmed_by_underscore(basename: str) -> str:
    """
    Heuristic to map BAM basenames to VCF sample IDs (trim at first underscore).
    - This is a dataset-specific heuristic; adjust if your naming differs.
    """
    clean = strip_ext(basename)
    if "_" in clean:
        return clean.split("_", 1)[0]
    return clean

def generate_sample_sheet_from_bams():
    """
    Enumerate BAM/CRAM files under BAMROOT, attempt to map them to VCF sample names
    using trimmed_by_underscore heuristic.
    - Writes a simple CSV sample sheet and also a renaming map (TSV) for inspection.
    - Returns (sheet_rows, basename_map, duplicates)
    """
    # First get samples from VCF (preferred via cyvcf2)
    vcf_samples = []
    try:
        v = VCF(str(VCF_PATH))
        vcf_samples = list(v.samples)
        v.close()
    except Exception:
        # fallback: bcftools query -l or parse header manually
        try:
            rc, out, err = run_cmd(["bcftools", "query", "-l", str(VCF_PATH)])
            if rc == 0:
                vcf_samples = [s.strip() for s in out.splitlines() if s.strip()]
        except Exception:
            try:
                opener = gzip.open if str(VCF_PATH).endswith(".gz") else open
                with opener(str(VCF_PATH), "rt") as fh:
                    for line in fh:
                        if line.startswith("#CHROM"):
                            parts = line.rstrip("\n").split("\t")
                            vcf_samples = parts[9:]
                            break
            except Exception as e:
                raise RuntimeError(f"Failed to read VCF samples: {e}")
    if not vcf_samples:
        raise RuntimeError("No samples found in VCF")
    vcf_set = set(vcf_samples)

    # Find BAM/CRAM files recursively under BAMROOT
    files = []
    for pat in ("*.bam", "*.cram"):
        files.extend(sorted(BAMROOT.rglob(pat)))
    seen = set()
    bams = []
    for p in files:
        try:
            rp = str(p.resolve())
        except Exception:
            rp = str(p)
        if rp not in seen:
            seen.add(rp)
            bams.append(Path(p))

    # Map trimmed basenames to vcf sample IDs if possible
    trim_to_bam = {}
    basename_map = {}
    duplicates = {}
    for b in bams:
        base = b.name
        tid = trimmed_by_underscore(base)
        if tid in vcf_set:
            if tid in trim_to_bam:
                # multiple BAMs map to same trimmed id (flag as duplicate)
                duplicates.setdefault(tid, []).append(str(b))
                basename_map[base] = ""
            else:
                trim_to_bam[tid] = str(b)
                basename_map[base] = tid
        else:
            basename_map[base] = ""  # unmatched
    sheet_rows = []
    for s in vcf_samples:
        if s in trim_to_bam:
            sheet_rows.append((s, trim_to_bam[s]))
    # Ensure parent directories exist
    SAMPLE_SHEET.parent.mkdir(parents=True, exist_ok=True)
    MAP_OUT.parent.mkdir(parents=True, exist_ok=True)
    # Write sample sheet CSV (sample_id,bam_path)
    with SAMPLE_SHEET.open("w", encoding="utf-8") as fh:
        fh.write("sample_id,bam_path\n")
        for sid, path in sheet_rows:
            fh.write(f"{sid},{path}\n")
    # Write renaming map for manual inspection; unmatched basenames have empty sample_id
    with MAP_OUT.open("w", encoding="utf-8") as fh:
        fh.write("bam_basename\tsample_id\n")
        for base, sid in basename_map.items():
            fh.write(f"{base}\t{sid}\n")
        if any(v == "" for v in basename_map.values()):
            fh.write("\n# Unmatched BAM basenames have empty sample_id entries for inspection\n")
    return sheet_rows, basename_map, duplicates

# -----------------------------------------------------------------------------
# Top-level per-sample worker that is picklable for ProcessPoolExecutor
# -----------------------------------------------------------------------------
def process_sample_task(sample_id, bam_path, vcf_path_str, prepared_ref_str, settings):
    """
    High-level worker executed per sample:
    - flagstat summary
    - BAM metrics (pysam or samtools)
    - optional featureCounts gene counts
    - mean MAPQ calculation
    - VAF computation (AD preferred)
    - genotype summary (GT)
    - mutational spectrum (6-class always; 96-context optional)
    - computes derived mutation metrics (Ti/Tv, CpG counts)
    Returns (record_dict, logs_list)
    """
    vcf_path = Path(vcf_path_str)
    prepared_ref = Path(prepared_ref_str) if prepared_ref_str else None
    rec = {"sample": sample_id, "bam": bam_path}
    logs = []
    # 1) flagstat
    try:
        fs = run_flagstat(bam_path)
        rec.update({
            "total_reads": fs.get("total"),
            "mapped_reads": fs.get("mapped"),
            "mapped_pct": fs.get("mapped_pct"),
            "primary_reads": fs.get("primary"),
            "secondary_reads": fs.get("secondary"),
            "duplicates": fs.get("duplicates"),
            "paired_in_seq": fs.get("paired"),
            "properly_paired": fs.get("properly_paired"),
            "properly_paired_pct": fs.get("properly_paired_pct"),
            "read1": fs.get("read1"),
            "read2": fs.get("read2"),
        })
        logs.append(f"flagstat total={fs.get('total')} mapped_pct={fs.get('mapped_pct')}")
    except Exception as e:
        rec["flagstat_error"] = str(e)
        logs.append(f"flagstat error: {e}")

    # 2) BAM metrics (insert sizes, softclip, MAPQ low fraction)
    if settings.get("compute_bam_metrics"):
        try:
            bm = compute_bam_metrics(bam_path, max_reads=settings.get("max_reads"), mapq_threshold=settings.get("mapq_threshold", 20))
            rec.update(bm)
            if "reads_examined" in bm:
                rec["reads_examined"] = int(bm.get("reads_examined") or 0)
            logs.append("bam metrics computed")
        except Exception as e:
            rec["bam_metrics_error"] = str(e)
            logs.append(f"bam metrics error: {e}")

    # 3) optional gene counts via featureCounts
    if settings.get("compute_gene_counts") and GTF:
        try:
            out_path, info = compute_gene_counts_featurecounts(bam_path, GTF, OUT_DIR, sample_id, threads=settings.get("threads", 1), stranded=settings.get("gene_counts_stranded", 0))
            if out_path:
                rec["gene_counts_path"] = str(out_path)
                if "genes_with_counts" in info:
                    rec["genes_with_counts"] = int(info["genes_with_counts"])
            else:
                rec["gene_counts_skipped_reason"] = info.get("skipped_reason", "unknown")
            logs.append("gene counts processed")
        except Exception as e:
            rec["gene_counts_error"] = str(e)
            logs.append(f"gene counts error: {e}")

    # 4) Mean MAPQ computation (small helper: compute_mean_mapq implemented below)
    if settings.get("compute_mapq"):
        try:
            mm = compute_mean_mapq(bam_path)
            rec["mean_mapq"] = mm
        except Exception as e:
            rec["mapq_error"] = str(e)

    # 5) Validate sample presence in VCF
    try:
        v = VCF(str(vcf_path))
        has_sample = sample_id in v.samples
        v.close()
    except Exception as e:
        rec["vcf_error"] = str(e)
        return rec, logs

    if not has_sample:
        rec["vcf_missing"] = True
        return rec, logs

    # 6) VAFs & diagnostics (AD preferred)
    try:
        vaf_list, counters = compute_sample_vafs(vcf_path, sample_id, min_depth=settings.get("min_depth", 10), max_variants=settings.get("max_variants"))
        stats = summarize_vafs(vaf_list)
        rec.update({
            "vaf_count": stats["count"],
            "vaf_mean": stats["mean"],
            "vaf_median": stats["median"],
            "vaf_stddev": stats["stddev"],
            "vaf_q25": stats["q25"],
            "vaf_q75": stats["q75"],
            "vaf_min": stats["min"],
            "vaf_max": stats["max"],
            "vaf_prop_ge_0.5": stats["prop_ge_0.5"],
            "vcf_total_variants_scanned": counters["total_seen"],
            "vcf_ad_missing": counters["ad_missing"],
            "vcf_ad_nonzero": counters.get("ad_nonzero", 0),
            "vcf_depth_filtered": counters["depth_filtered"],
            "vaf_list": vaf_list
        })
        logs.append(f"vaf variants considered: {stats['count']}")
    except Exception as e:
        rec["vaf_error"] = str(e)

    # 7) GT summary
    try:
        gt = compute_genotype_summary(vcf_path, sample_id, max_variants=settings.get("max_variants"))
        rec.update({
            "gt_total": gt["total"],
            "gt_hom_ref": gt["hom_ref"],
            "gt_het": gt["het"],
            "gt_hom_alt": gt["hom_alt"],
            "gt_missing": gt["missing"],
            "vcf_gt_nonref_count": gt.get("nonref_count", 0)
        })
        logs.append(f"genotype counts: total={gt['total']}")
    except Exception as e:
        rec["gt_error"] = str(e)

    # 8) Mutational spectrum (always produce 6-class; context optional)
    try:
        ms = compute_mutational_spectrum(vcf_path, sample_id, ref_fasta_path=prepared_ref if prepared_ref else None, max_variants=settings.get("max_variants"), min_depth=settings.get("min_depth", 10), min_vaf=settings.get("min_vaf", 0.0))
        if not ms.get("skipped"):
            rec["mut_six"] = ms.get("six", {})
            if ms.get("context96"):
                rec["mut_context96"] = ms.get("context96", {})
            if "variants_examined" in ms:
                rec["mutational_variants_examined"] = int(ms.get("variants_examined") or 0)
            logs.append("mutational spectrum computed")
            if sum(int(v or 0) for v in rec["mut_six"].values()) == 0:
                logs.append("warning: mut_six counts are all zero for this sample")
        else:
            rec["mut_spectrum_skipped_reason"] = ms.get("reason", "skipped")
    except Exception as e:
        rec["mut_spectrum_error"] = str(e)
        logs.append(f"mutational spectrum error: {e}")

    # 9) Derived mutation metrics: Ti/Tv, CpG enrichment etc.
    try:
        rec = compute_derived_mutation_metrics(rec)
    except Exception as e:
        rec["derived_mutation_error"] = str(e)

    return rec, logs

# -----------------------------------------------------------------------------
# compute_mean_mapq helper: fallback option using pysam or samtools view
# -----------------------------------------------------------------------------
def compute_mean_mapq(bam_path):
    """
    Compute mean mapping quality across alignments in BAM.
    - Use pysam if available; otherwise use samtools view (fast enough for many BAMs).
    - Excludes secondary/supplementary/unmapped reads.
    """
    if pysam is not None:
        try:
            total = 0; sum_mapq = 0.0
            with pysam.AlignmentFile(str(bam_path), "rb") as fh:
                for r in fh.fetch(until_eof=True):
                    if r.is_secondary or r.is_supplementary or r.is_unmapped:
                        continue
                    mq = r.mapping_quality if r.mapping_quality is not None else 0
                    sum_mapq += mq; total += 1
            return (sum_mapq / total) if total > 0 else None
        except Exception:
            pass
    if not shutil.which("samtools"):
        return None
    cmd = ["samtools", "view", "-F", "3328", str(bam_path)]
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    sum_mapq = 0.0; count = 0
    for line in p.stdout:
        parts = line.split("\t")
        if len(parts) > 4:
            try:
                mq = float(parts[4]); sum_mapq += mq; count += 1
            except Exception:
                continue
    p.stdout.close(); p.wait()
    return (sum_mapq / count) if count > 0 else None

# -----------------------------------------------------------------------------
# Derived mutation metrics (Ti/Tv, CpG-centric C->T fraction, etc.)
# -----------------------------------------------------------------------------
def compute_derived_mutation_metrics(rec):
    """
    Given a record dict that may contain 'mut_six' (6-class dict) and 'mut_context96'
    (dict mapping A[C>T]G -> count), compute Ti/Tv, fractions, and CpG counts.
    - Updates the record in-place and returns it.
    """
    six = rec.get("mut_six")
    if isinstance(six, str):
        try:
            import ast
            six = ast.literal_eval(six)
        except Exception:
            try:
                import json as _json
                six = _json.loads(six.replace("'", '"'))
            except Exception:
                six = {}
    has_six_counts = isinstance(six, dict) and any(int(six.get(k, 0) or 0) for k in ["C>A","C>G","C>T","T>A","T>C","T>G"])
    if not has_six_counts:
        # try derive from context96 if present
        ctx = rec.get("mut_context96")
        if isinstance(ctx, str):
            try:
                import ast
                ctx = ast.literal_eval(ctx)
            except Exception:
                try:
                    import json as _json
                    ctx = _json.loads(ctx.replace("'", '"'))
                except Exception:
                    ctx = {}
        if isinstance(ctx, dict) and ctx:
            derived = {"C>A":0,"C>G":0,"C>T":0,"T>A":0,"T>C":0,"T>G":0}
            for k, v in ctx.items():
                try:
                    cnt = int(v)
                except Exception:
                    try:
                        cnt = int(float(v))
                    except Exception:
                        cnt = 0
                # find the mutation string inside the context
                m = re.search(r"\[([ACGT]>[ACGT])\]", str(k))
                if not m:
                    m2 = re.search(r"([ACGT]>[ACGT])", str(k))
                    sub = m2.group(1) if m2 else None
                else:
                    sub = m.group(1)
                if sub in derived:
                    derived[sub] += cnt
            if any(v > 0 for v in derived.values()):
                six = derived
                rec["mut_six"] = six
                rec["mut_six_derived"] = True
                has_six_counts = True
    if has_six_counts:
        counts = {k: int(six.get(k, 0) or 0) for k in ["C>A","C>G","C>T","T>A","T>C","T>G"]}
        total = sum(counts.values())
        ti = counts.get("C>T", 0) + counts.get("T>C", 0)
        tv = total - ti
        rec["ti_count"] = ti
        rec["tv_count"] = tv
        rec["ti_tv_ratio"] = (ti / tv) if (tv and tv > 0) else None
        rec["ti_frac"] = (ti / total) if (total and total > 0) else None
        rec["pct_C_to_T"] = (counts.get("C>T", 0) / total) if (total and total > 0) else None
    else:
        rec.update({k: None for k in ["ti_count","tv_count","ti_tv_ratio","ti_frac","pct_C_to_T"]})
    # If context96 present, compute CpG-specific C->T counts
    ctx = rec.get("mut_context96")
    if isinstance(ctx, dict) and ctx:
        CpG_ct = 0
        CtoT_total = 0
        for k, v in ctx.items():
            try:
                cnt = int(v)
            except Exception:
                try:
                    cnt = int(float(v))
                except Exception:
                    cnt = 0
            if "[C>" in k:
                CtoT_total += cnt
                # simplistic check: last base 'G' indicates CpG on right context (A[C>T]G)
                if k.strip()[-1].upper() == "G":
                    CpG_ct += cnt
        rec["CpG_CtoT_count"] = CpG_ct
        rec["CpG_fraction_CtoT"] = (CpG_ct / CtoT_total) if (CtoT_total and CtoT_total > 0) else None
    else:
        rec["CpG_CtoT_count"] = None
        rec["CpG_fraction_CtoT"] = None
    return rec

# -----------------------------------------------------------------------------
# Utility group statistics functions (median, MAD etc.)
# -----------------------------------------------------------------------------
def _to_float(v):
    if v is None:
        return None
    try:
        return float(v)
    except Exception:
        try:
            return float_from(v)
        except Exception:
            return None

def _quantile(sorted_vals, p):
    n = len(sorted_vals)
    if n == 0:
        return None
    idx = int(p * (n - 1))
    idx = max(0, min(n - 1, idx))
    return sorted_vals[idx]

def _mad(vals):
    if not vals:
        return None
    med = median(vals)
    devs = [abs(x - med) for x in vals]
    mad = median(devs)
    return mad * 1.4826  # consistent scaling to estimate std

def compute_group_stats(records, group_key_fn, keys_to_aggregate):
    """
    Compute group-wise statistics for lists of numeric keys across records.
    Returns dict: group -> metric -> stats-dict
    stats-dict includes mean, std, median, q25, q75, iqr, mad, sem, n
    """
    groups = {}
    for rec in records:
        group = group_key_fn(rec)
        if group is None:
            group = "UNKNOWN"
        groups.setdefault(group, []).append(rec)
    out = {}
    for g, recs in groups.items():
        stats = {}
        for key in keys_to_aggregate:
            vals = []
            for r in recs:
                v = _to_float(r.get(key))
                if v is None:
                    continue
                if isinstance(v, float) and (math.isnan(v) or math.isinf(v)):
                    continue
                vals.append(v)
            if not vals:
                stats[key] = {"mean": None, "std": None, "median": None, "q25": None, "q75": None,
                              "iqr": None, "mad": None, "sem": None, "n": 0}
            else:
                vals_sorted = sorted(vals)
                n = len(vals_sorted)
                try:
                    m = mean(vals_sorted)
                except Exception:
                    m = None
                try:
                    s = pstdev(vals_sorted) if n > 1 else 0.0
                except Exception:
                    s = None
                med = median(vals_sorted)
                q25 = _quantile(vals_sorted, 0.25)
                q75 = _quantile(vals_sorted, 0.75)
                iqr = (q75 - q25) if (q25 is not None and q75 is not None) else None
                mad = _mad(vals_sorted)
                sem = (s / math.sqrt(n)) if (s is not None and n > 0) else None
                stats[key] = {"mean": m, "std": s, "median": med, "q25": q25, "q75": q75,
                              "iqr": iqr, "mad": mad, "sem": sem, "n": n}
        out[g] = stats
    return out

# -----------------------------------------------------------------------------
# MAIN: orchestrate sample discovery, processing, summary, and plotting
# -----------------------------------------------------------------------------
def main():
    # Ensure output directories exist
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    PLOTS_DIR.mkdir(parents=True, exist_ok=True)

    max_reads_cfg = SETTINGS.get("max_reads")
    max_vars_cfg = SETTINGS.get("max_variants")
    print(f"CONFIG: max_reads={max_reads_cfg!s}, max_variants={max_vars_cfg!s}")
    if max_reads_cfg is None or max_vars_cfg is None:
        print("WARNING: one or both of max_reads/max_variants are None — the script will perform full scans which may be slow.", file=sys.stderr)

    print("Preparing sample sheet and inputs (hard-coded)...")
    sheet_rows, basename_map, duplicates = generate_sample_sheet_from_bams()
    print(f"Wrote sample sheet to {SAMPLE_SHEET} ({len(sheet_rows)} matched samples)")
    if duplicates:
        print("WARNING: duplicates found for some trimmed IDs; see map file.")

    samples = sheet_rows
    if not samples:
        print("No samples matched. Exiting.", file=sys.stderr)
        return

    # prepare reference for context-aware spectrum if requested
    def prepare_reference_for_indexing_inline(ref_path):
        """
        Ensure the reference fasta is indexable by samtools/pysam.
        - If .fai exists, return the ref path
        - If compressed in a way that samtools cannot index, attempt to create bgzipped or uncompressed copy and index it
        - Returns path to a usable reference (may be a temp file)
        """
        ref = Path(ref_path)
        if not ref.exists():
            raise RuntimeError(f"Reference not found: {ref}")
        fai_guess = ref.parent / (ref.name + ".fai")
        if fai_guess.exists():
            return ref
        rc, out, err = run_cmd(["samtools", "faidx", str(ref)])
        if rc == 0:
            return ref
        err_l = (err or "").lower()
        # Some samtools versions cannot index gzipped files missing the right compression (e.g., plain gzip)
        if "cannot index files compressed with gzip" not in err_l and rc != 0:
            raise RuntimeError(f"samtools faidx failed for {ref}:\n{err}")
        # Attempt bgzip compression if bgzip available
        base_no_gz = ref.name[:-3] if ref.name.endswith(".gz") else ref.name
        target_bgz = ref.parent / (base_no_gz + ".bgz.fa.gz")
        target_uncompressed = ref.parent / (base_no_gz + ".uncompressed.fa")
        if shutil.which("bgzip"):
            try:
                with open(target_bgz, "wb") as outfh:
                    p = subprocess.Popen(["bgzip", "-c", str(ref)], stdout=outfh, stderr=subprocess.PIPE, text=True)
                    _, berr = p.communicate()
                    if p.returncode == 0:
                        rc2, out2, err2 = run_cmd(["samtools", "faidx", str(target_bgz)])
                        if rc2 == 0:
                            return target_bgz
            except Exception:
                pass
        # As a last resort, uncompress to a temporary plain fasta and index it
        with gzip.open(str(ref), "rb") as inf, open(str(target_uncompressed), "wb") as outf:
            shutil.copyfileobj(inf, outf)
        rc3, out3, err3 = run_cmd(["samtools", "faidx", str(target_uncompressed)])
        if rc3 == 0:
            return target_uncompressed
        raise RuntimeError("Failed to prepare reference for indexing")

    if REF_FA and REF_FA.exists():
        try:
            prepared_ref = prepare_reference_for_indexing_inline(REF_FA)
            print("Prepared reference for indexing:", prepared_ref)
            prepared_ref_str = str(prepared_ref)
        except Exception as e:
            print("Warning: reference preparation failed:", e, file=sys.stderr)
            prepared_ref_str = ""
    else:
        prepared_ref_str = ""

    records = []
    text_lines = ["Results summary", f"BAMROOT: {BAMROOT}", f"Sample sheet: {SAMPLE_SHEET}", f"VCF: {VCF_PATH}",
                  f"min_depth: {SETTINGS['min_depth']}", f"min_vaf: {SETTINGS['min_vaf']}", f"max_variants: {SETTINGS['max_variants']}", ""]

    tasks = [(sid, bam) for (sid, bam) in samples]

    # Parallel processing using ProcessPoolExecutor; keep the worker function picklable and avoid closures
    if SETTINGS["jobs"] > 1:
        with ProcessPoolExecutor(max_workers=SETTINGS["jobs"]) as exe:
            futures = {exe.submit(process_sample_task, sid, bam, str(VCF_PATH), prepared_ref_str, SETTINGS): sid for sid, bam in tasks}
            for fut in as_completed(futures):
                sid = futures[fut]
                try:
                    rec, logs = fut.result()
                    records.append(rec)
                    text_lines.append("="*60)
                    text_lines.append(f"Sample: {rec.get('sample')}")
                    text_lines.append(f"BAM: {rec.get('bam')}")
                    for L in logs:
                        text_lines.append(L)
                    print(f"[done] {sid}")
                except Exception as e:
                    print(f"[failed] {sid}: {e}", file=sys.stderr)
                    records.append({"sample": sid, "worker_error": str(e)})
    else:
        for sid, bam in tasks:
            rec, logs = process_sample_task(sid, bam, str(VCF_PATH), prepared_ref_str, SETTINGS)
            records.append(rec)
            text_lines.append("="*60)
            text_lines.append(f"Sample: {rec.get('sample')}")
            text_lines.append(f"BAM: {rec.get('bam')}")
            for L in logs:
                text_lines.append(L)

    # ---------------- compute group (family) statistics ----------------
    def sample_family_fn(rec):
        s = rec.get("sample")
        if not s:
            return "UNKNOWN"
        return trimmed_by_underscore(s)

    group_keys = [
        "total_reads","mapped_reads","mapped_pct","primary_reads","secondary_reads","duplicates",
        "paired_in_seq","properly_paired","properly_paired_pct","read1","read2",
        "insert_count","insert_mean","insert_median","insert_std","low_mapq_frac","pct_softclipped",
        "mean_mapq","mean_depth","vaf_count","vaf_mean","vaf_median","vaf_stddev","vaf_q25","vaf_q75",
        "vaf_min","vaf_max","vaf_prop_ge_0.5","vcf_total_variants_scanned","vcf_ad_missing","vcf_ad_nonzero","vcf_depth_filtered",
        "vcf_gt_nonref_count","gt_total","gt_hom_ref","gt_het","gt_hom_alt","gt_missing",
        "ti_count","tv_count","ti_tv_ratio","ti_frac","pct_C_to_T","CpG_CtoT_count","CpG_fraction_CtoT",
        "reads_examined","mutational_variants_examined","genes_with_counts"
    ]

    family_stats = compute_group_stats(records, sample_family_fn, group_keys)
    overall_stats_dict = compute_group_stats(records, lambda r: "ALL", group_keys)
    overall_stats = overall_stats_dict.get("ALL", {})

    # Append human-readable family & overall block
    text_lines.append("\n" + "="*60)
    text_lines.append("Group / Family summary (classic + robust stats)")
    text_lines.append("family\tn\tmetric\tmean±std\tmedian\tq25-q75\tIQR\tMAD\tSEM\tn_values")
    for fam in sorted(family_stats.keys()):
        fam_stats = family_stats[fam]
        for k in group_keys:
            s = fam_stats.get(k, None)
            if s is None:
                continue
            mean_v = s.get("mean"); std_v = s.get("std"); med = s.get("median")
            q25 = s.get("q25"); q75 = s.get("q75")
            iqr = s.get("iqr"); mad = s.get("mad"); sem = s.get("sem"); nn = s.get("n", 0)
            if mean_v is None and med is None:
                continue
            meanstd = f"{mean_v:.3g} ±{std_v:.3g}" if (mean_v is not None and std_v is not None) else (f"{mean_v:.3g}" if mean_v is not None else "")
            median_str = f"{med:.3g}" if med is not None else ""
            qrange = f"{q25:.3g}-{q75:.3g}" if (q25 is not None and q75 is not None) else ""
            iqr_s = f"{iqr:.3g}" if iqr is not None else ""
            mad_s = f"{mad:.3g}" if mad is not None else ""
            sem_s = f"{sem:.3g}" if sem is not None else ""
            text_lines.append(f"{fam}\t{nn}\t{k}\t{meanstd}\t{median_str}\t{qrange}\t{iqr_s}\t{mad_s}\t{sem_s}\t{nn}")

    text_lines.append("\nOverall (ALL groups):")
    for k in group_keys:
        s = overall_stats.get(k)
        if not s:
            continue
        mean_v = s.get("mean"); std_v = s.get("std"); med = s.get("median")
        q25 = s.get("q25"); q75 = s.get("q75")
        iqr = s.get("iqr"); mad = s.get("mad"); sem = s.get("sem"); nn = s.get("n", 0)
        if mean_v is None and med is None:
            continue
        meanstd = f"{mean_v:.3g} ±{std_v:.3g}" if (mean_v is not None and std_v is not None) else (f"{mean_v:.3g}" if mean_v is not None else "")
        median_str = f"{med:.3g}" if med is not None else ""
        qrange = f"{q25:.3g}-{q75:.3g}" if (q25 is not None and q75 is not None) else ""
        iqr_s = f"{iqr:.3g}" if iqr is not None else ""
        mad_s = f"{mad:.3g}" if mad is not None else ""
        sem_s = f"{sem:.3g}" if sem is not None else ""
        text_lines.append(f"ALL\t{nn}\t{k}\t{meanstd}\t{median_str}\t{qrange}\t{iqr_s}\t{mad_s}\t{sem_s}\t{nn}")
    text_lines.append("="*60 + "\n")

    # Machine-readable family_summary.tsv
    fam_tsv = OUT_DIR / "family_summary.tsv"
    with fam_tsv.open("w") as fh:
        cols = ["family","metric","n","mean","std","median","q25","q75","iqr","mad","sem"]
        fh.write("\t".join(cols) + "\n")
        for fam in sorted(family_stats.keys()):
            fam_stats = family_stats[fam]
            for k in group_keys:
                s = fam_stats.get(k)
                if s is None:
                    fh.write("\t".join([fam, k, "0"] + [""] * (len(cols)-3)) + "\n")
                else:
                    fh.write("\t".join([
                        fam,
                        k,
                        str(s.get("n",0)),
                        "" if s.get("mean") is None else f"{s['mean']:.6g}",
                        "" if s.get("std") is None else f"{s['std']:.6g}",
                        "" if s.get("median") is None else f"{s['median']:.6g}",
                        "" if s.get("q25") is None else f"{s['q25']:.6g}",
                        "" if s.get("q75") is None else f"{s['q75']:.6g}",
                        "" if s.get("iqr") is None else f"{s['iqr']:.6g}",
                        "" if s.get("mad") is None else f"{s['mad']:.6g}",
                        "" if s.get("sem") is None else f"{s['sem']:.6g}",
                    ]) + "\n")
    print(f"Wrote family summary TSV -> {fam_tsv}")

    # -----------------------------------------------------------------------------
    # Output: text, JSON, TSV. Complex fields are JSON-encoded in the TSV to make
    # downstream parsing robust (e.g., mut_six, mut_context96, vaf_list).
    # -----------------------------------------------------------------------------
    out_txt = OUT_DIR / "results_summary.txt"
    out_json = OUT_DIR / "results_summary.json"
    out_tsv = OUT_DIR / "results_summary.tsv"
    out_txt.write_text("\n".join(text_lines))
    try:
        out_json.write_text(json.dumps(records, indent=2))
    except Exception as e:
        print("Failed to write JSON:", e, file=sys.stderr)
    try:
        if records:
            # Define a preferred order for columns, then append any extras found in records
            keys_priority = ["sample","bam","total_reads","mapped_reads","mapped_pct","primary_reads","secondary_reads","duplicates",
                             "paired_in_seq","properly_paired","properly_paired_pct","read1","read2",
                             "insert_count","insert_mean","insert_median","insert_std","low_mapq_frac","pct_softclipped",
                             "mean_mapq","mean_depth","vaf_count","vaf_mean","vaf_median","vaf_stddev","vaf_q25","vaf_q75",
                             "vaf_min","vaf_max","vaf_prop_ge_0.5","vcf_total_variants_scanned","vcf_ad_missing","vcf_ad_nonzero","vcf_depth_filtered",
                             "vcf_gt_nonref_count","gt_total","gt_hom_ref","gt_het","gt_hom_alt","gt_missing","ti_count","tv_count","ti_tv_ratio",
                             "ti_frac","pct_C_to_T","CpG_CtoT_count","CpG_fraction_CtoT","reads_examined","mutational_variants_examined","genes_with_counts","mut_six","mut_context96","vaf_list"]
            all_keys = []
            for r in records:
                for k in r.keys():
                    if k not in all_keys:
                        all_keys.append(k)
            keys = [k for k in keys_priority if k in all_keys]
            keys += [k for k in all_keys if k not in keys]
            with open(out_tsv, "w") as fh:
                fh.write("\t".join(keys) + "\n")
                for rec in records:
                    row = []
                    for k in keys:
                        val = rec.get(k, "")
                        # JSON-encode lists/dicts to keep TSV parseable
                        if isinstance(val, (dict, list)):
                            try:
                                row.append(json.dumps(val, separators=(",", ":")))
                            except Exception:
                                row.append(str(val))
                        else:
                            row.append("" if val is None else str(val))
                    fh.write("\t".join(row) + "\n")
    except Exception:
        pass

    print(f"Wrote text -> {out_txt}, json -> {out_json}, tsv -> {out_tsv}")
    # Produce plots if requested
    if SETTINGS["make_plots"]:
        try:
            make_plots(records, PLOTS_DIR)
        except Exception as e:
            print("Plotting failed:", e, file=sys.stderr)

# -----------------------------------------------------------------------------
# If script is invoked directly, run main()
# -----------------------------------------------------------------------------
if __name__ == "__main__":
    main()

# -----------------------------------------------------------------------------
# POST-SCRIPT: Practical notes, troubleshooting, and suggested improvements
# -----------------------------------------------------------------------------
#
# 1) Speed & scaling:
#    - Full VCF and BAM scans are expensive. Use SETTINGS["max_variants"] and
#      SETTINGS["max_reads"] during development to run quick checks.
#    - For production, consider region-limited scanning (e.g., exome capture regions)
#      or using bcftools query to pre-filter VCF to PASS variants and SNPs only.
#    - For BAM scanning, tools like mosdepth or sambamba can compute coverage faster.
#
# 2) Reproducibility:
#    - Record exact package versions: cyvcf2, pysam, samtools, featureCounts, matplotlib, seaborn.
#      The script writes session info nowhere by default — consider adding that.
#    - Containerize (Docker or Singularity) the environment to ensure reproducible runs.
#
# 3) Robustness:
#    - There are many VCF variants of FORMAT usage. AD may be missing or may represent per-ALT counts.
#      This script sums per-ALT counts for VAF; verify this matches your assumed convention.
#    - GT encodings may be phased (|) — cyvcf2 returns numeric allele indices regardless of phase.
#    - Multi-allelic sites are handled conservatively: alt_count sums all non-ref AD elements when computing VAF.
#
# 4) Improvements you may want:
#    - Convert to a CLI-driven tool (argparse / click) rather than hard-coded paths.
#    - Replace print-based logging with Python logging module with loglevels and file handlers.
#    - Add unit tests:
#        * Small synthetic VCF (2-3 variants) and BAM fixtures (samtools viewable) to validate compute_sample_vafs, compute_mutational_spectrum, and compute_bam_metrics.
#    - Add a checkpointing mechanism to persist per-sample results as they are produced so long runs can be resumed.
#
# 5) Edge cases:
#    - If VCF uses unusual genotype encodings (e.g., phased haploid GATK-compat), adjust GT parsing logic accordingly.
#    - Reference fasta sequence names must match VCF/VCF chromosome names for context lookup (e.g., 'chr1' vs '1').
#    - samtools / pysam failure modes: "Can't index gzip" etc. The script tries bgzip/uncompress heuristics but test carefully.
#
# 6) Downstream:
#    - The TSV JSON-encodes complex columns (mut_six, mut_context96, vaf_list) to make it straightforward
#      to ingest in R or Python for further analysis (use json.loads).
#    - Consider exporting smaller summary CSVs for each metric separately if consumers need column-normalized tables.
#
# If you'd like, I can:
#  - Convert this script to accept arguments and produce a runnable CLI.
#  - Add robust logging and sessionInfo saving.
#  - Implement example unit tests and small VCF/BAM fixtures for CI.
#  - Optimize VCF scanning using bcftools + bcftools query where appropriate.
#
# End of annotated script.
