#!/bin/bash
#
# isec_annotation_annotated.sh
#
# Heavily annotated wrapper script that:
#  - Checks for a usable Python interpreter (a conda env in this user's setup)
#  - Verifies presence of a VEP-annotated VCF and its index
#  - Invokes an embedded Python script (here-document) that:
#      * Reads a simple positions file (chrom pos per line)
#      * Annotates each position by consulting:
#          - a VEP-annotated VCF (CSQ in INFO) for variant-level consequences
#          - a bgzipped/tabix-indexed GTF (fallback gene overlap if VEP absent)
#      * Writes a TSV summarising presence in VCF, REF/ALT/QUAL/INFO, nearby genes,
#        and parsed VEP fields (Consequence, Gene, SYMBOL)
#
# This annotated version documents intent, failure modes, and suggestions for
# improvements and robustness at each step. The script is left runnable as-is
# (do not forget to edit the hard-coded paths at the top to match your system).
#
# RUN:
#   Edit the HARD-CODED PATHS section below and run:
#     ./isec_annotation_annotated.sh
#
# NOTES:
#  - This wrapper uses a specific Python interpreter from the user's conda env.
#    You can change PYTHON to a different Python executable or simply 'python3'.
#  - The Python body uses pysam (pysam.VariantFile and pysam.TabixFile). Ensure
#    pysam is installed in the selected interpreter.
#  - The script expects BGZIP+TABIX indexed VCF and GTF files (.tbi present).
#    VEP annotations are expected in the CSQ INFO field; the script parses the
#    CSQ header to extract field order (necessary to map values to names).
#
# SECURITY / REPRODUCIBILITY:
#  - Avoid running arbitrary VCFs from untrusted sources through this script.
#  - For reproducibility, pin the Python environment (conda env yaml) or use
#    a container image that includes the expected pysam version.
#

# -------------------- USER-EDITABLE: interpreter & paths --------------------
# Path to the Python interpreter to use for the embedded Python script.
# Points at a user's conda env in the original context; replace with your env or `python3`.
PYTHON="/home/jgdaniel/miniconda3/envs/vep_env/bin/python"

# Quick check: make sure the interpreter is present and executable.
if [ ! -x "${PYTHON}" ]; then
  echo "ERROR: interpreter not found or not executable: ${PYTHON}" >&2
  exit 1
fi

# === HARD-CODED PATHS: change ONLY these if you intend to run the script ===
# POS_FILE: a text file of positions (chrom pos per line) that you want to annotate.
# Example format (two columns separated by whitespace or tab):
#   chr1    1234567
#   chr2    7654321
#
POS_FILE="/home/jgd/Documents/bioinformatics_working/2025-11-06_enu-proc/intersection/2026-03-02_isec_with_ck-data/per_family/au01/lists/nothrom_specific_noBG.tsv"

# VEP_VCF: VEP-annotated VCF (bgzipped) that contains CSQ annotations in INFO.
# This script uses the VEP-annotated VCF as the primary source for consequence info.
VEP_VCF="/home/jgd/Documents/bioinformatics_working/2025-11-06_enu-proc/vcf/gatk/gatk-ensembl/vep_annotated/gatk-ensembl_filtered.20260224.nodash.vep.vcf.gz"

# GTF_FILE: bgzipped GTF with tabix index (.tbi) used as a fallback to fetch overlapping genes.
GTF_FILE="/home/jgd/Documents/bioinformatics_working/ref/danio_rerio.grcz11.115.gtf.gz"

# OUT_FILE: where the annotated TSV will be written. Parent directory will be created if needed.
OUT_FILE="/home/jgd/Documents/bioinformatics_working/output/isec_annotation/au01_annotated.tsv"

# VEP_SPECIES / VEP_ASSEMBLY: informational only here (we don't invoke VEP inside this wrapper).
VEP_SPECIES="danio_rerio"
VEP_ASSEMBLY="GRCz11"

# VEP_FORKS: how many parallel jobs you might give to VEP if you were to call it.
# Not used in this wrapper but kept as a useful parameter reminder.
VEP_FORKS="${SLURM_CPUS_PER_TASK:-$(nproc)}"
# -------------------------------------------------------------------------------

echo "Using python: ${PYTHON}"

# Basic interpreter check: run a short Python snippet that imports pysam and prints versions
# This is a helpful early failure for missing pysam or misconfigured interpreters.
"${PYTHON}" - <<'PYTEST'
import sys, traceback
try:
    import pysam
    print("Python executable:", sys.executable)
    print("Python version:", sys.version.replace('\n',' '))
    print("pysam version:", getattr(pysam, "__version__", "unknown"))
except Exception:
    traceback.print_exc()
    raise SystemExit(2)
PYTEST

# Verify the VEP VCF and its index exist. We require both the bgzipped VCF and .tbi index.
if [ -f "${VEP_VCF}" ] && [ -f "${VEP_VCF}.tbi" ]; then
  echo "Found existing VEP-annotated VCF and index: ${VEP_VCF} ${VEP_VCF}.tbi"
else
  echo "ERROR: VEP-annotated VCF or index not found. Please ensure they exist."
  exit 1
fi

# Now invoke the Python annotation logic. We pass through any CLI args ($@) to the Python script.
# The Python code is embedded below; it is a self-contained script that will execute under the selected interpreter.
"${PYTHON}" - "$@" <<'PYCODE'
#!/usr/bin/env python3
"""
Embedded Python: annotate positions using a VEP-annotated VCF and a bgzipped GTF.

This code block is executed by the shell wrapper above. It is intentionally
self-contained and uses pysam for efficient indexed access to both the VCF and GTF.

Key behavior:
 - For each (chrom,pos) pair from a positions file it tries to fetch any VCF records
   at that position. If found, it collects REF/ALT/QUAL/INFO and attempts to extract
   the VEP CSQ entry (parsing the VCF header to discover the CSQ field order).
 - If no VCF record is found, or to supplement VEP output, it queries the GTF
   (tabix indexed) for overlapping gene features and collects gene names.
 - The output is a TSV with these columns:
     CHROM  POS  PRESENT(yes/no)  REF  ALT  QUAL  INFO  GENES  VEP_CONSEQUENCE  VEP_GENE  VEP_SYMBOL
 - The script takes care to:
     * validate input files and indexes,
     * create the output directory if needed and check write permissions,
     * be conservative when parsing complex INFO/CSQ fields,
     * tolerate missing fields and provide informative fallbacks.
"""

# Standard libraries
import os
import sys
import re

# Third-party: pysam (must be installed in the PYTHON environment)
try:
    import pysam
except Exception as e:
    sys.stderr.write("pysam is required: {}\n".format(e))
    raise

# ---- NOTE ----
# These constants are intentionally duplicated from the wrapper to keep the embedded script
# self-contained and easy to extract if you want to run it directly with Python.
# Edit them here if you prefer to run the python part directly instead of the wrapper.
POS_FILE = "/home/jgd/Documents/bioinformatics_working/2025-11-06_enu-proc/intersection/2026-03-02_isec_with_ck-data/per_family/ai08/lists/nothrom_specific_noBG.tsv"
VEP_VCF = "/home/jgd/Documents/bioinformatics_working/2025-11-06_enu-proc/vcf/gatk/gatk-ensembl/vep_annotated/gatk-ensembl_filtered.20260224.nodash.vep.vcf.gz"
GTF_FILE = "/home/jgd/Documents/bioinformatics_working/ref/danio_rerio.grcz11.115.gtf.gz"
OUT_FILE = "/home/jgd/Documents/bioinformatics_working/output/isec_annotation/ai08_annotated.tsv"
# ---------------------------------------------------------------------------------

def fail(msg, rc=1):
    """Print an error message to stderr and exit with rc."""
    sys.stderr.write(msg + "\n")
    sys.exit(rc)

# ---------------------------------------------------------------------------------
# Parsing helpers
# ---------------------------------------------------------------------------------
def parse_attrs(attr_string):
    """
    Quick GTF attribute parser.
    GTF attributes often look like: key "value"; key2 "value2";
    This parser extracts those key/value pairs into a dict.
    Notes:
      - This is a minimal parser sufficient for standard GTFs.
      - It does not fully implement all edge-cases of GTF quoting or escaped characters.
    """
    d = {}
    # split on ';' then split first whitespace to separate key and value
    for part in attr_string.strip().rstrip(';').split(';'):
        part = part.strip()
        if not part:
            continue
        # Many GTFs use 'key "value"' style; we split on first space
        if ' ' in part:
            k, v = part.split(' ', 1)
            # strip surrounding quotes if present
            d[k] = v.strip().strip('"')
    return d

def parse_positions(path):
    """
    Generator that yields (chrom, pos) for each non-empty, well-formed line
    in the POS_FILE. Accepts files with whitespace-separated columns where
    the first two columns are CHROM and POS.
    Skips lines that don't contain a valid integer position and prints a warning.
    """
    with open(path, 'r') as fh:
        for lineno, line in enumerate(fh, 1):
            s = line.strip()
            if not s:
                continue
            cols = s.split()
            if len(cols) < 2:
                sys.stderr.write(f"Skipping malformed line {lineno}: {s}\n")
                continue
            chrom = cols[0]
            try:
                pos = int(cols[1])
            except ValueError:
                sys.stderr.write(f"Skipping line with non-integer pos {lineno}: {s}\n")
                continue
            yield chrom, pos

def genes_at_pos(gtf_tb, chrom, pos):
    """
    Query the tabix-indexed GTF for features overlapping 'pos' on 'chrom'.
    Extract gene_name, gene_id, or transcript_id (in that preference order).
    Returns a sorted list of unique gene identifiers overlapping the coordinate.
    If the chromosome is not present in the GTF index, fetch() may raise ValueError;
    we handle it and return an empty list.
    """
    genes = set()
    try:
        # pysam.TabixFile.fetch uses 0-based start, end; rec lines are raw GTF lines
        for rec in gtf_tb.fetch(chrom, pos-1, pos):
            if not rec or rec.startswith('#'):
                continue
            cols = rec.split('\t')
            # Validate columns length
            if len(cols) < 9:
                continue
            ad = parse_attrs(cols[8])
            # Prefer human-friendly gene_name if present, else gene_id or transcript_id
            if 'gene_name' in ad:
                genes.add(ad['gene_name'])
            elif 'gene_id' in ad:
                genes.add(ad['gene_id'])
            elif 'transcript_id' in ad:
                genes.add(ad['transcript_id'])
    except ValueError:
        # chromosome not present in index: ignore and return empty list (caller treats as INTERGENIC)
        pass
    return sorted(genes)

def ensure_outdir(out_path):
    """
    Ensure parent directory exists and is writable.
    Creates the parent directory with mode 2775 (group writable) by default.
    Performs a write test to ensure we have permission to create files there.
    """
    parent = os.path.dirname(os.path.abspath(out_path))
    if not parent:
        return
    try:
        os.makedirs(parent, mode=0o2775, exist_ok=True)
    except Exception as e:
        fail(f"Failed to create output directory '{parent}': {e}")
    try:
        testfile = os.path.join(parent, f".write_test_{os.getpid()}")
        with open(testfile, "w") as fh:
            fh.write("ok")
        os.remove(testfile)
    except Exception as e:
        fail(f"No write permission in output directory '{parent}': {e}")

# ---------------------------------------------------------------------------------
# VEP CSQ parsing helpers (header-driven)
# ---------------------------------------------------------------------------------
def get_vep_csq_fields_from_header(vcf_variantfile):
    """
    VEP stores consequence annotations in INFO field CSQ using a pipe-delimited format.
    The VCF header contains a line like:
      ##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. Format: Allele|Consequence|IMPACT|SYMBOL|...">
    This function extracts the "Format: ..." list of field names (the order matters).
    Returns a list of field names (or None if not found).
    """
    header_text = str(vcf_variantfile.header)
    for line in header_text.splitlines():
        if line.startswith("##INFO=<ID=CSQ"):
            # Find the `Format: ...` segment inside the header
            m = re.search(r'Format:\s*([^">]+)', line)
            if m:
                fields = [f.strip() for f in m.group(1).split('|')]
                return fields
    return None

def parse_csq_entry(entry_str, csq_fields):
    """
    Given a single CSQ string (pipe-delimited) and the csq_fields list (from header),
    return a dict mapping field name -> value. If there are fewer parts than fields,
    pad with empty strings. If csq_fields is None, return empty dict.
    """
    if not entry_str or not csq_fields:
        return {}
    parts = entry_str.split('|')
    # pad if necessary
    if len(parts) < len(csq_fields):
        parts += [''] * (len(csq_fields) - len(parts))
    return dict(zip(csq_fields, parts))

# ---------------------------------------------------------------------------------
# Main annotation routine
# ---------------------------------------------------------------------------------
def annotate_positions_file(pos_file, vep_vcf, gtf_file, out_file):
    """
    Orchestrates annotation:
    - Ensures outputs dirs and input files exist
    - Opens VEP VCF (pysam.VariantFile) and GTF (pysam.TabixFile)
    - For each position: try to find matching VCF records; parse REF/ALT/QUAL/INFO; parse CSQ
      If no VCF records found, still attempt GTF overlap
    - Write a tab-delimited row per position with selected annotation fields
    """
    ensure_outdir(out_file)

    # Input existence and basic index checks — fail early with clear messages
    if not os.path.isfile(pos_file):
        fail(f"Positions file not found: {pos_file}")
    if not os.path.isfile(vep_vcf):
        fail(f"VEP-annotated VCF file not found: {vep_vcf}")
    if not os.path.isfile(vep_vcf + ".tbi"):
        fail(f"VEP-annotated VCF index not found: {vep_vcf}.tbi (tabix -p vcf)")
    if not os.path.isfile(gtf_file):
        fail(f"GTF file not found: {gtf_file}")
    if not os.path.isfile(gtf_file + ".tbi"):
        fail(f"GTF index not found: {gtf_file}.tbi (tabix -p gff)")

    # Open VCF with pysam.VariantFile for indexed access and header inspection
    try:
        vcf_in = pysam.VariantFile(vep_vcf)
    except Exception as e:
        fail(f"Failed to open VEP VCF/index '{vep_vcf}': {e}")

    # Open GTF as TabixFile for coordinate queries
    try:
        gtf_tb = pysam.TabixFile(gtf_file)
    except Exception as e:
        fail(f"Failed to open GTF/index '{gtf_file}': {e}")

    # Discover CSQ header fields (order of fields in VEP's CSQ pipe-delimited annotation)
    csq_fields = get_vep_csq_fields_from_header(vcf_in)

    # Open output TSV and write header
    try:
        with open(out_file, "w") as out:
            out.write("CHROM\tPOS\tPRESENT\tREF\tALT\tQUAL\tINFO\tGENES\tVEP_CONSEQUENCE\tVEP_GENE\tVEP_SYMBOL\n")
            # Iterate positions using the simple parser (generator) above
            for chrom, pos in parse_positions(pos_file):
                matched = []
                try:
                    # Use vcf_in.fetch to iterate VCF records overlapping the 0-based interval [pos-1, pos)
                    # pysam.VariantFile.fetch returns VariantRecord objects
                    for rec in vcf_in.fetch(chrom, pos-1, pos):
                        if rec is None:
                            continue
                        # rec.pos is 1-based position
                        try:
                            rpos = int(rec.pos)
                        except Exception:
                            continue
                        if rpos != pos:
                            # in some complex VCFs fetch may return variants on same contig but different pos
                            continue
                        # capture core fields; rec.alts can be a tuple of ALT alleles
                        ref = rec.ref or ""
                        alt = ",".join(rec.alts) if rec.alts else ""
                        qual = str(rec.qual) if rec.qual is not None else ""
                        # Compose a simplified INFO string for output
                        info_items = []
                        for k in rec.info.keys():
                            v = rec.info.get(k)
                            # rec.info values can be lists/tuples for repeated keys; convert to string
                            if isinstance(v, (list, tuple)):
                                vstr = ",".join(str(x) for x in v)
                            else:
                                vstr = str(v)
                            info_items.append(f"{k}={vstr}")
                        info = ";".join(info_items)
                        matched.append((ref, alt, qual, info, rec))
                except ValueError:
                    # Could be invalid chromosome name for this VCF (not present in index); treat as no match
                    matched = []

                # Find gene names overlapping this coordinate via GTF fetch as a fallback / supplement
                genes = genes_at_pos(gtf_tb, chrom, pos)
                genes_str = ';'.join(genes) if genes else 'INTERGENIC'

                # Initialize VEP-derived fields to empty strings; they will remain empty if not present
                vep_consequence = ""
                vep_gene = ""
                vep_symbol = ""

                if matched:
                    # Choose the first matched VCF record as the representative record for that coordinate.
                    # Note: there may be multiple records (e.g., multi-allelic, multiple variant calls). You can
                    # choose to iterate ALL matched recs and write multiple lines per position if you need per-ALT output.
                    rec_obj = matched[0][4]
                    csq_entries = rec_obj.info.get("CSQ")
                    # CSQ may be a tuple/list of entries (one per alt or multiple annotations).
                    csq_entry = csq_entries[0] if csq_entries else ""
                    csq = parse_csq_entry(csq_entry, csq_fields) if csq_fields else {}
                    # Common VEP keys:
                    # - Consequence (string like 'missense_variant')
                    # - Gene (Ensembl gene ID)
                    # - SYMBOL (HGNC symbol or gene symbol)
                    vep_consequence = csq.get("Consequence", "")
                    vep_gene = csq.get("Gene", "") or csq.get("Gene_ID", "")
                    vep_symbol = csq.get("SYMBOL", "") or csq.get("Gene_symbol", "") or csq.get("Gene_Name", "")

                # Write output row. If no matched VCF record, 'PRESENT' is 'no' and REF/ALT/QUAL/INFO empty.
                if not matched:
                    out.write(f"{chrom}\t{pos}\tno\t\t\t\t\t{genes_str}\t{vep_consequence}\t{vep_gene}\t{vep_symbol}\n")
                else:
                    # If multiple matched records exist, we concatenated their string fields with semicolons
                    refs = ";".join(m[0] for m in matched)
                    alts = ";".join(m[1] for m in matched)
                    quals = ";".join(m[2] for m in matched)
                    infos = ";".join(m[3] for m in matched)
                    out.write(f"{chrom}\t{pos}\tyes\t{refs}\t{alts}\t{quals}\t{infos}\t{genes_str}\t{vep_consequence}\t{vep_gene}\t{vep_symbol}\n")
    except Exception as e:
        fail(f"Failed to open/write output file '{out_file}': {e}")

    print(f"Annotation complete. Output written to: {out_file}")

# ---------------------------------------------------------------------------------
# Entrypoint for embedded script
# ---------------------------------------------------------------------------------
def main():
    annotate_positions_file(POS_FILE, VEP_VCF, GTF_FILE, OUT_FILE)

if __name__ == "__main__":
    main()
PYCODE

# -------------------- END OF WRAPPER --------------------
#
# POST-RUN ADVICE & TROUBLESHOOTING
# --------------------------------
# 1) If pysam import fails in the interpreter test at top:
#    - Activate the conda environment where pysam is installed or install pysam:
#         mamba install -c bioconda pysam
#    - Confirm 'PYTHON' points to Python inside that env.
#
# 2) If the script reports missing .tbi files:
#    - Create an index with:
#         tabix -p vcf path/to/file.vcf.gz
#      or for GTF:
#         tabix -p gff file.gtf.gz
#
# 3) Multi-allelic sites:
#    - This script currently aggregates matched records and concatenates REF/ALT/QUAL/INFO strings.
#      If you need per-ALT granularity (for AFs or allele-specific CSQ fields), prefer re-writing the
#      Python section to emit one output row per ALT by iterating rec.alts and per-ALT CSQ entries.
#
# 4) CSQ parsing caveats:
#    - VEP CSQ header varies by version. The script extracts the "Format:" list from the VCF header.
#      If your header is non-standard, you may need to inspect its line manually:
#         zcat vep.vcf.gz | grep -m1 '##INFO=<ID=CSQ'
#
# 5) Performance:
#    - This approach streams positions from a text file and performs a random-access fetch on the VCF
#      and GTF for each position (tabix-backed). For many thousands of positions this is fine; for
#      hundreds of thousands consider batching or using bcftools view to pre-extract records, or
#      using cyvcf2 which supports fast random access and C-level parsing.
#
# 6) Suggested enhancements:
#    - Add command-line flags (argparse) to the embedded Python so you can call it directly:
#         python annotate_isec.py --pos-file positions.txt --vcf vep.vcf.gz --gtf genes.gtf.gz --out annotations.tsv
#    - Output intermediate RDS/pickle objects for reproducibility: save recs, parsed CSQ fields, etc.
#    - Provide an option to expand multi-allelic sites into one row per ALT with allele-specific annotations.
#    - Add unit-tests around CSQ parsing using a small example VCF snippet embedded in tests.
#
# If you'd like, I can:
#  - Rewrite the Python annotation body into a standalone script (annotate_positions.py) with CLI,
#    better multi-allelic handling, and unit tests.
#  - Convert CSQ parsing to handle VEP's allele-specific fields more precisely.
#
# End of file.
