#!/bin/bash
#
# isec_annotation_heavily_annotated.sh
#
# Very heavily annotated wrapper that annotates genomic positions (intersection results)
# using an existing VEP-annotated VCF and a bgzipped/tabix-indexed GTF as fallback.
#
# This file contains:
#  - a shell wrapper that:
#      * verifies the chosen Python interpreter and pysam
#      * verifies VEP VCF + index and GTF + index are present
#      * invokes an embedded Python script (here-document) to perform the annotation work
#  - an embedded Python script (inside a here-doc) that:
#      * reads a position list (CHROM POS) line-by-line
#      * queries the VEP-annotated VCF for records at each coordinate and parses CSQ
#      * queries a tabix-indexed GTF for overlapping gene features if needed
#      * writes a summarized TSV for downstream use
#
# This file is intentionally verbose: every step, assumption, and common failure mode
# is explained inline so you can adapt the script safely to your environment.
#
# Edit only the HARD-CODED PATHS at the top of the file to point to your files, or
# convert the embedded Python into a standalone script that accepts CLI arguments.
#
# Run:
#   chmod +x isec_annotation_heavily_annotated.sh
#   ./isec_annotation_heavily_annotated.sh
#
# (If you prefer to run the Python portion directly, extract the here-doc block and run it
#  with your Python interpreter; the embedded script is self-contained.)
#
# -------------------------------------------------------------------------------------

# -------------------- SELECT PYTHON (edit if needed) ----------------------------
# Use an explicit interpreter so we run in the environment where pysam is installed.
# Replace this path with `python3` or another interpreter as needed.
PYTHON="/home/jgdaniel/miniconda3/envs/vep_env/bin/python"

# Quick sanity: ensure the interpreter exists and is executable.
if [ ! -x "${PYTHON}" ]; then
  echo "ERROR: interpreter not found or not executable: ${PYTHON}" >&2
  echo "  -> Update PYTHON variable near the top of this script to point to a valid python3 binary." >&2
  exit 1
fi

# -------------------- HARD-CODED PATHS (edit here only) -------------------------
# POS_FILE: text file with positions to annotate (CHROM POS per line).
# Example two-line file:
#   chr1    1234567
#   1       7654321     # supports either 'chr' prefix or not (script will not normalize)
#
POS_FILE="/home/jgd/Documents/bioinformatics_working/2025-11-06_enu-proc/intersection/2026-03-02_isec_with_ck-data/per_family/au01/lists/nothrom_specific_noBG.tsv"

# VEP-annotated VCF (bgzipped) already produced by your pipeline. The script relies on VEP's CSQ in INFO.
VEP_VCF="/home/jgd/Documents/bioinformatics_working/2025-11-06_enu-proc/vcf/gatk/gatk-ensembl/vep_annotated/gatk-ensembl_filtered.20260224.nodash.vep.vcf.gz"

# GTF (bgzipped) used as a fallback to determine overlapping genes at the coordinate.
# It must be tabix-indexed with -p gff (creates .tbi).
GTF_FILE="/home/jgd/Documents/bioinformatics_working/ref/danio_rerio.grcz11.115.gtf.gz"

# OUT_FILE: where to write annotated TSV. Parent directory will be created if missing.
OUT_FILE="/home/jgd/Documents/bioinformatics_working/output/isec_annotation/au01_annotated.tsv"

# Informational variables (not used in the script but useful for provenance)
VEP_SPECIES="danio_rerio"
VEP_ASSEMBLY="GRCz11"
VEP_FORKS="${SLURM_CPUS_PER_TASK:-$(nproc)}"
# -------------------------------------------------------------------------------

# Print which python we will attempt to run (handy for debugging path issues)
echo "Using python: ${PYTHON}"

# Run a tiny Python snippet with the selected interpreter to verify that pysam is importable
# and that the interpreter is the one you expect. The snippet prints the Python executable,
# version and pysam version. If this fails, the embedded script will also fail later.
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

# Verify VEP VCF and tabix index exist — we require both .vcf.gz and .vcf.gz.tbi
if [ -f "${VEP_VCF}" ] && [ -f "${VEP_VCF}.tbi" ]; then
  echo "Found existing VEP-annotated VCF and index: ${VEP_VCF} ${VEP_VCF}.tbi"
else
  echo "ERROR: VEP-annotated VCF or index not found: ${VEP_VCF} or ${VEP_VCF}.tbi" >&2
  echo "  -> Create the index with: tabix -p vcf ${VEP_VCF}" >&2
  exit 1
fi

# Verify GTF and index exist
if [ -f "${GTF_FILE}" ] && [ -f "${GTF_FILE}.tbi" ]; then
  echo "Found GTF and index: ${GTF_FILE} ${GTF_FILE}.tbi"
else
  echo "ERROR: GTF or GTF.tbi not found: ${GTF_FILE} or ${GTF_FILE}.tbi" >&2
  echo "  -> Create index with: tabix -p gff ${GTF_FILE}" >&2
  exit 1
fi

# -------------------------------------------------------------------------------------
# Embedded Python script: performs the annotation logic using pysam.
# The script is delimited by the here-doc; it is executed by the interpreter above.
# It is heavily commented so you can extract it and run it separately if preferred.
# -------------------------------------------------------------------------------------
"${PYTHON}" - "$@" <<'PYCODE'
#!/usr/bin/env python3
"""
Embedded annotator (very heavily documented).

Purpose
-------
Given:
 - POS_FILE: text file with coordinates (CHROM POS)
 - VEP_VCF: a bgzipped VCF annotated with VEP (CSQ in INFO)
 - GTF_FILE: bgzipped GTF with tabix index (.tbi)

Produce:
 - OUT_FILE: TSV summarizing VEP info and overlapping gene names for each input position.

Design considerations
---------------------
- Primary source of functional consequence / gene mapping is the VEP CSQ field in the VCF INFO.
  The script parses the VCF header to learn CSQ's pipe-delimited field order, then maps entries correctly.
- When VEP info is absent or insufficient, we consult the GTF (tabix) for overlapping gene features,
  extract gene_name/gene_id/transcript_id attributes, and return them as a fallback.
- The script writes one output row per input coordinate. It consolidates multiple VCF records at the
  same coordinate by concatenating simple text fields (REF/ALT/QUAL/INFO). If you need one row per ALT
  (per-allele), consider expanding matched VCF records into multiple output rows — instructions below.
- The script uses pysam for safe, indexed access to both VCF and GTF files (fast and robust).

Failure modes
-------------
- pysam not installed in the selected interpreter -> import error (the wrapper detects this earlier).
- VCF or GTF index (.tbi) missing -> script exits with an instructive message.
- CSQ header not present or non-standard -> CSQ parsing will be best-effort and may return empty values.
- Multi-allelic sites: current behavior concats REF/ALT values; not allele-specific.

Improvements you might want
---------------------------
- Convert this embedded script to a standalone CLI with argparse for flexibility.
- Add explicit multi-allelic handling: expand each ALT into its own output row and select the matching CSQ entry per allele.
- Output compressed TSV (.tsv.gz) for large outputs.
- Add more robust parsing for GTF attributes (quoted strings with embedded semicolons).
"""

# Standard imports
import os
import sys
import re

# pysam: required for VariantFile and TabixFile access
try:
    import pysam
except Exception as e:
    sys.stderr.write("pysam import failed in embedded script: {}\n".format(e))
    raise

# === The constants below are intentionally the same as in the wrapper for clarity ===
POS_FILE = "/home/jgd/Documents/bioinformatics_working/2025-11-06_enu-proc/intersection/2026-03-02_isec_with_ck-data/per_family/ai08/lists/nothrom_specific_noBG.tsv"
VEP_VCF = "/home/jgd/Documents/bioinformatics_working/2025-11-06_enu-proc/vcf/gatk/gatk-ensembl/vep_annotated/gatk-ensembl_filtered.20260224.nodash.vep.vcf.gz"
GTF_FILE = "/home/jgd/Documents/bioinformatics_working/ref/danio_rerio.grcz11.115.gtf.gz"
OUT_FILE = "/home/jgd/Documents/bioinformatics_working/output/isec_annotation/ai08_annotated.tsv"
# ====================================================================================

# Utility: exit with message
def fail(msg, rc=1):
    sys.stderr.write(msg + "\n")
    sys.exit(rc)

# -----------------------
# Parsing helpers
# -----------------------
def parse_attrs(attr_string):
    """
    Minimal parser for the GTF/GFF attributes field (column 9).
    Example attribute string:
      gene_id "ENSG000001234"; gene_name "BRCA1"; transcript_id "ENST00000380152";
    This function returns a dict: {'gene_id': 'ENSG000001234', 'gene_name':'BRCA1', ...}
    Caveats:
      - Not a full GTF parser: it assumes keys and quoted values separated by whitespace,
        and attributes separated by semicolons.
      - If your GTF uses slightly different conventions, adapt as necessary.
    """
    d = {}
    # split on semicolon; for each part, split on the first whitespace to get key and value
    for part in attr_string.strip().rstrip(';').split(';'):
        part = part.strip()
        if not part:
            continue
        if ' ' in part:
            k, v = part.split(' ', 1)
            d[k] = v.strip().strip('"')
    return d

def parse_positions(path):
    """
    Generator that yields (chrom, pos) from POS_FILE.
    Accepts rows with at least two whitespace-separated columns (chrom and pos).
    Skips malformed lines with a warning.
    """
    with open(path, 'r') as fh:
        for lineno, line in enumerate(fh, 1):
            s = line.strip()
            if not s:
                continue
            cols = s.split()
            if len(cols) < 2:
                sys.stderr.write(f"Warning: Skipping malformed line {lineno}: {s}\n")
                continue
            chrom = cols[0]
            try:
                pos = int(cols[1])
            except ValueError:
                sys.stderr.write(f"Warning: Skipping non-integer pos at line {lineno}: {s}\n")
                continue
            yield chrom, pos

def genes_at_pos(gtf_tb, chrom, pos):
    """
    Query the tabix-indexed GTF for features overlapping the given coordinate.
    - pysam.TabixFile.fetch uses 0-based start and end positions (start inclusive, end exclusive).
    - We request the interval [pos-1, pos) to get features overlapping the 1-based position pos.
    - From each GTF line we parse column 9 attributes and extract gene_name, gene_id, or transcript_id.
    - Returns sorted list of unique names (can be empty).
    """
    genes = set()
    try:
        for rec in gtf_tb.fetch(chrom, pos-1, pos):
            if not rec or rec.startswith('#'):
                continue
            cols = rec.split('\t')
            if len(cols) < 9:
                continue
            ad = parse_attrs(cols[8])
            if 'gene_name' in ad:
                genes.add(ad['gene_name'])
            elif 'gene_id' in ad:
                genes.add(ad['gene_id'])
            elif 'transcript_id' in ad:
                genes.add(ad['transcript_id'])
    except ValueError:
        # chromosome not present in the index (common if naming mismatch like '1' vs 'chr1')
        pass
    return sorted(genes)

def ensure_outdir(out_path):
    """
    Create parent directory for the output file and check we can write into it.
    Writes a short temporary test file to ensure filesystem permissions are sufficient.
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

# -----------------------
# VEP CSQ parsing helpers
# -----------------------
def get_vep_csq_fields_from_header(vcf_variantfile):
    """
    Read the VCF header text and attempt to extract the 'Format: ...' field for the CSQ INFO.
    VEP adds a header line similar to:
      ##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. Format: Allele|Consequence|IMPACT|SYMBOL|...">
    This function returns the list of CSQ field names in order, or None if not found.
    """
    header_text = str(vcf_variantfile.header)
    for line in header_text.splitlines():
        if line.startswith("##INFO=<ID=CSQ"):
            m = re.search(r'Format:\s*([^">]+)', line)
            if m:
                fields = [f.strip() for f in m.group(1).split('|')]
                return fields
    return None

def parse_csq_entry(entry_str, csq_fields):
    """
    Parse a single CSQ string (pipe-delimited) into a dict using csq_fields as keys.
    Ensures we pad the parts if some trailing fields are missing.
    If csq_fields is None, returns empty dict.
    """
    if not entry_str or not csq_fields:
        return {}
    parts = entry_str.split('|')
    if len(parts) < len(csq_fields):
        parts += [''] * (len(csq_fields) - len(parts))
    return dict(zip(csq_fields, parts))

# -----------------------
# Main annotation function
# -----------------------
def annotate_positions_file(pos_file, vep_vcf, gtf_file, out_file):
    """
    Orchestrate everything:
      1. Validate inputs and indexes
      2. Open VCF and GTF using pysam
      3. For each position in pos_file:
           - fetch VCF records at that position (if any)
           - parse REF/ALT/QUAL and reconstruct a simplified INFO string
           - parse CSQ entries if present (map to known CSQ keys parsed from header)
           - query the GTF for overlapping genes (fallback)
           - write one TSV row summarizing results
    """
    ensure_outdir(out_file)

    # Sanity checks for input files and their tabix indexes:
    if not os.path.isfile(pos_file):
        fail(f"Positions file not found: {pos_file}")
    if not os.path.isfile(vep_vcf):
        fail(f"VEP-annotated VCF file not found: {vep_vcf}")
    if not os.path.isfile(vep_vcf + ".tbi"):
        fail(f"VEP-annotated VCF index not found: {vep_vcf}.tbi (create with tabix -p vcf)")
    if not os.path.isfile(gtf_file):
        fail(f"GTF file not found: {gtf_file}")
    if not os.path.isfile(gtf_file + ".tbi"):
        fail(f"GTF index not found: {gtf_file}.tbi (create with tabix -p gff)")

    # Open VCF and GTF (these operations validate indexes)
    try:
        vcf_in = pysam.VariantFile(vep_vcf)
    except Exception as e:
        fail(f"Failed to open VEP VCF/index '{vep_vcf}': {e}")

    try:
        gtf_tb = pysam.TabixFile(gtf_file)
    except Exception as e:
        fail(f"Failed to open GTF/index '{gtf_file}': {e}")

    # Parse CSQ header to know the order of the pipe-delimited fields
    csq_fields = get_vep_csq_fields_from_header(vcf_in)
    if csq_fields is None:
        # Warn but continue: we can still report VCF presence and fallback to GTF
        sys.stderr.write("Warning: CSQ header not found in VCF. VEP fields will not be parsed.\n")

    # Open output and write header row
    try:
        with open(out_file, "w") as out:
            out.write("CHROM\tPOS\tPRESENT\tREF\tALT\tQUAL\tINFO\tGENES\tVEP_CONSEQUENCE\tVEP_GENE\tVEP_SYMBOL\n")
            for chrom, pos in parse_positions(pos_file):
                # Accumulate matched VCF records at this coordinate
                matched = []
                try:
                    # VariantFile.fetch uses 0-based start, end-exclusive ranges
                    for rec in vcf_in.fetch(chrom, pos-1, pos):
                        if rec is None:
                            continue
                        try:
                            rpos = int(rec.pos)
                        except Exception:
                            # If rec.pos is not parseable, skip the record
                            continue
                        # Ensure we only keep records exactly at the coordinate (some fetch results may be nearby)
                        if rpos != pos:
                            continue
                        ref = rec.ref or ""
                        alt = ",".join(rec.alts) if rec.alts else ""
                        qual = str(rec.qual) if rec.qual is not None else ""
                        # Reconstruct a compact INFO string: key=value (join list entries by comma)
                        info_items = []
                        for k in rec.info.keys():
                            v = rec.info.get(k)
                            if isinstance(v, (list, tuple)):
                                vstr = ",".join(str(x) for x in v)
                            else:
                                vstr = str(v)
                            info_items.append(f"{k}={vstr}")
                        info = ";".join(info_items)
                        matched.append((ref, alt, qual, info, rec))
                except ValueError:
                    # fetch may throw ValueError if chromosome is absent in VCF index
                    matched = []

                # GTF fallback: fetch overlapping genes at this coordinate
                genes = genes_at_pos(gtf_tb, chrom, pos)
                genes_str = ';'.join(genes) if genes else 'INTERGENIC'

                # Default empty VEP outputs
                vep_consequence = ""
                vep_gene = ""
                vep_symbol = ""

                if matched:
                    # Use the first matched VCF record as representative (simple behavior).
                    # For multi-allelic or multiple records you may want to expand rows (see below).
                    rec_obj = matched[0][4]
                    csq_entries = rec_obj.info.get("CSQ")
                    csq_entry = csq_entries[0] if csq_entries else ""
                    csq = parse_csq_entry(csq_entry, csq_fields) if csq_fields else {}
                    # Common CSQ keys: Consequence, Gene, SYMBOL, etc.
                    vep_consequence = csq.get("Consequence", "")
                    vep_gene = csq.get("Gene", "") or csq.get("Gene_ID", "")
                    vep_symbol = csq.get("SYMBOL", "") or csq.get("Gene_symbol", "") or csq.get("Gene_Name", "")

                # Compose line(s) and write. If matched is empty, mark PRESENT=no.
                if not matched:
                    out.write(f"{chrom}\t{pos}\tno\t\t\t\t\t{genes_str}\t{vep_consequence}\t{vep_gene}\t{vep_symbol}\n")
                else:
                    # If multiple matched records are present (rare), we concatenate simple fields with semicolons.
                    # This preserves some information but is not allele-disambiguated.
                    refs = ";".join(m[0] for m in matched)
                    alts = ";".join(m[1] for m in matched)
                    quals = ";".join(m[2] for m in matched)
                    infos = ";".join(m[3] for m in matched)
                    out.write(f"{chrom}\t{pos}\tyes\t{refs}\t{alts}\t{quals}\t{infos}\t{genes_str}\t{vep_consequence}\t{vep_gene}\t{vep_symbol}\n")
    except Exception as e:
        fail(f"Failed to open/write output file '{out_file}': {e}")

    print(f"Annotation complete. Output written to: {out_file}")

# Entrypoint
def main():
    annotate_positions_file(POS_FILE, VEP_VCF, GTF_FILE, OUT_FILE)

if __name__ == "__main__":
    main()
PYCODE

# -------------------------- EXTENSIVE ADDITIONAL NOTES --------------------------
#
# You requested "heavily annotate" — the wrapper and the embedded Python above include
# many comments. Here are extra, pragmatic recommendations and examples you will find
# useful when running/adapting this in the real world.
#
# 1) Common problems & fixes
#    - pysam import error
#        * Activate the conda environment where pysam is installed: `conda activate vep_env`
#        * Install pysam: `mamba install -c bioconda pysam` or `pip install pysam`
#        * Then point PYTHON to that environment's python binary.
#
#    - "VCF index not found" or "GTF index not found"
#        * Index a bgzipped VCF: `tabix -p vcf myfile.vcf.gz`
#        * Index a bgzipped GTF (GFF-like): `tabix -p gff myfile.gtf.gz`
#        * Confirm the index file location: index should be sibling file myfile.vcf.gz.tbi
#
#    - No CSQ in VCF header or missing CSQ content
#        * Inspect header: `zgrep '^##INFO=<ID=CSQ' my.vep.vcf.gz | sed -n '1p'`
#        * If header is missing, VEP may not have been run or the VCF may have CSQ under a different INFO tag.
#        * If CSQ present but entries are empty, check your VEP run or re-run VEP with the appropriate options.
#
#    - Chromosome naming mismatch between files (e.g., '1' vs 'chr1')
#        * If your VCF uses '1' and GTF uses 'chr1' you'll get no GTF overlaps. Normalize names ahead of time
#          (e.g., in the VCF using `bcftools annotate --rename-chrs` or by altering the GTF).
#
# 2) Multi-allelic sites — how to improve per-ALT annotation
#    - VCF can contain records with ALTs like "A,C" and AD entries like "10,3,2".
#    - The current script concatenates alt alleles and CSQ values and picks the first CSQ entry when present.
#    - If you need allele-specific annotations:
#       * For each VariantRecord rec in vcf_in.fetch(...):
#           for alt_index, alt in enumerate(rec.alts):
#               # AD array: rec.format('AD')[sample_index] or rec.samples[...]['AD'] if you parse samples
#               # map csq_entries to the matching allele when VEP outputs allele-specific CSQ (often prefixed by the allele)
#               # Write one output row per allele with the matching CSQ sub-entry.
#       * This is more code but produces unambiguous allele-level results.
#
# 3) Example VEP CSQ header line (for parsing)
#    - Typical header line looks like:
#      ##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from VEP. Format: Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|...">
#    - The function get_vep_csq_fields_from_header extracts the pipe-separated "Format: ..." fields.
#
# 4) About the GTF parsing
#    - parse_attrs() handles the common `key "value";` pattern, which is used by Ensembl GTFs.
#    - If your GTF has attributes separated differently (rare), adjust parse_attrs accordingly.
#    - For robust parsing of many GTF/GFF flavors, consider using libraries like `gffutils` or `pysam`’s built-in GFF parsing.
#
# 5) Performance & scaling
#    - For small numbers of positions (hundreds to tens of thousands) this approach is fine.
#    - For hundreds of thousands to millions of positions consider:
#        * Using bcftools to pre-extract the variants in a region or of interest:
#            bcftools view -R pos_file.txt -Oz input.vep.vcf.gz -o subset.vcf.gz
#            bcftools index subset.vcf.gz
#          Then parse subset.vcf.gz in one pass.
#        * Using cyvcf2 to vectorize VCF scanning; cyvcf2 is very fast for large queries and can yield allele arrays directly.
#        * Parallelizing across chromosomes or position chunks (careful to merge outputs).
#
# 6) Testing
#    - Create a tiny test VCF with one or two records and a tiny GTF; run the script on that test set to validate behavior.
#    - Example minimal VCF/GTF and expected output can be used for unit tests.
#
# 7) Reproducibility
#    - Record the pysam version, python version, and the VEP version used to generate the VCF.
#    - Optionally write `session` information or a small YAML manifest documenting input file checksums.
#
# 8) Converting to a CLI tool (recommended)
#    - Extract the embedded Python into `annotate_positions.py`, add argparse arguments:
#         --pos-file, --vcf, --gtf, --out, --expand-multiallelic, --csq-fields, --threads
#    - This makes the tool reusable and easier to test/CI.
#
# 9) Security note
#    - If the POS_FILE content comes from external sources, validate or sanitize it before using.
#
# 10) If you would like more:
#    - I can rewrite the Python body as a standalone script with:
#         * per-allele output option,
#         * better CSQ-to-allele matching,
#         * optional sample-level genotype/AD extraction,
#         * optional output compression and a summary report.
#    - Tell me which features you want and I will produce the new script.
#
# -------------------------------------------------------------------------------------

# End of heavily annotated wrapper
