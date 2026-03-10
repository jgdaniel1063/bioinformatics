#!/usr/bin/env python3
"""
tx2gene_generate_hardcoded.py

Hard-coded helper to produce a tx2gene mapping file (target_id -> gene_id).
This version removes CLI argument parsing and uses top-of-file placeholders you
should edit to point to your local files.

Usage:
  - Edit the three placeholder variables below:
      MODE = "GTF" or "FASTA"
      GTF_PATH = "/path/to/annotations.gtf.gz"
      FASTA_PATH = "/path/to/transcripts.fa.gz"
      OUT_PATH = "/path/to/tx2gene.tsv"
  - Then run:
      python tx2gene_generate_hardcoded.py

Notes:
  - Recommended: use MODE="GTF" and point GTF_PATH to a GTF/GFF3 (bgzipped or plain).
  - If only a transcript FASTA is available, set MODE="FASTA" and set FASTA_PATH.
  - The script removes transcript version suffixes by default (ENSEMBL-like ".1").
  - Output is a two-column TSV with header: target_id<TAB>gene_id

This file is intentionally verbose and annotated to explain assumptions, failure
modes, and how to adapt heuristics to non-standard file formats.

Author: Generated for you (hard-coded placeholders)
Date: 2026-03-10
"""
from __future__ import annotations
import gzip
import io
import os
import re
import sys
from collections import defaultdict
from typing import Dict, Tuple, Optional

# ============================
# USER: HARD-CODED PLACEHOLDERS
# Edit these paths to match your environment. Do NOT leave them pointing to
# non-existent locations when you run the script.
# ============================

# MODE: "GTF" prefers parsing a GTF/GFF3 file (recommended). "FASTA" will parse transcript FASTA headers heuristically.
MODE = "GTF"    # choose "GTF" or "FASTA"

# When MODE == "GTF", set GTF_PATH to your GTF/GFF3 file (can be gzipped).
GTF_PATH = "/path/to/annotations.gtf.gz"

# When MODE == "FASTA", set FASTA_PATH to your transcriptome FASTA (can be gzipped).
FASTA_PATH = "/path/to/transcripts.fa.gz"

# Output tx2gene mapping file (two-column TSV). Parent directory will be created.
OUT_PATH = "/path/to/tx2gene.tsv"

# Which attribute to prefer from GTF: "gene_id" (typical) or "gene_name" (human-readable symbol)
PREFER_GENE_FIELD = "gene_id"

# Remove transcript version suffix like ".1" from transcript IDs (True recommended for Ensembl)
REMOVE_VERSION = True

# If using FASTA, you can provide a custom regex to extract gene name from header
# Example: r'gene:([A-Za-z0-9_.:-]+)' (use a single capturing group)
FASTA_GENE_REGEX: Optional[str] = None  # set to string or None

# Minimum number of mappings required before script will write output (safety)
MIN_ROWS = 1

# ============================
# End of user-editable placeholders
# ============================

# ----------------------------
# Helpers: file opening (supports .gz)
# ----------------------------
def smart_open_read(path: str):
    """Open plain text or gzipped file for reading in text mode."""
    if path.endswith(".gz"):
        # io.TextIOWrapper ensures text mode with proper decoding
        return io.TextIOWrapper(gzip.open(path, "rb"))
    return open(path, "r", encoding="utf-8")

# ----------------------------
# GTF attribute parsing
# ----------------------------
_attr_kv_re = re.compile(r'\s*([^ "\t]+)\s+"([^"]+)"\s*;?')  # matches key "value";

def parse_gtf_attrs(attr_text: str) -> Dict[str,str]:
    """
    Parse GTF/GFF3-style attributes string into a dict.
    Supports typical 'key "value"; key2 "value2";' formats and falls back to key=value.
    """
    d: Dict[str,str] = {}
    for m in _attr_kv_re.finditer(attr_text):
        k, v = m.group(1), m.group(2)
        d[k] = v
    if d:
        return d
    # fallback (GFF3 style key=value;key2=value2)
    for part in re.split(r';\s*', attr_text.strip().rstrip(';')):
        if not part:
            continue
        if '=' in part:
            k, v = part.split('=', 1)
            d[k.strip()] = v.strip().strip('"')
    return d

# ----------------------------
# Extract tx->gene mapping from GTF
# ----------------------------
def tx2gene_from_gtf(gtf_path: str, prefer_gene_field: str = "gene_id", remove_version: bool = True
                    ) -> Tuple[Dict[str,str], Dict[str, set]]:
    """
    Parse a GTF and return mapping transcript_id -> gene_id (or gene_name if prefer specified).
    Also return conflicts dict for diagnostics: transcript_id -> set(gene_ids) when multiple observed.
    Notes:
      - Prefers lines with 'transcript' in column 3, but accepts any line with transcript_id attr.
      - remove_version: strip trailing .<number> from transcript IDs (common for Ensembl).
    """
    mapping: Dict[str,str] = {}
    conflicts_by_tid: Dict[str, set] = defaultdict(set)
    count_lines = 0
    count_mapped = 0
    if not os.path.exists(gtf_path):
        raise FileNotFoundError(f"GTF not found: {gtf_path}")
    with smart_open_read(gtf_path) as fh:
        for ln in fh:
            if not ln or ln.startswith('#'):
                continue
            cols = ln.rstrip("\n").split('\t')
            if len(cols) < 9:
                continue
            count_lines += 1
            feature = cols[2]
            attrs = parse_gtf_attrs(cols[8])
            # Common attribute keys
            tid = attrs.get("transcript_id") or attrs.get("transcript") or attrs.get("ID")
            gid = attrs.get(prefer_gene_field) or attrs.get("gene_id") or attrs.get("gene_name")
            if tid is None or gid is None:
                continue
            if remove_version:
                tid = re.sub(r"\.\d+$", "", tid)
            # track mapping and conflicts
            if tid in mapping:
                if mapping[tid] != gid:
                    conflicts_by_tid[tid].add(mapping[tid])
                    conflicts_by_tid[tid].add(gid)
            else:
                mapping[tid] = gid
                count_mapped += 1
    # Build a summary conflicts dict keyed by tid (for human inspection)
    real_conflicts = {tid: gids for tid, gids in conflicts_by_tid.items() if len(gids) > 1}
    sys.stderr.write(f"[INFO] Parsed GTF lines: {count_lines}, unique transcript mappings: {count_mapped}\n")
    if real_conflicts:
        sys.stderr.write(f"[WARNING] {len(real_conflicts)} transcripts had multiple gene mappings (examples below):\n")
        for i, (tid, gids) in enumerate(real_conflicts.items()):
            if i >= 10:
                break
            sys.stderr.write(f"  {tid} -> {sorted(gids)}\n")
    return mapping, real_conflicts

# ----------------------------
# Extract tx->gene mapping heuristically from FASTA headers
# ----------------------------
def tx2gene_from_fasta(fa_path: str, fasta_gene_regex: Optional[str] = None, remove_version: bool = True) -> Dict[str,str]:
    """
    Best-effort extraction of transcript_id -> gene_id from FASTA headers.
    Heuristics:
      - transcript id: first token after '>'
      - gene id: search header text for patterns (common examples)
      - if fasta_gene_regex is provided, it is tried first (must have exactly one capturing group)
    Limitations:
      - FASTA headers are not standardized; this may miss many mappings. GTF is preferred.
    """
    if not os.path.exists(fa_path):
        raise FileNotFoundError(f"FASTA not found: {fa_path}")
    mapping: Dict[str,str] = {}
    patterns = []
    if fasta_gene_regex:
        patterns.append(re.compile(fasta_gene_regex))
    # common heuristics
    patterns += [
        re.compile(r'gene:([A-Za-z0-9_.:-]+)'),
        re.compile(r'gene_id:?([A-Za-z0-9_.:-]+)'),
        re.compile(r'gene_name:?([A-Za-z0-9_.:-]+)'),
        re.compile(r'GN:([A-Za-z0-9_.:-]+)'),
        re.compile(r'gene=([A-Za-z0-9_.:-]+)')
    ]
    examples_shown = 0
    with smart_open_read(fa_path) as fh:
        for ln in fh:
            if not ln.startswith('>'):
                continue
            header = ln[1:].strip()
            tid = header.split()[0]
            if remove_version:
                tid = re.sub(r'\.\d+$', '', tid)
            gid = None
            for pat in patterns:
                m = pat.search(header)
                if m:
                    # choose last capturing group if pattern has multiple
                    gid = m.group(m.lastindex)
                    break
            if gid:
                mapping[tid] = gid
            else:
                # print a few examples of headers that failed to parse to help user craft a regex
                if examples_shown < 5:
                    sys.stderr.write(f"[NOTE] Could not extract gene id from FASTA header: '{header[:140]}'\n")
                    examples_shown += 1
    sys.stderr.write(f"[INFO] Parsed FASTA: {len(mapping)} transcript->gene mappings inferred (heuristic)\n")
    return mapping

# ----------------------------
# Write tx2gene TSV
# ----------------------------
def write_tx2gene_tsv(mapping: Dict[str,str], out_path: str, include_header: bool = True):
    """
    Write a two-column TSV of target_id\tgene_id. Parent directory is created if needed.
    """
    parent = os.path.dirname(os.path.abspath(out_path)) or "."
    os.makedirs(parent, exist_ok=True)
    with open(out_path, "w", encoding="utf-8") as oh:
        if include_header:
            oh.write("target_id\tgene_id\n")
        for tid, gid in mapping.items():
            oh.write(f"{tid}\t{gid}\n")
    sys.stderr.write(f"[INFO] Wrote tx2gene mapping with {len(mapping)} rows to: {out_path}\n")

# ----------------------------
# Main routine (hard-coded flow)
# ----------------------------
def main():
    print(f"[RUN] MODE={MODE}; GTF_PATH={GTF_PATH}; FASTA_PATH={FASTA_PATH}; OUT_PATH={OUT_PATH}", file=sys.stderr)

    mapping: Dict[str,str] = {}
    conflicts = {}
    if MODE.upper() == "GTF":
        # Preferred path: parse the GTF
        if not os.path.exists(GTF_PATH):
            sys.exit(f"[ERROR] GTF not found (change GTF_PATH placeholder): {GTF_PATH}")
        mapping, conflicts = tx2gene_from_gtf(GTF_PATH, prefer_gene_field=PREFER_GENE_FIELD, remove_version=REMOVE_VERSION)
    elif MODE.upper() == "FASTA":
        # Fallback: parse FASTA heuristically
        if not os.path.exists(FASTA_PATH):
            sys.exit(f"[ERROR] FASTA not found (change FASTA_PATH placeholder): {FASTA_PATH}")
        mapping = tx2gene_from_fasta(FASTA_PATH, fasta_gene_regex=FASTA_GENE_REGEX, remove_version=REMOVE_VERSION)
    else:
        sys.exit("[ERROR] MODE must be either 'GTF' or 'FASTA' — edit the script top placeholders accordingly.")

    if len(mapping) < MIN_ROWS:
        sys.exit(f"[ERROR] Generated mapping has only {len(mapping)} rows which is less than MIN_ROWS={MIN_ROWS}. Aborting.")

    # Optionally warn about mapping density; useful for FASTA heuristics
    sys.stderr.write(f"[INFO] Total mappings generated: {len(mapping)}\n")
    if conflicts:
        sys.stderr.write(f"[WARNING] Conflicts detected — see earlier messages. Affected transcripts: {len(conflicts)}\n")

    # Write out the tx2gene file
    write_tx2gene_tsv(mapping, OUT_PATH, include_header=True)

    # Print a small preview on stdout for convenience
    print("# Example rows (first 10):")
    i = 0
    for tid, gid in mapping.items():
        print(tid + "\t" + gid)
        i += 1
        if i >= 10:
            break
    print("# Done. Use the generated tx2gene in sleuth (target_mapping).", file=sys.stderr)

if __name__ == "__main__":
    main()
