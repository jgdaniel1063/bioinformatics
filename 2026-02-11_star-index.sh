#!/usr/bin/env bash
# build_star_index.sh - Simple STAR index builder
# Usage: ./build_star_index.sh
set -euo pipefail

# Edit these
GENOME_DIR="/media/jgd/jd-1tb/bioinformatics/reference/star_index"
REF_FA="/home/jgd/Documents/bioinformatics_working/reference/danRer11.fa.gz"
GTF_FILE="/home/jgd/Documents/bioinformatics_working/reference/danRer11.ensGene.gtf"  # optional, set to "" if none
SJDB_OVERHANG=149
THREADS=6

ts() { date -u +%Y-%m-%dT%H:%M:%SZ; }

if ! command -v STAR >/dev/null 2>&1; then
  echo "[$(ts)] ERROR: STAR not in PATH"
  exit 1
fi

if [[ -z "$REF_FA" ]]; then
  echo "[$(ts)] ERROR: Set REF_FA"
  exit 1
fi

mkdir -p "$GENOME_DIR"

# Uncompress if needed
REF_TO_USE="$REF_FA"
if [[ "$REF_FA" == *.gz ]]; then
  REF_TO_USE="$(mktemp --tmpdir star_ref_XXXXXX.fa)"
  echo "[$(ts)] Uncompressing $REF_FA -> $REF_TO_USE"
  gzip -dc "$REF_FA" > "$REF_TO_USE"
  if [[ ! -s "$REF_TO_USE" ]]; then
    echo "[$(ts)] ERROR: Uncompressed file is empty"
    rm -f "$REF_TO_USE"
    exit 1
  fi
  trap "rm -f '$REF_TO_USE'" EXIT
fi

# Build index
CMD=(STAR --runThreadN "$THREADS" --runMode genomeGenerate --genomeDir "$GENOME_DIR" --genomeFastaFiles "$REF_TO_USE" --sjdbOverhang "$SJDB_OVERHANG")
if [[ -n "$GTF_FILE" && -f "$GTF_FILE" ]]; then
  CMD+=(--sjdbGTFfile "$GTF_FILE")
fi

echo "[$(ts)] Running: ${CMD[*]}"
"${CMD[@]}"
echo "[$(ts)] Index built in $GENOME_DIR"