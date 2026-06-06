#!/usr/bin/env bash

# Hard-coded SnpEff annotation script

# Input VCF (change this to your actual file)
INPUT="/home/jgd/data/my_variants.vcf.gz"

# Output VCF (change this to wherever you want it saved)
OUTPUT="/home/jgd/data/my_variants.annotated.vcf.gz"

# SnpEff settings
DB="GRCz11.115"
SNPEFF_DIR="/home/jgd/snpEff"

echo "Annotating $INPUT using SnpEff database $DB..."

# Annotate effects only
zcat "$INPUT" \
  | java -jar $SNPEFF_DIR/snpEff.jar $DB \
  | bgzip > "$OUTPUT"

# Index the output VCF
tabix -p vcf "$OUTPUT"

echo "Done. Output written to $OUTPUT"
