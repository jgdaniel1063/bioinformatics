#!/bin/bash

# Simple STAR paired-end alignment for all files in input folder
# Memory-limited to prevent crashes

ulimit -v 40000000  # Increased to 40GB

INPUT_DIR="/media/jgd/jd-1tb/bioinformatics/2024-04-01_proc-enu_paper/trimming"
GENOME_DIR="/media/jgd/jd-1tb/bioinformatics/reference/star_index"
OUTPUT_DIR="/home/jgd/Documents/bioinformatics_working/output"
GTF_FILE="/home/jgd/Documents/bioinformatics_working/reference/danRer11.ensGene.gtf"

mkdir -p "$OUTPUT_DIR"

for r1 in "$INPUT_DIR"/*_r1_R1.fastq.gz; do
    [ -f "$r1" ] || continue
    r2="${r1/_r1_R1.fastq.gz/_r1_R2.fastq.gz}"
    [ -f "$r2" ] || { echo "Skipping $r1 - no pair"; continue; }
    
    sample=$(basename "$r1" | sed 's/_r1_R1\.fastq\.gz//')
    out_prefix="$OUTPUT_DIR/${sample}_aligned"
    
    echo "Aligning $sample"
    STAR --genomeDir "$GENOME_DIR" \
         --readFilesIn "$r1" "$r2" \
         --readFilesCommand zcat \
         --runThreadN 2 \
         --outFilterMultimapNmax 1 \
         --outSAMmapqUnique 60 \
         --outSAMtype BAM SortedByCoordinate \
         --outSAMunmapped Within \
         --outSAMattributes NH HI AS nM MD XS \
         --outSAMattrRGline "ID: SM: LB:_lib PL:ILLUMINA PU:.run" \
         --sjdbGTFfile "$GTF_FILE" \
         --outFileNamePrefix "$out_prefix/" \
         --limitBAMsortRAM 10000000000 \
         --genomeLoad NoSharedMemory
done

echo "Done"