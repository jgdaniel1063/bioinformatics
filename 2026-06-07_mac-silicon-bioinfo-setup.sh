#!/usr/bin/env bash

set -e

echo "=== Checking for Homebrew ==="
if ! command -v brew &> /dev/null; then
    echo "Homebrew not found. Installing..."
    /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
    echo 'eval "$(/opt/homebrew/bin/brew shellenv)"' >> ~/.zprofile
    eval "$(/opt/homebrew/bin/brew shellenv)"
else
    echo "Homebrew already installed."
fi

echo "=== Installing Miniconda (ARM64) ==="
brew install --cask miniconda
eval "$(/opt/homebrew/Caskroom/miniconda/base/bin/conda shell.zsh hook)"
conda init zsh

echo "=== Creating bioinformatics environment ==="
conda create -y -n bioinfo python=3.11
conda activate bioinfo

echo "=== Installing JupyterLab ==="
conda install -y -c conda-forge jupyterlab

echo "=== Installing Python scientific stack ==="
conda install -y -c conda-forge numpy pandas scipy matplotlib seaborn scikit-learn biopython

echo "=== Installing JupyterLab extensions ==="
pip install lckr-jupyterlab-variableinspector
conda install -y -c conda-forge jupyterlab-git
pip install jupyterlab_code_formatter black

echo "=== Installing R (ARM-native) and IRkernel ==="
conda install -y -c conda-forge r-base r-irkernel r-essentials

echo "=== Installing R bioinformatics packages ==="
Rscript - <<EOF
install.packages("BiocManager", repos="https://cloud.r-project.org")
BiocManager::install(c(
    "DESeq2",
    "tximport",
    "edgeR",
    "limma",
    "GenomicFeatures",
    "rtracklayer",
    "IsoformSwitchAnalyzeR"
))
IRkernel::installspec(user = TRUE)
EOF

echo "=== Installing core bioinformatics tools (ARM-native) ==="
conda install -y -c bioconda -c conda-forge \
    salmon \
    star \
    samtools \
    bcftools \
    bedtools \
    htslib \
    minimap2 \
    fastqc \
    multiqc \
    cutadapt \
    fastp \
    gffread \
    seqkit \
    hisat2 \
    bowtie2 \
    bwa \
    picard \
    subread \
    sra-tools

echo "=== Setup complete! ==="
echo "Activate environment with: conda activate bioinfo"
echo "Launch JupyterLab with: jupyter lab"
