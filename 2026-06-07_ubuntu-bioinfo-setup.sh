#!/usr/bin/env bash

set -e

echo "=== Installing Miniconda ==="
if [ ! -d "$HOME/miniconda3" ]; then
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
    bash miniconda.sh -b -p $HOME/miniconda3
    eval "$($HOME/miniconda3/bin/conda shell.bash hook)"
    conda init
else
    echo "Miniconda already installed."
    eval "$($HOME/miniconda3/bin/conda shell.bash hook)"
fi

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

echo "=== Installing R and IRkernel ==="
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

echo "=== Setup complete! ==="
echo "Activate environment with: conda activate bioinfo"
echo "Launch JupyterLab with: jupyter lab"
