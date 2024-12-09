#!/bin/bash

#Usage: ./setup_vOTUcluster.sh /path/to/db 4
# Check if the correct number of arguments is provided
# Create and activate the Conda environment
#mamba env create -f vOTUcluster
#mamba activate vOTUcluster

#clean cache
#conda clean --all

#instasll packages
#mamba install virsorter=2 --clobber
#mamba create -n vOTUcluster8.27 -c conda-forge -c bioconda "python=3.8" 
# 设置临时缓存目录，以防使用本地缓存
export CONDA_PKGS_DIRS=$(mktemp -d)

# 安装 dRep、ViralVerify、Genomad、CheckV
#pip3 install checkm-genome==1.2.2
echo "Installing dRep, ViralVerify, Genomad, and CheckV..."
mamba install -c conda-forge -c bioconda \
    genomad=1.8.0 --yes
mamba install -c conda-forge -c bioconda \
    checkm-genome=1.2.2 --yes
mamba install -c conda-forge -c bioconda \
    dRep=3.5.0 --yes
mamba install -c conda-forge -c bioconda \
    checkv=1.0.3 --yes

echo "Installing dRep, ViralVerify, Genomad, and CheckV..."
mamba install -c conda-forge -c bioconda viralverify=1.1 --yes

# 安装 scikit-learn、imbalanced-learn、pandas、seaborn、pyhmmer
echo "Installing scikit-learn, imbalanced-learn, pandas, seaborn, and pyhmmer..."
mamba install -c conda-forge -c bioconda \
    scikit-learn=0.22.1 imbalanced-learn pandas seaborn pyhmmer==0.10.14 --yes

# 安装 prodigal、screed、ruamel.yaml、snakemake、click、conda-package-handling
echo "Installing prodigal, screed, ruamel.yaml, snakemake, and related packages..."
mamba install -c conda-forge -c bioconda \
    prodigal screed ruamel.yaml "snakemake>=5.18,<=5.26" click "conda-package-handling<=1.9" --yes

# 安装 numpy、checkm-genome
echo "Installing numpy and checkm-genome..."
mamba install -c conda-forge -c bioconda numpy=1.23.5 --yes

# 安装 megahit、spades、fastp
echo "Installing megahit, spades, and fastp..."
mamba install -c conda-forge -c bioconda megahit spades fastp --yes
pip3 install bio

# 清理临时缓存目录
echo "Cleaning up temporary cache..."
rm -rf $CONDA_PKGS_DIRS


# Clone the repository and install the package
git clone https://github.com/liusihang/VirSorter2-pyhmmerAcc
cd VirSorter2-pyhmmerAcc
pip install -e .

# Run Python scripts
cd ..
python ./setupscript/Move2bin.py

# Get the current conda environment path
CONDA_ENV_PATH=$(conda info --base)/envs/$(basename "$CONDA_DEFAULT_ENV")

# Add execute permissions to all files in the current conda environment's bin folder
chmod +x "$CONDA_ENV_PATH/bin"

echo "All packages installed successfully!"
