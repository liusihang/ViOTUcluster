#!/usr/bin/env bash

# Set a temporary cache directory to avoid using the local cache
export CONDA_PKGS_DIRS=$(mktemp -d)

# Install dRep, ViralVerify, Genomad, and CheckV
echo "Installing dRep, ViralVerify, Genomad, and CheckV..."
mamba install -c conda-forge -c bioconda genomad=1.8.0 --yes
mamba install -c conda-forge -c bioconda checkm-genome=1.2.2 --yes
mamba install -c conda-forge -c bioconda dRep=3.5.0 --yes
mamba install -c conda-forge -c bioconda checkv=1.0.3 --yes

echo "Installing ViralVerify..."
mamba install -c conda-forge -c bioconda viralverify=1.1 --yes

# Install scikit-learn, imbalanced-learn, pandas, seaborn, and pyhmmer
echo "Installing scikit-learn, imbalanced-learn, pandas, seaborn, and pyhmmer..."
mamba install -c conda-forge -c bioconda \
    scikit-learn=0.22.1 imbalanced-learn pandas seaborn pyhmmer==0.10.14 --yes

# Install Prodigal, Screed, ruamel.yaml, Snakemake, Click, and conda-package-handling
echo "Installing Prodigal, Screed, ruamel.yaml, Snakemake, and related packages..."
mamba install -c conda-forge -c bioconda \
    prodigal screed ruamel.yaml "snakemake>=5.18,<=5.26" click "conda-package-handling<=1.9" --yes

# Install NumPy
echo "Installing NumPy..."
mamba install -c conda-forge -c bioconda numpy=1.23.5 --yes

# Install MEGAHIT, SPAdes, and fastp
echo "Installing MEGAHIT, SPAdes, and fastp..."
mamba install -c conda-forge -c bioconda megahit spades fastp --yes

# Install BioPython
echo "Installing BioPython..."
pip3 install bio

# Clean up the temporary cache directory
echo "Cleaning up temporary cache..."
rm -rf "$CONDA_PKGS_DIRS"

# Clone the VirSorter2-pyhmmerAcc repository and install the package
git clone https://github.com/liusihang/VirSorter2-pyhmmerAcc || { echo "Git clone failed"; exit 1; }
cd VirSorter2-pyhmmerAcc
pip install -e . || { echo "Pip install failed"; exit 1; }

# Run the Move2bin Python script
cd ..
python ./setupscript/Move2bin.py

# Get the current Conda environment path
CONDA_ENV_PATH=$(conda info --base)/envs/$(basename "$CONDA_DEFAULT_ENV")

# Add execute permissions to all files in the Conda environment's bin folder
chmod +x "$CONDA_ENV_PATH/bin/"*

echo "All packages installed successfully!"