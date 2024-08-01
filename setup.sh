#!/bin/bash

#Usage: ./setup_vOTUcluster.sh /path/to/db 4
# Check if the correct number of arguments is provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <db_directory> <num_threads>"
    exit 1
fi

DB_DIR=$1
NUM_THREADS=$2

# Create and activate the Conda environment
# mamba create -n vOTUcluster -c conda-forge -c bioconda "python=3.8" scikit-learn=0.22.1 imbalanced-learn pandas seaborn pyhmmer==0.10.14 prodigal screed ruamel.yaml "snakemake>=5.18,<=5.26" click "conda-package-handling<=1.9"
# mamba activate vOTUcluster

# Clone the repository and install the package
git clone https://github.com/liusihang/VirSorter2-pyhmmerAcc.git
cd VirSorter2-pyhmmerAcc
pip install -e .

# Install additional packages
mamba install -c conda-forge -c bioconda viralverify
mamba install -c conda-forge -c bioconda geNomad
mamba install -c conda-forge -c bioconda checkv
mamba install -c conda-forge -c bioconda drep

# Run Python scripts
python Move2bin.py
python ReplaceVf.py

# Setup the database for VirSorter
virsorter setup -d "$DB_DIR" -j "$NUM_THREADS"
checkv download_database "$DB_DIR"
genomad download-database "$DB_DIR"

