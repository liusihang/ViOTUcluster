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
mamba env create -f ./setupscript/vOTUcluster.yaml
mamba activate vOTUcluster

# Clone the repository and install the package
git clone https://gitclone.com/github.com/liusihang/VirSorter2-pyhmmerAcc.git
cd VirSorter2-pyhmmerAcc
pip install -e .

# Run Python scripts
python ./setupscript/Move2bin.py
python ./setupscript/ReplaceVf.py

# Setup the database for VirSorter
virsorter setup -d "$DB_DIR" -j "$NUM_THREADS"
checkv download_database "$DB_DIR"
genomad download-database "$DB_DIR"

