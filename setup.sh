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
#mamba env create -f vOTUcluster
#mamba activate vOTUcluster

#clean cache
conda clean --all

#instasll packages
#mamba install virsorter=2 --clobber
#mamba create -n vOTUcluster8.27 -c conda-forge -c bioconda "python=3.8" 
mamba install -c conda-forge -c bioconda dRep=3.5.0 viralverify=1.1 genomad=1.8.0 checkv=1.0.3 scikit-learn=0.22.1 imbalanced-learn pandas seaborn pyhmmer==0.10.14 prodigal screed ruamel.yaml "snakemake>=5.18,<=5.26" click "conda-package-handling<=1.9"
mamba install -c conda-forge -c bioconda numpy=1.23.5 checkm-genome=1.2.2
mamba install -c conda-forge -c bioconda megahit spades fastp

# Clone the repository and install the package
git clone https://gitclone.com/github.com/liusihang/VirSorter2-pyhmmerAcc
cd VirSorter2-pyhmmerAcc
pip install -e .

# Run Python scripts
cd ..
python ./setupscript/Move2bin.py

# Get the current conda environment path
CONDA_ENV_PATH=$(conda info --base)/envs/$(basename "$CONDA_DEFAULT_ENV")

# Add execute permissions to all files in the current conda environment's bin folder
chmod +x "$CONDA_ENV_PATH/bin/*"

# Setup the database for VirSorter
virsorter setup -d "$DB_DIR" -j "$NUM_THREADS"
checkv download_database "$DB_DIR"
genomad download-database "$DB_DIR"

# Download and unzip ViralVerify database
echo "Downloading ViralVerify database..."
curl -L -o "$DB_DIR/viralverify.zip" "https://figshare.com/ndownloader/files/17904323?private_link=f897d463b31a35ad7bf0"

echo "Unzipping ViralVerify database..."
unzip "$DB_DIR/viralverify.zip" -d "$DB_DIR/viralverify"