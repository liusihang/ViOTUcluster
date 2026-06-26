#!/usr/bin/env bash

set -euo pipefail

#echo "Conda environment activated: $(conda info --envs)"

# Accept command-line arguments
INPUT_FASTA=$1
OUTPUT=$2

# Check if input file exists
if [ ! -f "$INPUT_FASTA" ]; then
    echo "Error: Input FASTA file '$INPUT_FASTA' not found."
    exit 1
fi

# Create output directories if they don't exist
mkdir -p "$OUTPUT/split_files" "$OUTPUT/iPhop_results"
VIOTUCLUSTER_PYTHON=${VIOTUCLUSTER_PYTHON:-python}

echo -e "\n\n\n# Performing Host Prediction Analysis!!!\n\n\n"
pwd

# Get the base filename of the input file (without path)
BASE_INPUT_FASTA=$(basename "$INPUT_FASTA")

# Split the input fasta file into smaller files, each containing 3000 sequences
awk -v output_dir="$OUTPUT/split_files" -v base_name="$BASE_INPUT_FASTA" 'BEGIN {n_seq=0;} 
     /^>/ {
        if (n_seq % 3000 == 0) {
            file = sprintf("%s/%s_%d.fna", output_dir, base_name, n_seq);
        } 
        print >> file; 
        n_seq++; 
        next;
     } 
     { print >> file; }' "$INPUT_FASTA"

# Change to the split files directory
cd "$OUTPUT/split_files" || exit

# List all split fna files
ls *.fna > iPhop

# Run the Python script for iPhop prediction
export OUTPUT
"$VIOTUCLUSTER_PYTHON" -m ViOTUcluster.run_iphop

echo "All tasks completed successfully."

# Merge all annotation results into the output directory
awk 'FNR==1 && NR!=1{next;} {print}' ./*_iPhopResult/Host_prediction_to_genus_m90.csv > "$OUTPUT/Combined_Host_prediction_to_genus_m90.tsv"

echo "Host prediction complete. Results combined and saved to $OUTPUT/Combined_Host_prediction_to_genus_m90.tsv"

# Clean up temporary files
echo "Cleaning up temporary files..."
rm -rf "$OUTPUT/split_files"
rm -rf "$OUTPUT/iPhop_results"/*_iPhopResult

echo "Cleanup complete."
#conda deactivate
