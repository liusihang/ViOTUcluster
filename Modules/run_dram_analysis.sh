#!/usr/bin/env bash

#source activate DRAM

# Check command-line arguments
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input_fasta_file> <output_directory>"
    exit 1
fi

# Accept command-line arguments
INPUT_FASTA=$1
OUTPUT_DIR=$2

# Check if output file already exists
if [ -f "$OUTPUT_DIR/DRAM_annotations.tsv" ]; then
    echo "Output file '$OUTPUT_DIR/DRAM_annotations.tsv' already exists. Skipping analysis."
    exit 0
fi

# Check if input file exists
if [ ! -f "$INPUT_FASTA" ]; then
    echo "Error: Input FASTA file '$INPUT_FASTA' not found."
    exit 1
fi

# Create output directories if they don't exist
mkdir -p "$OUTPUT_DIR"
mkdir -p "$OUTPUT_DIR/split_files"

echo -e "\n\n\n# Performing DRAM analysis!!!\n\n\n"
pwd

# Get the base filename of the input file (without path)
BASE_INPUT_FASTA=$(basename "$INPUT_FASTA")

# Split the input fasta file into smaller files, each containing 1000 sequences
awk -v output_dir="$OUTPUT_DIR/split_files" -v base_name="$BASE_INPUT_FASTA" 'BEGIN {n_seq=0;} 
     /^>/ {
        if (n_seq % 1000 == 0) {
            file = sprintf("%s/%s_%d.fna", output_dir, base_name, n_seq);
        } 
        print >> file; 
        n_seq++; 
        next;
     } 
     { print >> file; }' "$INPUT_FASTA"

# Change to the split files directory
cd "$OUTPUT_DIR/split_files" || exit

# List all split fna files
ls *.fna > DRAM

# Run the Python script for DRAM annotation
python "${ScriptDir}/run_DRAM.py"

all_tasks_completed=false

# Monitor task completion
while [ "$all_tasks_completed" == "false" ]; do
    sleep 30
    all_tasks_completed=true

    # Iterate over all directories ending with _DRAMAnnot
    for dir in *_DRAMAnnot; do
        if [ ! -f "$dir/annotations.tsv" ]; then
            echo "DRAM annotation still in progress in $dir."
            all_tasks_completed=false
            break
        fi
    done

    # If not completed, wait another 30 seconds
    if [ "$all_tasks_completed" == "false" ]; then
        sleep 30
    fi
done

echo "All DRAM annotations completed."

# Merge all annotation results into the output directory
awk 'FNR==1 && NR!=1{next;} {print}' "$OUTPUT_DIR"/split_files/*_DRAMAnnot/annotations.tsv > "$OUTPUT_DIR/DRAM_annotations.tsv"

echo "Annotation complete. Results combined and saved to $OUTPUT_DIR/DRAM_annotations.tsv"

# Clean up temporary files
echo "Cleaning up temporary files..."
rm -rf "$OUTPUT_DIR/split_files"
rm -rf "$OUTPUT_DIR/DRAM_results"/*_DRAMAnnot

echo "Cleanup complete."
#conda deactivate