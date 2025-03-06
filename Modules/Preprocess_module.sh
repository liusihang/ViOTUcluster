#!/usr/bin/env bash

# Check if the correct number of arguments is provided
if [ $# -lt 3 ]; then
    echo "Usage: $0 <input_directory> <assembly_software> <output_directory>"
    echo "assembly_software options: megahit, metaspades"
    exit 1
fi

# Get input arguments
INPUT_DIR=$1
ASSEMBLY_SOFTWARE=$2
OUTPUT_DIR=$3
THREADS=$4

# Validate assembly software option
if [[ "$ASSEMBLY_SOFTWARE" != "megahit" && "$ASSEMBLY_SOFTWARE" != "metaspades" ]]; then
    echo "Error: Invalid assembly software. Please choose either 'megahit' or 'metaspades'."
    exit 1
fi

# Create output directories
Cleanreads="${OUTPUT_DIR}/Cleanreads"
ASSEMBLY_DIR="${OUTPUT_DIR}/Assembly"
CONTIGS_DIR="${OUTPUT_DIR}/Contigs"

mkdir -p "$Cleanreads" "$ASSEMBLY_DIR" "$CONTIGS_DIR"
if [ $? -ne 0 ]; then
    echo "Error: Failed to create output directories."
    exit 1
fi

python "${ScriptDir}/ContigsPreprocess.py" -i ${INPUT_DIR} -o ${OUTPUT_DIR} -c ${THREADS} -a ${ASSEMBLY_SOFTWARE}

echo "Processing completed. Contigs are saved in $CONTIGS_DIR."