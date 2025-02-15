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

# Iterate over all possible _R1 files (including fastq.gz, fq.gz, fastq, fq formats)
for R1_FILE in "${INPUT_DIR}"/*_R1.{fq.gz,fastq.gz,fq,fastq}; do
    # Skip if no files match the pattern
    [ -e "$R1_FILE" ] || continue
    
    # Find the corresponding _R2 file
    R2_FILE="${R1_FILE/_R1/_R2}"

    # Ensure R2 file exists
    if [ -f "$R2_FILE" ]; then
        PREFIX=$(basename "$R1_FILE" | sed 's/_R1.*//')

        # Check if fastp output files already exist
        if [ -f "${Cleanreads}/${PREFIX}_R1.fq.gz" ] && [ -f "${Cleanreads}/${PREFIX}_R2.fq.gz" ]; then
            echo "Fastp already completed for $PREFIX, skipping..."
        else
            echo "Running fastp for $PREFIX..."
            fastp -i "$R1_FILE" -I "$R2_FILE" -o "${Cleanreads}/${PREFIX}_R1.fq.gz" -O "${Cleanreads}/${PREFIX}_R2.fq.gz"
            if [ $? -ne 0 ]; then
                echo "Error: fastp processing failed for $PREFIX."
                exit 1
            fi
        fi

        # Select assembly software and check if contigs already exist
        if [ "$ASSEMBLY_SOFTWARE" == "megahit" ]; then
            if [ -f "${CONTIGS_DIR}/${PREFIX}.fa" ]; then
                echo "Megahit assembly already completed for $PREFIX, skipping..."
            else
                echo "Running megahit for $PREFIX..."
                megahit -1 "${Cleanreads}/${PREFIX}_R1.fq.gz" -2 "${Cleanreads}/${PREFIX}_R2.fq.gz" -o "${ASSEMBLY_DIR}/${PREFIX}_megahit_out" -t "${THREADS}"
                if [ $? -ne 0 ]; then
                    echo "Error: Megahit assembly failed for $PREFIX."
                    exit 1
                fi

                # Extract assembly results
                cp "${ASSEMBLY_DIR}/${PREFIX}_megahit_out/final.contigs.fa" "${CONTIGS_DIR}/${PREFIX}.fa"
                if [ $? -ne 0 ]; then
                    echo "Error: Failed to copy megahit contigs for $PREFIX."
                    exit 1
                fi
            fi

        elif [ "$ASSEMBLY_SOFTWARE" == "metaspades" ]; then
            if [ -f "${CONTIGS_DIR}/${PREFIX}.fa" ]; then
                echo "Metaspades assembly already completed for $PREFIX, skipping..."
            else
                echo "Running metaspades for $PREFIX..."
                metaspades.py -1 "${Cleanreads}/${PREFIX}_R1.fq.gz" -2 "${Cleanreads}/${PREFIX}_R2.fq.gz" -o "${ASSEMBLY_DIR}/${PREFIX}_spades_out" -t "${THREADS}"
                if [ $? -ne 0 ]; then
                    echo "Error: Metaspades assembly failed for $PREFIX."
                    exit 1
                fi

                # Extract assembly results
                cp "${ASSEMBLY_DIR}/${PREFIX}_spades_out/scaffolds.fasta" "${CONTIGS_DIR}/${PREFIX}.fa"
                if [ $? -ne 0 ]; then
                    echo "Error: Failed to copy metaspades contigs for $PREFIX."
                    exit 1
                fi
            fi
        fi
    else
        echo "Paired file for $R1_FILE not found. Skipping..."
    fi
done

echo "Processing completed. Contigs are saved in $CONTIGS_DIR."