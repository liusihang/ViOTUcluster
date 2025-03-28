#!/usr/bin/env bash

# Define the path to vOTU.Abundance.csv
ABUNDANCE_CSV="$OUTPUT_DIR/Summary/vOTU/vOTU.Abundance.csv"

# Check if vOTU.Abundance.csv already exists
if [ -f "$ABUNDANCE_CSV" ]; then
    echo "vOTU.Abundance.csv already exists, skipping..."
else
    # Create temporary directory for TPM calculation
    echo "Creating temporary directory for TPM calculation..."
    mkdir -p "$OUTPUT_DIR/Summary/Viralcontigs/TPMTemp"
    if [ $? -ne 0 ]; then
        echo "Error: Failed to create temporary directory for TPM calculation."
        exit 1
    fi

    # Build BWA index
    echo "Building BWA index..."
    bwa index -p "$OUTPUT_DIR/Summary/Viralcontigs/TPMTemp/TempIndex" "$OUTPUT_DIR/Summary/vOTU/vOTU.fasta"
    if [ $? -ne 0 ]; then
        echo "Error: Failed to build BWA index."
        exit 1
    fi

    for FILE in $FILES; do
        echo "Processing $FILE..."
        
        BASENAME=$(basename "$FILE" .fa)
        BASENAME=${BASENAME%.fasta}
        
        # Skip if coverage file already exists
        if [ -f "$OUTPUT_DIR/Summary/Viralcontigs/Temp/${BASENAME}_coverage.tsv" ]; then
            echo "Skipping $BASENAME as coverage file already exists."
            continue
        fi
        
        # Find Read1 and Read2 files
        Read1=$(find "${RAW_SEQ_DIR}" -type f -name "${BASENAME}_R1*" | head -n 1)
        Read2=$(find "${RAW_SEQ_DIR}" -type f -name "${BASENAME}_R2*" | head -n 1)
        
        if [ -z "$Read1" ] || [ -z "$Read2" ]; then
            echo "Error: Read1 or Read2 files not found for $BASENAME."
            exit 1
        fi

        # Align reads using BWA-MEM
        echo "Aligning reads, converting SAM to BAM, and sorting for ${BASENAME}..."
        bwa mem -t "${THREADS}" "$OUTPUT_DIR/Summary/Viralcontigs/TPMTemp/TempIndex" "${Read1}" "${Read2}" 2>> "${OUTPUT_DIR}/Log/Summary.log" | \
        sambamba view -S -f bam -t "${THREADS}" /dev/stdin 2>> "${OUTPUT_DIR}/Log/Summary.log" | \
        sambamba sort -t "${THREADS}" -o "$OUTPUT_DIR/Summary/Viralcontigs/TPMTemp/${BASENAME}_sorted_gene.bam" /dev/stdin 2>> "${OUTPUT_DIR}/Log/Summary.log"

        if [ $? -ne 0 ]; then
            echo "Error: Alignment, conversion, or sorting failed for ${BASENAME}."
            exit 1
        fi

        echo "Indexing sorted BAM file..."
        sambamba index -t "${THREADS}" "$OUTPUT_DIR/Summary/Viralcontigs/TPMTemp/${BASENAME}_sorted_gene.bam"
        if [ $? -ne 0 ]; then
            echo "Error: Failed to generate BAM index for ${BASENAME}."
            exit 1
        fi

        # Calculate coverage
        echo "Calculating coverage..."
        mkdir -p "$OUTPUT_DIR/Summary/Viralcontigs/TPMTemp/binsf"
        if [ $? -ne 0 ]; then
            echo "Error: Failed to create directory for coverage calculation."
            exit 1
        fi
        cp "$OUTPUT_DIR/Summary/vOTU/vOTU.fasta" "$OUTPUT_DIR/Summary/Viralcontigs/TPMTemp/binsf"
        checkm coverage -x fasta -m 0 -t "${THREADS}" "$OUTPUT_DIR/Summary/Viralcontigs/TPMTemp/binsf" "$OUTPUT_DIR/Summary/Viralcontigs/Temp/${BASENAME}_coverage.tsv" "$OUTPUT_DIR/Summary/Viralcontigs/TPMTemp/${BASENAME}_sorted_gene.bam"
        if [ $? -ne 0 ]; then
            echo "Error: Failed to calculate coverage for $BASENAME."
            exit 1
        fi

        # Clean up temporary files
        rm -r "$OUTPUT_DIR/Summary/Viralcontigs/TPMTemp/binsf"
    done

    # Run TPM calculation Python script
    echo "Running TPM calculation..."
    python "${ScriptDir}/TPM_caculate.py" "$OUTPUT_DIR/Summary/Viralcontigs/Temp" "$ABUNDANCE_CSV"
    if [ $? -ne 0 ]; then
        echo "Error: Failed to run TPM calculation."
        exit 1
    fi

    # Clean up temporary directories
    rm -r "$OUTPUT_DIR/Summary/Viralcontigs"
    echo "TPM calculation completed successfully."
fi