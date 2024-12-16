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

    # Perform Binning analysis
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
        echo "Aligning reads for $BASENAME..."
        bwa mem -t "${THREADS}" "$OUTPUT_DIR/Summary/Viralcontigs/TPMTemp/TempIndex" "${Read1}" "${Read2}" > "$OUTPUT_DIR/Summary/Viralcontigs/TPMTemp/${BASENAME}_gene.sam"
        if [ $? -ne 0 ]; then
            echo "Error: Failed to perform BWA alignment for $BASENAME."
            exit 1
        fi

        # Convert SAM to BAM
        echo "Converting SAM to BAM..."
        samtools view -bS --threads "${THREADS}" "$OUTPUT_DIR/Summary/Viralcontigs/TPMTemp/${BASENAME}_gene.sam" > "$OUTPUT_DIR/Summary/Viralcontigs/TPMTemp/${BASENAME}_gene.bam"
        if [ $? -ne 0 ]; then
            echo "Error: Failed to convert SAM to BAM for $BASENAME."
            exit 1
        fi

        # Sort BAM file by coordinates
        echo "Sorting BAM file..."
        samtools sort "$OUTPUT_DIR/Summary/Viralcontigs/TPMTemp/${BASENAME}_gene.bam" -o "$OUTPUT_DIR/Summary/Viralcontigs/TPMTemp/${BASENAME}_sorted_gene.bam" --threads "${THREADS}"
        if [ $? -ne 0 ]; then
            echo "Error: Failed to sort BAM file for $BASENAME."
            exit 1
        fi

        # Index the sorted BAM file
        echo "Indexing BAM file..."
        samtools index "$OUTPUT_DIR/Summary/Viralcontigs/TPMTemp/${BASENAME}_sorted_gene.bam"
        if [ $? -ne 0 ]; then
            echo "Error: Failed to generate BAM index for $BASENAME."
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
        checkm coverage -x fasta -m 20 -t "${THREADS}" "$OUTPUT_DIR/Summary/Viralcontigs/TPMTemp/binsf" "$OUTPUT_DIR/Summary/Viralcontigs/Temp/${BASENAME}_coverage.tsv" "$OUTPUT_DIR/Summary/Viralcontigs/TPMTemp/${BASENAME}_sorted_gene.bam"
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