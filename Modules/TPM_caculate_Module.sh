#!/usr/bin/env bash

set -e  # Exit immediately if any command exits with a non-zero status

# Create a temporary directory to store intermediate files for TPM calculation
mkdir -p "$OUTPUT_DIR/Summary/Viralcontigs/TPMTemp"
if [ $? -ne 0 ]; then
    echo "Error: Failed to create temporary directory for TPM calculation."
    exit 1
fi

# Build BWA index for vOTU sequences
bwa index -p "$OUTPUT_DIR/Summary/Viralcontigs/TPMTemp/TempIndex" "$OUTPUT_DIR/Summary/Viralcontigs/vOTU.fasta"
if [ $? -ne 0 ]; then
    echo "Error: Failed to build BWA index."
    exit 1
fi

# Perform Binning analysis
for FILE in $FILES; do
    echo "Processing $FILE"
    
    BASENAME=$(basename "$FILE" .fa)
    BASENAME=${BASENAME%.fasta}
    
    # 检查是否已经存在coverage文件，若存在则跳过
    if [ -f "$OUTPUT_DIR/Summary/Viralcontigs/Temp/${BASENAME}_coverage.tsv" ]; then
        echo "Skipping $BASENAME as coverage file already exists."
        continue
    fi
    
    # 检查并设置Read1和Read2的路径
    Read1=$(find "${RAW_SEQ_DIR}" -type f -name "${BASENAME}_R1*" | head -n 1)
    Read2=$(find "${RAW_SEQ_DIR}" -type f -name "${BASENAME}_R2*" | head -n 1)
    
    if [ -z "$Read1" ] || [ -z "$Read2" ]; then
        echo "Error: Read1 or Read2 files not found for $BASENAME."
        exit 1
    fi

    # Perform alignment using BWA-MEM
    bwa mem -t "${THREADS}" "$OUTPUT_DIR/Summary/Viralcontigs/TPMTemp/TempIndex" ${Read1} ${Read2} > "$OUTPUT_DIR/Summary/Viralcontigs/TPMTemp/${BASENAME}_gene.sam"
    if [ $? -ne 0 ]; then
        echo "Error: Failed to perform BWA alignment for $BASENAME."
        exit 1
    fi

    # Convert SAM file to BAM format
    samtools view -bS --threads "${THREADS}" "$OUTPUT_DIR/Summary/Viralcontigs/TPMTemp/${BASENAME}_gene.sam" > "$OUTPUT_DIR/Summary/Viralcontigs/TPMTemp/${BASENAME}_gene.bam"
    if [ $? -ne 0 ]; then
        echo "Error: Failed to convert SAM to BAM for $BASENAME."
        exit 1
    fi

    # Sort the BAM file by alignment position coordinates
    samtools sort "$OUTPUT_DIR/Summary/Viralcontigs/TPMTemp/${BASENAME}_gene.bam" -o "$OUTPUT_DIR/Summary/Viralcontigs/TPMTemp/${BASENAME}_sorted_gene.bam" --threads "${THREADS}"
    if [ $? -ne 0 ]; then
        echo "Error: Failed to sort BAM file for $BASENAME."
        exit 1
    fi

    # Generate index for the sorted BAM file
    samtools index "$OUTPUT_DIR/Summary/Viralcontigs/TPMTemp/${BASENAME}_sorted_gene.bam"
    if [ $? -ne 0 ]; then
        echo "Error: Failed to generate BAM index for $BASENAME."
        exit 1
    fi

    # Calculate coverage
    mkdir -p "$OUTPUT_DIR/Summary/Viralcontigs/TPMTemp/binsf"
    if [ $? -ne 0 ]; then
        echo "Error: Failed to create directory for coverage calculation."
        exit 1
    fi
    cp "$OUTPUT_DIR/Summary/Viralcontigs/vOTU.fasta" "$OUTPUT_DIR/Summary/Viralcontigs/TPMTemp/binsf"
    checkm coverage -x fasta -m 20 -t "${THREADS}" "$OUTPUT_DIR/Summary/Viralcontigs/TPMTemp/binsf" "$OUTPUT_DIR/Summary/Viralcontigs/Temp/${BASENAME}_coverage.tsv" "$OUTPUT_DIR/Summary/Viralcontigs/TPMTemp/${BASENAME}_sorted_gene.bam"
    if [ $? -ne 0 ]; then
        echo "Error: Failed to calculate coverage for $BASENAME."
        exit 1
    fi

    # Clean up intermediate files
    rm -r "$OUTPUT_DIR/Summary/Viralcontigs/TPMTemp/binsf"
    #rm -r "$OUTPUT_DIR/Summary/Viralcontigs/TPMTemp/TempIndex"
done

# Run the TPM calculation Python script
python ${ScriptDir}/TPM_caculate.py "$OUTPUT_DIR/Summary/Viralcontigs/Temp" "$OUTPUT_DIR/Summary/vOTU/vOTU.Abundance.csv"
if [ $? -ne 0 ]; then
    echo "Error: Failed to run TPM calculation."
    exit 1
fi
# Clean up intermediate files
rm -r "$OUTPUT_DIR/Summary/Viralcontigs/Temp"
rm -r "$OUTPUT_DIR/Summary/Viralcontigs/TPMTemp"

echo "TPM calculation completed successfully."