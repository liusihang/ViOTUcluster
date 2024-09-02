#!/usr/bin/env bash

# Create a temporary directory to store intermediate files for TPM calculation
mkdir -p "$OUTPUT_DIR/Summary/Viralcontigs/TPMTemp"

# Build BWA index for vOTU sequences
bwa index -p "$OUTPUT_DIR/Summary/Viralcontigs/TPMTemp/TempIndex" "$OUTPUT_DIR/Summary/Viralcontigs/vOTU.fasta"

# Perform Binning analysis
for FILE in $FILES; do
    echo "Processing $FILE"
    
    BASENAME=$(basename "$FILE" .fa)
    BASENAME=${BASENAME%.fasta}
    
    # Check and set the paths for Read1 and Read2
    PREFIX="${RAW_SEQ_DIR}/${BASENAME}"
    if [ -f "${PREFIX}_R1.fq" ] && [ -f "${PREFIX}_R2.fq" ]; then
        Read1="${PREFIX}_R1.fq"
        Read2="${PREFIX}_R2.fq"
    elif [ -f "${PREFIX}_R1.fastq" ] && [ -f "${PREFIX}_R2.fastq" ]; then
        Read1="${PREFIX}_R1.fastq"
        Read2="${PREFIX}_R2.fastq"
    elif [ -f "${PREFIX}_R1.fq.gz" ] && [ -f "${PREFIX}_R2.fq.gz" ]; then
        Read1="${PREFIX}_R1.fq.gz"
        Read2="${PREFIX}_R2.fq.gz"
    elif [ -f "${PREFIX}_R1.fastq.gz" ] && [ -f "${PREFIX}_R2.fastq.gz" ]; then
        Read1="${PREFIX}_R1.fastq.gz"
        Read2="${PREFIX}_R2.fastq.gz"
    else
        echo "Error: Paired-end files for $BASENAME not found in the expected formats (.fq, .fastq, .fq.gz, .fastq.gz)"
        exit 1
    fi
    
    # Perform alignment using BWA-MEM
    bwa mem -t 104 "$OUTPUT_DIR/Summary/Viralcontigs/TPMTemp/TempIndex" ${Read1} ${Read2} > "$OUTPUT_DIR/Summary/Viralcontigs/TPMTemp/${BASENAME}_gene.sam"
    
    # Convert SAM file to BAM format
    samtools view -bS --threads 104 "$OUTPUT_DIR/Summary/Viralcontigs/TPMTemp/${BASENAME}_gene.sam" > "$OUTPUT_DIR/Summary/Viralcontigs/TPMTemp/${BASENAME}_gene.bam"
    
    # Sort the BAM file by alignment position coordinates
    samtools sort "$OUTPUT_DIR/Summary/Viralcontigs/TPMTemp/${BASENAME}_gene.bam" -o "$OUTPUT_DIR/Summary/Viralcontigs/TPMTemp/${BASENAME}_sorted_gene.bam" --threads 104
    
    # Generate index for the sorted BAM file
    samtools index "$OUTPUT_DIR/Summary/Viralcontigs/TPMTemp/${BASENAME}_sorted_gene.bam"
    
    # Calculate coverage by creating necessary folders
    mkdir -p "$OUTPUT_DIR/Summary/Viralcontigs/TPMTemp/binsf"
    cp "$OUTPUT_DIR/Summary/Viralcontigs/vOTU.fasta" "$OUTPUT_DIR/Summary/Viralcontigs/TPMTemp/binsf"
    
    # Perform coverage calculation using CheckM
    checkm coverage -x fna -m 20 -t 104 "$OUTPUT_DIR/Summary/Viralcontigs/TPMTemp/binsf" "$OUTPUT_DIR/Summary/Viralcontigs/Temp/${BASENAME}_coverage.tsv" "$OUTPUT_DIR/Summary/Viralcontigs/TPMTemp/${BASENAME}_sorted_gene.bam"
    
    # Clean up intermediate files
    rm -r "$OUTPUT_DIR/Summary/Viralcontigs/TPMTemp/binsf"
    rm "$OUTPUT_DIR/Summary/Viralcontigs/TPMTemp/TempIndex"
done

# Run the TPM calculation Python script
python TPM_caculate.py "$OUTPUT_DIR/Summary/Viralcontigs/Temp" "$OUTPUT_DIR/Summary/Viralcontigs/Merged_vOTU_TPM.csv"