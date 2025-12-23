#!/usr/bin/env bash
# Merge all
echo "[🔄] Merging final sequences..."

# Define the path for the quality_summary.tsv file
QUALITY_SUMMARY="$OUTPUT_DIR/Summary/vOTU/vOTU_CheckRes/quality_summary.tsv"

# Define paths for input FASTA files
DREP_VIRAL_FASTA="$OUTPUT_DIR/Summary/dRepRes/DrepViralcontigs.fasta"
DREP_BINS_FASTA="$OUTPUT_DIR/Summary/dRepRes/DrepBins.fasta"

# Check if quality_summary.tsv already exists
if [ -f "$QUALITY_SUMMARY" ]; then
  echo "[✅] quality_summary.tsv already exists, skipping vOTU merging and CheckV analysis."
else
  # Create vOTU directory
  echo "[📁] Creating vOTU directory..."
  mkdir -p "$OUTPUT_DIR/Summary/vOTU"

  # Rename DrepViralcontigs.fasta file
  echo "[🔄] Renaming sequences..."
  python "${ScriptDir}/Rename.py" -i "$DREP_VIRAL_FASTA"

  # Merge DrepViralcontigs.fasta and DrepBins.fasta into vOTU.fasta
  echo "[🔄] Merging FASTA files..."
  cat "$DREP_VIRAL_FASTA" "$DREP_BINS_FASTA" > "$OUTPUT_DIR/Summary/vOTU/vOTU.fasta"

  # Run CheckV analysis
  echo "[🔄] Running CheckV analysis..."
  checkv end_to_end "$OUTPUT_DIR/Summary/vOTU/vOTU.fasta" "$OUTPUT_DIR/Summary/vOTU/vOTU_CheckRes" -t "${THREADS}" -d "$DATABASE/checkv-db-v1.5"
  if [ $? -ne 0 ]; then
    echo "[❌] Error: CheckV analysis failed."
    exit 1
  fi

  echo "[✅] CheckV analysis completed successfully."
fi

# Define the path for vOTU.Abundance.csv
ABUNDANCE_CSV="$OUTPUT_DIR/Summary/vOTU/vOTU.Abundance.csv"

# Check if vOTU.Abundance.csv already exists
if [ -f "$ABUNDANCE_CSV" ]; then
  echo "[✅] vOTU.Abundance.csv already exists, skipping TPM calculation."
else
  set -e

  # Create temporary directory for TPM calculation
  mkdir -p "$OUTPUT_DIR/Summary/Viralcontigs/TPMTemp"
  if [ $? -ne 0 ]; then
    echo "[❌] Error: Failed to create temporary directory for TPM calculation."
    exit 1
  fi

  # Check and generate sorted BAM file with index
  if [ ! -f "$OUTPUT_DIR/Summary/Viralcontigs/TPMTemp/TempIndex.sa" ]; then
      echo "[🔄] Building BWA index..."
      bwa index -b "100000000" -p "$OUTPUT_DIR/Summary/Viralcontigs/TPMTemp/TempIndex" "$OUTPUT_DIR/Summary/vOTU/vOTU.fasta"
      if [ $? -ne 0 ]; then
        echo "[❌] Error: Failed to build BWA index."
        exit 1
      fi
  else
      echo "[✅] BWA index already completed. Skipping..."
  fi

  # Function to process a single sample for BAM generation and coverage calculation
  process_sample_bam() {
    local FILE=$1
    echo "[🔄] Processing $FILE"
    
    local BASENAME=$(basename "$FILE" .fa)
    BASENAME=${BASENAME%.fasta}
    
    if [ -f "$OUTPUT_DIR/Summary/SeperateRes/Coverage/${BASENAME}_coverage.tsv" ]; then
      echo "[⏭️] Skipping $BASENAME as coverage file already exists."
      return 0
    fi
    
    local Read1=$(find "${RAW_SEQ_DIR}" -type f -name "${BASENAME}_R1*" | head -n 1)
    local Read2=$(find "${RAW_SEQ_DIR}" -type f -name "${BASENAME}_R2*" | head -n 1)
    
    if [ -z "$Read1" ] || [ -z "$Read2" ]; then
      echo "[❌] Error: Read1 or Read2 files not found for $BASENAME."
      return 1
    fi

    # Check and generate sorted BAM file with index
    if [ ! -f "$OUTPUT_DIR/Summary/Viralcontigs/TPMTemp/${BASENAME}_sorted_gene.bam" ]; then
        echo "[🔄] Running alignment, conversion, and sorting for ${BASENAME}..."
        bwa mem -t "${THREADS}" "$OUTPUT_DIR/Summary/Viralcontigs/TPMTemp/TempIndex" "${Read1}" "${Read2}" 2>> "${OUTPUT_DIR}/Log/Summary.log" | \
          sambamba view -S -f bam -t "${THREADS}" /dev/stdin | \
          sambamba sort -t "${THREADS}" -o "$OUTPUT_DIR/Summary/Viralcontigs/TPMTemp/${BASENAME}_sorted_gene.bam" /dev/stdin

        if [ $? -ne 0 ]; then
          echo "[❌] Error: Failed to complete alignment pipeline for $BASENAME."
          return 1
        fi

        sambamba index -t "${THREADS}" "$OUTPUT_DIR/Summary/Viralcontigs/TPMTemp/${BASENAME}_sorted_gene.bam"
        
        if [ $? -ne 0 ]; then
          echo "[❌] Error: Failed to generate BAM index for $BASENAME."
          return 1
        fi
    else
        echo "[✅] Alignment already completed for ${BASENAME}. Skipping..."
    fi

    # Create directory for coverage calculation
    mkdir -p "$OUTPUT_DIR/Summary/SeperateRes/Coverage"
    coverm contig --methods tpm --bam-files "$OUTPUT_DIR/Summary/Viralcontigs/TPMTemp/${BASENAME}_sorted_gene.bam" -t "${THREADS}" -o "$OUTPUT_DIR/Summary/SeperateRes/Coverage/${BASENAME}_coverage.tsv"
    if [ $? -ne 0 ]; then
      echo "[❌] Error: Failed to calculate coverage for $BASENAME."
      return 1
    fi
    
    echo "[✅] BAM processing completed for ${BASENAME}."
    return 0
  }
  
  export -f process_sample_bam
  export OUTPUT_DIR RAW_SEQ_DIR THREADS
  
  # Run BAM processing in parallel using TPM_tasks as concurrency limit
  TPM_tasks=${TPM_tasks:-4}
  echo "[🔄] Processing ${#FILES[@]} samples with parallel (max ${TPM_tasks} concurrent jobs)..."
  parallel -j ${TPM_tasks} process_sample_bam ::: $FILES

  # Run TPM calculation
  python ${ScriptDir}/TPM_caculate.py "$OUTPUT_DIR/Summary/SeperateRes/Coverage" "$ABUNDANCE_CSV"
  if [ $? -ne 0 ]; then
    echo "[❌] Error: Failed to run TPM calculation."
    exit 1
  fi

  #rm -r "$OUTPUT_DIR/Summary/Viralcontigs/Temp"
  #rm -r "$OUTPUT_DIR/Summary/Viralcontigs/TPMTemp"

  echo "[✅] TPM calculation completed successfully."
fi

# Define the path for vOTU.Taxonomy.csv
TAXONOMY_CSV="$OUTPUT_DIR/Summary/vOTU/vOTU.Taxonomy.csv"

# Check if vOTU.Taxonomy.csv already exists
if [ -f "$TAXONOMY_CSV" ]; then
  echo "[✅] vOTU.Taxonomy.csv already exists, skipping Taxonomy prediction."
else
  echo "[🔬] Starting taxonomy prediction..."
  genomad annotate "$OUTPUT_DIR/Summary/vOTU/vOTU.fasta" "$OUTPUT_DIR/Summary/vOTU/TaxAnnotate" $DATABASE/genomad_db -t "$THREADS"
  python ${ScriptDir}/format_taxonomy.py "$OUTPUT_DIR/Summary/vOTU/TaxAnnotate/vOTU_annotate/vOTU_taxonomy.tsv" "$TAXONOMY_CSV" "$ABUNDANCE_CSV"
  echo "[✅] Taxonomy prediction completed successfully."
fi

rm -r "$OUTPUT_DIR/Summary/temp"
rm -r "$OUTPUT_DIR/Summary/dRepRes"
rm -r "$OUTPUT_DIR/Summary/Viralcontigs"

echo "[✅] All files processed and combined successfully."
