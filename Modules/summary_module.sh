#!/usr/bin/env bash
# Merge all
echo "Merging final sequences..."
mkdir "$OUTPUT_DIR/Summary/vOTU"
python "${ScriptDir}/Rename.py" -i "$OUTPUT_DIR/Summary/temp/DrepViralcontigs.fasta"
cat "$OUTPUT_DIR/Summary/temp/DrepViralcontigs.fasta" "$OUTPUT_DIR/Summary/temp/DrepBins.fasta" > "$OUTPUT_DIR/Summary/vOTU/vOTU.fasta"
checkv end_to_end "$OUTPUT_DIR/Summary/vOTU/vOTU.fasta" "$OUTPUT_DIR/Summary/vOTU/vOTU_CheckRes" -t "${THREADS}" -d "$DATABASE/checkv-db-v1.5"

for folder in "$OUTPUT_DIR/SeprateFile/"*/ ; do
  if [ -d "$folder" ]; then
    folderName=$(basename "$folder")
    echo "Processing folder: $folderName"

    combined_dest="$OUTPUT_DIR/Summary/SeperateRes/${folderName}_ViralList.fasta"
    touch "$combined_dest"

    unbined_source="${folder}Binning/Summary/Finialfasta/unbined.fasta"
    if [ -f "$unbined_source" ]; then
      cat "$unbined_source" >> "$combined_dest"
      echo -e "\n" >> "$combined_dest"
      rm "$unbined_source"
    else
      echo "Warning: $unbined_source does not exist."
    fi

    bestbins_source="${folder}Binning/Summary/Finialfasta/Bestbins/"
    if [ -d "$bestbins_source" ]; then
      for fasta_file in "$bestbins_source"*.fasta; do
        if [ -f "$fasta_file" ]; then
          cat "$fasta_file" >> "$combined_dest"
          echo -e "\n" >> "$combined_dest"
          rm "$fasta_file"
        else
          echo "Warning: No .fasta files found in $bestbins_source"
        fi
      done
      rmdir "$bestbins_source" 2>/dev/null || true
    else
      echo "Warning: $bestbins_source does not exist."
    fi

    finialfasta_dir="${folder}Binning/Summary/Finialfasta/"
    rmdir "$finialfasta_dir" 2>/dev/null || true

    summary_dir="${folder}Binning/Summary/"
    rmdir "$summary_dir" 2>/dev/null || true

    binning_dir="${folder}Binning/"
    rmdir "$binning_dir" 2>/dev/null || true
  fi
done

ABUNDANCE_CSV="$OUTPUT_DIR/Summary/vOTU/vOTU.Abundance.csv"

if [ -f "$ABUNDANCE_CSV" ]; then
    echo "vOTU.Abundance.csv already exists, skipping TPM calculation."
else
    set -e

    mkdir -p "$OUTPUT_DIR/Summary/Viralcontigs/TPMTemp"
    if [ $? -ne 0 ]; then
        echo "Error: Failed to create temporary directory for TPM calculation."
        exit 1
    fi

    bwa index -p "$OUTPUT_DIR/Summary/Viralcontigs/TPMTemp/TempIndex" "$OUTPUT_DIR/Summary/vOTU/vOTU.fasta"
    if [ $? -ne 0 ]; then
        echo "Error: Failed to build BWA index."
        exit 1
    fi

    for FILE in $FILES; do
        echo "Processing $FILE"
        
        BASENAME=$(basename "$FILE" .fa)
        BASENAME=${BASENAME%.fasta}
        
        if [ -f "$OUTPUT_DIR/Summary/Viralcontigs/Temp/${BASENAME}_coverage.tsv" ]; then
            echo "Skipping $BASENAME as coverage file already exists."
            continue
        fi
        
        Read1=$(find "${RAW_SEQ_DIR}" -type f -name "${BASENAME}_R1*" | head -n 1)
        Read2=$(find "${RAW_SEQ_DIR}" -type f -name "${BASENAME}_R2*" | head -n 1)
        
        if [ -z "$Read1" ] || [ -z "$Read2" ]; then
            echo "Error: Read1 or Read2 files not found for $BASENAME."
            exit 1
        fi

        bwa mem -t "${THREADS}" "$OUTPUT_DIR/Summary/Viralcontigs/TPMTemp/TempIndex" ${Read1} ${Read2} > "$OUTPUT_DIR/Summary/Viralcontigs/TPMTemp/${BASENAME}_gene.sam"
        if [ $? -ne 0 ]; then
            echo "Error: Failed to perform BWA alignment for $BASENAME."
            exit 1
        fi

        samtools view -bS --threads "${THREADS}" "$OUTPUT_DIR/Summary/Viralcontigs/TPMTemp/${BASENAME}_gene.sam" > "$OUTPUT_DIR/Summary/Viralcontigs/TPMTemp/${BASENAME}_gene.bam"
        if [ $? -ne 0 ]; then
            echo "Error: Failed to convert SAM to BAM for $BASENAME."
            exit 1
        fi

        samtools sort "$OUTPUT_DIR/Summary/Viralcontigs/TPMTemp/${BASENAME}_gene.bam" -o "$OUTPUT_DIR/Summary/Viralcontigs/TPMTemp/${BASENAME}_sorted_gene.bam" --threads "${THREADS}"
        if [ $? -ne 0 ]; then
            echo "Error: Failed to sort BAM file for $BASENAME."
            exit 1
        fi

        samtools index "$OUTPUT_DIR/Summary/Viralcontigs/TPMTemp/${BASENAME}_sorted_gene.bam"
        if [ $? -ne 0 ]; then
            echo "Error: Failed to generate BAM index for $BASENAME."
            exit 1
        fi

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

        rm -r "$OUTPUT_DIR/Summary/Viralcontigs/TPMTemp/binsf"
    done

    python ${ScriptDir}/TPM_caculate.py "$OUTPUT_DIR/Summary/Viralcontigs/Temp" "$ABUNDANCE_CSV"
    if [ $? -ne 0 ]; then
        echo "Error: Failed to run TPM calculation."
        exit 1
    fi

    rm -r "$OUTPUT_DIR/Summary/Viralcontigs/Temp"
    rm -r "$OUTPUT_DIR/Summary/Viralcontigs/TPMTemp"

    echo "TPM calculation completed successfully."
fi

TAXONOMY_CSV="$OUTPUT_DIR/Summary/vOTU/vOTU.Taxonomy.csv"

if [ -f "$TAXONOMY_CSV" ]; then
    echo "vOTU.Taxonomy.csv already exists, skipping Taxonomy prediction."
else
    echo "Starting taxonomy prediction..."
    genomad annotate "$OUTPUT_DIR/Summary/vOTU/vOTU.fasta" "$OUTPUT_DIR/Summary/vOTU/TaxAnnotate" $DATABASE/genomad_db -t "$THREADS"
    python ${ScriptDir}/format_taxonomy.py "$OUTPUT_DIR/Summary/vOTU/TaxAnnotate/vOTU_annotate/vOTU_taxonomy.tsv" "$TAXONOMY_CSV"
    echo "Taxonomy prediction completed successfully."
fi

echo "All files processed and combined successfully."