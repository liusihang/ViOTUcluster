#!/usr/bin/env bash
# Merge all
echo "Merging final sequences..."
mkdir "$OUTPUT_DIR/Summary/vOTU"
python "${ScriptDir}/Rename.py" -i "$OUTPUT_DIR/Summary/temp/DrepViralcontigs.fasta"
cat "$OUTPUT_DIR/Summary/temp/DrepViralcontigs.fasta" "$OUTPUT_DIR/Summary/temp/DrepBins.fasta" > "$OUTPUT_DIR/Summary/vOTU/vOTU.fasta"
checkv end_to_end "$OUTPUT_DIR/Summary/vOTU/vOTU.fasta" "$OUTPUT_DIR/Summary/vOTU/vOTU_CheckRes" -t "${THREADS}" -d "$DATABASE/checkv-db-v1.5"

# 创建目标目录用于存放汇总后的 ViralList.fasta 文件
mkdir -p "$OUTPUT_DIR/Summary/SeparateRes"

# 遍历当前目录下的所有子目录
for folder in "$OUTPUT_DIR/Summary/SeprateFile/"*/ ; do
  if [ -d "$folder" ]; then
    # 获取文件夹的名字，去除尾部的斜线
    folderName=$(basename "$folder")
    echo "Processing folder: $folderName"

    # 初始化目标文件路径，用于汇总
    combined_dest="$OUTPUT_DIR/Summary/SeparateRes/${folderName}_ViralList.fasta"
    touch "$combined_dest"

    # 处理unbined.fasta文件
    unbined_source="${folder}Binning/Summary/Finialfasta/unbined.fasta"
    if [ -f "$unbined_source" ]; then
      cat "$unbined_source" >> "$combined_dest"
      echo -e "\n" >> "$combined_dest" # 添加一个换行符以确保格式正确
      # 删除源文件
      rm "$unbined_source"
    else
      echo "Warning: $unbined_source does not exist."
    fi

    # 处理Bestbins文件夹中的fasta文件
    bestbins_source="${folder}Binning/Summary/Finialfasta/Bestbins/"
    if [ -d "$bestbins_source" ]; then
      for fasta_file in "$bestbins_source"*.fasta; do
        if [ -f "$fasta_file" ]; then
          cat "$fasta_file" >> "$combined_dest"
          echo -e "\n" >> "$combined_dest" # 添加一个换行符以确保格式正确
          # 删除源文件
          rm "$fasta_file"
        else
          echo "Warning: No .fasta files found in $bestbins_source"
        fi
      done
      # 删除Bestbins目录，如果已空
      rmdir "$bestbins_source" 2>/dev/null || true
    else
      echo "Warning: $bestbins_source does not exist."
    fi

    # 删除Finialfasta目录，如果已空
    finialfasta_dir="${folder}Binning/Summary/Finialfasta/"
    rmdir "$finialfasta_dir" 2>/dev/null || true

    # 删除Summary目录，如果已空
    summary_dir="${folder}Binning/Summary/"
    rmdir "$summary_dir" 2>/dev/null || true

    # 删除Binning目录，如果已空
    binning_dir="${folder}Binning/"
    rmdir "$binning_dir" 2>/dev/null || true
  fi
done

set -e  # Exit immediately if any command exits with a non-zero status

# Create a temporary directory to store intermediate files for TPM calculation
mkdir -p "$OUTPUT_DIR/Summary/Viralcontigs/TPMTemp"
if [ $? -ne 0 ]; then
    echo "Error: Failed to create temporary directory for TPM calculation."
    exit 1
fi

# Build BWA index for vOTU sequences
bwa index -p "$OUTPUT_DIR/Summary/Viralcontigs/TPMTemp/TempIndex" "$OUTPUT_DIR/Summary/vOTU/vOTU.fasta"
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
    cp "$OUTPUT_DIR/Summary/vOTU/vOTU.fasta" "$OUTPUT_DIR/Summary/Viralcontigs/TPMTemp/binsf"
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


#Taxonmy prediction
genomad annotate "$OUTPUT_DIR/Summary/vOTU/vOTU.fasta" "$OUTPUT_DIR/Summary/vOTU/TaxAnnotate" $DATABASE/genomad_db -t "$THREADS"
python ${ScriptDir}/format_taxonomy.py "$OUTPUT_DIR/Summary/vOTU/TaxAnnotate/vOTU_annotate/vOTU_taxonomy.tsv" "$OUTPUT_DIR/Summary/vOTU/vOTU.Taxonomy.csv"
#rm -r "$OUTPUT_DIR/Summary/vOTU/TaxAnnotate"

echo "All files processed and combined successfully."