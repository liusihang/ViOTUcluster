#!/usr/bin/env bash

# 设置错误捕获机制：如果发生错误，脚本将停止执行
set -e
trap 'echo "An error occurred. Exiting..."; exit 1;' ERR

# 进行 Binning 分析
for FILE in $FILES; do
  echo "Processing $FILE"

  BASENAME=$(basename "$FILE" .fa)
  BASENAME=${BASENAME%.fasta}
  OUT_DIR="$OUTPUT_DIR/SeprateFile/${BASENAME}"

  # 匹配 BASENAME_R1* 的格式
  Read1=$(find "${RAW_SEQ_DIR}" -maxdepth 1 -type f -name "${BASENAME}_R1*" | head -n 1)
  Read2=$(find "${RAW_SEQ_DIR}" -maxdepth 1 -type f -name "${BASENAME}_R2*" | head -n 1)

  if [ -z "$Read1" ] || [ -z "$Read2" ]; then
    echo "Error: Paired-end files for $BASENAME not found in the expected formats."
    exit 1
  fi

  echo "Using Read1: $Read1"
  echo "Using Read2: $Read2"

  # 创建 Binning 文件夹
  if [ ! -d "$OUT_DIR/Binning" ]; then
    mkdir -p "$OUT_DIR/Binning"
  fi

  # 检查并生成 BAM 文件
  if [ ! -f "$OUT_DIR/Binning/alignment.bam" ]; then
    bwa index -p "$OUT_DIR/Binning/assembly_index" "$OUT_DIR/${BASENAME}_filtered.fasta" >> "${OUTPUT_DIR}/Log/binning_merge.log" 2>&1
    bwa mem -t "${THREADS}" "$OUT_DIR/Binning/assembly_index" $Read1 $Read2 > "$OUT_DIR/Binning/alignment.sam" 2>> "${OUTPUT_DIR}/Log/binning_merge.log"
    samtools view -S -b "$OUT_DIR/Binning/alignment.sam" > "$OUT_DIR/Binning/alignment.bam" 2>> "${OUTPUT_DIR}/Log/binning_merge.log"
  else
    echo "$OUT_DIR/Binning/alignment.bam already exists. Skipping alignment." >> "${OUTPUT_DIR}/Log/binning_merge.log"
  fi

  VRHYME_DIR="$OUT_DIR/Binning/vRhyme_results_${BASENAME}_filtered"
  LOG_FILE="$VRHYME_DIR/log_vRhyme_${BASENAME}_filtered.log"

  # 检查 vRhyme 结果是否已经存在且已完成
  if [ -f "$LOG_FILE" ] && grep -q "Writing finalized bin sequences to individual fasta files" "$LOG_FILE"; then
    echo "vRhyme results for $BASENAME already exist. Skipping vRhyme run."
  else
    echo "Running vRhyme for $BASENAME..."
    source $(conda info --base)/etc/profile.d/conda.sh
    conda activate vRhyme

    # 删除已有的 vRhyme 结果文件夹
    if [ -d "$VRHYME_DIR" ]; then
      echo "Deleting existing vRhyme results directory: $VRHYME_DIR"
      rm -rf "$VRHYME_DIR"
    fi

    vRhyme -i "$OUT_DIR/${BASENAME}_filtered.fasta" \
          -b "$OUT_DIR/Binning/alignment.bam" \
          -t "${THREADS_PER_FILE}" \
          -o "$VRHYME_DIR"

    conda deactivate
  fi

  # 判断是否需要重新组装
  if [ "$REASSEMBLE" = true ]; then
    ALL_BINS_FA="$OUT_DIR/Binning/summary_bins_contigs.fa"
    cat "$VRHYME_DIR/vRhyme_best_bins_fasta/"*.fasta > "$ALL_BINS_FA"

    if [ ! -d "$OUT_DIR/Binning/reads_for_reassembly" ]; then
      mkdir -p "$OUT_DIR/Binning/reads_for_reassembly"
    fi

    bwa index -p "$OUT_DIR/Binning/all_bins_index" "$ALL_BINS_FA"
    bwa mem -t "${THREADS}" "$OUT_DIR/Binning/all_bins_index" $Read1 $Read2 | python "${ScriptDir}/filter_reads_for_bin_reassembly.py" "$VRHYME_DIR/vRhyme_best_bins_fasta" "$OUT_DIR/Binning/reads_for_reassembly" $STRICT_MAX $PERMISSIVE_MAX

    for FASTQ_FILE in "$OUT_DIR/Binning/reads_for_reassembly/"*_1.fastq; do
      BIN_BASENAME=$(basename "$FASTQ_FILE" _1.fastq)
      OriginalBin=${BIN_BASENAME%%.*}
      EXTRACTED_DIR="$OUT_DIR/Binning/reassembled_bins"

      if [ ! -d "$EXTRACTED_DIR/${BIN_BASENAME}_tmp" ]; then
        mkdir -p "$EXTRACTED_DIR/${BIN_BASENAME}_tmp"
      fi

      spades.py -t "${THREADS}" --tmp "$EXTRACTED_DIR/${BIN_BASENAME}_tmp" --careful --untrusted-contigs "$VRHYME_DIR/vRhyme_best_bins_fasta/${OriginalBin}.fasta" -1 "$OUT_DIR/Binning/reads_for_reassembly/${BIN_BASENAME}_1.fastq" -2 "$OUT_DIR/Binning/reads_for_reassembly/${BIN_BASENAME}_2.fastq" -o "$OUT_DIR/Binning/reassembile/${BIN_BASENAME}"

      if [ ! -d "$EXTRACTED_DIR" ]; then
        mkdir -p "$EXTRACTED_DIR"
      fi
      cp "$OUT_DIR/Binning/reassembile/${BIN_BASENAME}/contigs.fasta" "$EXTRACTED_DIR/${BIN_BASENAME}.fasta"
    done

    if [ ! -d "$OUT_DIR/Binning/Summary" ]; then
      mkdir "$OUT_DIR/Binning/Summary"
    fi
    python "${ScriptDir}/concat_fasta_sequences.py" "$EXTRACTED_DIR" "$OUT_DIR/Binning/Summary/tempsummary.fasta"
  fi

  # 创建 bins 文件夹
  if [ ! -d "${OUTPUT_DIR}/Summary/SeperateRes/bins" ]; then
    mkdir -p "${OUTPUT_DIR}/Summary/SeperateRes/bins"
  fi

  # 定义 unbined 输出文件路径
  UNBINNED_FASTA="$OUTPUT_DIR/Summary/SeperateRes/unbined/${BASENAME}_unbined.fasta"

  # 如果 unbined fasta 文件已存在，则跳过以下所有步骤
  if [ -f "$UNBINNED_FASTA" ]; then
      echo "All processing steps already completed for $FILE, skipping..."
  else
      # 重命名并复制 vRhyme 结果
      for vRhymeFILE in "$VRHYME_DIR/vRhyme_best_bins_fasta/"*.fasta; do
          NEW_NAME=$(basename "$vRhymeFILE" | sed "s/^vRhyme_/${BASENAME}_/")
          NEW_PATH="$OUT_DIR/Binning/Summary/Finialfasta/Bestbins/$NEW_NAME"
          
          # 检查目标文件夹是否存在，如果不存在则创建
          if [ ! -d "${OUT_DIR}/Binning/Summary/Finialfasta/Bestbins" ]; then
              mkdir -p "${OUT_DIR}/Binning/Summary/Finialfasta/Bestbins"
          fi
          
          mv "$vRhymeFILE" "$NEW_PATH"
          cp "$NEW_PATH" "${OUTPUT_DIR}/Summary/SeperateRes/bins"
      done

      # 获取合并的 bins 和 unbined 序列
      python "${ScriptDir}/MergeBins.py" -i "$VRHYME_DIR/vRhyme_best_bins_fasta" -o "${OUTPUT_DIR}/Summary/SeperateRes/bins/${BASENAME}_bins.fasta"

      # 生成 unbined 序列
      if [ ! -d "$OUTPUT_DIR/Summary/SeperateRes/unbined" ]; then
          mkdir -p "$OUTPUT_DIR/Summary/SeperateRes/unbined"
      fi
      python "${ScriptDir}/unbined.py" -i "$VRHYME_DIR/vRhyme_best_bins_fasta" -r "$OUT_DIR/${BASENAME}_filtered.fasta" -o "$UNBINNED_FASTA"

      # 合并 bins 和 unbined 序列
      cat "${OUTPUT_DIR}/Summary/SeperateRes/bins/${BASENAME}_bins.fasta" "${OUTPUT_DIR}/Summary/SeperateRes/unbined/${BASENAME}_unbined.fasta" > "${OUTPUT_DIR}/Summary/SeperateRes/${BASENAME}_viralseqs.fasta"
      
      echo "Rebinning and reassembly complete for $FILE"
  fi
done