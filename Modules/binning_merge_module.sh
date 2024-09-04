#!/usr/bin/env bash

# 设置错误捕获机制：如果发生错误，脚本将停止执行
set -e
trap 'echo "An error occurred. Exiting..."; exit 1;' ERR

# 进行Binning分析
for FILE in $FILES; do
  echo "Processing $FILE"

  BASENAME=$(basename "$FILE" .fa)
  BASENAME=${BASENAME%.fasta}
  OUT_DIR="$OUTPUT_DIR/SeprateFile/${BASENAME}"

  # 匹配BASENAME_R1* 的格式
  Read1=$(find "${RAW_SEQ_DIR}" -type f -name "${BASENAME}_R1*" | head -n 1)
  Read2=$(find "${RAW_SEQ_DIR}" -type f -name "${BASENAME}_R2*" | head -n 1)

  if [ -z "$Read1" ] || [ -z "$Read2" ]; then
    echo "Error: Paired-end files for $BASENAME not found in the expected formats."
    exit 1
  fi

  echo "Using Read1: $Read1"
  echo "Using Read2: $Read2"

  mkdir -p "$OUT_DIR/Binning"

  if [ ! -f "$OUT_DIR/Binning/alignment.bam" ]; then
    # 创建BWA索引，并将索引文件输出到指定目录
    bwa index -p "$OUT_DIR/Binning/assembly_index" "$OUT_DIR/${BASENAME}_filtered.fasta"

    # 比对读数到装配的contigs上，并将结果输出到指定目录
    bwa mem -t 104 "$OUT_DIR/Binning/assembly_index" $Read1 $Read2 > "$OUT_DIR/Binning/alignment.sam"
    samtools view -S -b "$OUT_DIR/Binning/alignment.sam" > "$OUT_DIR/Binning/alignment.bam"
  else
    echo "Alignment already completed for $FILE, skipping..."
  fi

  # 设置conda环境和日志文件
  BASE_CONDA_PREFIX=$(conda info --base)
  conda_sh="$BASE_CONDA_PREFIX/etc/profile.d/conda.sh"
  CURRENT_ENV=$(basename "$CONDA_DEFAULT_ENV")
  Log_file="$OUT_DIR/Binning/vRhyme_results_${BASENAME}_filtered/log_vRhyme_${BASENAME}_filtered.log"

  if [ ! -f "$Log_file" ]; then
    # 激活vRhyme环境并运行vRhyme
    source ${conda_sh}
    conda activate vRhyme
    echo "Conda environment activated: $(conda info --envs)"
    which vRhyme

    vRhyme -i "$OUT_DIR/${BASENAME}_filtered.fasta" -b "$OUT_DIR/Binning/alignment.bam" -t 104 -o "$OUT_DIR/Binning/vRhyme_results_${BASENAME}_filtered"

    # 监控vRhyme任务是否完成
    all_tasks_completed=false
    while [ "$all_tasks_completed" == "false" ]; do
      sleep 30
      all_tasks_completed=true
      if ! grep -q "vRhyme binning complete" "$Log_file"; then
        all_tasks_completed=false
        echo "vRhyme command is still running in the background for $BASENAME"
      fi
      if [ "$all_tasks_completed" == "false" ]; then
        sleep 30
      fi
    done
  else
    echo "Viral binning already completed for $FILE, skipping..."
  fi

  conda activate "$CURRENT_ENV"

  # 判断是否需要重新组装
  if [ "$REASSEMBLE" = true ]; then
    ALL_BINS_FA="$OUT_DIR/Binning/summary_bins_contigs.fa"
    cat "$OUT_DIR/Binning/vRhyme_results_${BASENAME}_filtered/vRhyme_best_bins_fasta/"*.fasta > "$ALL_BINS_FA"

    bwa index -p "$OUT_DIR/Binning/all_bins_index" "$ALL_BINS_FA"

    mkdir -p "$OUT_DIR/Binning/reads_for_reassembly"
    bwa mem -t 104 "$OUT_DIR/Binning/all_bins_index" $Read1 $Read2 | python "${ScriptDir}/filter_reads_for_bin_reassembly.py" "$OUT_DIR/Binning/vRhyme_results_${BASENAME}_filtered/vRhyme_best_bins_fasta" "$OUT_DIR/Binning/reads_for_reassembly" $STRICT_MAX $PERMISSIVE_MAX

    for FASTQ_FILE in "$OUT_DIR/Binning/reads_for_reassembly/"*_1.fastq; do
      BIN_BASENAME=$(basename "$FASTQ_FILE" _1.fastq)
      OriginalBin=${BIN_BASENAME%%.*}
      EXTRACTED_DIR="$OUT_DIR/Binning/reassembled_bins"
      TMP_DIR="$EXTRACTED_DIR/${BIN_BASENAME}_tmp"
      mkdir -p "$TMP_DIR"

      spades.py -t 104 --tmp $TMP_DIR --careful --untrusted-contigs "$OUT_DIR/Binning/vRhyme_results_${BASENAME}_filtered/vRhyme_best_bins_fasta/${OriginalBin}.fasta" -1 "${OUT_DIR}/Binning/reads_for_reassembly/${BIN_BASENAME}_1.fastq" -2 "${OUT_DIR}/Binning/reads_for_reassembly/${BIN_BASENAME}_2.fastq" -o "$OUT_DIR/Binning/reassembile/${BIN_BASENAME}"

      mkdir -p "$EXTRACTED_DIR"
      cp "$OUT_DIR/Binning/reassembile/${BIN_BASENAME}/contigs.fasta" "$EXTRACTED_DIR/${BIN_BASENAME}.fasta"
      cp "$OUT_DIR/Binning/vRhyme_results_${BASENAME}_filtered/vRhyme_best_bins_fasta/${OriginalBin}.fasta" "$EXTRACTED_DIR/${OriginalBin}.origin.fasta"

      echo "Reassembly complete for ${BIN_BASENAME}"
      rm -rf "$TMP_DIR"
    done
        
    mkdir "$OUT_DIR/Binning/Summary"
    python "${ScriptDir}/concat_fasta_sequences.py" $EXTRACTED_DIR "$OUT_DIR/Binning/Summary/tempsummary.fasta"
    checkv end_to_end "$OUT_DIR/Binning/Summary/tempsummary.fasta" "$OUT_DIR/Binning/Summary/CheckRes" -t 100 -d "$DATABASE/checkv-db-v1.5"
    mkdir "$OUT_DIR/Binning/Summary/Finialbins"
    python "${ScriptDir}/select_best_bins.py" "$OUT_DIR/Binning/Summary/CheckRes/quality_summary.tsv" $EXTRACTED_DIR "$OUT_DIR/Binning/Summary/Finialfasta/Bestbins"

  else
    mkdir "$OUT_DIR/Binning/Summary"
    mkdir -p "$OUT_DIR/Binning/Summary/Finialfasta/Bestbins"
    cp "$OUT_DIR/Binning/vRhyme_results_${BASENAME}_filtered/vRhyme_best_bins_fasta/"*.fasta "$OUT_DIR/Binning/Summary/Finialfasta/Bestbins/"
  fi

  for vRhymeFILE in "$OUT_DIR/Binning/Summary/Finialfasta/Bestbins"/*.fasta; do
    NEW_NAME=$(basename "$vRhymeFILE" | sed 's/^vRhyme_//')
    mv "$vRhymeFILE" "$OUT_DIR/Binning/Summary/Finialfasta/Bestbins/$NEW_NAME"
  done

  python "${ScriptDir}/unbined.py" -i "$OUT_DIR/Binning/vRhyme_results_${BASENAME}_filtered/vRhyme_best_bins_fasta" -r "$OUT_DIR/${BASENAME}_filtered.fasta" -o "$OUT_DIR/Binning/Summary/Finialfasta/unbined.fasta"

  echo "Rebinning and reassembly complete for $FILE"
done