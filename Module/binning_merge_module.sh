#!/usr/bin/env bash

# 进行Binning分析
for FILE in $FILES; do
  echo "Processing $FILE"

  # MIN_CONTIG_LEN=300  # 最小contig长度
  STRICT_MAX=2  # 严格比对的最大SNP数
  PERMISSIVE_MAX=5  # 宽松比对的最大SNP数

  BASENAME=$(basename "$FILE" .fa)
  BASENAME=${BASENAME%.fasta}
  OUT_DIR="$OUTPUT_DIR/SeprateFile/${BASENAME}"
  # 检查并设置Read1和Read2的路径
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

  Log_file="$OUT_DIR/Binning/vRhyme_results_${BASENAME}_filtered/log_vRhyme_${BASENAME}_filtered.log"
  if [ ! -f "$Log_file" ]; then
    # vRhyme binning
    nohup bash -c "
      source ~/miniconda3/etc/profile.d/conda.sh;
      conda activate vRhyme;
      echo 'Conda environment activated:';
      conda info --envs;
      which vRhyme;
      vRhyme -i \"$OUT_DIR/${BASENAME}_filtered.fasta\" -b \"$OUT_DIR/Binning/alignment.bam\" -t 104 -o \"$OUT_DIR/Binning/vRhyme_results_${BASENAME}_filtered\";
    " > "$Log_file" 2>&1 &

    # 检测任务是否完成
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

  if [ "$REASSEMBLE" = true ]; then
    # 创建大的FA文件用于比对
    ALL_BINS_FA="$OUT_DIR/Binning/summary_bins_contigs.fa"
    cat "$OUT_DIR/Binning/vRhyme_results_${BASENAME}_filtered/vRhyme_best_bins_fasta/"*.fasta > "$ALL_BINS_FA"

    # 创建新的BWA索引，并将索引文件输出到指定目录
    bwa index -p "$OUT_DIR/Binning/all_bins_index" "$ALL_BINS_FA"

    # Aligning all reads back to entire assembly and splitting reads into individual fastq files based on their bin membership
    mkdir -p "$OUT_DIR/Binning/reads_for_reassembly"
    bwa mem -t 104 "$OUT_DIR/Binning/all_bins_index" $Read1 $Read2 \
      | python "${ScriptDir}/filter_reads_for_bin_reassembly.py" "$OUT_DIR/Binning/vRhyme_results_${BASENAME}_filtered/vRhyme_best_bins_fasta" "$OUT_DIR/Binning/reads_for_reassembly" $STRICT_MAX $PERMISSIVE_MAX

    # 重新组装vRhyme_results_${BASENAME}_filtered/vRhyme_best_bins_fasta文件夹里的所有fasta文件
    for FASTQ_FILE in "$OUT_DIR/Binning/reads_for_reassembly/"*_1.fastq; do
      BIN_BASENAME=$(basename "$FASTQ_FILE" _1.fastq)
      OriginalBin=${BIN_BASENAME%%.*}

      # 创建目录
      EXTRACTED_DIR="$OUT_DIR/Binning/reassembled_bins"
      TMP_DIR="$EXTRACTED_DIR/${BIN_BASENAME}_tmp"
      mkdir -p "$TMP_DIR"

      # 使用SPAdes进行bin的重新组装
      spades.py -t 104 --tmp $TMP_DIR --careful \
        --untrusted-contigs "$OUT_DIR/Binning/vRhyme_results_${BASENAME}_filtered/vRhyme_best_bins_fasta/${OriginalBin}.fasta" \
        -1 "${OUT_DIR}/Binning/reads_for_reassembly/${BIN_BASENAME}_1.fastq" \
        -2 "${OUT_DIR}/Binning/reads_for_reassembly/${BIN_BASENAME}_2.fastq" \
        -o "$OUT_DIR/Binning/reassembile/${BIN_BASENAME}"

      # 清理短的contigs
      # python rm_short_contigs.py "$OUT_DIR/Binning/reassembile/${BIN_BASENAME}/contigs.fasta" 500

      # 提取并重命名scaffolds.fasta
      mkdir -p "$EXTRACTED_DIR"
      cp "$OUT_DIR/Binning/reassembile/${BIN_BASENAME}/contigs.fasta" "$EXTRACTED_DIR/${BIN_BASENAME}.fasta"
      cp "$OUT_DIR/Binning/vRhyme_results_${BASENAME}_filtered/vRhyme_best_bins_fasta/${OriginalBin}.fasta" "$EXTRACTED_DIR/${OriginalBin}.origin.fasta"

      echo "Reassembly complete for ${BIN_BASENAME}"

      # 删除临时目录
      rm -rf "$TMP_DIR"
    done
        
    # Find best bin
    mkdir "$OUT_DIR/Binning/Summary"
    python "${ScriptDir}/concat_fasta_sequences.py" $EXTRACTED_DIR "$OUT_DIR/Binning/Summary/tempsummary.fasta"
    checkv end_to_end "$OUT_DIR/Binning/Summary/tempsummary.fasta" "$OUT_DIR/Binning/Summary/CheckRes" -t 100 -d "$DATABASE/checkv-db-v1.5"

    # Select best bin
    mkdir "$OUT_DIR/Binning/Summary/Finialbins"
    python "${ScriptDir}/select_best_bins.py" "$OUT_DIR/Binning/Summary/CheckRes/quality_summary.tsv" $EXTRACTED_DIR "$OUT_DIR/Binning/Summary/Finialfasta/Bestbins"

  else
    # 如果不重新组装，则直接将vRhyme_best_bins_fasta目录中的所有bin放入Bestbins文件夹中
    mkdir "$OUT_DIR/Binning/Summary"
    mkdir -p "$OUT_DIR/Binning/Summary/Finialfasta/Bestbins"
    cp "$OUT_DIR/Binning/vRhyme_results_${BASENAME}_filtered/vRhyme_best_bins_fasta/"*.fasta "$OUT_DIR/Binning/Summary/Finialfasta/Bestbins/"
  fi

  EXTRACTED_DIR="$OUT_DIR/Binning/reassembled_bins"
  # 去掉$OUT_DIR/Binning/Summary/Finialfasta/Bestbins文件夹内所有fasta文件名字里的vRhyme_
  for vRhymeFILE in "$OUT_DIR/Binning/Summary/Finialfasta/Bestbins"/*.fasta; do
    NEW_NAME=$(basename "$vRhymeFILE" | sed 's/^vRhyme_//')
    mv "$vRhymeFILE" "$OUT_DIR/Binning/Summary/Finialfasta/Bestbins/$NEW_NAME"
  done

  # 提取未归类的contigs
  python "${ScriptDir}/unbined.py" -i "$OUT_DIR/Binning/vRhyme_results_${BASENAME}_filtered/vRhyme_best_bins_fasta" -r "$OUT_DIR/${BASENAME}_filtered.fasta" -o "$OUT_DIR/Binning/Summary/Finialfasta/unbined.fasta"

  echo "Rebinning and reassembly complete for $FILE"
done
