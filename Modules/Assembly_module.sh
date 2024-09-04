#!/usr/bin/env bash

# 进行 Assembly
for FILE in $FILES; do
  echo "Processing $FILE"

  BASENAME=$(basename "$FILE" .fa)
  BASENAME=${BASENAME%.fasta}
  OUT_DIR="$OUTPUT_DIR/SeprateFile/${BASENAME}"
  # 检查并设置Read1和Read2的路径
  # 匹配BASENAME_R1* 的格式
  Read1=$(find "${RAW_SEQ_DIR}" -type f -name "${BASENAME}_R1*" | head -n 1)
  Read2=$(find "${RAW_SEQ_DIR}" -type f -name "${BASENAME}_R2*" | head -n 1)

  # Assembly
  echo -e "\n \n \n # 进行fastap处理交叉验证!!! \n \n \n"
  fastp -i in.R1.fq.gz -I in.R2.fq.gz -o out.R1.fq.gz -O out.R2.fq.gz

  spades.py -1 out.R1.fq.gz -2 out.R2.fq.gz -o output_spades #contigs scaffold
  megahit -1 out.R1.fq.gz -2 out.R2.fq.gz -o megahit_out --presets meta-large #extraction
done