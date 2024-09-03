#!/usr/bin/env bash

# 进行 Assembly
for FILE in $FILES; do
  echo "Processing $FILE"

  BASENAME=$(basename "$FILE" .fa)
  BASENAME=${BASENAME%.fasta}
  OUT_DIR="$OUTPUT_DIR/SeprateFile/${BASENAME}"
  PREDICTION_DIR="$OUT_DIR/RoughViralPrediction"
  Genomad_dir="$PREDICTION_DIR/genomadres"
  Viralverify_dir="$PREDICTION_DIR/viralverify"
  Virsorter_dir="$PREDICTION_DIR/virsorter2"
  Read1="$RAW_SEQ_DIR/${BASENAME}_R1.fq"
  Read2="$RAW_SEQ_DIR/${BASENAME}_R2.fq"

  # Assembly
  echo -e "\n \n \n # 进行fastap处理交叉验证!!! \n \n \n"
  fastp -i in.R1.fq.gz -I in.R2.fq.gz -o out.R1.fq.gz -O out.R2.fq.gz

  spades.py -1 out.R1.fq.gz -2 out.R2.fq.gz -o output_spades #contigs scaffold
  megahit -1 out.R1.fq.gz -2 out.R2.fq.gz -o megahit_out --presets meta-large #extraction
done