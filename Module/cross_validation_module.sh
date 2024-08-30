#!/usr/bin/env bash

# 进行CrossValid分析
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

  # CrossValid virus contigs
  echo -e "\n \n \n # 进行virus contigs交叉验证!!! \n \n \n"
  python "${ScriptDir}/CrossValid.py" "$Genomad_dir" "$Viralverify_dir" "$Virsorter_dir" "$BASENAME" "$OUT_DIR"

  # Sequence提取
  python "${ScriptDir}/FilterRawResSeqs.py" "$FILE" "$BASENAME" "$OUT_DIR"
done
