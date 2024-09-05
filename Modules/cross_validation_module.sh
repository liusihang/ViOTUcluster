#!/usr/bin/env bash

# 导出必要的变量
export OUTPUT_DIR DATABASE CONCENTRATION_TYPE ScriptDir

process_file() {
  local FILE=$1
  echo "Processing $FILE"

  BASENAME=$(basename "$FILE" .fa)
  BASENAME=${BASENAME%.fasta}
  OUT_DIR="$OUTPUT_DIR/SeprateFile/${BASENAME}"
  PREDICTION_DIR="$OUT_DIR/RoughViralPrediction"
  Genomad_dir="$PREDICTION_DIR/genomadres"
  Viralverify_dir="$PREDICTION_DIR/viralverify"
  Virsorter_dir="$PREDICTION_DIR/virsorter2"

  # 如果已经存在 quality_summary.tsv 文件，则跳过该文件
  if [ -f "$OUT_DIR/${BASENAME}_CheckRes/quality_summary.tsv" ]; then
    echo "Skipping $BASENAME as quality_summary.tsv already exists."
    return 0
  fi

  # CrossValid virus contigs
  echo -e "\n\n\n # Performing cross-validation for virus contigs!!! \n\n\n"
  if ! python "${ScriptDir}/CrossValid.py" "$Genomad_dir" "$Viralverify_dir" "$Virsorter_dir" "$BASENAME" "$OUT_DIR" "$CONCENTRATION_TYPE"; then
    echo "Error during cross-validation for $BASENAME. Exiting..."
    return 1
  fi

  # Extract Sequences
  if ! python "${ScriptDir}/FilterRawResSeqs.py" "$FILE" "$BASENAME" "$OUT_DIR"; then
    echo "Error during sequence extraction for $BASENAME. Exiting..."
    return 1
  fi

  # CheckVFilter
  mkdir -p "$OUT_DIR/${BASENAME}_CheckRes"
  if ! checkv end_to_end "$OUT_DIR/${BASENAME}_filtered.fasta" "$OUT_DIR/${BASENAME}_CheckRes" -t 104 -d "$DATABASE/checkv-db-v1.5"; then
    echo "Error during CheckVFilter for $BASENAME. Exiting..."
    return 1
  fi

  if ! python "${ScriptDir}/check_removal.py" "$OUT_DIR/${BASENAME}_CheckRes/quality_summary.tsv" "$OUT_DIR/${BASENAME}_filtered.fasta"; then
    echo "Error during check removal for $BASENAME. Exiting..."
    return 1
  fi

  echo "Processing for $BASENAME completed."
}

export -f process_file

# 并行运行
parallel process_file ::: $FILES

echo "CrossValid analysis completed."