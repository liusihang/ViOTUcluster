#!/usr/bin/env bash

# 导出必要的变量
export OUTPUT_DIR DATABASE CONCENTRATION_TYPE ScriptDir FILES THREADS

set -x  # 启用调试模式

process_file() {
  local FILE=$1
  echo "Processing file: $FILE"
  local BASENAME=$(basename "$FILE")
  # 提取扩展名
  local EXTENSION="${BASENAME##*.}"
  # 如果扩展名是 .fa 或 .fasta，则移除扩展名
  if [ "$EXTENSION" = "fa" ] || [ "$EXTENSION" = "fasta" ]; then
      BASENAME="${BASENAME%.*}"
  fi
  echo "Basename: $BASENAME"

  out_dir="$OUTPUT_DIR/SeprateFile/${BASENAME}"
  PREDICTION_DIR="$out_dir/RoughViralPrediction"
  genomad_dir="$PREDICTION_DIR/genomadres"
  viralverify_dir="$PREDICTION_DIR/viralverify"
  virsorter_dir="$PREDICTION_DIR/virsorter2"

  # 如果已经存在 quality_summary.tsv 文件，则跳过该文件
  if [ -f "$out_dir/${BASENAME}_filtered.fasta" ]; then
    echo "Skipping $BASENAME as quality_summary.tsv already exists."
    return 0
  fi

  # 交叉验证病毒片段
  echo -e "\n\n\n # Performing cross-validation for virus contigs!!! \n\n\n"
  echo "Running CrossValid.py..."
  if ! python "${ScriptDir}/CrossValidAnA.py" "$genomad_dir" "$viralverify_dir" "$virsorter_dir" "$BASENAME" "$out_dir" "$CONCENTRATION_TYPE"; then
    echo "Error during cross-validation for $BASENAME. Exiting..."
    return 1
  fi

  # 从原始结果中提取序列
  echo "Running FilterRawResSeqs.py..."
  if ! python "${ScriptDir}/FilterRawResSeqs.py" "$FILE" "$BASENAME" "$out_dir"; then
    echo "Error during sequence extraction for $BASENAME. Exiting..."
    return 1
  fi

  return 0

  # 创建用于 CheckV 结果的目录
  echo "Creating CheckV results directory..."
  mkdir -p "$out_dir/${BASENAME}_CheckRes"

  # 运行 CheckV 以评估病毒序列的质量
  echo "Running CheckV..."
  if ! checkv end_to_end "$out_dir/${BASENAME}_filtered.fasta" "$out_dir/${BASENAME}_CheckRes" -t "${THREADS}" -d "$DATABASE/checkv-db-v1.5"; then
    echo "Error during CheckV for $BASENAME. Exiting..."
    return 1
  fi

  # 运行 check_removal.py 根据质量报告移除低质量序列
  echo "Running check_removal.py..."
  if ! python "${ScriptDir}/check_removal.py" "$out_dir/${BASENAME}_CheckRes/quality_summary.tsv" "$out_dir/${BASENAME}_filtered.fasta"; then
    echo "Error during check removal for $BASENAME. Exiting..."
    return 1
  fi

  echo "Processing for $BASENAME completed."

  # 删除 ${BASENAME}_CheckRes 目录以释放空间
  echo "Removing CheckV results directory..."
  rm -rf "$out_dir/${BASENAME}_CheckRes"
}

export -f process_file

# 并行运行
echo "Starting parallel processing..."
parallel process_file ::: $FILES

echo "CrossValid analysis completed."