#!/usr/bin/env bash

# 设置错误捕获机制：如果发生错误，脚本将停止执行
set -e
trap 'echo "An error occurred. Exiting..."; exit 1;' ERR

# 确保 BASE_CONDA_PREFIX 被正确设置
BASE_CONDA_PREFIX=$(conda info --base)
conda_sh="$BASE_CONDA_PREFIX/etc/profile.d/conda.sh"

# 导出必要的变量和函数
export OUTPUT_DIR DATABASE Group CONCENTRATION_TYPE

# 处理文件的主函数
process_file() {
  local FILE=$1
  local BASENAME=$(basename "$FILE" .fa)
  BASENAME=${BASENAME%.fasta}

  local OUT_DIR="$OUTPUT_DIR/SeprateFile/${BASENAME}"
  mkdir -p "$OUT_DIR"

  local PREDICTION_DIR="$OUT_DIR/RoughViralPrediction"
  mkdir -p "$PREDICTION_DIR"

  local Genomad_dir="$PREDICTION_DIR/genomadres"
  mkdir -p "$Genomad_dir"
  local Viralverify_dir="$PREDICTION_DIR/viralverify"
  mkdir -p "$Viralverify_dir"

  echo "Processing $FILE"

  # 执行 viralverify 分析
  if [ ! -f "$Viralverify_dir/${BASENAME}_result_table.csv" ]; then
    echo "Running viralverify prediction..."
    source ${conda_sh}
    viralverify -f "$FILE" -o "$Viralverify_dir" --hmm "$DATABASE/ViralVerify/nbc_hmms.hmm" -t 104 > "$Viralverify_dir/viralverify.log" 2>&1
    echo "viralverify prediction completed!"
  else
    echo "viralverify prediction already completed for $FILE, skipping..."
  fi

  # 仅当是 concentration 时才运行 Virsorter2
  if [ "$CONCENTRATION_TYPE" == "concentration" ]; then
    local Virsorter_dir="$PREDICTION_DIR/virsorter2"
    mkdir -p "$Virsorter_dir"

    if [ ! -f "$Virsorter_dir/final-viral-score.tsv" ]; then
      echo "Running Virsorter2 prediction..."
      virsorter run -w "$Virsorter_dir" -i "$FILE" --include-groups "$Group" -j 104 all --min-score 0.5 --min-length 2000 --keep-original-seq -d "$DATABASE/db" > "$Virsorter_dir/virsorter.log" 2>&1
      echo "Virsorter2 prediction completed!"
    else
      echo "Virsorter2 prediction already completed for $FILE, skipping..."
    fi
  fi
}

export -f process_file
parallel process_file ::: $FILES

# 进行 Genomad 分析
for FILE in $RawFILES; do
  echo "Processing $FILE"
  BASENAME=$(basename "$FILE" .fa)
  BASENAME=${BASENAME%.fasta}

  OUT_DIR="$OUTPUT_DIR/SeprateFile/${BASENAME}"
  PREDICTION_DIR="$OUT_DIR/RoughViralPrediction"

  # genomad 预测
  echo -e "\n \n \n # 进行genomad预测!!! \n \n \n"
  Genomad_dir="$PREDICTION_DIR/genomadres"
  mkdir -p "$Genomad_dir"

  if [ ! -f "$Genomad_dir/${BASENAME}_summary/${BASENAME}_virus_summary.tsv" ]; then
    genomad end-to-end --enable-score-calibration "$FILE" "$Genomad_dir" "$DATABASE/genomad_db"
    echo -e "\n \n \n # genomad预测完成!!! \n \n \n"
  else
    echo "genomad prediction already completed for $FILE, skipping..."
  fi
done

# 检测所有 Virsorter2 任务是否完成 (仅当 concentration 时)
if [ "$CONCENTRATION_TYPE" == "concentration" ]; then
  all_tasks_completed=false
  while [ "$all_tasks_completed" == "false" ]; do
    all_tasks_completed=true
    for FILE in $FILES; do
      BASENAME=$(basename "$FILE" .fa)
      BASENAME=${BASENAME%.fasta}
      Virsorter_dir="$OUTPUT_DIR/SeprateFile/${BASENAME}/RoughViralPrediction/virsorter2"

      if [ ! -f "$Virsorter_dir/final-viral-score.tsv" ]; then
        all_tasks_completed=false
        echo "Virsorter2 still in processing"
        break
      fi
    done

    if [ "$all_tasks_completed" == "false" ]; then
      sleep 30
    fi
  done
fi

echo "All files have been processed."
