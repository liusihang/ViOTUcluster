#!/usr/bin/env bash

# 记录开始时间
start_time=$(date +%s)

# 初始化变量
INPUT_DIR=""
RAW_SEQ_DIR=""
OUTPUT_DIR=""
DATABASE=""
SAMPLETYPE=""
REASSEMBLE=false
VERSION="0.0.5"  # 定义版本号
CONCENTRATION_TYPE=""

# 设置日志文件路径
LOG_FILE=""

# 使用 getopts 处理短选项和手动处理长选项
while [[ $# -gt 0 ]]; do
  case $1 in
    -i)
      INPUT_DIR=$2
      shift 2
      ;;
    -r)
      RAW_SEQ_DIR=$2
      shift 2
      ;;
    -o)
      OUTPUT_DIR=$2
      LOG_FILE="${OUTPUT_DIR}/pipeline.log"  # 设置日志文件路径
      shift 2
      ;;
    -d)
      DATABASE=$2
      shift 2
      ;;
    -t) # 用 't' 代表 'type'
      SAMPLETYPE=$2
      shift 2
      ;;
    -v) # 显示版本信息
      echo "Version: $VERSION"
      exit 0
      ;;
    --non-con) # non-concentration
      CONCENTRATION_TYPE="non-concentration"
      shift
      ;;
    --con) # concentration
      CONCENTRATION_TYPE="concentration"
      shift
      ;;
    -h) # 显示帮助信息
      echo "Usage: $0 [options]"
      echo ""
      echo "Options:"
      echo "  -i <input_path_to_search>   Specify the input directory to search for FASTA files."
      echo "  -r <input_path_raw_seqs>    Specify the input directory to search for raw seqs files."
      echo "  -o <output_path>            Specify the output directory for the results."
      echo "  -d <database_path>          Specify the path to the database required for analysis."
      echo "  -t <sample_type>            Specify the sample type: DNA, RNA, or Mix."
      echo "  --non-con                   Specify non-concentration processing."
      echo "  --con                       Specify concentration processing."
      echo "  --reassemble                Enable reassembly of bins."
      echo "  -v                          Display the version of this script."
      echo "  -h                          Display this help and exit."
      exit 0
      ;;
    --reassemble)
      REASSEMBLE=true
      shift
      ;;
    *)
      echo "Unknown option: $1" >&2
      exit 1
      ;;
  esac
done

# 检查必要的输入是否已提供
if [ -z "$INPUT_DIR" ] || [ -z "$RAW_SEQ_DIR" ] || [ -z "$OUTPUT_DIR" ] || [ -z "$DATABASE" ] || [ -z "$SAMPLETYPE" ] || [ -z "$CONCENTRATION_TYPE" ]; then
    echo "Usage: $0 -i <input_path_to_search> -r <input_path_raw_seqs> -o <output_path> -d <database_path> -t <sample_type> --non-con/--con [--reassemble]"
    exit 1
fi

# 验证配对文件是否存在并且格式正确
for FILE in "${RAW_SEQ_DIR}"/*_R1.*; do
  BASENAME=$(basename "$FILE" | sed 's/_R1\..*//')
  PREFIX="${RAW_SEQ_DIR}/${BASENAME}"

  if [ -f "${PREFIX}_R1.fq" ] && [ -f "${PREFIX}_R2.fq" ]; then
    echo "Found paired files: ${PREFIX}_R1.fq and ${PREFIX}_R2.fq"
  elif [ -f "${PREFIX}_R1.fastq" ] && [ -f "${PREFIX}_R2.fastq" ]; then
    echo "Found paired files: ${PREFIX}_R1.fastq and ${PREFIX}_R2.fastq"
  elif [ -f "${PREFIX}_R1.fq.gz" ] && [ -f "${PREFIX}_R2.fq.gz" ]; then
    echo "Found paired files: ${PREFIX}_R1.fq.gz and ${PREFIX}_R2.fq.gz"
  elif [ -f "${PREFIX}_R1.fastq.gz" ] && [ -f "${PREFIX}_R2.fastq.gz" ]; then
    echo "Found paired files: ${PREFIX}_R1.fastq.gz and ${PREFIX}_R2.fastq.gz"
  else
    echo "Error: Paired-end files for ${BASENAME} not found in the expected formats (.fq, .fastq, .fq.gz, .fastq.gz)"
    exit 1
  fi
done

# 根据 CONCENTRATION_TYPE 进行不同的处理
if [ "$CONCENTRATION_TYPE" == "non-concentration" ]; then
    echo "Running non-concentration specific steps..."
    # 在这里添加 non-concentration 的特定处理步骤
elif [ "$CONCENTRATION_TYPE" == "concentration" ]; then
    echo "Running concentration specific steps..."
    # 在这里添加 concentration 的特定处理步骤
else
    echo "Error: Invalid concentration type."
    exit 1
fi

# 重定向所有输出到日志文件
exec > >(tee -a "$LOG_FILE") 2>&1

# 根据 SAMPLETYPE 设置 Group 变量
case $SAMPLETYPE in
  DNA)
    Group="dsDNAphage, NCLDV, ssDNA, lavidaviridae"
    ;;
  RNA)
    Group="RNA, lavidaviridae"
    ;;
  Mix)
    Group="dsDNAphage, NCLDV, RNA, ssDNA, lavidaviridae"
    ;;
  *)
    echo "Unknown sample type: $SAMPLETYPE"
    exit 1
    ;;
esac

# 输出选择的处理类型
echo "Processing with $CONCENTRATION_TYPE mode."

# 获取当前 Conda 环境的 bin 文件夹位置
if [ -z "$CONDA_PREFIX" ]; then
  echo "Conda environment is not activated."
  exit 1
fi
ScriptDir="${CONDA_PREFIX}/bin"

# 导出参数作为环境变量
export INPUT_DIR OUTPUT_DIR DATABASE SAMPLETYPE REASSEMBLE ScriptDir RAW_SEQ_DIR

# 过滤小于2000bp的序列
mkdir -p "${OUTPUT_DIR}/FilteredSeqs"
python "${ScriptDir}/filter_contigs.py" 100 "${INPUT_DIR}" "${OUTPUT_DIR}/FilteredSeqs"

# Find all .fa and .fasta files in the specified directory
FILES=$(find "${OUTPUT_DIR}/FilteredSeqs" -type f \( -name "*.fa" -o -name "*.fasta" \))
RawFILES=$(find "${INPUT_DIR}" -maxdepth 1 -type f \( -name "*.fa" -o -name "*.fasta" \))

# 导出必要的变量和函数
export OUTPUT_DIR DATABASE Group FILES RawFILES CONCENTRATION_TYPE ScriptDir RAW_SEQ_DIR

# 按顺序执行各个模块并记录时间和日志

mkdir -p "${OUTPUT_DIR}/Log"

echo "Starting viral prediction..."
#Viral prediction module
module_start_time=$(date +%s)
# 根据 CONCENTRATION_TYPE 进行不同的处理
viral_prediction_module.sh > "${OUTPUT_DIR}/Log/cross_validation.log" 2>&1
module_end_time=$(date +%s)
module_runtime=$((module_end_time - module_start_time))
echo "Viral prediction completed in ${module_runtime} seconds." >> "${OUTPUT_DIR}/Log/viral_prediction.log"


# Cross Validation module
echo "Starting Cross Validation module..."
# 检查日志文件是否存在以及是否包含关键字
LOG_FILE="${OUTPUT_DIR}/Log/cross_validation.log"
KEYWORD="Cross Validation completed"

if [ -f "$LOG_FILE" ] && grep -q "$KEYWORD" "$LOG_FILE"; then
    echo "Cross Validation module already completed, skipping..."
else
    module_start_time=$(date +%s)
    # 根据 CONCENTRATION_TYPE 进行不同的处理
    if [ "$CONCENTRATION_TYPE" == "non-concentration" ]; then
        cross_validation_module.sh --non-concentration > "$LOG_FILE" 2>&1
    else
        cross_validation_module.sh --concentration > "$LOG_FILE" 2>&1
    fi
    module_end_time=$(date +%s)
    module_runtime=$((module_end_time - module_start_time))
    echo "Cross Validation completed in ${module_runtime} seconds." >> "$LOG_FILE"
fi


#Binning and merge module
echo "Starting Binning and merge module..."
module_start_time=$(date +%s)
binning_merge_module.sh > "${OUTPUT_DIR}/Log/binning_merge.log" 2>&1
module_end_time=$(date +%s)
module_runtime=$((module_end_time - module_start_time))
echo "Binning completed in ${module_runtime} seconds." >> "${OUTPUT_DIR}/Log/binning_merge.log"

#Summary module
echo "Starting Summary module..."
module_start_time=$(date +%s)
summary_module.sh > "${OUTPUT_DIR}/Log/summary.log" 2>&1
module_end_time=$(date +%s)
module_runtime=$((module_end_time - module_start_time))
echo "Summary completed in ${module_runtime} seconds." >> "${OUTPUT_DIR}/Log/summary.log"

# dRep module
echo "Starting dRep module..."
LOG_FILE="${OUTPUT_DIR}/Log/drep.log"
# 检查日志文件是否存在且包含关键字
if [ -f "$LOG_FILE" ] && grep -q "Combined fasta files and quality summaries completed." "$LOG_FILE"; then
  echo "dRep module already completed, skipping..."
else
  module_start_time=$(date +%s)
  # 执行 dRep 模块
  drep_module.sh > "$LOG_FILE" 2>&1
  module_end_time=$(date +%s)
  module_runtime=$((module_end_time - module_start_time))
  echo "dRep completed in ${module_runtime} seconds." >> "$LOG_FILE"
fi

#TPM caculate module
LOG_FILE="${OUTPUT_DIR}/Log/TPM_caculate.log"
# 检查日志文件是否存在，并且是否包含 "TPM calculation completed successfully."
if [ -f "$LOG_FILE" ] && grep -q "TPM calculation completed successfully." "$LOG_FILE"; then
    echo "TPM calculation already completed, skipping..."
else
    echo "Starting TPM calculate module..."
    module_start_time=$(date +%s)
    # 执行TPM计算模块并将日志输出到文件
    TPM_caculate_Module.sh > "$LOG_FILE" 2>&1
    module_end_time=$(date +%s)
    module_runtime=$((module_end_time - module_start_time))
    # 将运行时间追加到日志文件中
    echo "TPM calculate completed in ${module_runtime} seconds." >> "$LOG_FILE"
fi


#DRAM module
echo "Starting DRAM module..."
module_start_time=$(date +%s)
mkdir -p "$OUTPUT_DIR/Summary/DRAM"
run_dram_analysis.sh "$OUTPUT_DIR/Summary/Viralcontigs/vOTU.fasta" "$OUTPUT_DIR/Summary/DRAM" > "${OUTPUT_DIR}/Log/DRAM.log" 2>&1
module_end_time=$(date +%s)
module_runtime=$((module_end_time - module_start_time))
echo "DRAM completed in ${module_runtime} seconds." >> "${OUTPUT_DIR}/Log/DRAM.log"

#iPhop module
echo "Starting iPhop module..."
module_start_time=$(date +%s)
mkdir -p "$OUTPUT_DIR/Summary/iPhop"
run_iphop_analysis.sh "$OUTPUT_DIR/Summary/Viralcontigs/vOTU.fasta" "$OUTPUT_DIR/Summary/iPhop" > "${OUTPUT_DIR}/Log/iPhop.log" 2>&1
module_end_time=$(date +%s)
module_runtime=$((module_end_time - module_start_time))
echo "iPhop completed in ${module_runtime} seconds." >> "${OUTPUT_DIR}/Log/iPhop.log"

# 记录结束时间并计算整个流程的总耗时
end_time=$(date +%s)
total_runtime=$((end_time - start_time))
echo "Total runtime: ${total_runtime} seconds" >> "${OUTPUT_DIR}/pipeline.log"
