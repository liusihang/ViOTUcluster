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
VERSION="0.0.3"  # 定义版本号

# 设置日志文件路径
LOG_FILE=""

# 使用 getopts 处理命令行选项
while getopts ":i:r:o:d:t:rvh" opt; do
  case $opt in
    i)
      INPUT_DIR=$OPTARG
      ;;
    r)
      RAW_SEQ_DIR=$OPTARG
      ;;
    o)
      OUTPUT_DIR=$OPTARG
      LOG_FILE="${OUTPUT_DIR}/pipeline.log"  # 设置日志文件路径
      ;;
    d)
      DATABASE=$OPTARG
      ;;
    t) # 用 't' 代表 'type'
      SAMPLETYPE=$OPTARG
      ;;
    v) # 显示版本信息
      echo "Version: $VERSION"
      exit 0
      ;;
    h) # 显示帮助信息
      echo "Usage: $0 [options]"
      echo ""
      echo "Options:"
      echo "  -i <input_path_to_search>   Specify the input directory to search for FASTA files."
      echo "  -r <input_path_raw_seqs>    Specify the input directory to search for raw seqs files."
      echo "  -o <output_path>            Specify the output directory for the results."
      echo "  -d <database_path>          Specify the path to the database required for analysis."
      echo "  -t <sample_type>            Specify the sample type: DNA, RNA, or Mix."
      echo "  --reassemble                Enable reassembly of bins."
      echo "  -v                          Display the version of this script."
      echo "  -h                          Display this help and exit."
      exit 0
      ;;
    -)
      case "${OPTARG}" in
        reassemble)
          REASSEMBLE=true
          ;;
        *)
          echo "Unknown option: --${OPTARG}" >&2
          exit 1
          ;;
      esac
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
    :)
      echo "Option -$OPTARG requires an argument." >&2
      exit 1
      ;;
  esac
done

# 检查必要的输入是否已提供
if [ -z "$INPUT_DIR" ] || [ -z "$RAW_SEQ_DIR" ] || [ -z "$OUTPUT_DIR" ] || [ -z "$DATABASE" ] || [ -z "$SAMPLETYPE" ]; then
    echo "Usage: $0 -i <input_path_to_search> -r <input_path_raw_seqs> -o <output_path> -d <database_path> -t <sample_type> [--reassemble]"
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

# 获取当前 Conda 环境的 bin 文件夹位置
if [ -z "$CONDA_PREFIX" ]; then
  echo "Conda environment is not activated."
  exit 1
fi
ScriptDir="${CONDA_PREFIX}/bin"


# 导出参数作为环境变量
export INPUT_DIR OUTPUT_DIR DATABASE SAMPLETYPE REASSEMBLE

# 获取当前 Conda 环境的 bin 文件夹位置
if [ -z "$CONDA_PREFIX" ]; then
  echo "Conda environment is not activated."
  exit 1
fi
ScriptDir="${CONDA_PREFIX}/bin"

# 过滤小于2000bp的序列
mkdir -p "${OUTPUT_DIR}/FilteredSeqs"
python "${ScriptDir}/filter_contigs.py" 100 "${INPUT_DIR}" "${OUTPUT_DIR}/FilteredSeqs"

# Find all .fa and .fasta files in the specified directory
FILES=$(find "${OUTPUT_DIR}/FilteredSeqs" -type f \( -name "*.fa" -o -name "*.fasta" \))
RawFILES=$(find "${INPUT_DIR}" -maxdepth 1 -type f \( -name "*.fa" -o -name "*.fasta" \))

# 导出必要的变量和函数
export OUTPUT_DIR DATABASE Group FILES RawFILES

# 按顺序执行各个模块并记录时间和日志

mkdir -p "${OUTPUT_DIR}/Log"

echo "Starting viral prediction..."
#Viral prediction module
module_start_time=$(date +%s)
./viral_prediction_module.sh > "${OUTPUT_DIR}/Log/viral_prediction.log" 2>&1
module_end_time=$(date +%s)
module_runtime=$((module_end_time - module_start_time))
echo "Viral prediction completed in ${module_runtime} seconds." >> "${OUTPUT_DIR}/Log/viral_prediction.log"

#Cross Validation module
echo "Starting Cross Validation module..."
module_start_time=$(date +%s)
./cross_validation_module.sh > "${OUTPUT_DIR}/Log/cross_validation.log" 2>&1
module_end_time=$(date +%s)
module_runtime=$((module_end_time - module_start_time))
echo "Cross Validation completed in ${module_runtime} seconds." >> "${OUTPUT_DIR}/Log/cross_validation.log"

#Binning and merge module
echo "Starting Binning and merge module..."
module_start_time=$(date +%s)
./binning_merge_module.sh > "${OUTPUT_DIR}/Log/binning_merge.log" 2>&1
module_end_time=$(date +%s)
module_runtime=$((module_end_time - module_start_time))
echo "Binning completed in ${module_runtime} seconds." >> "${OUTPUT_DIR}/Log/binning.log"

#Summary module
echo "Starting Summary module..."
module_start_time=$(date +%s)
./summary_module.sh > "${OUTPUT_DIR}/Log/summary.log" 2>&1
module_end_time=$(date +%s)
module_runtime=$((module_end_time - module_start_time))
echo "Summary completed in ${module_runtime} seconds." >> "${OUTPUT_DIR}/Log/summary.log"

#dRep module
echo "Starting dRep module..."
module_start_time=$(date +%s)
./drep_module.sh > "${OUTPUT_DIR}/Log/drep.log" 2>&1
module_end_time=$(date +%s)
module_runtime=$((module_end_time - module_start_time))
echo "dRep completed in ${module_runtime} seconds." >> "${OUTPUT_DIR}/Log/drep.log"

# 记录结束时间并计算整个流程的总耗时
end_time=$(date +%s)
total_runtime=$((end_time - start_time))
echo "Total runtime: ${total_runtime} seconds" >> "${OUTPUT_DIR}/pipeline.log"