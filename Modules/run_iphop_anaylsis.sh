#!/bin/bash

# 检查命令行参数
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input_fasta_file> <output_directory>"
    exit 1
fi

# 接受命令行参数
INPUT_FASTA=$1
OUTPUT_DIR=$2

# 检查输入文件是否存在
if [ ! -f "$INPUT_FASTA" ]; then
    echo "Error: Input FASTA file '$INPUT_FASTA' not found."
    exit 1
fi

# 创建输出目录（如果不存在）
mkdir -p "$OUTPUT_DIR"
mkdir -p "$OUTPUT_DIR/split_files"
mkdir -p "$OUTPUT_DIR/iPhop_results"

echo -e "\n\n\n# 进行Host prediction分析!!!\n\n\n"
pwd

# 获取输入文件的基础文件名（不带路径）
BASE_INPUT_FASTA=$(basename "$INPUT_FASTA")

# 分割输入的fasta文件，每个新文件包含3000条序列
awk -v output_dir="$OUTPUT_DIR/split_files" -v base_name="$BASE_INPUT_FASTA" 'BEGIN {n_seq=0;} 
     /^>/ {
        if (n_seq % 3000 == 0) {
            file = sprintf("%s/%s_%d.fna", output_dir, base_name, n_seq);
        } 
        print >> file; 
        n_seq++; 
        next;
     } 
     { print >> file; }' "$INPUT_FASTA"

# 切换到分割后的文件目录
cd "$OUTPUT_DIR/split_files" || exit

# 列出所有分割后的fna文件
ls *.fna > iPhop

# 设置conda环境和日志文件
BASE_CONDA_PREFIX=$(conda info --base)
conda_sh="$BASE_CONDA_PREFIX/etc/profile.d/conda.sh"
CURRENT_ENV=$(basename "$CONDA_DEFAULT_ENV")

# 激活vRhyme环境并运行vRhyme
source ${conda_sh}
conda activate iphop
echo "Conda environment activated: $(conda info --envs)"
which iphop

# 使用 GNU Parallel 并行执行 iPhop host prediction
cat iPhop | parallel -j 2 'fa_file={}; iphop predict --fa_file "${fa_file}" --db_dir /media/em/student/db/Virus/Aug_2023_pub_rw --out_dir "${fa_file}_iPhopResult" -t 10'

# 监控任务是否完成
all_tasks_completed=false
while [ "$all_tasks_completed" == "false" ]; do
    sleep 30
    all_tasks_completed=true
    if ls *_iPhopResult/Host_prediction_to_genus_m90.csv 1> /dev/null 2>&1; then
        echo "iPhop prediction still in progress."
        all_tasks_completed=false
    fi
    if [ "$all_tasks_completed" == "false" ]; then
        sleep 30
    fi
done

# 合并所有注释结果到输出目录
awk 'FNR==1 && NR!=1{next;} {print}' ./*_iPhopResult/Host_prediction_to_genus_m90.csv > "$OUTPUT_DIR/Combined_Host_prediction_to_genus_m90.tsv"

echo "Host prediction complete. Results combined and saved to $OUTPUT_DIR/Combined_Host_prediction_to_genus_m90.tsv"

# 删除中间产生的临时文件
echo "Cleaning up temporary files..."
rm -rf "$OUTPUT_DIR/split_files"
rm -rf "$OUTPUT_DIR/iPhop_results"/*_iPhopResult

echo "Cleanup complete."
conda activate "$CURRENT_ENV"