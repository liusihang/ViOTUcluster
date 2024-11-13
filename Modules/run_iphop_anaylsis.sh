#!/bin/bash

source activate iphop

echo "Conda environment activated: $(conda info --envs)"

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

# 使用 Python 控制 iPhop 预测
python "${ScriptDir}/run_iphop.py"


# 定义一个函数来检查内存使用情况
check_memory_usage() {
    # 获取总内存和已用内存（以 KB 为单位）
    total_memory=$(free -k | awk '/^Mem:/ {print $2}')
    used_memory=$(free -k | awk '/^Mem:/ {print $3}')
    
    # 计算内存使用百分比
    memory_usage=$(awk -v used="$used_memory" -v total="$total_memory" 'BEGIN { print (used / total) * 100 }')

    # 如果内存使用超过 99.5%，则返回 1，否则返回 0
    if awk 'BEGIN {exit ARGV[1] <= 99.5}' "$memory_usage"; then
        return 1  # 内存使用超过 99.5%
    else
        return 0  # 内存使用低于 99.5%
    fi
}

# 监控任务是否完成
all_tasks_completed=false
while [ "$all_tasks_completed" == "false" ]; do
    sleep 30
    all_tasks_completed=true

    # 检查内存使用情况
    if ! check_memory_usage; then
        echo "Memory usage exceeded 99.5%. Stopping all tasks."
        kill 0  # 终止所有子进程
        exit 1
    fi

    # 检查任务完成情况
    for dir in *_iPhopResult; do
        if [ ! -f "$dir/Host_prediction_to_genus_m90.csv" ]; then
            echo "DRAM annotation still in progress in $dir."
            all_tasks_completed=false
            break
        fi
    done

    # 如果尚未完成，则再等待 30 秒
    if [ "$all_tasks_completed" == "false" ]; then
        sleep 30
    fi
done

echo "All tasks completed successfully."

# 合并所有注释结果到输出目录
awk 'FNR==1 && NR!=1{next;} {print}' ./*_iPhopResult/Host_prediction_to_genus_m90.csv > "$OUTPUT_DIR/Combined_Host_prediction_to_genus_m90.tsv"

echo "Host prediction complete. Results combined and saved to $OUTPUT_DIR/Combined_Host_prediction_to_genus_m90.tsv"

# 删除中间产生的临时文件
echo "Cleaning up temporary files..."
rm -rf "$OUTPUT_DIR/split_files"
rm -rf "$OUTPUT_DIR/iPhop_results"/*_iPhopResult

echo "Cleanup complete."
conda deactivate