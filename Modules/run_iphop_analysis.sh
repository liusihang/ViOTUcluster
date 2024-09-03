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
mkdir -p "$OUTPUT_DIR/iphop_results"

echo -e "\n\n\n# 进行iPHOP分析!!!\n\n\n"
pwd

# 分割输入的fasta文件，每个新文件包含1000条序列
awk 'BEGIN {n_seq=0; file_idx=0;} 
     /^>/{
        if(n_seq % 1000 == 0 && n_seq != 0) {
            file_idx++;
        }
        file=sprintf("'"$OUTPUT_DIR"'/split_files/%s_%d.fna", FILENAME, file_idx);
        n_seq++;
     }
     {print >> file}' "$INPUT_FASTA"

# 切换到分割后的文件目录
cd "$OUTPUT_DIR/split_files" || { echo "Directory 'split_files' not found!"; exit 1; }

# 列出所有分割后的fna文件
ls *.fna > iPhop_list

# 使用 GNU Parallel 并行执行 iPHOP 注释
cat iPhop_list | parallel -j 2 iphop predict --fa_file {} --db_dir /media/em/student/db/Virus/Aug_2023_pub_rw --out_dir "$OUTPUT_DIR/iphop_results/{}" -t 10

# 等待所有并行任务完成
wait

# 合并所有注释结果到输出目录
awk 'FNR==1 && NR!=1 {next;} {print}' "$OUTPUT_DIR/iphop_results"/*/Host_prediction_to_genus_m90.csv > "$OUTPUT_DIR/Combined_Host_prediction_to_genus_m90.csv.tsv"

echo "Annotation complete. Results combined and saved to $OUTPUT_DIR/Combined_Host_prediction_to_genus_m90.csv.tsv"

# 删除中间产生的临时文件
echo "Cleaning up temporary files..."
rm -rf "$OUTPUT_DIR/split_files"
rm -rf "$OUTPUT_DIR/iphop_results"/*

echo "Cleanup complete."