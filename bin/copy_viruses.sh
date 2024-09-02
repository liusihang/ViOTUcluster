#!/bin/bash

# 检查输入参数
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input_directory> <output_directory>"
    exit 1
fi

# 输入和输出目录
INPUT_DIR="$1"
OUTPUT_DIR="$2"

# 创建输出目录（如果不存在）
mkdir -p "$OUTPUT_DIR"

# 找到所有以 _viral_predictionsList_CheckRes 结尾的文件夹
find "$INPUT_DIR" -type d -name "*_viral_predictionsList_CheckRes" | while read -r dir; do
    # 提取文件夹前面的部分
    folder_name=$(basename "$dir")
    prefix_name=${folder_name%_viral_predictionsList_CheckRes}

    # 源文件和目标文件
    src_file="$dir/viruses.fna"
    dest_file="$OUTPUT_DIR/${prefix_name}.fna"

    # 复制并重命名文件
    if [ -f "$src_file" ]; then
        cp "$src_file" "$dest_file"
        echo "Copied $src_file to $dest_file"
    else
        echo "File $src_file does not exist."
    fi
done
