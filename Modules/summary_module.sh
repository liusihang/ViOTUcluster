#!/usr/bin/env bash

# 创建目标目录用于存放所有组合后的fasta文件和重命名的quality summary文件
mkdir -p "$OUTPUT_DIR/Summary/unbined"
mkdir -p "$OUTPUT_DIR/Summary/bins"

# 遍历当前目录下的所有子目录
for folder in "$OUTPUT_DIR/SeprateFile/"*/ ; do
  if [ -d "$folder" ]; then
    # 获取文件夹的名字，去除尾部的斜线
    folderName=$(basename "$folder")
    echo "Processing folder: $folderName"

    # 处理unbined.fasta文件
    unbined_source="${folder}Binning/Summary/Finialfasta/unbined.fasta"
    if [ -f "$unbined_source" ]; then
      unbined_dest="$OUTPUT_DIR/Summary/unbined/${folderName}_unbined.fasta"
      cp "$unbined_source" "$unbined_dest"
      # 删除源文件
      rm "$unbined_source"
    else
      echo "Warning: $unbined_source does not exist."
    fi

    # 处理Bestbins文件夹中的fasta文件
    bestbins_source="${folder}Binning/Summary/Finialfasta/Bestbins/"
    if [ -d "$bestbins_source" ]; then
      for fasta_file in "$bestbins_source"*.fasta; do
        if [ -f "$fasta_file" ]; then
          fasta_filename=$(basename "$fasta_file")
          fasta_dest="$OUTPUT_DIR/Summary/bins/${folderName}_${fasta_filename}"
          cp "$fasta_file" "$fasta_dest"
          # 删除源文件
          rm "$fasta_file"
        else
          echo "Warning: No .fasta files found in $bestbins_source"
        fi
      done
      # 删除Bestbins目录，如果已空
      rmdir "$bestbins_source" 2>/dev/null || true
    else
      echo "Warning: $bestbins_source does not exist."
    fi

    # 删除Finialfasta目录，如果已空
    finialfasta_dir="${folder}Binning/Summary/Finialfasta/"
    rmdir "$finialfasta_dir" 2>/dev/null || true

    # 删除Summary目录，如果已空
    summary_dir="${folder}Binning/Summary/"
    rmdir "$summary_dir" 2>/dev/null || true

    # 删除Binning目录，如果已空
    binning_dir="${folder}Binning/"
    rmdir "$binning_dir" 2>/dev/null || true
  fi
done

# Clean up intermediate files
rm -r "$OUTPUT_DIR/Summary/temp"

echo "All files processed and cleaned up successfully."
