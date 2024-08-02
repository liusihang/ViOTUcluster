#!/usr/bin/env python
import os
import glob
import sys

# 指定文件夹路径和输出文件名称
if len(sys.argv) < 3:
    print("Usage: MergedClusterSeqs.py <folder_path> <output_path>")
    sys.exit(1)
folder_path = sys.argv[1]
output_path = sys.argv[2]


def merge_and_rename_fasta_files(folder_path, output_path):
    # 初始化序列计数器
    sequence_counter = 1
    if not os.path.exists(output_path):
        os.makedirs(output_path,exist_ok=True)
    output_file = os.path.join(output_path, 'merged_sequences.fasta')
    # 创建或打开输出文件
    with open(output_file, 'w') as outfile:
        # 遍历文件夹中的所有.fasta文件
        for fasta_file in glob.glob(os.path.join(folder_path, '*.fasta')):
            # 打开并读取fasta文件
            with open(fasta_file, 'r') as infile:
                for line in infile:
                    # 如果行以'>'开头，它标志着一个序列的开始
                    if line.startswith('>'):
                        # 为序列生成一个新的名称
                        new_name = f'>vOTU{sequence_counter}\n'
                        outfile.write(new_name)
                        sequence_counter += 1
                    else:
                        # 直接写入序列数据
                        outfile.write(line)

# 调用函数
merge_and_rename_fasta_files(folder_path,output_path)
