# -*- coding: utf-8 -*-
#!/usr/bin/env python

import sys
import os
from Bio import SeqIO

def filter_sequences(min_length, input_dir, output_dir):
    # 确保输出目录存在
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # 遍历输入目录中的所有文件
    for filename in os.listdir(input_dir):
        if filename.endswith(('.fasta', '.fna', '.fa')):
            input_file_path = os.path.join(input_dir, filename)
            output_file_path = os.path.join(output_dir, filename)
            
            # 过滤序列
            with open(input_file_path, 'r') as input_handle, open(output_file_path, 'w') as output_handle:
                sequences = SeqIO.parse(input_handle, 'fasta')
                filtered = (record for record in sequences if len(record.seq) >= min_length)
                SeqIO.write(filtered, output_handle, 'fasta')

if __name__ == '__main__':
    # 检查命令行参数
    if len(sys.argv) != 4:
        print("Usage: python filter_fasta.py <min_length> <input_dir> <output_dir>")
        sys.exit(1)

    min_len = int(sys.argv[1])
    in_dir = sys.argv[2]
    out_dir = sys.argv[3]
    
    filter_sequences(min_len, in_dir, out_dir)
