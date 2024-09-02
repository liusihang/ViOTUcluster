#!/usr/bin/env python
import sys
import pandas as pd
from Bio import SeqIO
import os

def filter_sequences(tsv_file, fasta_file):
    # 读取TSV文件
    data = pd.read_csv(tsv_file, sep='\t')

    # 筛选出viral_genes为0且(host_genes大于0或者provirus为Yes)的行
    filter_condition = (data['viral_genes'] == 0) & ((data['host_genes'] > 0) | (data['provirus'] == 'Yes'))
    filtered_ids = set(data[filter_condition].iloc[:, 0])  # 假设第1列是序列ID

    # 读取FASTA文件并过滤序列
    sequences = SeqIO.parse(fasta_file, "fasta")
    filtered_sequences = [seq for seq in sequences if seq.id not in filtered_ids]

    # 将剩余的序列写入临时FASTA文件
    temp_fasta = fasta_file + ".tmp"
    SeqIO.write(filtered_sequences, temp_fasta, "fasta")

    # 用临时文件替换原始FASTA文件
    os.replace(temp_fasta, fasta_file)

    print(f"Filtered {len(filtered_ids)} sequences. Output saved to {fasta_file}.")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python filter_fasta.py <input_tsv> <input_fasta>")
    else:
        tsv_file = sys.argv[1]
        fasta_file = sys.argv[2]
        filter_sequences(tsv_file, fasta_file)