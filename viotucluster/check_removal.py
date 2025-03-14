#!/usr/bin/env python
import sys
import pandas as pd
from Bio import SeqIO
import os

def filter_sequences(tsv_file, fasta_file):
    data = pd.read_csv(tsv_file, sep='\t')
    filter_condition = (data['viral_genes'] == 0) & ((data['host_genes'] > 0) | (data['provirus'] == 'Yes'))
    filtered_ids = set(data[filter_condition].iloc[:, 0])  # 假设第1列是序列ID

    sequences = SeqIO.parse(fasta_file, "fasta")
    filtered_sequences = [seq for seq in sequences if seq.id not in filtered_ids]
    temp_fasta = fasta_file + ".tmp"
    SeqIO.write(filtered_sequences, temp_fasta, "fasta")

    os.replace(temp_fasta, fasta_file)

    print(f"Filtered {len(filtered_ids)} sequences. Output saved to {fasta_file}.")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python filter_fasta.py <input_tsv> <input_fasta>")
    else:
        tsv_file = sys.argv[1]
        fasta_file = sys.argv[2]
        filter_sequences(tsv_file, fasta_file)