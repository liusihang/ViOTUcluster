#!/usr/bin/env python3
import os
import sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

def concatenate_fasta_sequences(input_dir, output_fasta):
    concatenated_sequences = []

    # 遍历输入文件夹中的每个fasta文件
    for filename in os.listdir(input_dir):
        if filename.endswith(".fa") or filename.endswith(".fasta"):
            file_path = os.path.join(input_dir, filename)
            sequence_name = os.path.splitext(filename)[0]  # 序列名称是fasta文件的名称
            sequences = []

            # 读取fasta文件中的所有序列
            for record in SeqIO.parse(file_path, "fasta"):
                sequences.append(str(record.seq))

            # 将所有序列连接成一个新序列
            concatenated_sequence = "".join(sequences)
            new_record = SeqRecord(Seq(concatenated_sequence), id=sequence_name, description="")
            concatenated_sequences.append(new_record)

    # 将所有新序列写入输出fasta文件
    with open(output_fasta, "w") as output_handle:
        SeqIO.write(concatenated_sequences, output_handle, "fasta")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python concat_fasta_sequences.py <input_directory> <output_fasta>")
        sys.exit(1)

    input_dir = sys.argv[1]
    output_fasta = sys.argv[2]

    concatenate_fasta_sequences(input_dir, output_fasta)
    print(f"Concatenated sequences saved to {output_fasta}")