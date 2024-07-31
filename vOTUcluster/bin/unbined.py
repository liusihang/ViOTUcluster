import os
import argparse
from Bio import SeqIO

def format_fasta(output_fasta):
    """Ensure the output FASTA file is formatted correctly."""
    with open(output_fasta, "r") as input_handle, open(output_fasta + ".tmp", "w") as output_handle:
        for record in SeqIO.parse(input_handle, "fasta"):
            output_handle.write(f">{record.id}\n{str(record.seq)}\n")
    os.replace(output_fasta + ".tmp", output_fasta)

def extract_real_names_from_folder(input_folder):
    """提取文件夹中所有fasta文件中的真正序列名"""
    real_names = set()
    for file_name in os.listdir(input_folder):
        if file_name.endswith(".fasta"):
            file_path = os.path.join(input_folder, file_name)
            for record in SeqIO.parse(file_path, "fasta"):
                real_name = record.id.split("__")[1]  # 提取形如k127_45641的真正序列名
                real_names.add(real_name)
    return real_names

def filter_sequences(reference_fasta, real_names, output_file):
    """过滤并生成新的fasta文件"""
    with open(output_file, "w") as out_fasta:
        for record in SeqIO.parse(reference_fasta, "fasta"):
            if record.id not in real_names:
                SeqIO.write(record, out_fasta, "fasta")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Filter sequences from fasta files")
    parser.add_argument('-i', '--input_folder', required=True, help='Path to the input folder containing fasta files')
    parser.add_argument('-r', '--reference_fasta', required=True, help='Path to the reference fasta file containing sequence names')
    parser.add_argument('-o', '--output_file', required=True, help='Path to the output fasta file')
    
    args = parser.parse_args()
    
    real_names = extract_real_names_from_folder(args.input_folder)
    filter_sequences(args.reference_fasta, real_names, args.output_file)
    format_fasta(args.output_file)
    
    print(f"Filtered sequences have been saved to {args.output_file}")