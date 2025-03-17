#!/usr/bin/env python
import os
import glob
import sys


def merge_and_rename_fasta_files(folder_path, output_path):
    """
    Merge all FASTA files in the specified folder into a single FASTA file with renamed sequence headers.

    Each sequence header is replaced with a unique identifier in the format '>vOTU<number>'.

    Parameters:
    - folder_path (str): Path to the folder containing input FASTA files.
    - output_path (str): Path to the directory where the merged FASTA file will be saved.
    """
    sequence_counter = 1
    os.makedirs(output_path, exist_ok=True)
    output_file = os.path.join(output_path, 'merged_sequences.fasta')

    with open(output_file, 'w') as outfile:
        for fasta_file in glob.glob(os.path.join(folder_path, '*.fasta')):
            with open(fasta_file, 'r') as infile:
                for line in infile:
                    if line.startswith('>'):
                        new_name = f'>vOTU{sequence_counter}\n'
                        outfile.write(new_name)
                        sequence_counter += 1
                    else:
                        outfile.write(line)


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: MergedClusterSeqs.py <folder_path> <output_path>")
        sys.exit(1)

    folder_path = sys.argv[1]
    output_path = sys.argv[2]

    merge_and_rename_fasta_files(folder_path, output_path)