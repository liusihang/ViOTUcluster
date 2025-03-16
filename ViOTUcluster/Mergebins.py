#!/usr/bin/env python

import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import argparse

def combine_sequences_in_folder(input_folder):
    """Combine sequences in each fasta file into one sequence, and assign id bin_1, bin_2, etc."""
    combined_sequences = []
    bin_number = 1
    for file_name in sorted(os.listdir(input_folder)):
        if file_name.endswith(".fasta"):
            file_path = os.path.join(input_folder, file_name)
            sequences = []
            for record in SeqIO.parse(file_path, "fasta"):
                sequences.append(str(record.seq))
            if sequences:
                combined_seq = ''.join(sequences)
                new_record = SeqRecord(Seq(combined_seq), id=f'bin_{bin_number}', description='')
                combined_sequences.append(new_record)
                bin_number += 1
    return combined_sequences

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Combine sequences in fasta files into one sequence per file.")
    parser.add_argument('-i', '--input_folder', required=True, help='Path to the input folder containing fasta files')
    parser.add_argument('-o', '--output_file', required=True, help='Path to the output fasta file')
    
    args = parser.parse_args()
    
    output_dir = os.path.dirname(args.output_file)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # Step 1: Combine sequences in input folder
    combined_sequences = combine_sequences_in_folder(args.input_folder)
    
    # Step 2: Write combined sequences to output fasta
    with open(args.output_file, "w") as output_handle:
        SeqIO.write(combined_sequences, output_handle, "fasta")
    
    print(f"Combined sequences have been saved to {args.output_file}")

#python script.py -i /path/to/input_folder -o /path/to/output_file.fasta