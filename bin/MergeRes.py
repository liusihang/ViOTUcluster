#!/usr/bin/env python

import os
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def format_fasta(output_fasta):
    """Ensure the output FASTA file is formatted correctly."""
    with open(output_fasta, "r") as input_handle, open(output_fasta + ".tmp", "w") as output_handle:
        for record in SeqIO.parse(input_handle, "fasta"):
            output_handle.write(f">{record.id}\n{str(record.seq)}\n")
    os.replace(output_fasta + ".tmp", output_fasta)

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

def extract_real_names_from_folder(input_folder):
    """Extract real sequence names from all fasta files in the folder."""
    real_names = set()
    for file_name in os.listdir(input_folder):
        if file_name.endswith(".fasta"):
            file_path = os.path.join(input_folder, file_name)
            for record in SeqIO.parse(file_path, "fasta"):
                real_name = record.id.split("__")[1]  # Extract the real sequence name like k127_45641
                real_names.add(real_name)
    return real_names

def filter_sequences(reference_fasta, real_names, output_file):
    """Filter sequences and generate a new fasta file."""
    filtered_records = []
    for record in SeqIO.parse(reference_fasta, "fasta"):
        if record.id not in real_names:
            filtered_records.append(record)
    SeqIO.write(filtered_records, output_file, "fasta")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Filter sequences from fasta files and combine sequences")
    parser.add_argument('-i', '--input_folder', required=True, help='Path to the input folder containing fasta files')
    parser.add_argument('-r', '--reference_fasta', required=True, help='Path to the reference fasta file containing sequence names')
    parser.add_argument('-o', '--output_file', required=True, help='Path to the output fasta file')
    
    args = parser.parse_args()
    
    # Step 1: Combine sequences in input folder
    combined_sequences = combine_sequences_in_folder(args.input_folder)
    
    # Step 2: Extract real names from input folder
    real_names = extract_real_names_from_folder(args.input_folder)
    
    # Step 3: Filter sequences from reference fasta
    filtered_sequences_file = args.output_file + ".filtered.tmp"
    filter_sequences(args.reference_fasta, real_names, filtered_sequences_file)
    
    # Step 4: Combine the combined_sequences and the filtered sequences into the final output
    final_sequences = combined_sequences + list(SeqIO.parse(filtered_sequences_file, "fasta"))
    SeqIO.write(final_sequences, args.output_file, "fasta")
    
    # Format the output fasta
    format_fasta(args.output_file)
    
    # Clean up temporary file
    os.remove(filtered_sequences_file)
    
    print(f"Combined sequences have been saved to {args.output_file}")