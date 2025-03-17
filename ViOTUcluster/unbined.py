#!/usr/bin/env python3

import os
import argparse
from Bio import SeqIO

def format_fasta(output_fasta):
    """
    Ensure the output FASTA file is formatted correctly by rewriting it.

    Parameters:
    - output_fasta (str): Path to the output FASTA file to be formatted.
    """
    with open(output_fasta, "r") as input_handle, open(output_fasta + ".tmp", "w") as output_handle:
        for record in SeqIO.parse(input_handle, "fasta"):
            output_handle.write(f">{record.id}\n{str(record.seq)}\n")
    os.replace(output_fasta + ".tmp", output_fasta)

def extract_real_names_from_folder(input_folder):
    """
    Extract the real sequence names from all FASTA files in the input folder.

    Assumes that sequence IDs contain a '__' separator and the real name is the second part,
    e.g., 'k127__45641'.

    Parameters:
    - input_folder (str): Path to the folder containing input FASTA files.

    Returns:
    - set: A set of real sequence names.
    """
    real_names = set()
    for file_name in os.listdir(input_folder):
        if file_name.endswith(".fasta"):
            file_path = os.path.join(input_folder, file_name)
            for record in SeqIO.parse(file_path, "fasta"):
                try:
                    real_name = record.id.split("__")[1]  # Extract the real sequence name, e.g., k127_45641
                    real_names.add(real_name)
                except IndexError:
                    print(f"Warning: Sequence ID '{record.id}' does not contain '__'. Skipping...")
    return real_names

def filter_sequences(reference_fasta, real_names, output_file):
    """
    Filter out sequences from the reference FASTA file that are present in the real_names set
    and write the remaining sequences to a new FASTA file.

    Parameters:
    - reference_fasta (str): Path to the reference FASTA file containing sequence names.
    - real_names (set): Set of real sequence names to be excluded.
    - output_file (str): Path to the output FASTA file with filtered sequences.
    """
    with open(output_file, "w") as out_fasta:
        for record in SeqIO.parse(reference_fasta, "fasta"):
            if record.id not in real_names:
                SeqIO.write(record, out_fasta, "fasta")

def main():
    """
    Main function to execute the sequence filtering and formatting process.

    Expects three command-line arguments:
    1. Input folder containing FASTA files.
    2. Reference FASTA file containing sequence names.
    3. Output FASTA file path for the filtered sequences.
    """
    parser = argparse.ArgumentParser(description="Filter sequences from FASTA files based on real sequence names.")
    parser.add_argument('-i', '--input_folder', required=True, help='Path to the input folder containing FASTA files')
    parser.add_argument('-r', '--reference_fasta', required=True, help='Path to the reference FASTA file containing sequence names')
    parser.add_argument('-o', '--output_file', required=True, help='Path to the output FASTA file')
    
    args = parser.parse_args()
    
    # Extract real sequence names from the input folder
    real_names = extract_real_names_from_folder(args.input_folder)
    
    # Filter sequences from the reference FASTA based on the extracted real names
    filter_sequences(args.reference_fasta, real_names, args.output_file)
    
    # Ensure the output FASTA file is correctly formatted
    format_fasta(args.output_file)
    
    print(f"Filtered sequences have been saved to {args.output_file}")

if __name__ == "__main__":
    main()