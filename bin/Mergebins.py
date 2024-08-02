#!/usr/bin/env python
import os
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def concatenate_fasta_sequences(input_folder, output_file):
    bin_counter = 1
    concatenated_sequences = []

    for file_name in os.listdir(input_folder):
        if file_name.endswith(".fa") or file_name.endswith(".fasta"):
            file_path = os.path.join(input_folder, file_name)
            concatenated_seq = ""

            for record in SeqIO.parse(file_path, "fasta"):
                concatenated_seq += str(record.seq)

            concatenated_record = SeqRecord(
                Seq(concatenated_seq),
                id=f"bin_{bin_counter}",
                description=""
            )
            concatenated_sequences.append(concatenated_record)
            bin_counter += 1

    with open(output_file, "w") as output_handle:
        SeqIO.write(concatenated_sequences, output_handle, "fasta")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Concatenate fasta sequences in each file and rename them.")
    parser.add_argument('-i', '--input_folder', required=True, help='Path to the input folder containing fasta files')
    parser.add_argument('-o', '--output_file', required=True, help='Path to the output fasta file')

    args = parser.parse_args()
    
    concatenate_fasta_sequences(args.input_folder, args.output_file)

#python concatenate_fasta.py -i path/to/your/fasta/folder -o Allbins.fasta
