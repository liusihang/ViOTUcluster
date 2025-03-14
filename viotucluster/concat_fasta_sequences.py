#!/usr/bin/env python3
import os
import sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

def concatenate_fasta_sequences(input_dir, output_fasta):
    # List to store the concatenated sequences
    concatenated_sequences = []

    # Iterate through each fasta file in the input directory
    for filename in os.listdir(input_dir):
        if filename.endswith(".fa") or filename.endswith(".fasta"):
            file_path = os.path.join(input_dir, filename)
            sequence_name = os.path.splitext(filename)[0]  # Use the filename (without extension) as the sequence name
            sequences = []

            # Read all sequences in the fasta file
            for record in SeqIO.parse(file_path, "fasta"):
                sequences.append(str(record.seq))

            # Concatenate all sequences into one single sequence
            concatenated_sequence = "".join(sequences)
            new_record = SeqRecord(Seq(concatenated_sequence), id=sequence_name, description="")
            concatenated_sequences.append(new_record)

    # Write the concatenated sequences to the output fasta file
    with open(output_fasta, "w") as output_handle:
        SeqIO.write(concatenated_sequences, output_handle, "fasta")

if __name__ == "__main__":
    # Ensure the correct number of arguments are provided
    if len(sys.argv) != 3:
        print("Usage: python concat_fasta_sequences.py <input_directory> <output_fasta>")
        sys.exit(1)

    input_dir = sys.argv[1]
    output_fasta = sys.argv[2]

    # Call the function to concatenate sequences
    concatenate_fasta_sequences(input_dir, output_fasta)
    print(f"Concatenated sequences saved to {output_fasta}")