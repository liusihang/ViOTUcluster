#!/usr/bin/env python

import argparse
from Bio import SeqIO

def rename_sequences(input_file):
    count = 1
    records = []
    for record in SeqIO.parse(input_file, "fasta"):
        record.id = f"vOTU{count}"
        record.description = ""
        records.append(record)
        count += 1

    # Write the renamed sequences back to the same input file
    with open(input_file, 'w') as output_handle:
        SeqIO.write(records, output_handle, "fasta")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Rename sequences in a FASTA file")
    parser.add_argument('-i', '--input_file', required=True, help='Path to the input FASTA file')

    args = parser.parse_args()

    rename_sequences(args.input_file)

    print(f"Renamed sequences have been saved to {args.input_file}")
