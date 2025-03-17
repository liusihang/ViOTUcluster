#!/usr/bin/env python

import csv
import os
import sys
from Bio import SeqIO

def format_fasta(output_fasta):
    """Ensure the output FASTA file is formatted correctly."""
    with open(output_fasta, "r") as input_handle, open(output_fasta + ".tmp", "w") as output_handle:
        for record in SeqIO.parse(input_handle, "fasta"):
            output_handle.write(f">{record.id}\n{str(record.seq)}\n")
    os.replace(output_fasta + ".tmp", output_fasta)

# Accept command-line arguments
fasta = sys.argv[1]
Inputname = sys.argv[2]
OUT_DIR = sys.argv[3]

# Read the list of IDs from the CSV file
csv_values = set()
csv_path = os.path.join(OUT_DIR, f"{Inputname}_viral_predictionsList.csv")
with open(csv_path, 'r', encoding='utf-8') as csvfile:
    reader = csv.reader(csvfile)
    for row in reader:
        if row:  # Check if the row is not empty
            csv_values.add(row[0])  # Assume the first column contains the IDs

# Read the FASTA file and find matching sequences
matches = []
for record in SeqIO.parse(fasta, "fasta"):
    record_id = record.id.split()[0]  # Extract the first part of the ID before any spaces
    if record_id in csv_values:
        matches.append(record)

# If matching sequences are found, save them to a new FASTA file
if matches:
    output_fasta_filename = os.path.join(OUT_DIR, f"{Inputname}_filtered.fasta")
    SeqIO.write(matches, output_fasta_filename, "fasta")

    # Ensure the output FASTA file is formatted correctly
    format_fasta(output_fasta_filename)

print("All file were be processed.")
