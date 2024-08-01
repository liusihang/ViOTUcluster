import csv
import os
import sys
from Bio import SeqIO

# Accept command-line arguments
fasta = sys.argv[1]
Inputname = sys.argv[2]
OUT_DIR = sys.argv[3]

# Read the list of IDs from the CSV file
csv_values = set()
csv_path = os.path.join(OUT_DIR, f"{Inputname}.csv")
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

print("所有文件处理完成。")

