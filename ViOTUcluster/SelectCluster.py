#!/usr/bin/env python3
import csv
import sys


fasta_file = sys.argv[1]
csv_file = sys.argv[2]
output_file = sys.argv[3]

csv.field_size_limit(1000000000)

def extract_sequences(fasta_file, csv_file, output_file):
    with open(csv_file, 'r') as csvf:
        csv_reader = csv.reader(csvf, delimiter= "\t")
        sequence_ids = [row[0] for row in csv_reader]

    with open(fasta_file, 'r') as fasta, open(output_file, 'w') as out:
        write_sequence = False
        for line in fasta:
            if line.startswith('>'):
                sequence_id = line[1:].strip().split()[0] 
                if sequence_id in sequence_ids:
                    write_sequence = True
                    out.write(line)
                else:
                    write_sequence = False
            else:
                if write_sequence:
                    out.write(line)

extract_sequences(fasta_file, csv_file, output_file)
