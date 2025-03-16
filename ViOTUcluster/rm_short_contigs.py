#!/usr/bin/env python3
import sys
from Bio import SeqIO

def filter_contigs(input_fasta, min_length):
    output_fasta = input_fasta 
    min_length = int(min_length)
    
    with open(input_fasta, "r") as input_handle, open(output_fasta, "w") as output_handle:
        for record in SeqIO.parse(input_handle, "fasta"):
            if len(record.seq) >= min_length:
                SeqIO.write(record, output_handle, "fasta")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python rm_short_contigs.py <input_fasta> <min_length>")
        sys.exit(1)
    
    input_fasta = sys.argv[1]
    min_length = sys.argv[2]
    
    filter_contigs(input_fasta, min_length)
    print(f"Filtered fasta file saved as {input_fasta}")