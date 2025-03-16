#!/usr/bin/env python3
# usage: 
# bwa mem -a assembly.fa reads_1.fastq reads_2.fastq | ./filter_reads_for_bin_reassembly.py original_bin_folder output_dir strict_snp_cutoff permissive_snp_cutoff
from __future__ import print_function
import sys
import os

strict_snp_cutoff = int(sys.argv[3])
permissive_snp_cutoff = int(sys.argv[4])

complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'a': 't', 't': 'a', 'c': 'g', 'g': 'c', 'N': 'N', 'n': 'n'}

def rev_comp(seq):
    rev_comp = ""
    for n in seq:
        rev_comp += complement[n]
    return rev_comp[::-1]

# Load bin contigs
print("Loading contig to bin mappings...")
contig_bins = {}
for bin_file in os.listdir(sys.argv[1]):
    if bin_file.endswith(".fa") or bin_file.endswith(".fasta"): 
        bin_name = ".".join(bin_file.split("/")[-1].split(".")[:-1])
        with open(os.path.join(sys.argv[1], bin_file)) as f:
            for line in f:
                if line[0] != ">": 
                    continue
                contig_bins[line[1:-1]] = bin_name

# Store the read names and what bins they belong in in these dictionaries
# strict stores only perfectly aligning reads and permissive stores any aligned reads

print("Parsing sam file and writing reads to appropriate files depending on what bin they aligned to...")
files = {}
opened_bins = {}
F_line = None

for line in sys.stdin:
    if line[0] == "@": 
        continue
    cut = line.strip().split("\t")
    binary_flag = bin(int(cut[1]))

    if binary_flag[-7] == "1":
        F_line = line
        continue
    elif binary_flag[-8] == "1":
        R_line = line

        # Get fields for forward and reverse reads
        F_cut = F_line.strip().split("\t")
        R_cut = R_line.strip().split("\t")

        # Skip non-aligned reads
        if F_cut[2] == "*" and R_cut[2] == "*": 
            continue

        # Make sure the R and F reads aligned to the same bin
        if F_cut[2] != R_cut[2]:
            if F_cut[2] not in contig_bins or R_cut[2] not in contig_bins: 
                continue
            bin1 = contig_bins[F_cut[2]]
            bin2 = contig_bins[R_cut[2]]
            if bin1 != bin2: 
                continue
            bin_name = bin1
        else:
            contig = F_cut[2]
            if contig not in contig_bins: 
                continue
            bin_name = contig_bins[contig]

        # Make sure the reads aligned again
        if "NM:i:" not in F_line and "NM:i:" not in R_line: 
            continue
        
        # Open the relevant output files
        if bin_name not in opened_bins:
            opened_bins[bin_name] = None
            files[os.path.join(sys.argv[2], bin_name + ".strict_1.fastq")] = open(os.path.join(sys.argv[2], bin_name + ".strict_1.fastq"), "w")
            files[os.path.join(sys.argv[2], bin_name + ".strict_2.fastq")] = open(os.path.join(sys.argv[2], bin_name + ".strict_2.fastq"), "w")
            files[os.path.join(sys.argv[2], bin_name + ".permissive_1.fastq")] = open(os.path.join(sys.argv[2], bin_name + ".permissive_1.fastq"), "w")
            files[os.path.join(sys.argv[2], bin_name + ".permissive_2.fastq")] = open(os.path.join(sys.argv[2], bin_name + ".permissive_2.fastq"), "w")

        # Count how many mismatches there are between the two reads
        cumulative_mismatches = 0
        for field in F_cut:
            if field.startswith("NM:i:"):
                cumulative_mismatches += int(field.split(":")[-1])
                break
        for field in R_cut:
            if field.startswith("NM:i:"):
                cumulative_mismatches += int(field.split(":")[-1])
                break

        # Determine alignment type from bitwise FLAG
        F_binary_flag = bin(int(F_cut[1]))
        R_binary_flag = bin(int(R_cut[1]))

        # If the reads are reversed, fix them
        if F_binary_flag[-5] == "1":
            F_cut[9] = rev_comp(F_cut[9])
            F_cut[10] = F_cut[10][::-1]
        if R_binary_flag[-5] == "1":
            R_cut[9] = rev_comp(R_cut[9])
            R_cut[10] = R_cut[10][::-1]

        # Strict assembly
        if cumulative_mismatches < strict_snp_cutoff:
            files[os.path.join(sys.argv[2], bin_name + ".strict_1.fastq")].write('@' + F_cut[0] + "/1" + "\n" + F_cut[9] + "\n+\n" + F_cut[10] + "\n")
            files[os.path.join(sys.argv[2], bin_name + ".strict_2.fastq")].write('@' + R_cut[0] + "/2" + "\n" + R_cut[9] + "\n+\n" + R_cut[10] + "\n")

        # Permissive assembly
        if cumulative_mismatches < permissive_snp_cutoff:
            files[os.path.join(sys.argv[2], bin_name + ".permissive_1.fastq")].write('@' + F_cut[0] + "/1" + "\n" + F_cut[9] + "\n+\n" + F_cut[10] + "\n")
            files[os.path.join(sys.argv[2], bin_name + ".permissive_2.fastq")].write('@' + R_cut[0] + "/2" + "\n" + R_cut[9] + "\n+\n" + R_cut[10] + "\n")

print("Closing files")
for f in files.values():
    f.close()

print("Finished splitting reads!")
