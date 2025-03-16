#!/usr/bin/env python3
import sys
import os
from multiprocessing import cpu_count

def split_fasta(input_file, output_dir):
    num_cpus = cpu_count()
    print(f"Detected {num_cpus} CPUs.")


    with open(input_file, "r") as f:
        fasta_data = f.readlines()

    sequences_per_file = len(fasta_data) // (4 * num_cpus)
    print(f"Each file will contain {sequences_per_file} sequences.")

    os.makedirs(output_dir, exist_ok=True)

    current_file = 1
    current_sequences = 0
    current_output = None
    for line in fasta_data:
        if line.startswith(">") and current_sequences % sequences_per_file == 0:
            if current_output:
                current_output.close()
            current_output = open(os.path.join(output_dir, f"output_{current_file}.fasta"), "w")
            current_file += 1
        current_output.write(line)
        current_sequences += 1
    if current_output:
        current_output.close()

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script.py input.fasta output_directory")
        sys.exit(1)
    input_file = sys.argv[1]
    output_dir = sys.argv[2]
    split_fasta(input_file, output_dir)
