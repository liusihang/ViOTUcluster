# -*- coding: utf-8 -*-
#!/usr/bin/env python3

import sys
import os
from Bio import SeqIO

def filter_sequences(min_length, input_dir, output_dir):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    for filename in os.listdir(input_dir):
        if filename.endswith(('.fasta', '.fa', '.fna')):
            input_file_path = os.path.join(input_dir, filename)
            output_file_path = os.path.join(output_dir, filename)

            with open(input_file_path, 'r') as input_handle, open(output_file_path, 'w') as output_handle:
                sequences = SeqIO.parse(input_handle, 'fasta')
                filtered_seqs = (seq for seq in sequences if len(seq.seq) >= min_length)
                SeqIO.write(filtered_seqs, output_handle, 'fasta')

    # Create a completion flag file
    with open(os.path.join(output_dir, 'Done'), 'w') as f:
        f.write('Filtering completed successfully.')

if __name__ == '__main__':
    import os
    import sys
    from Bio import SeqIO

    if len(sys.argv) != 4:
        print("Usage: python script.py <min_length> <input_dir> <output_dir>")
        sys.exit(1)

    min_len = int(sys.argv[1])
    input_dir = sys.argv[2]
    output_dir = sys.argv[3]

    done_file = os.path.join(output_dir, 'Done')

    if os.path.exists(done_file):
        print('Filtering already completed. Exiting.')
        sys.exit(0)
    else:
        filter_sequences(min_len, input_dir, output_dir)
        print('Filtering completed successfully.')