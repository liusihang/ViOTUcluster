#!/usr/bin/env python3
import collections
import os
import argparse
import pyhmmer
from pyhmmer.easel import SequenceFile
from pyhmmer.plan7 import HMMFile

def run_pyhmmer(faa_file, hmm_dbs, output_file, threads, evalue, bitscore):
    """
    Runs PyHMMER on given faa file with multiple databases
    :param faa_file: input .faa file path
    :param hmm_dbs: dictionary of HMM database file paths and their corresponding names
    :param output_file: output file path
    :param threads: number of threads
    :param evalue: evalue threshold for pyhmmer
    :param bitscore: bitscore threshold for pyhmmer
    :return: None
    """
    # Define result named tuple
    Result = collections.namedtuple("Result", ["protein", "db", "phrog", "bitscore", "evalue"])

    # Run hmmscan and get all results
    results = []
    for db_name, hmm_db_file in hmm_dbs.items():
        with pyhmmer.plan7.HMMFile(hmm_db_file) as hmms:  # Load HMMs
            with pyhmmer.easel.SequenceFile(faa_file, digital=True) as seqs:  # Load sequences
                for hits in pyhmmer.hmmer.hmmsearch(hmms, seqs, cpus=int(threads), E=float(evalue), T=bitscore):  # Run hmmscan
                    protein = hits.query_name.decode()  # Get protein from the hit
                    for hit in hits:
                        if hit.included:
                            # Include the hit to the result collection
                            results.append(Result(protein, db_name, hit.name.decode(), hit.score, hit.evalue))

    # Write results to output file
    with open(output_file, 'w') as out_f:
        #out_f.write("Query Sequence\tDatabase\tTarget HMM\tScore\tE-value\n")
        for result in results:
            out_f.write(f"{result.phrog}\t{result.db}\t{result.protein}\t{result.bitscore}\t{result.evalue}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run pyhmmer on given faa file with multiple databases")
    parser.add_argument('--faa_file', required=True, help='Input .faa file path')
    parser.add_argument('--db_dir', required=True, help='Database directory path')
    parser.add_argument('--output_file', required=True, help='Output file path')
    parser.add_argument('--threads', type=int, default=1, help='Number of threads to use')
    parser.add_argument('--evalue', type=float, default=1e-5, help='E-value threshold for pyhmmer')
    parser.add_argument('--bitscore', type=float, default=30, help='Bitscore threshold for pyhmmer')

    args = parser.parse_args()

    hmm_dbs = {
        "vir": f"{args.db_dir}/hmm/viral/combined.hmm",
        "arc": f"{args.db_dir}/hmm/pfam/Pfam-A-Archaea.hmm",
        "bac": f"{args.db_dir}/hmm/pfam/Pfam-A-Bacteria.hmm",
        "euk": f"{args.db_dir}/hmm/pfam/Pfam-A-Eukaryota.hmm",
        "mix": f"{args.db_dir}/hmm/pfam/Pfam-A-Mixed.hmm",
    }

    run_pyhmmer(args.faa_file, hmm_dbs, args.output_file, args.threads, args.evalue, args.bitscore)

#python run_pyhmmer.py --faa_file /path/to/your/file.faa --db_dir /path/to/your/db_dir --output_file /path/to/your/output.tsv --threads 100 --evalue 1e-5 --bitscore 30
