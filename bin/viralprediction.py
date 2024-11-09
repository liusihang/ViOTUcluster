#!/usr/bin/env python3

import os
import sys
import subprocess
import multiprocessing
from concurrent.futures import ThreadPoolExecutor, as_completed
import threading
import time
import glob
import signal

# Ensure necessary environment variables are set
required_env_vars = ['OUTPUT_DIR', 'DATABASE', 'Group', 'CONCENTRATION_TYPE', 'THREADS']
for var in required_env_vars:
    if var not in os.environ:
        print(f"Environment variable {var} is not set.")
        sys.exit(1)

# Get environment variables
OUTPUT_DIR = os.environ['OUTPUT_DIR']
DATABASE = os.environ['DATABASE']
Group = os.environ['Group']
CONCENTRATION_TYPE = os.environ['CONCENTRATION_TYPE']
THREADS = int(os.environ['THREADS'])  # Get the maximum number of threads

# Get file list from environment variable $FILES, or from default directory
FILES = os.environ.get('FILES')
if FILES:
    files_list = FILES.strip().split()
else:
    files_list = glob.glob(os.path.join(OUTPUT_DIR, 'FilteredSeqs', '*.fa')) + \
                 glob.glob(os.path.join(OUTPUT_DIR, 'FilteredSeqs', '*.fasta'))

if not files_list:
    print("No files to process.")
    sys.exit(1)

# Calculate the number of cores needed. 
CORES_TO_USE = THREADS 

# Get all available cores
all_cores = list(range(multiprocessing.cpu_count()))
# Get the first CORES_TO_USE cores
assigned_cores = all_cores[:CORES_TO_USE]
print(f"Assigning tasks to cores: {assigned_cores}")

# Function to process a single file
def process_file(file_path):
    basename = os.path.basename(file_path)
    if basename.endswith('.fasta'):
        basename = basename[:-6]  # 去掉 ".fasta"
    elif basename.endswith('.fa'):
        basename = basename[:-3]  # 去掉 ".fa"
    
    out_dir = os.path.join(OUTPUT_DIR, 'SeprateFile', basename)
    os.makedirs(out_dir, exist_ok=True)

    prediction_dir = os.path.join(out_dir, 'RoughViralPrediction')
    os.makedirs(prediction_dir, exist_ok=True)

    viralverify_dir = os.path.join(prediction_dir, 'viralverify')
    os.makedirs(viralverify_dir, exist_ok=True)

    virsorter_dir = os.path.join(prediction_dir, 'virsorter2')
    os.makedirs(virsorter_dir, exist_ok=True)

    genomad_dir = os.path.join(prediction_dir, 'genomadres')
    os.makedirs(genomad_dir, exist_ok=True)

    try:
        # Start ViralVerify task and bind to specified cores
        viralverify_result = os.path.join(viralverify_dir, f'{basename}_result_table.csv')
        if not os.path.isfile(viralverify_result):
            print(f"Starting Viralverify prediction for {file_path}.")
            process = subprocess.Popen(
                ['viralverify', '-f', file_path, '-o', viralverify_dir,
                '--hmm', os.path.join(DATABASE, 'ViralVerify', 'nbc_hmms.h3m'),
                '-t', str(THREADS)],
                stdout=subprocess.PIPE, stderr=subprocess.PIPE
            )
            
            # Bind process to specified cores if supported
            if hasattr(os, 'sched_setaffinity'):
                os.sched_setaffinity(process.pid, assigned_cores)

            process.wait()
            if process.returncode != 0:
                raise RuntimeError(f"Viralverify failed with exit code {process.returncode}")

            print(f"Viralverify prediction completed for {file_path}.")
        else:
            print(f"Viralverify prediction already completed for {file_path}, skipping...")
        
        # VirSorter2 task
        virsorter_result = os.path.join(virsorter_dir, 'final-viral-score.tsv')
        if CONCENTRATION_TYPE == 'non-concentration':
            print(f"Skipping Virsorter2 prediction for {file_path} due to non-concentration mode.")
        else:
            if not os.path.isfile(virsorter_result):
                print(f"Starting Virsorter2 prediction for {file_path}.")
                virsorter_cmd = [
                    'virsorter', 'run', '-w', virsorter_dir, '-i', file_path,
                    '--include-groups', Group, '-j', str(THREADS),
                    'all', '--min-score', '0.5', '--min-length', '300',
                    '--keep-original-seq', '-d', os.path.join(DATABASE, 'db')
                ]
                process = subprocess.Popen(virsorter_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

                # 将任务绑定到指定核心
                if hasattr(os, 'sched_setaffinity'):
                    os.sched_setaffinity(process.pid, assigned_cores)

                process.wait()
                if process.returncode != 0:
                    raise RuntimeError(f"Virsorter2 failed with exit code {process.returncode}")

                print(f"Virsorter2 prediction completed for {file_path}.")
            else:
                print(f"Virsorter2 prediction already completed for {file_path}, skipping...")

        # Genomad task
        genomad_result_dir = os.path.join(genomad_dir, f"{basename}_summary")
        os.makedirs(genomad_result_dir, exist_ok=True)

        if not os.path.isfile(os.path.join(genomad_result_dir, f"{basename}_virus_summary.tsv")):
            print(f"Starting Genomad prediction for {file_path}.")
            if CONCENTRATION_TYPE == 'non-concentration':
                genomad_cmd = [
                    'genomad', 'end-to-end', '--enable-score-calibration',
                    '--min-score', '0.7',
                    '--max-fdr', '0.05',
                    '--min-number-genes', '0',
                    '--min-plasmid-marker-enrichment', '0',
                    '--min-virus-marker-enrichment', '0',
                    '--min-plasmid-hallmarks', '1',
                    '--min-plasmid-hallmarks-short-seqs', '0',
                    '--min-virus-hallmarks', '0',
                    '--min-virus-hallmarks-short-seqs', '0',
                    '--max-uscg', '2',
                    file_path, genomad_dir, os.path.join(DATABASE, 'genomad_db'), '-t', str(THREADS)
                ]
            else:
                genomad_cmd = [
                    'genomad', 'end-to-end', '--enable-score-calibration',
                    '--min-score', '0.85',
                    '--max-fdr', '0.05',
                    '--min-number-genes', '1',
                    '--min-plasmid-marker-enrichment', '1.5',
                    '--min-virus-marker-enrichment', '1.5',
                    '--min-plasmid-hallmarks', '1',
                    '--min-plasmid-hallmarks-short-seqs', '1',
                    '--min-virus-hallmarks', '0',
                    '--min-virus-hallmarks-short-seqs', '1',
                    '--max-uscg', '2',
                    file_path, genomad_dir, os.path.join(DATABASE, 'genomad_db'), '-t', str(THREADS)
                ]
            process = subprocess.Popen(genomad_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

            # Bind Genomad task to specified cores
            if hasattr(os, 'sched_setaffinity'):
                os.sched_setaffinity(process.pid, assigned_cores)

            process.wait()
            if process.returncode != 0:
                raise RuntimeError(f"Genomad failed with exit code {process.returncode}")

            print(f"Genomad prediction completed for {file_path}.")
        else:
            print(f"{basename}: Genomad prediction already completed, skipping...")

        print(f"All predictions completed for {file_path}")

    except (subprocess.CalledProcessError, OSError) as e:
        print(f"An error occurred while processing {file_path}: {e}")
        raise

# Check if VirSorter2 and Genomad tasks are completed
def check_virsorter_completion():
    all_tasks_completed = False
    while not all_tasks_completed:
        all_tasks_completed = True
        for file_path in files_list:
            basename = os.path.basename(file_path)
            if basename.endswith('.fasta'):
                basename = basename[:-6]  #".fasta"
            elif basename.endswith('.fa'):
                basename = basename[:-3]  #".fa"
            virsorter_dir = os.path.join(
                OUTPUT_DIR, 'SeprateFile', basename, 'RoughViralPrediction', 'virsorter2'
            )

            virsorter_result = os.path.join(virsorter_dir, 'final-viral-score.tsv')
            if not os.path.isfile(virsorter_result):
                all_tasks_completed = False
                print("Virsorter2 still in processing")
                break

        if not all_tasks_completed:
            time.sleep(30)

# Main function
def main():
    # Handle termination signals
    def signal_handler(sig, frame):
        print("Process interrupted. Exiting gracefully...")
        sys.exit(0)

    signal.signal(signal.SIGINT, signal_handler)
    signal.signal(signal.SIGTERM, signal_handler)

    print(f"Total available threads: {THREADS}")
    print(f"Using {CORES_TO_USE} cores.")

    # Use ThreadPoolExecutor to start all tasks simultaneously
    with ThreadPoolExecutor(max_workers=len(files_list)) as executor:
        futures = []
        for file_path in files_list:
            # Submit task to thread pool
            future = executor.submit(process_file, file_path)
            futures.append(future)

        # Wait for all tasks to complete
        for future in as_completed(futures):
            try:
                future.result()
            except Exception as e:
                print(f"Task generated an exception: {e}")

    # Check if VirSorter2 and Genomad tasks are completed
    check_virsorter_completion()

    print("All files have been processed.")

if __name__ == "__main__":
    main()