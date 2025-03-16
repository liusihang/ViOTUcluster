#!/usr/bin/env python3

import os
import sys
import subprocess
import multiprocessing
from concurrent.futures import ThreadPoolExecutor, as_completed
import time
import glob
import signal

# Ensure necessary environment variables are set
required_env_vars = ['THREADS', 'OUTPUT_DIR']
for var in required_env_vars:
    if var not in os.environ:
        print(f"Environment variable {var} is not set.")
        sys.exit(1)

# Get environment variables
THREADS = int(os.environ['THREADS'])
OUTPUT_DIR = os.environ['OUTPUT_DIR']

# Get the number of cores available
all_cores = list(range(multiprocessing.cpu_count()))  # Get all core numbers of the system
CORES_TO_USE = THREADS
assigned_cores = all_cores[:CORES_TO_USE]  # Assign cores to be used
#print(f"Assigning tasks to cores: {assigned_cores}")

# Function to run DRAM annotation on a single file
def run_dram_annotation(fa_file):
    try:
        output_dir = f"{fa_file}_DRAMAnnot"
        conda_prefix = os.environ.get("CONDA_PREFIX")
        env_path = os.path.join(conda_prefix, "envs", "DRAM")
        dram_cmd = [
            'conda','run',
            '-p',env_path,
            'DRAM-v.py', 'annotate', 
            '-i', fa_file,
            '-o', output_dir,
            '--threads', str(THREADS)
        ]
        print(f"Running DRAM annotation for {fa_file}")
        process = subprocess.Popen(dram_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        # Bind process to assigned cores if supported
        if hasattr(os, 'sched_setaffinity'):
            os.sched_setaffinity(process.pid, assigned_cores)

        process.wait()

        if process.returncode != 0:
            raise RuntimeError(f"DRAM annotation failed for {fa_file} with exit code {process.returncode}")

        print(f"DRAM annotation completed for {fa_file}")
    except subprocess.CalledProcessError as e:
        print(f"An error occurred while processing {fa_file}: {e}")
        raise

# Function to monitor the completion of all DRAM tasks
def monitor_dram_tasks(files_list):
    all_tasks_completed = False
    while not all_tasks_completed:
        all_tasks_completed = True
        for fa_file in files_list:
            output_dir = f"{fa_file}_DRAMAnnot"
            result_file = os.path.join(output_dir, 'annotations.tsv')

            if not os.path.isfile(result_file):
                all_tasks_completed = False
                print(f"DRAM annotation still in progress for {fa_file}")
                break

        if not all_tasks_completed:
            time.sleep(60)

# Main function
def main():
    # Handle termination signals
    def signal_handler(sig, frame):
        print("Process interrupted. Exiting gracefully...")
        sys.exit(0)

    signal.signal(signal.SIGINT, signal_handler)
    signal.signal(signal.SIGTERM, signal_handler)

    # Read list of split files from the "DRAM" file
    with open(os.path.join(OUTPUT_DIR, 'split_files', 'DRAM')) as f:
        files_list = [line.strip() for line in f if line.strip()]

    if not files_list:
        print("No files to process.")
        sys.exit(1)

    print(f"Using {CORES_TO_USE} cores for DRAM annotation.")

    # Use ThreadPoolExecutor to process files in parallel
    with ThreadPoolExecutor(max_workers=min(THREADS, len(files_list))) as executor:
        futures = []
        for fa_file in files_list:
            # Submit task to thread pool
            future = executor.submit(run_dram_annotation, fa_file)
            futures.append(future)

        # Wait for all tasks to complete
        for future in as_completed(futures):
            try:
                future.result()
            except Exception as e:
                print(f"Task generated an exception: {e}")
if __name__ == "__main__":
    main()
