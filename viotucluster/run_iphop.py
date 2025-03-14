#!/usr/bin/env python3

import os
import sys
import subprocess
import multiprocessing
from concurrent.futures import ThreadPoolExecutor, as_completed
import time
import signal

# Ensure necessary environment variables are set
required_env_vars = ['DATABASE', 'THREADS_PER_FILE', 'OUTPUT_DIR']
for var in required_env_vars:
    if var not in os.environ:
        print(f"Environment variable {var} is not set.")
        sys.exit(1)

# Get environment variables
THREADS_PER_FILE = int(os.environ['THREADS_PER_FILE'])
OUTPUT_DIR = os.environ['OUTPUT_DIR']
DATABASE = os.environ['DATABASE']

# Get the number of cores available
all_cores = list(range(multiprocessing.cpu_count()))  # Get all core numbers of the system
CORES_TO_USE = THREADS_PER_FILE 
assigned_cores = all_cores[:CORES_TO_USE]  # Assign cores to be used
print(f"Assigning tasks to cores: {assigned_cores}")

# Function to process a single iPhop prediction
def run_iphop_prediction(fa_file):
    try:
        output_dir = f"{fa_file}_iPhopResult"
        iphop_cmd = [
            'iphop', 'predict',
            '--fa_file', fa_file,
            '--db_dir', os.path.join(DATABASE, 'Aug_2023_pub_rw'),
            '--out_dir', output_dir,
            '-t', "20"
        ]
        print(f"Running iPhop prediction for {fa_file}")
        process = subprocess.Popen(iphop_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        # Bind process to assigned cores if supported
        if hasattr(os, 'sched_setaffinity'):
            os.sched_setaffinity(process.pid, assigned_cores)

        process.wait()

        if process.returncode != 0:
            raise RuntimeError(f"iPhop prediction failed for {fa_file} with exit code {process.returncode}")

        print(f"iPhop prediction completed for {fa_file}")
    except subprocess.CalledProcessError as e:
        print(f"An error occurred while processing {fa_file}: {e}")
        raise

# Function to monitor the completion of all iPhop tasks
def monitor_iphop_tasks(files_list):
    all_tasks_completed = False
    while not all_tasks_completed:
        all_tasks_completed = True
        for fa_file in files_list:
            output_dir = f"{fa_file}_iPhopResult"
            result_file = os.path.join(output_dir, 'Host_prediction_to_genus_m90.csv')

            if not os.path.isfile(result_file):
                all_tasks_completed = False
                print(f"iPhop prediction still in progress for {fa_file}")
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

    # Read list of split files from the "iPhop" file
    with open(os.path.join(OUTPUT_DIR, 'split_files', 'iPhop')) as f:
        files_list = [line.strip() for line in f if line.strip()]

    if not files_list:
        print("No files to process.")
        sys.exit(1)

    print(f"Using {CORES_TO_USE} cores for iPhop prediction.")

    # Use ThreadPoolExecutor to process files in parallel
    with ThreadPoolExecutor(max_workers=min(THREADS_PER_FILE, len(files_list))) as executor:
        futures = []
        for fa_file in files_list:
            # Submit task to thread pool
            future = executor.submit(run_iphop_prediction, fa_file)
            futures.append(future)

        # Wait for all tasks to complete
        for future in as_completed(futures):
            try:
                future.result()
            except Exception as e:
                print(f"Task generated an exception: {e}")

    # Monitor the completion of iPhop predictions
    #monitor_iphop_tasks(files_list)

    print("All iPhop predictions have been processed.")

if __name__ == "__main__":
    main()
