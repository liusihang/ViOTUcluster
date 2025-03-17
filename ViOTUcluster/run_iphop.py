#!/usr/bin/env python3

import os
import sys
import subprocess
import multiprocessing
from concurrent.futures import ThreadPoolExecutor, as_completed
import time
import signal
import logging

# Ensure necessary environment variables are set
required_env_vars = ['DATABASE', 'THREADS', 'OUTPUT_DIR','OUTPUT']
for var in required_env_vars:
    if var not in os.environ:
        print(f"Environment variable {var} is not set.")
        sys.exit(1)

# Get environment variables
THREADS = int(os.environ['THREADS'])
OUTPUT_DIR = os.environ['OUTPUT_DIR']
OUTPUT = os.environ['OUTPUT']
DATABASE = os.environ['DATABASE']

# Setup logging
log_dir = os.path.join(OUTPUT_DIR, 'Log')
if not os.path.exists(log_dir):
    os.makedirs(log_dir)
log_file = os.path.join(log_dir, 'iPhop.log')
logging.basicConfig(
    filename=log_file,
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)

# Get the number of cores available
all_cores = list(range(multiprocessing.cpu_count()))  # Get all core numbers of the system
CORES_TO_USE = THREADS 
assigned_cores = all_cores[:CORES_TO_USE]  # Assign cores to be used
# logging.info(f"Assigning tasks to cores: {assigned_cores}")

# Function to process a single iPhop prediction
def run_iphop_prediction(fa_file):
    try:
        output_dir = f"{fa_file}_iPhopResult"
        conda_prefix = os.environ.get("CONDA_PREFIX")
        env_path = os.path.join(conda_prefix, "envs", "iphop")
        t_value = min(THREADS, 20)  # Set t value: max 20, otherwise THREADS
        iphop_cmd = [
            'conda', 'run',
            '-p', env_path,
            'iphop', 'predict',
            '--fa_file', fa_file,
            '--db_dir', os.path.join(DATABASE, 'Aug_2023_pub_rw'),
            '--out_dir', output_dir,
            '-t', str(t_value)
        ]
        
        print(f"Running iPhop prediction for {fa_file}")
        process = subprocess.Popen(iphop_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        # Bind process to assigned cores if supported
        if hasattr(os, 'sched_setaffinity'):
            os.sched_setaffinity(process.pid, assigned_cores)

        stdout, stderr = process.communicate()

        if process.returncode != 0:
            logging.error(f"iPhop prediction failed for {fa_file} with exit code {process.returncode}")
            logging.error(stderr.decode())
            raise RuntimeError(f"iPhop prediction failed for {fa_file} with exit code {process.returncode}")

        print(f"iPhop prediction completed for {fa_file}")
        logging.info(stdout.decode())
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
    with open(os.path.join(OUTPUT, 'split_files', 'iPhop')) as f:
        files_list = [line.strip() for line in f if line.strip()]

    if not files_list:
        print("No files to process.")
        sys.exit(1)

    logging.info(f"Using {CORES_TO_USE} cores for iPhop prediction.")

    # Use ThreadPoolExecutor to process files in parallel
    with ThreadPoolExecutor(max_workers=min(THREADS, len(files_list))) as executor:
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
                logging.error(f"Task generated an exception: {e}")

    # Monitor the completion of iPhop predictions
    # monitor_iphop_tasks(files_list)

    print("All iPhop predictions have been processed.")

if __name__ == "__main__":
    main()
