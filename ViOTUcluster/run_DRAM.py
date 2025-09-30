#!/usr/bin/env python3

import os
import sys
import subprocess
import multiprocessing
from concurrent.futures import ThreadPoolExecutor, as_completed
import time
import signal
import logging
from logging.handlers import RotatingFileHandler

# Ensure necessary environment variables are set
required_env_vars = ['THREADS', 'OUTPUT_DIR']
for var in required_env_vars:
    if var not in os.environ:
        print(f"Environment variable {var} is not set.")
        sys.exit(1)

# Get environment variables
THREADS = int(os.environ['THREADS'])
OUTPUT_DIR = os.environ['OUTPUT_DIR']
LOG_DIR = os.path.join(OUTPUT_DIR, 'Log')
os.makedirs(LOG_DIR, exist_ok=True)


def configure_logger() -> logging.Logger:
    logger = logging.getLogger("run_DRAM")
    if logger.handlers:
        return logger
    logger.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    log_file = os.path.join(LOG_DIR, 'run_DRAM.log')
    file_handler = RotatingFileHandler(log_file, maxBytes=5 * 1024 * 1024, backupCount=2)
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)

    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)
    return logger


LOGGER = configure_logger()

# Get the number of cores available
all_cores = list(range(multiprocessing.cpu_count()))  # Get all core numbers of the system
CORES_TO_USE = THREADS
assigned_cores = all_cores[:CORES_TO_USE]  # Assign cores to be used
#print(f"Assigning tasks to cores: {assigned_cores}")

LOGGER.info("run_DRAM initialized with THREADS=%s; log file stored under %s", THREADS, os.path.join(LOG_DIR, 'run_DRAM.log'))

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
        LOGGER.info("Running DRAM annotation for %s", fa_file)
        process = subprocess.Popen(dram_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        # Bind process to assigned cores if supported
        if hasattr(os, 'sched_setaffinity'):
            os.sched_setaffinity(process.pid, assigned_cores)

        stdout, stderr = process.communicate()

        if process.returncode != 0:
            stderr_text = stderr.decode(errors='ignore').strip()
            if stderr_text:
                LOGGER.error("DRAM stderr for %s:\n%s", fa_file, stderr_text)
            raise RuntimeError(f"DRAM annotation failed for {fa_file} with exit code {process.returncode}")
        stdout_text = stdout.decode(errors='ignore').strip()
        if stdout_text:
            LOGGER.debug("DRAM stdout for %s:\n%s", fa_file, stdout_text)
        LOGGER.info("DRAM annotation completed for %s", fa_file)
    except subprocess.CalledProcessError as e:
        LOGGER.error("Subprocess error while processing %s: %s", fa_file, e)
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
                LOGGER.info("DRAM annotation still in progress for %s", fa_file)
                break

        if not all_tasks_completed:
            time.sleep(60)

# Main function
def main():
    # Handle termination signals
    def signal_handler(sig, frame):
        LOGGER.warning("Process interrupted (signal %s). Exiting gracefully...", sig)
        sys.exit(0)

    signal.signal(signal.SIGINT, signal_handler)
    signal.signal(signal.SIGTERM, signal_handler)

    # Read list of split files from the "DRAM" file
    with open(os.path.join(OUTPUT_DIR, 'split_files', 'DRAM')) as f:
        files_list = [line.strip() for line in f if line.strip()]

    if not files_list:
        LOGGER.error("No files to process.")
        sys.exit(1)

    LOGGER.info("Using %s cores for DRAM annotation across %s files.", CORES_TO_USE, len(files_list))

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
                LOGGER.error("Task generated an exception: %s", e)
if __name__ == "__main__":
    main()
