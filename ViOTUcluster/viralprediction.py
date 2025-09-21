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

# ===== 1) Check =====
required_env_vars = ['OUTPUT_DIR', 'DATABASE', 'Group', 'CONCENTRATION_TYPE', 'THREADS']
for var in required_env_vars:
    if var not in os.environ:
        print(f"Environment variable {var} is not set.")
        sys.exit(1)

# ===== 2) Read Envs =====
OUTPUT_DIR = os.environ['OUTPUT_DIR']
DATABASE = os.environ['DATABASE']
Group = os.environ['Group']
CONCENTRATION_TYPE = os.environ['CONCENTRATION_TYPE']
THREADS = int(os.environ['THREADS'])  

# 全局并发上限
try:
    MAX_TASKS = int(os.environ.get('MAX_TASKS', max(1, multiprocessing.cpu_count())))
except ValueError:
    MAX_TASKS = max(1, multiprocessing.cpu_count())

# 全局并发闸
_TASK_SEM = threading.BoundedSemaphore(value=max(1, MAX_TASKS))

# ===== 3) File list =====
FILES = os.environ.get('FILES')
if FILES:
    files_list = FILES.strip().split()
else:
    files_list = glob.glob(os.path.join(OUTPUT_DIR, 'FilteredSeqs', '*.fa')) + \
                 glob.glob(os.path.join(OUTPUT_DIR, 'FilteredSeqs', '*.fasta'))

if not files_list:
    print("No files to process.")
    sys.exit(1)

# ===== 4) Multi-processing setup =====
CORES_TO_USE = THREADS
all_cores = list(range(multiprocessing.cpu_count()))
assigned_cores = all_cores[:CORES_TO_USE]
print(f"Assigning tasks to cores: {assigned_cores}")
print(f"Global max concurrent tasks (MAX_TASKS): {MAX_TASKS}")
print(f"Per-task threads passed to tools (THREADS): {THREADS}")

def run_command(cmd, cores=None):
    """Run an external command and bind it to specified cores (if supported),
       with a global concurrency gate (_TASK_SEM)."""
    _TASK_SEM.acquire()  # 全局并发闸
    try:
        process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        # 子进程亲和性
        if hasattr(os, 'sched_setaffinity') and cores:
            try:
                os.sched_setaffinity(process.pid, cores)
            except Exception as e:
                print(f"[Warning] Failed to set CPU affinity for {' '.join(cmd)}: {e}")
        stdout, stderr = process.communicate()
        if process.returncode != 0:
            raise RuntimeError(
                f"Command {' '.join(cmd)} failed with exit code {process.returncode}.\nError: {stderr.decode(errors='ignore')}"
            )
        return stdout.decode(errors='ignore')
    finally:
        _TASK_SEM.release()  # 释放并发令牌

def process_file(file_path):
    basename = os.path.basename(file_path).replace('.fasta', '').replace('.fa', '')
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

    tasks = []
    inner_workers = max(1, min(3, MAX_TASKS))
    with ThreadPoolExecutor(max_workers=inner_workers) as executor:
        # ViralVerify
        viralverify_result = os.path.join(viralverify_dir, f'{basename}_result_table.csv')
        if not os.path.isfile(viralverify_result):
            viralverify_cmd = [
                'viralverify', '-f', file_path, '-o', viralverify_dir,
                '--hmm', os.path.join(DATABASE, 'ViralVerify', 'nbc_hmms.h3m'),
                '-t', str(THREADS)
            ]
            tasks.append(executor.submit(run_command, viralverify_cmd, assigned_cores))
        else:
            print(f"Viralverify prediction already completed for {file_path}, skipping...")

        # VirSorter2
        virsorter_result = os.path.join(virsorter_dir, 'final-viral-score.tsv')
        if not os.path.isfile(virsorter_result):
            virsorter_cmd = [
                'virsorter', 'run', '-w', virsorter_dir, '-i', file_path,
                '--include-groups', Group, '-j', str(THREADS),
                'all', '--min-score', '0.5', '--min-length', '300',
                '--keep-original-seq', '-d', os.path.join(DATABASE, 'db')
            ]
            tasks.append(executor.submit(run_command, virsorter_cmd, assigned_cores))
        else:
            print(f"Virsorter2 prediction already completed for {file_path}, skipping...")

        # Genomad
        genomad_result_dir = os.path.join(genomad_dir, f"{basename}_summary")
        os.makedirs(genomad_result_dir, exist_ok=True)
        genomad_result = os.path.join(genomad_result_dir, f"{basename}_virus_summary.tsv")

        if not os.path.isfile(genomad_result):
            if CONCENTRATION_TYPE == "concentration":
                genomad_cmd = [
                    'genomad', 'end-to-end', '--enable-score-calibration',
                    file_path, genomad_dir, os.path.join(DATABASE, 'genomad_db'),
                    '-t', str(THREADS),
                    '--min-score', '0.7', '--max-fdr', '0.05',
                    '--min-number-genes', '0',
                    '--min-virus-marker-enrichment', '1.5',
                    '--min-plasmid-marker-enrichment', '0',
                    '--min-plasmid-hallmarks', '1',
                    '--min-plasmid-hallmarks-short-seqs', '0',
                    '--max-uscg', '2'
                ]
            else:
                genomad_cmd = [
                    'genomad', 'end-to-end', '--enable-score-calibration',
                    file_path, genomad_dir, os.path.join(DATABASE, 'genomad_db'),
                    '-t', str(THREADS),
                    '--min-score', '0.8', '--max-fdr', '0.05',
                    '--min-number-genes', '1',
                    '--min-virus-marker-enrichment', '0',
                    '--min-plasmid-marker-enrichment', '1.5',
                    '--min-plasmid-hallmarks', '1',
                    '--min-plasmid-hallmarks-short-seqs', '1',
                    '--max-uscg', '2'
                ]
            tasks.append(executor.submit(run_command, genomad_cmd, assigned_cores))
        else:
            print(f"Genomad prediction already completed for {file_path}, skipping...")

        #Wait for all tasks to complete and handle exceptions
        for future in as_completed(tasks):
            try:
                future.result()
            except Exception as e:
                print(f"An error occurred while processing {file_path}: {e}")

    print(f"All predictions completed for {file_path}")

def check_virsorter_completion():
    all_tasks_completed = False
    while not all_tasks_completed:
        all_tasks_completed = True
        for file_path in files_list:
            basename = os.path.basename(file_path)
            if basename.endswith('.fasta'):
                basename = basename[:-6]
            elif basename.endswith('.fa'):
                basename = basename[:-3]
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

def main():
    def signal_handler(sig, frame):
        print("Process interrupted. Exiting gracefully...")
        sys.exit(0)
    signal.signal(signal.SIGINT, signal_handler)
    signal.signal(signal.SIGTERM, signal_handler)

    print(f"Total available threads (logical CPUs): {multiprocessing.cpu_count()}")
    print(f"Using per-task threads (THREADS): {THREADS}")
    print(f"Global concurrency cap (MAX_TASKS): {MAX_TASKS}")
    print(f"Total files to process: {len(files_list)}")

    # 外层线程池：以文件为粒度的提交并发，限到 MAX_TASKS
    outer_workers = max(1, min(len(files_list), MAX_TASKS))
    with ThreadPoolExecutor(max_workers=outer_workers) as executor:
        futures = [executor.submit(process_file, fp) for fp in files_list]
        for future in as_completed(futures):
            try:
                future.result()
            except Exception as e:
                print(f"Task generated an exception: {e}")

    # check_virsorter_completion()
    # print("All files have been processed.")

if __name__ == "__main__":
    main()
