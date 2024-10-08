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

# 确保必要的环境变量已设置
required_env_vars = ['OUTPUT_DIR', 'DATABASE', 'Group', 'CONCENTRATION_TYPE', 'THREADS']
for var in required_env_vars:
    if var not in os.environ:
        print(f"Environment variable {var} is not set.")
        sys.exit(1)

# 获取环境变量
OUTPUT_DIR = os.environ['OUTPUT_DIR']
DATABASE = os.environ['DATABASE']
Group = os.environ['Group']
CONCENTRATION_TYPE = os.environ['CONCENTRATION_TYPE']
THREADS = int(os.environ['THREADS'])  # 获取最大线程数

# 从环境变量 $FILES 获取文件列表，或从默认目录获取
FILES = os.environ.get('FILES')
if FILES:
    files_list = FILES.strip().split()
else:
    # 假设文件位于默认目录；根据需要调整
    files_list = glob.glob(os.path.join(OUTPUT_DIR, 'FilteredSeqs', '*.fa')) + \
                 glob.glob(os.path.join(OUTPUT_DIR, 'FilteredSeqs', '*.fasta'))

if not files_list:
    print("No files to process.")
    sys.exit(1)

# 计算所需的核心数。例如，16 个线程可以使用 8 个核心（每个核心 2 个线程）
CORES_TO_USE = THREADS  # 例如 16 线程 -> 8 核心

# 获取可用的所有核心，假设系统总共有 16 个核心（您可以动态获取核心数）
all_cores = list(range(multiprocessing.cpu_count()))  # 获取系统所有核心编号

# 获取前 CORES_TO_USE 个核心
assigned_cores = all_cores[:CORES_TO_USE]
print(f"Assigning tasks to cores: {assigned_cores}")

# 处理单个文件的函数
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

    print(f"Processing {file_path}")

    try:
        # 启动外部任务并绑定到指定的核心
        process = subprocess.Popen(['viralverify', '-f', file_path, '-o', viralverify_dir,
                                    '--hmm', os.path.join(DATABASE, 'ViralVerify', 'nbc_hmms.hmm'),
                                    '-t', str(THREADS)])

        # 将进程绑定到指定核心
        os.sched_setaffinity(process.pid, assigned_cores)

        # 运行 Genomad 任务
        #logging.info(f"{basename}: Running Genomad end-to-end with up to {THREADS_PER_FILE} threads...")
        genomad_cmd = [
            'genomad', 'end-to-end', '--enable-score-calibration', file_path,
            genomad_dir, os.path.join(DATABASE, 'genomad_db'), '-t', str(THREADS)
        ]
        process = subprocess.Popen(genomad_cmd)

        # 绑定 Genomad 任务到指定核心
        os.sched_setaffinity(process.pid, assigned_cores)


        # VirSorter2 任务
        if CONCENTRATION_TYPE == "concentration":
            virsorter_result = os.path.join(virsorter_dir, 'final-viral-score.tsv')
            if not os.path.isfile(virsorter_result):
                print(f"{basename}: Running Virsorter2 with up to {THREADS} threads...")
                virsorter_cmd = [
                    'virsorter', 'run', '-w', virsorter_dir, '-i', file_path,
                    '--include-groups', Group, '-j', str(THREADS),
                    'all', '--min-score', '0.5', '--min-length', '300',
                    '--keep-original-seq', '-d', os.path.join(DATABASE, 'db')
                ]
                process = subprocess.Popen(virsorter_cmd)

                # 绑定 Virsorter2 任务到指定核心
                os.sched_setaffinity(process.pid, assigned_cores)
                #process.wait()

        print(f"All predictions completed for {file_path}")

    except Exception as e:
        print(f"An error occurred while processing {file_path}: {e}")
        raise
# 主函数
def main():
    # 处理终止信号
    def signal_handler(sig, frame):
        print("Process interrupted. Exiting gracefully...")
        sys.exit(0)

    signal.signal(signal.SIGINT, signal_handler)
    signal.signal(signal.SIGTERM, signal_handler)

    print(f"Total available threads: {THREADS}")
    print(f"Using {CORES_TO_USE} cores.")

    # 使用 ThreadPoolExecutor 同时启动所有任务
    with ThreadPoolExecutor(max_workers=len(files_list)) as executor:
        futures = []
        for file_path in files_list:
            # 提交任务到线程池
            future = executor.submit(process_file, file_path)
            futures.append(future)

        # 等待所有任务完成
        for future in as_completed(futures):
            try:
                future.result()
            except Exception as e:
                print(f"Task generated an exception: {e}")

if __name__ == "__main__":
    main()