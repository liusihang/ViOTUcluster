#!/usr/bin/env python3
import os
import glob
import shutil
import argparse
import subprocess
import multiprocessing
from multiprocessing import Pool

def find_fastq_pairs(input_dir):
    """
    Recursively search for all *_R1* FASTQ files in the input directory
    and find their corresponding *_R2* files.
    """
    patterns = ["*_R1.fq.gz", "*_R1.fastq.gz", "*_R1.fq", "*_R1.fastq"]
    r1_files = []
    for pattern in patterns:
        r1_files.extend(glob.glob(os.path.join(input_dir, "**", pattern), recursive=True))
    r1_files = sorted(r1_files)
    pairs = []
    for r1 in r1_files:
        name = os.path.basename(r1)
        idx = name.rfind("_R1")
        if idx == -1:
            continue  # Skip files that do not contain '_R1'
        # Construct the corresponding R2 filename
        r2_name = name[:idx] + "_R2" + name[idx+3:]
        r2_path = os.path.join(os.path.dirname(r1), r2_name)
        if os.path.exists(r2_path):
            sample = name[:idx]  # Extract sample name (portion before _R1)
            pairs.append((sample, r1, r2_path))
        else:
            print(f"[Warning] Found R1 file but no matching R2: {r1}")
    return pairs

def run_fastp(r1_path, r2_path, out_r1, out_r2):
    """
    Run fastp to clean the given R1/R2 FASTQ files.
    """
    try:
        cmd = ["fastp", "-i", r1_path, "-I", r2_path, "-o", out_r1, "-O", out_r2]
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            raise subprocess.CalledProcessError(result.returncode, cmd, result.stdout, result.stderr)
        print(f"[fastp] Completed cleaning for {os.path.basename(r1_path)}")
        return True
    except subprocess.CalledProcessError as e:
        print(f"[Error] fastp failed for {os.path.basename(r1_path)}: {e.stderr or e.stdout}")
        return False
    except Exception as ex:
        print(f"[Error] Unexpected exception in fastp for {os.path.basename(r1_path)}: {ex}")
        return False

def run_assembler(clean_r1, clean_r2, assembler, out_dir, threads):
    """
    Run the specified assembler (megahit or metaspades) on the cleaned reads.
    """
    try:
        if assembler.lower() == "megahit":
            cmd = ["megahit", "-1", clean_r1, "-2", clean_r2, "-o", out_dir, "-t", str(threads)]
        elif assembler.lower() == "metaspades":
            cmd = ["metaspades.py", "-1", clean_r1, "-2", clean_r2, "-o", out_dir, "-t", str(threads)]
        else:
            raise ValueError(f"Unsupported assembler: {assembler}")
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            raise subprocess.CalledProcessError(result.returncode, cmd, result.stdout, result.stderr)
        print(f"[{assembler}] Assembly completed for reads: {os.path.basename(clean_r1)}")
        return True
    except subprocess.CalledProcessError as e:
        print(f"[Error] {assembler} failed: {e.stderr or e.stdout}")
        return False
    except Exception as ex:
        print(f"[Error] Unexpected exception running {assembler}: {ex}")
        return False

def main():
    parser = argparse.ArgumentParser(description="Parallel fastp cleaning and assembly pipeline.")
    parser.add_argument("-i", "--input_dir", required=True, help="Directory containing raw FASTQ files")
    parser.add_argument("-o", "--output_dir", required=True, help="Output directory for cleaned reads and assembly results")
    parser.add_argument("-c", "--cores", type=int, default=multiprocessing.cpu_count(),
                        help="Number of CPU cores to use in parallel")
    parser.add_argument("-a", "--assembler", choices=["megahit", "metaspades"], required=True,
                        help="Assembler to use: megahit or metaspades")
    args = parser.parse_args()

    input_dir = args.input_dir
    output_dir = args.output_dir
    cores_to_use = args.cores
    assembler = args.assembler

    # Create output subdirectories with names matching the original script
    cleaned_dir = os.path.join(output_dir, "Cleanreads")
    assembly_dir_base = os.path.join(output_dir, "Assembly")
    final_contigs_dir = os.path.join(output_dir, "Contigs")
    os.makedirs(cleaned_dir, exist_ok=True)
    os.makedirs(assembly_dir_base, exist_ok=True)
    os.makedirs(final_contigs_dir, exist_ok=True)

    # Find paired FASTQ files
    pairs = find_fastq_pairs(input_dir)
    if not pairs:
        print("No FASTQ pairs found in input directory.")
        exit(0)
    print(f"Found {len(pairs)} paired-end read sets.")

    # Limit cores to available CPU count
    total_cores = multiprocessing.cpu_count()
    if cores_to_use > total_cores:
        cores_to_use = total_cores
    all_cores = list(range(total_cores))
    assigned_cores = all_cores[:cores_to_use]
    print(f"Assigning tasks to cores: {assigned_cores}")

    # Prepare fastp cleaning tasks
    fastp_tasks = []
    cleaned_skip = []  # Samples already cleaned
    for sample, r1, r2 in pairs:
        # Determine output file names while preserving the original extension
        name = os.path.basename(r1)
        ext = name[name.rfind("_R1") + 3:]  # e.g., ".fq.gz" or ".fastq"
        out_r1 = os.path.join(cleaned_dir, f"{sample}_R1{ext}")
        out_r2 = os.path.join(cleaned_dir, f"{sample}_R2{ext}")
        if os.path.exists(out_r1) and os.path.exists(out_r2):
            cleaned_skip.append(sample)
        else:
            fastp_tasks.append((sample, r1, r2, out_r1, out_r2))  # Include sample in the tuple

    # Run fastp tasks in parallel
    cleaned_success = []  # Samples cleaned successfully
    cleaned_fail = []     # Samples failed during cleaning
    if fastp_tasks:
        print(f"Running fastp on {len(fastp_tasks)} sample(s) in parallel...")
        pool = Pool(processes=min(len(fastp_tasks), cores_to_use))
        result_objs = [pool.apply_async(run_fastp, task[1:]) for task in fastp_tasks]  # Pass r1, r2, out_r1, out_r2 to run_fastp
        pool.close()
        pool.join()
        # Collect fastp results
        for (sample, r1, r2, out_r1, out_r2), res in zip(fastp_tasks, result_objs):
            try:
                status = res.get()
            except Exception as exc:
                status = False
                print(f"[Error] fastp process crashed for sample {sample}: {exc}")
            if status:
                cleaned_success.append(sample)
            else:
                cleaned_fail.append(sample)
    else:
        print("No fastp cleaning needed (all samples already cleaned).")

    # Consider skipped samples as already cleaned
    cleaned_success.extend(cleaned_skip)

    # Run assembly tasks for each sample
    assembly_success = []  # Samples assembled successfully
    assembly_fail = []     # Samples that failed during assembly
    assembly_skip = []     # Samples skipped for assembly
    for sample, r1, r2 in pairs:
        if sample not in cleaned_success:
            assembly_skip.append(sample)
            continue
        # Define cleaned read file paths
        name = os.path.basename(r1)
        ext = name[name.rfind("_R1") + 3:]
        clean_r1 = os.path.join(cleaned_dir, f"{sample}_R1{ext}")
        clean_r2 = os.path.join(cleaned_dir, f"{sample}_R2{ext}")
        if not (os.path.exists(clean_r1) and os.path.exists(clean_r2)):
            print(f"[Warning] Cleaned reads for {sample} not found, skipping assembly.")
            assembly_skip.append(sample)
            continue
        # Create a directory for the sample's assembly output
        asm_sample_dir = os.path.join(assembly_dir_base, sample) 
        #os.makedirs(asm_sample_dir, exist_ok=True)#No,prevent megahit from creating a new directory
        # Check if assembly results already exist
        contig_file = os.path.join(asm_sample_dir,
                                   "final.contigs.fa" if assembler.lower() == "megahit" else "contigs.fasta")
        if os.path.exists(contig_file) and os.path.getsize(contig_file) > 0:
            assembly_skip.append(sample)
        else:
            print(f"Running {assembler} assembly for sample {sample}...")
            success = run_assembler(clean_r1, clean_r2, assembler, asm_sample_dir, cores_to_use)
            if success:
                assembly_success.append(sample)
            else:
                assembly_fail.append(sample)
                continue
        # Copy the final contigs file to the final_contigs directory
        if os.path.exists(contig_file):
            final_name = f"{sample}.fasta"
            final_path = os.path.join(final_contigs_dir, final_name)
            try:
                shutil.copy(contig_file, final_path)
            except Exception as e:
                print(f"[Error] Failed to copy contigs for sample {sample}: {e}")

    # Print pipeline summary
    print("\n=== Pipeline Summary ===")
    print(f"Total samples found: {len(pairs)}")
    if cleaned_success:
        print(f"fastp completed for {len(cleaned_success)} sample(s): {', '.join(cleaned_success)}")
    if cleaned_skip:
        print(f"fastp skipped for {len(cleaned_skip)} sample(s) (already cleaned): {', '.join(cleaned_skip)}")
    if cleaned_fail:
        print(f"fastp failed for {len(cleaned_fail)} sample(s): {', '.join(cleaned_fail)}")
    if assembly_success:
        print(f"Assembly ({assembler}) completed for {len(assembly_success)} sample(s): {', '.join(assembly_success)}")
    if assembly_skip:
        print(f"Assembly skipped for {len(assembly_skip)} sample(s) (already assembled or not processed): {', '.join(assembly_skip)}")
    if assembly_fail:
        print(f"Assembly failed for {len(assembly_fail)} sample(s): {', '.join(assembly_fail)}")
    print("Final contigs are available in:", final_contigs_dir)

if __name__ == "__main__":
    main()