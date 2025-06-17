# -*- coding: utf-8 -*-
#!/usr/bin/env python3

import sys
import os
import shutil # For os.replace or a fallback
from Bio import SeqIO

def _filter_single_fasta(min_length, input_file_path, output_file_path):
    """
    Helper function to filter a single FASTA file.
    If input_file_path and output_file_path are the same, it performs an in-place update safely.
    """
    # Determine if this is an in-place operation for this specific file
    is_effectively_inplace = (os.path.abspath(input_file_path) == os.path.abspath(output_file_path))
    
    # If in-place, write to a temporary file first in the same directory as the target output_file_path
    # Otherwise, temp_output_file_path is just the final output_file_path
    if is_effectively_inplace:
        # Create temp file in the same directory as the target output_file_path
        # Using a slightly more robust temp name to avoid potential collisions if script is run multiple times quickly
        base, ext = os.path.splitext(output_file_path)
        temp_output_file_path = f"{base}.tmpfilter{os.getpid()}{ext}"
    else:
        temp_output_file_path = output_file_path 

    processed_count = 0
    success = False
    try:
        with open(input_file_path, 'r') as input_handle, open(temp_output_file_path, 'w') as output_handle:
            sequences = SeqIO.parse(input_handle, 'fasta')
            filtered_seqs = (seq for seq in sequences if len(seq.seq) >= min_length)
            processed_count = SeqIO.write(filtered_seqs, output_handle, 'fasta')
        success = True # Reached here means writing to temp/output was successful
        
        # If it was an in-place operation, and writing to temp succeeded, replace the original.
        if is_effectively_inplace:
            # os.replace is atomic on POSIX and Windows (for files on the same filesystem)
            os.replace(temp_output_file_path, output_file_path) 
            print(f"Filtered {processed_count} sequences. File '{output_file_path}' updated in-place.")
        else: # Standard output to a different file/location
            # If not in-place, temp_output_file_path is already the final output_file_path, so no move needed.
            print(f"Filtered {processed_count} sequences from '{input_file_path}' to '{output_file_path}'.")

    except FileNotFoundError:
        print(f"Error: Input file '{input_file_path}' not found.")
        if is_effectively_inplace and os.path.exists(temp_output_file_path):
            os.remove(temp_output_file_path) # Clean up temp if input not found (though unlikely to reach here for temp)
    except Exception as e:
        print(f"Error processing file '{input_file_path}': {e}")
        # Clean up the temporary file if it exists and an error occurred during processing or replacement
        if os.path.exists(temp_output_file_path) and (is_effectively_inplace or not success):
             # If it was in-place, temp_output_file_path is distinct and should be removed.
             # If it was not in-place (temp_output_file_path == output_file_path) AND writing failed (not success),
             # then output_file_path might be partially written or empty, so remove it.
            try:
                os.remove(temp_output_file_path)
            except OSError as rm_err:
                print(f"Warning: Could not remove temporary/output file '{temp_output_file_path}': {rm_err}")
    finally:
        # This finally block is mostly for the case where is_effectively_inplace is true
        # and an error happened AFTER successful write to temp_output_file_path but BEFORE os.replace
        # or if os.replace itself failed. os.replace should handle the temp file itself on success.
        # If not in_place, and an error happened, the partially written output_file_path
        # (which is temp_output_file_path in that branch) would have been handled by the except block.
        if is_effectively_inplace and not success and os.path.exists(temp_output_file_path):
            try:
                os.remove(temp_output_file_path) # Final cleanup attempt for temp file on error
            except OSError:
                pass


def filter_sequences_flexible(min_length, input_path, output_dir):
    """
    Filters sequences from FASTA files in an input directory or a single input FASTA file.

    Args:
        min_length (int): The minimum length for a sequence to be kept.
        input_path (str): Path to the input directory or a single FASTA file.
        output_dir (str): Path to the output directory where filtered files will be saved.
                          If this is the same as the input directory (or input file's directory),
                          files will be updated in-place.
    """
    # If output_dir doesn't exist AND it's not intended for in-place modification
    # (i.e. output_dir is different from input_path or input_path's parent)
    # then create it.
    # For in-place, the output_dir (which is the input_dir) must already exist.
    if not os.path.exists(output_dir):
        # Heuristic: if input_path is a dir and output_dir is different, create output_dir.
        # If input_path is a file and output_dir is different from its parent, create output_dir.
        is_input_dir = os.path.isdir(input_path)
        input_base_dir = input_path if is_input_dir else os.path.dirname(os.path.abspath(input_path))
        
        if os.path.abspath(input_base_dir) != os.path.abspath(output_dir):
            os.makedirs(output_dir)
            print(f"Created output directory: {output_dir}")
        # If output_dir is same as input_base_dir but doesn't exist, something is wrong with input.
        # This case should be caught by the input_path existence check in __main__.


    processed_any_fasta = False

    if os.path.isdir(input_path):
        print(f"Processing directory: {input_path}")
        for filename in os.listdir(input_path):
            if filename.lower().endswith(('.fasta', '.fa', '.fna')):
                # Avoid processing our own temporary files if script is interrupted and run again
                if filename.startswith(".tmpfilter") and filename.endswith(os.path.splitext(filename)[1]): # Basic check
                    continue

                current_input_file_path = os.path.join(input_path, filename)
                current_output_file_path = os.path.join(output_dir, filename) # This might be same as input
                _filter_single_fasta(min_length, current_input_file_path, current_output_file_path)
                processed_any_fasta = True
        if not processed_any_fasta:
             print(f"No FASTA files found in directory: {input_path}")

    elif os.path.isfile(input_path):
        print(f"Processing single file: {input_path}")
        if input_path.lower().endswith(('.fasta', '.fa', '.fna')):
            filename = os.path.basename(input_path)
            # current_output_file_path could be same as input_path if output_dir is dirname(input_path)
            current_output_file_path = os.path.join(output_dir, filename)
            _filter_single_fasta(min_length, input_path, current_output_file_path)
            processed_any_fasta = True
        else:
            print(f"Warning: Input file '{input_path}' does not appear to be a FASTA file "
                  "(expected .fasta, .fa, or .fna extension). Skipping.")
    else:
        # This case should ideally be caught by the check in __main__
        print(f"Error: Input path '{input_path}' is not a valid file or directory.")
        sys.exit(1)

    # Create a completion flag file.
    # Consider if 'Done' file is appropriate if output_dir was the same as input_dir (in-place).
    # For now, it will create 'Done' in the output_dir regardless.
    # If operating in-place, this means 'Done' is in the input_dir.
    if processed_any_fasta or os.path.isdir(input_path):
        done_file_path = os.path.join(output_dir, 'Done')
        # Check if we should actually write the Done file for in-place ops.
        # For simplicity, current logic will write it. User must remove it to re-filter.
        try:
            with open(done_file_path, 'w') as f:
                f.write('Filtering completed successfully for this output directory.')
            print(f"Created completion flag: {done_file_path}")
        except IOError as e:
            print(f"Warning: Could not write Done file '{done_file_path}': {e}")


if __name__ == '__main__':
    if len(sys.argv) != 4:
        print("Usage: python script.py <min_length> <input_path_or_dir> <output_dir>")
        print("  <input_path_or_dir>: Path to a directory containing FASTA files or a single FASTA file.")
        print("  <output_dir>: Directory where filtered FASTA files will be saved.")
        print("                To modify files in-place, set <output_dir> to be the same as <input_path_or_dir> (if input is a dir)")
        print("                or the directory containing the input file (if input is a single file).")
        sys.exit(1)

    try:
        min_len_arg = int(sys.argv[1])
        if min_len_arg < 0:
            print("Error: min_length must be a non-negative integer.")
            sys.exit(1)
    except ValueError:
        print("Error: min_length must be an integer.")
        sys.exit(1)

    input_path_arg = sys.argv[2]
    output_dir_arg = sys.argv[3]

    # Ensure input path exists
    if not os.path.exists(input_path_arg):
        print(f"Error: Input path '{input_path_arg}' does not exist.")
        sys.exit(1)

    # Resolve paths to be absolute for reliable comparison for 'Done' file logic
    abs_input_path = os.path.abspath(input_path_arg)
    abs_output_dir = os.path.abspath(output_dir_arg)

    # Determine the effective input base directory for comparison with output directory
    if os.path.isdir(abs_input_path):
        abs_input_base_dir = abs_input_path
    else: # input is a file
        abs_input_base_dir = os.path.dirname(abs_input_path)

    # 'Done' file logic:
    # If the output directory is different from the input source directory,
    # then the 'Done' file signifies completion for that *distinct output directory*.
    # If output_dir is the same as the input source (in-place update), the 'Done' file
    # will be in the source directory and prevent re-filtering unless removed.
    # This behavior is kept from the original script.
    done_file = os.path.join(abs_output_dir, 'Done') # Use absolute path for output_dir
    if os.path.exists(done_file):
        print(f"Filtering already completed for output directory '{output_dir_arg}' (found 'Done' file). Exiting.")
        sys.exit(0)
    else:
        # Ensure output directory exists if it's not the same as the input's base directory
        # This is particularly for when input_path_arg is a file, and output_dir_arg is a new dir.
        # If output_dir_arg is meant to be the input file's directory for in-place, it already exists.
        if abs_input_base_dir != abs_output_dir and not os.path.exists(abs_output_dir):
            try:
                os.makedirs(abs_output_dir)
                print(f"Created output directory: {output_dir_arg}")
            except OSError as e:
                print(f"Error: Could not create output directory '{output_dir_arg}': {e}")
                sys.exit(1)
        elif not os.path.isdir(abs_output_dir): # If output_dir exists but is not a directory
             print(f"Error: Specified output_dir '{output_dir_arg}' exists but is not a directory.")
             sys.exit(1)


        filter_sequences_flexible(min_len_arg, input_path_arg, output_dir_arg) # Pass original args
        print(f"Filtering process finished. Check '{output_dir_arg}' for results.")