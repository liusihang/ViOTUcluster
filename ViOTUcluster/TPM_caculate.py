#!/usr/bin/env python3
import pandas as pd
import os
import sys

def calculate_tpm(tsv_file):
    """
    Calculate TPM (Transcripts Per Million) from a TSV file.

    Parameters:
    - tsv_file (str): Path to the TSV file containing gene expression data.

    Returns:
    - pd.DataFrame: A DataFrame with 'Sequence Id' as the index and 'TPM' as the column.
      Returns an empty DataFrame if the input file is empty or an error occurs.
    """
    try:
        df = pd.read_csv(tsv_file, sep='\t')
        
        # Calculate RPKM (Reads Per Kilobase of transcript, per Million mapped reads)
        df['RPKM'] = df['Mapped reads'] / (df['Sequence length (bp)'] / 1000)
        
        # Calculate the total number of mapped reads in millions
        total_mapped_reads = df['Mapped reads'].sum() / 1_000_000
        
        # Normalize RPKM to get TPM
        df['RPKM'] = df['RPKM'] / total_mapped_reads
        sum_rpk = df['RPKM'].sum()
        df['TPM'] = df['RPKM'] / sum_rpk * 1_000_000
        
        # Prepare the result with 'Sequence Id' as the index and 'TPM' as the value
        result = df[['Sequence Id', 'TPM']].copy()
        result.set_index('Sequence Id', inplace=True)
        
        return result
    except pd.errors.EmptyDataError:
        print(f"Warning: The file {tsv_file} is empty and will be skipped.")
        return pd.DataFrame() 
    except Exception as e:
        print(f"Error processing file {tsv_file}: {e}")
        return pd.DataFrame() 

def merge_tpm_files(input_folder, merged_output_file):
    """
    Merge TPM values from multiple TSV files into a single CSV file.

    Parameters:
    - input_folder (str): Path to the folder containing input TSV files.
    - merged_output_file (str): Path to save the merged TPM CSV file.
    """
    merged_df = pd.DataFrame()
    
    # Iterate through all TSV files in the input folder
    for file_name in os.listdir(input_folder):
        if file_name.endswith('.tsv'):
            file_path = os.path.join(input_folder, file_name)
            tpm_df = calculate_tpm(file_path)
            
            # If TPM data is available, add it to the merged DataFrame
            if not tpm_df.empty:
                column_name = os.path.splitext(file_name)[0]
                merged_df[column_name] = tpm_df['TPM']
    
    # If merged DataFrame is not empty, process and save it
    if not merged_df.empty:
        # Sort columns alphabetically
        merged_df = merged_df.sort_index(axis=1)
        
        # Remove specific suffix from column names
        merged_df.columns = merged_df.columns.str.replace('_coverage_TPM', '', regex=False)
        
        # Save the merged TPM data to a CSV file
        merged_df.to_csv(merged_output_file)
        print(f"Merged TPM file saved to {merged_output_file}")
    else:
        print("No valid TPM data found to merge.")

def main():
    """
    Main function to execute the TPM calculation and merging process.
    
    Expects two command-line arguments:
    1. Input folder containing TSV files.
    2. Output CSV file path for the merged TPM data.
    """
    if len(sys.argv) != 3:
        print("Usage: python script.py <input_folder> <merged_output_file>")
        sys.exit(1)
    
    input_folder = sys.argv[1]
    merged_output_file = sys.argv[2]
    
    merge_tpm_files(input_folder, merged_output_file)

if __name__ == "__main__":
    main()