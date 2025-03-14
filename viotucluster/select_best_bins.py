#!/usr/bin/env python3
import os
import sys
import pandas as pd
import shutil

def select_best_bins(quality_summary_file, extracted_dir, final_bins_dir):
    """
    Select the best bins based on completeness and contamination metrics,
    and copy the corresponding FASTA files to the final_bins directory.

    Parameters:
    - quality_summary_file (str): Path to the quality_summary.tsv file.
    - extracted_dir (str): Directory containing extracted bin FASTA files.
    - final_bins_dir (str): Directory where the best bins will be copied.
    """
    # Read the quality_summary.tsv file
    df = pd.read_csv(quality_summary_file, sep='\t')

    # Create the final_bins directory if it doesn't exist
    os.makedirs(final_bins_dir, exist_ok=True)

    # Extract bin IDs by removing suffixes like .permissive, .origin, .strict
    df['bin_id'] = df['contig_id'].apply(lambda x: x.rsplit('.', 1)[0])
    unique_bins = df['bin_id'].unique()

    best_bins = []

    # Iterate over each unique bin ID
    for bin_id in unique_bins:
        bin_df = df[df['bin_id'] == bin_id]
        
        # Filter out bins with contamination > 1
        filtered_bin_df = bin_df[bin_df['contamination'] <= 1]

        if not filtered_bin_df.empty:
            # Select the bin with the highest completeness
            best_bin = filtered_bin_df.loc[filtered_bin_df['completeness'].idxmax()]
        else:
            # If all contamination values > 1, retain the origin bin
            origin_bin_df = bin_df[bin_df['contig_id'].str.endswith(".origin")]
            if not origin_bin_df.empty:
                best_bin = origin_bin_df.iloc[0]
            else:
                print(f"Warning: Origin bin not found for {bin_id}. Skipping...")
                continue

        best_bins.append(best_bin)

        # Copy the best bin's FASTA file to the final_bins directory
        best_bin_filename = best_bin['contig_id'] + ".fasta"
        src_path = os.path.join(extracted_dir, best_bin_filename)
        dst_path = os.path.join(final_bins_dir, best_bin_filename)

        if os.path.exists(src_path):
            shutil.copy(src_path, dst_path)
        else:
            print(f"Warning: {src_path} does not exist. Skipping...")

    # Save the best bins information to a summary TSV file
    best_bins_df = pd.DataFrame(best_bins)
    best_bins_df.to_csv(os.path.join(final_bins_dir, 'best_bins_summary.tsv'), sep='\t', index=False)

    print(f"Best bins copied to {final_bins_dir} and summary saved as best_bins_summary.tsv")

def main():
    if len(sys.argv) != 4:
        #print("Usage: python select_best_bins.py <quality_summary.tsv> <EXTRACTED_DIR> <final_bins_dir>")
        sys.exit(1)

    quality_summary_file = sys.argv[1]
    extracted_dir = sys.argv[2]
    final_bins_dir = sys.argv[3]

    select_best_bins(quality_summary_file, extracted_dir, final_bins_dir)

if __name__ == "__main__":
    main()