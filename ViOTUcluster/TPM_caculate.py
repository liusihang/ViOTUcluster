import os
import sys
import pandas as pd

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
            # Read TSV file with OTU as index and assuming TPM is in the second column
            try:
                tpm_df = pd.read_csv(file_path, sep='\t', index_col=0)
            except Exception as e:
                print(f"Error reading {file_name}: {e}")
                continue

            # Check if 'TPM' column exists, or使用第二列数据
            tpm_series = tpm_df.iloc[:, 0]

            # Rename the series using the file name (without extension)
            column_name = os.path.splitext(file_name)[0]
            tpm_series.name = column_name

            # 合并数据，按 OTU 名进行对齐
            merged_df = pd.concat([merged_df, tpm_series], axis=1)

    # If merged DataFrame is not empty, process and save it
    if not merged_df.empty:
        # Sort columns alphabetically
        merged_df = merged_df.sort_index(axis=1)

        # Remove '_coverage' from column names if present
        merged_df.columns = merged_df.columns.str.replace('_coverage', '', regex=False)

        # Set the index name to "OTU" so that the first column in the CSV is named OTU
        merged_df.index.name = "OTU"

        # Save the merged TPM data to a CSV file
        merged_df.to_csv(merged_output_file, float_format='%.10f')
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
