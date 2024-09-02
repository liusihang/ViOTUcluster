#!/usr/bin/env python
import pandas as pd
import os
import sys

def calculate_tpm(tsv_file):
    try:
        df = pd.read_csv(tsv_file, sep='\t')
        df['RPKM'] = df['Mapped reads'] / (df['Sequence length (bp)'] / 1000)
        total_mapped_reads = df['Mapped reads'].sum() / 1000000
        df['RPKM'] = df['RPKM'] / total_mapped_reads
        sum_rpk = df['RPKM'].sum()
        df['TPM'] = df['RPKM'] / sum_rpk * 1000000
        result = df[['Sequence Id', 'TPM']]
        result.set_index('Sequence Id', inplace=True)
        return result
    except pd.errors.EmptyDataError:
        print(f"Warning: The file {tsv_file} is empty and will be skipped.")
        return pd.DataFrame()  # 返回空数据框
    except Exception as e:
        print(f"Error processing file {tsv_file}: {e}")
        return pd.DataFrame()  # 返回空数据框

def merge_tpm_files(input_folder, merged_output_file):
    merged_df = pd.DataFrame()
    
    for file_name in os.listdir(input_folder):
        if file_name.endswith('.tsv'):
            file_path = os.path.join(input_folder, file_name)
            tpm_df = calculate_tpm(file_path)
            if not tpm_df.empty:
                column_name = os.path.splitext(file_name)[0]
                merged_df[column_name] = tpm_df['TPM']
    
    if not merged_df.empty:
        merged_df = merged_df.sort_index(axis=1)
        merged_df.columns = merged_df.columns.str.replace('_coverage_TPM', '')
        merged_df.to_csv(merged_output_file)
        print(f"Merged TPM file saved to {merged_output_file}")
    else:
        print("No valid TPM data found to merge.")

def main():
    if len(sys.argv) != 3:
        print("Usage: python script.py <input_folder> <merged_output_file>")
        sys.exit(1)
    
    input_folder = sys.argv[1]
    merged_output_file = sys.argv[2]
    
    merge_tpm_files(input_folder, merged_output_file)

if __name__ == "__main__":
    main()