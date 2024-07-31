#!/usr/bin/env python3
import os
import sys
import pandas as pd
import shutil

def select_best_bins(quality_summary_file, extracted_dir, final_bins_dir):
    # 读取quality_summary.tsv文件
    df = pd.read_csv(quality_summary_file, sep='\t')

    # 创建存放最佳bin的final_bins文件夹
    if not os.path.exists(final_bins_dir):
        os.makedirs(final_bins_dir)

    # 获取所有的bin ID（不包括后缀.permissive, .origin, .strict）
    df['bin_id'] = df['contig_id'].apply(lambda x: x.rsplit('.', 1)[0])
    unique_bins = df['bin_id'].unique()

    best_bins = []

    # 遍历每个唯一的bin ID
    for bin_id in unique_bins:
        bin_df = df[df['bin_id'] == bin_id]
        
        # 过滤掉contamination值大于1的bin
        filtered_bin_df = bin_df[bin_df['contamination'] <= 1]

        if not filtered_bin_df.empty:
            # 找到completeness值最大的bin
            best_bin = filtered_bin_df.loc[filtered_bin_df['completeness'].idxmax()]
        else:
            # 所有的contamination都大于1，则保留origin bin
            origin_bin_df = bin_df[bin_df['contig_id'].str.endswith(".origin")]
            if not origin_bin_df.empty:
                best_bin = origin_bin_df.iloc[0]
            else:
                print(f"Warning: Origin bin not found for {bin_id}. Skipping...")
                continue

        best_bins.append(best_bin)

        # 从EXTRACTED_DIR中复制最佳bin到final_bins文件夹中
        best_bin_filename = best_bin['contig_id'] + ".fasta"
        src_path = os.path.join(extracted_dir, best_bin_filename)
        dst_path = os.path.join(final_bins_dir, best_bin_filename)

        if os.path.exists(src_path):
            shutil.copy(src_path, dst_path)
        else:
            print(f"Warning: {src_path} does not exist. Skipping...")

    # 将最佳bin的信息保存到一个新的文件中
    best_bins_df = pd.DataFrame(best_bins)
    best_bins_df.to_csv(os.path.join(final_bins_dir, 'best_bins_summary.tsv'), sep='\t', index=False)

    print(f"Best bins copied to {final_bins_dir} and summary saved as best_bins_summary.tsv")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python select_best_bins.py <quality_summary.tsv> <EXTRACTED_DIR> <final_bins_dir>")
        sys.exit(1)

    quality_summary_file = sys.argv[1]
    extracted_dir = sys.argv[2]
    final_bins_dir = sys.argv[3]

    select_best_bins(quality_summary_file, extracted_dir, final_bins_dir)
