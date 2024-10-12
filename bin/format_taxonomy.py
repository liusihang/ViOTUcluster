import pandas as pd
import sys
import os

def format_taxonomy(input_csv, output_dir):
    # 读取CSV文件
    df = pd.read_csv(input_csv)

    # 选择 'seq_name' 和 'lineage' 列，复制它们
    df_phyloseq = df[['seq_name', 'lineage']].copy()

    # 将 lineage 列拆分为多个分类级别，假设最多有7个分类级别 (Domain -> Species)
    taxonomy_levels = ['Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
    df_phyloseq[taxonomy_levels] = df_phyloseq['lineage'].str.split(';', expand=True)

    # 遍历每一行，检查是否存在缺失值
    for idx, row in df_phyloseq.iterrows():
        for i in range(1, len(taxonomy_levels)):
            current_level = taxonomy_levels[i]
            previous_level = taxonomy_levels[i - 1]
            # 如果当前级别为空，则用上一级别加 '_unclassified' 填充，且后续所有级别依次填充
            if pd.isna(row[current_level]):
                for j in range(i, len(taxonomy_levels)):
                    df_phyloseq.at[idx, taxonomy_levels[j]] = df_phyloseq.at[idx, taxonomy_levels[j - 1]] + '_unclassified'
                break  # 后续级别已经填充完毕，跳出循环

    # 将 'seq_name' 列重命名为 'OTUname'
    df_phyloseq = df_phyloseq.rename(columns={'seq_name': 'OTUname'})

    # 确保输出目录存在
    os.makedirs(output_dir, exist_ok=True)

    # 生成输出文件路径
    output_file = os.path.join(output_dir, 'phyloseq_taxonomy_formatted.csv')

    # 保存为CSV文件
    df_phyloseq.to_csv(output_file, index=False)

    print(f"文件已保存为: {output_file}")


if __name__ == "__main__":
    # 检查命令行参数
    if len(sys.argv) != 3:
        print("Usage: python script_name.py <input_csv_file> <output_directory>")
        sys.exit(1)

    input_csv_file = sys.argv[1]
    output_directory = sys.argv[2]

    # 调用函数进行格式化
    format_taxonomy(input_csv_file, output_directory)