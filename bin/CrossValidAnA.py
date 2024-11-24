import os
import pandas as pd
import sys

# Accept command line arguments
Genomadpath = sys.argv[1]
Viralverifypath = sys.argv[2]
Virsorterpath = sys.argv[3]
Inputfile = sys.argv[4]
OUT_DIR = sys.argv[5]
CONCENTRATION_TYPE = sys.argv[6]

#print(f"{Inputfile}_viurs_summary.tsv")
#print(find_file(Viralverifypath, f"{Inputfile}_result_table.csv"))
#print(find_file(Genomadpath, f"{Inputfile}_virus_summary.tsv"))

# Find file
def find_file(directory, filename):
    for root, dirs, files in os.walk(directory):
        if filename in files:
            return os.path.join(root, filename)
    return None

# Read and filter Virsorter2 data
def read_and_filter_virsorter2(file_path):
    data = pd.read_csv(file_path, sep='\t')
    if CONCENTRATION_TYPE == "non-concentration":
        filtered_data = data[
            ((data['max_score'] >= 0.95) & (data['hallmark'] >= 1)) |
            ((data['hallmark'] == 0) & (data['max_score'] > 0.99) & (data['max_score_group'] != 'ssDNA')) |
            ((data['max_score_group'] == 'ssDNA') & (data['max_score'] > 0.995) & (data['hallmark'] == 0) )
        ]
    else:
        filtered_data = data[
            (data['max_score'] >= 0.90) & (data['hallmark'] >= 1)
        ]
    # 处理序列名
    filtered_data.iloc[:, 0] = filtered_data.iloc[:, 0].apply(lambda x: x.split('||')[0] if pd.notnull(x) else x)
    return filtered_data.iloc[:, 0]

# Read and filter Genomad data
def read_and_filter_genomad(Path):
    filename = Inputfile + "_virus_summary.tsv"
    found_path = find_file(Path, filename)
    data = pd.read_csv(found_path, sep='\t')
    if CONCENTRATION_TYPE == "non-concentration":
        filtered_data = data[
            ((data['virus_score'] > 0.975) & (data['n_hallmarks'] >= 1) & (data['fdr'] <= 0.05))|
            ((data['n_hallmarks'] == 0) & (data['virus_score'] > 0.995) & (data['fdr'] <= 0.05)) 
            ]
    else:
        filtered_data = data[(data['virus_score'] > 0.90) & (data['fdr'] <= 0.05)]
    filtered_data.iloc[:, 0] = filtered_data.iloc[:, 0].apply(lambda x: x.split('||')[0] if pd.notnull(x) else x)
    return filtered_data.iloc[:, 0]

# Read and filter ViralVerify data
def read_and_filter_viralverify(Path):
    filename = Inputfile + "_result_table.csv"
    found_path = find_file(Path, filename)
    data = pd.read_csv(found_path)
    filtered_data = data[data['Prediction'] == "Virus"]
    filtered_data.iloc[:, 0] = filtered_data.iloc[:, 0].apply(lambda x: x.split('||')[0] if pd.notnull(x) else x)
    return filtered_data.iloc[:, 0]

# 合并列表
def merge_lists(*args):
    combined_set = set().union(*[set(list_) for list_ in args])
    return pd.Series(list(combined_set))

# 查找共同元素
def find_common_elements(*args):
    sets = [set(list_) for list_ in args]
    common_elements = set.intersection(*sets)
    return pd.Series(list(common_elements))

# Read and filter non-viral sequences from ViralVerify data (Plasmid, Chromosome, Uncertain - plasmid or chromosomal)
def read_and_filter_viralverify_nonviral(Path):
    filename = Inputfile + "_result_table.csv"
    found_path = find_file(Path, filename)
    data = pd.read_csv(found_path)
    filtered_data = data[data['Prediction'].isin(["Plasmid", "Chromosome", "Uncertain - plasmid or chromosomal"])]
    filtered_data.iloc[:, 0] = filtered_data.iloc[:, 0].apply(lambda x: x.split('||')[0] if pd.notnull(x) else x)
    return set(filtered_data.iloc[:, 0])

# Read and filter plasmid data from Genomad data
def read_plasmid_data(Path):
    filename = Inputfile + "_plasmid_summary.tsv"
    found_path = find_file(Path, filename)
    if found_path:
        data = pd.read_csv(found_path, sep='\t')
        plasmid_sequences = data.iloc[:, 0].apply(lambda x: x.split('||')[0] if pd.notnull(x) else x)
        return set(plasmid_sequences)
    return set()

def merge_lists(*args):
    # Combine lists into a set for unique identifiers
    combined_set = set().union(*[set(list_) for list_ in args])
    return pd.Series(list(combined_set))

def find_common_elements(*args):
    sets = [set(list_) for list_ in args]
    common_elements = set.intersection(*sets)
    return pd.Series(list(common_elements))

# Filter based on CONCENTRATION_TYPE for Pass1
if CONCENTRATION_TYPE == "concentration":
    virsorter2_list1 = read_and_filter_virsorter2(os.path.join(Virsorterpath, "final-viral-score.tsv"))
    genomad_list1 = read_and_filter_genomad(os.path.join(Genomadpath))
    viralverify_list1 = read_and_filter_viralverify(os.path.join(Viralverifypath))

    # 保存各自的 pass1 list
    virsorter2_list1.to_csv(os.path.join(OUT_DIR, Inputfile + "_virsorter2_pass1_list.csv"), index=False)
    genomad_list1.to_csv(os.path.join(OUT_DIR, Inputfile + "_genomad_pass1_list.csv"), index=False)
    viralverify_list1.to_csv(os.path.join(OUT_DIR, Inputfile + "_viralverify_pass1_list.csv"), index=False)
    # 合并 Pass1 结果
    Pass1_list = merge_lists(virsorter2_list1, genomad_list1, viralverify_list1)
    Pass1_list.to_csv(os.path.join(OUT_DIR, Inputfile + "_noRemove_list.csv"), index=False)
    # Merge Pass1 results for concentration type
    AllPass_series = merge_lists(virsorter2_list1, genomad_list1, viralverify_list1)
else:  # Non-concentration type
    genomad_list1 = read_and_filter_genomad(os.path.join(Genomadpath))
    viralverify_list1 = read_and_filter_viralverify(os.path.join(Viralverifypath))
    virsorter2_list1 = read_and_filter_virsorter2(os.path.join(Virsorterpath, "final-viral-score.tsv"))

    # 保存各自的 pass1 list
    virsorter2_list1.to_csv(os.path.join(OUT_DIR, Inputfile + "_virsorter2_pass1_list.csv"), index=False)
    genomad_list1.to_csv(os.path.join(OUT_DIR, Inputfile + "_genomad_pass1_list.csv"), index=False)
    viralverify_list1.to_csv(os.path.join(OUT_DIR, Inputfile + "_viralverify_pass1_list.csv"), index=False)
    # 合并 Pass1 结果
    Pass1_list = merge_lists(virsorter2_list1, genomad_list1, viralverify_list1)
    Pass1_list.to_csv(os.path.join(OUT_DIR, Inputfile + "_noRemove_list.csv"), index=False)

    # Merge Pass1 results for non-concentration type
    AllPass_series = merge_lists(genomad_list1, viralverify_list1)

# Read and filter plasmid data from Genomad data
plasmid_sequences = read_plasmid_data(Genomadpath)

# Read and filter non-viral sequences (Plasmid, Chromosome, Uncertain - plasmid or chromosomal) from ViralVerify data
nonviral_sequences = read_and_filter_viralverify_nonviral(Viralverifypath)

# 将 set 转换为 pandas.Series 或 pandas.DataFrame
pd.Series(list(plasmid_sequences)).to_csv(
    os.path.join(OUT_DIR, Inputfile + "_GenomadPlsmid_list.csv"), index=False, header=False
)

# 对非病毒序列的输出同样修改
pd.Series(list(nonviral_sequences)).to_csv(
    os.path.join(OUT_DIR, Inputfile + "_VVFPlsmid_list.csv"), index=False, header=False
)

# Remove all sequences in plasmid_sequences and nonviral_sequences from AllPass_series
AllPass_series = AllPass_series[~AllPass_series.isin(plasmid_sequences | nonviral_sequences)]

# Create DataFrame
AllPass_df = pd.DataFrame(AllPass_series, columns=['Sequence Id'])

# Save to file
filename = Inputfile + "_viral_predictionsList.csv"
AllPass_df.to_csv(os.path.join(OUT_DIR, filename), index=False)


print("过滤和合并完成，文件已保存。")
