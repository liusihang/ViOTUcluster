import os
import pandas as pd
import sys

# Accept command-line arguments
Genomadpath = sys.argv[1]
Viralverifypath = sys.argv[2]
Virsorterpath = sys.argv[3]
Inputfile = sys.argv[4]
OUT_DIR = sys.argv[5]

#print(f"{Inputfile}_viurs_summary.tsv")
#print(find_file(Viralverifypath, f"{Inputfile}_result_table.csv"))
#print(find_file(Genomadpath, f"{Inputfile}_virus_summary.tsv"))

# 查找文件函数
def find_file(directory, filename):
    for root, dirs, files in os.walk(directory):
        if filename in files:
            return os.path.join(root, filename)
    return None

# 读取和过滤 Virsorter2 数据
def read_and_filter_virsorter2(file_path, pass_type='pass1'):
    data = pd.read_csv(file_path, sep='\t')
    if pass_type == 'pass1':
        filtered_data = data[(data['max_score'] > 0.95) | (data['hallmark'] >= 2)]
    else:
        filtered_data = data[(data['max_score'] <= 0.95) & (data['max_score'] > 0.6)]
    filtered_data.iloc[:, 0] = filtered_data.iloc[:, 0].apply(lambda x: x.split('||')[0] if pd.notnull(x) else x)
    return filtered_data.iloc[:, 0]

# 读取和过滤 Genomad 数据
def read_and_filter_genomad(Path, pass_type='pass1'):
    filename = Inputfile + "_virus_summary.tsv"
    found_path = find_file(Path, filename)
    data = pd.read_csv(found_path, sep='\t')
    if pass_type == 'pass1':
        filtered_data = data[(data['virus_score'] > 0.8) & (data['n_hallmarks'] >= 1) & (data['fdr'] <= 0.05)]
    else:
        filtered_data = data[(data['virus_score'] > 0.6) & (data['virus_score'] <= 0.8) & (data['fdr'] <= 0.05)]
    filtered_data.iloc[:, 0] = filtered_data.iloc[:, 0].apply(lambda x: x.split('||')[0] if pd.notnull(x) else x)
    return filtered_data.iloc[:, 0]

# 读取和过滤 ViralVerify 数据
def read_and_filter_viralverify(Path, pass_type='pass1'):
    filename = Inputfile + "_result_table.csv"
    found_path = find_file(Path, filename)
    data = pd.read_csv(found_path)
    if pass_type == 'pass1':
        filtered_data = data[data['Prediction'] == "Virus"]
    else:
        filtered_data = data[data['Prediction'] == "Uncertain - viral or bacterial"]
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

# Pass1 过滤
virsorter2_list1 = read_and_filter_virsorter2(os.path.join(Virsorterpath, "final-viral-score.tsv"), 'pass1')
genomad_list1 = read_and_filter_genomad(os.path.join(Genomadpath), 'pass1')
viralverify_list1 = read_and_filter_viralverify(os.path.join(Viralverifypath), 'pass1')

# 保存各自的 pass1 list
virsorter2_list1.to_csv(os.path.join(OUT_DIR, Inputfile + "_virsorter2_pass1_list.csv"), index=False)
genomad_list1.to_csv(os.path.join(OUT_DIR, Inputfile + "_genomad_pass1_list.csv"), index=False)
viralverify_list1.to_csv(os.path.join(OUT_DIR, Inputfile + "_viralverify_pass1_list.csv"), index=False)

# 合并 Pass1 结果
Pass1_list = merge_lists(virsorter2_list1, genomad_list1, viralverify_list1)
Pass1_list.to_csv(os.path.join(OUT_DIR, Inputfile + "_combined_pass1_list.csv"), index=False)

# Double-check 过滤
virsorter2_list2 = read_and_filter_virsorter2(os.path.join(Virsorterpath, "final-viral-score.tsv"), 'doublecheck')
genomad_list2 = read_and_filter_genomad(os.path.join(Genomadpath), 'doublecheck')
viralverify_list2 = read_and_filter_viralverify(os.path.join(Viralverifypath), 'doublecheck')

# 合并 double-check 结果
Pass2_list = find_common_elements(virsorter2_list2, genomad_list2, viralverify_list2)

# 合并所有 Pass 结果
AllPass_series = merge_lists(Pass1_list, Pass2_list)

# 创建 DataFrame
AllPass_df = pd.DataFrame(AllPass_series, columns=['Sequence Id'])

# 保存所有 Pass 结果
AllPass_df.to_csv(os.path.join(OUT_DIR, Inputfile + "_all_pass_list.csv"), index=False)

print("过滤和合并完成，文件已保存。")
