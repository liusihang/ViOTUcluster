#!/usr/bin/env python3
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

def find_file(directory, filename):
    """Find a file in the given directory and its subdirectories."""
    for root, dirs, files in os.walk(directory):
        if filename in files:
            return os.path.join(root, filename)
    return None

def read_and_filter_virsorter2(file_path):
    """Read and filter Virsorter2 data based on concentration type."""
    data = pd.read_csv(file_path, sep='\t')
    if CONCENTRATION_TYPE == "non-concentration":
        filtered_data = data[
            ((data['max_score'] >= 0.95) & (data['hallmark'] >= 1)) |
            ((data['hallmark'] == 0) & (data['max_score'] > 0.99) & (data['max_score_group'] != 'ssDNA')) |
            ((data['max_score_group'] == 'ssDNA') & (data['max_score'] > 0.995) & (data['hallmark'] == 0))
        ]
    else:
        filtered_data = data[
            (data['max_score'] >= 0.90) & (data['hallmark'] >= 1)
        ]
    # Process sequence names
    filtered_data.iloc[:, 0] = filtered_data.iloc[:, 0].apply(lambda x: x.split('||')[0] if pd.notnull(x) else x)
    return filtered_data.iloc[:, 0]

def read_and_filter_genomad(path):
    """Read and filter Genomad data based on concentration type."""
    filename = f"{Inputfile}_virus_summary.tsv"
    found_path = find_file(path, filename)
    data = pd.read_csv(found_path, sep='\t')
    if CONCENTRATION_TYPE == "non-concentration":
        filtered_data = data[
            ((data['virus_score'] > 0.8) & (data['n_hallmarks'] >= 1) & (data['fdr'] <= 0.05)) |
            ((data['n_hallmarks'] == 0) & (data['virus_score'] > 0.995) & (data['fdr'] <= 0.05))
        ]
    else:
        filtered_data = data[(data['virus_score'] > 0.7) & (data['fdr'] <= 0.05)]
    # Process sequence names
    filtered_data.iloc[:, 0] = filtered_data.iloc[:, 0].apply(lambda x: x.split('||')[0] if pd.notnull(x) else x)
    return filtered_data.iloc[:, 0]

def read_and_filter_viralverify(path):
    """Read and filter ViralVerify data to include only virus predictions."""
    filename = f"{Inputfile}_result_table.csv"
    found_path = find_file(path, filename)
    data = pd.read_csv(found_path)
    filtered_data = data[data['Prediction'] == "Virus"]
    # Process sequence names
    filtered_data.iloc[:, 0] = filtered_data.iloc[:, 0].apply(lambda x: x.split('||')[0] if pd.notnull(x) else x)
    return filtered_data.iloc[:, 0]

def read_and_filter_viralverify_nonviral(path):
    """Read and filter ViralVerify data to include non-viral predictions."""
    filename = f"{Inputfile}_result_table.csv"
    found_path = find_file(path, filename)
    data = pd.read_csv(found_path)
    filtered_data = data[data['Prediction'].isin(["Plasmid", "Chromosome", "Uncertain - plasmid or chromosomal"])]
    # Process sequence names
    filtered_data.iloc[:, 0] = filtered_data.iloc[:, 0].apply(lambda x: x.split('||')[0] if pd.notnull(x) else x)
    return set(filtered_data.iloc[:, 0])

def read_plasmid_data(path):
    """Read and extract plasmid sequences from Genomad data."""
    filename = f"{Inputfile}_plasmid_summary.tsv"
    found_path = find_file(path, filename)
    if found_path:
        data = pd.read_csv(found_path, sep='\t')
        plasmid_sequences = data.iloc[:, 0].apply(lambda x: x.split('||')[0] if pd.notnull(x) else x)
        return set(plasmid_sequences)
    return set()

def merge_lists(*args):
    """Merge multiple lists into a unique set."""
    combined_set = set().union(*[set(list_) for list_ in args])
    return pd.Series(list(combined_set))

def find_common_elements(*args):
    """Find common elements across multiple lists."""
    sets = [set(list_) for list_ in args]
    common_elements = set.intersection(*sets)
    return pd.Series(list(common_elements))

# Read plasmid and non-viral sequences
plasmid_sequences = read_plasmid_data(Genomadpath)
nonviral_sequences = read_and_filter_viralverify_nonviral(Viralverifypath)

# Filter based on concentration type
if CONCENTRATION_TYPE == "concentration":
    virsorter2_list1 = read_and_filter_virsorter2(os.path.join(Virsorterpath, "final-viral-score.tsv"))
    genomad_list1 = read_and_filter_genomad(Genomadpath)
    viralverify_list1 = read_and_filter_viralverify(Viralverifypath)

    # Merge Pass1 results
    Pass1_list = merge_lists(virsorter2_list1, genomad_list1, viralverify_list1)
    AllPass_series = merge_lists(virsorter2_list1, genomad_list1, viralverify_list1)
    # Remove plasmid sequences
    AllPass_series = AllPass_series[~AllPass_series.isin(plasmid_sequences)]
else:  # Non-concentration type
    genomad_list1 = read_and_filter_genomad(Genomadpath)
    viralverify_list1 = read_and_filter_viralverify(Viralverifypath)
    virsorter2_list1 = read_and_filter_virsorter2(os.path.join(Virsorterpath, "final-viral-score.tsv"))
    # Merge Pass1 results
    Pass1_list = merge_lists(virsorter2_list1, genomad_list1, viralverify_list1)
    AllPass_series = merge_lists(genomad_list1, viralverify_list1)
    # Remove plasmid and non-viral sequences
    AllPass_series = AllPass_series[~AllPass_series.isin(plasmid_sequences | nonviral_sequences)]

# Save the results to a CSV file
AllPass_df = pd.DataFrame(AllPass_series, columns=['Sequence Id'])
filename = f"{Inputfile}_viral_predictionsList.csv"
AllPass_df.to_csv(os.path.join(OUT_DIR, filename), index=False)