#! /usr/bin/env python
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

# Find file
def find_file(directory, filename):
    for root, dirs, files in os.walk(directory):
        if filename in files:
            return os.path.join(root, filename)
    return None

# Read and filter Virsorter2 data
def read_and_filter_virsorter2(file_path, pass_type='pass1'):
    data = pd.read_csv(file_path, sep='\t')
    if pass_type == 'pass1':
        filtered_data = data[(data['max_score'] > 0.95) | (data['hallmark'] >= 2)]
    else:
        filtered_data = data[(data['max_score'] <= 0.95) & (data['max_score'] > 0.6)]
    filtered_data.iloc[:, 0] = filtered_data.iloc[:, 0].apply(lambda x: x.split('||')[0] if pd.notnull(x) else x)
    return filtered_data.iloc[:, 0]

# Read and filter Genomad data
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

# Read and filter ViralVerify data
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
    virsorter2_list1 = read_and_filter_virsorter2(os.path.join(Virsorterpath, "final-viral-score.tsv"), 'pass1')
    genomad_list1 = read_and_filter_genomad(os.path.join(Genomadpath), 'pass1')
    viralverify_list1 = read_and_filter_viralverify(os.path.join(Viralverifypath), 'pass1')

    # Merge Pass1 results for concentration type
    AllPass_series = merge_lists(virsorter2_list1, genomad_list1, viralverify_list1)
    # Read and filter plasmid data from Genomad data
    plasmid_sequences = read_plasmid_data(Genomadpath)

    # Read and filter non-viral sequences (Plasmid, Chromosome, Uncertain - plasmid or chromosomal) from ViralVerify data
    nonviral_sequences = read_and_filter_viralverify_nonviral(Viralverifypath)

    # Remove all sequences in plasmid_sequences and nonviral_sequences from AllPass_series
    AllPass_series = AllPass_series[~AllPass_series.isin(plasmid_sequences)]

    # Create DataFrame
    AllPass_df = pd.DataFrame(AllPass_series, columns=['Sequence Id'])

    # Save to file
    filename = Inputfile + "_viral_predictionsList.csv"
    AllPass_df.to_csv(os.path.join(OUT_DIR, filename), index=False)
    
else:  # Non-concentration type
    virsorter2_list1 = read_and_filter_virsorter2(os.path.join(Virsorterpath, "final-viral-score.tsv"), 'pass1')
    genomad_list1 = read_and_filter_genomad(os.path.join(Genomadpath), 'pass1')
    viralverify_list1 = read_and_filter_viralverify(os.path.join(Viralverifypath), 'pass1')

    # Merge Pass1 results for concentration type
    AllPass_series = merge_lists(virsorter2_list1, genomad_list1, viralverify_list1)
    # Read and filter plasmid data from Genomad data
    plasmid_sequences = read_plasmid_data(Genomadpath)

    # Read and filter non-viral sequences (Plasmid, Chromosome, Uncertain - plasmid or chromosomal) from ViralVerify data
    nonviral_sequences = read_and_filter_viralverify_nonviral(Viralverifypath)

    # Remove all sequences in plasmid_sequences and nonviral_sequences from AllPass_series
    AllPass_series = AllPass_series[~AllPass_series.isin(plasmid_sequences | nonviral_sequences)]

    # Create DataFrame
    AllPass_df = pd.DataFrame(AllPass_series, columns=['Sequence Id'])

    # Save to file
    filename = Inputfile + "_viral_predictionsList.csv"
    AllPass_df.to_csv(os.path.join(OUT_DIR, filename), index=False)
