#!/usr/bin/env python3
import pandas as pd
import sys
import os

def format_taxonomy(input_csv, output_file):
    """
    Format taxonomy data by splitting lineage into hierarchical levels and handling missing values.

    Main modifications:
    1. Remove the 'lineage' column from the output file.
    2. For missing taxonomy levels, fill them using the format of the previous level + "__unclassified".
       - From Domain to Family levels: If missing, fill with the previous level + "__unclassified".
       - For Genus and Species levels:
         * If Genus is missing, fill with Family + "__unclassified". If Species is also missing, fill with Family + "__unclassified".
         * If Genus exists but Species is missing, fill with Genus + "__unclassified".
       
    Parameters:
    - input_csv (str): Path to the input CSV file.
    - output_file (str): Path to save the formatted CSV file.
    """
    df = pd.read_csv(input_csv, sep='\t')

    # Select the 'seq_name' and 'lineage' columns
    df_phyloseq = df[['seq_name', 'lineage']].copy()

    # Define taxonomy levels from Domain to Species
    taxonomy_levels = ['Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
    df_phyloseq[taxonomy_levels] = df_phyloseq['lineage'].str.split(';', expand=True)

    # Fill missing values for Domain to Family levels (index 1 to 4)
    for idx, row in df_phyloseq.iterrows():
        for i in range(1, len(taxonomy_levels) - 2):  # i = 1 to 4, for Phylum, Class, Order, Family
            current_level = taxonomy_levels[i]
            previous_level = taxonomy_levels[i - 1]
            if pd.isna(row[current_level]) or row[current_level] == "":
                df_phyloseq.at[idx, current_level] = row[previous_level] + '__unclassified'
        
        # Special handling for Genus and Species to ensure consistent filling
        family = df_phyloseq.at[idx, 'Family']
        genus = row['Genus']
        species = row['Species']
        
        if pd.isna(genus) or genus == "":
            # If Genus is missing, fill with Family + '__unclassified'
            df_phyloseq.at[idx, 'Genus'] = family + '__unclassified'
            # If Species is also missing, fill with Family + '__unclassified'
            if pd.isna(species) or species == "":
                df_phyloseq.at[idx, 'Species'] = family + '__unclassified'
        else:
            # If Genus exists but Species is missing, fill with Genus + '__unclassified'
            if pd.isna(species) or species == "":
                df_phyloseq.at[idx, 'Species'] = df_phyloseq.at[idx, 'Genus'] + '__unclassified'

    # Remove the unnecessary 'lineage' column
    df_phyloseq.drop(columns=['lineage'], inplace=True)
    
    # Rename 'seq_name' to 'OTU'
    df_phyloseq = df_phyloseq.rename(columns={'seq_name': 'OTU'})

    # Save the formatted DataFrame to a CSV file
    df_phyloseq.to_csv(output_file, index=False)
    print(f"File saved as: {output_file}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script_name.py <input_csv_file> <output_file>")
        sys.exit(1)

    input_csv_file = sys.argv[1]
    output_file = sys.argv[2]

    format_taxonomy(input_csv_file, output_file)
