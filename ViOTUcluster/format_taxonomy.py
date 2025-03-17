import pandas as pd
import sys
import os

def format_taxonomy(input_csv, output_file):
    """
    Format taxonomy data by splitting lineage into hierarchical levels and handling missing values.

    Parameters:
    - input_csv (str): Path to the input CSV file.
    - output_file (str): Path to save the formatted CSV file.
    """
    df = pd.read_csv(input_csv, sep='\t')

    # Select 'seq_name' and 'lineage' columns
    df_phyloseq = df[['seq_name', 'lineage']].copy()

    # Define taxonomy levels from Domain to Species
    taxonomy_levels = ['Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
    df_phyloseq[taxonomy_levels] = df_phyloseq['lineage'].str.split(';', expand=True)

    # Fill missing taxonomy levels with the previous level + '_unclassified'
    for idx, row in df_phyloseq.iterrows():
        for i in range(1, len(taxonomy_levels)):
            current_level = taxonomy_levels[i]
            previous_level = taxonomy_levels[i - 1]
            if pd.isna(row[current_level]):
                for j in range(i, len(taxonomy_levels)):
                    df_phyloseq.at[idx, taxonomy_levels[j]] = df_phyloseq.at[idx, taxonomy_levels[j - 1]] + '_unclassified'
                break

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
