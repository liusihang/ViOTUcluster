#!/usr/bin/env python3
import sys
from typing import Dict, Optional

import pandas as pd

LEVELS = ['Realm', 'Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
FILL_SUFFIX = '__unclassified'

ALIASES: Dict[str, str] = {
    'realm': 'Realm',
    'kingdom': 'Kingdom',
    'phylum': 'Phylum',
    'class': 'Class',
    'order': 'Order',
    'family': 'Family',
    'genus': 'Genus',
    'species': 'Species',
    'domain': 'Realm',
    'superkingdom': 'Realm',
}


def parse_lineage(lineage_str: str) -> Dict[str, Optional[str]]:
    """
    Parse a lineage string into a dictionary keyed by the standard LEVELS list.
    """
    result = {lvl: None for lvl in LEVELS}
    if pd.isna(lineage_str):
        return result

    s = str(lineage_str).strip()
    if not s:
        return result

    tokens = [token.strip() for token in s.split(';') if token.strip()]
    is_kv = any(('=' in token) or (':' in token) for token in tokens)

    if is_kv:
        for token in tokens:
            if '=' in token:
                key, value = token.split('=', 1)
            elif ':' in token:
                key, value = token.split(':', 1)
            else:
                continue
            key = key.strip().lower()
            value = value.strip()
            if not value:
                continue
            if key in ALIASES:
                result[ALIASES[key]] = value
    else:
        for i, level in enumerate(LEVELS):
            if i < len(tokens):
                value = tokens[i].strip()
                result[level] = value if value else None

    return result


def fill_row_hierarchy(row: pd.Series) -> pd.Series:
    """
    Fill missing levels by propagating the closest known parent and appending
    the default suffix, matching the previous behavior.
    """
    if pd.isna(row['Realm']) or row['Realm'] == '':
        row['Realm'] = 'unclassified'

    for i in range(1, LEVELS.index('Family') + 1):
        current = LEVELS[i]
        previous = LEVELS[i - 1]
        if pd.isna(row[current]) or row[current] == '':
            base = row[previous] if pd.notna(row[previous]) and row[previous] != '' else 'unclassified'
            row[current] = f"{base}{FILL_SUFFIX}"

    if pd.isna(row['Genus']) or row['Genus'] == '':
        family = row['Family'] if pd.notna(row['Family']) and row['Family'] != '' else 'unclassified'
        row['Genus'] = f"{family}{FILL_SUFFIX}"

    if pd.isna(row['Species']) or row['Species'] == '':
        genus = row['Genus'] if pd.notna(row['Genus']) and row['Genus'] != '' else None
        family = row['Family'] if pd.notna(row['Family']) and row['Family'] != '' else 'unclassified'
        row['Species'] = f"{genus}{FILL_SUFFIX}" if genus else f"{family}{FILL_SUFFIX}"

    return row


def format_taxonomy(input_csv: str, output_file: str, abundance_csv: Optional[str] = None) -> None:
    """
    Transform the genomad taxonomy table into OTU + 8 ranked columns and,
    optionally, ensure that every OTU present in the abundance table exists in
    the taxonomy file (adding Unclassified_virus placeholders when needed).
    """
    df = pd.read_csv(input_csv, sep='\t')

    if 'seq_name' not in df.columns or 'lineage' not in df.columns:
        raise ValueError("Input file must contain 'seq_name' and 'lineage' columns.")

    parsed = df['lineage'].apply(parse_lineage)
    tax_df = pd.DataFrame(list(parsed.values), columns=LEVELS)

    out = pd.concat([df[['seq_name']].reset_index(drop=True), tax_df], axis=1)
    out = out.apply(fill_row_hierarchy, axis=1)

    out = out.rename(columns={'seq_name': 'OTU'})
    out = out[['OTU'] + LEVELS]

    if abundance_csv:
        try:
            abundance_df = pd.read_csv(abundance_csv)
            abundance_ids = (
                abundance_df.iloc[:, 0]
                .dropna()
                .astype(str)
                .tolist()
            )
        except Exception as exc:
            print(f"[WARN] Failed to read abundance file {abundance_csv}: {exc}", file=sys.stderr)
            abundance_ids = []

        if abundance_ids:
            taxonomy_ids = set(out['OTU'].astype(str))
            missing_ids = [otu_id for otu_id in abundance_ids if otu_id not in taxonomy_ids]

            if missing_ids:
                filler = pd.DataFrame({'OTU': missing_ids})
                for level in LEVELS:
                    filler[level] = 'Unclassified_virus'
                out = pd.concat([out, filler], ignore_index=True)
                print(f"[INFO] Added {len(missing_ids)} OTUs missing from taxonomy.")

    out.to_csv(output_file, index=False, sep='\t')
    print(f"File saved as: {output_file}")


if __name__ == "__main__":
    if len(sys.argv) not in (3, 4):
        print("Usage: python script_name.py <input_tsv_file> <output_tsv_file> [abundance_csv]")
        sys.exit(1)

    abundance_file = sys.argv[3] if len(sys.argv) == 4 else None
    format_taxonomy(sys.argv[1], sys.argv[2], abundance_file)
