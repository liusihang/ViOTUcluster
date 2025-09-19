#!/usr/bin/env python3
import pandas as pd
import sys

LEVELS = ['Realm', 'Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
FILL_SUFFIX = '__unclassified'

ALIASES = {
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

def parse_lineage(lineage_str: str) -> dict:
    """
    将谱系字符串解析为 {level: name} 字典。
    """
    result = {lvl: None for lvl in LEVELS}
    if pd.isna(lineage_str):
        return result

    s = str(lineage_str).strip()
    if not s:
        return result

    tokens = [t.strip() for t in s.split(';') if t.strip()]
    is_kv = any(('=' in t) or (':' in t) for t in tokens)

    if is_kv:
        for t in tokens:
            if '=' in t:
                k, v = t.split('=', 1)
            elif ':' in t:
                k, v = t.split(':', 1)
            else:
                continue
            k = k.strip().lower()
            v = v.strip()
            if not v:
                continue
            if k in ALIASES:
                result[ALIASES[k]] = v
    else:
        for i, lvl in enumerate(LEVELS):
            if i < len(tokens):
                val = tokens[i].strip()
                result[lvl] = val if val else None

    return result

def fill_row_hierarchy(row: pd.Series) -> pd.Series:
    """
    缺失填补规则：
    - Realm 缺失→ 'unclassified'
    - Kingdom..Family：用上一层 + '__unclassified'
    - Genus 缺失→ Family + '__unclassified'
    - Species 缺失→ 若 Genus 有值，用 Genus + '__unclassified'，否则用 Family + '__unclassified'
    """
    if pd.isna(row['Realm']) or row['Realm'] == '':
        row['Realm'] = 'unclassified'

    for i in range(1, LEVELS.index('Family') + 1):
        cur = LEVELS[i]
        prev = LEVELS[i - 1]
        if pd.isna(row[cur]) or row[cur] == '':
            base = row[prev] if pd.notna(row[prev]) and row[prev] != '' else 'unclassified'
            row[cur] = f"{base}{FILL_SUFFIX}"

    if pd.isna(row['Genus']) or row['Genus'] == '':
        family = row['Family'] if pd.notna(row['Family']) and row['Family'] != '' else 'unclassified'
        row['Genus'] = f"{family}{FILL_SUFFIX}"

    if pd.isna(row['Species']) or row['Species'] == '':
        genus = row['Genus'] if pd.notna(row['Genus']) and row['Genus'] != '' else None
        family = row['Family'] if pd.notna(row['Family']) and row['Family'] != '' else 'unclassified'
        row['Species'] = f"{genus}{FILL_SUFFIX}" if genus else f"{family}{FILL_SUFFIX}"

    return row

def format_taxonomy(input_csv, output_file):
    """
    输入：含 `seq_name` 与 `lineage` 的 TSV。
    输出：仅包含 9 列（OTU + 8 层级）的 TSV：
        OTU, Realm, Kingdom, Phylum, Class, Order, Family, Genus, Species。
    """
    df = pd.read_csv(input_csv, sep='\t')

    if 'seq_name' not in df.columns or 'lineage' not in df.columns:
        raise ValueError("Input file must contain 'seq_name' and 'lineage' columns.")

    parsed = df['lineage'].apply(parse_lineage)
    tax_df = pd.DataFrame(list(parsed.values), columns=LEVELS)

    out = pd.concat([df[['seq_name']].reset_index(drop=True), tax_df], axis=1)
    out = out.apply(fill_row_hierarchy, axis=1)

    out = out.rename(columns={'seq_name': 'OTU'})
    out = out[['OTU'] + LEVELS]  # 强制只保留 9 列

    out.to_csv(output_file, index=False, sep='\t')
    print(f"File saved as: {output_file}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script_name.py <input_tsv_file> <output_tsv_file>")
        sys.exit(1)
    format_taxonomy(sys.argv[1], sys.argv[2])