#!/usr/bin/env python3
import pandas as pd
import sys
import os

LEVELS = ['Realm', 'Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
FILL_SUFFIX = '__unclassified'

# 兼容别名到标准层级（小写比较）
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
    支持两种格式：
    - 位置式: 以 ';' 分隔，依次对应 LEVELS
    - 键值式: 'level=name' 或 'level:name'，level 大小写不敏感
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
        # 位置式：按 LEVELS 顺序填入
        for i, lvl in enumerate(LEVELS):
            if i < len(tokens):
                val = tokens[i].strip()
                result[lvl] = val if val else None

    return result

def fill_row_hierarchy(row: pd.Series) -> pd.Series:
    """
    按规则填补缺失层级。
    """
    # Realm 顶层
    if pd.isna(row['Realm']) or row['Realm'] == '':
        row['Realm'] = 'unclassified'

    # Kingdom..Family：用上一层 + __unclassified
    for i in range(1, LEVELS.index('Family') + 1):
        cur = LEVELS[i]
        prev = LEVELS[i - 1]
        if pd.isna(row[cur]) or row[cur] == '':
            base = row[prev] if pd.notna(row[prev]) and row[prev] != '' else 'unclassified'
            row[cur] = f"{base}{FILL_SUFFIX}"

    # Genus
    if pd.isna(row['Genus']) or row['Genus'] == '':
        family = row['Family'] if pd.notna(row['Family']) and row['Family'] != '' else 'unclassified'
        row['Genus'] = f"{family}{FILL_SUFFIX}"

    # Species
    if pd.isna(row['Species']) or row['Species'] == '':
        genus = row['Genus'] if pd.notna(row['Genus']) and row['Genus'] != '' else None
        family = row['Family'] if pd.notna(row['Family']) and row['Family'] != '' else 'unclassified'
        if genus:
            row['Species'] = f"{genus}{FILL_SUFFIX}"
        else:
            row['Species'] = f"{family}{FILL_SUFFIX}"

    return row

def format_taxonomy(input_csv, output_file):
    """
    读取包含 seq_name 与 lineage 的表（输入使用制表符分隔），
    拆分并标准化为 Realm→Species 八列；去掉原 lineage 列；导出为 TSV。
    """
    df = pd.read_csv(input_csv, sep='\t')

    if 'seq_name' not in df.columns or 'lineage' not in df.columns:
        raise ValueError("Input file must contain 'seq_name' and 'lineage' columns.")

    # 逐行解析 lineage
    parsed = df['lineage'].apply(parse_lineage)
    tax_df = pd.DataFrame(list(parsed.values), columns=LEVELS)

    # 合并并按规则补全
    out = pd.concat([df[['seq_name']].reset_index(drop=True), tax_df], axis=1)
    out = out.apply(fill_row_hierarchy, axis=1)

    # 重命名并删除 lineage
    out = out.rename(columns={'seq_name': 'OTU'})：
    other_cols = [c for c in df.columns if c not in ('seq_name', 'lineage')]
    if other_cols:
        out = pd.concat([out, df[other_cols].reset_index(drop=True)], axis=1)

    # 导出为 TSV（与输入保持一致）
    out.to_csv(output_file, index=False, sep='\t')
    print(f"File saved as: {output_file}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script_name.py <input_tsv_file> <output_tsv_file>")
        sys.exit(1)
    format_taxonomy(sys.argv[1], sys.argv[2])