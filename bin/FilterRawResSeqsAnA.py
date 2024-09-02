import csv
import os
import sys
from Bio import SeqIO

# 接受命令行参数
fasta = sys.argv[1]
Inputname = sys.argv[2]
OUT_DIR = sys.argv[3]

# 从CSV文件中读取ID列表
csv_values = set()
csv_path = os.path.join(OUT_DIR, f"{Inputname}.csv")
with open(csv_path, 'r', encoding='utf-8') as csvfile:
    reader = csv.reader(csvfile)
    for row in reader:
        if row:  # 检查行是否不为空
            # 处理CSV中的ID
            csv_id = row[0].split('|')[0]  # 取'|'前面的部分
            csv_values.add(csv_id)

# 读取FASTA文件并找到匹配的序列
matches = []
for record in SeqIO.parse(fasta, "fasta"):
    # 全字匹配序列ID
    if record.id in csv_values:
        matches.append(record)

# 如果找到匹配的序列，保存到新的FASTA文件
if matches:
    output_fasta_filename = os.path.join(OUT_DIR, f"{Inputname}_filtered.fasta")
    SeqIO.write(matches, output_fasta_filename, "fasta")

print("所有文件处理完成。")