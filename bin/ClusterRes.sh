#!/bin/bash

# 给定文件夹路径
inputDir="/media/em/student/LSH/Virus/vOTUclusterExperiments/SummaryFasta"

# 输出总目录
outputDir="/media/em/student/LSH/Virus/vOTUclusterExperiments/SummaryFasta"
mkdir -p "${outputDir}/FinalRes"

# 循环处理每个FASTA文件
for fastaFile in "${inputDir}"/*.fasta; do
    # 获取不包含扩展名的文件名
    baseName=$(basename "$fastaFile" .fasta)

    # 为每个文件创建输出子文件夹
    subdir="${outputDir}/${baseName}/DeRep"

    # 合并聚类结果到总目录下的FinalRes
    python SelectCluster.py "$fastaFile" "${subdir}/merged_sequences_clusters.tsv" "${outputDir}/FinalRes/${baseName}_vOTU.fa"

    # 删除临时文件
    rm -f "${subdir}/temp_db.*"
    rm -f "${subdir}/merged_sequences_blast.tsv"
done


