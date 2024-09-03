#!/usr/bin/env bash

# 检查参数是否正确
if [ $# -lt 3 ]; then
    echo "Usage: $0 <input_directory> <assembly_software> <output_directory>"
    echo "assembly_software options: megahit, metaspades"
    exit 1
fi

# 获取输入参数
INPUT_DIR=$1
ASSEMBLY_SOFTWARE=$2
OUTPUT_DIR=$3

# 验证组装软件选项
if [[ "$ASSEMBLY_SOFTWARE" != "megahit" && "$ASSEMBLY_SOFTWARE" != "metaspades" ]]; then
    echo "Error: Invalid assembly software. Please choose either 'megahit' or 'metaspades'."
    exit 1
fi

# 创建输出目录
Cleanreads="${OUTPUT_DIR}/Cleanreads"
ASSEMBLY_DIR="${OUTPUT_DIR}/Assembly"
CONTIGS_DIR="${OUTPUT_DIR}/Contigs"
mkdir -p $Cleanreads
mkdir -p $ASSEMBLY_DIR
mkdir -p $CONTIGS_DIR

# 遍历所有配对的_R1.fq.gz和_R2.fq.gz文件
for R1_FILE in "${INPUT_DIR}"/*_R1.fq.gz; do
    R2_FILE=${R1_FILE/_R1.fq.gz/_R2.fq.gz}

    # 确保R2文件存在
    if [ -f "$R2_FILE" ]; then
        # 获取文件的前缀名字
        PREFIX=$(basename $R1_FILE _R1.fq.gz)

        # 使用fastp进行处理
        echo "Running fastp for $PREFIX..."
        fastp -i $R1_FILE -I $R2_FILE -o ${Cleanreads}/${PREFIX}_R1.fq.gz -O ${Cleanreads}/${PREFIX}_R2.fq.gz

        # 选择组装软件并进行组装
        if [ "$ASSEMBLY_SOFTWARE" == "megahit" ]; then
            echo "Running megahit for $PREFIX..."
            megahit -1 ${Cleanreads}/${PREFIX}_R1.fq.gz -2 ${Cleanreads}/${PREFIX}_R2.fq.gz -o ${ASSEMBLY_DIR}/${PREFIX}_megahit_out
            # 提取组装结果
            cp ${ASSEMBLY_DIR}/${PREFIX}_megahit_out/final.contigs.fa ${CONTIGS_DIR}/${PREFIX}.fa
        elif [ "$ASSEMBLY_SOFTWARE" == "metaspades" ]; then
            echo "Running metaspades for $PREFIX..."
            metaspades.py -1 ${Cleanreads}/${PREFIX}_R1.fq.gz -2 ${Cleanreads}/${PREFIX}_R2.fq.gz -o ${ASSEMBLY_DIR}/${PREFIX}_spades_out
            # 提取组装结果
            cp ${ASSEMBLY_DIR}/${PREFIX}_spades_out/scaffolds.fasta ${CONTIGS_DIR}/${PREFIX}.fa
        fi
    else
        echo "Paired file for $R1_FILE not found. Skipping..."
    fi
done

echo "Processing completed. Contigs are saved in $CONTIGS_DIR."