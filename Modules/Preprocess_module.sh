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

# 创建输出目录并检查是否成功
Cleanreads="${OUTPUT_DIR}/Cleanreads"
ASSEMBLY_DIR="${OUTPUT_DIR}/Assembly"
CONTIGS_DIR="${OUTPUT_DIR}/Contigs"

mkdir -p $Cleanreads $ASSEMBLY_DIR $CONTIGS_DIR
if [ $? -ne 0 ]; then
    echo "Error: Failed to create output directories."
    exit 1
fi

# 遍历所有可能的_R1文件（包括fastq.gz, fq.gz, fastq, fq格式）
for R1_FILE in "${INPUT_DIR}"/*_R1.{fq.gz,fastq.gz,fq,fastq}; do
    # 跳过不存在的模式
    [ -e "$R1_FILE" ] || continue
    
    # 自动查找对应的_R2文件
    R2_FILE="${R1_FILE/_R1/_R2}"

    # 确保R2文件存在
    if [ -f "$R2_FILE" ]; then
        # 获取文件的前缀名字
        PREFIX=$(basename "$R1_FILE" | sed 's/_R1.*//')

        # 检查fastp输出文件是否已存在
        if [ -f "${Cleanreads}/${PREFIX}_R1.fq.gz" ] && [ -f "${Cleanreads}/${PREFIX}_R2.fq.gz" ]; then
            echo "Fastp already completed for $PREFIX, skipping..."
        else
            # 使用fastp进行处理并检查是否成功
            echo "Running fastp for $PREFIX..."
            fastp -i "$R1_FILE" -I "$R2_FILE" -o "${Cleanreads}/${PREFIX}_R1.fq.gz" -O "${Cleanreads}/${PREFIX}_R2.fq.gz"
            if [ $? -ne 0 ]; then
                echo "Error: fastp processing failed for $PREFIX."
                exit 1
            fi
        fi

        # 选择组装软件并检查目标文件是否已存在
        if [ "$ASSEMBLY_SOFTWARE" == "megahit" ]; then
            if [ -f "${CONTIGS_DIR}/${PREFIX}.fa" ]; then
                echo "Megahit assembly already completed for $PREFIX, skipping..."
            else
                echo "Running megahit for $PREFIX..."
                megahit -1 "${Cleanreads}/${PREFIX}_R1.fq.gz" -2 "${Cleanreads}/${PREFIX}_R2.fq.gz" -o "${ASSEMBLY_DIR}/${PREFIX}_megahit_out"
                if [ $? -ne 0 ]; then
                    echo "Error: Megahit assembly failed for $PREFIX."
                    exit 1
                fi

                # 提取组装结果
                cp "${ASSEMBLY_DIR}/${PREFIX}_megahit_out/final.contigs.fa" "${CONTIGS_DIR}/${PREFIX}.fa"
                if [ $? -ne 0 ]; then
                    echo "Error: Failed to copy megahit contigs for $PREFIX."
                    exit 1
                fi
            fi

        elif [ "$ASSEMBLY_SOFTWARE" == "metaspades" ]; then
            if [ -f "${CONTIGS_DIR}/${PREFIX}.fa" ]; then
                echo "Metaspades assembly already completed for $PREFIX, skipping..."
            else
                echo "Running metaspades for $PREFIX..."
                metaspades.py -1 "${Cleanreads}/${PREFIX}_R1.fq.gz" -2 "${Cleanreads}/${PREFIX}_R2.fq.gz" -o "${ASSEMBLY_DIR}/${PREFIX}_spades_out"
                if [ $? -ne 0 ]; then
                    echo "Error: Metaspades assembly failed for $PREFIX."
                    exit 1
                fi

                # 提取组装结果
                cp "${ASSEMBLY_DIR}/${PREFIX}_spades_out/scaffolds.fasta" "${CONTIGS_DIR}/${PREFIX}.fa"
                if [ $? -ne 0 ]; then
                    echo "Error: Failed to copy metaspades contigs for $PREFIX."
                    exit 1
                fi
            fi
        fi
    else
        echo "Paired file for $R1_FILE not found. Skipping..."
    fi
done

echo "Processing completed. Contigs are saved in $CONTIGS_DIR."