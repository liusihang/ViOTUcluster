#!/usr/bin/env bash

# 定义 vOTU.Abundance.csv 的路径
ABUNDANCE_CSV="$OUTPUT_DIR/Summary/vOTU/vOTU.Abundance.csv"

# 检查 vOTU.Abundance.csv 是否已存在
if [ -f "$ABUNDANCE_CSV" ]; then
    echo "vOTU.Abundance.csv already exited, skip..."
else
    # 创建临时目录用于 TPM 计算
    echo "Creating temporary directory for TPM calculation..."
    mkdir -p "$OUTPUT_DIR/Summary/Viralcontigs/TPMTemp"
    if [ $? -ne 0 ]; then
        echo "Error: Failed to create temporary directory for TPM calculation."
        exit 1
    fi

    # 构建 BWA 索引
    echo "Building BWA index..."
    bwa index -p "$OUTPUT_DIR/Summary/Viralcontigs/TPMTemp/TempIndex" "$OUTPUT_DIR/Summary/vOTU/vOTU.fasta"
    if [ $? -ne 0 ]; then
        echo "Error: Failed to build BWA index."
        exit 1
    fi

    # 执行 Binning 分析
    for FILE in $FILES; do
        echo "Processing $FILE..."
        
        BASENAME=$(basename "$FILE" .fa)
        BASENAME=${BASENAME%.fasta}
        
        # 检查是否已经存在 coverage 文件，若存在则跳过
        if [ -f "$OUTPUT_DIR/Summary/Viralcontigs/Temp/${BASENAME}_coverage.tsv" ]; then
            echo "Skipping $BASENAME as coverage file already exists."
            continue
        fi
        
        # 检查并设置 Read1 和 Read2 的路径
        Read1=$(find "${RAW_SEQ_DIR}" -type f -name "${BASENAME}_R1*" | head -n 1)
        Read2=$(find "${RAW_SEQ_DIR}" -type f -name "${BASENAME}_R2*" | head -n 1)
        
        if [ -z "$Read1" ] || [ -z "$Read2" ]; then
            echo "Error: Read1 or Read2 files not found for $BASENAME."
            exit 1
        fi

        # 使用 BWA-MEM 进行比对
        echo "Aligning reads for $BASENAME..."
        bwa mem -t "${THREADS}" "$OUTPUT_DIR/Summary/Viralcontigs/TPMTemp/TempIndex" ${Read1} ${Read2} > "$OUTPUT_DIR/Summary/Viralcontigs/TPMTemp/${BASENAME}_gene.sam"
        if [ $? -ne 0 ]; then
            echo "Error: Failed to perform BWA alignment for $BASENAME."
            exit 1
        fi

        # 将 SAM 文件转换为 BAM 格式
        echo "Converting SAM to BAM..."
        samtools view -bS --threads "${THREADS}" "$OUTPUT_DIR/Summary/Viralcontigs/TPMTemp/${BASENAME}_gene.sam" > "$OUTPUT_DIR/Summary/Viralcontigs/TPMTemp/${BASENAME}_gene.bam"
        if [ $? -ne 0 ]; then
            echo "Error: Failed to convert SAM to BAM for $BASENAME."
            exit 1
        fi

        # 按照坐标对 BAM 文件进行排序
        echo "Sorting BAM file..."
        samtools sort "$OUTPUT_DIR/Summary/Viralcontigs/TPMTemp/${BASENAME}_gene.bam" -o "$OUTPUT_DIR/Summary/Viralcontigs/TPMTemp/${BASENAME}_sorted_gene.bam" --threads "${THREADS}"
        if [ $? -ne 0 ]; then
            echo "Error: Failed to sort BAM file for $BASENAME."
            exit 1
        fi

        # 为排序后的 BAM 文件生成索引
        echo "Indexing BAM file..."
        samtools index "$OUTPUT_DIR/Summary/Viralcontigs/TPMTemp/${BASENAME}_sorted_gene.bam"
        if [ $? -ne 0 ]; then
            echo "Error: Failed to generate BAM index for $BASENAME."
            exit 1
        fi

        # 计算覆盖度
        echo "Calculating coverage..."
        mkdir -p "$OUTPUT_DIR/Summary/Viralcontigs/TPMTemp/binsf"
        if [ $? -ne 0 ]; then
            echo "Error: Failed to create directory for coverage calculation."
            exit 1
        fi
        cp "$OUTPUT_DIR/Summary/vOTU/vOTU.fasta" "$OUTPUT_DIR/Summary/Viralcontigs/TPMTemp/binsf"
        checkm coverage -x fasta -m 20 -t "${THREADS}" "$OUTPUT_DIR/Summary/Viralcontigs/TPMTemp/binsf" "$OUTPUT_DIR/Summary/Viralcontigs/Temp/${BASENAME}_coverage.tsv" "$OUTPUT_DIR/Summary/Viralcontigs/TPMTemp/${BASENAME}_sorted_gene.bam"
        if [ $? -ne 0 ]; then
            echo "Error: Failed to calculate coverage for $BASENAME."
            exit 1
        fi

        # 清理中间文件
        rm -r "$OUTPUT_DIR/Summary/Viralcontigs/TPMTemp/binsf"
    done

    # 运行 TPM 计算 Python 脚本
    echo "Running TPM calculation..."
    python ${ScriptDir}/TPM_caculate.py "$OUTPUT_DIR/Summary/Viralcontigs/Temp" "$ABUNDANCE_CSV"
    if [ $? -ne 0 ]; then
        echo "Error: Failed to run TPM calculation."
        exit 1
    fi

    # 清理中间文件
    rm -r "$OUTPUT_DIR/Summary/Viralcontigs/TPMTemp"
    echo "TPM calculation completed successfully."
fi