#!/usr/bin/env bash

# 在发生错误时，停止执行并记录日志
set -e
trap 'echo "An error occurred in the dRep module. Exiting..."; exit 1;' ERR

# dRep for bins
echo "Starting dRep for bins..."
mkdir -p "$OUTPUT_DIR/Summary/Viralcontigs"
mkdir -p "$OUTPUT_DIR/Summary/dRepRes"
dRep dereplicate "$OUTPUT_DIR/Summary/dRepRes" -g ${OUTPUT_DIR}/Summary/bins/*.fasta --ignoreGenomeQuality -pa 0.8 -sa 0.95 -nc 0.85 -comW 0 -conW 0 -strW 0 -N50W 0 -sizeW 1 -centW 0
echo "dRep for bins completed."

# Concat fasta sequences
echo "Concatenating fasta sequences..."
python "${ScriptDir}/concat_fasta_sequences.py" "$OUTPUT_DIR/Summary/dRepRes/dereplicated_genomes" "$OUTPUT_DIR/Summary/Viralcontigs/DrepBins.fasta"
echo "Fasta concatenation completed."

# dRep for contigs
echo "Starting dRep for unbined contigs..."
mkdir -p "$OUTPUT_DIR/Summary/temp"
cat "$OUTPUT_DIR/Summary/unbined/"*.fasta > "$OUTPUT_DIR/Summary/temp/merged_sequences.fasta"
echo "Contigs merging completed."

# Cluster
newDir="$OUTPUT_DIR/Summary/temp"
echo "Clustering..."

# First, create a blast+ database
makeblastdb -in "${newDir}/merged_sequences.fasta" -dbtype nucl -out "${newDir}/temp_db"

# Perform all-vs-all blastn of sequences
blastn -query "${newDir}/merged_sequences.fasta" -db "${newDir}/temp_db" -outfmt "6 std qlen slen" -max_target_seqs 10000 -out "${newDir}/merged_sequences_blast.tsv" -num_threads "${THREADS}"

# Calculate pairwise ANI
python "${ScriptDir}/anicalc.py" -i "${newDir}/merged_sequences_blast.tsv" -o "${newDir}/merged_sequences_ani.tsv"

# Perform UCLUST-like clustering using the MIUVIG recommended parameters (95% ANI + 85% AF)
python "${ScriptDir}/aniclust.py" --fna "${newDir}/merged_sequences.fasta" --ani "${newDir}/merged_sequences_ani.tsv" --out "${newDir}/merged_sequences_clusters.tsv" --min_ani 95 --min_tcov 85 --min_qcov 0

# Delete temporary files
echo "Cleaning up temporary files..."
rm -f "${newDir}/temp_db.*"
rm -f "${newDir}/merged_sequences_blast.tsv"

# Merge cluster results
echo "Merging cluster results..."
python "${ScriptDir}/SelectCluster.py" "${newDir}/merged_sequences.fasta" "${newDir}/merged_sequences_clusters.tsv" "$OUTPUT_DIR/Summary/Viralcontigs/DrepViralcontigs.fasta"

# Merge all
echo "Merging final sequences..."
python "${ScriptDir}/Rename.py" -i "$OUTPUT_DIR/Summary/Viralcontigs/DrepViralcontigs.fasta"
cat "$OUTPUT_DIR/Summary/Viralcontigs/DrepViralcontigs.fasta" "$OUTPUT_DIR/Summary/Viralcontigs/DrepBins.fasta" > "$OUTPUT_DIR/Summary/Viralcontigs/vOTU.fasta"
checkv end_to_end "$OUTPUT_DIR/Summary/Viralcontigs/vOTU.fasta" "$OUTPUT_DIR/Summary/Viralcontigs/vOTU_CheckRes" -t "${THREADS}" -d "$DATABASE/checkv-db-v1.5"

echo "Combined fasta files and quality summaries completed."