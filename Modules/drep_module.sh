#!/usr/bin/env bash

# Stop execution and log errors on failure
set -e
trap 'echo "[❌] An error occurred in the dRep module. Exiting..."; exit 1;' ERR

echo "[📁] Creating required directories..."
mkdir -p "$OUTPUT_DIR/Summary/temp"
mkdir -p "$OUTPUT_DIR/Summary/Viralcontigs"
mkdir -p "$OUTPUT_DIR/Summary/dRepRes"

# Define the path for DrepBins.fasta
DREP_BINS_FASTA="$OUTPUT_DIR/Summary/dRepRes/DrepBins.fasta"

# Check if DrepBins.fasta already exists
if [ -f "$DREP_BINS_FASTA" ]; then
    echo "[✅] DrepBins.fasta already exists, skipping dRep and fasta concatenation steps."
else
    # Generate a list of genome files
    GENOME_LIST_FILE="${OUTPUT_DIR}/Summary/temp/genome_list.txt"
    find "${OUTPUT_DIR}/Summary/SeperateRes/bins" -name "*.fasta" > "$GENOME_LIST_FILE"
    echo "Genome list file generated at $GENOME_LIST_FILE"

    echo "[🔄] Starting dRep dereplication for bins..."
    dRep dereplicate "$OUTPUT_DIR/Summary/dRepRes" -g "$GENOME_LIST_FILE" --ignoreGenomeQuality --skip_plots -pa 0.8 -sa 0.95 -nc 0.85 -comW 0 -conW 0 -strW 0 -N50W 0 -sizeW 1 -centW 0 -l 3000
    echo "[✅] dRep dereplication completed."

    echo "Concatenating dereplicated fasta sequences..."
    python "${ScriptDir}/concat_fasta_sequences.py" "$OUTPUT_DIR/Summary/dRepRes/dereplicated_genomes" "$DREP_BINS_FASTA"
    echo "[✅] Fasta concatenation completed."
fi

# Define the path for DrepViralcontigs.fasta
DREP_VIRAL_FASTA="$OUTPUT_DIR/Summary/dRepRes/DrepViralcontigs.fasta"

# Check if DrepViralcontigs.fasta already exists
if [ -f "$DREP_VIRAL_FASTA" ]; then
    echo "[✅] DrepViralcontigs.fasta already exists, skipping dRep and clustering steps."
else
    echo "Starting dRep for unbined contigs..."
    cat "$OUTPUT_DIR/Summary/SeperateRes/unbined/"*_unbined.fasta > "$OUTPUT_DIR/Summary/temp/merged_sequences.fasta"
    echo "Contigs merging completed."

    newDir="$OUTPUT_DIR/Summary/temp"
    
    echo "[🔄] Filtering sequences shorter than 5000bp..."
    awk 'BEGIN {RS=">";FS="\n"} NR>1 {seq=""; for(i=2;i<=NF;i++) seq=seq $i; if(length(seq)>=5000) print ">" $1 "\n" seq}' \
        "${newDir}/merged_sequences.fasta" > "${newDir}/merged_sequences_filtered.fasta"
    cp "${newDir}/merged_sequences_filtered.fasta" "${newDir}/merged_sequences.fasta"
    echo "[✅] Filtering completed."

    echo "[🔄] Building BLAST database and running clustering..."
    makeblastdb -in "${newDir}/merged_sequences.fasta" -dbtype nucl -out "${newDir}/temp_db"

    blastn -query "${newDir}/merged_sequences.fasta" -db "${newDir}/temp_db" -outfmt "6 std qlen slen" \
        -max_target_seqs 10000 -out "${newDir}/merged_sequences_blast.tsv" -num_threads "${THREADS}"

    python "${ScriptDir}/anicalc.py" -i "${newDir}/merged_sequences_blast.tsv" -o "${newDir}/merged_sequences_ani.tsv"

    python "${ScriptDir}/aniclust.py" --fna "${newDir}/merged_sequences.fasta" --ani "${newDir}/merged_sequences_ani.tsv" \
        --out "${newDir}/merged_sequences_clusters.tsv" --min_ani 95 --min_tcov 85 --min_qcov 0

    echo "[🧹] Cleaning up temporary files..."
    rm -f "${newDir}/temp_db.*"
    rm -f "${newDir}/merged_sequences_blast.tsv"

    echo "Selecting representative sequences from clusters..."
    python "${ScriptDir}/SelectCluster.py" "${newDir}/merged_sequences.fasta" "${newDir}/merged_sequences_clusters.tsv" "$DREP_VIRAL_FASTA"
    echo "[✅] dRep and clustering for unbined contigs completed."
fi