#!/usr/bin/env bash

# Stop execution and log errors on failure
set -e
trap 'echo "[❌] An error occurred in the dRep module. Exiting..."; exit 1;' ERR

DISABLE_BINNING=${DISABLE_BINNING:-false}

echo "[📁] Creating required directories..."
mkdir -p "$OUTPUT_DIR/Summary/temp"
mkdir -p "$OUTPUT_DIR/Summary/Viralcontigs"
mkdir -p "$OUTPUT_DIR/Summary/dRepRes"

# Define the path for DrepBins.fasta
DREP_BINS_FASTA="$OUTPUT_DIR/Summary/dRepRes/DrepBins.fasta"

if [ "$DISABLE_BINNING" = true ]; then
    echo "[ℹ️] Binning disabled upstream; creating empty DrepBins.fasta placeholder."
    : > "$DREP_BINS_FASTA"
elif [ -f "$DREP_BINS_FASTA" ]; then
    echo "[✅] DrepBins.fasta already exists, skipping dRep and fasta concatenation steps."
else
    shopt -s nullglob
    BIN_FASTA_FILES=("$OUTPUT_DIR/Summary/SeperateRes/bins/"*.fasta)
    shopt -u nullglob

    if [ "${#BIN_FASTA_FILES[@]}" -eq 0 ]; then
        echo "[⚠️] No bin FASTA files detected; creating empty DrepBins.fasta."
        : > "$DREP_BINS_FASTA"
    else
        GENOME_LIST_FILE="${OUTPUT_DIR}/Summary/temp/genome_list.txt"
        printf "%s\n" "${BIN_FASTA_FILES[@]}" > "$GENOME_LIST_FILE"
        echo "Genome list file generated at $GENOME_LIST_FILE"

        echo "[🔄] Starting dRep dereplication for bins..."
        dRep dereplicate "$OUTPUT_DIR/Summary/dRepRes" -g "$GENOME_LIST_FILE" --ignoreGenomeQuality --skip_plots -pa 0.8 -sa 0.95 -nc 0.85 -comW 0 -conW 0 -strW 0 -N50W 0 -sizeW 1 -centW 0 -l 3000
        echo "[✅] dRep dereplication completed."

        echo "Concatenating dereplicated fasta sequences..."
        python "${ScriptDir}/concat_fasta_sequences.py" "$OUTPUT_DIR/Summary/dRepRes/dereplicated_genomes" "$DREP_BINS_FASTA"
        echo "[✅] Fasta concatenation completed."
    fi
fi

# Define the path for DrepViralcontigs.fasta
DREP_VIRAL_FASTA="$OUTPUT_DIR/Summary/dRepRes/DrepViralcontigs.fasta"

# Check if DrepViralcontigs.fasta already exists
if [ -f "$DREP_VIRAL_FASTA" ]; then
    echo "[✅] DrepViralcontigs.fasta already exists, skipping dRep and clustering steps."
else
    echo "Starting dRep for unbined contigs..."
    shopt -s nullglob
    UNBINNED_FILES=("$OUTPUT_DIR/Summary/SeperateRes/unbined/"*_unbined.fasta)
    shopt -u nullglob

    if [ "${#UNBINNED_FILES[@]}" -eq 0 ]; then
        echo "[⚠️] No unbinned contig FASTA files found; creating empty DrepViralcontigs.fasta."
        : > "$DREP_VIRAL_FASTA"
    else
        cat "${UNBINNED_FILES[@]}" > "$OUTPUT_DIR/Summary/temp/merged_sequences.fasta"
        echo "Contigs merging completed."

        newDir="$OUTPUT_DIR/Summary/temp"
        rm -f "${newDir}/Done" # Ensure script runs even if Done file exists from previous run
        echo "[🔄] Filtering sequences shorter than ${MIN_LENGTH}bp"
        python "${ScriptDir}/filter_contigs.py" "${MIN_LENGTH}" "${newDir}/merged_sequences.fasta" "$newDir"
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
fi
