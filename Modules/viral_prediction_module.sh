#!/usr/bin/env bash

# Set error handling mechanism: if an error occurs, the script will stop executing
set -e
trap 'echo "An error occurred. Exiting..."; exit 1;' ERR

# Make sure BASE_CONDA_PREFIX is correctly set
BASE_CONDA_PREFIX=$(conda info --base)
conda_sh="$BASE_CONDA_PREFIX/etc/profile.d/conda.sh"

# Export necessary variables and functions
export OUTPUT_DIR DATABASE Group CONCENTRATION_TYPE THREADS_PER_FILE FILES

# Function to run viralverify
run_viralverify() {
  local FILE=$1
  local BASENAME=$(basename "$FILE")
  local OUT_DIR="$OUTPUT_DIR/SeprateFile/${BASENAME%.*}/RoughViralPrediction/viralverify"
  mkdir -p "$OUT_DIR"

  if [ ! -f "$OUT_DIR/${BASENAME}_result_table.csv" ]; then
    echo "Running viralverify prediction for $FILE..."
    viralverify -f "$FILE" -o "$OUT_DIR" --hmm "$DATABASE/ViralVerify/nbc_hmms.hmm" \
    -t "$THREADS_PER_FILE" > "$OUT_DIR/viralverify.log" 2>&1
  else
    echo "viralverify prediction already completed for $FILE, skipping..."
  fi
}

# Function to run VirSorter2
run_virsorter() {
  local FILE=$1
  local BASENAME=$(basename "$FILE")
  local OUT_DIR="$OUTPUT_DIR/SeprateFile/${BASENAME%.*}/RoughViralPrediction/virsorter2"
  mkdir -p "$OUT_DIR"

  if [ ! -f "$OUT_DIR/final-viral-score.tsv" ]; then
    echo "Running VirSorter2 prediction for $FILE..."
    virsorter run -w "$OUT_DIR" -i "$FILE" --include-groups "$Group" \
    -j "$THREADS_PER_FILE" all --min-score 0.5 --min-length 300 \
    --keep-original-seq -d "$DATABASE/db" > "$OUT_DIR/virsorter.log" 2>&1
  else
    echo "VirSorter2 prediction already completed for $FILE, skipping..."
  fi
}

# Function to run Genomad
run_genomad() {
  local FILE=$1
  local BASENAME=$(basename "$FILE")
  local OUT_DIR="$OUTPUT_DIR/SeprateFile/${BASENAME%.*}/RoughViralPrediction/genomadres"
  mkdir -p "$OUT_DIR"

  if [ ! -f "$OUT_DIR/${BASENAME}_summary/${BASENAME}_virus_summary.tsv" ]; then
    echo "Running Genomad prediction for $FILE..."
    if [ "$CONCENTRATION_TYPE" == "concentration" ]; then
      genomad end-to-end --enable-score-calibration "$FILE" "$OUT_DIR" "$DATABASE/genomad_db" \
      -t "$THREADS_PER_FILE" --min-score 0.7 --max-fdr 0.05 --min-number-genes 0 \
      --min-virus-marker-enrichment 1.5 --min-plasmid-marker-enrichment 0 --min-plasmid-hallmarks 1 \
      --min-plasmid-hallmarks-short-seqs 0 --max-uscg 2 
    else
      genomad end-to-end --enable-score-calibration "$FILE" "$OUT_DIR" "$DATABASE/genomad_db" \
      -t "$THREADS_PER_FILE" --min-score 0.8 --max-fdr 0.05 --min-number-genes 1 \
      --min-virus-marker-enrichment 0 --min-plasmid-marker-enrichment 1.5 --min-plasmid-hallmarks 1 \
      --min-plasmid-hallmarks-short-seqs 1 --max-uscg 2 
    fi
  else
    echo "Genomad prediction already completed for $FILE, skipping..."
  fi
}

# Export the functions to be used by parallel
export -f run_viralverify
export -f run_virsorter
export -f run_genomad

# Main function to process files with all three predictions in parallel
process_file() {
  local FILE=$1
  run_viralverify "$FILE"
  run_virsorter "$FILE"
  run_genomad "$FILE"
}

export -f process_file

# Run the process_file function in parallel for all input files
niceload --load "$THREADS_PER_FILE" process_file ::: $FILES

# Check if all Virsorter2 tasks are completed (only when it is concentration)
all_tasks_completed=false
while [ "$all_tasks_completed" == "false" ]; do
  all_tasks_completed=true
  for FILE in $FILES; do
    BASENAME=$(basename "$FILE" .fa)
    BASENAME=${BASENAME%.fasta}
    Virsorter_dir="$OUTPUT_DIR/SeprateFile/${BASENAME}/RoughViralPrediction/virsorter2"

    if [ ! -f "$Virsorter_dir/final-viral-score.tsv" ]; then
      all_tasks_completed=false
      echo "Virsorter2 still in processing"
      break
    fi
  done

  if [ "$all_tasks_completed" == "false" ]; then
    sleep 30
  fi
done

echo "All files have been processed."