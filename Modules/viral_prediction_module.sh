#!/usr/bin/env bash

# Set error handling mechanism: if an error occurs, the script will stop executing
set -e
trap 'echo "An error occurred. Exiting..."; exit 1;' ERR

# Make sure BASE_CONDA_PREFIX is correctly set
BASE_CONDA_PREFIX=$(conda info --base)
conda_sh="$BASE_CONDA_PREFIX/etc/profile.d/conda.sh"

# Export necessary variables and functions
export OUTPUT_DIR DATABASE Group CONCENTRATION_TYPE THREADS_PER_FILE FILES

# Main function to process files
process_file() {
  local FILE=$1
  local BASENAME=$(basename "$FILE")
  # Extract the extension
  local EXTENSION="${BASENAME##*.}"
  # Remove the extension if it's .fa or .fasta
  if [ "$EXTENSION" = "fa" ] || [ "$EXTENSION" = "fasta" ]; then
      BASENAME="${BASENAME%.*}"
  fi

  local OUT_DIR="$OUTPUT_DIR/SeprateFile/${BASENAME}"
  mkdir -p "$OUT_DIR"

  local PREDICTION_DIR="$OUT_DIR/RoughViralPrediction"
  mkdir -p "$PREDICTION_DIR"

  local Genomad_dir="$PREDICTION_DIR/genomadres"
  mkdir -p "$Genomad_dir"
  local Viralverify_dir="$PREDICTION_DIR/viralverify"
  mkdir -p "$Viralverify_dir"

  echo "Processing $FILE"

  # Run viralverify analysis with different parameters based on CONCENTRATION_TYPE
  if [ ! -f "$Viralverify_dir/${BASENAME}_result_table.csv" ]; then
    echo "Running viralverify prediction..."
    source ${conda_sh}
    if [ "$CONCENTRATION_TYPE" == "concentration" ]; then
      # Parameters for 'concentration' type
      viralverify -f "$FILE" -o "$Viralverify_dir" --hmm "$DATABASE/ViralVerify/nbc_hmms.hmm" \
      -t "$THREADS_PER_FILE" > "$Viralverify_dir/viralverify.log" 2>&1 &
    else
      # Parameters for 'non-concentration' type
      viralverify -f "$FILE" -o "$Viralverify_dir" --hmm "$DATABASE/ViralVerify/nbc_hmms.hmm" \
      -t "$THREADS_PER_FILE" > "$Viralverify_dir/viralverify.log" 2>&1 &
    fi
  else
    echo "viralverify prediction already completed for $FILE, skipping..."
  fi

  # Run VirSorter2 with different parameters based on CONCENTRATION_TYPE
  local Virsorter_dir="$PREDICTION_DIR/virsorter2"
  mkdir -p "$Virsorter_dir"

  if [ ! -f "$Virsorter_dir/final-viral-score.tsv" ]; then
    if [ "$CONCENTRATION_TYPE" == "concentration" ]; then
      echo "Running VirSorter2 prediction in concentration mode..."
      virsorter run -w "$Virsorter_dir" -i "$FILE" --include-groups "$Group" \
      -j "$THREADS_PER_FILE" all --min-score 0.5 --min-length 300 \
      --keep-original-seq -d "$DATABASE/db" > "$Virsorter_dir/virsorter.log" 2>&1 &
    else
      echo "Running VirSorter2 prediction in non-concentration mode..."
      virsorter run -w "$Virsorter_dir" -i "$FILE" --include-groups "$Group" \
      -j "$THREADS_PER_FILE" all --min-score 0.5 --min-length 300 \
      --keep-original-seq -d "$DATABASE/db" > "$Virsorter_dir/virsorter.log" 2>&1 &
    fi
  else
    echo "VirSorter2 prediction already completed for $FILE, skipping..."
  fi

  echo "All predictions submitted for $FILE"
}

export -f process_file
parallel process_file ::: $FILES

# Wait for background processes to finish
wait

# Perform Genomad analysis based on CONCENTRATION_TYPE
for FILE in $FILES; do
  local BASENAME=$(basename "$FILE")
  # Extract the extension
  local EXTENSION="${BASENAME##*.}"
  # Remove the extension if it's .fa or .fasta
  if [ "$EXTENSION" = "fa" ] || [ "$EXTENSION" = "fasta" ]; then
      BASENAME="${BASENAME%.*}"
  fi

  OUT_DIR="$OUTPUT_DIR/SeprateFile/${BASENAME}"
  PREDICTION_DIR="$OUT_DIR/RoughViralPrediction"

  # Genomad prediction
  echo -e "\n \n \n # Performing Genomad prediction!!! \n \n \n"
  Genomad_dir="$PREDICTION_DIR/genomadres"
  mkdir -p "$Genomad_dir"

  if [ ! -f "$Genomad_dir/${BASENAME}_summary/${BASENAME}_virus_summary.tsv" ]; then
    if [ "$CONCENTRATION_TYPE" == "concentration" ]; then
      echo "Running Genomad in concentration mode..."
      genomad end-to-end --enable-score-calibration "$FILE" "$Genomad_dir" "$DATABASE/genomad_db" \
      -t "$THREADS_PER_FILE" --default --min-score 0.7 --max-fdr 0.05 --min-number-genes 0 \
      --min-plasmid-marker-enrichment 0 --min-plasmid-hallmarks 1 \
      --min-plasmid-hallmarks-short-seqs 0 --max-uscg 2 
    else
      echo "Running Genomad in non-concentration mode..."
      genomad end-to-end --enable-score-calibration "$FILE" "$Genomad_dir" "$DATABASE/genomad_db" \
      -t "$THREADS_PER_FILE" --conservative --min-score 0.8 --max-fdr 0.05 --min-number-genes 1 \
      --min-plasmid-marker-enrichment 1.5 --min-plasmid-hallmarks 1 \
      --min-plasmid-hallmarks-short-seqs 1 --max-uscg 2 
    fi
    echo -e "\n \n \n # Genomad prediction completed!!! \n \n \n"
  else
    echo "Genomad prediction already completed for $FILE, skipping..."
  fi
done

echo "All files have been processed."