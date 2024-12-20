#!/usr/bin/env bash

# Export necessary variables
export OUTPUT_DIR DATABASE CONCENTRATION_TYPE ScriptDir FILES THREADS

# Uncomment to enable debug mode
#set -x

process_file() {
  local FILE=$1
  echo "Processing file: $FILE"

  local EXTENSION="${FILE##*.}"
  BASENAME=$(basename "$FILE" .fa)
  BASENAME=$(basename "$BASENAME" .fasta)
  echo "Basename: $BASENAME"

  out_dir="$OUTPUT_DIR/SeprateFile/${BASENAME}"
  PREDICTION_DIR="$out_dir/RoughViralPrediction"
  genomad_dir="$PREDICTION_DIR/genomadres"
  viralverify_dir="$PREDICTION_DIR/viralverify"
  virsorter_dir="$PREDICTION_DIR/virsorter2"

  # Skip the file if the filtered fasta already exists
  if [ -f "$out_dir/${BASENAME}_filtered.fasta" ]; then
    echo "Skipping $BASENAME as quality_summary.tsv already exists."
    return 0
  fi

  # Perform cross-validation for virus contigs
  echo -e "\n\n\n# Performing cross-validation for virus contigs!!!\n\n\n"
  echo "Running CrossValidAnA.py..."
  if ! python "${ScriptDir}/CrossValid.py" "$genomad_dir" "$viralverify_dir" "$virsorter_dir" "$BASENAME" "$out_dir" "$CONCENTRATION_TYPE"; then
    echo "Error during cross-validation for $BASENAME. Exiting..."
    return 1
  fi

  # Extract sequences from raw results
  echo "Running FilterRawResSeqs.py..."
  if ! python "${ScriptDir}/FilterRawResSeqs.py" "$FILE" "$BASENAME" "$out_dir"; then
    echo "Error during sequence extraction for $BASENAME. Exiting..."
    return 1
  fi

  # Create directory for CheckV results
  echo "Creating CheckV results directory..."
  mkdir -p "$out_dir/${BASENAME}_CheckRes"

  # Run CheckV to evaluate the quality of viral sequences
  echo "Running CheckV..."
  if ! checkv end_to_end "$out_dir/${BASENAME}_filtered.fasta" "$out_dir/${BASENAME}_CheckRes" -t "${THREADS}" -d "$DATABASE/checkv-db-v1.5"; then
    echo "Error during CheckV for $BASENAME. Exiting..."
    return 1
  fi

  # Run check_removal.py to remove low-quality sequences based on the quality report
  echo "Running check_removal.py..."
  if ! python "${ScriptDir}/check_removal.py" "$out_dir/${BASENAME}_CheckRes/quality_summary.tsv" "$out_dir/${BASENAME}_filtered.fasta"; then
    echo "Error during check removal for $BASENAME. Exiting..."
    return 1
  fi

  echo "Processing for $BASENAME completed."

  # Remove the CheckV results directory to free up space
  echo "Removing CheckV results directory..."
  rm -rf "$out_dir/${BASENAME}_CheckRes"
}

export -f process_file

# Run in parallel
echo "Starting parallel processing..."
parallel process_file ::: $FILES

echo "CrossValid analysis completed."
