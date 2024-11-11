#!/usr/bin/env bash

# Set error handling mechanism: if an error occurs, the script will stop executing
set -e
trap 'echo "An error occurred. Exiting..."; exit 1;' ERR

# Make sure BASE_CONDA_PREFIX is correctly set
BASE_CONDA_PREFIX=$(conda info --base)
conda_sh="$BASE_CONDA_PREFIX/etc/profile.d/conda.sh"

# Export necessary variables and functions
export OUTPUT_DIR DATABASE Group CONCENTRATION_TYPE THREADS_PER_FILE FILES

# Main Function
python -c 'import os; print(os.environ.get("OUTPUT_DIR"))'
python "${ScriptDir}/viralprediction.py"

# Check if all tasks are completed
all_tasks_completed=false
previous_incomplete_count=-1

while [ "$all_tasks_completed" == "false" ]; do
  all_tasks_completed=true
  incomplete_count=0

  for FILE in $FILES; do
    BASENAME=$(basename "$FILE" .fa)
    BASENAME=${BASENAME%.fasta}
    Virsorter_dir="$OUTPUT_DIR/SeprateFile/${BASENAME}/RoughViralPrediction/virsorter2"
    Genomad_dir="$OUTPUT_DIR/SeprateFile/${BASENAME}/RoughViralPrediction/genomadres"
    ViralVerify_dir="$OUTPUT_DIR/SeprateFile/${BASENAME}/RoughViralPrediction/viralverify"

    # Define file paths
    virsorter_file="$Virsorter_dir/final-viral-score.tsv"
    genomad_file="$Genomad_dir/${BASENAME}_summary/${BASENAME}_virus.fna"
    viralverify_file="$ViralVerify_dir/${BASENAME}_result_table.csv"

    # Check completion status based on CONCENTRATION_TYPE
    if [ "$CONCENTRATION_TYPE" == "concentration" ]; then
      if [ ! -f "$virsorter_file" ] || [ ! -f "$genomad_file" ] || [ ! -f "$viralverify_file" ]; then
        all_tasks_completed=false
        ((incomplete_count++))
      fi
    else
      # For non-concentration, only check genomad and viralverify
      if [ ! -f "$genomad_file" ] || [ ! -f "$viralverify_file" ]; then
        all_tasks_completed=false
        ((incomplete_count++))
      fi
    fi
  done

  # Log only if the incomplete count has changed
  if [ "$incomplete_count" -ne "$previous_incomplete_count" ]; then
    echo "$(date '+%Y-%m-%d %H:%M:%S') - Still $incomplete_count files remaining"
    previous_incomplete_count=$incomplete_count
  fi

  # Sleep for 30 seconds if tasks are still in progress
  if [ "$all_tasks_completed" == "false" ]; then
    sleep 30
  fi
done

echo "$(date '+%Y-%m-%d %H:%M:%S') - All Viral predictions completed."