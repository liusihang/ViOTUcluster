#!/usr/bin/env bash

# Set error handling mechanism: if an error occurs, the script will stop executing
set -e
trap 'echo "An error occurred. Exiting..."; exit 1;' ERR

# Make sure BASE_CONDA_PREFIX is correctly set
BASE_CONDA_PREFIX=$(conda info --base)
conda_sh="$BASE_CONDA_PREFIX/etc/profile.d/conda.sh"

# Export necessary variables and functions
export OUTPUT_DIR DATABASE Group CONCENTRATION_TYPE THREADS_PER_FILE FILES

# Function
#python -c 'import os; print(os.environ.get("OUTPUT_DIR"))'
python "${ScriptDir}/viralprediction.py"

# Check if all Virsorter2 tasks are completed 
all_tasks_completed=false
while [ "$all_tasks_completed" == "false" ]; do
  all_tasks_completed=true
  for FILE in $FILES; do
    BASENAME=$(basename "$FILE" .fa)
    BASENAME=${BASENAME%.fasta}
    Virsorter_dir="$OUTPUT_DIR/SeprateFile/${BASENAME}/RoughViralPrediction/virsorter2"
    echo "$Virsorter_dir"

    if [ ! -f "$Virsorter_dir/final-viral-score.tsv" ]; then
      all_tasks_completed=false
      echo "$(date '+%Y-%m-%d %H:%M:%S') - Viral prediction still in processing"
      break
    fi
  done

  if [ "$all_tasks_completed" == "false" ]; then
    sleep 30
  fi
done


echo "$(date '+%Y-%m-%d %H:%M:%S') - All Viral predictions completed."