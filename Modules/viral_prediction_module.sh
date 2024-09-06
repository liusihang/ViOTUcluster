#!/usr/bin/env bash

# Set error handling mechanism: if an error occurs, the script will stop executing
set -e
trap 'echo "An error occurred. Exiting..."; exit 1;' ERR

# Make sure BASE_CONDA_PREFIX is correctly set
BASE_CONDA_PREFIX=$(conda info --base)
conda_sh="$BASE_CONDA_PREFIX/etc/profile.d/conda.sh"

# Export necessary variables and functions
export OUTPUT_DIR DATABASE Group CONCENTRATION_TYPE

# Main function to process files
process_file() {
  local FILE=$1
  local BASENAME=$(basename "$FILE" .fa)
  BASENAME=${BASENAME%.fasta}

  local OUT_DIR="$OUTPUT_DIR/SeprateFile/${BASENAME}"
  mkdir -p "$OUT_DIR"

  local PREDICTION_DIR="$OUT_DIR/RoughViralPrediction"
  mkdir -p "$PREDICTION_DIR"

  local Genomad_dir="$PREDICTION_DIR/genomadres"
  mkdir -p "$Genomad_dir"
  local Viralverify_dir="$PREDICTION_DIR/viralverify"
  mkdir -p "$Viralverify_dir"

  echo "Processing $FILE"

  # Run viralverify analysis
  if [ ! -f "$Viralverify_dir/${BASENAME}_result_table.csv" ]; then
    echo "Running viralverify prediction..."
    source ${conda_sh}
    viralverify -f "$FILE" -o "$Viralverify_dir" --hmm "$DATABASE/ViralVerify/nbc_hmms.hmm" -t 104 > "$Viralverify_dir/viralverify.log" 2>&1 &
    VVERIFY_PID=$!  # Get the process ID of viralverify
  else
    echo "viralverify prediction already completed for $FILE, skipping..."
  fi

  # Run Virsorter2 only when it is concentration
  if [ "$CONCENTRATION_TYPE" == "concentration" ]; then
    local Virsorter_dir="$PREDICTION_DIR/virsorter2"
    mkdir -p "$Virsorter_dir"

    if [ ! -f "$Virsorter_dir/final-viral-score.tsv" ]; then
      echo "Running Virsorter2 prediction..."
      virsorter run -w "$Virsorter_dir" -i "$FILE" --include-groups "$Group" -j 104 all --min-score 0.5 --min-length 2000 --keep-original-seq -d "$DATABASE/db" > "$Virsorter_dir/virsorter.log" 2>&1 &
      VSORTER_PID=$!  # Get the process ID of virsorter2
    else
      echo "Virsorter2 prediction already completed for $FILE, skipping..."
    fi
  fi

  # Wait for all background tasks to complete
  wait $VVERIFY_PID
  if [ "$CONCENTRATION_TYPE" == "concentration" ]; then
    wait $VSORTER_PID
  fi

  echo "All predictions completed for $FILE"
}

export -f process_file
parallel process_file ::: $FILES

# Perform Genomad analysis
for FILE in $FILES; do
  echo "Processing $FILE"
  BASENAME=$(basename "$FILE" .fa)
  BASENAME=${BASENAME%.fasta}

  OUT_DIR="$OUTPUT_DIR/SeprateFile/${BASENAME}"
  PREDICTION_DIR="$OUT_DIR/RoughViralPrediction"

  # Genomad prediction
  echo -e "\n \n \n # Performing genomad prediction!!! \n \n \n"
  Genomad_dir="$PREDICTION_DIR/genomadres"
  mkdir -p "$Genomad_dir"

  if [ ! -f "$Genomad_dir/${BASENAME}_summary/${BASENAME}_virus_summary.tsv" ]; then
    genomad end-to-end --enable-score-calibration "$FILE" "$Genomad_dir" "$DATABASE/genomad_db"
    echo -e "\n \n \n # Genomad prediction completed!!! \n \n \n"
  else
    echo "genomad prediction already completed for $FILE, skipping..."
  fi
done

# Check if all Virsorter2 tasks are completed (only when it is concentration)
if [ "$CONCENTRATION_TYPE" == "concentration" ]; then
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
fi

echo "All files have been processed."
