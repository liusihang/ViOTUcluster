#!/usr/bin/env bash

#!/usr/bin/env bash

# 确保脚本在出现错误时停止
set -e
trap 'echo "An error occurred in script execution. Exiting..."; exit 1;' ERR

# Perform CrossValid analysis
for FILE in $FILES; do
  echo "Processing $FILE"

  BASENAME=$(basename "$FILE" .fa)
  BASENAME=${BASENAME%.fasta}
  OUT_DIR="$OUTPUT_DIR/SeprateFile/${BASENAME}"
  PREDICTION_DIR="$OUT_DIR/RoughViralPrediction"
  Genomad_dir="$PREDICTION_DIR/genomadres"
  Viralverify_dir="$PREDICTION_DIR/viralverify"
  Virsorter_dir="$PREDICTION_DIR/virsorter2"

  # CrossValid virus contigs
  echo -e "\n \n \n # Performing cross-validation for virus contigs!!! \n \n \n"
  if ! python "${ScriptDir}/CrossValid.py" "$Genomad_dir" "$Viralverify_dir" "$Virsorter_dir" "$BASENAME" "$OUT_DIR" "$CONCENTRATION_TYPE"; then
    echo "Error during cross-validation. Exiting..."
    exit 1
  fi

  # Extract Sequences
  if ! python "${ScriptDir}/FilterRawResSeqs.py" "$FILE" "$BASENAME" "$OUT_DIR"; then
    echo "Error during sequence extraction. Exiting..."
    exit 1
  fi

  # CheckVFilter
  mkdir -p "$OUT_DIR/${BASENAME}_CheckRes"
  if ! checkv end_to_end "$OUT_DIR/${BASENAME}_filtered.fasta" "$OUT_DIR/${BASENAME}_CheckRes" -t 104 -d "$DATABASE/checkv-db-v1.5"; then
    echo "Error during CheckVFilter. Exiting..."
    exit 1
  fi
  
  if ! python "${ScriptDir}/check_removal.py" "$OUT_DIR/${BASENAME}_CheckRes/quality_summary.tsv" "$OUT_DIR/${BASENAME}_filtered.fasta"; then
    echo "Error during check removal. Exiting..."
    exit 1
  fi
done

echo "CrossValid analysis completed."
