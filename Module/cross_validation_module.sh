#!/usr/bin/env bash

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
  python "${ScriptDir}/CrossValid.py" "$Genomad_dir" "$Viralverify_dir" "$Virsorter_dir" "$BASENAME" "$OUT_DIR" "$CONCENTRATION_TYPE"

  # Extract Sequences
  python "${ScriptDir}/FilterRawResSeqs.py" "$FILE" "$BASENAME" "$OUT_DIR"

  # CheckVFilter
  mkdir -p "$OUT_DIR/${BASENAME}_CheckRes"
  checkv end_to_end "$OUT_DIR/${BASENAME}_filtered.fasta" "$OUT_DIR/${BASENAME}_CheckRes" -t 104 -d "$DATABASE/checkv-db-v1.5"
  python "${ScriptDir}/check_removal.py" "$OUT_DIR/${BASENAME}_CheckRes/quality_summary.tsv" "$OUT_DIR/${BASENAME}_filtered.fasta"
done