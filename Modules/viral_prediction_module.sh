#!/usr/bin/env bash

# Set error handling mechanism: if an error occurs, the script will stop executing
set -euo pipefail
trap 'echo "[❌] An error occurred. Exiting..."; exit 1;' ERR

export OUTPUT_DIR DATABASE Group CONCENTRATION_TYPE THREADS_PER_FILE FILES MAX_PredictionTASKS
VIOTUCLUSTER_PYTHON=${VIOTUCLUSTER_PYTHON:-python}

# Run main viral prediction script
echo "[🔄] Starting viral prediction pipeline..."
"$VIOTUCLUSTER_PYTHON" -m ViOTUcluster.viralprediction

echo "[✅] $(date '+%Y-%m-%d %H:%M:%S') - All viral predictions completed successfully!"
