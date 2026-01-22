#!/bin/bash
# Run the python script for outlier detection
# Description: Plots MeanDP vs SMinAC and identifies SMinAC outliers using robust Z-score

SCRIPT_PATH="/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/tuning.rv/scripts/plot_zscore_outliers.py"
PYTHON_EXE="/home/b/b37974/anaconda3/envs/cteph_geno_pro/bin/python"

# Parameters
INPUT_FILE="/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/tuning.rv/results/00.qc_metrics/minac2/sample_metrics.txt.gz"
OUTPUT_DIR="/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/tuning.rv/results/tmp"
THRESHOLD=3.0
DIRECTION="both" # Options: both, lower, upper

echo "Running script: $SCRIPT_PATH using python: $PYTHON_EXE"
echo "Input: $INPUT_FILE"
echo "Output: $OUTPUT_DIR"
echo "Threshold: $THRESHOLD"
echo "Direction: $DIRECTION"

$PYTHON_EXE "$SCRIPT_PATH" \
    --input "$INPUT_FILE" \
    --output-dir "$OUTPUT_DIR" \
    --threshold "$THRESHOLD" \
    --direction "$DIRECTION"
