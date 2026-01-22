#!/bin/bash
# Wrapper script for Sensitivity Data Preparation
# Args: <input_vcf> <summary_file> <out_prefix> <threads> <script_dir>

INPUT_VCF=$1
SUMMARY_FILE=$2
OUT_PREFIX=$3
THREADS=$4
SCRIPT_DIR=$5

source activate cteph_geno_pro

echo "Running Sensitivity Extraction..."
python ${SCRIPT_DIR}/sensitivity_main.py \
    --vcf ${INPUT_VCF} \
    --summary ${SUMMARY_FILE} \
    --out-prefix ${OUT_PREFIX} \
    --threads ${THREADS}

echo "Done."
