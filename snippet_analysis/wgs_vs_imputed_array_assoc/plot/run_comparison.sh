#!/bin/bash

#########################################################################
# Script: run_comparison.sh
# Purpose: Compare WGS and Array association results with customizable
#          color metric (MAF, R2, etc.) and flexible visualization options
# Author: ZHAO TIE
# Date: December 2025
#########################################################################

#SBATCH --job-name=wgs_array_compare
#SBATCH --output=wgs_array_comparison_%j.out
#SBATCH --error=wgs_array_comparison_%j.err
#SBATCH --rsc p=1:t=48:c=24:m=109704M

set -e  # Exit on error
set -u  # Exit on undefined variable

# Print job information
echo "=========================================="
echo "Job ID: ${SLURM_JOB_ID}"
echo "Job Name: ${SLURM_JOB_NAME}"
echo "Node: ${SLURM_NODELIST}"
echo "Start Time: $(date)"
echo "=========================================="
echo ""

# Activate conda environment
echo "Activating conda environment: cteph_geno_pro"
# source ~/miniconda3/etc/profile.d/conda.sh
source activate cteph_geno_pro

# Define paths
BASE_DIR="/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/snippet_analysis/wgs_vs_imputed_array_assoc"
SCRIPT_DIR="${BASE_DIR}/scrpts"
RESULT_DIR="${BASE_DIR}/results/02.assoc_result"
ARRAY_QC_DIR="/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/snippet_analysis/wgs_vs_imputed_array/results/05.variant_stats_final"
OUTPUT_DIR="${BASE_DIR}/plot"

# Input files
WGS_ASSOC="${RESULT_DIR}/wgs/wgs_vs_array.wgs.sex.10pc.additive.PHENO1.glm.logistic"
ARRAY_ASSOC="${RESULT_DIR}/array/wgs_vs_array.array.sex.10pc.additive.PHENO1.glm.logistic"
VQC_FILE="${ARRAY_QC_DIR}/cteph_agp3k.imputed_array.vmiss_qc.variant_qc_final.tsv"

# Check if input files exist
echo "Checking input files..."
for file in "${WGS_ASSOC}" "${ARRAY_ASSOC}" "${VQC_FILE}"; do
    if [ ! -f "$file" ]; then
        echo "ERROR: File not found: $file"
        exit 1
    fi
    echo "  ✓ Found: $file"
done

# Create output directory if it doesn't exist
mkdir -p "${OUTPUT_DIR}"

# Run the comparison script
echo ""
echo "Running WGS vs Array comparison..."
echo "==============================================================================="

python3 "${SCRIPT_DIR}/wgs_vs_array_cli.py" \
    --wgs "${WGS_ASSOC}" \
    --array "${ARRAY_ASSOC}" \
    --vqc "${VQC_FILE}" \
    --color-col "MAF_ALL" \
    # --reverse-colormap \
    --output-dir "${OUTPUT_DIR}" \
    --sig-level 5e-8 \
    --top-n 20 \
    --prefix "wgs_array_comparison"

echo ""
echo "==============================================================================="
echo "Analysis completed successfully!"
echo ""
echo "Output files:"
echo "  - ${OUTPUT_DIR}/wgs_array_comparison_MAF_ALL.png"
echo "  - ${OUTPUT_DIR}/wgs_array_comparison_top20_discordant.tsv"
echo "  - ${OUTPUT_DIR}/wgs_array_comparison_*.log"
echo ""
echo "=========================================="
echo "Job completed at: $(date)"
echo "=========================================="
