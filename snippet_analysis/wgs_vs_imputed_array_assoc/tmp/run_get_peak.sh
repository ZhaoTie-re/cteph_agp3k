#!/bin/bash

#########################################################################
# Script: run_get_peak.sh
# Purpose: Extract lead variants and peak regions from GWAS results
# Author: ZHAO TIE
# Date: December 2025
#########################################################################

#SBATCH --job-name=get_peak
#SBATCH --output=get_peak_%j.out
#SBATCH --error=get_peak_%j.err
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
source activate cteph_geno_pro

# Define paths
BASE_DIR="/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/snippet_analysis/wgs_vs_imputed_array_assoc"
SCRIPT_DIR="${BASE_DIR}/scrpts"
RESULT_DIR="${BASE_DIR}/results/02.assoc_result"
OUTPUT_DIR="${BASE_DIR}/tmp/peaks"

# Input file
ASSOC_FILE="${RESULT_DIR}/array/wgs_vs_array.array.sex.10pc.additive.PHENO1.glm.logistic"

# Parameters
SIG_LEVEL=5e-8
LEAD_RANGE=500000
WINDOW=500000
PREFIX="array_peaks"

# Check if input file exists
echo "Checking input file..."
if [ ! -f "$ASSOC_FILE" ]; then
    echo "ERROR: File not found: $ASSOC_FILE"
    exit 1
fi
echo "  ✓ Found: $ASSOC_FILE"

# Create output directory if it doesn't exist
mkdir -p "${OUTPUT_DIR}"

# Run the peak extraction script
echo ""
echo "Running lead variant and peak region extraction..."
echo "==============================================================================="

python3 "${SCRIPT_DIR}/get_peak_cli.py" \
    --assoc "${ASSOC_FILE}" \
    --sig-level ${SIG_LEVEL} \
    --lead-range ${LEAD_RANGE} \
    --window ${WINDOW} \
    --output-dir "${OUTPUT_DIR}" \
    --prefix "${PREFIX}"

echo ""
echo "==============================================================================="
echo "Analysis completed successfully!"
echo ""
echo "Output files:"
echo "  - ${OUTPUT_DIR}/${PREFIX}_lead_variants.tsv"
echo "  - ${OUTPUT_DIR}/${PREFIX}_peak_regions.tsv"
echo "  - ${OUTPUT_DIR}/${PREFIX}_peak*_variants.txt"
echo ""
echo "=========================================="
echo "Job completed at: $(date)"
echo "=========================================="
