#!/bin/bash

# Configuration
VCF_PATH="/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/analysis/assoc_rvtest.rev1/results/02.snpeff_annotate/cteph_agp3k.rare.mac2.rm_samples.all.nochr.chrprefix.snpeff.vcf.gz"
OUT_DIR="/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/analysis/assoc_rvtest.rev1/results/tmp/check_impact_effect"

# Ensure output directory exists
mkdir -p "${OUT_DIR}"

# Function to check for bcftools
check_dependency() {
    if ! command -v bcftools &> /dev/null; then
        echo "Error: bcftools is not installed or not in PATH."
        echo "Trying to interpret 'source activate'..."
        # Try to activate the environment seen in the nextflow file if possible, or assume user runs this in correct env
        if [ -f ~/.bashrc ]; then source ~/.bashrc; fi
        if command -v conda &> /dev/null; then
             source activate cteph_geno_pro
        fi
    fi
}

check_dependency

echo "Processing VCF: ${VCF_PATH}"
echo "Saving results to: ${OUT_DIR}"

# 1. Extract 5 examples: Effect=sequence_feature AND Impact=LOW
echo "Extracting LOW impact sequence_details..."
OUTPUT_LOW="${OUT_DIR}/examples_sequence_feature_LOW.vcf"
# Get header first
bcftools view -h "${VCF_PATH}" > "${OUTPUT_LOW}"
# Append 5 records
bcftools view -H -i 'INFO/effect="sequence_feature" && INFO/impact="LOW"' "${VCF_PATH}" | head -n 5 >> "${OUTPUT_LOW}"

# 2. Extract 5 examples: Effect=sequence_feature AND Impact=MODERATE
echo "Extracting MODERATE impact sequence_details..."
OUTPUT_MOD="${OUT_DIR}/examples_sequence_feature_MODERATE.vcf"
# Get header first
bcftools view -h "${VCF_PATH}" > "${OUTPUT_MOD}"
# Append 5 records
bcftools view -H -i 'INFO/effect="sequence_feature" && INFO/impact="MODERATE"' "${VCF_PATH}" | head -n 5 >> "${OUTPUT_MOD}"

echo "Done."
echo "examples_sequence_feature_LOW.vcf created."
echo "examples_sequence_feature_MODERATE.vcf created."
