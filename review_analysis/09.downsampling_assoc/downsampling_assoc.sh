#!/bin/bash

# Script for downsampling association analysis
# This script extracts common samples and variants between WGS and array data,
# then performs association analysis

set -e

# Define paths
PLINK2="/home/b/b37974/plink2"
WGS_PREFIX="/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/wgs/19.tommo_panel_filter/cteph_agp3k.lowfreq_common"
ARRAY_PREFIX="/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/snippet_analysis/02.array_check_impute/03.extract_common/cteph_agp3k.array.common"
COV_FILE="/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/wgs/20.cov_pheno_prepare/cteph_agp3k.bbj.projection.cov_df.no_age.csv"
PHENO_FILE="/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/wgs/20.cov_pheno_prepare/cteph_agp3k.bbj.projection.pheno_df.csv"

# Output directory
OUT_DIR="/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/review_analysis/09.downsampling_assoc"
OUT_PREFIX="cteph_agp3k.downsampled"

echo "=================================================="
echo "Starting downsampling association analysis"
echo "=================================================="
echo ""

# Step 1: Extract variant IDs from array data (based on column 2 of .bim file)
echo "Step 1: Extracting variant IDs from array data..."
awk '{print $2}' ${ARRAY_PREFIX}.bim > ${OUT_DIR}/array_variants.txt
ARRAY_VARIANT_COUNT=$(wc -l < ${OUT_DIR}/array_variants.txt)
echo "Array variants: ${ARRAY_VARIANT_COUNT}"
echo ""

# Step 2: Extract sample IDs from array data
echo "Step 2: Extracting sample IDs from array data..."
awk '{print $1, $2}' ${ARRAY_PREFIX}.fam > ${OUT_DIR}/array_samples.txt
ARRAY_SAMPLE_COUNT=$(wc -l < ${OUT_DIR}/array_samples.txt)
echo "Array samples: ${ARRAY_SAMPLE_COUNT}"
echo ""

# Step 3: Extract common variants and samples from WGS data
echo "Step 3: Extracting common variants and samples from WGS data..."
${PLINK2} \
  --bfile ${WGS_PREFIX} \
  --extract ${OUT_DIR}/array_variants.txt \
  --keep ${OUT_DIR}/array_samples.txt \
  --make-bed \
  --out ${OUT_DIR}/${OUT_PREFIX} \
  --threads 16

echo ""

# Get final counts
FINAL_VARIANT_COUNT=$(wc -l < ${OUT_DIR}/${OUT_PREFIX}.bim)
FINAL_SAMPLE_COUNT=$(wc -l < ${OUT_DIR}/${OUT_PREFIX}.fam)

echo "=================================================="
echo "Subset Summary:"
echo "  Variants extracted: ${FINAL_VARIANT_COUNT}"
echo "  Samples extracted: ${FINAL_SAMPLE_COUNT}"
echo "=================================================="
echo ""

# Step 4: Perform association analysis
echo "Step 4: Performing association analysis..."
${PLINK2} \
  --bfile ${OUT_DIR}/${OUT_PREFIX} \
  --ci 0.95 \
  --covar ${COV_FILE} \
  --covar-name SEX,PC1_AVG-PC10_AVG \
  --glm omit-ref no-firth hide-covar \
  --out ${OUT_DIR}/${OUT_PREFIX}.sex.10pc.additive \
  --pheno ${PHENO_FILE} \
  --pheno-name PHENO1 \
  --threads 16

echo ""
echo "=================================================="
echo "Association analysis completed!"
echo "Output files:"
echo "  - ${OUT_DIR}/${OUT_PREFIX}.bed/bim/fam (extracted genotypes)"
echo "  - ${OUT_DIR}/${OUT_PREFIX}.sex.10pc.additive.PHENO1.glm.linear (results)"
echo "=================================================="
