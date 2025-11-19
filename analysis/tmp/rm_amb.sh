#!/bin/bash

# Script to run plink2 GWAS analysis while removing ambiguous samples
# Remove samples: PHOM0518, PHOM0582

# Set plink2 path
PLINK2="/home/b/b37974/plink2"

# Define input/output paths
BFILE="/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/wgs/19.tommo_panel_filter/cteph_agp3k.lowfreq_common"
COVAR="/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/wgs/20.cov_pheno_prepare/cteph_agp3k.bbj.projection.cov_df.no_age.csv"
PHENO="/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/wgs/20.cov_pheno_prepare/cteph_agp3k.bbj.projection.pheno_df.csv"
OUTPUT_PREFIX="/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/analysis/tmp/cteph_agp3k_rm_amb"

# Create a file with samples to remove
REMOVE_FILE="/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/analysis/tmp/remove_samples.txt"

echo "Creating remove sample file: $REMOVE_FILE"
# Create file with sample IDs to remove (FID IID format, FID=IID)
cat > $REMOVE_FILE << EOF
PHOM0518	PHOM0518
PHOM0582	PHOM0582
EOF

echo "Starting plink2 GWAS analysis..."
echo "Removing samples listed in: $REMOVE_FILE"
echo "Output prefix: $OUTPUT_PREFIX"

# Run plink2 with sample removal
$PLINK2 \
  --bfile $BFILE \
  --remove $REMOVE_FILE \
  --ci 0.95 \
  --covar $COVAR \
  --covar-name SEX,PC1_AVG-PC10_AVG \
  --glm omit-ref no-firth hide-covar \
  --out $OUTPUT_PREFIX \
  --pheno $PHENO \
  --pheno-name PHENO1 \
  --threads 16

echo "Analysis completed. Check output files with prefix: $OUTPUT_PREFIX"