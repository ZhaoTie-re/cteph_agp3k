#!/bin/bash
#SBATCH -p gr10478b
#SBATCH -t 48:0:0
#SBATCH --rsc p=1:t=16:c=8:m=36568M
#SBATCH --job-name=association

# Initialize environment
export PATH=/home/b/b37974/:$PATH

BED_FILE_PR="/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/wgs/18.clean_gt/cteph_agp3k.clean"
COVAR_FILE="/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/wgs/18.clean_gt/covariance.txt"
PHENO_FILE="/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/wgs/18.clean_gt/pheno.txt"

# Run association analysis
# sex + 10 PCs

plink2 \
    --bfile ${BED_FILE_PR} \
    --pheno ${PHENO_FILE} \
    --pheno-name PHENO1 \
    --covar ${COVAR_FILE} \
    --covar-name SEX, PC1_AVG-PC10_AVG \
    --glm no-firth hide-covar \
    --out cteph_agp3k.sex_pcs \
    --ci 0.95 \
    --threads 16

