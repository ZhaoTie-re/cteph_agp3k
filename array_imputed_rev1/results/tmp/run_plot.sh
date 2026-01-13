#!/bin/bash
#SBATCH --job-name=plot_gt_ds
#SBATCH --output=plot_gt_ds.log
#SBATCH --error=plot_gt_ds.err
#SBATCH --time=1:00:00

source activate cteph_geno_pro

INPUT_LIST="/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/array_imputed_rev1/results/tmp/gt_ds_analysis/check.ls"
VCF_DIR="/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/array_imputed_rev1/results/01.normalized"
OUT_DIR="/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/array_imputed_rev1/results/tmp/gt_ds_analysis"
SCRIPT="/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/array_imputed_rev1/scripts/plot_gt_vs_ds.py"
TABIX="/home/b/b37974/htslib-1.9/tabix"

mkdir -p ${OUT_DIR}

python3 ${SCRIPT} \
    --variant-list ${INPUT_LIST} \
    --vcf-dir ${VCF_DIR} \
    --output-dir ${OUT_DIR} \
    --tabix ${TABIX}
