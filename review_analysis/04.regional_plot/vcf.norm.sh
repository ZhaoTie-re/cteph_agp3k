#!/bin/bash
#SBATCH -p gr10478b
#SBATCH -t 168:0:0
#SBATCH --rsc p=1:t=64:c=32:m=146272M
#SBATCH --job-name=vcf_norm

set -euo pipefail

# ------------------------------
# Environment
# ------------------------------
export PATH=/home/b/b37974/:$PATH

vcf='/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/review_analysis/04.regional_plot/EAS.ALL.split_norm_af.1kg_30x.hg38.vcf.gz'
nagasaki_pipeline='/LARGE0/gr10478/b37974/Pulmonary_Hypertension/nagasaki_pipeline/data'

# Threads: prefer SLURM_CPUS_PER_TASK if available
THREADS=${SLURM_CPUS_PER_TASK:-64}

# Output
OUT_VCF='EAS.ALL.split_norm_af.1kg_30x.hg38.norm.setid.vcf.gz'

mkdir -p $(dirname "${OUT_VCF}")

echo "[INFO] Start normalization -> set-id -> index (threads=${THREADS})"

# ------------------------------
# Accelerated pipeline:
# 1) bcftools norm outputs uncompressed BCF (-Ou) to avoid extra disk I/O/compression
# 2) pipe to bcftools annotate to set ID as CHROM:POS:REF:ALT, compress once (-Oz)
# ------------------------------

time \
  bcftools norm \
    --multiallelics -any \
    --fasta-ref "${nagasaki_pipeline}/hs38DH.fa" \
    --check-ref s \
    --threads "${THREADS}" \
    -Ou "${vcf}" \
  | bcftools annotate \
      --set-id '%CHROM:%POS:%REF:%ALT' \
      --threads "${THREADS}" \
      -Oz -o "${OUT_VCF}"

# Index the final VCF (tbi)
echo "[INFO] Indexing ${OUT_VCF} (tbi)"
 time bcftools index -t --threads "${THREADS}" "${OUT_VCF}"

echo "[DONE] Output: ${OUT_VCF} and ${OUT_VCF}.tbi"