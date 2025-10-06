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
SORTED_VCF="${OUT_VCF%.vcf.gz}.sorted.vcf.gz"
FINAL_VCF="${OUT_VCF%.vcf.gz}.sorted.filtered.vcf.gz"

mkdir -p $(dirname "${OUT_VCF}")
# workspace for bcftools sort temp files
mkdir -p ./tmp_sort

# Fast path: if final outputs already exist, clean intermediates and exit
if [ -f "${FINAL_VCF}" ] && [ -f "${FINAL_VCF}.tbi" ]; then
  echo "[INFO] Found existing FINAL VCF: ${FINAL_VCF}; cleaning intermediates and exiting."
  # Remove intermediates if they exist
  rm -f "${OUT_VCF}" "${OUT_VCF}.tbi" "${SORTED_VCF}" "${SORTED_VCF}.tbi" || true
  exit 0
fi

# If sorted is already present, skip normalization/annotation
SKIP_NORM=0
if [ -f "${SORTED_VCF}" ]; then
  echo "[INFO] Detected existing sorted VCF: ${SORTED_VCF}; will skip normalization/annotation."
  SKIP_NORM=1
fi

if [ "${SKIP_NORM}" -eq 0 ]; then
  echo "[INFO] Start normalization -> set-id (threads=${THREADS})"
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

  echo "[INFO] Sorting ${OUT_VCF} -> ${SORTED_VCF}"
  time bcftools sort -T ./tmp_sort -Oz -o "${SORTED_VCF}" "${OUT_VCF}"
  echo "[INFO] Indexing ${SORTED_VCF} (tbi)"
  time bcftools index -t --threads "${THREADS}" "${SORTED_VCF}"
fi

# --- Filtering on sorted+indexed VCF: remove SV (symbolic ALT) and REF/ALT==N ---
if [ ! -f "${SORTED_VCF}" ]; then
  echo "[FATAL] ${SORTED_VCF} not found. Cannot proceed to filtering." >&2
  exit 1
fi

echo "[INFO] Filtering (keep snps,indels; drop symbolic ALT and REF/ALT=N): ${SORTED_VCF} -> ${FINAL_VCF}"
time \
  bcftools view -v snps,indels "${SORTED_VCF}" \
  | bcftools view -e 'ALT ~ "^<" || REF = "N" || ALT = "N"' \
  -Oz -o "${FINAL_VCF}"

echo "[INFO] Indexing FINAL ${FINAL_VCF} (tbi)"
 time bcftools index -t --threads "${THREADS}" "${FINAL_VCF}"

# --- Cleanup intermediates: keep only FINAL files ---
rm -f "${OUT_VCF}" "${OUT_VCF}.tbi" "${SORTED_VCF}" "${SORTED_VCF}.tbi" || true

echo "[DONE] Output: ${FINAL_VCF} and ${FINAL_VCF}.tbi"