#!/bin/bash
# 02_sample_qc.sh
# Usage: 02_sample_qc.sh <plink2_path> <prefix>

PLINK2_PATH=$1
PREFIX=$2

# 计算样本缺失率
${PLINK2_PATH} \
    --bfile ${PREFIX} \
    --missing sample-only \
    --out ${PREFIX} \
    --threads 4

# 过滤call-rate < 99%的样本 (即缺失率 > 0.01)
${PLINK2_PATH} \
    --bfile ${PREFIX} \
    --mind 0.01 \
    --make-bed \
    --out ${PREFIX}.sqc \
    --threads 4
