#!/bin/bash
# 01_sample_select.sh
# Usage: 01_sample_select.sh <plink2_path> <sample_list> <prefix>

PLINK2_PATH=$1
SAMPLE_LIST=$2
PREFIX=$3

# 格式化样本列表：添加列名并复制为两列（FID和IID）
echo -e "#FID\tIID" > sample_list.formatted.txt
awk '{print $1"\t"$1}' ${SAMPLE_LIST} >> sample_list.formatted.txt

# 选择样本
${PLINK2_PATH} \
    --bfile ${PREFIX} \
    --keep sample_list.formatted.txt \
    --make-bed \
    --out ${PREFIX}.selected \
    --threads 4
