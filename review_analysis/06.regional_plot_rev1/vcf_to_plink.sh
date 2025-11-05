#!/bin/bash
#SBATCH --job-name=vcf_to_plink
#SBATCH --output=vcf_to_plink_%j.log
#SBATCH --error=vcf_to_plink_%j.err
#SBATCH -p gr10478b
#SBATCH -t 168:0:0
#SBATCH --rsc p=1:t=64:c=32:m=146272M

# Script to convert VCF.gz to plink format (bed/bim/fam)

# 输入输出文件（直接指定）
INPUT_VCF="/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/review_analysis/04.regional_plot/EAS.ALL.split_norm_af.1kg_30x.hg38.norm.setid.sorted.filtered.vcf.gz"
OUTPUT_PREFIX="eas_all"

# plink2路径
PLINK2="/home/b/b37974/plink2"

# 检查输入文件是否存在
if [ ! -f "$INPUT_VCF" ]; then
    echo "Error: Input VCF file not found: $INPUT_VCF"
    exit 1
fi

# 检查plink2是否存在
if [ ! -f "$PLINK2" ]; then
    echo "Error: plink2 not found: $PLINK2"
    exit 1
fi

echo "Converting VCF to plink format..."
echo "Input VCF: $INPUT_VCF"
echo "Output prefix: $OUTPUT_PREFIX"
echo "Using plink2: $PLINK2"
echo ""

# 执行转换
$PLINK2 \
    --vcf "$INPUT_VCF" \
    --make-bed \
    --out "$OUTPUT_PREFIX" \
    --threads 16

# 检查是否成功
if [ $? -eq 0 ]; then
    echo ""
    echo "Conversion completed successfully!"
    echo "Output files:"
    echo "  - ${OUTPUT_PREFIX}.bed"
    echo "  - ${OUTPUT_PREFIX}.bim"
    echo "  - ${OUTPUT_PREFIX}.fam"
    
    # 显示文件大小
    if [ -f "${OUTPUT_PREFIX}.bed" ]; then
        ls -lh "${OUTPUT_PREFIX}.bed" "${OUTPUT_PREFIX}.bim" "${OUTPUT_PREFIX}.fam"
    fi
else
    echo ""
    echo "Error: Conversion failed!"
    exit 1
fi
