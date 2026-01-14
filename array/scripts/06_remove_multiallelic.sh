#!/bin/bash
# 06_remove_multiallelic.sh
# Usage: 06_remove_multiallelic.sh <plink2_path> <prefix>

PLINK2_PATH=$1
PREFIX=$2

echo "=== Multiallelic Variants Check Report ===" > multiallelic_check.log
echo "Date: $(date)" >> multiallelic_check.log
echo "" >> multiallelic_check.log

# 统计输入文件信息
input_variant_count=$(wc -l < ${PREFIX}.bim)
input_sample_count=$(wc -l < ${PREFIX}.fam)

printf "Input file: %s\n" "${PREFIX}" >> multiallelic_check.log
printf "Input variants: %'d\n" ${input_variant_count} >> multiallelic_check.log
printf "Input samples: %'d\n" ${input_sample_count} >> multiallelic_check.log
echo "" >> multiallelic_check.log

# 识别multiallelic variants（相同chr:pos的变体）
echo "Identifying multiallelic variants (same chr:pos)..." >> multiallelic_check.log

# 提取chr和pos（bim文件的第1列和第4列），统计重复
awk '{print $1":"$4}' ${PREFIX}.bim | sort | uniq -c | awk '$1 > 1 {print $2}' > multiallelic_sites.txt

multiallelic_site_count=$(wc -l < multiallelic_sites.txt)
printf "Multiallelic sites found: %'d\n" ${multiallelic_site_count} >> multiallelic_check.log
echo "" >> multiallelic_check.log

if [ ${multiallelic_site_count} -gt 0 ]; then
    # 提取所有multiallelic位点的详细信息
    echo "Extracting details of multiallelic variants..." >> multiallelic_check.log
    
    # 创建multiallelic variant ID列表
    touch multiallelic_variants.txt
    while IFS= read -r site; do
        chr=$(echo "$site" | cut -d: -f1)
        pos=$(echo "$site" | cut -d: -f2)
        # 找到所有该位点的变体ID
        awk -v chr="$chr" -v pos="$pos" '$1==chr && $4==pos {print $2}' ${PREFIX}.bim >> multiallelic_variants.txt
    done < multiallelic_sites.txt
    
    multiallelic_variant_count=$(wc -l < multiallelic_variants.txt)
    printf "Total variants at multiallelic sites: %'d\n" ${multiallelic_variant_count} >> multiallelic_check.log
    echo "" >> multiallelic_check.log
    
    # 显示前10个multiallelic位点的示例
    echo "Example multiallelic sites (first 10):" >> multiallelic_check.log
    head -10 multiallelic_sites.txt | while IFS= read -r site; do
        chr=$(echo "$site" | cut -d: -f1)
        pos=$(echo "$site" | cut -d: -f2)
        echo "  Site: $site" >> multiallelic_check.log
        awk -v chr="$chr" -v pos="$pos" '$1==chr && $4==pos {printf "    - ID: %s, ALT: %s, REF: %s\n", $2, $5, $6}' ${PREFIX}.bim >> multiallelic_check.log
    done
    echo "" >> multiallelic_check.log
    
    # 使用PLINK去除multiallelic variants
    echo "Removing multiallelic variants with PLINK2..." >> multiallelic_check.log
    
    ${PLINK2_PATH} \
        --bfile ${PREFIX} \
        --exclude multiallelic_variants.txt \
        --make-bed \
        --out ${PREFIX}.biallelic \
        --threads 8
    
    # 统计过滤后的结果
    output_variant_count=$(wc -l < ${PREFIX}.biallelic.bim)
    removed_count=$((input_variant_count - output_variant_count))
    
    printf "Variants after removal: %'d\n" ${output_variant_count} >> multiallelic_check.log
    printf "Variants removed: %'d\n" ${removed_count} >> multiallelic_check.log
    removal_rate=$(awk "BEGIN {printf \"%.2f%%\", (${removed_count}/${input_variant_count})*100}")
    echo "Removal rate: ${removal_rate}" >> multiallelic_check.log
    
else
    echo "No multiallelic variants detected!" >> multiallelic_check.log
    echo "All variants are biallelic. Copying files without modification..." >> multiallelic_check.log
    
    # 如果没有multiallelic variants，直接复制文件
    cp ${PREFIX}.bed ${PREFIX}.biallelic.bed
    cp ${PREFIX}.bim ${PREFIX}.biallelic.bim
    cp ${PREFIX}.fam ${PREFIX}.biallelic.fam
    
    printf "Output variants: %'d (unchanged)\n" ${input_variant_count} >> multiallelic_check.log
fi

echo "" >> multiallelic_check.log
echo "======================================" >> multiallelic_check.log
echo "SUMMARY" >> multiallelic_check.log
echo "======================================" >> multiallelic_check.log

final_variant_count=$(wc -l < ${PREFIX}.biallelic.bim)
final_sample_count=$(wc -l < ${PREFIX}.biallelic.fam)

printf "Input variants: %'d\n" ${input_variant_count} >> multiallelic_check.log
printf "Output variants: %'d\n" ${final_variant_count} >> multiallelic_check.log
printf "Samples (unchanged): %'d\n" ${final_sample_count} >> multiallelic_check.log
echo "" >> multiallelic_check.log

echo "Output files:" >> multiallelic_check.log
echo "  - ${PREFIX}.biallelic.bed (binary genotype file, biallelic only)" >> multiallelic_check.log
echo "  - ${PREFIX}.biallelic.bim (variant information, biallelic only)" >> multiallelic_check.log
echo "  - ${PREFIX}.biallelic.fam (sample information, unchanged)" >> multiallelic_check.log
echo "  - multiallelic_check.log (this log file)" >> multiallelic_check.log
if [ ${multiallelic_site_count} -gt 0 ]; then
    echo "  - multiallelic_variants.txt (list of removed variant IDs)" >> multiallelic_check.log
fi

# 清理临时文件
rm -f multiallelic_sites.txt
