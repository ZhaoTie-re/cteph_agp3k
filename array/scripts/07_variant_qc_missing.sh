#!/bin/bash
# 07_variant_qc_missing.sh
# Usage: 07_variant_qc_missing.sh <plink2_path> <prefix>

PLINK2_PATH=$1
PREFIX=$2

echo "=== Variant Missing Rate QC Report ===" > variant_missing_qc.log
echo "Date: $(date)" >> variant_missing_qc.log
echo "" >> variant_missing_qc.log

# 统计输入文件信息
input_variant_count=$(wc -l < ${PREFIX}.bim)
input_sample_count=$(wc -l < ${PREFIX}.fam)

printf "Input file: %s\n" "${PREFIX}" >> variant_missing_qc.log
printf "Input variants: %'d\n" ${input_variant_count} >> variant_missing_qc.log
printf "Input samples: %'d\n" ${input_sample_count} >> variant_missing_qc.log
echo "" >> variant_missing_qc.log

# 统计样本性别分布（用于说明性染色体过滤策略）
echo "Sample sex distribution:" >> variant_missing_qc.log
male_count=$(awk '$5==1' ${PREFIX}.fam | wc -l)
female_count=$(awk '$5==2' ${PREFIX}.fam | wc -l)
unknown_count=$(awk '$5==0' ${PREFIX}.fam | wc -l)
printf "  Males: %'d\n" ${male_count} >> variant_missing_qc.log
printf "  Females: %'d\n" ${female_count} >> variant_missing_qc.log
printf "  Unknown: %'d\n" ${unknown_count} >> variant_missing_qc.log
echo "" >> variant_missing_qc.log

# 统计每条染色体的变异数
echo "Variants per chromosome (before filtering):" >> variant_missing_qc.log
cut -f1 ${PREFIX}.bim | sort | uniq -c | awk '{printf "  Chr %s: %'"'"'d variants\n", $2, $1}' >> variant_missing_qc.log
echo "" >> variant_missing_qc.log

# 计算变异缺失率
echo "Step 1: Calculating variant missing rates..." >> variant_missing_qc.log
${PLINK2_PATH} \
    --bfile ${PREFIX} \
    --missing variant-only \
    --out ${PREFIX} \
    --threads 8

# 分析vmiss文件
total_variants=$(tail -n +2 ${PREFIX}.vmiss | wc -l)
echo "  Total variants analyzed: ${total_variants}" >> variant_missing_qc.log
echo "" >> variant_missing_qc.log

# 统计不同缺失率范围的变异数
echo "Missing rate distribution:" >> variant_missing_qc.log
awk 'NR>1 {
    missing_rate = $5;
    if (missing_rate == 0) perfect++;
    else if (missing_rate <= 0.001) very_low++;
    else if (missing_rate <= 0.005) low++;
    else if (missing_rate <= 0.01) moderate++;
    else if (missing_rate <= 0.05) high++;
    else very_high++;
}
END {
    printf "  Missing rate = 0:          %'"'"'d (%.2f%%)\n", perfect, (perfect/NR)*100;
    printf "  0 < missing rate ≤ 0.001:  %'"'"'d (%.2f%%)\n", very_low, (very_low/NR)*100;
    printf "  0.001 < missing rate ≤ 0.005: %'"'"'d (%.2f%%)\n", low, (low/NR)*100;
    printf "  0.005 < missing rate ≤ 0.01: %'"'"'d (%.2f%%)\n", moderate, (moderate/NR)*100;
    printf "  0.01 < missing rate ≤ 0.05: %'"'"'d (%.2f%%)\n", high, (high/NR)*100;
    printf "  Missing rate > 0.05:       %'"'"'d (%.2f%%)\n", very_high, (very_high/NR)*100;
}' ${PREFIX}.vmiss >> variant_missing_qc.log
echo "" >> variant_missing_qc.log

# 统计将被过滤的变异数（missing rate > 0.01，注意是严格大于）
# PLINK2 --geno 使用 > 而非 >=，即保留 missing_rate <= threshold 的变异
to_remove=$(awk 'NR>1 && $5 > 0.01' ${PREFIX}.vmiss | wc -l)
printf "Variants with missing rate > 0.01: %'d (%.2f%%)\n" ${to_remove} $(awk "BEGIN {printf \"%.2f\", (${to_remove}/${total_variants})*100}") >> variant_missing_qc.log
echo "" >> variant_missing_qc.log

# PLINK2性染色体过滤说明
echo "Note on sex chromosome filtering:" >> variant_missing_qc.log
echo "  PLINK2 --geno filter handles sex chromosomes intelligently:" >> variant_missing_qc.log
echo "  - Autosomes (chr1-22): Missing rate calculated across all samples" >> variant_missing_qc.log
echo "  - X chromosome: Missing rate calculated only for female samples" >> variant_missing_qc.log
echo "  - Y chromosome: Missing rate calculated only for male samples" >> variant_missing_qc.log
echo "  - MT (mitochondrial): Missing rate calculated across all samples" >> variant_missing_qc.log
echo "  This prevents incorrect filtering due to sex-specific biology." >> variant_missing_qc.log
echo "" >> variant_missing_qc.log

# 过滤缺失率 > 0.01 的变异
echo "Step 2: Filtering variants with missing rate > 0.01..." >> variant_missing_qc.log
${PLINK2_PATH} \
    --bfile ${PREFIX} \
    --geno 0.01 \
    --make-bed \
    --out ${PREFIX}.vmiss_qc \
    --threads 8

# 统计过滤后的结果
output_variant_count=$(wc -l < ${PREFIX}.vmiss_qc.bim)
removed_count=$((input_variant_count - output_variant_count))

printf "Variants after filtering: %'d\n" ${output_variant_count} >> variant_missing_qc.log
printf "Variants removed: %'d\n" ${removed_count} >> variant_missing_qc.log
removal_rate=$(awk "BEGIN {printf \"%.2f%%\", (${removed_count}/${input_variant_count})*100}")
echo "Removal rate: ${removal_rate}" >> variant_missing_qc.log
retention_rate=$(awk "BEGIN {printf \"%.2f%%\", (${output_variant_count}/${input_variant_count})*100}")
echo "Retention rate: ${retention_rate}" >> variant_missing_qc.log
echo "" >> variant_missing_qc.log

# 统计每条染色体过滤后的变异数
echo "Variants per chromosome (after filtering):" >> variant_missing_qc.log
cut -f1 ${PREFIX}.vmiss_qc.bim | sort | uniq -c | awk '{printf "  Chr %s: %'"'"'d variants\n", $2, $1}' >> variant_missing_qc.log
echo "" >> variant_missing_qc.log

# 详细统计每条染色体的过滤情况
echo "Filtering details by chromosome:" >> variant_missing_qc.log
for chr in $(cut -f1 ${PREFIX}.bim | sort -u); do
    before=$(awk -v c="$chr" '$1==c' ${PREFIX}.bim | wc -l)
    after=$(awk -v c="$chr" '$1==c' ${PREFIX}.vmiss_qc.bim | wc -l)
    removed=$((before - after))
    if [ ${before} -gt 0 ]; then
        pct=$(awk "BEGIN {printf \"%.2f\", (${removed}/${before})*100}")
        printf "  Chr %s: %'d → %'d (removed: %'d, %.2f%%)\n" "$chr" ${before} ${after} ${removed} ${pct} >> variant_missing_qc.log
    fi
done
echo "" >> variant_missing_qc.log

# 最终统计汇总
echo "======================================" >> variant_missing_qc.log
echo "SUMMARY" >> variant_missing_qc.log
echo "======================================" >> variant_missing_qc.log
printf "Input variants: %'d\n" ${input_variant_count} >> variant_missing_qc.log
printf "Output variants: %'d\n" ${output_variant_count} >> variant_missing_qc.log
printf "Variants removed: %'d\n" ${removed_count} >> variant_missing_qc.log
printf "Samples (unchanged): %'d\n" ${input_sample_count} >> variant_missing_qc.log
echo "" >> variant_missing_qc.log
echo "Quality control threshold:" >> variant_missing_qc.log
echo "  Maximum allowed missing rate: 1% (0.01, strict inequality)" >> variant_missing_qc.log
echo "  Filter behavior: PLINK2 --geno uses > (not >=)" >> variant_missing_qc.log
echo "    - Variants with missing_rate > 0.01 are EXCLUDED" >> variant_missing_qc.log
echo "    - Variants with missing_rate <= 0.01 are RETAINED" >> variant_missing_qc.log
echo "  Rationale: Variants with >1% missing genotypes may indicate" >> variant_missing_qc.log
echo "             technical issues or poor probe performance" >> variant_missing_qc.log
echo "" >> variant_missing_qc.log


echo "Output files:" >> variant_missing_qc.log
echo "  - ${PREFIX}.vmiss_qc.bed (binary genotype file, QC-passed)" >> variant_missing_qc.log
echo "  - ${PREFIX}.vmiss_qc.bim (variant information, QC-passed)" >> variant_missing_qc.log
echo "  - ${PREFIX}.vmiss_qc.fam (sample information, unchanged)" >> variant_missing_qc.log
echo "  - ${PREFIX}.vmiss (variant missing rate statistics)" >> variant_missing_qc.log
echo "  - ${PREFIX}.vmiss_qc.log (PLINK2 log file)" >> variant_missing_qc.log
echo "  - variant_missing_qc.log (this report)" >> variant_missing_qc.log
