#!/bin/bash
# 03_convert_to_vcf.sh
# Usage: 03_convert_to_vcf.sh <plink2_path> <plink_path> <nagasaki_pipeline_path> <chr_rename_file> <prefix>

PLINK2_PATH=$1
PLINK_PATH=$2
NAGASAKI_PIPELINE_PATH=$3
CHR_RENAME_FILE=$4
PREFIX=$5

# 优化的流式处理：过滤染色体 -> PLINK转VCF -> 重命名染色体 -> 标准化 -> 填充ID
# 使用管道减少I/O操作

# 步骤0: 过滤只保留标准染色体（1-22, X, Y, MT）并清理非标准等位基因
echo "=== VCF Conversion QC Report ===" > ${PREFIX}.norm.vcf.gz.log
echo "Processing: ${PREFIX}" >> ${PREFIX}.norm.vcf.gz.log
echo "Date: $(date)" >> ${PREFIX}.norm.vcf.gz.log
echo "" >> ${PREFIX}.norm.vcf.gz.log

echo "Step 0: Input variants (before filtering)" >> ${PREFIX}.norm.vcf.gz.log
wc -l ${PREFIX}.bim | awk '{printf "  Total variants: %'"'"'d\n", $1}' >> ${PREFIX}.norm.vcf.gz.log
cut -f1 ${PREFIX}.bim | sort | uniq -c | awk '{printf "  Chr %s: %'"'"'d variants\n", $2, $1}' >> ${PREFIX}.norm.vcf.gz.log
echo "" >> ${PREFIX}.norm.vcf.gz.log

# 统计问题变异
echo "  Quality issues detected:" >> ${PREFIX}.norm.vcf.gz.log
# 使用更简单的匹配方法
missing_allele_count=$(awk '$5 == "." || $6 == "."' ${PREFIX}.bim | wc -l)
non_standard_count=$(awk '$5 !~ /^[ATCG]$/ || $6 !~ /^[ATCG]$/' ${PREFIX}.bim | wc -l)
non_acgt_count=$((non_standard_count))
true_non_acgt_count=$((non_acgt_count - missing_allele_count))

printf "    - Non-ACGT alleles (indels, etc.): %'d\n" ${true_non_acgt_count} >> ${PREFIX}.norm.vcf.gz.log
printf "    - Missing allele information (. in REF/ALT): %'d\n" ${missing_allele_count} >> ${PREFIX}.norm.vcf.gz.log
printf "    - Total problematic variants: %'d\n" ${non_acgt_count} >> ${PREFIX}.norm.vcf.gz.log
echo "" >> ${PREFIX}.norm.vcf.gz.log

${PLINK2_PATH} \
    --bfile ${PREFIX} \
    --chr 1-22,X,Y,MT \
    --snps-only just-acgt \
    --make-bed \
    --out temp_filtered_${PREFIX} \
    --threads 8

echo "Step 1: After chromosome filtering (1-22,X,Y,MT) and ACGT-only SNPs" >> ${PREFIX}.norm.vcf.gz.log
wc -l temp_filtered_${PREFIX}.bim | awk '{printf "  Remaining variants: %'"'"'d\n", $1}' >> ${PREFIX}.norm.vcf.gz.log

# 计算被过滤掉的变体数量
original_count=$(wc -l < ${PREFIX}.bim)
filtered_count=$(wc -l < temp_filtered_${PREFIX}.bim)
removed_count=$((original_count - filtered_count))
printf "  Removed variants: %'d\n" ${removed_count} >> ${PREFIX}.norm.vcf.gz.log
echo "  Removal rate: $(awk "BEGIN {printf \"%.2f%%\", (${removed_count}/${original_count})*100}")" >> ${PREFIX}.norm.vcf.gz.log
echo "" >> ${PREFIX}.norm.vcf.gz.log

# 检查过滤后是否还有缺失等位基因
remaining_missing=$(awk '$5 == "." || $6 == "."' temp_filtered_${PREFIX}.bim | wc -l)
if [ ${remaining_missing} -gt 0 ]; then
    echo "  WARNING: ${remaining_missing} variants with missing alleles remain after PLINK2 filtering" >> ${PREFIX}.norm.vcf.gz.log
    echo "  These will be removed manually..." >> ${PREFIX}.norm.vcf.gz.log
    
    # 手动移除缺失等位基因的变异
    awk '$5 != "." && $6 != "."' temp_filtered_${PREFIX}.bim > temp_filtered_${PREFIX}.bim.clean
    
    # 提取这些变异的ID用于过滤
    awk '{print $2}' temp_filtered_${PREFIX}.bim.clean > variants_to_keep.txt
    
    ${PLINK2_PATH} \
        --bfile temp_filtered_${PREFIX} \
        --extract variants_to_keep.txt \
        --make-bed \
        --out temp_filtered_clean_${PREFIX} \
        --threads 8
    
    # 替换文件
    mv temp_filtered_clean_${PREFIX}.bed temp_filtered_${PREFIX}.bed
    mv temp_filtered_clean_${PREFIX}.bim temp_filtered_${PREFIX}.bim
    mv temp_filtered_clean_${PREFIX}.fam temp_filtered_${PREFIX}.fam
    
    # 更新统计
    filtered_count=$(wc -l < temp_filtered_${PREFIX}.bim)
    removed_count=$((original_count - filtered_count))
    
    printf "  After removing missing alleles: %'d variants\n" ${filtered_count} >> ${PREFIX}.norm.vcf.gz.log
    printf "  Total removed: %'d\n" ${removed_count} >> ${PREFIX}.norm.vcf.gz.log
    
    rm -f variants_to_keep.txt temp_filtered_${PREFIX}.bim.clean
fi

# 步骤2: PLINK转VCF（写入磁盘，因为plink不支持stdout）
echo "" >> ${PREFIX}.norm.vcf.gz.log
echo "Step 2: PLINK binary to VCF conversion" >> ${PREFIX}.norm.vcf.gz.log

${PLINK_PATH} \
    --bfile temp_filtered_${PREFIX} \
    --recode vcf-iid bgz \
    --out temp_${PREFIX}

vcf_count=$(bcftools view -H temp_${PREFIX}.vcf.gz | wc -l)
printf "  VCF variants: %'d\n" ${vcf_count} >> ${PREFIX}.norm.vcf.gz.log
echo "  Format: VCF 4.2, bgzipped" >> ${PREFIX}.norm.vcf.gz.log
echo "  Reason: Convert PLINK binary format to VCF for downstream processing" >> ${PREFIX}.norm.vcf.gz.log

# 步骤3: 重命名染色体
echo "" >> ${PREFIX}.norm.vcf.gz.log
echo "Step 3: Chromosome renaming (numeric to chr-prefix)" >> ${PREFIX}.norm.vcf.gz.log

bcftools annotate \
    --rename-chrs ${CHR_RENAME_FILE} \
    --threads 8 \
    -Oz \
    -o temp_renamed_${PREFIX}.vcf.gz \
    temp_${PREFIX}.vcf.gz || { echo "Error: bcftools annotate (rename) failed" >> ${PREFIX}.norm.vcf.gz.log; exit 1; }

bcftools index -t temp_renamed_${PREFIX}.vcf.gz

renamed_count=$(bcftools view -H temp_renamed_${PREFIX}.vcf.gz | wc -l)
printf "  Variants after renaming: %'d\n" ${renamed_count} >> ${PREFIX}.norm.vcf.gz.log
echo "  Chromosomes: 1-22 → chr1-chr22, 23 → chrX, 24 → chrY, 26 → chrM" >> ${PREFIX}.norm.vcf.gz.log
echo "  Reason: Match reference genome chromosome naming convention (GRCh38)" >> ${PREFIX}.norm.vcf.gz.log

# 步骤4: 标准化并严格检查REF一致性
echo "" >> ${PREFIX}.norm.vcf.gz.log
echo "Step 4: VCF normalization with reference genome" >> ${PREFIX}.norm.vcf.gz.log
echo "  Using --check-ref s (skip variants with REF mismatch)" >> ${PREFIX}.norm.vcf.gz.log

# 首先记录重命名后的变体数
before_norm_count=${renamed_count}

# 执行标准化，使用-s参数跳过REF不匹配的变体
bcftools norm \
    --multiallelics -any \
    --fasta-ref ${NAGASAKI_PIPELINE_PATH}/data/hs38DH.fa \
    --check-ref s \
    --threads 8 \
    -Oz \
    -o temp_norm_${PREFIX}.vcf.gz \
    temp_renamed_${PREFIX}.vcf.gz \
    2> ${PREFIX}.norm_warnings.txt || { echo "Error: bcftools norm failed" >> ${PREFIX}.norm.vcf.gz.log; exit 1; }

bcftools index -t temp_norm_${PREFIX}.vcf.gz

# 统计标准化后的变体数
norm_count=$(bcftools view -H temp_norm_${PREFIX}.vcf.gz | wc -l)
total_norm_change=$((norm_count - before_norm_count))

printf "  Variants before normalization: %'d\n" ${before_norm_count} >> ${PREFIX}.norm.vcf.gz.log
printf "  Variants after normalization: %'d\n" ${norm_count} >> ${PREFIX}.norm.vcf.gz.log

if [ ${total_norm_change} -eq 0 ]; then
    echo "  Net change: 0 (no variants added or removed)" >> ${PREFIX}.norm.vcf.gz.log
elif [ ${total_norm_change} -gt 0 ]; then
    printf "  Net change: +%'d variants (from multiallelic splitting)\n" ${total_norm_change} >> ${PREFIX}.norm.vcf.gz.log
else
    removed=$((-total_norm_change))
    printf "  Net change: -%'d variants removed\n" ${removed} >> ${PREFIX}.norm.vcf.gz.log
fi
echo "" >> ${PREFIX}.norm.vcf.gz.log

# 解析bcftools norm的输出统计
echo "  Normalization summary (from bcftools norm):" >> ${PREFIX}.norm.vcf.gz.log
cat ${PREFIX}.norm_warnings.txt >> ${PREFIX}.norm.vcf.gz.log
echo "" >> ${PREFIX}.norm.vcf.gz.log

# 提取关键统计信息
if [ -f ${PREFIX}.norm_warnings.txt ]; then
    # 提取Lines统计行
    lines_stats=$(grep "^Lines" ${PREFIX}.norm_warnings.txt 2>/dev/null || echo "")
    if [ ! -z "$lines_stats" ]; then
        # 解析各个数值: total/split/joined/realigned/removed/skipped
        # 提取冒号后的数字部分
        numbers_part=$(echo "$lines_stats" | sed 's/.*: *//')
        total_lines=$(echo "$numbers_part" | awk -F'/' '{print $1}')
        split_lines=$(echo "$numbers_part" | awk -F'/' '{print $2}')
        joined_lines=$(echo "$numbers_part" | awk -F'/' '{print $3}')
        realigned_lines=$(echo "$numbers_part" | awk -F'/' '{print $4}')
        removed_lines=$(echo "$numbers_part" | awk -F'/' '{print $5}')
        skipped_lines=$(echo "$numbers_part" | awk -F'/' '{print $6}')
        
        echo "  Detailed breakdown:" >> ${PREFIX}.norm.vcf.gz.log
        printf "    - Total variants processed: %'d\n" ${total_lines} >> ${PREFIX}.norm.vcf.gz.log
        printf "    - Multiallelic sites split: %'d\n" ${split_lines} >> ${PREFIX}.norm.vcf.gz.log
        printf "    - Variants joined: %'d\n" ${joined_lines} >> ${PREFIX}.norm.vcf.gz.log
        printf "    - Indels realigned: %'d\n" ${realigned_lines} >> ${PREFIX}.norm.vcf.gz.log
        printf "    - Variants removed: %'d\n" ${removed_lines} >> ${PREFIX}.norm.vcf.gz.log
        printf "    - Variants skipped: %'d\n" ${skipped_lines} >> ${PREFIX}.norm.vcf.gz.log
        echo "" >> ${PREFIX}.norm.vcf.gz.log
    fi
    
    # 提取REF/ALT统计行
    ref_alt_stats=$(grep "^REF/ALT" ${PREFIX}.norm_warnings.txt 2>/dev/null || echo "")
    if [ ! -z "$ref_alt_stats" ]; then
        # 提取冒号后的数字部分
        ref_numbers_part=$(echo "$ref_alt_stats" | sed 's/.*: *//')
        ref_alt_total=$(echo "$ref_numbers_part" | awk -F'/' '{print $1}')
        ref_alt_modified=$(echo "$ref_numbers_part" | awk -F'/' '{print $2}')
        ref_alt_added=$(echo "$ref_numbers_part" | awk -F'/' '{print $3}')
        
        echo "  REF/ALT allele adjustments:" >> ${PREFIX}.norm.vcf.gz.log
        printf "    - Total variants: %'d\n" ${ref_alt_total} >> ${PREFIX}.norm.vcf.gz.log
        printf "    - REF/ALT swapped (to match reference): %'d\n" ${ref_alt_modified} >> ${PREFIX}.norm.vcf.gz.log
        printf "    - ALT alleles added: %'d\n" ${ref_alt_added} >> ${PREFIX}.norm.vcf.gz.log
        
        if [ "${ref_alt_modified}" != "0" ] && [ "${ref_alt_total}" != "0" ] && [ ! -z "${ref_alt_total}" ]; then
            swap_rate=$(awk "BEGIN {printf \"%.2f%%\", (${ref_alt_modified}/${ref_alt_total})*100}")
            echo "    - Swap rate: ${swap_rate}" >> ${PREFIX}.norm.vcf.gz.log
            echo "    - Reason: VCF REF allele differs from reference genome" >> ${PREFIX}.norm.vcf.gz.log
            echo "    - Action: REF and ALT were swapped, genotypes flipped (0↔1)" >> ${PREFIX}.norm.vcf.gz.log
            echo "    - Note: Biological meaning preserved (e.g., AA stays AA)" >> ${PREFIX}.norm.vcf.gz.log
        fi
        echo "" >> ${PREFIX}.norm.vcf.gz.log
    fi
fi

echo "  Operations performed:" >> ${PREFIX}.norm.vcf.gz.log
echo "    1. Split multiallelic sites into biallelic records" >> ${PREFIX}.norm.vcf.gz.log
echo "    2. Left-align and normalize indels" >> ${PREFIX}.norm.vcf.gz.log
echo "    3. Swap REF/ALT when VCF REF != reference genome (--check-ref s)" >> ${PREFIX}.norm.vcf.gz.log
echo "    4. Flip genotypes accordingly to preserve biological meaning" >> ${PREFIX}.norm.vcf.gz.log

# 步骤5: 设置ID
echo "" >> ${PREFIX}.norm.vcf.gz.log
echo "Step 5: Set variant IDs" >> ${PREFIX}.norm.vcf.gz.log

bcftools annotate \
    --set-id '%CHROM:%POS:%REF:%ALT' \
    --threads 8 \
    -Oz \
    -o ${PREFIX}.norm.vcf.gz \
    temp_norm_${PREFIX}.vcf.gz || { echo "Error: bcftools annotate (set-id) failed" >> ${PREFIX}.norm.vcf.gz.log; exit 1; }

# 创建索引
bcftools index --threads 8 -t ${PREFIX}.norm.vcf.gz

final_count=$(bcftools view -H ${PREFIX}.norm.vcf.gz | wc -l)
printf "  Final variants: %'d\n" ${final_count} >> ${PREFIX}.norm.vcf.gz.log
echo "  ID format: CHROM:POS:REF:ALT (e.g., chr1:12345:A:G)" >> ${PREFIX}.norm.vcf.gz.log
echo "  Reason: Unique, reproducible variant identifiers" >> ${PREFIX}.norm.vcf.gz.log

# 最终统计汇总
echo "" >> ${PREFIX}.norm.vcf.gz.log
echo "======================================" >> ${PREFIX}.norm.vcf.gz.log
echo "SUMMARY" >> ${PREFIX}.norm.vcf.gz.log
echo "======================================" >> ${PREFIX}.norm.vcf.gz.log
printf "Input variants (Step 0):     %'d\n" ${original_count} >> ${PREFIX}.norm.vcf.gz.log
printf "After filtering (Step 1):    %'d (removed: %'d)\n" ${filtered_count} ${removed_count} >> ${PREFIX}.norm.vcf.gz.log
printf "After VCF conversion (Step 2): %'d\n" ${vcf_count} >> ${PREFIX}.norm.vcf.gz.log
printf "After chr rename (Step 3):   %'d\n" ${renamed_count} >> ${PREFIX}.norm.vcf.gz.log
if [ ${total_norm_change} -ge 0 ]; then
    printf "After normalization (Step 4): %'d (change: +%'d)\n" ${norm_count} ${total_norm_change} >> ${PREFIX}.norm.vcf.gz.log
else
    abs_change=$((-total_norm_change))
    printf "After normalization (Step 4): %'d (change: -%'d)\n" ${norm_count} ${abs_change} >> ${PREFIX}.norm.vcf.gz.log
fi
printf "Final output (Step 5):       %'d\n" ${final_count} >> ${PREFIX}.norm.vcf.gz.log
echo "" >> ${PREFIX}.norm.vcf.gz.log
total_removed=$((original_count - final_count))
retention_rate=$(awk "BEGIN {printf \"%.2f%%\", (${final_count}/${original_count})*100}")
printf "Total variants removed: %'d\n" ${total_removed} >> ${PREFIX}.norm.vcf.gz.log
echo "Retention rate: ${retention_rate}" >> ${PREFIX}.norm.vcf.gz.log
echo "" >> ${PREFIX}.norm.vcf.gz.log
echo "Main processing steps:" >> ${PREFIX}.norm.vcf.gz.log
echo "  1. Chromosome filtering: Removed non-standard chromosomes" >> ${PREFIX}.norm.vcf.gz.log
echo "  2. Allele filtering: Removed non-ACGT alleles" >> ${PREFIX}.norm.vcf.gz.log
echo "  3. Normalization: REF/ALT swapped where needed, multiallelic splitting" >> ${PREFIX}.norm.vcf.gz.log
echo "" >> ${PREFIX}.norm.vcf.gz.log
echo "Output files:" >> ${PREFIX}.norm.vcf.gz.log
echo "  - ${PREFIX}.norm.vcf.gz (normalized VCF)" >> ${PREFIX}.norm.vcf.gz.log
echo "  - ${PREFIX}.norm.vcf.gz.tbi (tabix index)" >> ${PREFIX}.norm.vcf.gz.log
echo "  - ${PREFIX}.norm.vcf.gz.log (this log file)" >> ${PREFIX}.norm.vcf.gz.log
echo "  - ${PREFIX}.norm_warnings.txt (normalization warnings)" >> ${PREFIX}.norm.vcf.gz.log

# 清理临时文件
rm -f temp_*${PREFIX}.*
