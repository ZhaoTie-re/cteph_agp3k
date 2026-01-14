#!/bin/bash
# 05_merge_plink.sh
# Usage: 05_merge_plink.sh <plink2_path> <plink_path>

PLINK2_PATH=$1
PLINK_PATH=$2

# 提取所有.bim文件的路径
bim_files=$(ls -1 *.bim)

echo "=== PLINK File Merging Report ===" > merge.log
echo "Date: $(date)" >> merge.log
echo "" >> merge.log

# 统计每个文件的变体数
echo "Input files:" >> merge.log
file_count=0
for bim in ${bim_files}; do
    file_count=$((file_count + 1))
    prefix=${bim%.bim}
    variant_count=$(wc -l < ${bim})
    sample_count=$(wc -l < ${prefix}.fam)
    printf "  %d. %s: %'d variants, %'d samples\n" ${file_count} ${prefix} ${variant_count} ${sample_count} >> merge.log
done
echo "" >> merge.log

# 提取每个.bim文件的第2列（variant ID）
echo "Extracting variant IDs from all files..." >> merge.log
for bim in ${bim_files}; do
    awk '{print $2}' ${bim} | sort > ${bim}.variants.txt
done

# 找到所有文件共有的变体（交集）
echo "Finding common variants across all files..." >> merge.log

# 获取第一个文件作为起点
first_bim=$(echo "${bim_files}" | head -1)
cp ${first_bim}.variants.txt common_variants.txt

# 依次与其他文件求交集
for bim in ${bim_files}; do
    if [ "${bim}" != "${first_bim}" ]; then
        comm -12 common_variants.txt ${bim}.variants.txt > temp_common.txt
        mv temp_common.txt common_variants.txt
    fi
done

common_count=$(wc -l < common_variants.txt)
printf "Common variants found: %'d\n" ${common_count} >> merge.log
echo "" >> merge.log

# 为每个PLINK文件提取共同变体
echo "Extracting common variants from each file..." >> merge.log
first_filtered=""
file_num=0
for bim in ${bim_files}; do
    prefix=${bim%.bim}
    
    ${PLINK2_PATH} \
        --bfile ${prefix} \
        --extract common_variants.txt \
        --make-bed \
        --out ${prefix}.common \
        --threads 4
    
    # 验证提取后的变体数
    extracted_count=$(wc -l < ${prefix}.common.bim)
    printf "  %s: %'d variants extracted\n" ${prefix} ${extracted_count} >> merge.log
    
    # 记录第一个文件，其他文件添加到merge列表
    file_num=$((file_num + 1))
    if [ ${file_num} -eq 1 ]; then
        first_filtered=${prefix}.common
    else
        echo "${prefix}.common" >> merge_list.txt
    fi
done
echo "" >> merge.log

# 检查文件数量并执行相应操作
if [ ${file_num} -gt 1 ]; then
    # 多个文件，执行合并
    echo "Merging ${file_num} PLINK files using PLINK 1.9..." >> merge.log
    echo "Base file: ${first_filtered}" >> merge.log
    echo "Files to merge with base:" >> merge.log
    cat merge_list.txt >> merge.log
    echo "" >> merge.log
    
    ${PLINK_PATH} \
        --bfile ${first_filtered} \
        --merge-list merge_list.txt \
        --keep-allele-order \
        --allow-extra-chr \
        --make-bed \
        --out cteph_agp3k.ajsa.sqc.norm.plink1
    
    echo "Merge completed successfully with PLINK 1.9" >> merge.log
    echo "" >> merge.log
    
    # 使用PLINK2转换为PLINK2格式
    echo "Converting merged file to PLINK2 format..." >> merge.log
    ${PLINK2_PATH} \
        --bfile cteph_agp3k.ajsa.sqc.norm.plink1 \
        --make-bed \
        --out cteph_agp3k.ajsa.sqc.norm \
        --threads 8
    
    echo "Conversion to PLINK2 format completed" >> merge.log
    
    # 清理PLINK1中间文件
    rm -f cteph_agp3k.ajsa.sqc.norm.plink1.bed
    rm -f cteph_agp3k.ajsa.sqc.norm.plink1.bim
    rm -f cteph_agp3k.ajsa.sqc.norm.plink1.fam
    rm -f cteph_agp3k.ajsa.sqc.norm.plink1.log
    
elif [ ${file_num} -eq 1 ]; then
    # 只有一个文件，使用PLINK2转换
    echo "Only one file found, converting to PLINK2 format..." >> merge.log
    ${PLINK2_PATH} \
        --bfile ${first_filtered} \
        --make-bed \
        --out cteph_agp3k.ajsa.sqc.norm \
        --threads 8
    
    echo "Conversion to PLINK2 format completed" >> merge.log
else
    echo "ERROR: No files found to merge!" >> merge.log
    exit 1
fi

# 最终统计
echo "" >> merge.log
echo "======================================" >> merge.log
echo "MERGE SUMMARY" >> merge.log
echo "======================================" >> merge.log

final_variant_count=$(wc -l < cteph_agp3k.ajsa.sqc.norm.bim)
final_sample_count=$(wc -l < cteph_agp3k.ajsa.sqc.norm.fam)

printf "Input files: %'d\n" ${file_count} >> merge.log
printf "Common variants: %'d\n" ${common_count} >> merge.log
printf "Final variants in merged file: %'d\n" ${final_variant_count} >> merge.log
printf "Total samples in merged file: %'d\n" ${final_sample_count} >> merge.log
echo "" >> merge.log

# 统计合并后的性别分布
echo "Sex distribution in merged file:" >> merge.log
male_count=$(awk '$5==1' cteph_agp3k.ajsa.sqc.norm.fam | wc -l)
female_count=$(awk '$5==2' cteph_agp3k.ajsa.sqc.norm.fam | wc -l)
unknown_count=$(awk '$5==0' cteph_agp3k.ajsa.sqc.norm.fam | wc -l)
printf "  Males: %'d\n" ${male_count} >> merge.log
printf "  Females: %'d\n" ${female_count} >> merge.log
printf "  Unknown: %'d\n" ${unknown_count} >> merge.log
echo "" >> merge.log

# 统计合并后的表型分布（注意：此时表型尚未更新）
echo "Phenotype distribution in merged file (before update):" >> merge.log
case_count=$(awk '$6==2' cteph_agp3k.ajsa.sqc.norm.fam | wc -l)
ctrl_count=$(awk '$6==1' cteph_agp3k.ajsa.sqc.norm.fam | wc -l)
missing_pheno=$(awk '$6==-9 || $6==0' cteph_agp3k.ajsa.sqc.norm.fam | wc -l)
printf "  Cases (pheno=2): %'d\n" ${case_count} >> merge.log
printf "  Controls (pheno=1): %'d\n" ${ctrl_count} >> merge.log
printf "  Missing/Unknown: %'d\n" ${missing_pheno} >> merge.log
echo "  Note: Phenotype will be updated in the next process step" >> merge.log
echo "" >> merge.log

echo "Output files:" >> merge.log
echo "  - cteph_agp3k.ajsa.sqc.norm.bed (merged binary genotype file)" >> merge.log
echo "  - cteph_agp3k.ajsa.sqc.norm.bim (merged variant information)" >> merge.log
echo "  - cteph_agp3k.ajsa.sqc.norm.fam (merged sample information)" >> merge.log
echo "  - merge.log (this log file)" >> merge.log

# 清理临时文件
rm -f *.common.* *.variants.txt common_variants.txt merge_list.txt
