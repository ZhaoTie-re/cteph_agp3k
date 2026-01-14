#!/bin/bash
# 10_variant_qc_filter.sh
# Usage: 10_variant_qc_filter.sh <plink2_path> <region> <prefix> <tsv> <maf_mode> <maf_threshold> <hwe_mode> <hwe_ctrl_threshold> <hwe_case_threshold>

PLINK2_PATH=$1
REGION=$2
PREFIX=$3
TSV=$4
MAF_MODE=$5
MAF_THRESHOLD=$6
HWE_MODE=$7
HWE_CTRL_THRESHOLD=$8
HWE_CASE_THRESHOLD=$9

OUTPUT_PREFIX="cteph_agp3k.array.sqc.vqc"

echo "=== Variant Filtering by MAF and HWE Report ===" > variant_filter.log
echo "Date: $(date)" >> variant_filter.log
echo "Region: ${REGION}" >> variant_filter.log
echo "" >> variant_filter.log

# 统计输入文件信息
input_variant_count=$(wc -l < ${PREFIX}.bim)
input_sample_count=$(wc -l < ${PREFIX}.fam)

printf "Input file: %s\n" "${PREFIX}" >> variant_filter.log
printf "Total variants: %'d\n" ${input_variant_count} >> variant_filter.log
printf "Total samples: %'d\n" ${input_sample_count} >> variant_filter.log
echo "" >> variant_filter.log

# 过滤条件
echo "Filtering criteria:" >> variant_filter.log
echo "  MAF filter mode: ${MAF_MODE}" >> variant_filter.log
echo "  MAF threshold: MAF_${MAF_MODE} < ${MAF_THRESHOLD}" >> variant_filter.log
echo "  HWE filter mode: ${HWE_MODE}" >> variant_filter.log
if [[ "${HWE_MODE}" == *"CTRL"* ]]; then
    echo "  HWE CTRL threshold: HWE_CTRL < ${HWE_CTRL_THRESHOLD}" >> variant_filter.log
fi
if [[ "${HWE_MODE}" == *"CASE"* ]]; then
    echo "  HWE CASE threshold: HWE_CASE < ${HWE_CASE_THRESHOLD}" >> variant_filter.log
fi
echo "" >> variant_filter.log

# 使用Python提取需要移除的变异ID
echo "Step 1: Identifying variants to remove..." >> variant_filter.log

python3 << EOF 2>&1 | tee -a variant_filter.log
import pandas as pd
import sys

print("Reading variant QC summary file...")

try:
    # 读取TSV文件
    df = pd.read_csv("${TSV}", sep="\t")
    print(f"  Total variants in TSV: {len(df)}")
    
    # 检查所需的列是否存在
    maf_col = "MAF_${MAF_MODE}"
    hwe_mode = "${HWE_MODE}"
    
    if maf_col not in df.columns:
        print(f"Error: Column '{maf_col}' not found in TSV file")
        print(f"Available columns: {', '.join(df.columns)}")
        sys.exit(1)
    
    # 检查HWE相关列
    use_hwe_ctrl = "CTRL" in hwe_mode
    use_hwe_case = "CASE" in hwe_mode
    
    if use_hwe_ctrl and "HWE_CTRL" not in df.columns:
        print(f"Error: Column 'HWE_CTRL' not found in TSV file")
        print(f"Available columns: {', '.join(df.columns)}")
        sys.exit(1)
    
    if use_hwe_case and "HWE_CASE" not in df.columns:
        print(f"Error: Column 'HWE_CASE' not found in TSV file")
        print(f"Available columns: {', '.join(df.columns)}")
        sys.exit(1)
    
    # 应用过滤条件
    maf_threshold = float("${MAF_THRESHOLD}")
    hwe_ctrl_threshold = float("${HWE_CTRL_THRESHOLD}")
    hwe_case_threshold = float("${HWE_CASE_THRESHOLD}")
    
    # MAF过滤
    maf_fail = df[maf_col] < maf_threshold
    
    # HWE过滤（根据模式）
    hwe_fail = pd.Series([False] * len(df), index=df.index)
    hwe_ctrl_fail = pd.Series([False] * len(df), index=df.index)
    hwe_case_fail = pd.Series([False] * len(df), index=df.index)
    
    if use_hwe_ctrl:
        hwe_ctrl_fail = df["HWE_CTRL"] < hwe_ctrl_threshold
        hwe_fail = hwe_fail | hwe_ctrl_fail
    
    if use_hwe_case:
        hwe_case_fail = df["HWE_CASE"] < hwe_case_threshold
        hwe_fail = hwe_fail | hwe_case_fail
    
    # 合并所有过滤条件
    to_remove = df[maf_fail | hwe_fail]
    
    print(f"\nFiltering results:")
    print(f"  Variants failing MAF filter ({maf_col} < {maf_threshold}): {maf_fail.sum()} ({maf_fail.sum()/len(df)*100:.2f}%)")
    
    if use_hwe_ctrl:
        print(f"  Variants failing HWE_CTRL filter (< {hwe_ctrl_threshold}): {hwe_ctrl_fail.sum()} ({hwe_ctrl_fail.sum()/len(df)*100:.2f}%)")
    
    if use_hwe_case:
        print(f"  Variants failing HWE_CASE filter (< {hwe_case_threshold}): {hwe_case_fail.sum()} ({hwe_case_fail.sum()/len(df)*100:.2f}%)")
    
    print(f"  Variants failing any filter (to remove): {len(to_remove)} ({len(to_remove)/len(df)*100:.2f}%)")
    print(f"  Variants passing all filters (to keep): {len(df) - len(to_remove)} ({(len(df)-len(to_remove))/len(df)*100:.2f}%)")
    
    # 保存需要移除的变异ID（不包含列名）
    to_remove["VARIANT_ID"].to_csv("${PREFIX}.variants_to_remove.txt", index=False, header=False)
    
    print(f"\nVariant IDs to remove saved to: ${PREFIX}.variants_to_remove.txt")
    
    # 统计失败组合
    print(f"\nFailure breakdown:")
    
    if use_hwe_ctrl and use_hwe_case:
        # 三个过滤器都启用
        maf_only = maf_fail & ~hwe_ctrl_fail & ~hwe_case_fail
        hwe_ctrl_only = ~maf_fail & hwe_ctrl_fail & ~hwe_case_fail
        hwe_case_only = ~maf_fail & ~hwe_ctrl_fail & hwe_case_fail
        maf_hwe_ctrl = maf_fail & hwe_ctrl_fail & ~hwe_case_fail
        maf_hwe_case = maf_fail & ~hwe_ctrl_fail & hwe_case_fail
        hwe_both = ~maf_fail & hwe_ctrl_fail & hwe_case_fail
        all_three = maf_fail & hwe_ctrl_fail & hwe_case_fail
        
        print(f"  MAF only: {maf_only.sum()} variants ({maf_only.sum()/len(df)*100:.2f}%)")
        print(f"  HWE_CTRL only: {hwe_ctrl_only.sum()} variants ({hwe_ctrl_only.sum()/len(df)*100:.2f}%)")
        print(f"  HWE_CASE only: {hwe_case_only.sum()} variants ({hwe_case_only.sum()/len(df)*100:.2f}%)")
        print(f"  MAF + HWE_CTRL: {maf_hwe_ctrl.sum()} variants ({maf_hwe_ctrl.sum()/len(df)*100:.2f}%)")
        print(f"  MAF + HWE_CASE: {maf_hwe_case.sum()} variants ({maf_hwe_case.sum()/len(df)*100:.2f}%)")
        print(f"  HWE_CTRL + HWE_CASE: {hwe_both.sum()} variants ({hwe_both.sum()/len(df)*100:.2f}%)")
        print(f"  All three filters: {all_three.sum()} variants ({all_three.sum()/len(df)*100:.2f}%)")
    
    elif use_hwe_ctrl:
        # 只有MAF和HWE_CTRL
        maf_only = maf_fail & ~hwe_ctrl_fail
        hwe_ctrl_only = ~maf_fail & hwe_ctrl_fail
        both = maf_fail & hwe_ctrl_fail
        
        print(f"  MAF only: {maf_only.sum()} variants ({maf_only.sum()/len(df)*100:.2f}%)")
        print(f"  HWE_CTRL only: {hwe_ctrl_only.sum()} variants ({hwe_ctrl_only.sum()/len(df)*100:.2f}%)")
        print(f"  Both MAF + HWE_CTRL: {both.sum()} variants ({both.sum()/len(df)*100:.2f}%)")
    
    elif use_hwe_case:
        # 只有MAF和HWE_CASE
        maf_only = maf_fail & ~hwe_case_fail
        hwe_case_only = ~maf_fail & hwe_case_fail
        both = maf_fail & hwe_case_fail
        
        print(f"  MAF only: {maf_only.sum()} variants ({maf_only.sum()/len(df)*100:.2f}%)")
        print(f"  HWE_CASE only: {hwe_case_only.sum()} variants ({hwe_case_only.sum()/len(df)*100:.2f}%)")
        print(f"  Both MAF + HWE_CASE: {both.sum()} variants ({both.sum()/len(df)*100:.2f}%)")
    
    else:
        # 只有MAF过滤
        print(f"  MAF only: {maf_fail.sum()} variants ({maf_fail.sum()/len(df)*100:.2f}%)")
    
except Exception as e:
    print(f"Error processing TSV file: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)
EOF

echo "" >> variant_filter.log

# 统计要移除的变异数
variants_to_remove=$(wc -l < ${PREFIX}.variants_to_remove.txt)
printf "Variants identified for removal: %'d\n" ${variants_to_remove} >> variant_filter.log
echo "" >> variant_filter.log

# 使用PLINK2移除这些变异
if [ ${variants_to_remove} -gt 0 ]; then
    echo "Step 2: Removing filtered variants with PLINK2..." >> variant_filter.log
    
    ${PLINK2_PATH} \
        --bfile ${PREFIX} \
        --exclude ${PREFIX}.variants_to_remove.txt \
        --make-bed \
        --out ${OUTPUT_PREFIX} \
        --threads 8
    
    # 统计过滤后的结果
    output_variant_count=$(wc -l < ${OUTPUT_PREFIX}.bim)
    removed_count=$((input_variant_count - output_variant_count))
    
    echo "" >> variant_filter.log
    printf "Variants after filtering: %'d\n" ${output_variant_count} >> variant_filter.log
    printf "Variants removed: %'d\n" ${removed_count} >> variant_filter.log
    removal_rate=$(awk "BEGIN {printf \"%.2f%%\", (${removed_count}/${input_variant_count})*100}")
    echo "Removal rate: ${removal_rate}" >> variant_filter.log
    retention_rate=$(awk "BEGIN {printf \"%.2f%%\", (${output_variant_count}/${input_variant_count})*100}")
    echo "Retention rate: ${retention_rate}" >> variant_filter.log
    
else
    echo "Step 2: No variants to remove (all passed filters)" >> variant_filter.log
    echo "  Copying input files to output..." >> variant_filter.log
    
    cp ${PREFIX}.bed ${OUTPUT_PREFIX}.bed
    cp ${PREFIX}.bim ${OUTPUT_PREFIX}.bim
    cp ${PREFIX}.fam ${OUTPUT_PREFIX}.fam
    
    output_variant_count=${input_variant_count}
    echo "  All variants retained: ${output_variant_count}" >> variant_filter.log
fi

echo "" >> variant_filter.log

# ===== 最终统计汇总 =====
echo "======================================" >> variant_filter.log
echo "SUMMARY" >> variant_filter.log
echo "======================================" >> variant_filter.log

printf "Region: %s\n" "${REGION}" >> variant_filter.log
printf "Input variants: %'d\n" ${input_variant_count} >> variant_filter.log
printf "Output variants: %'d\n" ${output_variant_count} >> variant_filter.log

if [ ${variants_to_remove} -gt 0 ]; then
    removed_count=$((input_variant_count - output_variant_count))
    printf "Variants removed: %'d\n" ${removed_count} >> variant_filter.log
    removal_rate=$(awk "BEGIN {printf \"%.2f%%\", (${removed_count}/${input_variant_count})*100}")
    retention_rate=$(awk "BEGIN {printf \"%.2f%%\", (${output_variant_count}/${input_variant_count})*100}")
    echo "  Removal rate: ${removal_rate}" >> variant_filter.log
    echo "  Retention rate: ${retention_rate}" >> variant_filter.log
else
    echo "Variants removed: 0 (all passed filters)" >> variant_filter.log
fi

printf "Samples (unchanged): %'d\n" ${input_sample_count} >> variant_filter.log
echo "" >> variant_filter.log

echo "Quality control filters applied:" >> variant_filter.log
echo "  1. MAF filter: MAF_${MAF_MODE} < ${MAF_THRESHOLD}" >> variant_filter.log
echo "     Rationale: Variants with very low minor allele frequency may have" >> variant_filter.log
echo "                insufficient statistical power for association testing" >> variant_filter.log
echo "" >> variant_filter.log

filter_num=2
if [[ "${HWE_MODE}" == *"CTRL"* ]]; then
    echo "  ${filter_num}. HWE_CTRL filter: HWE_CTRL < ${HWE_CTRL_THRESHOLD}" >> variant_filter.log
    echo "     Rationale: Significant deviation from Hardy-Weinberg equilibrium in" >> variant_filter.log
    echo "                controls may indicate genotyping errors or population stratification" >> variant_filter.log
    echo "     Note: Using controls avoids removing disease-associated variants" >> variant_filter.log
    echo "" >> variant_filter.log
    filter_num=$((filter_num + 1))
fi

if [[ "${HWE_MODE}" == *"CASE"* ]]; then
    echo "  ${filter_num}. HWE_CASE filter: HWE_CASE < ${HWE_CASE_THRESHOLD}" >> variant_filter.log
    echo "     Rationale: Extreme deviation from HWE in cases may indicate" >> variant_filter.log
    echo "                genotyping errors or technical artifacts" >> variant_filter.log
    echo "" >> variant_filter.log
fi

echo "Filter logic: Variants are REMOVED if they fail ANY criterion" >> variant_filter.log
echo "              (i.e., OR logic among all active filters)" >> variant_filter.log
echo "" >> variant_filter.log
echo "Active filter modes:" >> variant_filter.log
echo "  - MAF mode: ${MAF_MODE} (filtering based on MAF_${MAF_MODE})" >> variant_filter.log
echo "  - HWE mode: ${HWE_MODE} (filtering based on ${HWE_MODE} samples)" >> variant_filter.log
echo "" >> variant_filter.log

echo "Output files:" >> variant_filter.log
echo "  - ${OUTPUT_PREFIX}.bed (filtered binary genotype file)" >> variant_filter.log
echo "  - ${OUTPUT_PREFIX}.bim (filtered variant information)" >> variant_filter.log
echo "  - ${OUTPUT_PREFIX}.fam (sample information, unchanged)" >> variant_filter.log
echo "  - ${PREFIX}.variants_to_remove.txt (list of removed variant IDs)" >> variant_filter.log
echo "  - variant_filter.log (this report)" >> variant_filter.log
