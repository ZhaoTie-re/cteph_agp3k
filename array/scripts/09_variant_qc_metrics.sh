#!/bin/bash
# 09_variant_qc_metrics.sh
# Usage: 09_variant_qc_metrics.sh <plink2_path> <region> <prefix>

PLINK2_PATH=$1
REGION=$2
PREFIX=$3

echo "=== Variant QC Statistics Calculation Report ===" > variant_qc_calculation.log
echo "Date: $(date)" >> variant_qc_calculation.log
echo "Region: ${REGION}" >> variant_qc_calculation.log
echo "" >> variant_qc_calculation.log

# 统计输入文件信息
input_variant_count=$(wc -l < ${PREFIX}.bim)
input_sample_count=$(wc -l < ${PREFIX}.fam)
case_count=$(awk '$6==2' ${PREFIX}.fam | wc -l)
ctrl_count=$(awk '$6==1' ${PREFIX}.fam | wc -l)
missing_pheno_count=$(awk '$6==-9 || $6==0' ${PREFIX}.fam | wc -l)

printf "Input file: %s\n" "${PREFIX}" >> variant_qc_calculation.log
printf "Total variants: %'d\n" ${input_variant_count} >> variant_qc_calculation.log
printf "Total samples: %'d\n" ${input_sample_count} >> variant_qc_calculation.log
printf "  Cases (pheno=2): %'d\n" ${case_count} >> variant_qc_calculation.log
printf "  Controls (pheno=1): %'d\n" ${ctrl_count} >> variant_qc_calculation.log
printf "  Missing phenotype: %'d\n" ${missing_pheno_count} >> variant_qc_calculation.log
echo "" >> variant_qc_calculation.log

# ===== 步骤1: 计算AAF (ALL) =====
echo "Step 1: Calculating ALT allele frequency (AAF) for ALL samples..." >> variant_qc_calculation.log
${PLINK2_PATH} \
    --bfile ${PREFIX} \
    --freq \
    --out ${PREFIX}.all \
    --threads 8

echo "  AAF (ALL) calculation completed" >> variant_qc_calculation.log
echo "" >> variant_qc_calculation.log

# ===== 步骤2: 计算AAF (CASE) =====
echo "Step 2: Calculating ALT allele frequency (AAF) for CASE samples..." >> variant_qc_calculation.log
${PLINK2_PATH} \
    --bfile ${PREFIX} \
    --keep <(awk '$6==2 {print $1, $2}' ${PREFIX}.fam) \
    --freq \
    --out ${PREFIX}.case \
    --threads 8

echo "  AAF (CASE) calculation completed" >> variant_qc_calculation.log
echo "" >> variant_qc_calculation.log

# ===== 步骤3: 计算AAF (CTRL) =====
echo "Step 3: Calculating ALT allele frequency (AAF) for CTRL samples..." >> variant_qc_calculation.log
${PLINK2_PATH} \
    --bfile ${PREFIX} \
    --keep <(awk '$6==1 {print $1, $2}' ${PREFIX}.fam) \
    --freq \
    --out ${PREFIX}.ctrl \
    --threads 8

echo "  AAF (CTRL) calculation completed" >> variant_qc_calculation.log
echo "" >> variant_qc_calculation.log

# ===== 步骤4: 计算HWE (ALL) =====
echo "Step 4: Calculating Hardy-Weinberg equilibrium for ALL samples..." >> variant_qc_calculation.log
${PLINK2_PATH} \
    --bfile ${PREFIX} \
    --hardy \
    --out ${PREFIX}.all \
    --threads 8

echo "  HWE (ALL) calculation completed" >> variant_qc_calculation.log
echo "" >> variant_qc_calculation.log

# ===== 步骤5: 计算HWE (CASE) =====
echo "Step 5: Calculating Hardy-Weinberg equilibrium for CASE samples..." >> variant_qc_calculation.log
${PLINK2_PATH} \
    --bfile ${PREFIX} \
    --keep <(awk '$6==2 {print $1, $2}' ${PREFIX}.fam) \
    --hardy \
    --out ${PREFIX}.case \
    --threads 8

echo "  HWE (CASE) calculation completed" >> variant_qc_calculation.log
echo "" >> variant_qc_calculation.log

# ===== 步骤6: 计算HWE (CTRL) =====
echo "Step 6: Calculating Hardy-Weinberg equilibrium for CTRL samples..." >> variant_qc_calculation.log
${PLINK2_PATH} \
    --bfile ${PREFIX} \
    --keep <(awk '$6==1 {print $1, $2}' ${PREFIX}.fam) \
    --hardy \
    --out ${PREFIX}.ctrl \
    --threads 8

echo "  HWE (CTRL) calculation completed" >> variant_qc_calculation.log
echo "" >> variant_qc_calculation.log

# ===== 步骤7: 合并所有统计结果 =====
echo "Step 7: Merging all statistics into summary file..." >> variant_qc_calculation.log

python3 << EOF
import pandas as pd
import sys

print("Reading PLINK2 output files...")

try:
    # 读取AAF文件 (PLINK2 --freq 输出格式)
    # 列: #CHROM ID REF ALT ALT_FREQS OBS_CT
    afreq_all = pd.read_csv("${PREFIX}.all.afreq", sep="\t")
    afreq_case = pd.read_csv("${PREFIX}.case.afreq", sep="\t")
    afreq_ctrl = pd.read_csv("${PREFIX}.ctrl.afreq", sep="\t")
    
    print(f"  AAF ALL: {len(afreq_all)} variants")
    print(f"  AAF CASE: {len(afreq_case)} variants")
    print(f"  AAF CTRL: {len(afreq_ctrl)} variants")
    
    # 读取HWE文件 (PLINK2 --hardy 输出格式)
    # 列: #CHROM ID A1 AX HETX HOMX1 HOMX2 HET_A1 HOMXAX P
    hardy_all = pd.read_csv("${PREFIX}.all.hardy", sep="\t")
    hardy_case = pd.read_csv("${PREFIX}.case.hardy", sep="\t")
    hardy_ctrl = pd.read_csv("${PREFIX}.ctrl.hardy", sep="\t")
    
    print(f"  HWE ALL: {len(hardy_all)} variants")
    print(f"  HWE CASE: {len(hardy_case)} variants")
    print(f"  HWE CTRL: {len(hardy_ctrl)} variants")
    
    # 提取需要的列并重命名
    # AAF = ALT allele frequency (即ALT_FREQS列)
    afreq_all_sub = afreq_all[["ID", "ALT_FREQS"]].rename(columns={"ALT_FREQS": "AAF_ALL"})
    afreq_case_sub = afreq_case[["ID", "ALT_FREQS"]].rename(columns={"ALT_FREQS": "AAF_CASE"})
    afreq_ctrl_sub = afreq_ctrl[["ID", "ALT_FREQS"]].rename(columns={"ALT_FREQS": "AAF_CTRL"})
    
    # HWE = P value (即P列)
    hardy_all_sub = hardy_all[["ID", "P"]].rename(columns={"P": "HWE_ALL"})
    hardy_case_sub = hardy_case[["ID", "P"]].rename(columns={"P": "HWE_CASE"})
    hardy_ctrl_sub = hardy_ctrl[["ID", "P"]].rename(columns={"P": "HWE_CTRL"})
    
    # 合并所有AAF数据
    merged = afreq_all_sub.merge(afreq_case_sub, on="ID", how="outer")
    merged = merged.merge(afreq_ctrl_sub, on="ID", how="outer")
    
    # 计算MAF (Minor Allele Frequency)
    # MAF = min(AAF, 1-AAF)
    merged["MAF_ALL"] = merged["AAF_ALL"].apply(lambda x: min(x, 1-x) if pd.notna(x) else None)
    merged["MAF_CASE"] = merged["AAF_CASE"].apply(lambda x: min(x, 1-x) if pd.notna(x) else None)
    merged["MAF_CTRL"] = merged["AAF_CTRL"].apply(lambda x: min(x, 1-x) if pd.notna(x) else None)
    
    # 合并HWE数据
    merged = merged.merge(hardy_all_sub, on="ID", how="outer")
    merged = merged.merge(hardy_case_sub, on="ID", how="outer")
    merged = merged.merge(hardy_ctrl_sub, on="ID", how="outer")
    
    # 重新排列列顺序
    merged = merged[["ID", "AAF_ALL", "AAF_CASE", "AAF_CTRL", 
                     "MAF_ALL", "MAF_CASE", "MAF_CTRL",
                     "HWE_ALL", "HWE_CASE", "HWE_CTRL"]]
    
    # 重命名ID列为VARIANT_ID
    merged = merged.rename(columns={"ID": "VARIANT_ID"})
    
    # 解析VARIANT_ID格式 (CHROM:POS:REF:ALT) 用于排序
    print("\nParsing VARIANT_ID for sorting...")
    merged[["CHROM", "POS", "REF", "ALT"]] = merged["VARIANT_ID"].str.split(":", expand=True)
    
    # 转换POS为整数
    merged["POS"] = merged["POS"].astype(int)
    
    # 提取染色体编号用于排序 (去除chr前缀)
    # chr1 -> 1, chr2 -> 2, ..., chr22 -> 22, chrX -> 23, chrY -> 24, chrM -> 25
    def get_chr_sort_key(chrom):
        chrom_clean = chrom.replace("chr", "").replace("Chr", "").replace("CHR", "")
        if chrom_clean.isdigit():
            return int(chrom_clean)
        elif chrom_clean.upper() == "X":
            return 23
        elif chrom_clean.upper() == "Y":
            return 24
        elif chrom_clean.upper() in ["M", "MT"]:
            return 25
        elif chrom_clean.upper() == "PAR1":
            return 26
        elif chrom_clean.upper() == "PAR2":
            return 27
        else:
            return 99  # 其他未知染色体排在最后
    
    merged["CHR_SORT_KEY"] = merged["CHROM"].apply(get_chr_sort_key)
    
    # 按照染色体编号(主升序)和位置(次升序)排序
    merged = merged.sort_values(by=["CHR_SORT_KEY", "POS"])
    
    # 删除临时排序列
    merged = merged.drop(columns=["CHROM", "POS", "REF", "ALT", "CHR_SORT_KEY"])
    
    print(f"  Sorted by chromosome (chr1-chr22) and position")
    
    # 保存为TSV文件
    output_file = "${PREFIX}.variant_qc_sum.tsv"
    merged.to_csv(output_file, sep="\t", index=False, na_rep="NA")
    
    print(f"  Summary file created: {output_file}")
    print(f"  Total variants in summary: {len(merged)}")
    
    # 统计一些基本信息
    print("\nStatistics summary:")
    print(f"  AAF_ALL range: [{merged['AAF_ALL'].min():.4f}, {merged['AAF_ALL'].max():.4f}]")
    print(f"  MAF_ALL range: [{merged['MAF_ALL'].min():.4f}, {merged['MAF_ALL'].max():.4f}]")
    print(f"  HWE_ALL range: [{merged['HWE_ALL'].min():.4e}, {merged['HWE_ALL'].max():.4e}]")
    
    # 统计MAF和AAF不同的变异数 (ALT allele是major allele的情况)
    # MAF != AAF 意味着 AAF > 0.5，即ALT allele是major allele
    maf_aaf_diff_all = (merged['MAF_ALL'] != merged['AAF_ALL']).sum()
    maf_aaf_diff_case = (merged['MAF_CASE'] != merged['AAF_CASE']).sum()
    maf_aaf_diff_ctrl = (merged['MAF_CTRL'] != merged['AAF_CTRL']).sum()
    
    print(f"\nVariants where MAF != AAF (ALT allele is major allele, AAF > 0.5):")
    print(f"  ALL: {maf_aaf_diff_all} variants ({maf_aaf_diff_all/len(merged)*100:.2f}%)")
    print(f"  CASE: {maf_aaf_diff_case} variants ({maf_aaf_diff_case/len(merged)*100:.2f}%)")
    print(f"  CTRL: {maf_aaf_diff_ctrl} variants ({maf_aaf_diff_ctrl/len(merged)*100:.2f}%)")
    print("  Note: MAF = min(AAF, 1-AAF), so MAF != AAF when ALT is the major allele")
    
    # 统计HWE显著偏离的变异数 (P < 1e-6)
    hwe_fail_all = (merged['HWE_ALL'] < 1e-6).sum()
    hwe_fail_case = (merged['HWE_CASE'] < 1e-6).sum()
    hwe_fail_ctrl = (merged['HWE_CTRL'] < 1e-6).sum()
    
    print(f"\nHWE violations (P < 1e-6):")
    print(f"  ALL: {hwe_fail_all} variants ({hwe_fail_all/len(merged)*100:.2f}%)")
    print(f"  CASE: {hwe_fail_case} variants ({hwe_fail_case/len(merged)*100:.2f}%)")
    print(f"  CTRL: {hwe_fail_ctrl} variants ({hwe_fail_ctrl/len(merged)*100:.2f}%)")
    
except Exception as e:
    print(f"Error processing files: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)
EOF

echo "  Statistics merged successfully" >> variant_qc_calculation.log
echo "" >> variant_qc_calculation.log

# ===== 最终统计汇总 =====
echo "======================================" >> variant_qc_calculation.log
echo "SUMMARY" >> variant_qc_calculation.log
echo "======================================" >> variant_qc_calculation.log

printf "Region: %s\n" "${REGION}" >> variant_qc_calculation.log
printf "Total variants: %'d\n" ${input_variant_count} >> variant_qc_calculation.log
printf "Samples used:\n" >> variant_qc_calculation.log
printf "  ALL: %'d\n" ${input_sample_count} >> variant_qc_calculation.log
printf "  CASE: %'d\n" ${case_count} >> variant_qc_calculation.log
printf "  CTRL: %'d\n" ${ctrl_count} >> variant_qc_calculation.log
echo "" >> variant_qc_calculation.log

echo "Statistics calculated:" >> variant_qc_calculation.log
echo "  1. AAF (ALT allele frequency): Frequency of the ALT allele (bim column 5)" >> variant_qc_calculation.log
echo "  2. MAF (Minor allele frequency): min(AAF, 1-AAF)" >> variant_qc_calculation.log
echo "  3. HWE (Hardy-Weinberg equilibrium): P-value for deviation from HWE" >> variant_qc_calculation.log
echo "" >> variant_qc_calculation.log
echo "Important notes:" >> variant_qc_calculation.log
echo "  - AAF = Frequency of ALT allele (can be >0.5 if ALT is the major allele)" >> variant_qc_calculation.log
echo "  - MAF = Frequency of the less common allele (always ≤0.5)" >> variant_qc_calculation.log
echo "  - When AAF > 0.5: ALT is the major allele, MAF = 1 - AAF" >> variant_qc_calculation.log
echo "  - When AAF ≤ 0.5: ALT is the minor allele, MAF = AAF" >> variant_qc_calculation.log
echo "  - MAF != AAF indicates ALT allele is more common than REF allele" >> variant_qc_calculation.log
echo "" >> variant_qc_calculation.log
echo "Stratification:" >> variant_qc_calculation.log
echo "  - ALL: All samples included" >> variant_qc_calculation.log
echo "  - CASE: Only case samples (phenotype = 2)" >> variant_qc_calculation.log
echo "  - CTRL: Only control samples (phenotype = 1)" >> variant_qc_calculation.log
echo "" >> variant_qc_calculation.log

echo "Output files:" >> variant_qc_calculation.log
echo "  - ${PREFIX}.variant_qc_sum.tsv (main summary file)" >> variant_qc_calculation.log
echo "      Columns: VARIANT_ID, AAF_ALL, AAF_CASE, AAF_CTRL," >> variant_qc_calculation.log
echo "               MAF_ALL, MAF_CASE, MAF_CTRL," >> variant_qc_calculation.log
echo "               HWE_ALL, HWE_CASE, HWE_CTRL" >> variant_qc_calculation.log
echo "  - ${PREFIX}.all.afreq (PLINK2 frequency output, ALL)" >> variant_qc_calculation.log
echo "  - ${PREFIX}.case.afreq (PLINK2 frequency output, CASE)" >> variant_qc_calculation.log
echo "  - ${PREFIX}.ctrl.afreq (PLINK2 frequency output, CTRL)" >> variant_qc_calculation.log
echo "  - ${PREFIX}.all.hardy (PLINK2 HWE output, ALL)" >> variant_qc_calculation.log
echo "  - ${PREFIX}.case.hardy (PLINK2 HWE output, CASE)" >> variant_qc_calculation.log
echo "  - ${PREFIX}.ctrl.hardy (PLINK2 HWE output, CTRL)" >> variant_qc_calculation.log
echo "  - variant_qc_calculation.log (this report)" >> variant_qc_calculation.log
