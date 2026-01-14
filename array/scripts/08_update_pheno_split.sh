#!/bin/bash
# 08_update_pheno_split.sh
# Usage: 08_update_pheno_split.sh <plink2_path> <sample_info> <id_col> <outcome_col> <case_val> <prefix>

PLINK2_PATH=$1
SAMPLE_INFO=$2
ID_COL=$3
OUTCOME_COL=$4
CASE_VAL=$5
PREFIX=$6

echo "=== Phenotype Update and Chromosome Split Report ===" > pheno_update.log
echo "Date: $(date)" >> pheno_update.log
echo "" >> pheno_update.log

# ===== 步骤1: 读取Excel并创建表型更新文件 =====
echo "Step 1: Reading phenotype information from Excel..." >> pheno_update.log

python3 << EOF
import pandas as pd
import sys

print("Reading Excel file for phenotype information...")
try:
    df = pd.read_excel("${SAMPLE_INFO}")
    print(f"Total records in Excel: {len(df)}")
    
    # 检查必需的列是否存在
    if "${ID_COL}" not in df.columns:
        print(f"Error: Required column '${ID_COL}' not found in Excel file")
        print(f"Available columns: {', '.join(df.columns)}")
        sys.exit(1)
    
    if "${OUTCOME_COL}" not in df.columns:
        print(f"Error: Required column '${OUTCOME_COL}' not found in Excel file")
        print(f"Available columns: {', '.join(df.columns)}")
        sys.exit(1)
    
    # 提取ID和分组列
    pheno_data = df[["${ID_COL}", "${OUTCOME_COL}"]].copy()
    
    # 转换表型编码: case->2, control->1, missing/other->-9
    pheno_data['Pheno_Code'] = pheno_data["${OUTCOME_COL}"].apply(
        lambda x: '2' if str(x).strip().upper() == "${CASE_VAL}".upper() else 
                  ('1' if pd.notna(x) and str(x).strip() != '' else '-9')
    )
    
    # 统计表型分布
    case_count = (pheno_data['Pheno_Code'] == '2').sum()
    ctrl_count = (pheno_data['Pheno_Code'] == '1').sum()
    missing_count = (pheno_data['Pheno_Code'] == '-9').sum()
    
    print(f"Phenotype distribution:")
    print(f"  Cases (${CASE_VAL}): {case_count}")
    print(f"  Controls (other non-missing): {ctrl_count}")
    print(f"  Missing/Unknown: {missing_count}")
    
    # 创建PLINK格式的表型更新文件: FID IID Phenotype
    with open("${PREFIX}.update_pheno.txt", "w") as f:
        for _, row in pheno_data.iterrows():
            sample_id = str(row["${ID_COL}"])
            pheno_code = row['Pheno_Code']
            f.write(f"{sample_id}\t{sample_id}\t{pheno_code}\n")
    
    print(f"Phenotype file created: ${PREFIX}.update_pheno.txt")
    
except Exception as e:
    print(f"Error reading Excel file: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)
EOF

# 记录Python脚本输出到log
echo "  Python output logged above" >> pheno_update.log
echo "" >> pheno_update.log

# ===== 步骤2: 更新PLINK fam文件的表型信息 =====
echo "Step 2: Updating phenotype information in PLINK files..." >> pheno_update.log

${PLINK2_PATH} \
    --bfile ${PREFIX} \
    --pheno ${PREFIX}.update_pheno.txt \
    --make-bed \
    --out ${PREFIX}.pheno \
    --threads 8

# 验证表型更新
echo "  Phenotype update completed" >> pheno_update.log
case_in_fam=$(awk '$6==2' ${PREFIX}.pheno.fam | wc -l)
ctrl_in_fam=$(awk '$6==1' ${PREFIX}.pheno.fam | wc -l)
missing_in_fam=$(awk '$6==-9 || $6==0' ${PREFIX}.pheno.fam | wc -l)

printf "  Cases in updated FAM: %'d\n" ${case_in_fam} >> pheno_update.log
printf "  Controls in updated FAM: %'d\n" ${ctrl_in_fam} >> pheno_update.log
printf "  Missing in updated FAM: %'d\n" ${missing_in_fam} >> pheno_update.log
echo "" >> pheno_update.log

# ===== 步骤3: 分析染色体分布 =====
echo "Step 3: Analyzing chromosome distribution..." >> pheno_update.log

# 统计每条染色体的变异数，显示时添加"Chr"前缀
echo "Chromosome distribution in input file:" >> pheno_update.log
cut -f1 ${PREFIX}.pheno.bim | sort | uniq -c | \
    awk '{printf "  Chr %s: %'"'"'d variants\n", $2, $1}' >> pheno_update.log
echo "" >> pheno_update.log

# ===== 步骤4: 拆分染色体区域 =====
echo "Step 4: Splitting genotype data by chromosomal regions..." >> pheno_update.log
echo "" > chr_split_summary.log

# 数据格式: 染色体标记为 1-22, X, Y, MT, PAR1, PAR2 (无chr前缀)
echo "  Chromosome format: numeric without 'chr' prefix (1,2,...,22,X,Y,MT,PAR1,PAR2)" >> pheno_update.log
echo "" >> pheno_update.log

# 检查是否存在常染色体 (1-22)
autosome_count=$(awk '$1 ~ /^([1-9]|1[0-9]|2[0-2])$/' ${PREFIX}.pheno.bim | wc -l)
if [ ${autosome_count} -gt 0 ]; then
    echo "  Extracting autosomes (1-22)..." >> pheno_update.log
    ${PLINK2_PATH} \
        --bfile ${PREFIX}.pheno \
        --chr 1-22 \
        --make-bed \
        --out ${PREFIX}.autosomes \
        --threads 8 2>&1 | tee -a chr_split_summary.log
    
    if [ -f ${PREFIX}.autosomes.bim ]; then
        variants=$(wc -l < ${PREFIX}.autosomes.bim)
        printf "    Autosomes: %'d variants extracted\n" ${variants} >> pheno_update.log
    fi
else
    echo "    Autosomes: No variants found (skipped)" >> pheno_update.log
fi

# 检查是否存在X染色体 (标记为 "X"，不含PAR区域)
# 注意: 如果数据中PAR1和PAR2是独立标记的，X染色体应该已经排除了PAR区域
chrx_count=$(awk '$1 == "X"' ${PREFIX}.pheno.bim | wc -l)
if [ ${chrx_count} -gt 0 ]; then
    echo "  Extracting chrX (non-PAR, marked as 'X')..." >> pheno_update.log
    ${PLINK2_PATH} \
        --bfile ${PREFIX}.pheno \
        --chr X \
        --make-bed \
        --out ${PREFIX}.chrX \
        --threads 8 2>&1 | tee -a chr_split_summary.log
    
    if [ -f ${PREFIX}.chrX.bim ]; then
        variants=$(wc -l < ${PREFIX}.chrX.bim)
        printf "    chrX (non-PAR): %'d variants extracted\n" ${variants} >> pheno_update.log
    fi
else
    echo "    chrX: No variants found (skipped)" >> pheno_update.log
fi

# 检查是否存在Y染色体 (标记为 "Y"，不含PAR区域)
# 注意: 如果数据中PAR1和PAR2是独立标记的，Y染色体应该已经排除了PAR区域
chry_count=$(awk '$1 == "Y"' ${PREFIX}.pheno.bim | wc -l)
if [ ${chry_count} -gt 0 ]; then
    echo "  Extracting chrY (non-PAR, marked as 'Y')..." >> pheno_update.log
    ${PLINK2_PATH} \
        --bfile ${PREFIX}.pheno \
        --chr Y \
        --make-bed \
        --out ${PREFIX}.chrY \
        --threads  8 2>&1 | tee -a chr_split_summary.log
    
    if [ -f ${PREFIX}.chrY.bim ]; then
        variants=$(wc -l < ${PREFIX}.chrY.bim)
        printf "    chrY (non-PAR): %'d variants extracted\n" ${variants} >> pheno_update.log
    fi
else
    echo "    chrY: No variants found (skipped)" >> pheno_update.log
fi

# 检查是否存在PAR1区域 (直接标记为 "PAR1")
par1_count=$(awk '$1 == "PAR1"' ${PREFIX}.pheno.bim | wc -l)
if [ ${par1_count} -gt 0 ]; then
    echo "  Extracting PAR1 region (marked as 'PAR1')..." >> pheno_update.log
    ${PLINK2_PATH} \
        --bfile ${PREFIX}.pheno \
        --chr PAR1 \
        --make-bed \
        --out ${PREFIX}.PAR1 \
        --threads 8 2>&1 | tee -a chr_split_summary.log
    
    if [ -f ${PREFIX}.PAR1.bim ]; then
        variants=$(wc -l < ${PREFIX}.PAR1.bim)
        printf "    PAR1: %'d variants extracted\n" ${variants} >> pheno_update.log
    fi
else
    echo "    PAR1: No variants found (skipped)" >> pheno_update.log
fi

# 检查是否存在PAR2区域 (直接标记为 "PAR2")
par2_count=$(awk '$1 == "PAR2"' ${PREFIX}.pheno.bim | wc -l)
if [ ${par2_count} -gt 0 ]; then
    echo "  Extracting PAR2 region (marked as 'PAR2')..." >> pheno_update.log
    ${PLINK2_PATH} \
        --bfile ${PREFIX}.pheno \
        --chr PAR2 \
        --make-bed \
        --out ${PREFIX}.PAR2 \
        --threads 8 2>&1 | tee -a chr_split_summary.log
    
    if [ -f ${PREFIX}.PAR2.bim ]; then
        variants=$(wc -l < ${PREFIX}.PAR2.bim)
        printf "    PAR2: %'d variants extracted\n" ${variants} >> pheno_update.log
    fi
else
    echo "    PAR2: No variants found (skipped)" >> pheno_update.log
fi

# 检查是否存在线粒体染色体 (标记为 "MT")
chrm_count=$(awk '$1 == "MT"' ${PREFIX}.pheno.bim | wc -l)
if [ ${chrm_count} -gt 0 ]; then
    echo "  Extracting chrM (mitochondrial)..." >> pheno_update.log
    ${PLINK2_PATH} \
        --bfile ${PREFIX}.pheno \
        --chr chrM \
        --make-bed \
        --out ${PREFIX}.chrM \
        --threads 8 2>&1 | tee -a chr_split_summary.log
    
    if [ -f ${PREFIX}.chrM.bim ]; then
        variants=$(wc -l < ${PREFIX}.chrM.bim)
        printf "    chrM: %'d variants extracted\n" ${variants} >> pheno_update.log
    fi
else
    echo "    chrM: No variants found (skipped)" >> pheno_update.log
fi

echo "" >> pheno_update.log

# ===== 最终统计汇总 =====
echo "======================================" >> pheno_update.log
echo "SUMMARY" >> pheno_update.log
echo "======================================" >> pheno_update.log

input_variants=$(wc -l < ${PREFIX}.pheno.bim)
input_samples=$(wc -l < ${PREFIX}.pheno.fam)

printf "Input data:\n" >> pheno_update.log
printf "  Total variants: %'d\n" ${input_variants} >> pheno_update.log
printf "  Total samples: %'d\n" ${input_samples} >> pheno_update.log
printf "  Cases: %'d\n" ${case_in_fam} >> pheno_update.log
printf "  Controls: %'d\n" ${ctrl_in_fam} >> pheno_update.log
printf "  Missing: %'d\n" ${missing_in_fam} >> pheno_update.log
echo "" >> pheno_update.log

printf "Chromosomal regions created:\n" >> pheno_update.log
for region in autosomes chrX chrY PAR1 PAR2 chrM; do
    if [ -f ${PREFIX}.${region}.bim ]; then
        variants=$(wc -l < ${PREFIX}.${region}.bim)
        printf "  %s: %'d variants\n" ${region} ${variants} >> pheno_update.log
    else
        printf "  %s: Not created (no variants)\n" ${region} >> pheno_update.log
    fi
done
echo "" >> pheno_update.log

echo "Output files:" >> pheno_update.log
echo "  - ${PREFIX}.update_pheno.txt (phenotype update file)" >> pheno_update.log
echo "  - ${PREFIX}.pheno.{bed,bim,fam} (updated phenotype files)" >> pheno_update.log
echo "  - ${PREFIX}.<region>.{bed,bim,fam} (split chromosome files)" >> pheno_update.log
echo "  - pheno_update.log (this report)" >> pheno_update.log
echo "  - chr_split_summary.log (PLINK2 split logs)" >> pheno_update.log
echo "" >> pheno_update.log

echo "Note: Only regions with variants present in the input data are created." >> pheno_update.log
echo "Missing regions indicate no variants were found for that chromosomal location." >> pheno_update.log

# 清理临时文件
rm -f ${PREFIX}.pheno.bed ${PREFIX}.pheno.bim ${PREFIX}.pheno.fam
