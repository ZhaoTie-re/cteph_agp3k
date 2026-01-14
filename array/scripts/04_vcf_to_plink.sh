#!/bin/bash
# 04_vcf_to_plink.sh
# Usage: 04_vcf_to_plink.sh <plink2_path> <sample_info> <id_col> <sex_col> <vcf_file> <prefix>

PLINK2_PATH=$1
SAMPLE_INFO=$2
ID_COL=$3
SEX_COL=$4
VCF_FILE=$5
PREFIX=$6

# 执行Python脚本读取Excel并创建性别更新文件
python3 << EOF
import pandas as pd
import sys

# 读取Excel文件中的样本信息
print("Reading sample information from Excel file...")
try:
    df = pd.read_excel("${SAMPLE_INFO}")
    print(f"Total records in Excel: {len(df)}")
    
    # 检查必需的列是否存在
    if "${ID_COL}" not in df.columns or "${SEX_COL}" not in df.columns:
        print(f"Error: Required columns '${ID_COL}' and/or '${SEX_COL}' not found in Excel file")
        print(f"Available columns: {', '.join(df.columns)}")
        sys.exit(1)
    
    # 提取ID和Sex列
    sex_data = df[["${ID_COL}", "${SEX_COL}"]].copy()
    
    # 转换性别编码: F->2 (female), M->1 (male)
    sex_mapping = {'F': '2', 'M': '1', 'f': '2', 'm': '1'}
    sex_data['Sex_Code'] = sex_data["${SEX_COL}"].map(sex_mapping)
    
    # 检查是否有未识别的性别代码
    unknown_sex = sex_data[sex_data['Sex_Code'].isna()]
    if len(unknown_sex) > 0:
        print(f"Warning: {len(unknown_sex)} samples with unknown sex codes:")
        print(unknown_sex["${SEX_COL}"].value_counts())
        print("These will be set as unknown (0)")
        sex_data['Sex_Code'] = sex_data['Sex_Code'].fillna('0')
    
    # 创建PLINK格式的性别更新文件: FID IID Sex
    # FID和IID都使用样本ID
    with open("${PREFIX}.update_sex.txt", "w") as f:
        for _, row in sex_data.iterrows():
            sample_id = str(row["${ID_COL}"])
            sex_code = row['Sex_Code']
            f.write(f"{sample_id}\t{sample_id}\t{sex_code}\n")
    
    print(f"Sex information written for {len(sex_data)} samples")
    print(f"  Females (F): {(sex_data['Sex_Code'] == '2').sum()}")
    print(f"  Males (M): {(sex_data['Sex_Code'] == '1').sum()}")
    print(f"  Unknown: {(sex_data['Sex_Code'] == '0').sum()}")
    
except Exception as e:
    print(f"Error reading Excel file: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)
EOF

# 转换VCF到PLINK格式
# 使用--double-id确保FID和IID相同
# --split-par b38: 处理X染色体伪常染色体区域（PAR）
# --update-sex: 使用从Excel读取的性别信息
${PLINK2_PATH} \
    --vcf ${VCF_FILE} \
    --double-id \
    --split-par b38 \
    --update-sex ${PREFIX}.update_sex.txt \
    --make-bed \
    --out ${PREFIX} \
    --threads 8

# 验证转换结果
echo "VCF to PLINK conversion completed" > ${PREFIX}.conversion.log
echo "Date: $(date)" >> ${PREFIX}.conversion.log
echo "" >> ${PREFIX}.conversion.log

# 统计变异数
variant_count=$(wc -l < ${PREFIX}.bim)
sample_count=$(wc -l < ${PREFIX}.fam)

printf "Variants: %'d\n" ${variant_count} >> ${PREFIX}.conversion.log
printf "Samples: %'d\n" ${sample_count} >> ${PREFIX}.conversion.log
echo "" >> ${PREFIX}.conversion.log

# 统计性别信息
echo "Sex distribution:" >> ${PREFIX}.conversion.log
male_count=$(awk '$5==1' ${PREFIX}.fam | wc -l)
female_count=$(awk '$5==2' ${PREFIX}.fam | wc -l)
unknown_count=$(awk '$5==0' ${PREFIX}.fam | wc -l)
printf "  Males: %'d\n" ${male_count} >> ${PREFIX}.conversion.log
printf "  Females: %'d\n" ${female_count} >> ${PREFIX}.conversion.log
printf "  Unknown: %'d\n" ${unknown_count} >> ${PREFIX}.conversion.log
echo "" >> ${PREFIX}.conversion.log

# 显示前几行样本ID以确认FID=IID和性别
echo "Sample ID format (first 5 samples):" >> ${PREFIX}.conversion.log
head -5 ${PREFIX}.fam | awk '{sex=$5; if(sex==1) sex_str="Male"; else if(sex==2) sex_str="Female"; else sex_str="Unknown"; print "  FID: " $1 "  IID: " $2 "  Sex: " sex_str " (" $5 ")"}' >> ${PREFIX}.conversion.log
echo "" >> ${PREFIX}.conversion.log

echo "Output files:" >> ${PREFIX}.conversion.log
echo "  - ${PREFIX}.bed (binary genotype file)" >> ${PREFIX}.conversion.log
echo "  - ${PREFIX}.bim (variant information)" >> ${PREFIX}.conversion.log
echo "  - ${PREFIX}.fam (sample information, FID=IID, with sex)" >> ${PREFIX}.conversion.log
echo "  - ${PREFIX}.update_sex.txt (sex information from Excel)" >> ${PREFIX}.conversion.log
echo "" >> ${PREFIX}.conversion.log
echo "Sex encoding: 1=Male, 2=Female, 0=Unknown" >> ${PREFIX}.conversion.log
