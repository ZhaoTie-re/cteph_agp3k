#!/bin/bash

# 脚本功能：根据bed文件前缀提取对应fam文件中的样本ID
# 输出文件保存在当前工作目录

# 设置bed文件前缀（请根据实际情况修改）
bed_prefix="/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/wgs/19.tommo_panel_filter/cteph_agp3k.lowfreq_common"

# 设置输出文件名（保存在当前工作目录）
output_file="$(basename ${bed_prefix}).samples.txt"

# 检查fam文件是否存在
fam_file="${bed_prefix}.fam"

if [ ! -f "$fam_file" ]; then
    echo "错误: 找不到fam文件: $fam_file"
    exit 1
fi

# 提取样本ID (fam文件第2列为样本ID)
echo "正在从 $fam_file 提取样本ID..."
awk '{print $2}' "$fam_file" > "$output_file"

# 统计样本数量
sample_count=$(wc -l < "$output_file")

echo "完成! 已从 $fam_file 提取 $sample_count 个样本ID"
echo "结果保存在: $output_file"

# 显示前几行预览
echo ""
echo "样本ID预览 (前10行):"
head -n 10 "$output_file"

if [ $sample_count -gt 10 ]; then
    echo "..."
    echo "(总共 $sample_count 个样本)"
fi