#!/bin/zsh
#SBATCH --job-name=tuning_vmiss_calculation
#SBATCH --output=tuning_vmiss_%j.out
#SBATCH --error=tuning_vmiss_%j.err
#SBATCH -p gr10478b
#SBATCH -t 168:0:0
#SBATCH --rsc p=1:t=16:c=8:m=36568M

# 设置工作目录
cd /LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/tuning.vmiss

# 激活conda环境
source $(conda info --base)/etc/profile.d/conda.sh
conda activate cteph_geno_pro

# 定义bed文件前缀并导出为环境变量
export BED_PREFIX="/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/wgs/14.rm_maf0_vmiss/cteph_agp3k.sqc.rm_maf0_vmiss1"
# 定义info文件路径并导出为环境变量
export INFO_PATH="/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/info/cteph_agp3k_jhrpv4.rev1.xlsx"   

# 创建输出目录（如果不存在）
OUTPUT_DIR="./vmiss_results"
mkdir -p ${OUTPUT_DIR}   

# 运行Python脚本
python3 << 'EOF'
import pandas as pd
import numpy as np
import sys
import importlib

# 添加scripts目录到路径
sys.path.insert(0, '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/tuning.vmiss/scripts')

# 导入或重新加载模块
if 'tuning_vmiss_tools' in sys.modules:
    importlib.reload(sys.modules['tuning_vmiss_tools'])
else:
    import tuning_vmiss_tools

from tuning_vmiss_tools import calculate_variant_metrics

# 读取环境变量中的bed_prefix
import os
bed_prefix = os.environ.get('BED_PREFIX')
info_path = os.environ.get('INFO_PATH')

# 调用函数计算变体指标
# 返回输出文件路径，不会将大量数据加载到内存中
output_file = calculate_variant_metrics(
    bed_prefix=bed_prefix,
    info_path=info_path,           # 指定info文件路径
    output_dir='./vmiss_results',  # 指定输出目录
    maf_group='ctrl',              # 使用对照组(AGP3K)计算MAF
    pheno_column='OUTCOME2',       # 表型列（默认值，AGP3K=ctrl, CTEPH=case）
    n_threads=16                    # 使用8个线程加速
)

# 显示结果文件路径
print(f"\n✓ 计算完成！结果文件: {output_file}")

# 可选：读取并预览前几行（仅用于检查，不会加载整个文件）
print(f"\n前10行预览:")
df_preview = pd.read_csv(output_file, sep='\t', nrows=10)
print(df_preview)

# 显示统计信息
print(f"\n文件列名: {df_preview.columns.tolist()}")
print(f"\n说明：")
print(f"  - VMISS: 所有样本的变体缺失率")
print(f"  - VMISS_15X: 15X深度样本的变体缺失率")
print(f"  - VMISS_30X: 30X深度样本的变体缺失率")
print(f"  - MAF: 对照组(AGP3K)的次要等位基因频率")
EOF

echo "任务完成！"
