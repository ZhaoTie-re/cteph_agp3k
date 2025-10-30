#!/bin/bash
#SBATCH --job-name=cteph_plot
#SBATCH --output=cteph_plot_%j.out
#SBATCH --error=cteph_plot_%j.err
#SBATCH --rsc p=1:t=16:c=8:m=36568M
#SBATCH --time=01:00:00

cd /LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/snippet_analysis/02.array_check_impute

echo "绘图作业开始: $(date)"
echo "作业ID: $SLURM_JOB_ID"

# 激活conda环境
source ~/.bashrc
conda activate cteph_geno_pro

# 创建绘图脚本
cat > plot_script.py << 'EOF'
import sys, importlib
sys.path.insert(0, "/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/snippet_analysis/scripts")

import array_check_prepare_tools
importlib.reload(array_check_prepare_tools)

import matplotlib
matplotlib.use('Agg')

from array_check_prepare_tools import plot_wgs_array_beta_comparison

# 使用和notebook中相同的参数
plot = plot_wgs_array_beta_comparison(
    merged_sumstat_file='/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/snippet_analysis/02.array_check_impute/cteph_agp3k.array_wgs_compration.merged_sumstat.tsv',
    output_prefix="cteph_agp3k.array_wgs_compration",
    save_pdf=False,
)

print("绘图完成！")
print(plot.get('summary_text', ''))
EOF

python plot_script.py

echo "绘图作业完成: $(date)"
ls -lah cteph_agp3k.array_wgs_compration.*
