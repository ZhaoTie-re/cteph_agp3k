# %%
import pandas as pd
import numpy as np
import subprocess
import argparse
import os

os.environ["PATH"] = "/home/b/b37974/:" + os.environ["PATH"]

parser = argparse.ArgumentParser(description="Sample QC for heterozygosity and missingness")
parser.add_argument("--bed_prefix", type=str, help="Path to the bed prefix")
parser.add_argument("--threads", type=int, default=8, help="Number of threads")
parser.add_argument("--miss_threshold", type=float, default=0.01, help="Missingness threshold")
parser.add_argument("--het_threshold", type=str, default='5sd', help="Heterozygosity threshold (e.g., '5sd')")
args = parser.parse_args()


bed_prefix = args.bed_prefix
sample_qc_prefix = "sample.qc"


sample_miss_cmd = [
    "/home/b/b37974/plink2",
    "--bfile", bed_prefix,
    "--missing",
    "--out", f"{sample_qc_prefix}.miss",
    "--threads", str(args.threads)
]

sample_het_cmd = [
    "/home/b/b37974/plink2",
    "--bfile", bed_prefix,
    "--het",
    "--out", f"{sample_qc_prefix}.het",
    "--threads", str(args.threads)
]

subprocess.run(sample_miss_cmd, check=True)
subprocess.run(sample_het_cmd, check=True)

# %%
heterozygosity_df = pd.read_csv(f"{sample_qc_prefix}.het.het", delim_whitespace=True)[["#FID", "IID", "F"]]
missing_df = pd.read_csv(f"{sample_qc_prefix}.miss.smiss", delim_whitespace=True)[["#FID", "IID", "F_MISS"]]

merged_df = pd.merge(heterozygosity_df, missing_df, on=["#FID", "IID"])

# %%
s_miss_threshold = args.miss_threshold
s_het_threshold = args.het_threshold
threshold_value = int(s_het_threshold[:-2])

# %%
import matplotlib.pyplot as plt
import seaborn as sns

plt.style.use('default')

fig, ax = plt.subplots(2, 1, figsize=(13, 5), 
                       gridspec_kw={'height_ratios': [1, 4]}, 
                       sharex=True)

# KDE plot
sns.kdeplot(merged_df['F_MISS'], ax=ax[0], color='navy')

# 均值和标准差
mean_f = merged_df['F'].mean()
std_f = merged_df['F'].std()
f_upper = mean_f + threshold_value * std_f
f_lower = mean_f - threshold_value * std_f

# 分组条件
cond_f_outlier = (merged_df['F'] > f_upper) | (merged_df['F'] < f_lower)
cond_miss_outlier = merged_df['F_MISS'] > s_miss_threshold
only_f_outlier = cond_f_outlier & (~cond_miss_outlier)
only_miss_outlier = (~cond_f_outlier) & cond_miss_outlier
both_outlier = cond_f_outlier & cond_miss_outlier
normal = ~(cond_f_outlier | cond_miss_outlier)

# 分别绘制不同颜色的点
ax[1].scatter(merged_df.loc[normal, 'F_MISS'], merged_df.loc[normal, 'F'], 
              c='steelblue', s=20, label=f'Pass (n={normal.sum()})', alpha=0.9)
ax[1].scatter(merged_df.loc[only_f_outlier, 'F_MISS'], merged_df.loc[only_f_outlier, 'F'], 
              c='green', s=20, label=f'Only F Outlier (n={only_f_outlier.sum()})', alpha=0.9)
ax[1].scatter(merged_df.loc[only_miss_outlier, 'F_MISS'], merged_df.loc[only_miss_outlier, 'F'], 
              c='darkorange', s=20, label=f'Only Missingness Outlier (n={only_miss_outlier.sum()})', alpha=0.9)
ax[1].scatter(merged_df.loc[both_outlier, 'F_MISS'], merged_df.loc[both_outlier, 'F'], 
              c='firebrick', s=20, label=f'Both Outliers (n={both_outlier.sum()})', alpha=0.9)

# 阈值线：水平（F）和垂直（Missing）
ax[1].axvline(s_miss_threshold, color='darkred', linestyle='--', lw=1)
ax[1].axhline(f_upper, color='darkred', linestyle='--', lw=1)
ax[1].axhline(f_lower, color='darkred', linestyle='--', lw=1)

# 坐标轴
ax[1].set_xscale('log')
ax[1].set_xlabel('Missing Rate (log scale)', fontsize=11)
ax[1].set_ylabel('Heterozygosity Rate (F)', fontsize=11)
ax[1].grid(True, linestyle='--', linewidth=0.5)

# 美化上图
ax[0].set_ylabel('')
ax[0].set_yticks([])
ax[0].spines['right'].set_visible(False)
ax[0].spines['top'].set_visible(False)

# 图例在外部右上
ax[1].legend(loc='upper left', bbox_to_anchor=(1.05, 1), borderaxespad=0, fontsize=10, frameon=False)

# 总标题
fig.suptitle(
    f'Heterozygosity vs. Missing Rate across Samples\n'
    f'(Missing threshold: {s_miss_threshold}, F threshold: ±{threshold_value}SD)',
    fontsize=12,
    x=0.32
)

plt.tight_layout(rect=[0, 0, 0.85, 0.95])  # 留出 legend 空间
# plt.show()
plt.savefig('sample_qc.pdf', dpi=300, bbox_inches='tight')

# %%
miss_sample = merged_df[
    (merged_df['F_MISS'] > s_miss_threshold)
][['#FID', 'IID']]

het_sample = merged_df[
    (merged_df['F'] > mean_f + threshold_value * std_f) | 
    (merged_df['F'] < mean_f - threshold_value * std_f)
][['#FID', 'IID']]

sample_out = merged_df[
    (merged_df['F_MISS'] > s_miss_threshold) | 
    (merged_df['F'] > mean_f + threshold_value * std_f) | 
    (merged_df['F'] < mean_f - threshold_value * std_f)
][['#FID', 'IID']]

sample_pass_qc = merged_df[
    (merged_df['F_MISS'] <= s_miss_threshold) & 
    (merged_df['F'] <= mean_f + threshold_value * std_f) & 
    (merged_df['F'] >= mean_f - threshold_value * std_f)
]['IID'].tolist()

pd.Series(sample_pass_qc).to_csv(
    "miss.het.pass.sample.txt", 
    index=False, header=False
)
