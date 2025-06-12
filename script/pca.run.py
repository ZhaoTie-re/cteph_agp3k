# %%
import pandas as pd
import numpy as np
import subprocess
import os
import argparse

# Set up environment
os.environ["PATH"] = "/home/b/b37974/:" + os.environ["PATH"]
# Argument parser for command line execution
parser = argparse.ArgumentParser(description="PCA and visualization for plink2 data")
parser.add_argument("--hwe_pass_ids", type=str, help="Path to the file with IDs passing HWE QC")
parser.add_argument("--check_maf", type=str, help="Path to the file with MAF check results")
parser.add_argument("--maf_selection", type=float, default=0.01, help="MAF threshold for filtering")
parser.add_argument("--bed_prefix", type=str, help="Path to the bed prefix for plink2 operations")
parser.add_argument("--high_ld", type=str, help="Path to the high LD regions file (plinkQC provided)")
parser.add_argument("--threads", type=int, default=8, help="Number of threads for plink2 operations")
parser.add_argument("--case_prefix", type=str, default="PHOM", help="Prefix for case samples (used in plots)")
parser.add_argument("--case_name", type=str, default="CTEPH", help="Name for case samples (used in plots)")
parser.add_argument("--control_name", type=str, default="AGP3K", help="Name for control samples (used in plots)")
args = parser.parse_args()

hwe_pass_ids = args.hwe_pass_ids
check_maf = args.check_maf
maf_selection = args.maf_selection

# 逐行读取check_maf文件，筛选MAF大于等于maf_selection的ID，避免一次性加载全部数据
selected_ids = set()
with open(check_maf, 'r') as f:
    header = next(f)
    for line in f:
        parts = line.rstrip('\n').split('\t')
        if len(parts) < 2:
            continue
        try:
            maf = float(parts[1])
        except ValueError:
            continue
        if maf >= maf_selection:
            selected_ids.add(parts[0])

# 读取hwe_pass_ids文件
with open(hwe_pass_ids, 'r') as f:
    hwe_ids = set(line.strip() for line in f if line.strip())

# 取交集
intersect_ids = set(selected_ids) & hwe_ids

# 保存到新文件
output_file = "pass.hwe.maf.txt"
with open(output_file, 'w') as f:
    for id_ in sorted(intersect_ids):
        f.write(f"{id_}\n")

# %%

bed_prefix = args.bed_prefix
output_prefix = f"cteph_agp3k.s_qc.gt_qc.v_qc.hwe.pass.hwe.maf{maf_selection}"

cmd = [
    "/home/b/b37974/plink2",
    "--bfile", bed_prefix,
    "--extract", "pass.hwe.maf.txt",
    "--snps-only", "just-acgt",
    "--make-bed",
    "--out", output_prefix,
    "--threads", str(args.threads)
]

subprocess.run(cmd, check=True)


# %%
high_ld = args.high_ld

cmd_ld = [
    "/home/b/b37974/plink2",
    "--bfile", output_prefix,
    "--exclude", "range", high_ld,
    "--indep-pairwise", "50", "5", "0.2",
    "--make-bed",
    "--out", f"{output_prefix}.no_high_ld",
    "--threads", str(args.threads)
]

subprocess.run(cmd_ld, check=True)

# %%
cmd_pca = [
    "/home/b/b37974/plink2",
    "--bfile", f"{output_prefix}.no_high_ld",
    "--extract", f"{output_prefix}.no_high_ld.prune.in",
    "--pca", "allele-wts", "10",
    "--out", f"{output_prefix}.no_high_ld.prune.pca",
    "--threads", str(args.threads)
]
subprocess.run(cmd_pca, check=True)

# %%
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages

# 读取数据
eigenvec_file = f"{output_prefix}.no_high_ld.prune.pca.eigenvec"
eigenval_file = f"{output_prefix}.no_high_ld.prune.pca.eigenval"
pca_df = pd.read_csv(eigenvec_file, delim_whitespace=True)
eigenval_df = pd.read_csv(eigenval_file, delim_whitespace=True, header=None)

# 分组标记
case_prefix = args.case_prefix
case_name = args.case_name
control_name = args.control_name
pca_df['Group'] = pca_df['IID'].apply(lambda x: f'{case_name}' if str(x).startswith(f'{case_prefix}') else f'{control_name}')

plt.style.use('default') 
# 设置美化主题
sns.set_theme(style="whitegrid", context="talk", font_scale=1.2)

# 创建PDF文件
with PdfPages('pca_pairwise_plots.pdf') as pdf:
    for i in range(0, 10, 2):  # PC1 vs PC2, PC3 vs PC4, ..., PC9 vs PC10
        pc_x = f"PC{i+1}"
        pc_y = f"PC{i+2}"
        pc_x_var = eigenval_df.iloc[i, 0] / eigenval_df[0].sum() * 100
        pc_y_var = eigenval_df.iloc[i+1, 0] / eigenval_df[0].sum() * 100

        # 创建图像
        plt.figure(figsize=(10, 8))
        ax = sns.scatterplot(
            data=pca_df,
            x=pc_x, y=pc_y,
            hue='Group',
            palette={f'{control_name}': '#4DBBD5', f'{case_name}': '#E64B35'},
            s=40, alpha=0.85,
            edgecolor='black', linewidth=0.4
        )

        # 设置标签和标题
        plt.xlabel(f"{pc_x} ({pc_x_var:.2f}%)", fontsize=16)
        plt.ylabel(f"{pc_y} ({pc_y_var:.2f}%)", fontsize=16)
        plt.title(f"PCA: {pc_x} vs {pc_y}", fontsize=18, weight='bold')

        # 图例设置
        plt.legend(
            title=None,
            loc='upper left',
            bbox_to_anchor=(1.02, 1),
            borderaxespad=0,
            frameon=False
        )

        # 辅助线与网格
        plt.grid(True, linestyle='--', linewidth=0.5, color='gray', alpha=0.4)
        plt.axhline(0, color='gray', linewidth=0.6, linestyle='--', alpha=0.5)
        plt.axvline(0, color='gray', linewidth=0.6, linestyle='--', alpha=0.5)
        
        plt.tight_layout()
        pdf.savefig()  # 保存当前页
        plt.close()    # 关闭当前图，释放内存
