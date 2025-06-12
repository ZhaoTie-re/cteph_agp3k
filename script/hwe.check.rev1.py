# %%
import os
import pandas as pd
import numpy as np
import seaborn as sns
import subprocess
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import argparse

# 设置参数
parser = argparse.ArgumentParser(description="HWE check for case-control data")
parser.add_argument("--bed_prefix", type=str, required=True, help="Path to the bed prefix")
parser.add_argument("--threads", type=int, default=4, help="Number of threads to use")
parser.add_argument("--plt_maf", type=float, default=0.01, help="MAF threshold for plotting (used for plotting)")
parser.add_argument("--plt_hwe_control", type=float, default=1e-6, help="HWE threshold for control group (used for plotting)")
parser.add_argument("--plt_hwe_case", type=float, default=1e-10, help="HWE threshold for case group (used for plotting)")
parser.add_argument("--mode", type=str, choices=['only_control', 'case_control'], default='case_control', help="Mode for HWE filtering: 'only_control' or 'case_control'")
parser.add_argument("--control_hwe_threshold", type=float, default=1e-6, help="Control group HWE threshold for filtering")
parser.add_argument("--case_hwe_threshold", type=float, default=1e-10, help="Case group HWE threshold for filtering")
parser.add_argument("--out", type=str, help="Output prefix for the hwe filtered plink file")
args = parser.parse_args()

os.environ["PATH"] = "/home/b/b37974/:" + os.environ["PATH"]

bed_prefix = args.bed_prefix
hwe_qc_prefix = "hwe_check"

maf_thresfold = args.plt_maf

aaf_output = f"{hwe_qc_prefix}.aaf"
subprocess.run([
    "/home/b/b37974/plink2",
    "--bfile", bed_prefix,
    "--freq",
    "--keep-if", "PHENO1==1", 
    "--out", aaf_output, 
    "--threads", str(args.threads)
])

case_hardy_output = f"{hwe_qc_prefix}.case"
subprocess.run([
    "/home/b/b37974/plink2",
    "--bfile", bed_prefix,
    "--hardy",
    "--keep-if", "PHENO1==2", 
    "--out", case_hardy_output,
    "--threads", str(args.threads)
])

control_hardy_output = f"{hwe_qc_prefix}.control"
subprocess.run([
    "/home/b/b37974/plink2",
    "--bfile", bed_prefix,
    "--hardy",
    "--keep-if", "PHENO1==1", 
    "--out", control_hardy_output,
    "--threads", str(args.threads)
])


# %%
# 参数设置
columns_afreq = ['ID', 'ALT_FREQS']
columns_hardy = ['ID', 'P']
chunksize = 100000  # 可适当加大以减少IO

# 读取 afreq 文件并计算 MAF
afreq_chunks = []
for chunk in pd.read_csv(f"{aaf_output}.afreq", delim_whitespace=True,
                         usecols=columns_afreq, chunksize=chunksize,
                         dtype={'ID': str, 'ALT_FREQS': float}):
    chunk['MAF'] = chunk['ALT_FREQS'].clip(upper=1).apply(lambda x: min(x, 1 - x))
    afreq_chunks.append(chunk[['ID', 'MAF']])
afreq_df = pd.concat(afreq_chunks, ignore_index=True)
afreq_df.to_csv(f"{hwe_qc_prefix}.maf", index=False, sep='\t')

# 读取 case hardy
case_chunks = []
for chunk in pd.read_csv(f"{case_hardy_output}.hardy", delim_whitespace=True,
                         usecols=columns_hardy, chunksize=chunksize,
                         dtype={'ID': str, 'P': float}):
    chunk = chunk.rename(columns={'P': 'CASE_HWE'})
    case_chunks.append(chunk)
case_df = pd.concat(case_chunks, ignore_index=True)

# 读取 control hardy
control_chunks = []
for chunk in pd.read_csv(f"{control_hardy_output}.hardy", delim_whitespace=True,
                         usecols=columns_hardy, chunksize=chunksize,
                         dtype={'ID': str, 'P': float}):
    chunk = chunk.rename(columns={'P': 'CONTROL_HWE'})
    control_chunks.append(chunk)
control_df = pd.concat(control_chunks, ignore_index=True)

# 合并三个表（按ID）
merged_df = afreq_df.merge(case_df, on='ID', how='left').merge(control_df, on='ID', how='left')

# %%
# 阈值设置(用于绘图)
threshold_control = args.plt_hwe_control
threshold_case = args.plt_hwe_case

def plot_hwe_comparison(plot_data, maf_label, ax_main, ax_top, ax_right):
    x = plot_data['CONTROL_HWE']
    y = plot_data['CASE_HWE']

    # 四象限计数
    q1 = ((x >= threshold_control) & (y >= threshold_case)).sum()
    q2 = ((x < threshold_control) & (y >= threshold_case)).sum()
    q3 = ((x < threshold_control) & (y < threshold_case)).sum()
    q4 = ((x >= threshold_control) & (y < threshold_case)).sum()

    # 主图
    ax_main.scatter(x, y, alpha=0.3, c='grey', s=20)
    ax_main.axvline(x=threshold_control, color='red', linestyle='--', label=f'Control HWE threshold = {threshold_control}')
    ax_main.axhline(y=threshold_case, color='blue', linestyle='--', label=f'Case HWE threshold = {threshold_case}')
    ax_main.set_xscale('log')
    ax_main.set_yscale('log')
    ax_main.set_xlabel('CONTROL_HWE_P')
    ax_main.set_ylabel('CASE_HWE_P')
    ax_main.grid(True)

    # KDE 顶部
    sns.kdeplot(x, ax=ax_top, fill=True, color='gray', linewidth=1.5)
    ax_top.set_xscale('log')
    ax_top.axvline(x=threshold_control, color='red', linestyle='--')
    ax_top.set_xlabel('')
    ax_top.set_ylabel('')
    ax_top.tick_params(axis='both', which='both', bottom=False, top=False, left=False, right=False, labelbottom=False, labelleft=False)

    # KDE 右侧
    sns.kdeplot(y, ax=ax_right, fill=True, color='gray', linewidth=1.5, vertical=True)
    ax_right.set_yscale('log')
    ax_right.axhline(y=threshold_case, color='blue', linestyle='--')
    ax_right.set_xlabel('')
    ax_right.set_ylabel('')
    ax_right.tick_params(axis='both', which='both', bottom=False, top=False, left=False, right=False, labelbottom=False, labelleft=False)

    # 插图
    top_ylim = ax_top.get_ylim()
    right_xlim = ax_right.get_xlim()
    inset_ax = inset_axes(ax_main, width="25%", height="25%", loc='upper right',
                          bbox_to_anchor=(0.28, 0.28, 1, 1), bbox_transform=ax_main.transAxes, borderpad=0)
    inset_ax.set_xlim(right_xlim)
    inset_ax.set_ylim(top_ylim)
    for spine in inset_ax.spines.values():
        spine.set_visible(True)
    inset_ax.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)
    xmid = (right_xlim[0] + right_xlim[1]) / 2
    ymid = (top_ylim[0] + top_ylim[1]) / 2
    inset_ax.axvline(x=xmid, linestyle='--', color='red', linewidth=1.5)
    inset_ax.axhline(y=ymid, linestyle='--', color='blue', linewidth=1.5)
    x_min, x_max = inset_ax.get_xlim()
    y_min, y_max = inset_ax.get_ylim()
    inset_ax.text((xmid + x_max) / 2, (ymid + y_max) / 2, f'{q1:,}', ha='center', va='center', fontsize=6.5)
    inset_ax.text((x_min + xmid) / 2, (ymid + y_max) / 2, f'{q2:,}', ha='center', va='center', fontsize=6.5)
    inset_ax.text((x_min + xmid) / 2, (y_min + ymid) / 2, f'{q3:,}', ha='center', va='center', fontsize=6.5)
    inset_ax.text((xmid + x_max) / 2, (y_min + ymid) / 2, f'{q4:,}', ha='center', va='center', fontsize=6.5)
    ax_main.text(0.2, 0.2, maf_label, transform=ax_main.transAxes,
                 fontsize=12, fontweight='bold', va='top', ha='left',
                 bbox=dict(facecolor='white', edgecolor='none', alpha=0.7))

# 图形布局：两列子图
fig = plt.figure(figsize=(16, 8))
outer_gs = gridspec.GridSpec(1, 2, width_ratios=[1, 1], wspace=0.2)

legend_handles = []
legend_labels = []

for i, (maf_filter, label) in enumerate([(merged_df['MAF'] >= maf_thresfold, f'MAF ≥ {maf_thresfold}'),
                                         (merged_df['MAF'] < maf_thresfold, f'MAF < {maf_thresfold}')]):
    plot_data = merged_df[maf_filter]
    inner_gs = gridspec.GridSpecFromSubplotSpec(2, 2, subplot_spec=outer_gs[i],
                                                width_ratios=[4, 1], height_ratios=[1, 4],
                                                wspace=0.05, hspace=0.05)
    ax_main = plt.subplot(inner_gs[1, 0])
    ax_top = plt.subplot(inner_gs[0, 0], sharex=ax_main)
    ax_right = plt.subplot(inner_gs[1, 1], sharey=ax_main)

    # 绘图
    plot_hwe_comparison(plot_data, label, ax_main, ax_top, ax_right)

    # 只收集一次 legend（第一张图即可）
    if i == 0:
        legend_handles, legend_labels = ax_main.get_legend_handles_labels()

# 添加总标题和共享图例
fig.suptitle('HWE P by MAF Threshold', fontsize=16, y=0.99)
fig.legend(legend_handles, legend_labels, loc='upper center', bbox_to_anchor=(0.5, 0.95),
           ncol=2, frameon=True, fontsize=10)
plt.tight_layout(rect=[0, 0, 1, 0.95])  # 为 legend 留出空间
# plt.show()
# save to pdf and png
fig.savefig(f"{hwe_qc_prefix}.maf_hwe.pdf", bbox_inches='tight', dpi=300)
fig.savefig(f"{hwe_qc_prefix}.maf_hwe.png", bbox_inches='tight', dpi=300)

# %%
def get_pass_hwe_ids(merged_df, mode, control_hwe_threshold=None, case_hwe_threshold=None):
    """
    根据mode和阈值筛选通过HWE的ID。

    参数:
        merged_df: 包含ID, CONTROL_HWE, CASE_HWE等列的DataFrame
        mode: 'only_control' 或 'case_control'
        control_hwe_threshold: 控制组HWE阈值 (float)
        case_hwe_threshold: 病例组HWE阈值 (float, 仅case_control模式需要)

    返回:
        pass_hwe: 满足条件的ID组成的pd.Series
    """
    if mode == 'only_control':
        if control_hwe_threshold is None:
            raise ValueError("control_hwe_threshold must be provided for only_control mode.")
        pass_hwe = merged_df.loc[merged_df['CONTROL_HWE'] >= control_hwe_threshold, 'ID']
    elif mode == 'case_control':
        if control_hwe_threshold is None or case_hwe_threshold is None:
            raise ValueError("Both control_hwe_threshold and case_hwe_threshold must be provided for case_control mode.")
        pass_hwe = merged_df.loc[
            (merged_df['CONTROL_HWE'] >= control_hwe_threshold) &
            (merged_df['CASE_HWE'] >= case_hwe_threshold),
            'ID'
        ]
    else:
        raise ValueError("mode must be 'only_control' or 'case_control'.")
    return pass_hwe

# 下面的参数用于筛选通过HWE测试的ID
mode = args.mode

if mode == 'only_control':
    pass_ids = get_pass_hwe_ids(merged_df, mode, control_hwe_threshold=args.control_hwe_threshold)
elif mode == 'case_control':
    pass_ids = get_pass_hwe_ids(merged_df, mode, control_hwe_threshold=args.control_hwe_threshold, case_hwe_threshold=args.case_hwe_threshold)
else:
    raise ValueError("Unsupported mode.")

print(f"The number of IDs passing HWE test in {mode} mode: {len(pass_ids)}")
pass_ids.to_csv("pass_hwe_ids.txt", index=False, header=False)

# %%
# Define the output prefix for the subset
subset_bed_prefix = args.out

# Run plink2 to extract the subset
subprocess.run([
    "/home/b/b37974/plink2",
    "--bfile", bed_prefix,
    "--extract", "pass_hwe_ids.txt",
    "--make-bed",
    "--out", subset_bed_prefix, 
    "--threads", str(args.threads)
])


