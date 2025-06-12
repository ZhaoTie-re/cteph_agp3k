# %%
import pandas as pd
import subprocess
import os
import argparse

# Set up environment
os.environ["PATH"] = "/home/b/b37974/:" + os.environ["PATH"]

# Argument parser for command line execution
parser = argparse.ArgumentParser(description="Perform PCA and visualization for BBJ data")
parser.add_argument("--bbj_bed_prefix", type=str, help="Path to the BBJ bed prefix for plink2 operations (after bbj_prepare.py)")
parser.add_argument("--output_prefix", type=str, default="bbj.b38.auto.sqc.vqc.norm", help="Output prefix for the processed data")
parser.add_argument("--threads", type=int, default=4, help="Number of threads for plink2 operations")
parser.add_argument("--high_ld", type=str, default="/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/info/high-LD-regions-hg38-GRCh38.txt", help="Path to the high LD regions file (plinkQC provided)")
args = parser.parse_args()

def run_command(command):
    """Run a shell command (tuple/list or string) and return the output."""
    if isinstance(command, (tuple, list)):
        result = subprocess.run(command, capture_output=True, text=True)
    else:
        result = subprocess.run(command, shell=True, capture_output=True, text=True)

    if result.returncode != 0:
        raise RuntimeError(f"Command failed:\n{command}\n{result.stderr}")
    return result.stdout.strip()

bbj_bed_prefix = args.bbj_bed_prefix
output_prefix = args.output_prefix

fam_file = f"{bbj_bed_prefix}.fam"
output_file = "update_fam.txt"

fam_df = pd.read_csv(fam_file, delim_whitespace=True, header=None, dtype=str)

update_df = pd.DataFrame({
    "#OLD_FID": fam_df[0].astype(str),
    "OLD_IID": fam_df[1].astype(str),
    "NEW_FID": "bbj_" + fam_df[0].astype(str),
    "NEW_IID": "bbj_" + fam_df[1].astype(str),
})

update_df.to_csv(output_file, sep="\t", header=True, index=False)

command_update_fam = [
    "/home/b/b37974/plink2",
    "--bfile", bbj_bed_prefix,
    "--update-ids", output_file,
    "--make-bed",
    "--out", output_prefix,
    "--threads", str(args.threads)  # Convert threads to string for command
]

run_command(command_update_fam)

# %%

high_ld = args.high_ld

command_ld = [
    "/home/b/b37974/plink2",
    "--bfile", output_prefix,
    "--snps-only", "just-acgt",
    "--exclude", "range", high_ld,
    "--indep-pairwise", "50", "5", "0.2",
    "--make-bed",
    "--out", f"{output_prefix}.no_high_ld",
    "--threads", str(args.threads)  # Convert threads to string for command
]

run_command(command_ld)

# %%
command_pca = [
    "/home/b/b37974/plink2",
    "--bfile", f"{output_prefix}.no_high_ld",
    "--extract", f"{output_prefix}.no_high_ld.prune.in",
    "--freq", "counts", 
    "--pca", "allele-wts", "approx", 
    "--out", f"{output_prefix}.no_high_ld.prune.pca",
    "--threads", str(args.threads)  # Convert threads to string for command
]

subprocess.run(command_pca)

# %%
# import matplotlib.pyplot as plt
# from matplotlib.lines import Line2D
# import seaborn as sns
# from matplotlib.backends.backend_pdf import PdfPages

# # 读取数据
# eigenvec_file = f"{output_prefix}.no_high_ld.prune.pca.eigenvec"
# eigenval_file = f"{output_prefix}.no_high_ld.prune.pca.eigenval"
# pca_df = pd.read_csv(eigenvec_file, delim_whitespace=True)
# eigenval_df = pd.read_csv(eigenval_file, delim_whitespace=True, header=None)

# plt.style.use('default') 
# # 设置美化主题
# sns.set_theme(style="whitegrid", context="talk", font_scale=1.2)

# with PdfPages("bbj_pca_pairwise_plots.pdf") as pdf:
#     for i in range(0, 10, 2):  # PC1 vs PC2, ..., PC9 vs PC10
#         pc_x = f"PC{i+1}"
#         pc_y = f"PC{i+2}"
#         pc_x_var = eigenval_df.iloc[i, 0] / eigenval_df[0].sum() * 100
#         pc_y_var = eigenval_df.iloc[i+1, 0] / eigenval_df[0].sum() * 100

#         plt.figure(figsize=(10, 8))
#         ax = sns.scatterplot(
#             data=pca_df,
#             x=pc_x, y=pc_y,
#             color="#C0C0C0",
#             s=40, alpha=0.85,
#             edgecolor='black', linewidth=0.4
#         )

#         # 添加点状图例
#         bbj_dot = Line2D(
#             [0], [0],
#             marker='o',
#             color='w',  # 图例中的背景色
#             label='BBJ',
#             markerfacecolor="#C0C0C0",
#             markeredgecolor='black',
#             markersize=8,
#             linewidth=0
#         )
#         ax.legend(
#             handles=[bbj_dot],
#             loc='upper left',
#             bbox_to_anchor=(1.02, 1),
#             frameon=False
#         )

#         plt.xlabel(f"{pc_x} ({pc_x_var:.2f}%)", fontsize=16)
#         plt.ylabel(f"{pc_y} ({pc_y_var:.2f}%)", fontsize=16)
#         plt.title(f"PCA: {pc_x} vs {pc_y}", fontsize=18, weight='bold')

#         plt.grid(True, linestyle='--', linewidth=0.5, color='gray', alpha=0.4)
#         plt.axhline(0, color='gray', linewidth=0.6, linestyle='--', alpha=0.5)
#         plt.axvline(0, color='gray', linewidth=0.6, linestyle='--', alpha=0.5)

#         plt.tight_layout()
#         pdf.savefig(dpi=300)  # Save each plot to the PDF
#         plt.close()

# %%
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.pyplot import imread

# 读取数据
eigenvec_file = f"{output_prefix}.no_high_ld.prune.pca.eigenvec"
eigenval_file = f"{output_prefix}.no_high_ld.prune.pca.eigenval"
pca_df = pd.read_csv(eigenvec_file, delim_whitespace=True)
eigenval_df = pd.read_csv(eigenval_file, delim_whitespace=True, header=None)

plt.style.use('default')
sns.set_theme(style="whitegrid", context="talk", font_scale=1.2)

with PdfPages("bbj_pca_pairwise_plots.pdf") as pdf:
    for i in range(0, 10, 2):  # PC1 vs PC2, ..., PC9 vs PC10
        pc_x = f"PC{i+1}"
        pc_y = f"PC{i+2}"
        pc_x_var = eigenval_df.iloc[i, 0] / eigenval_df[0].sum() * 100
        pc_y_var = eigenval_df.iloc[i+1, 0] / eigenval_df[0].sum() * 100

        plt.figure(figsize=(10, 8))
        ax = sns.scatterplot(
            data=pca_df,
            x=pc_x, y=pc_y,
            color="#C0C0C0",
            s=40, alpha=0.85,
            edgecolor='black', linewidth=0.4
        )

        # 添加点状图例
        bbj_dot = Line2D(
            [0], [0],
            marker='o',
            color='w',
            label='BBJ',
            markerfacecolor="#C0C0C0",
            markeredgecolor='black',
            markersize=8,
            linewidth=0
        )
        ax.legend(
            handles=[bbj_dot],
            loc='upper left',
            bbox_to_anchor=(1.02, 1),
            frameon=False
        )

        plt.xlabel(f"{pc_x} ({pc_x_var:.2f}%)", fontsize=16)
        plt.ylabel(f"{pc_y} ({pc_y_var:.2f}%)", fontsize=16)
        plt.title(f"PCA: {pc_x} vs {pc_y}", fontsize=18, weight='bold')

        plt.grid(True, linestyle='--', linewidth=0.5, color='gray', alpha=0.4)
        plt.axhline(0, color='gray', linewidth=0.6, linestyle='--', alpha=0.5)
        plt.axvline(0, color='gray', linewidth=0.6, linestyle='--', alpha=0.5)

        plt.tight_layout()

        # 保存位图 PNG，再嵌入 PDF
        tmp_png = f"tmp_pca_plot_{i}_{i+1}.png"
        plt.savefig(tmp_png, dpi=300)
        plt.close()

        fig = plt.figure(figsize=(10, 8))
        img = imread(tmp_png)
        plt.imshow(img)
        plt.axis('off')
        pdf.savefig(fig)
        plt.close()
        os.remove(tmp_png)

