# %%
# usr/bin/env python3
# -*- coding: utf-8 -*-
# /LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/script/variant.parameter.py
import pandas as pd
import numpy as np
# import pysam
import matplotlib.pyplot as plt
import subprocess
import argparse

# setting up the argument parser
parser = argparse.ArgumentParser(description="Plot the distribution of Mapping Quality (MQ) and Variant Quality Score Recalibration (VQSLOD) for a given chromosome")
parser.add_argument("--chrom", type=str, help="Chromosome number")
parser.add_argument("--vcf_path", type=str, help="Path to the VCF file")
parser.add_argument("--MQ", type=float, default=58.75, help="Threshold for Mapping Quality (MQ)")
parser.add_argument("--VQSLOD", type=float, default=10, help="Threshold for Variant Quality Score Recalibration (VQSLOD)")
parser.add_argument("--output_pdf", type=str, default="output.pdf", help="Output PDF file name")
args = parser.parse_args()

chrom = args.chrom
vcf_path= args.vcf_path
# vcf_file = pysam.VariantFile(vcf_path)
# samples = list(vcf_file.header.samples)
# print(f"Number of samples: {len(samples)}")

# %%

command = [
    "bcftools",
    "query",
    "-f", "%CHROM\t%POS\t%REF\t%ALT\t%INFO/MQ\t%INFO/VQSLOD\n",
    vcf_path
]

output_file = f"{chrom}.MQ.VQSLOD.txt"

with open(output_file, 'w') as file:
    result = subprocess.run(command, stdout=file, text=True)

if result.returncode == 0:
    print("Command Success:", output_file)
else:
    print("Command Failed:", result.returncode)

# %%
MQ_VQSLOD_df = pd.read_csv(output_file, sep="\s+", header=None)
# Rename columns
MQ_VQSLOD_df.columns = ["CHROM", "POS", "REF", "ALT", "MQ", "VQSLOD"]
# Convert columns to numeric to avoid errors
MQ_VQSLOD_df["MQ"] = pd.to_numeric(MQ_VQSLOD_df["MQ"], errors='coerce')
MQ_VQSLOD_df["VQSLOD"] = pd.to_numeric(MQ_VQSLOD_df["VQSLOD"], errors='coerce')

# %%

plt.style.use('default')

# set parameters
MQ_threshold = args.MQ
VQSLOD_threshold = args.VQSLOD

# subplots with 2 rows and 1 column
fig, axes = plt.subplots(2, 1, figsize=(10, 8))
fig.subplots_adjust(hspace=0.4) # adjust the space between the plots

# 1st plot: MQ distribution``
axes[0].hist(MQ_VQSLOD_df["MQ"], bins=100, edgecolor='black', linewidth=0.5, color='steelblue', alpha=0.7)
axes[0].axvline(x=MQ_threshold, color='crimson', linestyle='--', linewidth=2, label=f'MQ = {MQ_threshold}')

num_sites_gt_thres = (MQ_VQSLOD_df["MQ"] > MQ_threshold).sum()
mq_min, mq_max = MQ_VQSLOD_df["MQ"].min(), MQ_VQSLOD_df["MQ"].max()
total_sites = MQ_VQSLOD_df.shape[0]
percentage_gt_thres = (num_sites_gt_thres / total_sites) * 100

axes[0].text(mq_max + 2.5, axes[0].get_ylim()[1] * 0.9, f'Variants MQ>{MQ_threshold}: {num_sites_gt_thres:,}({percentage_gt_thres:.2f}%)', color='crimson', fontsize=12, horizontalalignment='left')
axes[0].text(mq_max + 2.5, axes[0].get_ylim()[1] * 0.8, f'Total Variants: {total_sites:,}', color='darkgreen', fontsize=12, horizontalalignment='left')
axes[0].text(mq_max + 2.5, axes[0].get_ylim()[1] * 0.7, f'MQ Range: {mq_min:.2f} ~ {mq_max:.2f}', color='navy', fontsize=12, horizontalalignment='left')

axes[0].set_xlabel('Mapping Quality (MQ)', fontsize=14)
axes[0].set_ylabel('Frequency', fontsize=14)
axes[0].set_title(f'Distribution of {chrom} MQ', fontsize=16)
axes[0].legend(fontsize=12, loc='upper left')

# 2nd plot: VQSLOD distribution
axes[1].hist(MQ_VQSLOD_df["VQSLOD"], bins=100, edgecolor='black', linewidth=0.5, color='steelblue', alpha=0.7)
axes[1].axvline(x=VQSLOD_threshold, color='crimson', linestyle='--', linewidth=2, label=f'VQSLOD = {VQSLOD_threshold}')

num_sites_vq_thres = (MQ_VQSLOD_df["VQSLOD"] > VQSLOD_threshold).sum()
vq_min, vq_max = MQ_VQSLOD_df["VQSLOD"].min(), MQ_VQSLOD_df["VQSLOD"].max()
percentage_vq_thres = (num_sites_vq_thres / total_sites) * 100

axes[1].text(vq_max + 2.5, axes[1].get_ylim()[1] * 0.9, f'Variants VQSLOD>{VQSLOD_threshold}: {num_sites_vq_thres:,}({percentage_vq_thres:.2f}%)', color='crimson', fontsize=12, horizontalalignment='left')
axes[1].text(vq_max + 2.5, axes[1].get_ylim()[1] * 0.8, f'Total Variants: {total_sites:,}', color='darkgreen', fontsize=12, horizontalalignment='left')
axes[1].text(vq_max + 2.5, axes[1].get_ylim()[1] * 0.7, f'VQSLOD Range: {vq_min:.2f} - {vq_max:.2f}', color='navy', fontsize=12, horizontalalignment='left')

axes[1].set_xlabel('Variant Quality Score Recalibration (VQSLOD)', fontsize=14)
axes[1].set_ylabel('Frequency', fontsize=14)
axes[1].set_title(f'Distribution of {chrom} VQSLOD', fontsize=16)
axes[1].legend(fontsize=12, loc='upper left')

# plt.show()
# save the plot into pdf
output_pdf = args.output_pdf
plt.savefig(output_pdf, dpi=300, bbox_inches='tight')
print(f"Plot saved: {output_pdf}")
