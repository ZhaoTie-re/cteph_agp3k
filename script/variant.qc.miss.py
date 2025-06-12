# %%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import argparse
import os
import subprocess

os.environ["PATH"] = "/home/b/b37974/:" + os.environ["PATH"]

parser = argparse.ArgumentParser(description="Variant QC for missingness and MAF")
parser.add_argument("--bed_prefix", type=str, help="Path to the bed prefix")
parser.add_argument("--threads", type=int, default=4, help="Number of threads")
parser.add_argument("--maf_threshold", type=float, default=0.01, help="MAF threshold for classification")
parser.add_argument("--vmiss_threshold", type=float, default=0.04, help="Variant missingness threshold for filtering")
parser.add_argument("--out", type=str, help="Output prefix for the filtered plink file")
args = parser.parse_args()

bed_prefix = args.bed_prefix
variant_qc_prefix = "variant.qc"

maf_threshold = args.maf_threshold
vmiss_threshold = args.vmiss_threshold
threads = args.threads

# Calculate variant missingness using plink2
missingness_output = f"{variant_qc_prefix}.missingness"
subprocess.run([
    "/home/b/b37974/plink2",
    "--bfile", bed_prefix,
    "--missing", "variant-only",
    "--out", missingness_output, 
    "--threads", str(threads)
])

# Calculate Alt Allele Frequency (AAF) for controls 
aaf_output = f"{variant_qc_prefix}.aaf"
subprocess.run([
    "/home/b/b37974/plink2",
    "--bfile", bed_prefix,
    "--freq",
    "--keep-if", "PHENO1==1",
    "--out", aaf_output, 
    "--threads", str(threads)
])

# %%
# Stream-friendly loading
aaf_cols = ["#CHROM", "ID", "REF", "ALT", "ALT_FREQS"]
vmiss_cols = ["#CHROM", "ID", "F_MISS"]

aaf_df = pd.read_csv(f"{aaf_output}.afreq", delim_whitespace=True, usecols=aaf_cols)
vmiss_df = pd.read_csv(f"{missingness_output}.vmiss", delim_whitespace=True, usecols=vmiss_cols)

# %%
# Efficient MAF computation
aaf_df.dropna(subset=['ALT_FREQS'], inplace=True)
aaf_df['MAF'] = aaf_df['ALT_FREQS'].apply(lambda x: min(x, 1 - x))

# Merge necessary columns only
merged_df = pd.merge(aaf_df, vmiss_df, on=["#CHROM", "ID"], how="inner")

# %%
mono_variants = merged_df[merged_df["MAF"] == 0]["ID"]
merged_df = merged_df[~merged_df["ID"].isin(mono_variants)]

# Filter data based on MAF threshold
low_maf = merged_df[merged_df['MAF'] < maf_threshold]
high_maf = merged_df[merged_df['MAF'] >= maf_threshold]

# Calculate proportions and counts
def calc_stats(df):
    total = len(df)
    passed = len(df[df['F_MISS'] < vmiss_threshold])
    ratio = passed / total if total > 0 else 0
    return total, passed, ratio

low_total, low_passed, low_ratio = calc_stats(low_maf)
high_total, high_passed, high_ratio = calc_stats(high_maf)
all_total, all_passed, all_ratio = calc_stats(merged_df)

# Create subplots
fig, axes = plt.subplots(1, 3, figsize=(18, 6), sharey=False)

# Common plot function
def plot_histogram(ax, data, title, color, ratio, total, passed):
    ax.hist(data['F_MISS'], bins=50, color=color, alpha=0.7, edgecolor='black')
    ax.axvline(x=vmiss_threshold, color='red', linestyle='--', linewidth=2, label=f'Threshold = {vmiss_threshold}')
    ax.legend(loc='upper right', fontsize=10)
    ax.set_title(title, fontsize=14)
    ax.set_xlabel('Miss Rate', fontsize=12)
    ax.set_ylabel('Frequency', fontsize=12)
    ax.text(0.5, 0.85,
            f"Total: {total:,}\nPassed: {passed:,}\nRatio: {ratio:.2%}",
            transform=ax.transAxes,
            fontsize=11,
            ha='center',
            va='top',
            bbox=dict(boxstyle="round,pad=0.4", facecolor="white", edgecolor=color))
    
# Set the main title
fig.suptitle('Missing Rate across Variants', fontsize=16)

# Plot each
plot_histogram(axes[0], low_maf, f'Variants MAF < {maf_threshold}', 'blue', low_ratio, low_total, low_passed)
plot_histogram(axes[1], high_maf, f'Variants MAF ≥ {maf_threshold}', 'green', high_ratio, high_total, high_passed)
plot_histogram(axes[2], merged_df, 'All Variants', 'orange', all_ratio, all_total, all_passed)

# Style and layout
plt.tight_layout()
plt.style.use('default')
plt.subplots_adjust(top=0.85)

# Save the figure to pdf
output_pdf = f"{variant_qc_prefix}.missing_rate.pdf"
plt.savefig(output_pdf, format='pdf', dpi=300)
# Save the figure to png
output_png = f"{variant_qc_prefix}.missing_rate.png"
plt.savefig(output_png, format='png', dpi=300)
# plt.close(fig)

# %%
mono_variants.to_csv("mono_variants.txt", index=False, header=False)
poly_passed = merged_df[merged_df["F_MISS"] < vmiss_threshold]["ID"]
poly_passed.to_csv("passed_smiss.txt", index=False, header=False)

# %%
# Define the output prefix for the subset
subset_bed_prefix = args.out

# Run plink2 to extract the subset
subprocess.run([
    "/home/b/b37974/plink2",
    "--bfile", bed_prefix,
    "--extract", "passed_smiss.txt",
    "--make-bed",
    "--out", subset_bed_prefix, 
    "--threads", str(threads)
])

