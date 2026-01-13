#!/usr/bin/env python3
import argparse
import logging
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

def setup_logging(log_file=None):
    handlers = [logging.StreamHandler(sys.stdout)]
    if log_file:
        handlers.append(logging.FileHandler(log_file))
    
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=handlers
    )
    return logging.getLogger(__name__)

def read_glm_file(filepath, logger):
    logger.info(f"Loading: {filepath}")
    df = pd.read_csv(filepath, sep='\t')
    logger.info(f"  Original variants: {len(df):,}")
    return df[['ID', 'P', 'BETA', 'SE', '#CHROM', 'POS', 'REF', 'ALT']]

def clean_sumstats(data, dataset_name, logger):
    # Basic cleaning
    data = data.copy()
    for col in ['P', 'BETA', 'SE']:
        data[col] = pd.to_numeric(data[col], errors='coerce')
    
    data = data.dropna(subset=['P', 'BETA', 'SE'])
    logger.info(f"  {dataset_name}: Clean variants: {len(data):,}")
    return data

def merge_datasets(data1, data2, name1, name2, logger):
    logger.info(f"Merging {name1} and {name2}")
    
    # Rename columns
    cols = ['ID', 'P', 'BETA', 'SE', '#CHROM', 'POS', 'REF', 'ALT']
    d1 = data1[cols].rename(columns={
        'P': f'P_{name1}', 'BETA': f'BETA_{name1}', 'SE': f'SE_{name1}',
        '#CHROM': 'CHROM', 'POS': 'POS', 'REF': 'REF', 'ALT': 'ALT'
    })
    
    # For data2, we only need ID + stats to merge
    d2 = data2[['ID', 'P', 'BETA', 'SE']].rename(columns={
        'P': f'P_{name2}', 'BETA': f'BETA_{name2}', 'SE': f'SE_{name2}'
    })
    
    merged = pd.merge(d1, d2, on='ID', how='inner')
    logger.info(f"  Common variants: {len(merged):,}")
    return merged

def plot_manhattan(df, outlier_ids, title, filename, logger):
    logger.info(f"Generating Manhattan plot: {title}")
    
    # Prepare data
    df = df.copy()
    
    # Ensure CHROM is numeric for sorting/plotting
    # Handle X, Y, MT if present
    if df['#CHROM'].dtype == 'object':
        df['#CHROM'] = df['#CHROM'].replace({'X': 23, 'Y': 24, 'XY': 25, 'MT': 26, 'M': 26})
        # Remove 'chr' prefix if present
        df['#CHROM'] = df['#CHROM'].astype(str).str.replace('chr', '')
        df['#CHROM'] = pd.to_numeric(df['#CHROM'], errors='coerce')
    
    df = df.dropna(subset=['#CHROM', 'POS', 'P'])
    df['#CHROM'] = df['#CHROM'].astype(int)
    df = df.sort_values(by=['#CHROM', 'POS'])
    
    # Calculate offsets for x-axis
    chrom_ends = df.groupby('#CHROM')['POS'].max()
    chrom_offsets = {}
    current_offset = 0
    img_chroms = sorted(df['#CHROM'].unique())
    
    for chrom in img_chroms:
        chrom_offsets[chrom] = current_offset
        current_offset += chrom_ends.loc[chrom]
    
    # Vectorized offset calculation
    df['offset'] = df['#CHROM'].map(chrom_offsets)
    df['X_POS'] = df['POS'] + df['offset']
    
    df['-log10P'] = -np.log10(df['P'])
    
    # Separate outliers and normal points
    is_outlier = df['ID'].isin(outlier_ids)
    
    normal_points = df[~is_outlier]
    outlier_points = df[is_outlier]
    
    # Plotting
    plt.figure(figsize=(14, 6))
    
    # Plot background points alternating colors
    # Use simpler colors
    colors = ['#AAAAAA', '#888888'] 
    
    for i, chrom in enumerate(img_chroms):
        chrom_mask = (normal_points['#CHROM'] == chrom)
        chrom_data = normal_points[chrom_mask]
        plt.scatter(chrom_data['X_POS'], chrom_data['-log10P'], 
                    c=colors[i % 2], s=2, rasterized=True, linewidths=0)
    
    # Plot outliers
    if not outlier_points.empty:
        plt.scatter(outlier_points['X_POS'], outlier_points['-log10P'], 
                    c='red', s=25, marker='D', rasterized=True, 
                    label='Comparision Outliers', linewidths=0.5, edgecolors='black', zorder=10)

    # X axis labels
    x_ticks = []
    x_labels = []
    for chrom in img_chroms:
        mid_point = chrom_offsets[chrom] + chrom_ends[chrom] / 2
        x_ticks.append(mid_point)
        x_labels.append(str(chrom))
    
    # Only show some labels if too many chromosomes
    if len(img_chroms) > 20:
         # Show 1..22, X
        keep_labels = set(list(range(1, 23)) + [23, 24, 25, 26])
        x_ticks = [t for t, l in zip(x_ticks, x_labels) if int(l) in keep_labels]
        x_labels = [l for l in x_labels if int(l) in keep_labels]

    plt.xticks(x_ticks, x_labels, fontsize=9)
    plt.xlabel('Chromosome')
    plt.ylabel('-log10(P)')
    plt.title(title)
    
    # Add legend if outliers exist
    if not outlier_points.empty:
        plt.legend(loc='upper right')
    
    # Threshold line
    plt.axhline(y=-np.log10(5e-8), color='blue', linestyle='--', lw=0.8, alpha=0.5)
    
    # plt.tight_layout() # tight_layout was causing issues/slowness with many points
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.close()
    logger.info(f"Saved Manhattan plot to {filename}")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--wgs', required=True)
    parser.add_argument('--array', required=True)
    parser.add_argument('--out-prefix', required=True)
    parser.add_argument('--threshold', type=float, default=1.5, help='Beta difference threshold')
    
    args = parser.parse_args()
    logger = setup_logging()
    
    # Load
    wgs = read_glm_file(args.wgs, logger)
    array = read_glm_file(args.array, logger)
    
    wgs_clean = clean_sumstats(wgs, 'WGS', logger)
    array_clean = clean_sumstats(array, 'Imputed Array', logger)
    
    merged = merge_datasets(array_clean, wgs_clean, 'Imputed Array', 'WGS', logger)
    
    # Calculate difference
    merged['diff'] = merged['BETA_WGS'] - merged['BETA_Imputed Array']
    merged['abs_diff'] = merged['diff'].abs()
    
    # Filter One: Large difference
    outliers = merged[merged['abs_diff'] > args.threshold].copy()
    logger.info(f"Found {len(outliers)} outliers with |diff| > {args.threshold}")
    
    # Save outliers
    out_file = f"{args.out_prefix}.outliers.tsv"
    outliers.to_csv(out_file, sep='\t', index=False)
    logger.info(f"Saved outliers to {out_file}")
    
    # Plot Scatter
    plt.figure(figsize=(8, 8))
    
    # Plot all points
    # Downsample for scatter plot background to save space if needed
    plt.scatter(merged['BETA_Imputed Array'], merged['BETA_WGS'], c='blue', alpha=0.1, s=5, 
                label='Common Variants', rasterized=True)
    
    # Highlight outliers
    plt.scatter(outliers['BETA_Imputed Array'], outliers['BETA_WGS'], c='red', alpha=0.8, s=15, 
                label=f'Outliers (|diff| > {args.threshold})', zorder=10, rasterized=True)
    
    # Add diagonal
    lims = [
        np.min([plt.xlim(), plt.ylim()]),  # min of both axes
        np.max([plt.xlim(), plt.ylim()]),  # max of both axes
    ]
    plt.plot(lims, lims, 'k--', alpha=0.75, zorder=0)
    
    plt.xlabel('Imputed Array Beta')
    plt.ylabel('WGS Beta')
    plt.title(f'Beta Comparison (Outliers highlighted)')
    plt.legend()
    
    plt.grid(True, alpha=0.3)
    plt.savefig(f"{args.out_prefix}.plot.png", dpi=300)
    logger.info(f"Saved plot to {args.out_prefix}.plot.png")
    
    # Plot Manhattan
    outlier_ids = set(outliers['ID'])
    
    plot_manhattan(wgs_clean, outlier_ids, "WGS Manhattan", f"{args.out_prefix}.wgs_manhattan.png", logger)
    plot_manhattan(array_clean, outlier_ids, "Imputed Array Manhattan", f"{args.out_prefix}.array_manhattan.png", logger)

if __name__ == "__main__":
    main()
