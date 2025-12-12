#!/usr/bin/env python3
"""
WGS vs Array Association Results Comparison Tool

This script compares association results between WGS and Array data,
generating publication-quality comparison plots colored by MAF.

Author: ZHAO TIE
Date: December 2025
"""

import argparse
import sys
import os
from datetime import datetime
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import logging


def setup_logging(output_dir):
    """Setup logging to both file and console"""
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_file = os.path.join(output_dir, f"wgs_array_comparison_{timestamp}.log")
    
    # Create logger
    logger = logging.getLogger('WGS_Array_Comparison')
    logger.setLevel(logging.INFO)
    
    # File handler
    fh = logging.FileHandler(log_file)
    fh.setLevel(logging.INFO)
    
    # Console handler
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(logging.INFO)
    
    # Formatter
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    fh.setFormatter(formatter)
    ch.setFormatter(formatter)
    
    logger.addHandler(fh)
    logger.addHandler(ch)
    
    return logger, log_file


def parse_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description='Compare WGS and Array association results with MAF coloring',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example usage:
    python wgs_vs_array_cli.py \\
        --wgs wgs.glm.logistic \\
        --array array.glm.logistic \\
        --vqc variant_qc_sum.tsv \\
        --maf-col MAF_ALL \\
        --output-dir ./results \\
        --sig-level 5e-8
        """
    )
    
    parser.add_argument('--wgs', required=True,
                        help='Path to WGS association results file (PLINK2 .glm.logistic format)')
    parser.add_argument('--array', required=True,
                        help='Path to Array association results file (PLINK2 .glm.logistic format)')
    parser.add_argument('--vqc', required=True,
                        help='Path to variant QC summary file containing color metric data')
    parser.add_argument('--color-col', default='MAF_ALL',
                        help='Column name for coloring points (e.g., MAF_ALL, IMPUTED_R2, etc. Default: MAF_ALL)')
    parser.add_argument('--reverse-colormap', action='store_true',
                        help='Reverse colormap for continuous variables (high values = dark colors)')
    parser.add_argument('--output-dir', default='.',
                        help='Output directory for plots and logs (default: current directory)')
    parser.add_argument('--sig-level', type=float, default=5e-8,
                        help='Significance level threshold (default: 5e-8)')
    parser.add_argument('--top-n', type=int, default=10,
                        help='Number of top discordant SNPs to display (default: 10)')
    parser.add_argument('--prefix', default='wgs_array_comparison',
                        help='Prefix for output files (default: wgs_array_comparison)')
    
    return parser.parse_args()


def load_association_data(wgs_path, array_path, logger):
    """Load WGS and Array association data"""
    logger.info("=" * 70)
    logger.info("STEP 1: Loading Association Data")
    logger.info("=" * 70)
    
    logger.info(f"Reading WGS data from: {wgs_path}")
    wgs_df = pd.read_csv(wgs_path, sep='\t')
    logger.info(f"  WGS data shape: {wgs_df.shape}")
    
    logger.info(f"Reading Array data from: {array_path}")
    array_df = pd.read_csv(array_path, sep='\t')
    logger.info(f"  Array data shape: {array_df.shape}")
    
    # Extract relevant columns
    wgs_data = wgs_df[['ID', 'P', 'BETA', 'SE']].copy()
    array_data = array_df[['ID', 'P', 'BETA', 'SE']].copy()
    
    # Rename columns
    wgs_data.columns = ['SNPID', 'P_wgs', 'BETA_wgs', 'SE_wgs']
    array_data.columns = ['SNPID', 'P_array', 'BETA_array', 'SE_array']
    
    wgs_original = len(wgs_data)
    array_original = len(array_data)
    
    logger.info(f"WGS original variants:   {wgs_original:>8,}")
    logger.info(f"Array original variants: {array_original:>8,}")
    
    return wgs_data, array_data, wgs_original, array_original


def filter_invalid_pvalues(wgs_data, array_data, logger):
    """Filter variants with invalid P-values"""
    logger.info("\n" + "=" * 70)
    logger.info("STEP 2: Filtering Invalid P-values")
    logger.info("=" * 70)
    
    # Replace '.' with NaN and convert to numeric
    wgs_data['P_wgs'] = pd.to_numeric(wgs_data['P_wgs'].replace('.', np.nan), errors='coerce')
    array_data['P_array'] = pd.to_numeric(array_data['P_array'].replace('.', np.nan), errors='coerce')
    
    # Identify invalid P-values
    wgs_p_invalid = wgs_data['P_wgs'].isna() | np.isinf(wgs_data['P_wgs'])
    array_p_invalid = array_data['P_array'].isna() | np.isinf(array_data['P_array'])
    
    wgs_p_removed = wgs_p_invalid.sum()
    array_p_removed = array_p_invalid.sum()
    
    logger.info(f"WGS variants with invalid P (., NA, inf):   {wgs_p_removed:>8,}")
    logger.info(f"Array variants with invalid P (., NA, inf): {array_p_removed:>8,}")
    
    # Remove invalid variants
    wgs_data = wgs_data[~wgs_p_invalid].copy()
    array_data = array_data[~array_p_invalid].copy()
    
    logger.info(f"WGS variants remaining:   {len(wgs_data):>8,}")
    logger.info(f"Array variants remaining: {len(array_data):>8,}")
    
    return wgs_data, array_data, wgs_p_removed, array_p_removed


def filter_invalid_beta_se(wgs_data, array_data, logger):
    """Filter variants with invalid BETA/SE values"""
    logger.info("\n" + "=" * 70)
    logger.info("STEP 3: Filtering Invalid BETA/SE Values")
    logger.info("=" * 70)
    
    # Convert BETA and SE to numeric
    for col in ['BETA_wgs', 'SE_wgs']:
        wgs_data[col] = pd.to_numeric(wgs_data[col].replace('.', np.nan), errors='coerce')
    
    for col in ['BETA_array', 'SE_array']:
        array_data[col] = pd.to_numeric(array_data[col].replace('.', np.nan), errors='coerce')
    
    # Identify invalid BETA/SE
    wgs_beta_se_invalid = (
        wgs_data['BETA_wgs'].isna() | np.isinf(wgs_data['BETA_wgs']) |
        wgs_data['SE_wgs'].isna() | np.isinf(wgs_data['SE_wgs'])
    )
    
    array_beta_se_invalid = (
        array_data['BETA_array'].isna() | np.isinf(array_data['BETA_array']) |
        array_data['SE_array'].isna() | np.isinf(array_data['SE_array'])
    )
    
    wgs_beta_se_removed = wgs_beta_se_invalid.sum()
    array_beta_se_removed = array_beta_se_invalid.sum()
    
    logger.info(f"WGS variants with invalid BETA/SE (., NA, inf):   {wgs_beta_se_removed:>8,}")
    logger.info(f"Array variants with invalid BETA/SE (., NA, inf): {array_beta_se_removed:>8,}")
    
    # Remove invalid variants
    wgs_data = wgs_data[~wgs_beta_se_invalid].copy()
    array_data = array_data[~array_beta_se_invalid].copy()
    
    logger.info(f"WGS variants remaining:   {len(wgs_data):>8,}")
    logger.info(f"Array variants remaining: {len(array_data):>8,}")
    
    return wgs_data, array_data, wgs_beta_se_removed, array_beta_se_removed


def merge_and_transform(wgs_data, array_data, logger):
    """Merge datasets and calculate -log10(P)"""
    logger.info("\n" + "=" * 70)
    logger.info("STEP 4: Merging Common Variants")
    logger.info("=" * 70)
    
    # Merge on common SNPIDs
    merged_data = pd.merge(wgs_data, array_data, on='SNPID', how='inner')
    logger.info(f"Common variants after merge: {len(merged_data):>8,}")
    
    # Calculate -log10(P)
    merged_data['neglog10P_wgs'] = -np.log10(merged_data['P_wgs'])
    merged_data['neglog10P_array'] = -np.log10(merged_data['P_array'])
    
    # Check for infinite values after transformation
    inf_check = (
        np.isinf(merged_data['neglog10P_wgs']) | 
        np.isinf(merged_data['neglog10P_array'])
    )
    inf_removed = inf_check.sum()
    
    if inf_removed > 0:
        logger.info(f"Variants with inf after -log10(P) transformation: {inf_removed:>8,}")
        merged_data = merged_data[~inf_check].copy()
        logger.info(f"Final common variants: {len(merged_data):>8,}")
    
    return merged_data


def clean_column_name_for_legend(col_name):
    """Clean column name for display in legend by removing suffixes"""
    # Remove common suffixes: _CASE, _CTRL, _ALL, and trailing underscores
    import re
    cleaned = re.sub(r'_(CASE|CTRL|ALL)$', '', col_name)
    cleaned = re.sub(r'_+$', '', cleaned)  # Remove trailing underscores
    return cleaned


def load_color_data(vqc_path, color_col, merged_data, logger):
    """Load color metric data with memory optimization"""
    logger.info("\n" + "=" * 70)
    logger.info("STEP 5: Loading Color Metric Data")
    logger.info("=" * 70)
    
    logger.info(f"Reading variant QC file: {vqc_path}")
    
    # Read header first
    with open(vqc_path, 'r') as f:
        header = f.readline().strip().split('\t')
    
    logger.info(f"  Available columns: {', '.join(header)}")
    
    # Determine variant ID column
    if 'VARIANT_ID' in header:
        id_col = 'VARIANT_ID'
    elif 'SNPID' in header:
        id_col = 'SNPID'
    else:
        id_col = header[0]
        logger.warning(f"  Using first column '{id_col}' as variant ID")
    
    # Check if color column exists
    if color_col not in header:
        logger.error(f"Color column '{color_col}' not found in VQC file!")
        logger.error(f"Available columns: {', '.join(header)}")
        raise ValueError(f"Color column '{color_col}' not found")
    
    # Load only necessary columns (don't enforce dtype yet)
    logger.info(f"  Loading columns: {id_col}, {color_col}")
    vqc_subset = pd.read_csv(
        vqc_path,
        sep='\t',
        usecols=[id_col, color_col],
        dtype={id_col: str}
    )
    
    logger.info(f"  Loaded {len(vqc_subset):,} variants from VQC file")
    
    # Rename to standard column name
    vqc_subset.columns = ['SNPID', 'COLOR_METRIC']
    
    # Detect if color column is categorical or continuous
    try:
        # Try converting to numeric
        vqc_subset['COLOR_METRIC'] = pd.to_numeric(vqc_subset['COLOR_METRIC'], errors='coerce')
        is_categorical = False
        logger.info(f"  Detected {color_col} as continuous variable")
    except:
        is_categorical = True
        logger.info(f"  Detected {color_col} as categorical variable")
    
    # If most values couldn't be converted to numeric, treat as categorical
    if not is_categorical and vqc_subset['COLOR_METRIC'].isna().sum() > len(vqc_subset) * 0.5:
        # Reload as string
        vqc_subset = pd.read_csv(
            vqc_path,
            sep='\t',
            usecols=[id_col, color_col],
            dtype={id_col: str, color_col: str}
        )
        vqc_subset.columns = ['SNPID', 'COLOR_METRIC']
        is_categorical = True
        logger.info(f"  Re-detected {color_col} as categorical variable (high NA rate)")
    
    # Apply -log10 transformation for HWE p-values
    if not is_categorical and 'HWE' in color_col.upper():
        logger.info(f"  Applying -log10 transformation for HWE p-values")
        # Remove invalid values (0, negative, NaN)
        valid_mask = (vqc_subset['COLOR_METRIC'] > 0) & vqc_subset['COLOR_METRIC'].notna()
        n_invalid = (~valid_mask).sum()
        if n_invalid > 0:
            logger.info(f"  Removing {n_invalid:,} variants with invalid HWE values (≤0 or NA)")
        vqc_subset.loc[valid_mask, 'COLOR_METRIC'] = -np.log10(vqc_subset.loc[valid_mask, 'COLOR_METRIC'])
        vqc_subset.loc[~valid_mask, 'COLOR_METRIC'] = np.nan
        logger.info(f"  HWE transformed to -log10(HWE)")
    
    # Merge with association data
    logger.info("  Merging color metric data with association results...")
    merged_with_color = pd.merge(merged_data, vqc_subset, on='SNPID', how='left')
    
    # Clean up
    del vqc_subset
    del merged_data
    import gc
    gc.collect()
    
    # Statistics
    color_available = merged_with_color['COLOR_METRIC'].notna().sum()
    color_missing = merged_with_color['COLOR_METRIC'].isna().sum()
    
    logger.info(f"  SNPs with {color_col} data:    {color_available:>8,} ({color_available/len(merged_with_color)*100:.2f}%)")
    logger.info(f"  SNPs without {color_col} data: {color_missing:>8,} ({color_missing/len(merged_with_color)*100:.2f}%)")
    
    if color_available > 0:
        if is_categorical:
            categories = merged_with_color['COLOR_METRIC'].value_counts()
            logger.info(f"  {color_col} categories: {len(categories)}")
            for cat, count in categories.items():
                logger.info(f"    {cat}: {count:,} ({count/color_available*100:.2f}%)")
        else:
            logger.info(f"  {color_col} range:  {merged_with_color['COLOR_METRIC'].min():.4f} - {merged_with_color['COLOR_METRIC'].max():.4f}")
            logger.info(f"  {color_col} mean:   {merged_with_color['COLOR_METRIC'].mean():.4f}")
            logger.info(f"  {color_col} median: {merged_with_color['COLOR_METRIC'].median():.4f}")
    
    return merged_with_color, color_col, is_categorical


def generate_plots(merged_with_color, color_col, is_categorical, reverse_cmap, sig_level, output_prefix, logger):
    """Generate comparison plots"""
    logger.info("\n" + "=" * 70)
    logger.info("STEP 6: Generating Plots")
    logger.info("=" * 70)
    
    # Clean column name for legend display
    legend_label = clean_column_name_for_legend(color_col)
    
    # Add -log10 prefix for HWE columns
    if 'HWE' in color_col.upper() and not is_categorical:
        legend_label = f"-log$_{{10}}$({legend_label})"
    
    # Set publication-quality parameters
    plt.rcParams['font.family'] = 'Arial'
    plt.rcParams['font.size'] = 10
    plt.rcParams['axes.linewidth'] = 1
    
    # Create figure
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    
    # ========== Left panel: P-value comparison ==========
    max_val_p = max(merged_with_color['neglog10P_wgs'].max(), merged_with_color['neglog10P_array'].max())
    padding_p = max_val_p * 0.1
    limit_p = max_val_p + padding_p
    
    # Prepare color mapping
    if is_categorical:
        # Categorical coloring
        categories = merged_with_color['COLOR_METRIC'].dropna().unique()
        n_categories = len(categories)
        logger.info(f"  Using categorical color scheme with {n_categories} categories")
        
        # Create high-contrast color map for categories
        if n_categories <= 8:
            # Set1: Maximum contrast for up to 8 categories
            cmap = plt.cm.get_cmap('Set1', n_categories)
        elif n_categories <= 12:
            # Paired: High contrast paired colors for up to 12 categories
            cmap = plt.cm.get_cmap('Paired', n_categories)
        else:
            # tab20: For more categories, use tab20
            cmap = plt.cm.get_cmap('tab20', min(n_categories, 20))
        
        # Map categories to integers
        cat_to_int = {cat: i for i, cat in enumerate(sorted(categories))}
        color_values = merged_with_color['COLOR_METRIC'].map(cat_to_int)
        
        scatter1 = ax1.scatter(
            merged_with_color['neglog10P_wgs'], merged_with_color['neglog10P_array'],
            c=color_values, cmap=cmap,
            alpha=0.5, s=15, edgecolors='none', rasterized=True, zorder=3,
            vmin=-0.5, vmax=n_categories-0.5
        )
    else:
        # Continuous coloring
        color_min = merged_with_color['COLOR_METRIC'].min()
        color_max = merged_with_color['COLOR_METRIC'].max()
        cmap_name = 'plasma_r' if reverse_cmap else 'plasma'
        logger.info(f"  Using continuous color scheme (range: {color_min:.4f} - {color_max:.4f})")
        logger.info(f"  Colormap: {cmap_name} ({'reversed' if reverse_cmap else 'normal'})")
        
        scatter1 = ax1.scatter(
            merged_with_color['neglog10P_wgs'], merged_with_color['neglog10P_array'],
            c=merged_with_color['COLOR_METRIC'], cmap=cmap_name,
            alpha=0.5, s=15, edgecolors='none', rasterized=True, zorder=3,
            vmin=color_min, vmax=color_max
        )
    
    ax1.set_xlim(-padding_p * 0.5, limit_p)
    ax1.set_ylim(-padding_p * 0.5, limit_p)
    
    # Diagonal line
    xlim = ax1.get_xlim()
    ylim = ax1.get_ylim()
    diag_min = max(xlim[0], ylim[0])
    diag_max = min(xlim[1], ylim[1])
    ax1.plot([diag_min, diag_max], [diag_min, diag_max],
            color='#E63946', linestyle='--', lw=1.2, alpha=0.8, zorder=2)
    
    # Significance threshold
    sig_threshold = -np.log10(sig_level)
    ax1.axhline(y=sig_threshold, color='gray', linestyle=':', alpha=0.6, linewidth=1, zorder=1)
    ax1.axvline(x=sig_threshold, color='gray', linestyle=':', alpha=0.6, linewidth=1, zorder=1)
    
    ax1.set_xlabel('WGS -log$_{10}$(P)', fontsize=11)
    ax1.set_ylabel('Imputed Array -log$_{10}$(P)', fontsize=11)
    ax1.set_title(f'P-value Comparison (Colored by {legend_label})', fontsize=12, fontweight='bold', pad=10)
    ax1.grid(True, alpha=0.15, linestyle='-', linewidth=0.5, zorder=0)
    ax1.tick_params(axis='both', which='major', labelsize=10)
    
    # Correlation
    corr_p = np.corrcoef(merged_with_color['neglog10P_wgs'], merged_with_color['neglog10P_array'])[0, 1]
    ax1.text(0.98, 0.02, f'r = {corr_p:.3f}',
            transform=ax1.transAxes, fontsize=10,
            verticalalignment='bottom', horizontalalignment='right',
            bbox=dict(boxstyle='round,pad=0.4', facecolor='white',
                     edgecolor='gray', linewidth=0.8, alpha=0.9))
    
    # ========== Right panel: BETA comparison ==========
    max_val_beta = max(abs(merged_with_color['BETA_wgs']).max(), abs(merged_with_color['BETA_array']).max())
    padding_beta = max_val_beta * 0.15
    limit_beta = max_val_beta + padding_beta
    
    if is_categorical:
        scatter2 = ax2.scatter(
            merged_with_color['BETA_wgs'], merged_with_color['BETA_array'],
            c=color_values, cmap=cmap,
            alpha=0.45, s=12, edgecolors='none', rasterized=True, zorder=3,
            vmin=-0.5, vmax=n_categories-0.5
        )
    else:
        scatter2 = ax2.scatter(
            merged_with_color['BETA_wgs'], merged_with_color['BETA_array'],
            c=merged_with_color['COLOR_METRIC'], cmap=cmap_name,
            alpha=0.45, s=12, edgecolors='none', rasterized=True, zorder=3,
            vmin=color_min, vmax=color_max
        )
    
    ax2.set_xlim(-limit_beta, limit_beta)
    ax2.set_ylim(-limit_beta, limit_beta)
    
    # Diagonal line
    xlim = ax2.get_xlim()
    ylim = ax2.get_ylim()
    diag_min = max(xlim[0], ylim[0])
    diag_max = min(xlim[1], ylim[1])
    ax2.plot([diag_min, diag_max], [diag_min, diag_max],
            color='#E63946', linestyle='--', lw=1.2, alpha=0.8, zorder=2)
    
    # Zero lines
    ax2.axhline(y=0, color='black', linestyle='-', alpha=0.25, linewidth=0.8, zorder=1)
    ax2.axvline(x=0, color='black', linestyle='-', alpha=0.25, linewidth=0.8, zorder=1)
    
    ax2.set_xlabel('WGS β', fontsize=11)
    ax2.set_ylabel('Imputed Array β', fontsize=11)
    ax2.set_title(f'Effect Size (β) Comparison (Colored by {legend_label})', fontsize=12, fontweight='bold', pad=10)
    ax2.grid(True, alpha=0.15, linestyle='-', linewidth=0.5, zorder=0)
    ax2.tick_params(axis='both', which='major', labelsize=10)
    
    # Correlation
    corr_beta = np.corrcoef(merged_with_color['BETA_wgs'], merged_with_color['BETA_array'])[0, 1]
    ax2.text(0.98, 0.02, f'r = {corr_beta:.3f}',
            transform=ax2.transAxes, fontsize=10,
            verticalalignment='bottom', horizontalalignment='right',
            bbox=dict(boxstyle='round,pad=0.4', facecolor='white',
                     edgecolor='gray', linewidth=0.8, alpha=0.9))
    
    # Set equal aspect ratio
    ax1.set_aspect('equal', adjustable='box')
    ax2.set_aspect('equal', adjustable='box')
    
    # Adjust for colorbar
    plt.subplots_adjust(left=0.07, right=0.86, bottom=0.12, top=0.92, wspace=0.30)
    
    # Add colorbar
    cbar_ax = fig.add_axes([0.88, 0.20, 0.018, 0.60])
    cbar = fig.colorbar(scatter2, cax=cbar_ax)
    cbar.set_label(legend_label, fontsize=11.5, labelpad=15, rotation=270, va='bottom')
    cbar.ax.tick_params(labelsize=10, length=5, width=1)
    cbar.outline.set_linewidth(1)
    
    if is_categorical:
        # Categorical colorbar with category labels
        cbar.set_ticks([i for i in range(n_categories)])
        cbar.set_ticklabels([cat for cat in sorted(categories)])
    else:
        # Continuous colorbar with automatic ticks
        cbar.locator = ticker.MaxNLocator(nbins=6)
        cbar.update_ticks()
    
    # Save plot as PNG (with color column name in filename)
    png_path = f"{output_prefix}_{color_col}.png"
    
    logger.info(f"  Saving PNG: {png_path}")
    plt.savefig(png_path, dpi=300, bbox_inches='tight')
    
    plt.close()
    
    # Reset matplotlib parameters
    plt.rcParams.update(plt.rcParamsDefault)
    
    logger.info(f"  P-value correlation: {corr_p:.4f}")
    logger.info(f"  BETA correlation:    {corr_beta:.4f}")
    
    return corr_p, corr_beta


def print_summary(wgs_original, array_original, wgs_p_removed, array_p_removed,
                 wgs_beta_se_removed, array_beta_se_removed, wgs_final, array_final,
                 common_variants, merged_with_color, corr_p, corr_beta, top_n, output_prefix, logger):
    """Print comprehensive summary"""
    logger.info("\n" + "=" * 70)
    logger.info("FILTERING SUMMARY")
    logger.info("=" * 70)
    
    logger.info("\nWGS Variants:")
    logger.info(f"  Original:                  {wgs_original:>8,}")
    logger.info(f"  Removed (invalid P):       {wgs_p_removed:>8,} ({wgs_p_removed/wgs_original*100:>5.2f}%)")
    logger.info(f"  Removed (invalid BETA/SE): {wgs_beta_se_removed:>8,} ({wgs_beta_se_removed/(wgs_original-wgs_p_removed)*100:>5.2f}%)")
    logger.info(f"  Final clean variants:      {wgs_final:>8,} ({wgs_final/wgs_original*100:>5.2f}%)")
    
    logger.info("\nArray Variants:")
    logger.info(f"  Original:                  {array_original:>8,}")
    logger.info(f"  Removed (invalid P):       {array_p_removed:>8,} ({array_p_removed/array_original*100:>5.2f}%)")
    logger.info(f"  Removed (invalid BETA/SE): {array_beta_se_removed:>8,} ({array_beta_se_removed/(array_original-array_p_removed)*100:>5.2f}%)")
    logger.info(f"  Final clean variants:      {array_final:>8,} ({array_final/array_original*100:>5.2f}%)")
    
    logger.info("\nCommon Variants:")
    logger.info(f"  Overlapping variants:      {common_variants:>8,}")
    logger.info(f"  WGS coverage:              {common_variants/wgs_final*100:>5.2f}% of clean WGS variants")
    logger.info(f"  Array coverage:            {common_variants/array_final*100:>5.2f}% of clean Array variants")
    
    logger.info("\n" + "=" * 70)
    logger.info("DATA RANGES")
    logger.info("=" * 70)
    logger.info(f"P-value range (WGS):   {merged_with_color['P_wgs'].min():.2e} - {merged_with_color['P_wgs'].max():.2e}")
    logger.info(f"P-value range (Array): {merged_with_color['P_array'].min():.2e} - {merged_with_color['P_array'].max():.2e}")
    logger.info(f"BETA range (WGS):      {merged_with_color['BETA_wgs'].min():>7.4f} - {merged_with_color['BETA_wgs'].max():>7.4f}")
    logger.info(f"BETA range (Array):    {merged_with_color['BETA_array'].min():>7.4f} - {merged_with_color['BETA_array'].max():>7.4f}")
    
    logger.info("\n" + "=" * 70)
    logger.info("CORRELATION STATISTICS")
    logger.info("=" * 70)
    logger.info(f"-log10(P) Pearson correlation: {corr_p:.4f}")
    logger.info(f"BETA Pearson correlation:      {corr_beta:.4f}")
    
    # Top discordant SNPs
    logger.info("\n" + "=" * 70)
    logger.info(f"TOP {top_n} SNPS WITH LARGEST P-VALUE DIFFERENCE")
    logger.info("=" * 70)
    
    merged_with_color['neglogP_diff'] = np.abs(merged_with_color['neglog10P_wgs'] - merged_with_color['neglog10P_array'])
    merged_with_color['BETA_diff'] = np.abs(merged_with_color['BETA_wgs'] - merged_with_color['BETA_array'])
    
    top_discordant = merged_with_color.nlargest(top_n, 'neglogP_diff')[
        ['SNPID', 'P_wgs', 'P_array', 'BETA_wgs', 'BETA_array', 'neglogP_diff', 'BETA_diff']
    ]
    
    # Save to TSV file
    top_discordant_file = f"{output_prefix}_top{top_n}_discordant.tsv"
    top_discordant.to_csv(top_discordant_file, sep='\t', index=False, float_format='%.6g')
    logger.info(f"\nTop {top_n} discordant SNPs saved to: {top_discordant_file}")
    
    logger.info("\n" + top_discordant.to_string(index=False))


def main():
    """Main execution function"""
    # Parse arguments
    args = parse_arguments()
    
    # Create output directory if needed
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Setup logging
    logger, log_file = setup_logging(args.output_dir)
    
    logger.info("=" * 70)
    logger.info("WGS vs Array Association Comparison")
    logger.info("=" * 70)
    logger.info(f"WGS file:     {args.wgs}")
    logger.info(f"Array file:   {args.array}")
    logger.info(f"VQC file:     {args.vqc}")
    logger.info(f"Color column: {args.color_col}")
    logger.info(f"Output dir:   {args.output_dir}")
    logger.info(f"Sig level:    {args.sig_level}")
    logger.info(f"Top N SNPs:   {args.top_n}")
    logger.info(f"Log file:     {log_file}")
    
    try:
        # Step 1: Load data
        wgs_data, array_data, wgs_original, array_original = load_association_data(
            args.wgs, args.array, logger
        )
        
        # Step 2: Filter P-values
        wgs_data, array_data, wgs_p_removed, array_p_removed = filter_invalid_pvalues(
            wgs_data, array_data, logger
        )
        
        # Step 3: Filter BETA/SE
        wgs_data, array_data, wgs_beta_se_removed, array_beta_se_removed = filter_invalid_beta_se(
            wgs_data, array_data, logger
        )
        
        wgs_final = len(wgs_data)
        array_final = len(array_data)
        
        # Step 4: Merge and transform
        merged_data = merge_and_transform(wgs_data, array_data, logger)
        common_variants = len(merged_data)
        
        # Step 5: Load color metric data
        merged_with_color, color_col, is_categorical = load_color_data(args.vqc, args.color_col, merged_data, logger)
        
        # Step 6: Generate plots
        output_prefix = os.path.join(args.output_dir, args.prefix)
        corr_p, corr_beta = generate_plots(merged_with_color, color_col, is_categorical, args.reverse_colormap, args.sig_level, output_prefix, logger)
        
        # Print summary
        print_summary(
            wgs_original, array_original, wgs_p_removed, array_p_removed,
            wgs_beta_se_removed, array_beta_se_removed, wgs_final, array_final,
            common_variants, merged_with_color, corr_p, corr_beta, args.top_n, output_prefix, logger
        )
        
        logger.info("\n" + "=" * 70)
        logger.info("ANALYSIS COMPLETED SUCCESSFULLY")
        logger.info("=" * 70)
        
        return 0
        
    except Exception as e:
        logger.error(f"\nERROR: {str(e)}")
        logger.exception("Full traceback:")
        return 1


if __name__ == '__main__':
    sys.exit(main())
