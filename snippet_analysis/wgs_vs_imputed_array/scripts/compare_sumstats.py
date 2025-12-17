#!/usr/bin/env python3
"""
Compare two GWAS summary statistics with optional variant-level annotation coloring.
"""

import argparse
import logging
import sys
import gzip
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from pathlib import Path
from multiprocessing import Pool, cpu_count
from functools import partial
import subprocess
import os


def setup_logging(log_file):
    """Setup logging configuration."""
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler(sys.stdout)
        ]
    )
    return logging.getLogger(__name__)


def read_glm_file(filepath, logger):
    """Read PLINK2 .glm.logistic file in chunks for memory efficiency."""
    logger.info(f"Loading: {filepath}")
    
    # Read the file
    df = pd.read_csv(filepath, sep='\t')
    
    logger.info(f"  Original variants: {len(df):,}")
    
    # Extract necessary columns
    data = df[['ID', 'P', 'BETA', 'SE']].copy()
    
    return data


def clean_sumstats(data, dataset_name, logger):
    """Clean summary statistics by removing missing/invalid values."""
    original_count = len(data)
    
    # Replace '.' with NaN
    for col in data.columns:
        if col != 'ID':
            data[col] = data[col].replace('.', np.nan)
            data[col] = pd.to_numeric(data[col], errors='coerce')
    
    # Remove invalid P-values
    p_invalid = data['P'].isna() | np.isinf(data['P'])
    p_removed = p_invalid.sum()
    data = data[~p_invalid].copy()
    
    logger.info(f"  {dataset_name}: Removed {p_removed:,} variants with invalid P-values")
    
    # Remove invalid BETA/SE
    beta_se_invalid = (
        data['BETA'].isna() | np.isinf(data['BETA']) |
        data['SE'].isna() | np.isinf(data['SE'])
    )
    beta_se_removed = beta_se_invalid.sum()
    data = data[~beta_se_invalid].copy()
    
    logger.info(f"  {dataset_name}: Removed {beta_se_removed:,} variants with invalid BETA/SE")
    logger.info(f"  {dataset_name}: Clean variants: {len(data):,} ({len(data)/original_count*100:.2f}%)")
    
    return data


def read_variant_stats(filepath, variant_ids, columns_to_read, tabix_path, logger):
    """Read variant statistics for specific variants using tabix for fast lookup."""
    logger.info(f"Loading variant statistics from: {filepath}")
    logger.info(f"  Columns to extract: {', '.join(columns_to_read)}")
    logger.info(f"  Using tabix: {tabix_path}")
    
    # Check if tabix index exists
    tbi_file = f"{filepath}.tbi"
    if not os.path.exists(tbi_file):
        logger.warning(f"  Tabix index not found: {tbi_file}")
        logger.info(f"  Falling back to sequential reading")
        return read_variant_stats_sequential(filepath, variant_ids, columns_to_read, logger)
    
    # Check if tabix executable exists
    if not os.path.exists(tabix_path):
        logger.warning(f"  Tabix executable not found: {tabix_path}")
        logger.info(f"  Falling back to sequential reading")
        return read_variant_stats_sequential(filepath, variant_ids, columns_to_read, logger)
    
    # Convert variant IDs to dictionary for fast lookup
    # Assume variant ID format: chr:pos:ref:alt or chr_pos_ref_alt
    variant_dict = {}
    for var_id in variant_ids:
        # Parse variant ID to get chromosome and position
        parts = var_id.replace('_', ':').split(':')
        if len(parts) >= 2:
            chrom = parts[0].replace('chr', '')
            try:
                pos = int(parts[1])
                variant_dict[var_id] = (chrom, pos)
            except ValueError:
                continue
    
    logger.info(f"  Parsed {len(variant_dict)} variant positions")
    
    # Group variants by chromosome for efficient querying
    chrom_variants = {}
    for var_id, (chrom, pos) in variant_dict.items():
        if chrom not in chrom_variants:
            chrom_variants[chrom] = []
        chrom_variants[chrom].append((pos, var_id))
    
    # Sort positions within each chromosome
    for chrom in chrom_variants:
        chrom_variants[chrom].sort()
    
    logger.info(f"  Variants span {len(chrom_variants)} chromosomes")
    
    # Read header to get column indices
    with gzip.open(filepath, 'rt') as f:
        header_line = f.readline().strip()
        all_columns = header_line.split('\t')
        
        # Get indices for columns we need
        col_indices = {}
        for col in ['ID'] + columns_to_read:
            if col in all_columns:
                col_indices[col] = all_columns.index(col)
            else:
                logger.warning(f"  Column {col} not found in file")
        
        if 'ID' not in col_indices:
            logger.error("  ID column not found!")
            return pd.DataFrame(columns=['ID'] + columns_to_read)
    
    # Query tabix for each chromosome region
    results = []
    for chrom, variants in chrom_variants.items():
        if not variants:
            continue
        
        # Get min and max positions for this chromosome
        min_pos = variants[0][0]
        max_pos = variants[-1][0]
        
        # Create a set of variant IDs for this chromosome for fast lookup
        chrom_var_ids = set([var_id for _, var_id in variants])
        
        # Try both with and without chr prefix
        query_success = False
        for chrom_fmt in [chrom, f"chr{chrom}"]:
            try:
                # Run tabix command
                cmd = [tabix_path, filepath, f"{chrom_fmt}:{min_pos}-{max_pos}"]
                result = subprocess.run(cmd, capture_output=True, text=True, check=True)
                
                # Process output lines
                for line in result.stdout.strip().split('\n'):
                    if not line:
                        continue
                    
                    fields = line.split('\t')
                    if len(fields) <= col_indices['ID']:
                        continue
                    
                    var_id = fields[col_indices['ID']]
                    
                    if var_id in chrom_var_ids:
                        # Extract required columns
                        row_data = {'ID': var_id}
                        for col in columns_to_read:
                            if col in col_indices and col_indices[col] < len(fields):
                                row_data[col] = fields[col_indices[col]]
                        results.append(row_data)
                
                query_success = True
                break
                
            except subprocess.CalledProcessError:
                continue
            except Exception as e:
                logger.warning(f"  Error querying chromosome {chrom_fmt}: {e}")
                continue
        
        if not query_success:
            logger.warning(f"  Could not query chromosome {chrom}")
    
    if results:
        variant_stats = pd.DataFrame(results)
        logger.info(f"  Loaded statistics for {len(variant_stats):,} variants using tabix")
    else:
        logger.warning("  No matching variants found")
        variant_stats = pd.DataFrame(columns=['ID'] + columns_to_read)
    
    return variant_stats


def read_variant_stats_sequential(filepath, variant_ids, columns_to_read, logger):
    """Fallback: Read variant statistics sequentially (original method)."""
    logger.info(f"  Using sequential reading method")
    
    chunk_size = 100000
    chunks = []
    
    with gzip.open(filepath, 'rt') as f:
        for chunk in pd.read_csv(f, sep='\t', chunksize=chunk_size):
            chunk_filtered = chunk[chunk['ID'].isin(variant_ids)]
            if len(chunk_filtered) > 0:
                chunks.append(chunk_filtered[['ID'] + columns_to_read])
    
    if chunks:
        variant_stats = pd.concat(chunks, ignore_index=True)
        logger.info(f"  Loaded statistics for {len(variant_stats):,} variants")
    else:
        logger.warning("  No matching variants found")
        variant_stats = pd.DataFrame(columns=['ID'] + columns_to_read)
    
    return variant_stats


def merge_datasets(data1, data2, dataset1_name, dataset2_name, logger):
    """Merge two cleaned datasets on common variant IDs."""
    logger.info(f"Merging {dataset1_name} and {dataset2_name}")
    
    # Rename columns
    data1.columns = ['SNPID', f'P_{dataset1_name}', f'BETA_{dataset1_name}', f'SE_{dataset1_name}']
    data2.columns = ['SNPID', f'P_{dataset2_name}', f'BETA_{dataset2_name}', f'SE_{dataset2_name}']
    
    # Merge on SNPID
    merged = pd.merge(data1, data2, on='SNPID', how='inner')
    
    # Calculate -log10(P)
    merged[f'neglog10P_{dataset1_name}'] = -np.log10(merged[f'P_{dataset1_name}'])
    merged[f'neglog10P_{dataset2_name}'] = -np.log10(merged[f'P_{dataset2_name}'])
    
    # Check for infinite values after transformation
    inf_check = (
        np.isinf(merged[f'neglog10P_{dataset1_name}']) | 
        np.isinf(merged[f'neglog10P_{dataset2_name}'])
    )
    inf_removed = inf_check.sum()
    
    if inf_removed > 0:
        logger.warning(f"  Removed {inf_removed:,} variants with inf after -log10(P) transformation")
        merged = merged[~inf_check].copy()
    
    logger.info(f"  Common variants: {len(merged):,}")
    
    return merged


def plot_basic_comparison(merged, dataset1_name, dataset2_name, output_file, sig_level, logger):
    """Plot basic P-value and BETA comparison without coloring."""
    logger.info(f"Generating basic comparison plot")
    
    # Set publication-quality parameters
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = ['DejaVu Sans', 'Bitstream Vera Sans', 'Computer Modern Sans Serif', 'Lucida Grande', 'Verdana', 'Geneva', 'Lucid', 'Arial', 'Helvetica', 'Avant Garde', 'sans-serif']
    plt.rcParams['font.size'] = 10
    plt.rcParams['axes.linewidth'] = 1
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5.5))
    
    # === Left panel: -log10(P) comparison ===
    max_val_p = max(merged[f'neglog10P_{dataset1_name}'].max(), 
                    merged[f'neglog10P_{dataset2_name}'].max())
    padding_p = max_val_p * 0.1
    limit_p = max_val_p + padding_p
    
    ax1.scatter(merged[f'neglog10P_{dataset1_name}'], merged[f'neglog10P_{dataset2_name}'],
               alpha=0.5, s=10, c='#2E86AB', edgecolors='none', rasterized=True, zorder=3)
    
    ax1.set_xlim(-padding_p * 0.5, limit_p)
    ax1.set_ylim(-padding_p * 0.5, limit_p)
    
    # Diagonal line
    diag_min = 0
    diag_max = limit_p
    ax1.plot([diag_min, diag_max], [diag_min, diag_max],
            color='#E63946', linestyle='--', lw=1.2, alpha=0.8, zorder=2)
    
    # Significance threshold
    sig_threshold = -np.log10(sig_level)
    ax1.axhline(y=sig_threshold, color='gray', linestyle=':', alpha=0.6, linewidth=1, zorder=1)
    ax1.axvline(x=sig_threshold, color='gray', linestyle=':', alpha=0.6, linewidth=1, zorder=1)
    
    ax1.set_xlabel(f'{dataset1_name} -log$_{{10}}$(P)', fontsize=11)
    ax1.set_ylabel(f'{dataset2_name} -log$_{{10}}$(P)', fontsize=11)
    ax1.set_title('P-value Comparison', fontsize=12, fontweight='bold', pad=10)
    ax1.grid(True, alpha=0.15, linestyle='-', linewidth=0.5, zorder=0)
    ax1.tick_params(axis='both', which='major', labelsize=10)
    
    corr_p = np.corrcoef(merged[f'neglog10P_{dataset1_name}'], 
                         merged[f'neglog10P_{dataset2_name}'])[0, 1]
    ax1.text(0.98, 0.02, f'r = {corr_p:.3f}',
            transform=ax1.transAxes, fontsize=10,
            verticalalignment='bottom', horizontalalignment='right',
            bbox=dict(boxstyle='round,pad=0.4', facecolor='white',
                     edgecolor='gray', linewidth=0.8, alpha=0.9))
    
    # === Right panel: BETA comparison ===
    max_val_beta = max(abs(merged[f'BETA_{dataset1_name}']).max(),
                       abs(merged[f'BETA_{dataset2_name}']).max())
    padding_beta = max_val_beta * 0.15
    limit_beta = max_val_beta + padding_beta
    
    ax2.scatter(merged[f'BETA_{dataset1_name}'], merged[f'BETA_{dataset2_name}'],
               alpha=0.35, s=8, c='#2E86AB', edgecolors='none', rasterized=True, zorder=3)
    
    ax2.set_xlim(-limit_beta, limit_beta)
    ax2.set_ylim(-limit_beta, limit_beta)
    
    # Diagonal line
    diag_min = -limit_beta
    diag_max = limit_beta
    ax2.plot([diag_min, diag_max], [diag_min, diag_max],
            color='#E63946', linestyle='--', lw=1.2, alpha=0.8, zorder=2)
    
    # Zero lines
    ax2.axhline(y=0, color='black', linestyle='-', alpha=0.25, linewidth=0.8, zorder=1)
    ax2.axvline(x=0, color='black', linestyle='-', alpha=0.25, linewidth=0.8, zorder=1)
    
    ax2.set_xlabel(f'{dataset1_name} β', fontsize=11)
    ax2.set_ylabel(f'{dataset2_name} β', fontsize=11)
    ax2.set_title('Effect Size (β) Comparison', fontsize=12, fontweight='bold', pad=10)
    ax2.grid(True, alpha=0.15, linestyle='-', linewidth=0.5, zorder=0)
    ax2.tick_params(axis='both', which='major', labelsize=10)
    
    corr_beta = np.corrcoef(merged[f'BETA_{dataset1_name}'],
                           merged[f'BETA_{dataset2_name}'])[0, 1]
    ax2.text(0.98, 0.02, f'r = {corr_beta:.3f}',
            transform=ax2.transAxes, fontsize=10,
            verticalalignment='bottom', horizontalalignment='right',
            bbox=dict(boxstyle='round,pad=0.4', facecolor='white',
                     edgecolor='gray', linewidth=0.8, alpha=0.9))
    
    ax1.set_aspect('equal', adjustable='box')
    ax2.set_aspect('equal', adjustable='box')
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    logger.info(f"  Saved: {output_file}")
    logger.info(f"  -log10(P) correlation: {corr_p:.4f}")
    logger.info(f"  BETA correlation: {corr_beta:.4f}")
    
    # Reset
    plt.rcParams.update(plt.rcParamsDefault)


def plot_colored_comparison(merged, dataset1_name, dataset2_name, color_column, 
                           output_file, sig_level, logger):
    """Plot comparison with coloring by a specific annotation column."""
    logger.info(f"Generating plot colored by: {color_column}")
    
    # Intelligent categorical vs continuous detection
    is_hwe = 'HWE' in color_column.upper()
    
    # Try to convert to numeric if possible
    if merged[color_column].dtype == 'object':
        try:
            numeric_values = pd.to_numeric(merged[color_column], errors='coerce')
            if numeric_values.notna().sum() / len(numeric_values) > 0.9:  # If >90% can be converted
                merged[color_column] = numeric_values
                logger.info(f"  Converted {color_column} from object to numeric")
        except:
            pass
    
    # Determine if categorical or continuous
    n_unique = merged[color_column].nunique()
    is_numeric = np.issubdtype(merged[color_column].dtype, np.number)
    
    # Categorical: object/category type OR numeric with <=20 unique values
    # Continuous: numeric type with >20 unique values
    if is_numeric and n_unique > 20:
        is_categorical = False
        logger.info(f"  Detected continuous variable (n_unique={n_unique})")
    elif not is_numeric or n_unique <= 20:
        is_categorical = True
        logger.info(f"  Detected categorical variable (n_unique={n_unique})")
    else:
        is_categorical = False
        logger.info(f"  Detected continuous variable")
    
    # Transform HWE values or prepare color values
    if is_hwe and is_numeric:
        logger.info(f"  Applying -log10 transformation for HWE values")
        color_values = -np.log10(merged[color_column].clip(lower=1e-300))
        color_label = f'-log$_{{10}}$({color_column})'
    else:
        color_values = merged[color_column]
        color_label = color_column
    
    # Set colormap and normalization
    if is_categorical:
        logger.info(f"  Detected categorical variable")
        unique_vals = color_values.unique()
        n_categories = len(unique_vals)
        logger.info(f"  Categories: {n_categories} ({', '.join(map(str, unique_vals[:10]))}...)")
        
        # Create discrete colormap
        if n_categories <= 10:
            cmap = plt.cm.tab10
        elif n_categories <= 20:
            cmap = plt.cm.tab20
        else:
            cmap = plt.cm.gist_rainbow
        
        # Map categories to integers
        category_map = {cat: i for i, cat in enumerate(unique_vals)}
        color_values_numeric = color_values.map(category_map)
        
        # Create boundaries for discrete colorbar
        bounds = np.arange(n_categories + 1) - 0.5
        norm = mcolors.BoundaryNorm(bounds, cmap.N)
    else:
        logger.info(f"  Detected continuous variable")
        cmap = 'plasma'
        norm = None
    
    # Set publication-quality parameters
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = ['DejaVu Sans', 'Bitstream Vera Sans', 'Computer Modern Sans Serif', 'Lucida Grande', 'Verdana', 'Geneva', 'Lucid', 'Arial', 'Helvetica', 'Avant Garde', 'sans-serif']
    plt.rcParams['font.size'] = 10
    plt.rcParams['axes.linewidth'] = 1
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    
    # === Left panel: -log10(P) comparison ===
    max_val_p = max(merged[f'neglog10P_{dataset1_name}'].max(),
                    merged[f'neglog10P_{dataset2_name}'].max())
    padding_p = max_val_p * 0.1
    limit_p = max_val_p + padding_p
    
    if is_categorical:
        scatter1 = ax1.scatter(merged[f'neglog10P_{dataset1_name}'], merged[f'neglog10P_{dataset2_name}'],
                   c=color_values_numeric, cmap=cmap, norm=norm,
                   alpha=0.5, s=15, edgecolors='none', rasterized=True, zorder=3)
    else:
        scatter1 = ax1.scatter(merged[f'neglog10P_{dataset1_name}'], merged[f'neglog10P_{dataset2_name}'],
                   c=color_values, cmap=cmap, norm=norm,
                   alpha=0.5, s=15, edgecolors='none', rasterized=True, zorder=3)
    
    ax1.set_xlim(-padding_p * 0.5, limit_p)
    ax1.set_ylim(-padding_p * 0.5, limit_p)
    
    # Diagonal line
    diag_max = limit_p
    ax1.plot([0, diag_max], [0, diag_max],
            color='#E63946', linestyle='--', lw=1.2, alpha=0.8, zorder=2)
    
    # Significance threshold
    sig_threshold = -np.log10(sig_level)
    ax1.axhline(y=sig_threshold, color='gray', linestyle=':', alpha=0.6, linewidth=1, zorder=1)
    ax1.axvline(x=sig_threshold, color='gray', linestyle=':', alpha=0.6, linewidth=1, zorder=1)
    
    ax1.set_xlabel(f'{dataset1_name} -log$_{{10}}$(P)', fontsize=11)
    ax1.set_ylabel(f'{dataset2_name} -log$_{{10}}$(P)', fontsize=11)
    ax1.set_title(f'P-value Comparison (Colored by {color_column})', fontsize=12, fontweight='bold', pad=10)
    ax1.grid(True, alpha=0.15, linestyle='-', linewidth=0.5, zorder=0)
    ax1.tick_params(axis='both', which='major', labelsize=10)
    
    corr_p = np.corrcoef(merged[f'neglog10P_{dataset1_name}'],
                         merged[f'neglog10P_{dataset2_name}'])[0, 1]
    ax1.text(0.98, 0.02, f'r = {corr_p:.3f}',
            transform=ax1.transAxes, fontsize=10,
            verticalalignment='bottom', horizontalalignment='right',
            bbox=dict(boxstyle='round,pad=0.4', facecolor='white',
                     edgecolor='gray', linewidth=0.8, alpha=0.9))
    
    # === Right panel: BETA comparison ===
    max_val_beta = max(abs(merged[f'BETA_{dataset1_name}']).max(),
                       abs(merged[f'BETA_{dataset2_name}']).max())
    padding_beta = max_val_beta * 0.15
    limit_beta = max_val_beta + padding_beta
    
    if is_categorical:
        scatter2 = ax2.scatter(merged[f'BETA_{dataset1_name}'], merged[f'BETA_{dataset2_name}'],
                   c=color_values_numeric, cmap=cmap, norm=norm,
                   alpha=0.45, s=12, edgecolors='none', rasterized=True, zorder=3)
    else:
        scatter2 = ax2.scatter(merged[f'BETA_{dataset1_name}'], merged[f'BETA_{dataset2_name}'],
                   c=color_values, cmap=cmap, norm=norm,
                   alpha=0.45, s=12, edgecolors='none', rasterized=True, zorder=3)
    
    ax2.set_xlim(-limit_beta, limit_beta)
    ax2.set_ylim(-limit_beta, limit_beta)
    
    # Diagonal line
    ax2.plot([-limit_beta, limit_beta], [-limit_beta, limit_beta],
            color='#E63946', linestyle='--', lw=1.2, alpha=0.8, zorder=2)
    
    # Zero lines
    ax2.axhline(y=0, color='black', linestyle='-', alpha=0.25, linewidth=0.8, zorder=1)
    ax2.axvline(x=0, color='black', linestyle='-', alpha=0.25, linewidth=0.8, zorder=1)
    
    ax2.set_xlabel(f'{dataset1_name} β', fontsize=11)
    ax2.set_ylabel(f'{dataset2_name} β', fontsize=11)
    ax2.set_title(f'Effect Size (β) Comparison (Colored by {color_column})', fontsize=12, fontweight='bold', pad=10)
    ax2.grid(True, alpha=0.15, linestyle='-', linewidth=0.5, zorder=0)
    ax2.tick_params(axis='both', which='major', labelsize=10)
    
    corr_beta = np.corrcoef(merged[f'BETA_{dataset1_name}'],
                           merged[f'BETA_{dataset2_name}'])[0, 1]
    ax2.text(0.98, 0.02, f'r = {corr_beta:.3f}',
            transform=ax2.transAxes, fontsize=10,
            verticalalignment='bottom', horizontalalignment='right',
            bbox=dict(boxstyle='round,pad=0.4', facecolor='white',
                     edgecolor='gray', linewidth=0.8, alpha=0.9))
    
    ax1.set_aspect('equal', adjustable='box')
    ax2.set_aspect('equal', adjustable='box')
    
    # Adjust for colorbar
    plt.subplots_adjust(left=0.07, right=0.86, bottom=0.12, top=0.92, wspace=0.30)
    
    # Add colorbar
    cbar_ax = fig.add_axes([0.88, 0.20, 0.018, 0.60])
    cbar = fig.colorbar(scatter2, cax=cbar_ax)
    cbar.set_label(color_label, fontsize=11.5, labelpad=15, rotation=270, va='bottom')
    cbar.ax.tick_params(labelsize=10, length=5, width=1)
    cbar.outline.set_linewidth(1)
    
    # Set colorbar ticks for categorical
    if is_categorical and n_categories <= 20:
        cbar.set_ticks(range(n_categories))
        cbar.ax.set_yticklabels(list(category_map.keys()))
    
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    logger.info(f"  Saved: {output_file}")
    
    # Reset
    plt.rcParams.update(plt.rcParamsDefault)


def plot_colored_comparison_worker(args):
    """Worker function for parallel plot generation."""
    color_col, merged_with_stats, dataset1_name, dataset2_name, output_prefix, sig_level, log_file = args
    
    # Setup logger for this worker
    logger = logging.getLogger(f"{__name__}.{color_col}")
    logger.setLevel(logging.INFO)
    handler = logging.FileHandler(log_file)
    handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
    logger.addHandler(handler)
    
    if color_col in merged_with_stats.columns:
        valid_count = merged_with_stats[color_col].notna().sum()
        logger.info(f"Processing column: {color_col}")
        logger.info(f"  Variants with valid {color_col}: {valid_count:,} / {len(merged_with_stats):,}")
        
        if valid_count > 0:
            merged_filtered = merged_with_stats[merged_with_stats[color_col].notna()].copy()
            output_file = f"{output_prefix}.{color_col}.png"
            
            try:
                plot_colored_comparison(merged_filtered, dataset1_name, dataset2_name,
                                      color_col, output_file, sig_level, logger)
                return (color_col, True, None)
            except Exception as e:
                logger.error(f"Error plotting {color_col}: {e}")
                return (color_col, False, str(e))
        else:
            logger.warning(f"Skipping {color_col}: No valid values")
            return (color_col, False, "No valid values")
    else:
        logger.warning(f"Column {color_col} not found")
        return (color_col, False, "Column not found")


def main():
    parser = argparse.ArgumentParser(
        description='Compare two GWAS summary statistics with optional variant annotation coloring'
    )
    parser.add_argument('--sumstat1', required=True, help='First summary statistics file (.glm.logistic)')
    parser.add_argument('--sumstat2', required=True, help='Second summary statistics file (.glm.logistic)')
    parser.add_argument('--name1', required=True, help='Name for first dataset')
    parser.add_argument('--name2', required=True, help='Name for second dataset')
    parser.add_argument('--variant-stats', required=True, help='Variant statistics file (compressed TSV.gz)')
    parser.add_argument('--color-columns', nargs='+', help='Columns to use for coloring (e.g., MAF_ALL IMPUTED_R2)')
    parser.add_argument('--output-prefix', required=True, help='Output file prefix')
    parser.add_argument('--sig-level', type=float, default=5e-8, help='Significance threshold (default: 5e-8)')
    parser.add_argument('--log-file', required=True, help='Log file path')
    parser.add_argument('--threads', type=int, default=None, help='Number of threads for parallel processing')
    parser.add_argument('--tabix', required=True, help='Path to tabix executable')
    
    args = parser.parse_args()
    
    # Setup logging
    logger = setup_logging(args.log_file)
    logger.info("="*70)
    logger.info("GWAS Summary Statistics Comparison")
    logger.info("="*70)
    logger.info(f"Dataset 1: {args.name1}")
    logger.info(f"Dataset 2: {args.name2}")
    logger.info(f"Significance level: {args.sig_level}")
    
    # Determine number of threads
    n_threads = args.threads if args.threads else max(1, cpu_count() - 1)
    logger.info(f"Using {n_threads} threads for parallel processing")
    
    # Read and clean data
    logger.info("\n" + "="*70)
    logger.info("STEP 1: Loading and Cleaning Data")
    logger.info("="*70)
    
    data1 = read_glm_file(args.sumstat1, logger)
    data2 = read_glm_file(args.sumstat2, logger)
    
    data1_clean = clean_sumstats(data1, args.name1, logger)
    data2_clean = clean_sumstats(data2, args.name2, logger)
    
    # Merge datasets
    logger.info("\n" + "="*70)
    logger.info("STEP 2: Merging Datasets")
    logger.info("="*70)
    
    merged = merge_datasets(data1_clean, data2_clean, args.name1, args.name2, logger)
    
    # Generate basic comparison plot
    logger.info("\n" + "="*70)
    logger.info("STEP 3: Generating Basic Comparison Plot")
    logger.info("="*70)
    
    basic_plot = f"{args.output_prefix}.basic.png"
    plot_basic_comparison(merged, args.name1, args.name2, basic_plot, args.sig_level, logger)
    
    # Generate colored plots if requested
    if args.color_columns:
        logger.info("\n" + "="*70)
        logger.info("STEP 4: Loading Variant Statistics and Generating Colored Plots")
        logger.info("="*70)
        
        # Read variant statistics using tabix
        variant_stats = read_variant_stats(args.variant_stats, set(merged['SNPID']), 
                                          args.color_columns, args.tabix, logger)
        
        # Merge with variant stats
        merged_with_stats = pd.merge(merged, variant_stats, left_on='SNPID', right_on='ID', how='left')
        logger.info(f"Merged with variant statistics: {len(merged_with_stats):,} variants")
        
        # Prepare arguments for parallel processing
        plot_args = [
            (col, merged_with_stats.copy(), args.name1, args.name2, 
             args.output_prefix, args.sig_level, f"{args.output_prefix}.{col}.log")
            for col in args.color_columns
        ]
        
        # Generate plots in parallel
        logger.info(f"Generating {len(args.color_columns)} plots in parallel with {n_threads} threads")
        
        if n_threads > 1 and len(args.color_columns) > 1:
            with Pool(processes=min(n_threads, len(args.color_columns))) as pool:
                results = pool.map(plot_colored_comparison_worker, plot_args)
        else:
            # Single-threaded execution
            results = [plot_colored_comparison_worker(arg) for arg in plot_args]
        
        # Log results
        logger.info("\nPlot generation results:")
        for col, success, error in results:
            if success:
                logger.info(f"  ✓ {col}: Success")
            else:
                logger.warning(f"  ✗ {col}: Failed - {error}")
    
    logger.info("\n" + "="*70)
    logger.info("Comparison Complete!")
    logger.info("="*70)


if __name__ == '__main__':
    main()
