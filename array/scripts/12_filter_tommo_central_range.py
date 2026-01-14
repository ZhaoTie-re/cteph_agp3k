#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script: 12_filter_tommo_central_range.py
========================================

Description:
    1. Filter variants with TOMMO_FILTER == 'PASS'.
    2. Calculate DIFF = AAF_CTRL - TOMMO_AAF.
    3. Identify variants within the 95% Central Range of DIFF (2.5th to 97.5th percentile).
    4. Generate ID list for these variants.
    5. Plot distributions and scatter plots (4 panels).
    6. Generate a summary report.

Author: GitHub Copilot
"""

import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
from datetime import datetime
from matplotlib.ticker import FuncFormatter

# --- Matplotlib Configuration for Professional/Academic Style ---
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'DejaVu Sans', 'Helvetica'],
    'font.size': 10,
    'axes.titlesize': 11,
    'axes.labelsize': 10,
    'xtick.labelsize': 9,
    'ytick.labelsize': 9,
    'legend.fontsize': 9,
    'figure.titlesize': 12,
    'axes.linewidth': 0.8,
    'grid.linewidth': 0.4,
    'grid.color': '#CCCCCC',
    'savefig.bbox': 'tight',
    'savefig.dpi': 300
})

def parse_args():
    parser = argparse.ArgumentParser(description="Filter ToMMo variants by 95% Central Range of DIFF")
    parser.add_argument("--input_tsv", required=True, help="Input variant_qc_with_tommo.tsv.gz")
    parser.add_argument("--output_prefix", required=True, help="Output prefix for generated files")
    return parser.parse_args()

def _log(msg, log_lines):
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    formatted_msg = f"[{timestamp}] {msg}"
    print(formatted_msg, file=sys.stderr)
    log_lines.append(formatted_msg)

def main():
    args = parse_args()
    log_lines = []
    
    _log(f"Script started: {sys.argv[0]}", log_lines)
    _log(f"Input file: {args.input_tsv}", log_lines)
    _log(f"Output prefix: {args.output_prefix}", log_lines)

    # 1. Load Data
    # We need columns: ID, AAF_CTRL, TOMMO_AAF, TOMMO_FILTER
    _log("Loading data...", log_lines)
    try:
        # Use chunking if file is massive, but for plotting we usually need all data.
        # Assuming memory is sufficient for these columns.
        df = pd.read_csv(
            args.input_tsv, 
            sep='\t', 
            usecols=['ID', 'AAF_CTRL', 'TOMMO_AAF', 'TOMMO_FILTER'],
            dtype={'ID': 'string', 'TOMMO_FILTER': 'string'}
        )
    except Exception as e:
        _log(f"Error reading input file: {e}", log_lines)
        sys.exit(1)

    total_variants = len(df)
    _log(f"Total variants loaded: {total_variants:,}", log_lines)

    # 2. Filter TOMMO_FILTER == 'PASS'
    df['TOMMO_FILTER'] = df['TOMMO_FILTER'].str.upper()
    df_pass = df[df['TOMMO_FILTER'] == 'PASS'].copy()
    pass_count = len(df_pass)
    
    if total_variants > 0:
        pass_pct = (pass_count / total_variants) * 100
    else:
        pass_pct = 0.0
        
    _log(f"Variants with TOMMO_FILTER == 'PASS': {pass_count:,} ({pass_pct:.2f}%)", log_lines)

    if pass_count == 0:
        _log("No PASS variants found. Exiting.", log_lines)
        # Create empty outputs to avoid pipeline crash
        open(f"{args.output_prefix}.tommo_pass_95pct.ids.txt", 'w').close()
        open(f"{args.output_prefix}.tommo_pass_95pct.report.txt", 'w').write("\n".join(log_lines))
        sys.exit(0)

    # 3. Calculate DIFF
    # Ensure numeric
    df_pass['AAF_CTRL'] = pd.to_numeric(df_pass['AAF_CTRL'], errors='coerce')
    df_pass['TOMMO_AAF'] = pd.to_numeric(df_pass['TOMMO_AAF'], errors='coerce')
    
    # Drop NaNs in AAFs
    df_pass = df_pass.dropna(subset=['AAF_CTRL', 'TOMMO_AAF'])
    valid_pass_count = len(df_pass)
    _log(f"PASS variants with valid AAFs: {valid_pass_count:,}", log_lines)

    if valid_pass_count == 0:
        _log("No valid AAF data in PASS variants. Exiting.", log_lines)
        sys.exit(0)

    df_pass['DIFF'] = df_pass['AAF_CTRL'] - df_pass['TOMMO_AAF']

    # 4. Calculate 95% Central Range
    lower_q = df_pass['DIFF'].quantile(0.025)
    upper_q = df_pass['DIFF'].quantile(0.975)
    
    _log(f"95% Central Range for DIFF (2.5% - 97.5%): [{lower_q:.6f}, {upper_q:.6f}]", log_lines)

    # 5. Extract Variants in Range
    mask_central = (df_pass['DIFF'] >= lower_q) & (df_pass['DIFF'] <= upper_q)
    df_central = df_pass[mask_central].copy()
    central_count = len(df_central)
    central_pct = (central_count / valid_pass_count) * 100
    _log(f"Variants in 95% Central Range: {central_count:,} ({central_pct:.2f}% of valid PASS)", log_lines)

    # 6. Output ID List
    out_ids_file = f"{args.output_prefix}.tommo_pass_95pct.ids.txt"
    df_central['ID'].to_csv(out_ids_file, index=False, header=False)
    _log(f"Saved ID list to: {out_ids_file}", log_lines)

    # 7. Plotting
    _log("Generating plots...", log_lines)
    fig, axes = plt.subplots(2, 2, figsize=(12, 12))
    
    # Colors
    color_hist = '#4c72b0'  # Muted blue
    color_scatter = '#4c72b0'
    color_central = '#55a868' # Muted green
    line_color = '#c44e52'    # Muted red

    # Helper for scatter plots
    def plot_scatter(ax, data, color, title):
        # No subsampling as requested, using rasterized=True for performance in PDF
        plot_data = data
        subtitle = f"(N={len(data):,})"
            
        # Adjusted alpha=0.5 and s=4 to make outliers more visible while maintaining density visualization
        ax.scatter(plot_data['TOMMO_AAF'], plot_data['AAF_CTRL'], 
                   alpha=0.5, s=4, color=color, edgecolors='none', rasterized=True)
        
        # Diagonal line on top (zorder=10)
        ax.plot([0, 1], [0, 1], color='black', linestyle='--', linewidth=1, zorder=10)
        
        ax.set_title(f"{title}\n{subtitle}")
        ax.set_xlabel("ToMMo AAF")
        ax.set_ylabel("Control AAF")
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.set_aspect('equal', adjustable='box')
        
        # Academic style: remove top and right spines
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        
        ax.grid(True, linestyle=':', alpha=0.6)

    # Helper for histograms
    def plot_hist(ax, data, color, title, show_lines=False):
        mean_val = np.mean(data)
        std_val = np.std(data)
        
        # Use black edges for histogram bars for better definition
        ax.hist(data, bins=100, color=color, alpha=0.7, density=True, edgecolor='black', linewidth=0.3)
        
        if show_lines:
            ax.axvline(lower_q, color=line_color, linestyle='--', linewidth=1.5, label=f'2.5%: {lower_q:.4f}')
            ax.axvline(upper_q, color=line_color, linestyle='--', linewidth=1.5, label=f'97.5%: {upper_q:.4f}')
            ax.legend(loc='upper right', frameon=False)
            
        ax.set_title(f"{title}\n(Mean={mean_val:.4f}, SD={std_val:.4f})")
        ax.set_xlabel("DIFF (AAF_CTRL - TOMMO_AAF)")
        ax.set_ylabel("Density")
        
        # Academic style: remove top and right spines
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        
        ax.grid(True, linestyle=':', alpha=0.6)

    # Panel A: Distribution of DIFF (All PASS)
    plot_hist(axes[0, 0], df_pass['DIFF'], color_hist, "A. DIFF Distribution (All PASS)", show_lines=True)

    # Panel B: Scatter AAF_CTRL vs TOMMO_AAF (All PASS)
    plot_scatter(axes[0, 1], df_pass, color_scatter, "B. AAF Comparison (All PASS)")

    # Panel C: Distribution of DIFF (Central 95%)
    plot_hist(axes[1, 0], df_central['DIFF'], color_central, "C. DIFF Distribution (Central 95%)", show_lines=False)

    # Panel D: Scatter AAF_CTRL vs TOMMO_AAF (Central 95%)
    plot_scatter(axes[1, 1], df_central, color_central, "D. AAF Comparison (Central 95%)")

    plt.tight_layout()
    out_plot_file = f"{args.output_prefix}.tommo_pass_95pct.plot.pdf"
    plt.savefig(out_plot_file)
    _log(f"Saved plot to: {out_plot_file}", log_lines)

    # 8. Write Report
    out_report_file = f"{args.output_prefix}.tommo_pass_95pct.report.txt"
    with open(out_report_file, 'w') as f:
        f.write("========================================================\n")
        f.write("          ToMMo 95% Central Range Filter Report         \n")
        f.write("========================================================\n")
        f.write(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"Input File: {args.input_tsv}\n")
        f.write(f"Output Prefix: {args.output_prefix}\n")
        f.write("\n")
        f.write("--------------------------------------------------------\n")
        f.write("Filtering Statistics\n")
        f.write("--------------------------------------------------------\n")
        f.write(f"1. Total Variants Loaded:              {total_variants:>12,}\n")
        f.write(f"2. Variants with TOMMO_FILTER='PASS':  {pass_count:>12,} ({pass_pct:.2f}%)\n")
        f.write(f"3. Valid AAF Data (Non-NaN):           {valid_pass_count:>12,}\n")
        f.write("\n")
        f.write("--------------------------------------------------------\n")
        f.write("95% Central Range Calculation (DIFF = AAF_CTRL - TOMMO_AAF)\n")
        f.write("--------------------------------------------------------\n")
        f.write(f"Lower Bound (2.5th percentile):  {lower_q:.6f}\n")
        f.write(f"Upper Bound (97.5th percentile): {upper_q:.6f}\n")
        f.write("\n")
        f.write("--------------------------------------------------------\n")
        f.write("Final Filtering Results\n")
        f.write("--------------------------------------------------------\n")
        f.write(f"Retained Variants (Inside Range):      {central_count:>12,} ({central_pct:.2f}% of valid PASS)\n")
        f.write(f"Excluded Variants (Outside Range):     {valid_pass_count - central_count:>12,}\n")
        f.write("\n")
        f.write("--------------------------------------------------------\n")
        f.write("Execution Log\n")
        f.write("--------------------------------------------------------\n")
        f.write("\n".join(log_lines))
        f.write("\n")
    
    print("Done.")

if __name__ == "__main__":
    main()
