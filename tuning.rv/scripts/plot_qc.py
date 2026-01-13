import matplotlib
matplotlib.use('Agg') # Non-interactive backend
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import argparse
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import FuncFormatter
import gzip
import sys

def setup_style():
    # Set academic style
    sns.set_theme(style="ticks", context="paper", font_scale=1.4)
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['axes.titlesize'] = 16
    plt.rcParams['axes.labelsize'] = 14
    plt.rcParams['legend.fontsize'] = 12
    plt.rcParams['xtick.labelsize'] = 12
    plt.rcParams['ytick.labelsize'] = 12
    plt.rcParams['axes.grid'] = True
    plt.rcParams['grid.alpha'] = 0.3
    plt.rcParams['grid.linestyle'] = '--'

def load_data(args):
    try:
        # Load sample metrics
        df = pd.read_csv(args.sample_metrics, sep='\t')
        
        # Get variant count
        variant_count = 0
        with gzip.open(args.variant_metrics, 'rt') as f:
            for i, _ in enumerate(f):
                variant_count = i
        # Correction for header if needed, assuming variant_metrics has header
        # If line count logic matches qc_summary.py: count = total_lines - 1 (header)
        # i is index of last line (0-based) = total_lines - 1. So i is actually the count excluding header if there is 1 header line.
        # e.g. 1 line (header) -> i=0. count=0. Correct.
        
        # Load stats
        stats_df = pd.read_csv(args.qc_stats, sep='\t')
        stats = stats_df.iloc[0]
        
        return df, variant_count, stats
    except Exception as e:
        print(f"Error loading data: {e}", file=sys.stderr)
        sys.exit(1)

def get_annotation_text(stats, variant_count, sample_count, args, page_type):
    min_ac = args.min_ac
    
    # Use formatted text blocks
    lines = []
    
    # Header
    lines.append(r"$\bf{STUDY\ DESIGN}$")
    lines.append(f"Samples ($N$): {sample_count:,}")
    lines.append(f"Variants ($V$): {variant_count:,}")
    lines.append(f"MinAC Cutoff: {min_ac}")
    lines.append("")
    
    if page_type == 'depth':
        try:
            mwu_p = float(stats['MWU_P_TargetDP'])
            wass_dist = float(stats['Wasserstein_Dist_TargetDP'])
            
            lines.append(r"$\bf{DISTRIBUTION\ TESTS}$")
            lines.append(r"Mann-Whitney $U$:")
            lines.append(f"  $p = {mwu_p:.2e}$")
            lines.append("")
            lines.append(r"Wasserstein Dist ($W_1$):")
            lines.append(f"  $W_1 = {wass_dist:.4f}$")
        except:
            lines.append("Stats unavailable")
        
    elif page_type == 'bias':
        try:
            rho = float(stats['Spearman_Rho_SMinAC_MeanDP'])
            rho_p = float(stats['Spearman_P_SMinAC_MeanDP'])
            beta = float(stats['Poisson_MeanDP_Beta'])
            ci_lower = float(stats['Poisson_MeanDP_95CI_Lower'])
            ci_upper = float(stats['Poisson_MeanDP_95CI_Upper'])
            pois_p = float(stats['Poisson_MeanDP_P'])
            
            lines.append(r"$\bf{BIAS\ METRICS}$")
            lines.append("Spearman Rank Correlation:")
            lines.append(rf"  $\rho = {rho:.3f}$")
            lines.append(rf"  $p = {rho_p:.2e}$")
            lines.append("")
            lines.append("Poisson Regression Model:")
            lines.append(r"  $S_{MinAC} \sim Group + Z(D_{mean})$")
            lines.append(rf"  $\beta_{{depth}} = {beta:.3f}$ (95% CI: {ci_lower:.3f}, {ci_upper:.3f})")
            lines.append(rf"  $p = {pois_p:.2e}$")
        except:
            lines.append("Stats unavailable")
    
    return "\n".join(lines)

def plot_page1_depth(pdf, df, stats, variant_count, args):
    # Adjust width ratios for better balance
    fig, axes = plt.subplots(1, 3, figsize=(18, 7), gridspec_kw={'width_ratios': [1, 1, 0.6], 'wspace': 0.3})
    
    # Identify key columns
    target_col = 'TargetDP' 
    if 'SMinAC' in df.columns:
        plot_col = 'SMinAC'
        xlabel_text = "Sample Minor Allele Burden ($S_{MinAC}$)"
        title_text_dist = "Distribution"
        title_text_cdf = "Cumulative Distribution"
    else:
        plot_col = 'MeanDP'
        xlabel_text = "Mean Sequencing Depth ($D_{mean}$)"
        title_text_dist = "Depth Distribution"
        title_text_cdf = "Depth Cumulative Distribution"
    
    if target_col not in df.columns:
        return

    df[target_col] = df[target_col].astype(str)
    
    # Use high-contrast categorical palette
    unique_targets = sorted(df[target_col].unique())
    palette = sns.color_palette("Set2", len(unique_targets))
    
    stat_type = args.hist_stat
    y_label = "Density" if stat_type == 'density' else "Sample Count"
    
    # Plot 1: Histogram OR KDE
    if stat_type == 'density':
        sns.kdeplot(data=df, x=plot_col, hue=target_col, fill=True, ax=axes[0], 
                    palette=palette, alpha=0.2, linewidth=1.5, common_norm=False)
    else:
        sns.histplot(data=df, x=plot_col, hue=target_col, kde=False, ax=axes[0], 
                     element="step", stat="count", palette=palette, alpha=0.3, linewidth=1.5)
        # Apply comma formatting to Y axis for counts
        axes[0].yaxis.set_major_formatter(FuncFormatter(lambda x, p: format(int(x), ',')))
        
    axes[0].set_title(title_text_dist, fontweight='bold')
    axes[0].set_xlabel(xlabel_text)
    axes[0].set_ylabel(y_label)
    sns.move_legend(axes[0], "upper right")

    # Plot 2: ECDF
    sns.ecdfplot(data=df, x=plot_col, hue=target_col, ax=axes[1], 
                 palette=palette, linewidth=3, alpha=0.9)
    axes[1].set_title(title_text_cdf, fontweight='bold')
    axes[1].set_xlabel(xlabel_text)
    axes[1].set_ylabel("Cumulative Proportion")
    axes[1].grid(True, which='both', linestyle='--', linewidth=0.5, alpha=0.5)

    # Apply comma formatting to X axis if plotting SMinAC (counts)
    if plot_col == 'SMinAC':
        axes[0].xaxis.set_major_formatter(FuncFormatter(lambda x, p: format(int(x), ',')))
        axes[1].xaxis.set_major_formatter(FuncFormatter(lambda x, p: format(int(x), ',')))
    
    # Plot 3: Stats Text
    annot_text = get_annotation_text(stats, variant_count, len(df), args, 'depth')
    axes[2].axis('off')
    
    # Refined text layout to match plot content area
    axes[2].text(0.0, 0.95, annot_text, transform=axes[2].transAxes, 
                 verticalalignment='top', fontsize=12, family='serif', linespacing=1.5)
    
    # Adjust layout
    plt.subplots_adjust(top=0.9)
    
    pdf.savefig(fig, bbox_inches='tight')
    plt.close()

def plot_page2_bias(pdf, df, stats, variant_count, args):
    fig, axes = plt.subplots(1, 3, figsize=(18, 7), gridspec_kw={'width_ratios': [1, 1, 0.6], 'wspace': 0.3})
    
    x_col = 'MeanDP'
    y_col = 'SMinAC'
    target_col = 'TargetDP'
    group_col = 'Group'
    
    df[target_col] = df[target_col].astype(str)
    
    unique_targets = sorted(df[target_col].unique())
    palette_target = sns.color_palette("Set2", len(unique_targets))
    
    # Plot 1: Standard Scatter
    sns.scatterplot(data=df, x=x_col, y=y_col, hue=target_col, alpha=0.7, s=40, edgecolor='w', linewidth=0.5, 
                    ax=axes[0], palette=palette_target)
    axes[0].set_title(f"Burden vs. Depth\n(by Target Depth)", fontweight='bold')
    axes[0].set_xlabel("Mean Depth ($D_{mean}$)")
    axes[0].set_ylabel("Sample Minor Allele Burden ($S_{MinAC}$)")
    axes[0].yaxis.set_major_formatter(FuncFormatter(lambda x, p: format(int(x), ',')))
    
    # Plot 2: Biological Scatter
    if group_col in df.columns:
        # Custom palette for Case/Control (Red/Blue ideally)
        groups = sorted(df[group_col].dropna().unique())
        # Use a distinct palette if Case/Control
        if set(groups) <= {'Case', 'Control', 'PH', 'Control'}:
             palette_group = {"Case": "#d62728", "PH": "#d62728", "Control": "#1f77b4"} # Tab10 Red/Blue
             # Fallback if other keys
             palette_group = {k: palette_group.get(k, "#7f7f7f") for k in groups}
        else:
             palette_group = "tab10"
             
        sns.scatterplot(data=df, x=x_col, y=y_col, hue=group_col, alpha=0.7, s=40, edgecolor='w', linewidth=0.5,
                        ax=axes[1], palette=palette_group)
    else:
        sns.scatterplot(data=df, x=x_col, y=y_col, alpha=0.7, s=40, ax=axes[1], color="#1f77b4")
        
    axes[1].set_title(f"Burden vs. Depth\n(by Phenotype Group)", fontweight='bold')
    axes[1].set_xlabel("Mean Depth ($D_{mean}$)")
    axes[1].set_ylabel("Sample Minor Allele Burden ($S_{MinAC}$)")
    axes[1].yaxis.set_major_formatter(FuncFormatter(lambda x, p: format(int(x), ',')))
    
    # Plot 3: Stats Text
    annot_text = get_annotation_text(stats, variant_count, len(df), args, 'bias')
    axes[2].axis('off')
    
    # Refined text layout to match plot content area
    axes[2].text(0.0, 0.95, annot_text, transform=axes[2].transAxes, 
                 verticalalignment='top', fontsize=12, family='serif', linespacing=1.5)

    plt.subplots_adjust(top=0.9)

    pdf.savefig(fig, bbox_inches='tight')
    plt.close()

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--sample-metrics', required=True)
    parser.add_argument('--variant-metrics', required=True)
    parser.add_argument('--qc-stats', required=True)
    parser.add_argument('--min-ac', required=True)
    parser.add_argument('--out-pdf', required=True)
    parser.add_argument('--hist-stat', choices=['count', 'density'], default='density', help='Statistic to use for histogram (count or density)')
    args = parser.parse_args()
    
    try:
        setup_style()
        df, variant_count, stats = load_data(args)
        
        with PdfPages(args.out_pdf) as pdf:
            plot_page1_depth(pdf, df, stats, variant_count, args)
            plot_page2_bias(pdf, df, stats, variant_count, args)
            
        print(f"Generated plots: {args.out_pdf}")
    except Exception as e:
        print(f"Error generating plots: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()
