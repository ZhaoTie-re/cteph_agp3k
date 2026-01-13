import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import argparse
import sys
from matplotlib.ticker import FuncFormatter, MaxNLocator

def setup_style():
    # Use built-in matplotlib style closely matching seaborn-ticks
    # 'seaborn-v0_8-ticks' is available in newer matplotlib, or manually set
    plt.style.use('fast') 
    
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.size'] = 11
    plt.rcParams['axes.titlesize'] = 14
    plt.rcParams['axes.labelsize'] = 12
    plt.rcParams['legend.fontsize'] = 10
    plt.rcParams['xtick.labelsize'] = 10
    plt.rcParams['ytick.labelsize'] = 10
    plt.rcParams['axes.grid'] = True
    plt.rcParams['grid.alpha'] = 0.4
    plt.rcParams['grid.linestyle'] = '--'
    
    # Tick direction out (like ticks style)
    plt.rcParams['xtick.direction'] = 'out'
    plt.rcParams['ytick.direction'] = 'out'
    
    # Remove top and right spines globally (optional, or do it per axis)
    plt.rcParams['axes.spines.top'] = False
    plt.rcParams['axes.spines.right'] = False

def format_func(x, p):
    return format(int(x), ',')

def plot_trends(args):
    try:
        df = pd.read_csv(args.qc_stats, sep='\t')
        
        # 1. Filter MinAC=0 unless requested
        if not args.include_zero:
            df = df[df['MinAC_Threshold'] != 0]

        df = df.sort_values('MinAC_Threshold') # Sort just in case

        if df.empty:
            print("No data to plot (check --include-zero or input file).", file=sys.stderr)
            # Create a dummy pdf to prevent pipeline failure
            fig, ax = plt.subplots(figsize=(4, 4))
            ax.text(0.5, 0.5, "No data available", ha='center', va='center')
            plt.savefig(args.out_pdf)
            return
        
        # Determine significance logic
        df['Spearman_Sig'] = df['Spearman_P_SMinAC_MeanDP'].fillna(1.0) < 0.05
        
        if args.sample_n:
            sample_n = args.sample_n
            threshold = 0.05 / sample_n
            df['Poisson_Sig'] = df['Poisson_MeanDP_P'].fillna(1.0) < threshold
        elif 'Sample_Count' in df.columns:
            sample_n = int(df['Sample_Count'].iloc[0])
            threshold = 0.05 / sample_n
            df['Poisson_Sig'] = df['Poisson_MeanDP_P'].fillna(1.0) < threshold
        else:
            sample_n = "N"
            threshold = 0.05
            print("Warning: Sample_Count column missing and --sample-n not provided. Using 0.05 threshold.", file=sys.stderr)
            df['Poisson_Sig'] = df['Poisson_MeanDP_P'].fillna(1.0) < threshold

        # Setup Figure
        fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(9, 14), sharex=True)
        
        x = df['MinAC_Threshold']
        
        # --- Plot 1: Variant Count ---
        ax1.plot(x, df['Variant_Count'], color='#333333', linewidth=2, zorder=1)
        ax1.scatter(x, df['Variant_Count'], color='#333333', s=40, zorder=2)
        
        ax1.set_ylabel("Variant Count")
        # Academic style usually puts labels A, B, C or no title, but we keep titles for clarity here
        ax1.set_title("A. Variant Count vs. MinAC Threshold", loc='left', fontweight='bold')
        ax1.yaxis.set_major_formatter(FuncFormatter(format_func))
        ax1.grid(True, linestyle=':', alpha=0.6)

        # --- Plot 2: Spearman Rho ---
        ax2.plot(x, df['Spearman_Rho_SMinAC_MeanDP'], color='#1f77b4', linewidth=1.5, alpha=0.7, zorder=1)
        
        sig_mask = df['Spearman_Sig']
        non_sig_mask = ~df['Spearman_Sig']
        
        # Significant (Hollow)
        ax2.scatter(x[sig_mask], df.loc[sig_mask, 'Spearman_Rho_SMinAC_MeanDP'], 
                    facecolors='white', edgecolors='#1f77b4', s=70, label=r'$P < 0.05$', zorder=2, linewidth=1.5)
        # Not Significant (Solid)
        ax2.scatter(x[non_sig_mask], df.loc[non_sig_mask, 'Spearman_Rho_SMinAC_MeanDP'], 
                    color='#1f77b4', s=70, label=r'$P \geq 0.05$', zorder=2)
            
        ax2.set_ylabel(r"Spearman's $\rho$ ($S_{MinAC}$ vs $D_{mean}$)")
        ax2.set_title("B. Correlation Bias vs. MinAC Threshold", loc='left', fontweight='bold')
        ax2.axhline(0, color='gray', linestyle='-', linewidth=1.0)
        # Place legend outside to avoid occlusion
        ax2.legend(loc='upper left', bbox_to_anchor=(1.01, 1), borderaxespad=0, frameon=True, edgecolor='gray')
        ax2.grid(True, linestyle=':', alpha=0.6)

        # --- Plot 3: Poisson Beta ---
        # CI Band
        ax3.fill_between(x, df['Poisson_MeanDP_95CI_Lower'], df['Poisson_MeanDP_95CI_Upper'], 
                         color='#d62728', alpha=0.15, label='95% CI')
        ax3.plot(x, df['Poisson_MeanDP_Beta'], color='#d62728', linewidth=1.5, alpha=0.7, zorder=1)
        
        # Labels for legend
        if isinstance(sample_n, int):
            thresh_sci = f"{threshold:.1e}" # e.g. 5.0e-05
            thresh_text = f"0.05 / {sample_n} ≈ {thresh_sci}"
            p_sig_label = fr'$P < {thresh_sci}$'
            p_nonsig_label = fr'$P \geq {thresh_sci}$'
        else:
            thresh_text = "0.05 / N"
            p_sig_label = r'$P < 0.05/N$'
            p_nonsig_label = r'$P \geq 0.05/N$'

        pois_sig_mask = df['Poisson_Sig']
        pois_non_sig_mask = ~df['Poisson_Sig']
        
        ax3.scatter(x[pois_sig_mask], df.loc[pois_sig_mask, 'Poisson_MeanDP_Beta'], 
                    facecolors='white', edgecolors='#d62728', s=70, label=p_sig_label, zorder=2, linewidth=1.5)
        ax3.scatter(x[pois_non_sig_mask], df.loc[pois_non_sig_mask, 'Poisson_MeanDP_Beta'], 
                    color='#d62728', s=70, label=p_nonsig_label, zorder=2)
        
        ax3.set_ylabel(r"Poisson $\beta_{depth}$ ($S_{MinAC} \sim D_{mean}$)")
        ax3.set_title("C. Depth Effect Size vs. MinAC Threshold", loc='left', fontweight='bold')
        ax3.set_xlabel("MinAC Threshold")
        ax3.axhline(0, color='gray', linestyle='-', linewidth=1.0)
        
        # Add N annotation
        ax3.text(0.02, 0.05, f"Significance Threshold:\n{thresh_text}", 
                 transform=ax3.transAxes, fontsize=10, 
                 bbox=dict(facecolor='white', alpha=0.9, edgecolor='lightgray', boxstyle='round,pad=0.5'))

        # Only show CI and P-value points in legend
        handles, labels = ax3.get_legend_handles_labels()
        # Place legend outside to avoid occlusion
        ax3.legend(loc='upper left', bbox_to_anchor=(1.01, 1), borderaxespad=0, frameon=True, edgecolor='gray')
        ax3.grid(True, linestyle=':', alpha=0.6)
        
        # 3. Force Integer Ticks and show all points
        # Set ticks exactly at the MinAC threshold values present in data
        all_x_values = df['MinAC_Threshold'].unique()
        ax3.set_xticks(all_x_values)
        ax3.set_xticklabels([str(int(val)) for val in all_x_values])
        
        # Despine handled by rcParams or manual
        # for ax in [ax1, ax2, ax3]:
        #    sns.despine(ax=ax) # Removed seaborn dependency

        plt.tight_layout()
        plt.savefig(args.out_pdf)
        print(f"Generated trend plots: {args.out_pdf}")
        
    except Exception as e:
        print(f"Error plotting trends: {e}", file=sys.stderr)
        sys.exit(1)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--qc-stats', required=True)
    parser.add_argument('--out-pdf', required=True)
    parser.add_argument('--include-zero', action='store_true', help='Include MinAC=0 in plots')
    parser.add_argument('--sample-n', type=int, help='Explicit sample size for significance correction')
    args = parser.parse_args()
    
    setup_style()
    plot_trends(args)

if __name__ == "__main__":
    main()
