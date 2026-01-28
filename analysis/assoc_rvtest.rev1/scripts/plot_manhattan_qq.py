import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import argparse
import sys
from scipy import stats

try:
    from adjustText import adjust_text
except ImportError:
    adjust_text = None

def parse_args():
    parser = argparse.ArgumentParser(description="Generate Manhattan and QQ plots for gene-based tests")
    parser.add_argument("--input", required=True, help="Input association file (must contain CHR, START/POS, Pvalue)")
    parser.add_argument("--output-prefix", required=True, help="Prefix for output plot files")
    parser.add_argument("--title", required=False, default="Gene-based Test Results", help="Title for the plots")
    return parser.parse_args()

def calculate_lambda_gc(pvals):
    """Calculate Genomic Control Lambda (lambda GC) using gwaslab logic"""
    # Logic based on gwaslab: https://github.com/Cloufield/gwaslab/
    # Uses median P-value converted to Chi2 via ISF
    
    if len(pvals) == 0:
        return np.nan
        
    # 1. Calculate Median P-value (ignoring NaNs)
    median_p = np.nanmedian(pvals)
    
    # 2. Convert Median P to Chi-squared (df=1)
    # using inverse survival function (isf) which is more precise for small P
    obs_median_chi2 = stats.chi2.isf(median_p, df=1)
    
    # 3. Expected Median Chi-squared under null (median of chi2(1))
    exp_median_chi2 = stats.chi2.ppf(0.5, df=1)
    
    lambda_gc = obs_median_chi2 / exp_median_chi2
    return lambda_gc

def map_chromosome(chr_val):
    """Map chromosome string to integer for sorting"""
    s = str(chr_val).strip().lower().replace('chr', '')
    if s == 'x': return 23
    if s == 'y': return 24
    if s == 'xy': return 25
    if s == 'm' or s == 'mt': return 26
    
    if s.isdigit():
        return int(s)
    else:
        # Unknown/Unmapped contigs -> End
        return 99

def main():
    args = parse_args()
    
    print(f"Reading data from {args.input}...")
    try:
        # Try tab first as per previous script output
        df = pd.read_csv(args.input, sep='\t')
        if df.shape[1] < 2:
            df = pd.read_csv(args.input, delim_whitespace=True)
    except Exception as e:
        print(f"Error reading file: {e}")
        sys.exit(1)

    # Normalize columns
    df.columns = [c.lower() for c in df.columns]
    
    # Identify required columns
    # Previous step added 'chr', 'start'
    p_col = next((c for c in ['pvalue', 'p.value', 'p_value', 'p'] if c in df.columns), None)
    chr_col = next((c for c in ['chr', 'chrom'] if c in df.columns), None)
    pos_col = next((c for c in ['start', 'pos', 'position'] if c in df.columns), None)
    gene_col = next((c for c in ['gene', 'genename', 'symbol'] if c in df.columns), None)
    
    if not p_col or not chr_col or not pos_col:
        print(f"Error: Missing columns. Found: {df.columns.tolist()}. Need P, CHR, START.")
        sys.exit(1)

    # Filter Valid Data
    df = df.dropna(subset=[p_col, chr_col, pos_col]).copy()
    # Filter P-values out of valid range (0, 1]
    # Small epsilon for P=0 to avoid -inf log
    df[p_col] = df[p_col].clip(lower=1e-300, upper=1.0)
    
    if df.empty:
        print("No valid data rows found for plotting.")
        sys.exit(0)

    # -------------------------------------------------------------------------
    # Prepare Data for Plotting
    # -------------------------------------------------------------------------
    
    # Map Chromosomes
    df['CHR_NUM'] = df[chr_col].apply(map_chromosome)
    df = df[df['CHR_NUM'] < 99] # Filter out weird contigs
    df = df.sort_values(by=['CHR_NUM', pos_col])
    
    # Calculate -log10 P
    df['LOG10_P'] = -np.log10(df[p_col])
    
    # Calculate GC Lambda
    lambda_val = calculate_lambda_gc(df[p_col].values)
    
    # -------------------------------------------------------------------------
    # Plotting Setup
    # -------------------------------------------------------------------------
    
    # -------------------------------------------------------------------------
    # Plotting Setup (Publication Quality)
    # -------------------------------------------------------------------------
    
    # Set publication style params manually to ensure consistency
    # Using 'Arial' or 'Helvetica' style logic via sans-serif
    plt.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans', 'Liberation Sans', 'sans-serif'],
        'font.size': 18,              # Even larger base font
        'axes.labelsize': 24,         # Even larger label font
        'axes.titlesize': 28,         # Even larger title font
        'xtick.labelsize': 18,        # Larger tick font
        'ytick.labelsize': 18,
        'figure.dpi': 400,            
        'axes.linewidth': 2.5,        # Thicker axes
        'axes.spines.top': False,
        'axes.spines.right': False,
        'legend.fontsize': 27,        # Increased 1.5x (was 18)
        'legend.frameon': False,
        'xtick.direction': 'out',
        'ytick.direction': 'out',
        'xtick.major.size': 10,
        'ytick.major.size': 10
    })
    
    # Wide figure, balanced aspect ratio
    # Increased width to (30, 10) and adjusted width_ratios to [1.8, 1] to ensure QQ plot is height-constrained (fills full height) 
    # rather than width-constrained, guaranteeing Y-axis alignment with Manhattan plot.
    fig = plt.figure(figsize=(30, 10), facecolor='white')
    # GridSpec: Manhattan gets more space (1.8 : 1)
    gs = fig.add_gridspec(1, 2, width_ratios=[1.8, 1], wspace=0.15)
    
    ax_man = fig.add_subplot(gs[0])
    ax_qq = fig.add_subplot(gs[1])
    
    # -------------------------------------------------------------------------
    # Calculation of Global Y Limit (Shared)
    # -------------------------------------------------------------------------
    
    # We want Y axes to align. Determine the max Y needed.
    max_logp_val = df['LOG10_P'].max()
    
    # Calculate Bonferroni
    n_tests = len(df)
    bonferroni_thresh = -np.log10(0.05 / n_tests)
    
    # Determine Ceiling
    # If explicit significant hits exist, go higher. If not, at least show threshold.
    global_ylim = max(max_logp_val, bonferroni_thresh) * 1.15
    # Ensure a minimum height for visual aesthetics (e.g. 8)
    global_ylim = max(global_ylim, 8.0)

    # -------------------------------------------------------------------------
    # Manhattan Plot (Left)
    # -------------------------------------------------------------------------
    
    chromosomes = sorted(df['CHR_NUM'].unique())
    colors = ['#4D4D4D', '#A6A6A6'] # Professional Grey Scale

    x_labels = []
    x_ticks = []
    
    # Add horizontal grid for readability (Subtle)
    ax_man.grid(axis='y', linestyle='-', linewidth=0.5, color='#E0E0E0', alpha=1.0, zorder=0)

    # Pre-calculate global offsets
    chr_offset_map = {}
    current_offset = 0
    
    for chrom in chromosomes:
        c_data = df[df['CHR_NUM'] == chrom]
        if c_data.empty: continue
        
        min_pos = c_data[pos_col].min()
        max_pos = c_data[pos_col].max()
        c_len = max_pos - min_pos
        
        chr_offset_map[chrom] = (current_offset, min_pos)
        
        mid_pt = current_offset + (c_len / 2)
        x_ticks.append(mid_pt)
        
        label = str(chrom)
        if chrom == 23: label = 'X'
        elif chrom == 24: label = 'Y'
        elif chrom == 25: label = 'XY'
        elif chrom == 26: label = 'MT'
        x_labels.append(label)
        
        current_offset += c_len + 1 # Buffer
        
    last_x = current_offset 

    # Plot Background Points (All) - zorder=1 to stay behind threshold lines
    for i, chrom in enumerate(chromosomes):
        if chrom not in chr_offset_map: continue
        
        c_data = df[df['CHR_NUM'] == chrom]
        offset, min_p = chr_offset_map[chrom]
        
        x_glob = offset + (c_data[pos_col] - min_p)
        
        # Consistent small dots
        ax_man.scatter(x_glob, c_data['LOG10_P'], 
                       color=colors[i % 2], s=27, alpha=1.0, linewidth=0, zorder=2)
            
    # Draw Bonferroni Threshold
    # Using a dark line for threshold
    ax_man.axhline(bonferroni_thresh, color='#CC0000', linestyle='--', linewidth=1.5, alpha=1.0, zorder=3,
                   label=r'Bonferroni ($P < %.1e$)' % (0.05/n_tests))
    
    # Highlight Significant Hits (Bonferroni) - RED
    sig_bonf_hits = df[df['LOG10_P'] >= bonferroni_thresh].copy()
    
    if not sig_bonf_hits.empty:
        hit_x = []
        hit_y = []
        for _, row in sig_bonf_hits.iterrows():
            c = row['CHR_NUM']
            if c not in chr_offset_map: continue
            
            offset, min_p = chr_offset_map[c]
            x = offset + (row[pos_col] - min_p)
            
            hit_x.append(x)
            hit_y.append(row['LOG10_P'])
            
        # Plot highlight points - Distinct Red
        ax_man.scatter(hit_x, hit_y, color='#CC0000', s=83, alpha=1.0, linewidth=0.5, edgecolor='black', zorder=4, 
                       label='Significant')

        # Annotate Significant Genes
        if gene_col:
            texts = []
            for _, row in sig_bonf_hits.iterrows():
                c = row['CHR_NUM']
                if c not in chr_offset_map: continue
                
                offset, min_p = chr_offset_map[c]
                x = offset + (row[pos_col] - min_p)
                y = row['LOG10_P']
                label = str(row[gene_col])
                
                # Add text
                t = ax_man.text(x, y, label, 
                                fontstyle='italic', fontsize=24, fontweight='bold', # Increased 1.5x (was 16)
                                ha='center', va='bottom', zorder=10)
                texts.append(t)
            
            # Use adjust_text if installed to prevent overlap
            if adjust_text:
                adjust_text(texts, ax=ax_man, 
                            arrowprops=dict(arrowstyle="-", color='black', lw=0.5, alpha=0.8),
                            expand_points=(1.5, 1.5))
    
    # Formatting Axes
    ax_man.set_xticks(x_ticks)
    # Stagger labels
    staggered_labels = [l if i % 2 == 0 else f"\n{l}" for i, l in enumerate(x_labels)]
    ax_man.set_xticklabels(staggered_labels, fontsize=18)
    # Limit x-axis range
    ax_man.set_xlim(-current_offset*0.015, current_offset*1.015)
    
    # Apply Shared Y-Limit
    ax_man.set_ylim(0, global_ylim)
        
    ax_man.set_xlabel('Chromosome', fontsize=24, fontweight='bold', labelpad=14)
    ax_man.set_ylabel(r'$-\log_{10}(P)$', fontsize=24, fontweight='bold', labelpad=14)
    
    # Title Processing
    title_text = args.title
    if "RVTest:" in args.title:
         if "burden" in args.title.lower():
             test_type = "Burden Test (CMC)"
         elif "skato" in args.title.lower():
             test_type = "SKAT-O Test"
         else:
             test_type = "Rare Variant Association"
         title_text = test_type

    ax_man.set_title(title_text, fontweight='bold', fontsize=28, pad=24)
    
    # Legend (Manhattan)
    # Filter duplicates just in case
    h_man, l_man = ax_man.get_legend_handles_labels()
    by_label_man = dict(zip(l_man, h_man))
    ax_man.legend(by_label_man.values(), by_label_man.keys(), loc='upper right', 
                  frameon=True, fancybox=False, edgecolor='black', fontsize=27, borderpad=0.8) # Increased 1.5x (was 18)


    # -------------------------------------------------------------------------
    # QQ Plot (Right) - Professional Square & Aligned
    # -------------------------------------------------------------------------
    
    # Force Y-axis to match Manhattan
    ax_qq.set_ylim(0, global_ylim)
    
    # Force X-axis to match Y-axis (to keep range consistent as requested)
    ax_qq.set_xlim(0, global_ylim)
    
    # Ensure square aspect ratio
    # adjustable='box' changes the physical box dimensions (shrinking width since we allocated extra)
    # anchor='W' keeps it left-aligned to minimize gap with Manhattan plot
    ax_qq.set_aspect('equal', adjustable='box', anchor='W') 
    
    # Grid
    ax_qq.grid(True, linestyle='-', linewidth=0.5, color='#E0E0E0', alpha=1.0)
    
    p_sorted = np.sort(df[p_col].values)
    observed_logp = -np.log10(p_sorted)
    
    n_points = len(df)
    pp = (np.arange(1, n_points + 1) - 0.5) / n_points
    expected_logp = -np.log10(pp)
    
    # Scatter points - Dark Blue Grey
    ax_qq.scatter(expected_logp, observed_logp, c='#2C3E50', s=30, alpha=0.8, linewidth=0, zorder=2, label='Observed')
    
    # Identity Line - Red dashed
    # Line goes from 0 to global_ylim
    ax_qq.plot([0, global_ylim], [0, global_ylim], color='#CC0000', linestyle='--', linewidth=1.5, zorder=3, label='Expected')
    
    # Confidence Interval
    index = np.arange(1, n_points + 1)
    lower_p = stats.beta.ppf(0.025, index, n_points - index + 1)
    upper_p = stats.beta.ppf(0.975, index, n_points - index + 1)
    lower_log = -np.log10(upper_p)
    upper_log = -np.log10(lower_p)
    
    ax_qq.fill_between(expected_logp, lower_log, upper_log, color='#B0BEC5', alpha=0.4, zorder=1, label='95% CI')
    
    # Text Box for Lambda and N (moved to Bottom Right) - functioning as the main legend info
    stats_text = f"$\lambda_{{GC}} = {lambda_val:.3f}$\n$N_{{genes}} = {n_points:,}$"
    
    ax_qq.text(0.95, 0.05, stats_text, 
               transform=ax_qq.transAxes, fontsize=27, # Increased 1.5x (was 18)
               verticalalignment='bottom', horizontalalignment='right',
               bbox=dict(boxstyle='square,pad=0.5', facecolor='white', alpha=1.0, edgecolor='black', linewidth=1.5))

    ax_qq.set_xlabel(r'Expected $-\log_{10}(P)$', fontsize=24, fontweight='bold', labelpad=14)
    # Hide Y label if it's redundant? No, keep it for clarity.
    ax_qq.set_ylabel(r'Observed $-\log_{10}(P)$', fontsize=24, fontweight='bold', labelpad=14)
    ax_qq.set_title("Q-Q Plot", fontweight='bold', fontsize=28, pad=24)
    
    # Standard legend removed per user request: "Only keep one lambda and gene number legend"

    plt.tight_layout()
    # Use bbox_inches='tight' to ensure large labels (like "Chromosome") are not cut off
    plt.savefig(f"{args.output_prefix}.png", dpi=300, bbox_inches='tight')
    plt.savefig(f"{args.output_prefix}.pdf", dpi=300, bbox_inches='tight')
    print(f"Saved plots to {args.output_prefix}.png/pdf")

if __name__ == "__main__":
    main()
