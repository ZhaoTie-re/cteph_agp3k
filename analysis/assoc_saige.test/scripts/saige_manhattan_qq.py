import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
import argparse
import sys
import multiprocessing
from scipy import stats

# Try to import adjustText for label placement, but don't fail if missing
try:
    from adjustText import adjust_text
except ImportError:
    adjust_text = None

def parse_args():
    parser = argparse.ArgumentParser(description="Generate Manhattan and QQ plots for SAIGE association results")
    parser.add_argument("--input", required=True, help="Input association file (must contain CHR, POS, p.value)")
    parser.add_argument("--output-prefix", required=True, help="Prefix for output plot files")
    parser.add_argument("--title", required=False, default="SAIGE Association Results", help="Title for the plots")
    return parser.parse_args()

def calculate_lambda_gc(pvals):
    """Calculate Genomic Control Lambda (lambda GC)"""
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
        return 99

def plot_single_pvalue(df_pvals, df_chr_pos, p_col, chr_col, pos_col, title, output_filename, total_variants):
    """
    Generate and save Manhattan and QQ plots for a specific p-value column.
    
    Args:
        df_pvals (pd.Series or np.ndarray): P-values for the column.
        df_chr_pos (pd.DataFrame): DataFrame containing CHR and POS columns.
        p_col (str): Name of the p-value column.
        chr_col (str): Name of the chromosome column.
        pos_col (str): Name of the position column.
        title (str): Title for the plot.
        output_filename (str): Output filename.
        total_variants (int): Total number of variants (for expected p-value calculation).
    """
    print(f"[Process {multiprocessing.current_process().name}] Plotting for column: {p_col}")
    
    # Create a local DataFrame for plotting to keep memory usage low within the process
    # We reconstruct the necessary dataframe from the inputs.
    # Note: inputs might be views or copies, let's treat them carefully.
    
    # Ensure inputs are 1D arrays or Series
    p_values = np.array(df_pvals)
    chr_values = np.array(df_chr_pos[chr_col])
    pos_values = np.array(df_chr_pos[pos_col])
    
    # Filter Valid Data (NaN P-values are useless for plotting)
    mask = ~np.isnan(p_values)
    p_values = p_values[mask]
    chr_values = chr_values[mask]
    pos_values = pos_values[mask]

    if len(p_values) == 0:
        print(f"No valid data rows found for column {p_col}.")
        return

    # Create a lightweight DataFrame
    df = pd.DataFrame({
        'CHR': chr_values,
        'POS': pos_values,
        'P': p_values
    })
    
    # Optimizing memory: Convert types if possible
    # CHR is often string or int, POS is int. P is float.
    
    # Clip P-values to avoid log(0)
    df['P'] = df['P'].clip(lower=1e-300, upper=1.0)
    
    # -------------------------------------------------------------------------
    # Prepare Data
    # -------------------------------------------------------------------------
    
    # Map Chromosomes
    df['CHR_NUM'] = df['CHR'].apply(map_chromosome)
    df = df[df['CHR_NUM'] < 99] # Filter out weird contigs
    df = df.sort_values(by=['CHR_NUM', 'POS'])
    
    # Calculate -log10 P
    df['LOG10_P'] = -np.log10(df['P'])
    
    # Calculate GC Lambda
    lambda_val = calculate_lambda_gc(df['P'].values)
    print(f"Lambda GC ({p_col}): {lambda_val}")
    
    # -------------------------------------------------------------------------
    # Plotting Setup
    # -------------------------------------------------------------------------
    
    # Set publication style params
    plt.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans', 'Liberation Sans', 'sans-serif'],
        'font.size': 24,
        'axes.labelsize': 32,
        'axes.titlesize': 36,
        'xtick.labelsize': 24,
        'ytick.labelsize': 24,
        'figure.dpi': 300,
        'axes.linewidth': 3.0,            
        'axes.spines.top': False,
        'axes.spines.right': False,
        'legend.fontsize': 24,
        'legend.frameon': False,
        'xtick.direction': 'out',
        'ytick.direction': 'out',
        'mathtext.fontset': 'custom',
        'mathtext.rm': 'Arial',
        'mathtext.it': 'Arial:italic',
        'mathtext.bf': 'Arial:bold'
    })
    
    fig = plt.figure(figsize=(28, 10), facecolor='white') 
    
    # Margins
    left_margin = 0.08
    right_margin = 0.02
    bottom_margin = 0.15
    top_margin = 0.90
    gap = 0.08 
    
    plot_height = top_margin - bottom_margin
    qq_width = plot_height * (10 / 28) 
    man_width = 1.0 - left_margin - right_margin - gap - qq_width
    
    rect_man = [left_margin, bottom_margin, man_width, plot_height]
    rect_qq = [left_margin + man_width + gap, bottom_margin, qq_width, plot_height]
    
    ax_man = fig.add_axes(rect_man)
    ax_qq = fig.add_axes(rect_qq)

    # Global Y Limit
    max_logp_val = df['LOG10_P'].max()
    n_tests = total_variants
    bonferroni_thresh = -np.log10(0.05 / n_tests) if n_tests > 0 else 8.0
    
    global_ylim = max(max_logp_val, bonferroni_thresh) * 1.30
    global_ylim = max(global_ylim, 8.0) 

    # -------------------------------------------------------------------------
    # Manhattan Plot
    # -------------------------------------------------------------------------
    
    chromosomes = sorted(df['CHR_NUM'].unique())
    colors = ['#1F4E79', '#8DB3E2'] 

    x_labels = []
    x_ticks = []
    
    chr_offset_map = {}
    current_offset = 0
    
    for chrom in chromosomes:
        c_data = df[df['CHR_NUM'] == chrom]
        if c_data.empty: continue
        
        min_pos = c_data['POS'].min()
        max_pos = c_data['POS'].max()
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
        
    for i, chrom in enumerate(chromosomes):
        if chrom not in chr_offset_map: continue
        c_data = df[df['CHR_NUM'] == chrom]
        offset, min_p = chr_offset_map[chrom]
        x_glob = offset + (c_data['POS'] - min_p)
        
        # Optimization: Downsample non-significant points for plotting speed/file size
        # Keep all significant hits (e.g., p < 1e-3) and downsample others
        
        # Simple threshold for "interesting" points to keep all of
        keep_thresh_logp = 2.0  # p < 0.01
        
        high_sig = c_data[c_data['LOG10_P'] >= keep_thresh_logp]
        low_sig = c_data[c_data['LOG10_P'] < keep_thresh_logp]
        
        # Plot high sig fully
        ax_man.scatter(offset + (high_sig['POS'] - min_p), high_sig['LOG10_P'], 
                       color=colors[i % 2], s=40, alpha=0.9, linewidth=0, zorder=2, rasterized=True)
                       
        # Plot low sig downsampled (e.g., 10%)
        if not low_sig.empty:
            # Random sample 10%
            if len(low_sig) > 10000:
                 low_sig_sampled = low_sig.sample(frac=0.1, random_state=42)
            else:
                 low_sig_sampled = low_sig 
                 
            ax_man.scatter(offset + (low_sig_sampled['POS'] - min_p), low_sig_sampled['LOG10_P'],
                           color=colors[i % 2], s=40, alpha=0.9, linewidth=0, zorder=1, rasterized=True)
            
    gws_thresh = -np.log10(5e-8)
    ax_man.axhline(gws_thresh, color='#D32F2F', linestyle='--', linewidth=2.0, alpha=0.8, zorder=3,
                   label=r'Genome-wide ($P < 5 \times 10^{-8}$)')

    bonferroni_p = 0.05/n_tests
    bonferroni_exp = int(np.floor(np.log10(bonferroni_p)))
    bonferroni_base = bonferroni_p / (10**bonferroni_exp)
    label_bonf = r'Bonferroni ($P < %.2f \times 10^{%d}$)' % (bonferroni_base, bonferroni_exp)
    
    ax_man.axhline(bonferroni_thresh, color='#7B1FA2', linestyle=':', linewidth=2.0, alpha=0.8, zorder=3,
                   label=label_bonf)
    
    sig_hits = df[df['LOG10_P'] >= gws_thresh].copy()
    if not sig_hits.empty:
        hit_x = []
        hit_y = []
        for _, row in sig_hits.iterrows():
            c = row['CHR_NUM']
            if c not in chr_offset_map: continue
            offset, min_p = chr_offset_map[c]
            x = offset + (row['POS'] - min_p)
            hit_x.append(x)
            hit_y.append(row['LOG10_P'])
        ax_man.scatter(hit_x, hit_y, color='#CC0000', s=40, alpha=1.0, linewidth=0.5, edgecolor='black', zorder=4)
        
    ax_man.set_xticks(x_ticks)
    staggered_labels = [l if i % 2 == 0 else f"\n{l}" for i, l in enumerate(x_labels)]
    ax_man.set_xticklabels(staggered_labels, fontsize=20, fontweight='bold')
    ax_man.tick_params(axis='x', length=0, pad=12) 
    
    ax_man.set_xlim(-current_offset*0.015, current_offset*1.015)
    ax_man.set_ylim(0, global_ylim)
        
    ax_man.set_xlabel('Chromosome', fontsize=28, fontweight='bold', labelpad=18)
    ax_man.set_ylabel(r'$-\log_{10}(P)$', fontsize=28, fontweight='bold', labelpad=18)
    ax_man.set_title(title, fontweight='bold', fontsize=32, pad=25, loc='left')
    
    h_man, l_man = ax_man.get_legend_handles_labels()
    by_label_man = dict(zip(l_man, h_man))
    leg = ax_man.legend(by_label_man.values(), by_label_man.keys(), loc='upper left', 
                  frameon=True, fancybox=True, edgecolor='#BDBDBD', framealpha=0.9,
                  handlelength=2.0, handletextpad=0.8, borderaxespad=1.0,
                  fontsize=24)
    leg.get_frame().set_linewidth(0.5) 
    
    # -------------------------------------------------------------------------
    # QQ Plot
    # -------------------------------------------------------------------------
    
    ax_qq.set_ylim(0, global_ylim)
    ax_qq.set_xlim(0, global_ylim)
    ax_qq.set_aspect('equal')
    ax_qq.grid(True, linestyle='-', linewidth=0.5, color='#E0E0E0', alpha=1.0)
    
    # Using df['P'] for QQ
    p_sorted = np.sort(df['P'].values)
    observed_logp = -np.log10(p_sorted)
    
    pp = (np.arange(1, len(p_sorted) + 1) - 0.5) / len(p_sorted)
    expected_logp = -np.log10(pp)
    
    try:
        if len(p_sorted) > 10000:
            head_indices = np.arange(1, 1001)
            tail_indices = np.unique(np.geomspace(1001, len(p_sorted), num=5000).astype(int))
            sample_indices = np.concatenate([head_indices, tail_indices])
        else:
            sample_indices = np.arange(1, len(p_sorted) + 1)
            
        ci_exp = -np.log10((sample_indices - 0.5) / len(p_sorted))
        lower_p = stats.beta.ppf(0.025, sample_indices, len(p_sorted) - sample_indices + 1)
        upper_p = stats.beta.ppf(0.975, sample_indices, len(p_sorted) - sample_indices + 1)
        ci_upper_y = -np.log10(lower_p)
        ci_lower_y = -np.log10(upper_p)
        ax_qq.fill_between(ci_exp, ci_lower_y, ci_upper_y, color='#B0BEC5', alpha=0.4, zorder=1, label='95% CI')
    except Exception as e:
        print(f"Warning: Could not calculate CI: {e}")

    colors_qq = ["#B0C4DE", "#4682B4", "#1F4E79", "#000080"]
    nodes_qq = [0.0, 0.3, 0.6, 1.0]
    cmap_qq = mcolors.LinearSegmentedColormap.from_list("saige_qq", list(zip(nodes_qq, colors_qq)))
    norm_qq = plt.Normalize(vmin=0, vmax=max(8, observed_logp.max()))
    
    ax_qq.scatter(expected_logp, observed_logp, c=observed_logp, cmap=cmap_qq, norm=norm_qq, 
                  s=40, alpha=1.0, linewidth=0, zorder=2, label='Observed', rasterized=True)
    
    ax_qq.plot([0, global_ylim], [0, global_ylim], color='#D50000', linestyle='--', linewidth=2.0, zorder=3, label='Expected')
        
    stats_text = f"$\\lambda_{{GC}} = {lambda_val:.4f}$"
    ax_qq.text(0.05, 0.95, stats_text, 
               transform=ax_qq.transAxes, fontsize=24, fontweight='bold', 
               verticalalignment='top', horizontalalignment='left',
               bbox=dict(boxstyle='round,pad=0.4', facecolor='white', alpha=0.9, edgecolor='#BDBDBD', linewidth=1.0))
    
    ax_qq.set_xlabel(r'Expected $-\log_{10}(P)$', fontsize=28, fontweight='bold', labelpad=18)
    ax_qq.set_ylabel(r'Observed $-\log_{10}(P)$', fontsize=28, fontweight='bold', labelpad=18)
    ax_qq.set_title("Q-Q Plot", fontweight='bold', fontsize=30, pad=25, loc='left')
    
    plt.savefig(output_filename, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved plotting results to {output_filename}")


def main():
    args = parse_args()
    
    print(f"Reading data from {args.input}...")
    
    # Optimization 1: Read only necessary columns to save memory
    # First, read header to identify columns
    try:
        # Try reading just the header
        try:
             header_df = pd.read_csv(args.input, sep='\t', nrows=0)
        except:
             header_df = pd.read_csv(args.input, delim_whitespace=True, nrows=0)
             
        if header_df.shape[1] < 2:
            header_df = pd.read_csv(args.input, delim_whitespace=True, nrows=0)
            
        columns = header_df.columns.tolist()
        col_map = {c.lower(): c for c in columns}
        
    except Exception as e:
        print(f"Error reading file header: {e}")
        sys.exit(1)

    # Identify required coordinate columns
    chr_col_key = next((c for c in ['chr', 'chrom'] if c in col_map), None)
    pos_col_key = next((c for c in ['pos', 'position', 'start'] if c in col_map), None)
    
    if not chr_col_key or not pos_col_key:
        print(f"Error: Missing CHR or POS columns. Found: {columns}.")
        sys.exit(1)
        
    chr_col = col_map[chr_col_key]
    pos_col = col_map[pos_col_key]
    print(f"Using columns: CHR='{chr_col}', POS='{pos_col}'")

    # Target P-value columns
    cols_to_plot = []
    
    # Check for 'p.value' (SPA)
    if 'p.value' in columns:
        cols_to_plot.append(('p.value', f"{args.title} (SPA)", f"{args.output_prefix}.SPA.png"))
    elif 'p.value' in col_map: 
         cols_to_plot.append((col_map['p.value'], f"{args.title} (SPA)", f"{args.output_prefix}.SPA.png"))

    # Check for 'p.value.NA' (No SPA)
    if 'p.value.NA' in columns:
        cols_to_plot.append(('p.value.NA', f"{args.title} (No SPA)", f"{args.output_prefix}.noSPA.png"))
    elif 'p.value.na' in col_map:
        cols_to_plot.append((col_map['p.value.na'], f"{args.title} (No SPA)", f"{args.output_prefix}.noSPA.png"))
        
    if not cols_to_plot:
        print("Error: Could not find 'p.value' or 'p.value.NA' (or similar) in input file.")
        
        # Fallback
        p_col_key = next((c for c in ['pvalue', 'p_value', 'p'] if c in col_map), None)
        if p_col_key:
             p_col = col_map[p_col_key]
             print(f"Falling back to generic column '{p_col}'")
             cols_to_plot.append((p_col, args.title, f"{args.output_prefix}.png"))
        else:
            sys.exit(1)

    # Prepare usecols list
    usecols = [chr_col, pos_col] + [c[0] for c in cols_to_plot]
    
    # Read entire file with specific columns
    # Optimization: Specify dtype to save memory
    dtype_map = {chr_col: str, pos_col: np.int32} # CHR as string (to handle X/Y), POS as int
    for c in [c[0] for c in cols_to_plot]:
        dtype_map[c] = np.float32 # P-values as float32 is enough precision for plotting usually
        
    print(f"Reading columns: {usecols}")
    try:
        try:
            df = pd.read_csv(args.input, sep='\t', usecols=usecols, dtype=dtype_map)
        except:
             df = pd.read_csv(args.input, delim_whitespace=True, usecols=usecols, dtype=dtype_map)
    except Exception as e:
        print(f"Error reading file body: {e}")
        sys.exit(1)
        
    # Process Parallel
    processes = []
    total_len = len(df)
    
    # Shared coordinate data (pass slice to processes or let them copy)
    # Since multiprocessing pickles arguments, passing large dataframe copies is bad.
    # However, each process needs different P-column but Same CHR/POS.
    # We can pass specific columns as Series/Arrays to the function to avoid pickling the whole DF multiple times if it had unused cols.
    
    # We strip DF to only what's needed for each process
    
    for p_col, title, out_fn in cols_to_plot:
        p = multiprocessing.Process(
            target=plot_single_pvalue, 
            args=(df[p_col], df[[chr_col, pos_col]], p_col, chr_col, pos_col, title, out_fn, total_len)
        )
        processes.append(p)
        p.start()
        
    for p in processes:
        p.join()

if __name__ == "__main__":
    multiprocessing.set_start_method('fork', force=True) # Ensure efficient forking on Linux
    main()
