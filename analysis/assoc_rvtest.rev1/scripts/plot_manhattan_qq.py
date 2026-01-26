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
    
    # Set publication style params manually to ensure consistency
    plt.rcParams.update({
        'font.family': 'sans-serif',
        'font.size': 10,
        'axes.labelsize': 12,
        'axes.titlesize': 14,
        'xtick.labelsize': 10,
        'ytick.labelsize': 10,
        'figure.dpi': 300,
        'axes.spines.top': False,
        'axes.spines.right': False
    })
    
    fig = plt.figure(figsize=(15, 6))
    gs = fig.add_gridspec(1, 2, width_ratios=[2, 1], wspace=0.2)
    
    ax_man = fig.add_subplot(gs[0])
    ax_qq = fig.add_subplot(gs[1])
    
    # -------------------------------------------------------------------------
    # Manhattan Plot (Left)
    # -------------------------------------------------------------------------
    
    chromosomes = sorted(df['CHR_NUM'].unique())
    colors = ['#4A4A4A', '#808080'] # Dark Grey / Light Grey classic academic
    # Alternatively: Blue/Navy ['#1f77b4', '#aec7e8'] or similar
    colors = ['#2C3E50', '#7F8C8D'] # Slate / Concrete

    x_labels = []
    x_ticks = []
    
    last_x = 0
    for i, chrom in enumerate(chromosomes):
        c_data = df[df['CHR_NUM'] == chrom]
        
        # Relative position
        # Gene based: use simple index or physical position?
        # Physical position is better for correct spacing
        r_pos = c_data[pos_col].values
        if len(r_pos) > 0:
            min_pos = r_pos.min()
            rel_pos = r_pos - min_pos
            
            # Scatter
            ax_man.scatter(last_x + rel_pos, c_data['LOG10_P'], 
                           color=colors[i % 2], s=15, alpha=0.9, linewidth=0)
            
            # Ticks
            mid_pt = last_x + (rel_pos.max() / 2)
            x_ticks.append(mid_pt)
            
            # Label
            label = str(chrom)
            if chrom == 23: label = 'X'
            elif chrom == 24: label = 'Y'
            elif chrom == 25: label = 'XY'
            elif chrom == 26: label = 'MT'
            x_labels.append(label)
            
            # Update last_x
            # Add a buffer between chromosomes equal to e.g. 5% of this chr length or fixed
            # Just simple concatenation of max pos
            last_x += rel_pos.max()
            
            # Add small gap strictly in index space?
            # For linear genome view, we usually just concat. 
            # To distinguish visually, we add a small offset
            last_x += 1 # small numeric offset if using index, usually big for pos
            # Actually, using index for plotting might be safer if positions are sparse gene centers
            # But let's stick to 'concatenated physical coordinates' approximation:
            # Shift by largest length? No, usually just add the max.
            
    # Draw Thresholds
    # Bonferroni
    n_tests = len(df)
    bonferroni_thresh = -np.log10(0.05 / n_tests)
    ax_man.axhline(bonferroni_thresh, color='#E74C3C', linestyle='--', linewidth=1.2, alpha=0.9,
                   label=f'Bonferroni (P=${0.05/n_tests:.1e}$)')
    
    # FDR (if available in file, usually q-value < 0.05)
    # Highlight significant hits with color instead of line/circle
    fdr_col = next((c for c in ['fdr', 'qval', 'q_value'] if c in df.columns), None)
    sig_fdr_hits = pd.DataFrame() # Define for later annotation logic
    
    if fdr_col:
        # Select hits that pass FDR but fail Bonferroni (to give them distinct color)
        # Hits passing Bonferroni are already obviously high
        # Or just highlight ALL FDR < 0.05 hits
        sig_fdr_hits = df[df[fdr_col] < 0.05]
        
        if not sig_fdr_hits.empty:
            # We want to re-plot these points with a distinct color/style
            # Need to recalculate their coordinates
            
            # Reconstruct chr_offset_map logic
            chr_offset_map = {}
            temp_offset = 0
            for c in chromosomes:
                 c_max = df[df['CHR_NUM'] == c][pos_col].max()
                 c_min = df[df['CHR_NUM'] == c][pos_col].min()
                 c_len = c_max - c_min
                 chr_offset_map[c] = temp_offset
                 temp_offset += c_len + 1 
            
            fdr_x = []
            fdr_y = []
            for idx, row in sig_fdr_hits.iterrows():
                g_chrom = row['CHR_NUM']
                if g_chrom not in chr_offset_map: continue
                g_rpos = row[pos_col]
                c_min = df[df['CHR_NUM'] == g_chrom][pos_col].min()
                x = chr_offset_map[g_chrom] + (g_rpos - c_min)
                fdr_x.append(x)
                fdr_y.append(row['LOG10_P'])
            
            # Plot highlight points
            ax_man.scatter(fdr_x, fdr_y, color='#E74C3C', s=25, alpha=1.0, zorder=3, 
                           label='FDR < 0.05')

    ax_man.set_xticks(x_ticks)
    # Stagger labels to show all chromosomes without overlap
    # Indices 0, 2, 4... (Chr 1, 3, 5...) on top row
    # Indices 1, 3, 5... (Chr 2, 4, 6...) on bottom row (prefixed with newline)
    staggered_labels = [l if i % 2 == 0 else f"\n{l}" for i, l in enumerate(x_labels)]
    ax_man.set_xticklabels(staggered_labels, fontsize=9, rotation=0)
        
    ax_man.set_xlabel('Chromosome', fontsize=12, fontweight='bold')
    ax_man.set_xlabel('Chromosome', fontsize=12, fontweight='bold')
    ax_man.set_ylabel(r'$-\log_{10}(P)$', fontsize=12, fontweight='bold')
    # Format Title
    # Parse title like "RVTest: impact_moderate_high - burden" to format nicely
    # Expected: "RVTest: <filter> - <method>"
    # If title structure matches, we can make it prettier.
    # Otherwise use as is.
    
    clean_title = args.title
    if "RVTest:" in args.title:
        try:
            # Example: RVTest: impact_moderate_high - burden
            parts = args.title.split(':')[-1].strip().split('-')
            if len(parts) >= 2:
                filt = parts[0].strip().replace('impact_', '')
                meth = parts[1].strip()
                # "Burden Test (Moderate & High Impact)"
                filt_clean = filt.replace('_', ' ').title()
                meth_clean = meth.title()
                if meth_clean.lower() == 'skato': meth_clean = 'SKAT-O'
                clean_title = f"{meth_clean} Test ({filt_clean} Impact)"
        except:
            pass

    # Increase title padding to make room for legend
    ax_man.set_title(clean_title, fontweight='bold', fontsize=13, pad=30)
    
    # Legend - placed above the plot area but below the title
    # Using lower center at y=1.0 puts it just above the top spine
    ax_man.legend(loc='lower center', bbox_to_anchor=(0.5, 1.0), ncol=2, 
                  borderaxespad=0, frameon=False, fontsize=9)

    # -------------------------------------------------------------------------
    # Annotate Top Genes
    # -------------------------------------------------------------------------
    if gene_col:
        # Priority 1: Bonferroni
        # Priority 2: FDR
        # Logic: Union of both sets
        
        sig_bonf = df[df['LOG10_P'] >= bonferroni_thresh]
        
        # Combine valid hits
        anno_frames = [sig_bonf]
        if fdr_col and not sig_fdr_hits.empty:
            anno_frames.append(sig_fdr_hits)
            
        if anno_frames:
            # Concat and drop duplicates
            anno_points = pd.concat(anno_frames).drop_duplicates()
        else:
            anno_points = pd.DataFrame()
        
        # If we have points to annotate
        if len(anno_points) > 0:
             texts = []
             # If too many points, limit to top 20
             if len(anno_points) > 20:
                 anno_points = anno_points.nlargest(20, 'LOG10_P')
                 
             # Re-map offsets if not already done (in case FDR highlight block was skipped)
             if 'chr_offset_map' not in locals():
                chr_offset_map = {}
                temp_offset = 0
                for c in chromosomes:
                     c_max = df[df['CHR_NUM'] == c][pos_col].max()
                     c_min = df[df['CHR_NUM'] == c][pos_col].min()
                     c_len = c_max - c_min
                     chr_offset_map[c] = temp_offset
                     temp_offset += c_len + 1

             for idx, row in anno_points.iterrows():
                 g_chrom = row['CHR_NUM']
                 if g_chrom not in chr_offset_map: continue

                 g_rpos = row[pos_col]
                 c_min = df[df['CHR_NUM'] == g_chrom][pos_col].min()
                 
                 x_coord = chr_offset_map[g_chrom] + (g_rpos - c_min)
                 y_coord = row['LOG10_P']
                 
                 gene_name = row[gene_col]
                 
                 # Italic gene name
                 t = ax_man.text(x_coord, y_coord, gene_name, fontsize=8.5, fontweight='bold', style='italic', color='#444444', 
                                 bbox=dict(facecolor='white', alpha=0.7, edgecolor='none', pad=1))
                 texts.append(t)
            
             if adjust_text:
                 # Improve repel settings to avoid overlap more aggressively
                 adjust_text(texts, ax=ax_man, 
                             arrowprops=dict(arrowstyle='-', color='#666666', lw=0.6),
                             force_points=0.3, force_text=1.0, 
                             expand_points=(1.2, 1.2), expand_text=(1.2, 1.2))
    
    # Remove logic that hides every 2nd label to ensure all chromosomes are shown
    # if len(x_labels) > 15: ... removed

    # Add numeric formatting to y-axis (though log scale usually doesn't need commas, but just in case)
    # ax_man.get_yaxis().set_major_formatter(plt.FuncFormatter(lambda x, loc: "{:,}".format(int(x))))
    
    # Clean spines
    ax_man.spines['top'].set_visible(False)
    ax_man.spines['right'].set_visible(False)
    ax_man.spines['left'].set_linewidth(0.8)
    ax_man.spines['bottom'].set_linewidth(0.8)


    # -------------------------------------------------------------------------
    # QQ Plot (Right)
    # -------------------------------------------------------------------------

    
    n_points = len(df)
    observed = np.sort(df['LOG10_P'].values)
    expected = -np.log10(np.arange(1, n_points + 1) / (n_points + 1))
    
    # Sort expected descending to match observed (which is sorted ascending? No, sort produces ascending)
    # -log10(P) large means small P.
    # small P -> large -log10.
    # sort(P) -> small to large.
    # -log(sort(P)) -> large to small.
    # But usually we plot expected vs observed 0 to max.
    # Let's match:
    # Expected: uniform 0..1 distributed -> sort -> k/(n+1)
    # -log10(k/(n+1)) -> decreases as k increases.
    # So we sort Observed ascending (small values first? No QQ usually is Large vs Large at the top right)
    
    # Standard way:
    # Expected: -log10( (N - i + 0.5) / N ) or similar
    # Observed: -log10( sort(P) )
    
    # Let's re-sort P low to high
    p_sorted = np.sort(df[p_col].values)
    observed_logp = -np.log10(p_sorted)
    
    # Expected P low to high: 1/(n+1), 2/(n+1) ... 
    # But for log scale:
    # i=1 (smallest P) -> 1/(n+1) -> largest log
    # i=n (largest P) -> n/(n+1) ~ 1 -> 0 log
    
    pp_expected = (np.arange(1, n_points + 1) - 0.5) / n_points # Hazen plotting position or simply i/(n+1)
    expected_logp = -np.log10(pp_expected)
    
    ax_qq.scatter(expected_logp, observed_logp, c='#2C3E50', s=15, alpha=0.7, linewidth=0)
    
    # Identity Line
    max_val = max(np.max(expected_logp), np.max(observed_logp))
    ax_qq.plot([0, max_val], [0, max_val], color='#E74C3C', linestyle='--')
    
    # Confidence Interval (95%)
    # Beta distribution based CI
    # For large N, calculating Beta ppf for every point is slow.
    # We can calculate it for the expected line points.
    
    # Only draw CI if N < 100000 or downsample?
    # Drawing shading is fast enough for ~20k genes.
    
    # k goes from 1 to n (smallest P to largest P)
    # k=1 corresponds to observed_logp[0] (largest value) if we sorted P ascending?
    # No, observed_logp = -log10(sort(P)). P[0] is smallest. -log(P[0]) is largest.
    # expected_logp = -log10( (i-0.5)/n ). i=1 is smallest fraction -> largest log.
    
    # The order matches.
    
    # CI for the i-th order statistic of Uniform(0,1) is Beta(i, n-i+1)
    # i=1..n
    index = np.arange(1, n_points + 1)
    
    # Lower/Upper bounds for P-value
    # alpha=0.95 -> 0.025 and 0.975
    # Warning: beta.ppf might range check.
    
    # Optimization: Calculate CI only for a subset of points for the ribbon
    # or just plotting distinct points.
    
    if n_points < 50000:
        cl = 0.95
        # Beta parameters: a=i, b=n-i+1
        # P_upper = beta.ppf(1-(1-cl)/2, i, n-i+1) which is small P
        # P_lower = beta.ppf((1-cl)/2, i, n-i+1) which is small P (wait)
        
        # P-values:
        # upper_p = stats.beta.ppf(0.975, index, n_points - index + 1)
        # lower_p = stats.beta.ppf(0.025, index, n_points - index + 1)
        
        # Log transformed:
        # lower_log = -np.log10(upper_p) # Because P is small, -log is big
        # upper_log = -np.log10(lower_p)
        
        # But for QQ plot, we usually shade around the diagonal.
        # Actually expected_logp IS the theoretical mean.
        # The CI is around expected_logp.
        
        # Let's use the approximate CI for the null distribution
        # -log10(L) and -log10(U)
        
        lower_p_bound = stats.beta.ppf(0.025, index, n_points - index + 1)
        upper_p_bound = stats.beta.ppf(0.975, index, n_points - index + 1)
        
        lower_log_bound = -np.log10(upper_p_bound) # Map upper P to lower log? No. P=0.9 -> log~0. P=0.01 -> log=2. Upper P corresponds to lower Log.
        upper_log_bound = -np.log10(lower_p_bound) 
        
        # Fill between expects matched x-axis. expected_logp matches index order.
        ax_qq.fill_between(expected_logp, lower_log_bound, upper_log_bound, color='gray', alpha=0.2, label='95% CI')

    # Add Lambda GC Label and N genes
    stats_text = f"$\lambda_{{GC}} = {lambda_val:.3f}$\n$N_{{genes}} = {n_points:,}$"
    ax_qq.text(0.05, 0.95, stats_text, 
               transform=ax_qq.transAxes, fontsize=11, fontweight='bold',
               verticalalignment='top',
               bbox=dict(boxstyle='round', facecolor='white', alpha=0.9, edgecolor='#CCCCCC'))

    ax_qq.set_xlabel(r'Expected $-\log_{10}(P)$')
    ax_qq.set_ylabel(r'Observed $-\log_{10}(P)$')
    ax_qq.set_title("Q-Q Plot")
    
    plt.tight_layout()
    plt.savefig(f"{args.output_prefix}.png", dpi=300)
    plt.savefig(f"{args.output_prefix}.pdf", dpi=300)
    print(f"Saved plots to {args.output_prefix}.png/pdf")

if __name__ == "__main__":
    main()
