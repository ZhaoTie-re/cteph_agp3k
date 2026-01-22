import matplotlib
matplotlib.use('Agg') # Force non-interactive backend for scripts
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns
import argparse
import os
import matplotlib.gridspec as gridspec
import numpy as np
import gc

# Set publication-quality style features
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']
plt.rcParams['svg.fonttype'] = 'none'
sns.set_context("talk", font_scale=1.0) # 'talk' context is better for PPT

# Standard SnpEff-like colors
IMPACT_COLORS = {
    'HIGH': '#d62728',      # Red
    'MODERATE': '#ff7f0e',  # Orange
    'LOW': '#2ca02c',       # Green
    'MODIFIER': '#1f77b4'   # Blue
}
DEFAULT_COLOR = '#7f7f7f'

def format_number(n):
    return f'{int(n):,}'

def setup_page_layout(fig, impacts, include_donut=False):
    """
    Creates a GridSpec layout for the page.
    Returns: gs_left (list of subplotspecs), ax_donut (or None)
    """
    n_rows = len(impacts)
    
    # Single column layout for bars - Donut is now on separate page
    gs_left_col = gridspec.GridSpec(n_rows, 1, hspace=0.6)
    return gs_left_col, None

def plot_impact_bars(ax, impact, subset):
    color = IMPACT_COLORS.get(impact, DEFAULT_COLOR)
    
    # Use raw labels as requested (no beautification)
    labels = subset['Effect'].tolist()
    counts = subset['Count'].tolist()
    
    bars = ax.barh(labels, counts, color=color, alpha=0.9, height=0.6, zorder=3)
    
    max_val = max(counts) if counts else 0
    ax.set_xlim(0, max_val * 1.15) 
    
    total_in_impact = sum(counts)
    for bar in bars:
        width = bar.get_width()
        count_str = format_number(width)
        pct_val = (width / total_in_impact * 100) if total_in_impact > 0 else 0
        
        text_x = width + (max_val * 0.01)
        # Bold text for PPT readability
        ax.text(text_x, bar.get_y() + bar.get_height()/2, 
                f"{count_str} ({pct_val:.1f}%)", 
                va='center', ha='left', fontsize=12, fontweight='bold', color='#404040')
    
    ax.set_title(f"{impact}", loc='left', fontsize=16, fontweight='bold', color=color, pad=15)
    sns.despine(ax=ax, left=True, bottom=False)
    ax.yaxis.set_tick_params(length=0) 
    ax.tick_params(axis='y', labelsize=12) # Larger labels for PPT
    ax.xaxis.grid(True, linestyle='--', alpha=0.5, zorder=0)

def plot_donut_page(fig, agg_df, present_impacts):
    """
    Plots the donut chart on a full separate page with optimized callout positioning
    """
    ax = fig.add_subplot(111)
    
    impact_counts = agg_df.groupby('Impact')['Count'].sum()
    impact_counts = impact_counts.reindex(present_impacts)
    impact_colors_list = [IMPACT_COLORS.get(i, DEFAULT_COLOR) for i in impact_counts.index]
    
    wedges, texts = ax.pie(
        impact_counts, 
        labels=None, 
        startangle=90, 
        colors=impact_colors_list,
        counterclock=False, 
        wedgeprops=dict(width=0.4, edgecolor='white', linewidth=3),
    )

    # Center Text
    total_vars = impact_counts.sum()
    ax.text(0, 0, f"Total\n{format_number(total_vars)}", ha='center', va='center', fontsize=26, fontweight='bold', color='#333333')
    
    ax.set_title("Variant Distribution by Impact", fontsize=30, fontweight='bold', pad=50)
    ax.axis('equal')
    # Manually set limits larger than 'equal' would implies, to accommodate callout labels
    # Pie radius is 1. Anchors are at 1.5. Text extends further.
    ax.set_xlim(-3.5, 3.5) 
    ax.set_ylim(-2.5, 2.5)

    # --- Annotation Layout with Collision Avoidance ---
    annotations = []
    for i, p in enumerate(wedges):
        # Calculate angle and grid coords on unit circle
        ang = (p.theta2 - p.theta1)/2. + p.theta1
        y = np.sin(np.deg2rad(ang))
        x = np.cos(np.deg2rad(ang))
        
        imp_name = impact_counts.index[i]
        color = IMPACT_COLORS.get(imp_name, DEFAULT_COLOR)
        count = impact_counts.values[i]
        pct = count / total_vars
        
        # Formatting label
        label_txt = f"{imp_name}\n{format_number(count)} ({pct:.1%})"
        
        annotations.append({
            'x': x, 'y': y, 
            'text': label_txt, 
            'color': color, 
            'angle': ang
        })

    # Split into Left and Right groups
    right_side = sorted([a for a in annotations if a['x'] >= 0], key=lambda item: item['y'], reverse=True)
    left_side = sorted([a for a in annotations if a['x'] < 0], key=lambda item: item['y'], reverse=True)
    
    def spread_y_positions(items, anchor_x):
        if not items: return
        ys = [item['y'] for item in items]
        min_dist = 0.35 # Sufficient spacing for larger font
        
        # Iterative spacing relaxation
        for _ in range(50):
            changed = False
            for i in range(len(ys) - 1):
                # distance between current and next (since sorted desc, ys[i] > ys[i+1])
                diff = ys[i] - ys[i+1]
                if diff < min_dist:
                    center = (ys[i] + ys[i+1]) / 2
                    shift = (min_dist - diff) / 2 + 0.01
                    ys[i] = center + shift
                    ys[i+1] = center - shift
                    changed = True
            if not changed:
                break
        
        for idx, item in enumerate(items):
            item['label_pos'] = (anchor_x, ys[idx])

    # Apply spacing logic: Push labels further out to 1.5 radius
    spread_y_positions(right_side, 1.5)
    spread_y_positions(left_side, -1.5)
    
    # Render Annotations
    for item in right_side + left_side:
        lx, ly = item['label_pos']
        ha = 'left' if lx > 0 else 'right'
        
        # "angle" connection style: horizontal line from text -> angled line to point
        connection = f"angle,angleA=0,angleB={item['angle']}"
        
        ax.annotate(item['text'], 
                    xy=(item['x'], item['y']), 
                    xytext=(lx, ly),
                    horizontalalignment=ha,
                    verticalalignment='center',
                    fontsize=18,          # Reasonably larger
                    fontweight='bold', 
                    color=item['color'],  # Match specific impact color
                    arrowprops=dict(arrowstyle="-", 
                                    color=item['color'], 
                                    lw=2,
                                    connectionstyle=connection))

def main():
    parser = argparse.ArgumentParser(description="Plot SnpEff stats (PPT Quality)")
    parser.add_argument("--input", required=True, help="Input TSV file")
    parser.add_argument("--output", required=True, help="Output plot prefix")
    args = parser.parse_args()

    try:
        df = pd.read_csv(args.input, sep='\t')
    except Exception as e:
        print(f"Error reading file: {e}")
        return

    required_cols = ['Impact', 'Effect', 'Count']
    if not all(col in df.columns for col in required_cols):
        print(f"Missing columns. Found: {df.columns}")
        return

    df = df[df['Impact'] != 'OTHER']
    agg_df = df.groupby(['Impact', 'Effect'])['Count'].sum().reset_index()

    impact_order = ['HIGH', 'MODERATE', 'LOW', 'MODIFIER']
    # Filter only impacts present in data
    present_impacts = [i for i in impact_order if i in agg_df['Impact'].unique()]
    
    if len(present_impacts) == 0:
        print("No data to plot found.")
        return

    # --- Limit Categories Logic ---
    MAX_EFFECTS = 15
    impact_data_map = {}
    
    # Track max label length to adjust margins dynamically
    max_label_len = 0

    for impact in present_impacts:
        subset = agg_df[agg_df['Impact'] == impact].sort_values('Count', ascending=True)
        
        if len(subset) > MAX_EFFECTS:
            top_subset = subset.tail(MAX_EFFECTS).copy()
            rest_df = subset.iloc[:-MAX_EFFECTS]
            rest_count = rest_df['Count'].sum()
            rest_types = len(rest_df)
            
            rest_row = pd.DataFrame({
                'Impact': [impact], 
                'Effect': [f'Other ({rest_types} types)'], 
                'Count': [rest_count]
            })
            plot_subset = pd.concat([rest_row, top_subset])
        else:
            plot_subset = subset
            
        impact_data_map[impact] = plot_subset
        
        # Check label lengths in this subset
        if not plot_subset.empty:
            current_max = plot_subset['Effect'].str.len().max()
            if current_max > max_label_len:
                max_label_len = current_max

    # --- Output PDF with Multiple Pages ---
    output_pdf = f"{args.output}.pdf"
    
    # --- Layout Logic for Bar Charts ---
    # Width increased to 20 inches to accommodate very long labels
    FIG_WIDTH = 20
    FIG_HEIGHT = 12
    
    # Calculate margin as a fraction of figure width
    # Estimate: char width ~ 0.12 inch (for fontsize 12-14)
    # Total label width = max_chars * 0.12
    # Margin Fraction = (Total Label Width) / FIG_WIDTH
    # + Buffer (0.05)
    
    estimated_label_width_inches = max_label_len * 0.13
    calculated_left_margin = (estimated_label_width_inches / FIG_WIDTH) + 0.05
    
    # Cap margin between 0.35 and 0.75
    calculated_left_margin = max(0.35, min(0.75, calculated_left_margin))
    
    print(f"Max label chars: {max_label_len}. Est width: {estimated_label_width_inches:.1f}in. Margin: {calculated_left_margin:.2f}")

    with PdfPages(output_pdf) as pdf:
        
        # Page 1: Donut Chart with Callouts
        # Use large figure size and wide limits to prevent text clipping
        fig_donut = plt.figure(figsize=(20, 14))
        plot_donut_page(fig_donut, agg_df, present_impacts)
        # Maximimize axes area, letting set_xlim control the "zoom"
        plt.subplots_adjust(left=0.02, right=0.98, top=0.92, bottom=0.08)
        pdf.savefig(fig_donut)
        plt.close(fig_donut)

        # Page 2: HIGH and MODERATE
        page1_impacts = [i for i in ['HIGH', 'MODERATE'] if i in present_impacts]
        
        if page1_impacts:
            fig1 = plt.figure(figsize=(FIG_WIDTH, FIG_HEIGHT)) 
            gs_left, _ = setup_page_layout(fig1, page1_impacts)
            
            for idx, impact in enumerate(page1_impacts):
                ax = fig1.add_subplot(gs_left[idx])
                plot_impact_bars(ax, impact, impact_data_map[impact])
                if idx == len(page1_impacts) - 1:
                    ax.set_xlabel("Count", fontweight='bold', fontsize=14)
            
            plt.subplots_adjust(left=calculated_left_margin, right=0.95, top=0.92, bottom=0.08)
            pdf.savefig(fig1)
            plt.close(fig1)

        # Page 3: LOW and MODIFIER
        page2_impacts = [i for i in ['LOW', 'MODIFIER'] if i in present_impacts]
        
        if page2_impacts:
            fig2 = plt.figure(figsize=(FIG_WIDTH, FIG_HEIGHT))
            gs_left, _ = setup_page_layout(fig2, page2_impacts)
            
            for idx, impact in enumerate(page2_impacts):
                ax = fig2.add_subplot(gs_left[idx])
                plot_impact_bars(ax, impact, impact_data_map[impact])
                if idx == len(page2_impacts) - 1:
                    ax.set_xlabel("Count", fontweight='bold', fontsize=14)
            
            plt.subplots_adjust(left=calculated_left_margin, right=0.95, top=0.92, bottom=0.08)
            pdf.savefig(fig2)
            plt.close(fig2)
            
    print(f"PDF saved to {output_pdf}")
    gc.collect()

if __name__ == "__main__":
    main()