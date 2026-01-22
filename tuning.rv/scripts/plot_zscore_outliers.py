import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import os
import argparse
import sys

def parse_args():
    parser = argparse.ArgumentParser(description="Plot MeanDP vs SMinAC and identify robust Z-score outliers.")
    parser.add_argument("--input", "-i", required=True, help="Path to sample metrics file.")
    parser.add_argument("--output-dir", "-o", required=True, help="Directory to save outputs.")
    parser.add_argument("--threshold", "-t", type=float, default=3.0, help="Z-score threshold for outlier detection (default: 3.0).")
    parser.add_argument("--direction", "-d", choices=['both', 'lower', 'upper'], default='both', 
                        help="Direction of outliers to detect: 'lower' (<-T), 'upper' (>T), or 'both' (default).")
    return parser.parse_args()

def main():
    args = parse_args()
    
    input_file = args.input
    output_dir = args.output_dir
    z_thresh = args.threshold
    direction = args.direction

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    print(f"Reading data from {input_file}")
    try:
        df = pd.read_csv(input_file, sep='\t')
    except Exception as e:
        print(f"Error reading file {input_file}: {e}")
        # Try delim_whitespace just in case
        print("Retrying with delim_whitespace=True")
        df = pd.read_csv(input_file, delim_whitespace=True)

    # 2. Robust Z-score for SMinAC
    print("Calculating Robust Z-score for SMinAC...")
    col = 'SMinAC'
    outliers = pd.DataFrame()
    filter_desc = "None"
    
    if col in df.columns:
        median = df[col].median()
        mad = np.median(np.abs(df[col] - median))
        
        # Avoid division by zero
        if mad == 0:
            print("Warning: MAD is 0, cannot calculate robust Z-score properly.")
            df['SMinAC_RobustZ'] = 0
        else:
            # Constant 0.6745 makes MAD consistent with sigma for normal distribution
            k = 0.6745
            df['SMinAC_RobustZ'] = k * (df[col] - median) / mad

        # 3. Identify outliers based on direction and threshold
        print(f"Identifying outliers with threshold Robust Z > {z_thresh} or Robust Z < -{z_thresh} (direction: {direction})")

        if direction == 'lower':
            outliers = df[df['SMinAC_RobustZ'] < -z_thresh]
            filter_desc = f"Robust Z < -{z_thresh}"
        elif direction == 'upper':
            outliers = df[df['SMinAC_RobustZ'] > z_thresh]
            filter_desc = f"Robust Z > {z_thresh}"
        else: # both
            outliers = df[(df['SMinAC_RobustZ'] < -z_thresh) | (df['SMinAC_RobustZ'] > z_thresh)]
            filter_desc = f"|Robust Z| > {z_thresh}"

        print(f"Found {len(outliers)} outliers.")
        
        outliers_path = os.path.join(output_dir, 'SMinAC_outliers.tsv')
        outliers.to_csv(outliers_path, sep='\t', index=False)
        print(f"Outliers saved to {outliers_path}")
        
    else:
        print(f"Error: Column '{col}' not found in dataframe.")

    # 1. Scatter Plot MeanDP vs SMinAC (Color by Group)
    print("Generating scatter plot MeanDP vs SMinAC...")
    
    # Set professional style
    sns.set_theme(style="whitegrid", context="paper", font_scale=1.2)
    
    # Force square figure
    plt.figure(figsize=(8, 8))
    
    # Main scatter plot
    if 'Group' in df.columns:
        ax = sns.scatterplot(data=df, x='MeanDP', y='SMinAC', hue='Group', alpha=0.6, edgecolor=None, s=30)
    else:
        ax = sns.scatterplot(data=df, x='MeanDP', y='SMinAC', alpha=0.6, edgecolor=None, s=30)
        print("Warning: 'Group' column not found, plotting without hue.")

    # Highlight outliers
    if col in df.columns and not outliers.empty:
        plt.scatter(outliers['MeanDP'], outliers['SMinAC'], color='red', s=80, marker='x', linewidth=2, label=f'Outliers ({filter_desc})')
        # Ensure legend handles are correct
        # Get current handles and labels
        handles, labels = ax.get_legend_handles_labels()
        # Add the outlier handle if not automatically added (matplotlib/seaborn interaction can vary)
        plt.legend(title='Group' if 'Group' in df.columns else None, loc='best', frameon=True, framealpha=0.9)

    plt.title(f'MeanDP vs SMinAC\n(Outliers: {filter_desc})', fontsize=14, fontweight='bold', pad=15)
    plt.xlabel(r'Mean Depth ($D_{mean}$)', fontsize=12, fontweight='bold')
    plt.ylabel(r'Sample Minor Allele Burden ($S_{MinAC}$)', fontsize=12, fontweight='bold')
    
    # Enforce square aspect ratio of the plot area if desired, distinct from figure size
    # But user likely just means the image file itself is square.
    # If we want data aspect ratio to optionally be square we could use ax.set_aspect('equal'), 
    # but MeanDP and SMinAC have different units/ranges (25 vs 25000), so 'equal' aspect ratio for data would flatten the plot.
    # So we stick to square figure size.
    
    plt.tight_layout()
    
    plot_path = os.path.join(output_dir, 'MeanDP_vs_SMinAC.png')
    plt.savefig(plot_path, dpi=300) # Higher DPI for professional quality
    plt.close()
    print(f"Plot saved to {plot_path}")

if __name__ == "__main__":
    main()
