#!/usr/bin/env python3
"""
Generate distribution plots for MQ and VQSLOD from extracted TSV files
Uses chunked reading and numpy for memory efficiency with large datasets
"""

import sys
import argparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import time
from pathlib import Path

# Add script directory to path for imports
sys.path.insert(0, str(Path(__file__).parent))
from logger_utils import setup_logging

def read_metrics_from_tsv(tsv_file, chunk_size=500000):
    """
    Read MQ and VQSLOD values from TSV file using chunked reading
    to prevent memory overflow with very large datasets
    Chunk size of 500K provides ~10 chunks for 5M variants per chromosome
    Returns: (mq_values, vqslod_values) as numpy arrays
    """
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] Reading TSV file: {tsv_file}")
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] Using chunked reading (chunk size: {chunk_size:,})")
    
    mq_chunks = []
    vqslod_chunks = []
    
    start_time = time.time()
    total_variants = 0
    
    with open(tsv_file, 'r') as f:
        # Skip header
        header = f.readline().strip()
        if header != "MQ\tVQSLOD":
            print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] Warning: Unexpected header: {header}")
            print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}]   Expected: MQ\\tVQSLOD")
        
        print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] Reading variant metrics in chunks...")
        
        mq_buffer = []
        vqslod_buffer = []
        
        for line_num, line in enumerate(f, start=1):
            line = line.strip()
            if not line:
                continue
            
            parts = line.split('\t')
            if len(parts) != 2:
                continue
            
            try:
                mq_str = parts[0]
                vqslod_str = parts[1]
                
                if ',' in mq_str:
                    valid_vals = []
                    for x in mq_str.split(','):
                        try:
                            valid_vals.append(float(x))
                        except ValueError:
                            continue
                    mq_val = max(valid_vals) if valid_vals else float('nan')
                else:
                    mq_val = float(mq_str)
                    
                if ',' in vqslod_str:
                    valid_vals = []
                    for x in vqslod_str.split(','):
                        try:
                            valid_vals.append(float(x))
                        except ValueError:
                            continue
                    vqslod_val = max(valid_vals) if valid_vals else float('nan')
                else:
                    vqslod_val = float(vqslod_str)
                    
                mq_buffer.append(mq_val)
                vqslod_buffer.append(vqslod_val)
            except ValueError:
                continue
            
            # Process chunk when buffer is full
            if len(mq_buffer) >= chunk_size:
                mq_chunks.append(np.array(mq_buffer, dtype=np.float64))
                vqslod_chunks.append(np.array(vqslod_buffer, dtype=np.float64))
                total_variants += len(mq_buffer)
                
                elapsed = time.time() - start_time
                rate = total_variants / elapsed
                print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}]   Processed {total_variants:,} variants ({rate:.0f} var/s)")
                
                mq_buffer = []
                vqslod_buffer = []
        
        # Process remaining data
        if mq_buffer:
            mq_chunks.append(np.array(mq_buffer, dtype=np.float64))
            vqslod_chunks.append(np.array(vqslod_buffer, dtype=np.float64))
            total_variants += len(mq_buffer)
    
    elapsed = time.time() - start_time
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] Reading completed in {elapsed:.2f}s")
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}]   Total variants: {total_variants:,}")
    
    # Handle case where no data was found
    if not mq_chunks:
        print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] Warning: No valid data found in TSV file")
        return np.array([], dtype=np.float64), np.array([], dtype=np.float64)
    
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] Concatenating {len(mq_chunks)} chunks...")
    
    # Concatenate all chunks efficiently
    mq_values = np.concatenate(mq_chunks) if len(mq_chunks) > 1 else mq_chunks[0]
    vqslod_values = np.concatenate(vqslod_chunks) if len(vqslod_chunks) > 1 else vqslod_chunks[0]
    
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] Data ready for plotting")
    
    # Filter out invalid values (inf, -inf, nan)
    mq_valid = np.isfinite(mq_values)
    vqslod_valid = np.isfinite(vqslod_values)
    
    if not np.all(mq_valid):
        n_invalid = np.sum(~mq_valid)
        print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] Warning: Filtered {n_invalid:,} invalid MQ values (inf/nan)")
        mq_values = mq_values[mq_valid]
    
    if not np.all(vqslod_valid):
        n_invalid = np.sum(~vqslod_valid)
        print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] Warning: Filtered {n_invalid:,} invalid VQSLOD values (inf/nan)")
        vqslod_values = vqslod_values[vqslod_valid]
    
    return mq_values, vqslod_values

def plot_distribution(values, threshold, title, xlabel, output_file):
    """Create publication-quality distribution plot with threshold line"""
    import time
    
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] Generating plot: {output_file}")
    
    if len(values) == 0:
        print(f"Warning: No values found for {title}")
        # Create empty plot
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.text(0.5, 0.5, 'No data available', ha='center', va='center', fontsize=14)
        ax.set_title(title, fontsize=16, fontweight='bold')
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        return
    
    # Calculate statistics using numpy for efficiency
    total_count = len(values)
    pass_count = np.sum(values > threshold)  # Vectorized comparison
    pass_percentage = (pass_count / total_count * 100) if total_count > 0 else 0
    
    mean_val = np.mean(values)
    median_val = np.median(values)
    std_val = np.std(values)
    
    # Set publication-quality style
    plt.style.use('seaborn-v0_8-paper')
    fig, ax = plt.subplots(figsize=(10, 7))
    
    # Create histogram with optimized bins
    bins = min(100, int(np.sqrt(len(values))) * 5)
    n, bins_edges, patches = ax.hist(values, bins=bins, alpha=0.75, color='#2E86AB', 
                                      edgecolor='#1A1A1A', linewidth=0.5, density=False)
    
    # Add threshold line with annotation
    ymax = ax.get_ylim()[1]
    ax.axvline(x=threshold, color='#E63946', linestyle='--', linewidth=2.5, 
               label=f'Threshold = {threshold}', zorder=10)
    
    # Add shaded region for passing variants
    ax.axvspan(threshold, values.max(), alpha=0.15, color='#06A77D', 
               label='Pass region', zorder=1)
    
    # Set labels and title with proper font sizes
    ax.set_xlabel(xlabel, fontsize=13, fontweight='bold')
    ax.set_ylabel('Frequency', fontsize=13, fontweight='bold')
    ax.set_title(title, fontsize=15, fontweight='bold', pad=40)
    
    # Improve tick labels
    ax.tick_params(axis='both', which='major', labelsize=11, width=1.2)
    ax.tick_params(axis='both', which='minor', width=0.8)
    
    # Add grid for better readability
    ax.grid(True, alpha=0.3, linestyle='--', linewidth=0.8, which='major')
    ax.set_axisbelow(True)
    
    # Add legend outside plot area, below title (centered), without shadow
    legend = ax.legend(fontsize=11, loc='lower center', bbox_to_anchor=(0.5, 1.02),
                       ncol=2, frameon=True, fancybox=False, shadow=False, framealpha=0.95)
    legend.get_frame().set_edgecolor('#1A1A1A')
    legend.get_frame().set_linewidth(1.2)
    
    # Set spine properties
    for spine in ax.spines.values():
        spine.set_linewidth(1.2)
        spine.set_edgecolor('#1A1A1A')
    
    # Tight layout
    plt.tight_layout()
    
    # Save with high DPI for publication
    plt.savefig(output_file, dpi=300, bbox_inches='tight', facecolor='white', 
                edgecolor='none')
    plt.close()
    
    import time
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] Saved publication-quality plot: {output_file}")

def main():
    parser = argparse.ArgumentParser(
        description='Generate MQ and VQSLOD distribution plots from extracted TSV files',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example:
  %(prog)s chr1.metrics.tsv chr1 58.75 10 chr1_output
        """
    )
    
    parser.add_argument('tsv_file', 
                       help='Input TSV file with MQ and VQSLOD values')
    parser.add_argument('chromosome', 
                       help='Chromosome name for labeling (e.g., chr1, chr22, chrX)')
    parser.add_argument('mq_threshold', type=float,
                       help='Mapping Quality (MQ) threshold for filtering')
    parser.add_argument('vqslod_threshold', type=float,
                       help='VQSLOD score threshold for filtering')
    parser.add_argument('output_prefix', 
                       help='Output prefix for plot files (e.g., chr1)')
    
    args = parser.parse_args()
    
    # Setup dual logging (stdout + file)
    log_file = setup_logging(f"{args.output_prefix}.plot.log")
    
    print(f"\n{'='*70}")
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] Starting plot generation for {args.chromosome}")
    print(f"{'='*70}")
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] Input TSV: {args.tsv_file}")
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] MQ threshold: {args.mq_threshold}, VQSLOD threshold: {args.vqslod_threshold}")
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] Log file: {log_file}")
    print(f"{'='*70}\n")
    
    # Read values from TSV
    mq_values, vqslod_values = read_metrics_from_tsv(args.tsv_file)
    
    print(f"\n[{time.strftime('%Y-%m-%d %H:%M:%S')}] Data summary:")
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}]   MQ values: {len(mq_values):,}")
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}]   VQSLOD values: {len(vqslod_values):,}\n")
    
    # Create MQ distribution plot
    mq_output = f"{args.output_prefix}.MQ_distribution.png"
    plot_distribution(mq_values, args.mq_threshold, 
                     f'MQ Distribution - {args.chromosome}',
                     'Mapping Quality (MQ)', mq_output)
    
    # Create VQSLOD distribution plot
    vqslod_output = f"{args.output_prefix}.VQSLOD_distribution.png"
    plot_distribution(vqslod_values, args.vqslod_threshold,
                     f'VQSLOD Distribution - {args.chromosome}',
                     'VQSLOD Score', vqslod_output)
    
    print(f"\n{'='*70}")
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] Plot generation completed for {args.chromosome}")
    print(f"{'='*70}\n")

if __name__ == "__main__":
    main()
