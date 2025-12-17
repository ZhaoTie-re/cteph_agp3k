#!/usr/bin/env python3
"""
Filter variants based on quality metrics with streaming processing.
Memory-efficient for large files (42M+ variants).
"""

import argparse
import json
import pandas as pd
import sys
import os
import time
import warnings
from datetime import datetime
from collections import defaultdict
import matplotlib.pyplot as plt
import seaborn as sns
from upsetplot import UpSet, from_memberships
import gzip

# Suppress FutureWarnings from upsetplot library (third-party dependency)
warnings.filterwarnings('ignore', category=FutureWarning, module='upsetplot')
warnings.filterwarnings('ignore', category=FutureWarning, message='.*fillna.*inplace.*')
warnings.filterwarnings('ignore', category=FutureWarning, message='.*Downcasting object dtype.*')


def log(message, level='INFO'):
    """Print log message with timestamp."""
    timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    print(f"[{timestamp}] [{level}] {message}", file=sys.stderr, flush=True)


def parse_filter_expression(expr):
    """
    Parse filter expression like 'VMISS_CASE >= 0.02' into components.
    Returns: (metric, operator, threshold)
    """
    import re
    
    # Match pattern: METRIC OPERATOR VALUE
    pattern = r'^\s*(\w+)\s*(>=|<=|>|<|==|!=)\s*([0-9.eE+-]+)\s*$'
    match = re.match(pattern, expr.strip())
    
    if not match:
        raise ValueError(f"Invalid filter expression: '{expr}'. Expected format: 'METRIC OPERATOR VALUE'")
    
    metric = match.group(1)
    operator = match.group(2)
    threshold = float(match.group(3))
    
    return metric, operator, threshold


def load_filter_config(config_file):
    """Load and validate filter configuration."""
    log(f"Loading filter configuration from {config_file}")
    
    with open(config_file, 'r') as f:
        config = json.load(f)
    
    filter_expressions = config.get('filters', [])
    
    if isinstance(filter_expressions, dict):
        # Old format compatibility
        active_filters = {k: v for k, v in filter_expressions.items() if v is not None}
    else:
        # New format: list of expressions
        active_filters = {}
        for expr in filter_expressions:
            metric, operator, threshold = parse_filter_expression(expr)
            active_filters[metric] = {'operator': operator, 'threshold': threshold}
    
    log(f"  ✓ Loaded {len(active_filters)} active filters")
    for metric, params in active_filters.items():
        log(f"    - {metric} {params['operator']} {params['threshold']}")
    
    return active_filters


def apply_filter(value, operator, threshold):
    """
    Apply comparison operator between value and threshold.
    Handles missing values (None/NaN).
    """
    if pd.isna(value) or value is None:
        return False
    
    try:
        value = float(value)
        threshold = float(threshold)
    except (ValueError, TypeError):
        return False
    
    if operator == '<':
        return value < threshold
    elif operator == '<=':
        return value <= threshold
    elif operator == '>':
        return value > threshold
    elif operator == '>=':
        return value >= threshold
    elif operator == '==':
        return value == threshold
    elif operator == '!=':
        return value != threshold
    else:
        raise ValueError(f"Unknown operator: {operator}")


def stream_filter_variants(input_file, filters, output_prefix, chunk_size=500000):
    """
    Stream through variant stats file, apply filters, and track which variants fail.
    Memory-efficient: processes in chunks, never loads entire file.
    
    Returns:
        dict: Statistics about filtered variants
    """
    log(f"Starting streaming variant filtering...")
    log(f"  Input file: {input_file}")
    log(f"  Chunk size: {chunk_size:,} variants")
    
    start_time = time.time()
    
    # Track filtering statistics
    total_variants = 0
    filtered_variants = set()  # IDs of variants to remove
    filter_counts = defaultdict(int)  # Count per filter
    variant_filter_membership = defaultdict(set)  # Which filters each variant fails
    
    # Open output file for excluded variant IDs
    exclude_file = f"{output_prefix}.exclude.txt"
    
    # Read file in chunks
    log("Processing variants in chunks...")
    chunk_num = 0
    
    try:
        # Use decompression based on file extension
        if input_file.endswith('.gz'):
            reader = pd.read_csv(input_file, sep='\t', compression='gzip', 
                               chunksize=chunk_size, dtype={'#CHROM': str})
        else:
            reader = pd.read_csv(input_file, sep='\t', chunksize=chunk_size, 
                               dtype={'#CHROM': str})
        
        for chunk_df in reader:
            chunk_num += 1
            chunk_start = time.time()
            chunk_variants = len(chunk_df)
            total_variants += chunk_variants
            
            # Check each filter
            for metric, params in filters.items():
                if metric not in chunk_df.columns:
                    log(f"  ⚠ WARNING: Metric '{metric}' not found in data", level='WARNING')
                    continue
                
                operator = params['operator']
                threshold = params['threshold']
                
                # Apply filter to this chunk
                mask = chunk_df[metric].apply(
                    lambda x: apply_filter(x, operator, threshold)
                )
                
                failed_ids = chunk_df.loc[mask, 'ID'].tolist()
                filter_counts[metric] += len(failed_ids)
                
                # Track membership for upset plot
                for vid in failed_ids:
                    filtered_variants.add(vid)
                    variant_filter_membership[vid].add(metric)
            
            chunk_elapsed = time.time() - chunk_start
            log(f"  Chunk {chunk_num}: Processed {chunk_variants:,} variants "
                f"in {chunk_elapsed:.2f}s (cumulative: {total_variants:,})")
    
    except Exception as e:
        log(f"Error during streaming: {str(e)}", level='ERROR')
        raise
    
    # Write excluded variants to file with natural sorting
    log(f"Writing excluded variant IDs to {exclude_file}...")
    
    # Natural sort: chr number first, then position
    def natural_sort_key(variant_id):
        """
        Sort variants by chromosome (natural order) and position.
        Handles format: chr:pos:ref:alt or chrN:pos:ref:alt
        """
        try:
            parts = variant_id.split(':')
            if len(parts) >= 2:
                chrom = parts[0]
                pos = int(parts[1])
                
                # Extract numeric part of chromosome
                # Handle both '1' and 'chr1' formats
                chrom_clean = chrom.replace('chr', '')
                
                # Try to convert to int, fallback to 99 for non-numeric (X, Y, MT, etc.)
                try:
                    chrom_num = int(chrom_clean)
                except ValueError:
                    # Non-numeric chromosomes go to the end
                    chrom_num = 99
                
                return (chrom_num, pos)
            else:
                # Fallback for malformed IDs
                return (99, 0)
        except:
            return (99, 0)
    
    with open(exclude_file, 'w') as f:
        for vid in sorted(filtered_variants, key=natural_sort_key):
            f.write(f"{vid}\n")
    
    elapsed = time.time() - start_time
    log(f"  ✓ Filtering complete in {elapsed:.2f}s ({elapsed/60:.2f} minutes)")
    log(f"  Total variants processed: {total_variants:,}")
    log(f"  Variants to exclude: {len(filtered_variants):,} ({len(filtered_variants)/total_variants*100:.2f}%)")
    
    # Prepare statistics
    stats = {
        'total_variants': total_variants,
        'excluded_variants': len(filtered_variants),
        'excluded_percentage': len(filtered_variants) / total_variants * 100 if total_variants > 0 else 0,
        'filter_counts': dict(filter_counts),
        'variant_filter_membership': variant_filter_membership,
        'processing_time': elapsed
    }
    
    return stats


def create_upset_plot(variant_filter_membership, output_file, filters):
    """
    Create publication-quality UpSet plot showing overlap between different filters.
    """
    log("Creating publication-quality UpSet plot...")
    
    if not variant_filter_membership:
        log("  No variants to plot (all passed filters)")
        return
    
    try:
        # Set publication-quality style
        plt.rcParams['font.family'] = 'sans-serif'
        plt.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans', 'Liberation Sans']
        plt.rcParams['font.size'] = 12
        plt.rcParams['axes.labelsize'] = 13
        plt.rcParams['axes.titlesize'] = 14
        plt.rcParams['xtick.labelsize'] = 11
        plt.rcParams['ytick.labelsize'] = 11
        plt.rcParams['legend.fontsize'] = 11
        plt.rcParams['figure.titlesize'] = 16
        plt.rcParams['axes.linewidth'] = 1.5
        plt.rcParams['xtick.major.width'] = 1.5
        plt.rcParams['ytick.major.width'] = 1.5
        
        # Convert to format suitable for upsetplot
        memberships = list(variant_filter_membership.values())
        
        # Create UpSet data (automatically preserves category names from memberships)
        upset_data = from_memberships(memberships)
        
        # Create figure with optimal layout (wider for better spacing)
        fig = plt.figure(figsize=(20, 11))
        
        # Create UpSet plot with enhanced styling
        # 关键：不使用 sort_categories_by 参数，让库自动处理类别名称
        upset = UpSet(
            upset_data,
            subset_size='count',
            show_counts=True,
            show_percentages=True,
            element_size=55,
            intersection_plot_elements=12,
            facecolor='#2E86AB',
            shading_color='#A9A9A9',
            with_lines=True,
            sort_by='cardinality'
        )
        
        upset.plot(fig=fig)
        
        # 简单格式化：只加粗Y轴类别标签
        for ax in fig.get_axes():
            if hasattr(ax, 'yaxis'):
                for label in ax.get_yticklabels():
                    label.set_fontweight('bold')
                    label.set_fontsize(11)
        
        # Enhance title with better spacing
        fig.suptitle(
            'Quality Control Filter Overlap Analysis',
            fontsize=18,
            fontweight='bold',
            y=0.96,
            x=0.5,
            ha='center'
        )
        
        # Add subtitle with filter count (formatted with comma)
        n_filters = len(filters)
        n_variants = len(variant_filter_membership)
        fig.text(
            0.5, 0.925,
            f'{n_filters} filters applied • {n_variants:,} variants excluded',
            ha='center',
            fontsize=13,
            color='#555555',
            style='italic'
        )
        
        # Add filter descriptions in a compact box (left side, next to "Intersection size" label)
        filter_desc = ["Applied Filters:"]
        filter_desc.append("─" * 45)
        for i, (metric, params) in enumerate(filters.items(), 1):
            threshold_str = f"{params['threshold']:.2e}" if params['threshold'] < 0.01 else f"{params['threshold']:.4f}"
            filter_desc.append(f"{i}. {metric:<15} {params['operator']:<3} {threshold_str}")
        
        desc_text = "\n".join(filter_desc)
        
        # Position at left side, vertically centered with main plot
        fig.text(
            0.005, 0.50,
            desc_text,
            fontsize=10,
            verticalalignment='center',
            horizontalalignment='left',
            bbox=dict(
                boxstyle='round,pad=0.7',
                facecolor='#F8F9FA',
                alpha=0.95,
                edgecolor='#2E86AB',
                linewidth=2
            ),
            family='monospace',
            transform=fig.transFigure
        )
        
        # Adjust layout with proper margins to accommodate the filter box on the left
        plt.subplots_adjust(
            left=0.20,  # Space on left for filter box
            right=0.97,
            top=0.88,   # Space at top for title and subtitle
            bottom=0.08,
            hspace=0.5,
            wspace=0.4
        )
        
        # Save with high quality settings
        plt.savefig(
            output_file,
            dpi=600,
            bbox_inches='tight',
            facecolor='white',
            edgecolor='none',
            format='png',
            pad_inches=0.3
        )
        
        # Also save as PDF for publication
        pdf_file = output_file.replace('.png', '.pdf')
        plt.savefig(
            pdf_file,
            bbox_inches='tight',
            facecolor='white',
            edgecolor='none',
            format='pdf',
            pad_inches=0.3
        )
        
        plt.close()
        
        log(f"  ✓ UpSet plot saved to {output_file} (PNG, 600 DPI)")
        log(f"  ✓ UpSet plot saved to {pdf_file} (PDF, vector format)")
        
    except Exception as e:
        log(f"  ⚠ Failed to create UpSet plot: {str(e)}", level='WARNING')
        import traceback
        log(f"  Traceback: {traceback.format_exc()}", level='WARNING')


def create_summary_report(stats, filters, output_file):
    """
    Create detailed text summary of filtering results.
    """
    log(f"Creating summary report: {output_file}")
    
    with open(output_file, 'w') as f:
        f.write("="*80 + "\n")
        f.write("VARIANT FILTERING SUMMARY REPORT\n")
        f.write("="*80 + "\n\n")
        
        f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        
        # Overall statistics
        f.write("OVERALL STATISTICS\n")
        f.write("-"*80 + "\n")
        f.write(f"Total variants processed:    {stats['total_variants']:>15,}\n")
        f.write(f"Variants excluded:           {stats['excluded_variants']:>15,} "
                f"({stats['excluded_percentage']:.2f}%)\n")
        f.write(f"Variants retained:           {stats['total_variants'] - stats['excluded_variants']:>15,} "
                f"({100 - stats['excluded_percentage']:.2f}%)\n")
        f.write(f"Processing time:             {stats['processing_time']:>15.2f}s "
                f"({stats['processing_time']/60:.2f} min)\n\n")
        
        # Filter-specific statistics
        f.write("FILTER-SPECIFIC COUNTS (with OR logic)\n")
        f.write("-"*80 + "\n")
        f.write(f"{'Filter Metric':<20} {'Operator':<10} {'Threshold':<15} {'Failed Count':<15} {'% of Total':<10}\n")
        f.write("-"*80 + "\n")
        
        for metric, params in filters.items():
            count = stats['filter_counts'].get(metric, 0)
            pct = count / stats['total_variants'] * 100 if stats['total_variants'] > 0 else 0
            f.write(f"{metric:<20} {params['operator']:<10} {params['threshold']:<15.2e} "
                   f"{count:<15,} {pct:<10.2f}%\n")
        
        f.write("\n")
        f.write("="*80 + "\n")
        f.write("NOTE: Multiple filters use OR logic - variants failing ANY filter are excluded.\n")
        f.write("      The sum of filter counts may exceed total excluded variants due to overlap.\n")
        f.write("      See UpSet plot for visualization of filter overlaps.\n")
        f.write("="*80 + "\n")
    
    log(f"  ✓ Summary report saved")


def main():
    parser = argparse.ArgumentParser(
        description='Filter variants based on quality metrics (streaming mode)',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    parser.add_argument('--input', required=True,
                       help='Input variant stats TSV.gz file (annotated)')
    parser.add_argument('--config', required=True,
                       help='JSON configuration file with filter thresholds')
    parser.add_argument('--output-prefix', required=True,
                       help='Output prefix for filtered variant list and plots')
    parser.add_argument('--chunk-size', type=int, default=500000,
                       help='Number of variants to process per chunk (default: 500k)')
    
    args = parser.parse_args()
    
    log("="*80)
    log("Variant Quality Filtering Pipeline (Streaming Mode)")
    log("="*80)
    log(f"Input file: {args.input}")
    log(f"Config file: {args.config}")
    log(f"Output prefix: {args.output_prefix}")
    log("="*80)
    
    overall_start = time.time()
    
    # Load filter configuration
    filters = load_filter_config(args.config)
    
    if not filters:
        log("⚠ No active filters found in configuration. No variants will be excluded.", level='WARNING')
        # Create empty exclude file
        with open(f"{args.output_prefix}.exclude.txt", 'w') as f:
            pass
        return
    
    # Apply filters with streaming
    stats = stream_filter_variants(args.input, filters, args.output_prefix, args.chunk_size)
    
    # Create visualizations and reports
    log("\nGenerating visualizations and reports...")
    
    # UpSet plot
    upset_file = f"{args.output_prefix}.upset.png"
    create_upset_plot(stats['variant_filter_membership'], upset_file, filters)
    
    # Summary report
    summary_file = f"{args.output_prefix}.summary.txt"
    create_summary_report(stats, filters, summary_file)
    
    overall_elapsed = time.time() - overall_start
    log("="*80)
    log(f"Pipeline completed successfully!")
    log(f"Total time: {overall_elapsed:.2f}s ({overall_elapsed/60:.2f} minutes)")
    log(f"Output files:")
    log(f"  - Exclude list:      {args.output_prefix}.exclude.txt")
    log(f"  - UpSet plot (PNG):  {args.output_prefix}.upset.png (600 DPI)")
    log(f"  - UpSet plot (PDF):  {args.output_prefix}.upset.pdf (vector)")
    log(f"  - Summary report:    {args.output_prefix}.summary.txt")
    log("="*80)


if __name__ == '__main__':
    main()
