#!/usr/bin/env python3
"""
Lead Variant and Peak Region Extraction Tool

This script identifies lead variants (most significant SNPs) from GWAS results
and extracts variant IDs within a specified range around each lead variant.

Author: ZHAO TIE
Date: December 2025
"""

import argparse
import pandas as pd
import numpy as np
import sys
import os
from datetime import datetime


def parse_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description='Identify lead variants and extract peak regions from GWAS results',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example usage:
    python get_peak_cli.py \\
        --assoc results.glm.logistic \\
        --sig-level 5e-8 \\
        --lead-range 500000 \\
        --output-dir ./peaks \\
        --prefix peaks
        """
    )
    
    parser.add_argument('--assoc', required=True,
                        help='Path to association results file (PLINK2 .glm.logistic format)')
    parser.add_argument('--sig-level', type=float, default=5e-8,
                        help='Genome-wide significance threshold (default: 5e-8)')
    parser.add_argument('--lead-range', type=int, default=500000,
                        help='Range (bp) around lead variant to extract (default: 500000)')
    parser.add_argument('--window', type=int, default=500000,
                        help='Window size (bp) for merging nearby lead variants. Only the variant with minimum P-value within this window is kept as lead (default: 500000)')
    parser.add_argument('--output-dir', default='.',
                        help='Output directory (default: current directory)')
    parser.add_argument('--prefix', default='peaks',
                        help='Prefix for output files (default: peaks)')
    
    return parser.parse_args()


def load_association_data(assoc_path):
    """Load association results"""
    print(f"Loading association data from: {assoc_path}")
    
    # Read association file
    df = pd.read_csv(assoc_path, sep='\t')
    print(f"  Total variants: {len(df):,}")
    
    # Check required columns
    required_cols = ['ID', '#CHROM', 'POS', 'P']
    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        raise ValueError(f"Missing required columns: {missing_cols}")
    
    # Extract relevant columns
    data = df[['ID', '#CHROM', 'POS', 'P']].copy()
    data.columns = ['VARIANT_ID', 'CHR', 'POS', 'P']
    
    # Convert P-values to numeric, handling '.' and NA
    data['P'] = pd.to_numeric(data['P'].replace('.', np.nan), errors='coerce')
    
    # Remove invalid P-values
    valid_mask = data['P'].notna() & (data['P'] > 0) & (data['P'] <= 1)
    invalid_count = (~valid_mask).sum()
    if invalid_count > 0:
        print(f"  Removing {invalid_count:,} variants with invalid P-values")
        data = data[valid_mask].copy()
    
    print(f"  Valid variants: {len(data):,}")
    
    return data


def identify_lead_variants(data, sig_level, window):
    """Identify lead variants from significant associations
    
    Within each window, only the variant with minimum P-value is kept as lead.
    """
    print(f"\nIdentifying lead variants (P < {sig_level:.2e})...")
    print(f"  Window size for merging: {window:,} bp")
    
    # Filter significant variants
    sig_data = data[data['P'] < sig_level].copy()
    print(f"  Significant variants: {len(sig_data):,}")
    
    if len(sig_data) == 0:
        print("  No significant variants found!")
        return pd.DataFrame()
    
    # Sort by chromosome and position
    sig_data = sig_data.sort_values(['CHR', 'POS']).reset_index(drop=True)
    
    # Identify independent peaks (lead variants)
    lead_variants = []
    
    for chrom in sig_data['CHR'].unique():
        chrom_data = sig_data[sig_data['CHR'] == chrom].copy()
        
        while len(chrom_data) > 0:
            # Find variant with minimum P-value
            min_idx = chrom_data['P'].idxmin()
            lead_var = chrom_data.loc[min_idx]
            
            lead_variants.append({
                'CHR': lead_var['CHR'],
                'POS': lead_var['POS'],
                'VARIANT_ID': lead_var['VARIANT_ID'],
                'P': lead_var['P'],
                'REGION_START': max(1, lead_var['POS'] - window // 2),
                'REGION_END': lead_var['POS'] + window // 2
            })
            
            # Remove all variants within window (keeping only minimum P-value as lead)
            chrom_data = chrom_data[
                (chrom_data['POS'] < lead_var['POS'] - window // 2) |
                (chrom_data['POS'] > lead_var['POS'] + window // 2)
            ]
    
    lead_df = pd.DataFrame(lead_variants)
    lead_df = lead_df.sort_values(['CHR', 'POS']).reset_index(drop=True)
    
    print(f"  Identified {len(lead_df)} independent lead variants")
    
    return lead_df


def extract_peak_regions(data, lead_df, lead_range):
    """Extract variant IDs within lead_range/2 around each lead variant"""
    print(f"\nExtracting variants within ±{lead_range//2:,} bp of lead variants...")
    
    half_range = lead_range // 2
    peak_regions = []
    
    for idx, lead in lead_df.iterrows():
        region_start = max(1, lead['POS'] - half_range)
        region_end = lead['POS'] + half_range
        
        # Extract variants in this region
        region_variants = data[
            (data['CHR'] == lead['CHR']) &
            (data['POS'] >= region_start) &
            (data['POS'] <= region_end)
        ].copy()
        
        region_variants['LEAD_VARIANT'] = lead['VARIANT_ID']
        region_variants['LEAD_CHR'] = lead['CHR']
        region_variants['LEAD_POS'] = lead['POS']
        region_variants['LEAD_P'] = lead['P']
        region_variants['PEAK_ID'] = idx + 1
        # Mark if this variant is the lead variant
        region_variants['IS_LEAD'] = region_variants['VARIANT_ID'] == lead['VARIANT_ID']
        
        peak_regions.append(region_variants)
        
        print(f"  Peak {idx+1}: chr{lead['CHR']}:{region_start:,}-{region_end:,} "
              f"(lead: {lead['VARIANT_ID']}, P={lead['P']:.2e}, n={len(region_variants):,} variants)")
    
    if peak_regions:
        all_regions = pd.concat(peak_regions, ignore_index=True)
        print(f"\nTotal variants extracted: {len(all_regions):,}")
        return all_regions
    else:
        return pd.DataFrame()


def save_results(lead_df, peak_regions, output_dir, prefix, args):
    """Save results to files"""
    os.makedirs(output_dir, exist_ok=True)
    
    # Create timestamp
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    
    # Save lead variants with metadata header
    lead_file = os.path.join(output_dir, f"{prefix}_lead_variants.tsv")
    with open(lead_file, 'w') as f:
        # Write metadata as comments
        f.write(f"# Lead Variant Extraction Results\n")
        f.write(f"# Generated: {timestamp}\n")
        f.write(f"# Association file: {args.assoc}\n")
        f.write(f"# Significance threshold: {args.sig_level:.2e}\n")
        f.write(f"# Lead range: {args.lead_range:,} bp (±{args.lead_range//2:,} bp)\n")
        f.write(f"# Window size: {args.window:,} bp (±{args.window//2:,} bp)\n")
        f.write(f"# Total lead variants: {len(lead_df)}\n")
        f.write(f"#\n")
    
    # Append data
    lead_df.to_csv(lead_file, sep='\t', index=False, float_format='%.6g', mode='a')
    print(f"\nLead variants saved to: {lead_file}")
    
    if len(peak_regions) > 0:
        # Reorder columns to put IS_LEAD and PEAK_ID first for better visibility
        cols = ['PEAK_ID', 'IS_LEAD', 'VARIANT_ID', 'CHR', 'POS', 'P', 
                'LEAD_VARIANT', 'LEAD_CHR', 'LEAD_POS', 'LEAD_P']
        peak_regions = peak_regions[cols]
        
        # Save all peak regions with metadata header
        regions_file = os.path.join(output_dir, f"{prefix}_peak_regions.tsv")
        with open(regions_file, 'w') as f:
            # Write metadata as comments
            f.write(f"# Peak Regions Extraction Results\n")
            f.write(f"# Generated: {timestamp}\n")
            f.write(f"# Association file: {args.assoc}\n")
            f.write(f"# Significance threshold: {args.sig_level:.2e}\n")
            f.write(f"# Lead range: {args.lead_range:,} bp (±{args.lead_range//2:,} bp)\n")
            f.write(f"# Window size: {args.window:,} bp (±{args.window//2:,} bp)\n")
            f.write(f"# Total peaks: {peak_regions['PEAK_ID'].nunique()}\n")
            f.write(f"# Total variants: {len(peak_regions):,}\n")
            f.write(f"#\n")
            f.write(f"# Column descriptions:\n")
            f.write(f"#   PEAK_ID: Peak number (independent locus)\n")
            f.write(f"#   IS_LEAD: TRUE if this variant is the lead variant for this peak\n")
            f.write(f"#   VARIANT_ID: Variant identifier\n")
            f.write(f"#   CHR: Chromosome\n")
            f.write(f"#   POS: Position (bp)\n")
            f.write(f"#   P: P-value\n")
            f.write(f"#   LEAD_VARIANT: Lead variant ID for this peak\n")
            f.write(f"#   LEAD_CHR: Lead variant chromosome\n")
            f.write(f"#   LEAD_POS: Lead variant position\n")
            f.write(f"#   LEAD_P: Lead variant P-value\n")
            f.write(f"#\n")
        
        # Append data
        peak_regions.to_csv(regions_file, sep='\t', index=False, float_format='%.6g', mode='a')
        print(f"Peak regions saved to: {regions_file}")
        
        # Save variant IDs only (one file per peak)
        for peak_id in peak_regions['PEAK_ID'].unique():
            peak_data = peak_regions[peak_regions['PEAK_ID'] == peak_id]
            lead_var = peak_data['LEAD_VARIANT'].iloc[0]
            
            ids_file = os.path.join(output_dir, f"{prefix}_peak{peak_id}_variants.txt")
            peak_data['VARIANT_ID'].to_csv(ids_file, index=False, header=False)
            print(f"Peak {peak_id} variant IDs saved to: {ids_file} (lead: {lead_var})")


def main():
    """Main execution function"""
    print("=" * 70)
    print("Lead Variant and Peak Region Extraction Tool")
    print("=" * 70)
    
    # Parse arguments
    args = parse_arguments()
    
    print(f"\nParameters:")
    print(f"  Association file:    {args.assoc}")
    print(f"  Significance level:  {args.sig_level:.2e}")
    print(f"  Lead range:          ±{args.lead_range//2:,} bp ({args.lead_range:,} bp total)")
    print(f"  Window size:         {args.window:,} bp (±{args.window//2:,} bp)")
    print(f"  Output directory:    {args.output_dir}")
    print(f"  Output prefix:       {args.prefix}")
    
    try:
        # Load data
        data = load_association_data(args.assoc)
        
        # Identify lead variants
        lead_df = identify_lead_variants(
            data, 
            args.sig_level, 
            args.window
        )
        
        if len(lead_df) == 0:
            print("\nNo significant peaks found. Exiting.")
            return 0
        
        # Extract peak regions
        peak_regions = extract_peak_regions(data, lead_df, args.lead_range)
        
        # Save results
        save_results(lead_df, peak_regions, args.output_dir, args.prefix, args)
        
        print("\n" + "=" * 70)
        print("Analysis completed successfully!")
        print("=" * 70)
        
        return 0
        
    except Exception as e:
        print(f"\nERROR: {str(e)}")
        import traceback
        traceback.print_exc()
        return 1


if __name__ == '__main__':
    sys.exit(main())