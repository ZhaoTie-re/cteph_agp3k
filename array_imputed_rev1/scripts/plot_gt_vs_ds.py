#!/usr/bin/env python3
import argparse
import os
import sys
import subprocess
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

def setup_plotting_style():
    # Use a clean, professional style manually since specific styles might not be installed
    plt.rcParams['figure.figsize'] = (10, 8)
    plt.rcParams['font.size'] = 12
    plt.rcParams['axes.grid'] = True
    plt.rcParams['axes.grid.axis'] = 'both'
    plt.rcParams['grid.alpha'] = 0.5
    plt.rcParams['axes.axisbelow'] = True # Grid behind plot elements

def parse_vcf_line(line):
    parts = line.strip().split('\t')
    chrom = parts[0]
    pos = parts[1]
    var_id = parts[2]
    ref = parts[3]
    alt = parts[4]
    info = parts[7]
    format_str = parts[8]
    samples = parts[9:]
    
    # Check what formats are available
    format_keys = format_str.split(':')
    try:
        gt_idx = format_keys.index('GT')
    except ValueError:
        return None # No GT?
    
    ds_idx = format_keys.index('DS') if 'DS' in format_keys else -1
    gp_idx = format_keys.index('GP') if 'GP' in format_keys else -1
    hds_idx = format_keys.index('HDS') if 'HDS' in format_keys else -1

    gt_values = []
    ds_values = []
    max_gp_values = [] 
    hds_consistent = [] 
    
    for sample in samples:
        sample_parts = sample.split(':')
        
        # Get GT
        if gt_idx < len(sample_parts):
            gt_str = sample_parts[gt_idx]
            if '.' in gt_str or gt_str == './.':
                gt_val = np.nan
            else:
                # Handle | and /
                alleles = gt_str.replace('|', '/').split('/')
                try:
                    gt_val = sum(int(a) for a in alleles)
                except ValueError:
                    gt_val = np.nan
        else:
            gt_val = np.nan
        
        # Get DS
        ds_val = np.nan
        if ds_idx != -1 and ds_idx < len(sample_parts):
            try:
                ds_val = float(sample_parts[ds_idx])
            except (ValueError, IndexError):
                pass

        # Get HDS
        hds_ok = True
        if hds_idx != -1 and hds_idx < len(sample_parts) and not np.isnan(ds_val):
            try:
                hds_str = sample_parts[hds_idx]
                hds_vals = [float(x) for x in hds_str.split(',')]
                if abs(sum(hds_vals) - ds_val) > 0.01:
                    hds_ok = False
            except (ValueError, IndexError):
                pass
            
        # Get GP
        gp_max = np.nan
        if gp_idx != -1 and gp_idx < len(sample_parts):
            try:
                gp_str = sample_parts[gp_idx]
                gps = [float(x) for x in gp_str.split(',')]
                if len(gps) == 3:
                    gp_max = max(gps)
            except (ValueError, IndexError):
                pass
                
        gt_values.append(gt_val)
        ds_values.append(ds_val)
        max_gp_values.append(gp_max)
        hds_consistent.append(hds_ok)
        
    return {
        'chrom': chrom,
        'pos': pos,
        'id': var_id,
        'ref': ref,
        'alt': alt,
        'GT': gt_values,
        'DS': ds_values,
        'MaxGP': max_gp_values,
        'HDS_OK': hds_consistent
    }

def process_variant(variant_str, vcf_dir, output_dir, tabix_path):
    try:
        parts = variant_str.split(':')
        chrom = parts[0]
        pos = parts[1]
        ref = parts[2]
        alt = parts[3]
    except (ValueError, IndexError):
        print(f"Skipping malformed variant string: {variant_str}")
        return

    vcf_file = os.path.join(vcf_dir, f"{chrom}.normalized.vcf.gz")
    if not os.path.exists(vcf_file):
        print(f"VCF file not found: {vcf_file}")
        return

    region = f"{chrom}:{pos}-{pos}"
    cmd = [tabix_path, vcf_file, region]
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        lines = result.stdout.strip().split('\n')
    except subprocess.CalledProcessError as e:
        print(f"Error running tabix for {variant_str}: {e}")
        return

    found = False
    for line in lines:
        if not line: continue
        data = parse_vcf_line(line)
        if not data: continue
        
        # Verify exact match
        if data['chrom'] == chrom and str(data['pos']) == str(pos) and data['ref'] == ref and data['alt'] == alt:
            found = True
            plot_gt_ds(data, output_dir, variant_str)
            break
            
    if not found:
        print(f"Variant {variant_str} not found")

def plot_gt_ds(data, output_dir, variant_str):
    df = pd.DataFrame({
        'GT': data['GT'], 
        'DS': data['DS'],
        'MaxGP': data['MaxGP'],
        'HDS_OK': data['HDS_OK']
    })
    
    df = df.dropna(subset=['DS'])
    
    # Use a professional layout with extra space on right for legends
    fig, ax = plt.subplots(figsize=(14, 8))
    
    # Adjust layout to make room on the right
    plt.subplots_adjust(right=0.75)
    
    # Jitter to avoid overlapping points
    jitter_strength = 0.08
    if 'GT' in df.columns and df['GT'].notna().any():
        x_data = df['GT'] + np.random.normal(0, jitter_strength, size=len(df))
        
        # === Background Zones (PLINK2 Exclusion) ===
        # Exclude regions where |DS - Integer| > 0.1
        ax.axhspan(0.1, 0.9, color='#e74c3c', alpha=0.1, zorder=0, label='PLINK2 Excluded Zone (|DS-k|>0.1)')
        ax.axhspan(1.1, 1.9, color='#e74c3c', alpha=0.1, zorder=0)

        # === Scatter with Professional Color Scale ===
        # Use RdYlGn (Red=Low Conf, Green=High Conf)
        if df['MaxGP'].notna().any():
            # Ensure the color scale emphasizes the high end (0.9-1.0)
            sc = ax.scatter(x_data, df['DS'], c=df['MaxGP'], cmap='RdYlGn', 
                             alpha=0.85, s=30, edgecolor='black', linewidth=0.3, zorder=3, 
                             vmin=0.5, vmax=1.0)
            
            # Place colorbar outside to the right, strictly in bottom half
            cax = fig.add_axes([0.78, 0.12, 0.02, 0.35]) # [left, bottom, width, height]
            cbar = plt.colorbar(sc, cax=cax)
            cbar.set_label('Genotype Confidence (Max GP)', fontsize=11, weight='bold')
            cbar.ax.tick_params(labelsize=9)
        else:
            ax.scatter(x_data, df['DS'], alpha=0.5, color='royalblue', s=30, edgecolor='black', linewidth=0.3, zorder=3)
            
        # === Formatting ===
        ax.set_xlabel("Genotype (GT) Hard Call", fontsize=13, weight='bold')
        ax.set_xticks([0, 1, 2])
        ax.set_xticklabels(['0/0 (HomRef)', '0/1 (Het)', '1/1 (HomAlt)'], fontsize=11)
        
        ax.set_ylabel("Dosage (DS)", fontsize=13, weight='bold')
        ax.set_yticks(np.arange(0, 2.25, 0.25))
        ax.set_yticklabels([f"{x:.2f}" for x in np.arange(0, 2.25, 0.25)], fontsize=11)
        
        ax.set_title(f"Genotype Quality & Dosage Discordance Analysis\n{variant_str}", fontsize=15, weight='bold', pad=20)
        
        # Perfect Concordance Line
        ax.plot([-0.3, 2.3], [-0.3, 2.3], 'k--', alpha=0.5, lw=1.5, zorder=2, label='Perfect Concordance (DS=GT)')
        
        ax.set_xlim(-0.5, 2.5)
        ax.set_ylim(-0.1, 2.1)
        
        # === Stats Box ===
        stats_text = []
        for gt in [0, 1, 2]:
            sub = df[df['GT'] == gt]
            if len(sub) > 0:
                # Count PLINK excluded
                excluded = sub[abs(sub['DS'] - gt) > 0.1]
                stats_text.append(f"GT {gt}: Total={len(sub)}, Rej={len(excluded)} ({len(excluded)/len(sub)*100:.1f}%)")
        
        hds_inconsistent = len(df[~df['HDS_OK']])
        if hds_inconsistent > 0:
            stats_text.append(f"HDS Error: {hds_inconsistent}")
            
        text_str = '\n'.join(stats_text)
        
        # Place Stats Box to the right, outside the plot
        # Move down to y=0.70 to avoid overlap with upper Legend
        props = dict(boxstyle='round,pad=0.6', facecolor='white', alpha=1.0, edgecolor='#cccccc')
        fig.text(0.77, 0.70, text_str, 
                 verticalalignment='top', fontsize=10, fontfamily='monospace',
                 bbox=props)

        # Legend strictly at top right
        ax.legend(loc='upper left', bbox_to_anchor=(1.02, 1), borderaxespad=0, framealpha=0.95)

    # plt.tight_layout() # Conflict with manual add_axes and subplots_adjust
    
    safe_var_name = variant_str.replace(':', '_')
    out_file = os.path.join(output_dir, f"{safe_var_name}.GT_vs_DS_Analysis.png")
    plt.savefig(out_file, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Plot saved to {out_file}")

def main():
    parser = argparse.ArgumentParser(description="Visualize GT vs DS for specific variants")
    parser.add_argument("--variant-list", required=True, help="File containing list of variants (chr:pos:ref:alt)")
    parser.add_argument("--vcf-dir", required=True, help="Directory containing normalized VCFs")
    parser.add_argument("--output-dir", required=True, help="Output directory for plots")
    parser.add_argument("--tabix", default="/usr/bin/tabix", help="Path to tabix")
    
    args = parser.parse_args()
    
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)
        
    with open(args.variant_list, 'r') as f:
        variants = [line.strip() for line in f if line.strip()]
        
    for var in variants:
        print(f"Processing {var}...")
        process_variant(var, args.vcf_dir, args.output_dir, args.tabix)

if __name__ == "__main__":
    setup_plotting_style()
    main()
