#!/usr/bin/env python3
"""
Find common samples and variants across WGS, Array Genotype, and Array Dosage datasets.

This script reads sample IDs and variant IDs from three different datasets:
- WGS: bed/bim/fam format
- Array GT: pgen/pvar/psam format (hard-called genotypes)
- Array DS: pgen/pvar/psam format (with dosage)

Outputs:
- common_samples.txt: FID IID format for PLINK2 --keep
- common_variants.txt: variant IDs for PLINK2 --extract
- intersection_summary.txt: detailed statistics
- *.sample_counts.txt: sample counts per dataset
- *.variant_counts.txt: variant counts per dataset
"""

import argparse
import sys
from pathlib import Path
from datetime import datetime


def read_fam_samples(fam_file):
    """Read sample IDs from .fam file (PLINK1 format)."""
    samples = set()
    with open(fam_file, 'r') as f:
        for line in f:
            fields = line.strip().split()
            if len(fields) >= 2:
                fid, iid = fields[0], fields[1]
                samples.add((fid, iid))
    return samples


def read_psam_samples(psam_file):
    """Read sample IDs from .psam file (PLINK2 format)."""
    samples = set()
    with open(psam_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) >= 2 and fields[0] != 'FID':
                fid, iid = fields[0], fields[1]
                samples.add((fid, iid))
    return samples


def read_bim_variants(bim_file):
    """Read variant IDs from .bim file (PLINK1 format)."""
    variants = set()
    with open(bim_file, 'r') as f:
        for line in f:
            fields = line.strip().split()
            if len(fields) >= 2:
                variant_id = fields[1]
                variants.add(variant_id)
    return variants


def read_pvar_variants(pvar_file):
    """Read variant IDs from .pvar file (PLINK2 format)."""
    variants = set()
    with open(pvar_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) >= 3:
                variant_id = fields[2]  # ID column
                variants.add(variant_id)
    return variants


def write_samples(samples, output_file):
    """Write samples in FID IID format for PLINK2 --keep."""
    with open(output_file, 'w') as f:
        for fid, iid in sorted(samples):
            f.write(f"{fid}\t{iid}\n")


def write_variants(variants, output_file):
    """Write variants for PLINK2 --extract."""
    with open(output_file, 'w') as f:
        for variant_id in sorted(variants):
            f.write(f"{variant_id}\n")


def main():
    parser = argparse.ArgumentParser(
        description='Find common samples and variants across WGS and Array datasets',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    parser.add_argument('--wgs-prefix', required=True,
                       help='Prefix for WGS files (bed/bim/fam)')
    parser.add_argument('--array-gt-prefix', required=True,
                       help='Prefix for Array Genotype files (pgen/pvar/psam)')
    parser.add_argument('--array-ds-prefix', required=True,
                       help='Prefix for Array Dosage files (pgen/pvar/psam)')
    parser.add_argument('--output-dir', required=True,
                       help='Output directory for results')
    parser.add_argument('--summary', default='intersection_summary.txt',
                       help='Summary statistics file')
    
    args = parser.parse_args()
    
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print("=" * 80)
    print("Finding Common Samples and Variants")
    print("=" * 80)
    print(f"Start time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print()
    
    # Read samples from each dataset
    print("Reading samples...")
    wgs_samples = read_fam_samples(f"{args.wgs_prefix}.fam")
    array_gt_samples = read_psam_samples(f"{args.array_gt_prefix}.psam")
    array_ds_samples = read_psam_samples(f"{args.array_ds_prefix}.psam")
    
    print(f"  WGS samples:             {len(wgs_samples):,}")
    print(f"  Array GT samples:        {len(array_gt_samples):,}")
    print(f"  Array DS samples:        {len(array_ds_samples):,}")
    print()
    
    # Write individual sample counts
    with open(output_dir / "wgs.sample_counts.txt", 'w') as f:
        f.write(f"WGS samples: {len(wgs_samples)}\n")
    with open(output_dir / "array_gt.sample_counts.txt", 'w') as f:
        f.write(f"Array GT samples: {len(array_gt_samples)}\n")
    with open(output_dir / "array_ds.sample_counts.txt", 'w') as f:
        f.write(f"Array DS samples: {len(array_ds_samples)}\n")
    
    # Read variants from each dataset
    print("Reading variants...")
    wgs_variants = read_bim_variants(f"{args.wgs_prefix}.bim")
    array_gt_variants = read_pvar_variants(f"{args.array_gt_prefix}.pvar")
    array_ds_variants = read_pvar_variants(f"{args.array_ds_prefix}.pvar")
    
    print(f"  WGS variants:            {len(wgs_variants):,}")
    print(f"  Array GT variants:       {len(array_gt_variants):,}")
    print(f"  Array DS variants:       {len(array_ds_variants):,}")
    print()
    
    # Write individual variant counts
    with open(output_dir / "wgs.variant_counts.txt", 'w') as f:
        f.write(f"WGS variants: {len(wgs_variants)}\n")
    with open(output_dir / "array_gt.variant_counts.txt", 'w') as f:
        f.write(f"Array GT variants: {len(array_gt_variants)}\n")
    with open(output_dir / "array_ds.variant_counts.txt", 'w') as f:
        f.write(f"Array DS variants: {len(array_ds_variants)}\n")
    
    # Find common samples (intersection of all three)
    print("Finding common samples...")
    common_samples = wgs_samples & array_gt_samples & array_ds_samples
    print(f"  Common samples:          {len(common_samples):,}")
    print(f"  Retention rate:          {len(common_samples)/len(wgs_samples)*100:.2f}% (vs WGS)")
    print()
    
    # Find common variants (intersection of all three)
    print("Finding common variants...")
    common_variants = wgs_variants & array_gt_variants & array_ds_variants
    print(f"  Common variants:         {len(common_variants):,}")
    print(f"  Retention rate:          {len(common_variants)/len(wgs_variants)*100:.2f}% (vs WGS)")
    print()
    
    # Check if we have any common samples/variants
    if len(common_samples) == 0:
        print("ERROR: No common samples found across all datasets!", file=sys.stderr)
        sys.exit(1)
    
    if len(common_variants) == 0:
        print("ERROR: No common variants found across all datasets!", file=sys.stderr)
        sys.exit(1)
    
    # Calculate pairwise intersections
    wgs_array_gt_samples = wgs_samples & array_gt_samples
    wgs_array_ds_samples = wgs_samples & array_ds_samples
    array_gt_ds_samples = array_gt_samples & array_ds_samples
    
    wgs_array_gt_variants = wgs_variants & array_gt_variants
    wgs_array_ds_variants = wgs_variants & array_ds_variants
    array_gt_ds_variants = array_gt_variants & array_ds_variants
    
    # Write common samples and variants
    print("Writing output files...")
    write_samples(common_samples, output_dir / "common_samples.txt")
    write_variants(common_variants, output_dir / "common_variants.txt")
    print(f"  ✓ {output_dir / 'common_samples.txt'}")
    print(f"  ✓ {output_dir / 'common_variants.txt'}")
    print()
    
    # Write comprehensive summary report
    summary_file = output_dir / args.summary
    with open(summary_file, 'w') as f:
        f.write("="*100 + "\n")
        f.write(" "*20 + "DATASET INTERSECTION ANALYSIS REPORT\n")
        f.write(" "*10 + "WGS vs Imputed Array Genotype (GT) vs Imputed Array Dosage (DS)\n")
        f.write("="*100 + "\n\n")
        
        f.write(f"Report generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"Analysis performed by: find_common_samples_variants.py\n\n")
        
        # Section 1: Input Datasets
        f.write("-"*100 + "\n")
        f.write("1. INPUT DATASETS\n")
        f.write("-"*100 + "\n\n")
        f.write(f"  [WGS]              {args.wgs_prefix}\n")
        f.write(f"                     Format: PLINK1 binary (bed/bim/fam)\n")
        f.write(f"                     Data type: Whole Genome Sequencing\n\n")
        f.write(f"  [Imputed Array GT] {args.array_gt_prefix}\n")
        f.write(f"                     Format: PLINK2 (pgen/pvar/psam)\n")
        f.write(f"                     Data type: Imputed array - Hard-called genotypes (0/1/2)\n\n")
        f.write(f"  [Imputed Array DS] {args.array_ds_prefix}\n")
        f.write(f"                     Format: PLINK2 (pgen/pvar/psam)\n")
        f.write(f"                     Data type: Imputed array - Dosage information (0.0-2.0)\n\n")
        
        # Section 2: Sample Statistics
        f.write("-"*100 + "\n")
        f.write("2. SAMPLE STATISTICS\n")
        f.write("-"*100 + "\n\n")
        f.write(f"  {'Dataset':<20} {'Total Samples':>15} {'Overlap with WGS':>20} {'Overlap %':>12}\n")
        f.write(f"  {'-'*20} {'-'*15} {'-'*20} {'-'*12}\n")
        f.write(f"  {'WGS':<20} {len(wgs_samples):>15,} {len(wgs_samples):>20,} {'100.00%':>12}\n")
        f.write(f"  {'Imputed Array GT':<20} {len(array_gt_samples):>15,} {len(wgs_array_gt_samples):>20,} {len(wgs_array_gt_samples)/len(wgs_samples)*100:>11.2f}%\n")
        f.write(f"  {'Imputed Array DS':<20} {len(array_ds_samples):>15,} {len(wgs_array_ds_samples):>20,} {len(wgs_array_ds_samples)/len(wgs_samples)*100:>11.2f}%\n")
        f.write("\n")
        f.write(f"  {'Three-way intersection:':<50} {len(common_samples):>15,}\n")
        f.write(f"  {'Retention rate (vs WGS baseline):':<50} {len(common_samples)/len(wgs_samples)*100:>14.2f}%\n\n")
        
        # Venn diagram style for samples
        f.write("  Pairwise Sample Overlaps:\n")
        f.write(f"    WGS ∩ Imputed Array GT:      {len(wgs_array_gt_samples):>10,} samples\n")
        f.write(f"    WGS ∩ Imputed Array DS:      {len(wgs_array_ds_samples):>10,} samples\n")
        f.write(f"    Array GT ∩ Array DS:         {len(array_gt_ds_samples):>10,} samples\n")
        f.write(f"    All three (final common):    {len(common_samples):>10,} samples\n\n")
        
        # Section 3: Variant Statistics
        f.write("-"*100 + "\n")
        f.write("3. VARIANT STATISTICS\n")
        f.write("-"*100 + "\n\n")
        f.write(f"  {'Dataset':<20} {'Total Variants':>15} {'Overlap with WGS':>20} {'Overlap %':>12}\n")
        f.write(f"  {'-'*20} {'-'*15} {'-'*20} {'-'*12}\n")
        f.write(f"  {'WGS':<20} {len(wgs_variants):>15,} {len(wgs_variants):>20,} {'100.00%':>12}\n")
        f.write(f"  {'Imputed Array GT':<20} {len(array_gt_variants):>15,} {len(wgs_array_gt_variants):>20,} {len(wgs_array_gt_variants)/len(wgs_variants)*100:>11.2f}%\n")
        f.write(f"  {'Imputed Array DS':<20} {len(array_ds_variants):>15,} {len(wgs_array_ds_variants):>20,} {len(wgs_array_ds_variants)/len(wgs_variants)*100:>11.2f}%\n")
        f.write("\n")
        f.write(f"  {'Three-way intersection:':<50} {len(common_variants):>15,}\n")
        f.write(f"  {'Retention rate (vs WGS baseline):':<50} {len(common_variants)/len(wgs_variants)*100:>14.2f}%\n\n")
        
        # Venn diagram style for variants
        f.write("  Pairwise Variant Overlaps:\n")
        f.write(f"    WGS ∩ Imputed Array GT:      {len(wgs_array_gt_variants):>10,} variants\n")
        f.write(f"    WGS ∩ Imputed Array DS:      {len(wgs_array_ds_variants):>10,} variants\n")
        f.write(f"    Array GT ∩ Array DS:         {len(array_gt_ds_variants):>10,} variants\n")
        f.write(f"    All three (final common):    {len(common_variants):>10,} variants\n\n")
        
        # Section 4: Output Files
        f.write("-"*100 + "\n")
        f.write("4. OUTPUT FILES\n")
        f.write("-"*100 + "\n\n")
        f.write(f"  Common samples list:   {output_dir / 'common_samples.txt'}\n")
        f.write(f"                         Format: FID TAB IID (compatible with PLINK2 --keep)\n")
        f.write(f"                         Records: {len(common_samples):,}\n\n")
        f.write(f"  Common variants list:  {output_dir / 'common_variants.txt'}\n")
        f.write(f"                         Format: Variant ID per line (compatible with PLINK2 --extract)\n")
        f.write(f"                         Records: {len(common_variants):,}\n\n")
        
        # Section 5: Quality Metrics
        f.write("-"*100 + "\n")
        f.write("5. DATA QUALITY METRICS\n")
        f.write("-"*100 + "\n\n")
        
        sample_concordance = len(common_samples) / min(len(wgs_samples), len(array_gt_samples), len(array_ds_samples)) * 100
        variant_concordance = len(common_variants) / min(len(wgs_variants), len(array_gt_variants), len(array_ds_variants)) * 100
        
        f.write(f"  Sample Concordance:    {sample_concordance:>6.2f}%  (common / smallest dataset)\n")
        f.write(f"  Variant Concordance:   {variant_concordance:>6.2f}%  (common / smallest dataset)\n\n")
        
        # Interpretation
        f.write("  Interpretation:\n")
        if sample_concordance > 95:
            f.write("    ✓ Excellent sample overlap across all three datasets\n")
        elif sample_concordance > 85:
            f.write("    ⚠ Good sample overlap, minor discrepancies detected\n")
        else:
            f.write("    ⚠ Moderate sample overlap, investigate discrepancies\n")
            
        if variant_concordance > 90:
            f.write("    ✓ Excellent variant overlap across all three datasets\n")
        elif variant_concordance > 75:
            f.write("    ⚠ Good variant overlap, some variants unique to individual datasets\n")
        else:
            f.write("    ⚠ Moderate variant overlap, substantial dataset-specific variants present\n")
        
        f.write("\n")
        f.write("="*100 + "\n")
        f.write(" "*40 + "END OF REPORT\n")
        f.write("="*100 + "\n")
    
    print(f"  ✓ {summary_file}")
    print()
    print("=" * 80)
    print("✓ Analysis complete!")
    print("=" * 80)


if __name__ == '__main__':
    main()
