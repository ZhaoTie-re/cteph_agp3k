#!/usr/bin/env python3
"""
Merge statistics from all chromosomes into a single summary file
"""

import sys
import os
import argparse
from pathlib import Path

# Add script directory to path for imports
sys.path.insert(0, str(Path(__file__).parent))
from logger_utils import setup_logging

def parse_stats_file(file_path):
    """Parse a single chromosome stats file"""
    with open(file_path, 'r') as f:
        lines = f.readlines()
        
        # Skip comment lines and find data line
        data_line = None
        for line in lines:
            if not line.startswith('#') and not line.startswith('Chromosome'):
                data_line = line.strip()
                break
        
        if not data_line:
            return None
        
        # Parse data (format: chr, total_raw, invalid, invalid_pct, total_valid, valid_pct, mq_pass, mq_pass_pct, vqslod_pass, vqslod_pass_pct, both_pass, both_pass_pct)
        data = data_line.split('\t')
        if len(data) >= 12:
            return {
                'chromosome': data[0],
                'total_raw': data[1],
                'invalid': data[2],
                'total_valid': data[4],
                'mq_pass': data[6],
                'vqslod_pass': data[8],
                'both_pass': data[10]
            }
    return None

def parse_number(num_str):
    """Parse comma-formatted number string to integer"""
    return int(num_str.replace(',', ''))

def format_number(num):
    """Format number with comma as thousands separator"""
    return f"{num:,}"

def main():
    parser = argparse.ArgumentParser(
        description='Merge variant statistics from all chromosomes into a summary file',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example:
  %(prog)s all_summary.txt chr1.stats.txt chr2.stats.txt ... chr22.stats.txt
        """
    )
    
    parser.add_argument('output_file', 
                       help='Output merged summary file path')
    parser.add_argument('stats_files', nargs='+',
                       help='Input statistics files from individual chromosomes')
    
    args = parser.parse_args()
    
    # Setup dual logging (stdout + file)
    log_file = setup_logging(args.output_file)
    
    print(f"Log file: {log_file}")
    print(f"Merging {len(args.stats_files)} statistics files...")
    
    # Collect all chromosome stats
    all_stats = []
    total_raw = 0
    total_invalid = 0
    total_valid = 0
    total_mq_pass = 0
    total_vqslod_pass = 0
    total_both_pass = 0
    
    for stats_file in sorted(args.stats_files):
        if not os.path.exists(stats_file):
            print(f"Warning: {stats_file} does not exist, skipping...")
            continue
        
        stats = parse_stats_file(stats_file)
        if stats:
            all_stats.append(stats)
            
            # Add to totals
            total_raw += parse_number(stats['total_raw'])
            total_invalid += parse_number(stats['invalid'])
            total_valid += parse_number(stats['total_valid'])
            total_mq_pass += parse_number(stats['mq_pass'])
            total_vqslod_pass += parse_number(stats['vqslod_pass'])
            total_both_pass += parse_number(stats['both_pass'])
            
            print(f"  Added {stats['chromosome']}: {stats['total_raw']} variants ({stats['invalid']} invalid)")
    
    # Sort chromosomes naturally: chr1-22, PAR, chrX, chrY
    def chr_sort_key(stat):
        chrom = stat['chromosome']
        if chrom == 'PAR':
            return (23, 0)
        elif chrom.startswith('chr'):
            chrom_num = chrom[3:]
            if chrom_num == 'X':
                return (24, 0)
            elif chrom_num == 'Y':
                return (25, 0)
            else:
                try:
                    return (int(chrom_num), 0)
                except ValueError:
                    return (98, chrom_num)
        return (100, chrom)
    
    all_stats.sort(key=chr_sort_key)
    
    # Write merged results
    with open(args.output_file, 'w') as out:
        # Write header with metadata
        out.write("#" + "="*100 + "\n")
        out.write("# Merged Variant Quality Control Statistics\n")
        out.write("#" + "="*100 + "\n")
        out.write(f"# Generated: {__import__('time').strftime('%Y-%m-%d %H:%M:%S')}\n")
        out.write(f"# Total Chromosomes: {len(all_stats)}\n")
        out.write("#\n")
        out.write("# Column Definitions:\n")
        out.write("#   Total_Raw     : Total number of variants in VCF file\n")
        out.write("#   Invalid       : Variants with inf/nan values in MQ or VQSLOD\n")
        out.write("#   Invalid_Pct   : Percentage of invalid variants (Invalid/Total_Raw * 100)\n")
        out.write("#   Total_Valid   : Variants with finite MQ and VQSLOD values\n")
        out.write("#   Valid_Pct     : Percentage of valid variants (Total_Valid/Total_Raw * 100)\n")
        out.write("#   MQ_Pass       : Valid variants with MQ > threshold\n")
        out.write("#   MQ_Pass_Pct   : Percentage passing MQ (MQ_Pass/Total_Raw * 100)\n")
        out.write("#   VQSLOD_Pass   : Valid variants with VQSLOD > threshold\n")
        out.write("#   VQSLOD_Pass_Pct: Percentage passing VQSLOD (VQSLOD_Pass/Total_Raw * 100)\n")
        out.write("#   Both_Pass     : Valid variants passing both MQ and VQSLOD thresholds\n")
        out.write("#   Both_Pass_Pct : Percentage passing both (Both_Pass/Total_Raw * 100)\n")
        out.write("#\n")
        out.write("# Note: All percentages are calculated based on Total_Raw variants\n")
        out.write("#" + "="*100 + "\n")
        out.write("\n")
        
        # Column header
        out.write("Chromosome\tTotal_Raw\tInvalid\tInvalid_Pct\tTotal_Valid\tValid_Pct\tMQ_Pass\tMQ_Pass_Pct\tVQSLOD_Pass\tVQSLOD_Pass_Pct\tBoth_Pass\tBoth_Pass_Pct\n")
        
        # Write each chromosome with percentages
        for stats in all_stats:
            raw = parse_number(stats['total_raw'])
            invalid = parse_number(stats['invalid'])
            valid = parse_number(stats['total_valid'])
            mq = parse_number(stats['mq_pass'])
            vqslod = parse_number(stats['vqslod_pass'])
            both = parse_number(stats['both_pass'])
            
            # Calculate percentages based on total_raw
            invalid_pct = (invalid / raw * 100) if raw > 0 else 0
            valid_pct = (valid / raw * 100) if raw > 0 else 0
            mq_pct = (mq / raw * 100) if raw > 0 else 0
            vqslod_pct = (vqslod / raw * 100) if raw > 0 else 0
            both_pct = (both / raw * 100) if raw > 0 else 0
            
            out.write(f"{stats['chromosome']}\t{stats['total_raw']}\t")
            out.write(f"{stats['invalid']}\t{invalid_pct:.2f}%\t")
            out.write(f"{stats['total_valid']}\t{valid_pct:.2f}%\t")
            out.write(f"{stats['mq_pass']}\t{mq_pct:.2f}%\t")
            out.write(f"{stats['vqslod_pass']}\t{vqslod_pct:.2f}%\t")
            out.write(f"{stats['both_pass']}\t{both_pct:.2f}%\n")
        
        # Write separator and totals
        out.write("\n")
        out.write("#" + "-"*100 + "\n")
        out.write("# SUMMARY STATISTICS\n")
        out.write("#" + "-"*100 + "\n")
        
        # Calculate total percentages based on total_raw
        invalid_pct = (total_invalid / total_raw * 100) if total_raw > 0 else 0
        valid_pct = (total_valid / total_raw * 100) if total_raw > 0 else 0
        mq_pct = (total_mq_pass / total_raw * 100) if total_raw > 0 else 0
        vqslod_pct = (total_vqslod_pass / total_raw * 100) if total_raw > 0 else 0
        both_pct = (total_both_pass / total_raw * 100) if total_raw > 0 else 0
        
        out.write(f"TOTAL\t{format_number(total_raw)}\t")
        out.write(f"{format_number(total_invalid)}\t{invalid_pct:.2f}%\t")
        out.write(f"{format_number(total_valid)}\t{valid_pct:.2f}%\t")
        out.write(f"{format_number(total_mq_pass)}\t{mq_pct:.2f}%\t")
        out.write(f"{format_number(total_vqslod_pass)}\t{vqslod_pct:.2f}%\t")
        out.write(f"{format_number(total_both_pass)}\t{both_pct:.2f}%\n")
        
        # Write summary notes
        out.write("\n")
        out.write("#" + "="*100 + "\n")
        out.write("# INTERPRETATION NOTES\n")
        out.write("#" + "="*100 + "\n")
        if total_raw > 0:
            out.write(f"# Total variants across all chromosomes: {format_number(total_raw)}\n")
            out.write(f"# Data quality: {invalid_pct:.2f}% invalid, {valid_pct:.2f}% valid\n")
            out.write(f"# Quality filtering results (% of total raw variants):\n")
            out.write(f"#   - MQ filter pass rate: {mq_pct:.2f}%\n")
            out.write(f"#   - VQSLOD filter pass rate: {vqslod_pct:.2f}%\n")
            out.write(f"#   - Combined filters pass rate: {both_pct:.2f}%\n")
            out.write("#" + "="*100 + "\n")
    
    print(f"\nMerged statistics saved to: {args.output_file}")
    print(f"\nSummary:")
    print(f"  Total raw variants: {format_number(total_raw)}")
    
    if total_raw > 0:
        print(f"  Invalid variants: {format_number(total_invalid)} ({total_invalid/total_raw*100:.2f}%)")
        print(f"  Valid variants: {format_number(total_valid)} ({total_valid/total_raw*100:.2f}%)")
        print(f"\n  Quality filtering pass rates (% of total raw):")
        print(f"    MQ pass: {format_number(total_mq_pass)} ({total_mq_pass/total_raw*100:.2f}%)")
        print(f"    VQSLOD pass: {format_number(total_vqslod_pass)} ({total_vqslod_pass/total_raw*100:.2f}%)")
        print(f"    Both pass: {format_number(total_both_pass)} ({total_both_pass/total_raw*100:.2f}%)")
    else:
        print(f"  Warning: No variants found across all chromosomes")

if __name__ == "__main__":
    main()
