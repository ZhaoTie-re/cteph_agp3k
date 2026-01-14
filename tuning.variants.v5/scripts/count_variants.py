#!/usr/bin/env python3
"""
Count variants passing different thresholds from extracted TSV files
Reads MQ and VQSLOD values from TSV format for fast in-memory counting
Uses chunked processing with numpy for memory efficiency
"""

import sys
import argparse
import os
import time
import numpy
from pathlib import Path

# Add script directory to path for imports
sys.path.insert(0, str(Path(__file__).parent))
from logger_utils import setup_logging

def format_number(num):
    """Format number with comma as thousands separator"""
    return f"{num:,}"

def count_variants_from_tsv(tsv_file, mq_threshold, vqslod_threshold, chunk_size=500000):
    """
    Count variants from TSV file using chunked processing with numpy
    for better performance and memory efficiency
    Chunk size of 500K provides ~10 chunks for 5M variants per chromosome
    Returns: (total_raw, total_invalid, total_valid, mq_pass, vqslod_pass, both_pass)
    """
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] Reading TSV file: {tsv_file}")
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] Using chunked processing (chunk size: {format_number(chunk_size)})")
    
    total_raw = 0
    total_invalid = 0
    total_valid = 0
    mq_pass = 0
    vqslod_pass = 0
    both_pass = 0
    value_error_count = 0
    
    start_time = time.time()
    
    with open(tsv_file, 'r') as f:
        # Skip header
        header = f.readline().strip()
        if header != "MQ\tVQSLOD":
            print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] Warning: Unexpected header: {header}")
            print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}]   Expected: MQ\\tVQSLOD")
        
        print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] Counting variants in chunks...")
        
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
                
                # Handle potential multi-value fields (comma separated) by taking the MAX value
                # This aligns with bcftools logic: "expression is true if any of the values satisfies the condition"
                if ',' in mq_str:
                    # Filter out non-numeric values (like '.') before taking max
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
                    # Filter out non-numeric values (like '.') before taking max
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
                value_error_count += 1
                continue
            
            # Process chunk when buffer is full
            if len(mq_buffer) >= chunk_size:
                # Use numpy for vectorized operations (much faster)
                mq_arr = numpy.array(mq_buffer, dtype=numpy.float64)
                vqslod_arr = numpy.array(vqslod_buffer, dtype=numpy.float64)
                
                chunk_raw = len(mq_arr)
                total_raw += chunk_raw
                
                # Filter out invalid values (inf, -inf, nan)
                valid_mask = numpy.isfinite(mq_arr) & numpy.isfinite(vqslod_arr)
                chunk_invalid = numpy.sum(~valid_mask)
                total_invalid += chunk_invalid
                
                mq_arr = mq_arr[valid_mask]
                vqslod_arr = vqslod_arr[valid_mask]
                
                chunk_valid = len(mq_arr)
                total_valid += chunk_valid
                chunk_mq_pass = numpy.sum(mq_arr > mq_threshold)
                chunk_vqslod_pass = numpy.sum(vqslod_arr > vqslod_threshold)
                chunk_both_pass = numpy.sum((mq_arr > mq_threshold) & (vqslod_arr > vqslod_threshold))
                
                mq_pass += chunk_mq_pass
                vqslod_pass += chunk_vqslod_pass
                both_pass += chunk_both_pass
                
                elapsed = time.time() - start_time
                rate = total_raw / elapsed
                print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}]   Processed {format_number(total_raw)} variants ({rate:.0f} var/s)")
                
                mq_buffer = []
                vqslod_buffer = []
        
        # Process remaining data
        if mq_buffer:
            mq_arr = numpy.array(mq_buffer, dtype=numpy.float64)
            vqslod_arr = numpy.array(vqslod_buffer, dtype=numpy.float64)
            
            chunk_raw = len(mq_arr)
            total_raw += chunk_raw
            
            # Filter out invalid values (inf, -inf, nan)
            valid_mask = numpy.isfinite(mq_arr) & numpy.isfinite(vqslod_arr)
            chunk_invalid = numpy.sum(~valid_mask)
            total_invalid += chunk_invalid
            
            if chunk_invalid > 0:
                print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}]   Found {format_number(chunk_invalid)} invalid values (inf/nan) in final chunk")
            
            mq_arr = mq_arr[valid_mask]
            vqslod_arr = vqslod_arr[valid_mask]
            
            chunk_valid = len(mq_arr)
            total_valid += chunk_valid
            chunk_mq_pass = numpy.sum(mq_arr > mq_threshold)
            chunk_vqslod_pass = numpy.sum(vqslod_arr > vqslod_threshold)
            chunk_both_pass = numpy.sum((mq_arr > mq_threshold) & (vqslod_arr > vqslod_threshold))
            
            mq_pass += int(chunk_mq_pass)
            vqslod_pass += int(chunk_vqslod_pass)
            both_pass += int(chunk_both_pass)
    
    # Add skipped lines to totals
    if value_error_count > 0:
        print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}]   Found {format_number(value_error_count)} lines with parsing errors (ValueError)")
        total_raw += value_error_count
        total_invalid += value_error_count

    elapsed = time.time() - start_time
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] Counting completed in {elapsed:.2f}s")
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}]   Total raw variants: {format_number(total_raw)}")
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}]   Invalid variants: {format_number(total_invalid)} ({total_invalid/total_raw*100:.2f}%)" if total_raw > 0 else "")
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}]   Valid variants: {format_number(total_valid)}")
    
    # Handle edge case where no variants found
    if total_valid == 0:
        print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] Warning: No valid variants found in TSV file")
    
    return total_raw, total_invalid, total_valid, mq_pass, vqslod_pass, both_pass

def main():
    parser = argparse.ArgumentParser(
        description='Count variants passing MQ and VQSLOD thresholds from TSV files',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example:
  %(prog)s chr1.metrics.tsv chr1 58.75 10 chr1.stats.txt
        """
    )
    
    parser.add_argument('tsv_file', 
                       help='Input TSV file with MQ and VQSLOD values')
    parser.add_argument('chromosome', 
                       help='Chromosome name for labeling (e.g., chr1, chr22, chrX)')
    parser.add_argument('mq_threshold', type=float,
                       help='Mapping Quality (MQ) threshold (variants with MQ > threshold pass)')
    parser.add_argument('vqslod_threshold', type=float,
                       help='VQSLOD score threshold (variants with VQSLOD > threshold pass)')
    parser.add_argument('output_file', 
                       help='Output statistics file path')
    
    args = parser.parse_args()
    
    # Setup dual logging (stdout + file)
    log_file = setup_logging(args.output_file)
    
    print(f"\n{'='*70}")
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] Starting variant counting for {args.chromosome}")
    print(f"{'='*70}")
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] Input TSV: {args.tsv_file}")
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] MQ threshold: {args.mq_threshold}, VQSLOD threshold: {args.vqslod_threshold}")
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] Log file: {log_file}")
    print(f"{'='*70}\n")
    
    # Count variants from TSV
    total_raw, total_invalid, total_valid, mq_pass, vqslod_pass, both_pass = count_variants_from_tsv(
        args.tsv_file, args.mq_threshold, args.vqslod_threshold
    )
    
    # Write results
    print(f"\n[{time.strftime('%Y-%m-%d %H:%M:%S')}] Writing results to {args.output_file}...")
    with open(args.output_file, 'w') as out:
        # Write header with column descriptions
        out.write(f"# Variant Quality Control Statistics\n")
        out.write(f"# MQ Threshold: {args.mq_threshold}, VQSLOD Threshold: {args.vqslod_threshold}\n")
        out.write(f"# All percentages are calculated based on Total_Raw variants\n")
        out.write(f"#\n")
        out.write(f"Chromosome\tTotal_Raw\tInvalid\tInvalid_Pct\tTotal_Valid\tValid_Pct\tMQ_Pass\tMQ_Pass_Pct\tVQSLOD_Pass\tVQSLOD_Pass_Pct\tBoth_Pass\tBoth_Pass_Pct\n")
        
        # Calculate percentages based on total_raw
        invalid_pct = (total_invalid / total_raw * 100) if total_raw > 0 else 0
        valid_pct = (total_valid / total_raw * 100) if total_raw > 0 else 0
        mq_pass_pct = (mq_pass / total_raw * 100) if total_raw > 0 else 0
        vqslod_pass_pct = (vqslod_pass / total_raw * 100) if total_raw > 0 else 0
        both_pass_pct = (both_pass / total_raw * 100) if total_raw > 0 else 0
        
        out.write(f"{args.chromosome}\t{format_number(total_raw)}\t")
        out.write(f"{format_number(total_invalid)}\t{invalid_pct:.2f}%\t")
        out.write(f"{format_number(total_valid)}\t{valid_pct:.2f}%\t")
        out.write(f"{format_number(mq_pass)}\t{mq_pass_pct:.2f}%\t")
        out.write(f"{format_number(vqslod_pass)}\t{vqslod_pass_pct:.2f}%\t")
        out.write(f"{format_number(both_pass)}\t{both_pass_pct:.2f}%\n")
    
    print(f"\n{'='*70}")
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] Summary for {args.chromosome}:")
    print(f"{'='*70}")
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}]   Total raw variants: {format_number(total_raw)}")
    if total_raw > 0:
        print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}]   Invalid variants: {format_number(total_invalid)} ({total_invalid/total_raw*100:.2f}%)")
        print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}]   Valid variants: {format_number(total_valid)} ({total_valid/total_raw*100:.2f}%)")
        print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}]   MQ > {args.mq_threshold}: {format_number(mq_pass)} ({mq_pass/total_raw*100:.2f}%)")
        print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}]   VQSLOD > {args.vqslod_threshold}: {format_number(vqslod_pass)} ({vqslod_pass/total_raw*100:.2f}%)")
        print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}]   Both pass: {format_number(both_pass)} ({both_pass/total_raw*100:.2f}%)")
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}]   Results saved to: {args.output_file}")
    print(f"{'='*70}\n")

if __name__ == "__main__":
    main()
