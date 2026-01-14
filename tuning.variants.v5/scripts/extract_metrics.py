#!/usr/bin/env python3
"""
Extract MQ and VQSLOD values from VCF files
Saves to TSV format for downstream analysis
"""

import sys
import subprocess
import argparse
import time
import os
from pathlib import Path

# Add script directory to path for imports
sys.path.insert(0, str(Path(__file__).parent))
from logger_utils import setup_logging

def extract_metrics_bcftools(vcf_file, bcftools_path, threads, output_file):
    """Extract MQ and VQSLOD values using bcftools and save to TSV"""
    
    # Use bcftools view with --threads to decompress in parallel, pipe to query
    # Output VCF format (-Ov) so query can read from stdin properly
    cmd = f"{bcftools_path} view --threads {threads} -Ov {vcf_file} | {bcftools_path} query -f '%INFO/MQ\\t%INFO/VQSLOD\\n'"
    
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] Extracting MQ and VQSLOD values...")
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] Using {threads} threads")
    start_time = time.time()
    
    try:
        # Stream output directly to file for better memory efficiency
        # Use larger buffer (64MB) and disable line buffering for maximum I/O speed
        with open(output_file, 'w', buffering=67108864) as out:  # 64MB buffer
            # Write header
            out.write("MQ\tVQSLOD\n")
            
            # Run bcftools with streaming output (no memory overhead)
            process = subprocess.Popen(
                cmd, 
                shell=True, 
                stdout=subprocess.PIPE, 
                stderr=subprocess.PIPE,
                text=True,
                bufsize=67108864  # 64MB buffer for pipe
            )
            
            # Ensure stdout is available
            if process.stdout is None:
                raise RuntimeError("Failed to capture stdout from bcftools process")
            
            line_count = 0
            batch_buffer = []
            batch_size = 100000  # Write every 100K lines for balanced performance
            last_log_time = time.time()
            last_write_time = time.time()
            log_interval = 10  # Log every 10 seconds
            write_interval = 30  # Force write every 30 seconds
            total_writes = 0
            
            print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] Processing variants (batch: {batch_size:,} lines or 30s, buffer: 64MB)...")
            
            # Use iter() with readline() for robust line iteration
            for line in iter(process.stdout.readline, ''):
                line = line.strip()
                if not line:
                    continue
                
                line_count += 1
                batch_buffer.append(line + '\n')
                
                current_time = time.time()
                
                # Write condition: either batch full OR 30 seconds elapsed
                should_write = (len(batch_buffer) >= batch_size or 
                               (batch_buffer and current_time - last_write_time >= write_interval))
                
                if should_write:
                    write_start = time.time()
                    out.writelines(batch_buffer)
                    out.flush()  # Force immediate write to disk
                    write_time = time.time() - write_start
                    total_writes += 1
                    batch_count = len(batch_buffer)
                    batch_buffer = []
                    last_write_time = current_time
                    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}]   ✓ Batch #{total_writes}: {batch_count:,} lines → Disk in {write_time:.2f}s (Total: {line_count:,})")
                
                # Time-based progress logging (every 10 seconds)
                if current_time - last_log_time >= log_interval:
                    elapsed = current_time - start_time
                    rate = line_count / elapsed
                    buffer_size = len(batch_buffer)
                    time_since_write = current_time - last_write_time
                    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}]   → {line_count:,} variants | {rate:.0f} var/s | Buffer: {buffer_size:,} | Last write: {time_since_write:.0f}s ago")
                    last_log_time = current_time
            
            # Write remaining buffered lines
            if batch_buffer:
                write_start = time.time()
                out.writelines(batch_buffer)
                out.flush()
                write_time = time.time() - write_start
                total_writes += 1
                print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}]   ✓ Final batch #{total_writes}: {len(batch_buffer):,} lines → Disk in {write_time:.2f}s")
            
            # Ensure all data is written with final flush
            out.flush()
            os.fsync(out.fileno())  # Force OS to write to physical disk
            print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}]   ✓ All data synchronized to disk ({total_writes} batches)")
            
            # Wait for process to complete and check return code
            process.wait()
            if process.returncode != 0:
                stderr_output = process.stderr.read() if process.stderr else "No error output"
                raise subprocess.CalledProcessError(process.returncode, cmd, stderr=stderr_output)
        
        elapsed = time.time() - start_time
        rate = line_count / elapsed if elapsed > 0 else 0
        print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] Extraction completed in {elapsed:.2f} seconds")
        print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] Total variants processed: {line_count:,}")
        print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] Average processing rate: {rate:.0f} variants/second")
        
        if line_count == 0:
            print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] Warning: No variants found in VCF file")
        
        print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] Output saved to: {output_file}")
        
        return line_count
    
    except subprocess.CalledProcessError as e:
        print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] Error running bcftools: {e}")
        print(f"stderr: {e.stderr}")
        sys.exit(1)
    except Exception as e:
        print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] Unexpected error: {e}")
        sys.exit(1)

def main():
    parser = argparse.ArgumentParser(
        description='Extract MQ and VQSLOD metrics from VCF files',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example:
  %(prog)s input.vcf.gz chr1 /path/to/bcftools 4 chr1.metrics.tsv
        """
    )
    
    parser.add_argument('vcf_file', 
                       help='Input VCF file (can be gzipped)')
    parser.add_argument('chromosome', 
                       help='Chromosome name for labeling (e.g., chr1, chr22, chrX)')
    parser.add_argument('bcftools_path', 
                       help='Path to bcftools executable')
    parser.add_argument('threads', type=int,
                       help='Number of threads for bcftools processing')
    parser.add_argument('output_file', 
                       help='Output TSV file with MQ and VQSLOD values')
    
    args = parser.parse_args()
    
    # Setup logging
    log_file = setup_logging(args.output_file)
    
    print(f"\n{'='*70}")
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] Starting metric extraction for {args.chromosome}")
    print(f"{'='*70}")
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] Input VCF: {args.vcf_file}")
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] Output TSV: {args.output_file}")
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] Log file: {log_file}")
    print(f"{'='*70}\n")
    
    # Extract metrics
    count = extract_metrics_bcftools(
        args.vcf_file, args.bcftools_path, args.threads, args.output_file
    )
    
    print(f"\n{'='*70}")
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] Metric extraction completed for {args.chromosome}")
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] Extracted {count:,} variant records")
    print(f"{'='*70}\n")

if __name__ == "__main__":
    main()
