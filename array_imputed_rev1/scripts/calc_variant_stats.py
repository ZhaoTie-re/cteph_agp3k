#!/usr/bin/env python3
"""
Calculate variant-level statistics for all samples, cases, and controls.

Statistics calculated:
- VMISS: Variant missingness rate (fraction of samples missing this variant)
- AAF: Alternative allele frequency
- MAF: Minor allele frequency
- HWE: Hardy-Weinberg equilibrium p-value

Output format: TSV file with columns:
#CHROM POS ID REF ALT VMISS_ALL VMISS_CASE VMISS_CTRL 
AAF_ALL AAF_CASE AAF_CTRL MAF_ALL MAF_CASE MAF_CTRL HWE_ALL HWE_CASE HWE_CTRL

Uses chunking and multiprocessing for memory efficiency with large datasets.
"""

import argparse
import subprocess
import pandas as pd
import gzip
import sys
import os
import time
import tempfile
import shutil
from datetime import datetime
from multiprocessing import Pool, cpu_count, Manager
from functools import partial
from threading import Lock


def log(message, level='INFO'):
    """Print log message with timestamp."""
    timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    print(f"[{timestamp}] [{level}] {message}", file=sys.stderr, flush=True)


def update_progress(chrom, status, progress_dict, lock, total_chroms):
    """Thread-safe progress update."""
    with lock:
        progress_dict[chrom] = status
        completed = sum(1 for s in progress_dict.values() if s == 'COMPLETED')
        in_progress = sum(1 for s in progress_dict.values() if s == 'PROCESSING')
        pending = total_chroms - completed - in_progress
        
        log(f"Progress: [{completed}/{total_chroms}] completed | "
            f"{in_progress} in progress | {pending} pending | "
            f"Current: chr{chrom} -> {status}")


def get_chromosomes(pvar_file):
    """Extract unique chromosomes from pvar file."""
    log("Extracting chromosomes from pvar file...")
    chromosomes = set()
    with open(pvar_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            chrom = line.split('\t')[0]
            chromosomes.add(chrom)
    
    chrom_list = sorted(chromosomes, key=lambda x: (x.replace('chr', '').zfill(2) if x.startswith('chr') else x))
    log(f"  ✓ Found {len(chrom_list)} chromosomes: {', '.join(chrom_list)}")
    return chrom_list


def process_chromosome(chrom, pfile, plink2_path, pheno_col, temp_dir, progress_dict, lock, total_chroms):
    """Process statistics for a single chromosome."""
    update_progress(chrom, 'PROCESSING', progress_dict, lock, total_chroms)
    
    log(f"[CHR{chrom}] Starting processing...")
    start_time = time.time()
    
    # Create temp prefix for this chromosome
    temp_prefix = os.path.join(temp_dir, f'chr{chrom}')
    
    try:
        # Extract this chromosome
        extract_cmd = [
            plink2_path,
            '--pfile', pfile,
            '--chr', str(chrom),
            '--make-pgen',
            '--out', temp_prefix
        ]
        log(f"[CHR{chrom}] Step 1/4: Extracting chromosome data...", level='DEBUG')
        subprocess.run(extract_cmd, check=True, capture_output=True, text=True)
        
        # Calculate stats for ALL, CASE, CTRL
        results = {}
        group_step = 2
        for group, filter_args in [
            ('ALL', []),
            ('CASE', ['--keep-if', f'{pheno_col} == 2']),
            ('CTRL', ['--keep-if', f'{pheno_col} == 1'])
        ]:
            log(f"[CHR{chrom}] Step {group_step}/4: Calculating statistics for {group}...", level='DEBUG')
            group_prefix = f'{temp_prefix}_{group}'
            
            base_cmd = [plink2_path, '--pfile', temp_prefix] + filter_args
            
            # Frequency
            subprocess.run(base_cmd + ['--freq', '--out', group_prefix], 
                          check=True, capture_output=True, text=True)
            
            # Missingness
            subprocess.run(base_cmd + ['--missing', 'variant-only', '--out', group_prefix],
                          check=True, capture_output=True, text=True)
            
            # Hardy-Weinberg
            subprocess.run(base_cmd + ['--hardy', '--out', group_prefix],
                          check=True, capture_output=True, text=True)
            
            # Read results with chunking to save memory
            freq_df = pd.read_csv(f'{group_prefix}.afreq', sep='\s+', 
                                 usecols=['#CHROM', 'ID', 'REF', 'ALT', 'ALT_FREQS'])
            miss_df = pd.read_csv(f'{group_prefix}.vmiss', sep='\s+',
                                 usecols=['ID', 'F_MISS'])
            hardy_df = pd.read_csv(f'{group_prefix}.hardy', sep='\s+',
                                  usecols=['ID', 'P'])
            
            # Merge
            merged = freq_df.merge(miss_df, on='ID', how='left')
            merged = merged.merge(hardy_df, on='ID', how='left')
            results[group] = merged
            
            group_step += 1
        
        # Combine all groups
        log(f"[CHR{chrom}] Step 5/5: Combining and formatting results...", level='DEBUG')
        final = results['ALL'][['#CHROM', 'ID', 'REF', 'ALT']].copy()
        
        # Extract POS from ID (format: CHROM:POS:REF:ALT)
        final['POS'] = final['ID'].str.split(':', expand=True)[1]
        final['POS'] = pd.to_numeric(final['POS'], errors='coerce')
        
        # Add statistics
        for stat in ['F_MISS', 'ALT_FREQS', 'P']:
            for group in ['ALL', 'CASE', 'CTRL']:
                col_name_map = {
                    'F_MISS': f'VMISS_{group}',
                    'ALT_FREQS': f'AAF_{group}',
                    'P': f'HWE_{group}'
                }
                final[col_name_map[stat]] = results[group][stat]
        
        # Calculate MAF
        for group in ['ALL', 'CASE', 'CTRL']:
            aaf_col = f'AAF_{group}'
            maf_col = f'MAF_{group}'
            final[maf_col] = final[aaf_col].apply(lambda x: min(x, 1-x) if pd.notna(x) else None)
        
        # Reorder columns
        final = final[['#CHROM', 'POS', 'ID', 'REF', 'ALT',
                      'VMISS_ALL', 'VMISS_CASE', 'VMISS_CTRL',
                      'AAF_ALL', 'AAF_CASE', 'AAF_CTRL',
                      'MAF_ALL', 'MAF_CASE', 'MAF_CTRL',
                      'HWE_ALL', 'HWE_CASE', 'HWE_CTRL']]
        
        # Save to temporary file
        output_file = os.path.join(temp_dir, f'result_chr{chrom}.tsv.gz')
        final.to_csv(output_file, sep='\t', index=False, compression='gzip')
        
        elapsed = time.time() - start_time
        log(f"[CHR{chrom}] ✓ COMPLETED: {len(final):,} variants processed in {elapsed:.2f}s ({len(final)/elapsed:.0f} variants/sec)")
        
        update_progress(chrom, 'COMPLETED', progress_dict, lock, total_chroms)
        
        return output_file
        
    except Exception as e:
        log(f"[CHR{chrom}] ✗ ERROR: {str(e)}", level='ERROR')
        update_progress(chrom, 'FAILED', progress_dict, lock, total_chroms)
        raise
    finally:
        # Clean up temporary files for this chromosome
        cleaned = 0
        for f in os.listdir(temp_dir):
            if f.startswith(f'chr{chrom}') and not f.startswith(f'result_chr{chrom}'):
                try:
                    os.remove(os.path.join(temp_dir, f))
                    cleaned += 1
                except:
                    pass
        if cleaned > 0:
            log(f"[CHR{chrom}] Cleaned up {cleaned} temporary files", level='DEBUG')


def main():
    parser = argparse.ArgumentParser(
        description='Calculate variant statistics for all samples, cases, and controls (optimized for large datasets)',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    parser.add_argument('--pfile', required=True,
                        help='Input PLINK2 pfile prefix')
    parser.add_argument('--plink2', required=True,
                        help='Path to PLINK2 executable')
    parser.add_argument('--tabix', required=True,
                        help='Path to tabix executable')
    parser.add_argument('--pheno-col', default='PHENO1',
                        help='Phenotype column name in .psam file')
    parser.add_argument('--threads', type=int, default=None,
                        help='Number of threads/processes to use (default: all available CPUs)')
    parser.add_argument('--output', required=True,
                        help='Output TSV.gz file path')
    parser.add_argument('--skip-order-check', action='store_true',
                        help='Skip checking variant order against pvar file (faster but may not match pvar order)')
    
    args = parser.parse_args()
    
    # Determine number of processes
    n_processes = args.threads if args.threads else cpu_count()
    log("="*80)
    log("Variant Statistics Calculation Pipeline (Parallel & Memory-Efficient)")
    log("="*80)
    log(f"Input pfile: {args.pfile}")
    log(f"Output file: {args.output}")
    log(f"Phenotype column: {args.pheno_col}")
    log(f"Parallel processes: {n_processes}")
    log("="*80)
    
    overall_start = time.time()
    
    # Create temporary directory
    temp_dir = tempfile.mkdtemp(prefix='variant_stats_')
    log(f"Created temporary directory: {temp_dir}")
    
    try:
        # Get list of chromosomes
        log("\n>>> PHASE 1/3: Identify chromosomes")
        chromosomes = get_chromosomes(f'{args.pfile}.pvar')
        
        # Process chromosomes in parallel
        log(f"\n>>> PHASE 2/3: Process {len(chromosomes)} chromosomes in parallel")
        log(f"Using {n_processes} parallel processes...")
        log("="*80)
        
        # Create shared progress tracking
        manager = Manager()
        progress_dict = manager.dict()
        lock = manager.Lock()
        
        # Initialize progress for all chromosomes
        for chrom in chromosomes:
            progress_dict[chrom] = 'PENDING'
        
        process_func = partial(
            process_chromosome,
            pfile=args.pfile,
            plink2_path=args.plink2,
            pheno_col=args.pheno_col,
            temp_dir=temp_dir,
            progress_dict=progress_dict,
            lock=lock,
            total_chroms=len(chromosomes)
        )
        
        phase2_start = time.time()
        with Pool(processes=n_processes) as pool:
            result_files = pool.map(process_func, chromosomes)
        
        phase2_elapsed = time.time() - phase2_start
        log("="*80)
        log(f"  ✓ All {len(chromosomes)} chromosomes processed successfully")
        log(f"  Total processing time: {phase2_elapsed:.2f}s ({phase2_elapsed/60:.2f} minutes)")
        log(f"  Average time per chromosome: {phase2_elapsed/len(chromosomes):.2f}s")
        log(f"  Speedup factor: {phase2_elapsed/(phase2_elapsed/len(chromosomes)*len(chromosomes)/n_processes):.1f}x")
        
        # Merge results
        log("\n>>> PHASE 3/3: Merge chromosome results and create final output")
        log("Concatenating chromosome results...")
        
        # Read and concatenate all chromosome results
        chunk_dfs = []
        total_variants = 0
        for i, result_file in enumerate(result_files, 1):
            if os.path.exists(result_file):
                df = pd.read_csv(result_file, sep='\t', compression='gzip')
                chunk_dfs.append(df)
                total_variants += len(df)
                chrom_name = os.path.basename(result_file).replace('result_', '').replace('.tsv.gz', '')
                log(f"  [{i}/{len(result_files)}] Loaded {chrom_name}: {len(df):,} variants (cumulative: {total_variants:,})")
        
        if not chunk_dfs:
            raise ValueError("No results to merge!")
        
        # Concatenate all chunks
        log(f"\nConcatenating {len(chunk_dfs)} chromosome datasets...")
        concat_start = time.time()
        final_df = pd.concat(chunk_dfs, ignore_index=True)
        concat_elapsed = time.time() - concat_start
        log(f"  ✓ Concatenation completed in {concat_elapsed:.2f}s")
        log(f"  Total variants: {len(final_df):,}")
        
        # Sort by original pvar order (if not skipped)
        if not args.skip_order_check:
            log("\nReading original pvar file for correct variant order...")
            pvar_ids = []
            with open(f'{args.pfile}.pvar', 'r') as f:
                for line in f:
                    if line.startswith('#'):
                        continue
                    parts = line.strip().split('\t')
                    if len(parts) >= 3:
                        pvar_ids.append(parts[2])  # ID is the 3rd column
            log(f"  ✓ Read {len(pvar_ids):,} variant IDs from pvar file")
            
            log("\nSorting by original pvar order...")
            sort_start = time.time()
            # Create a mapping of ID to original order
            id_order = {vid: i for i, vid in enumerate(pvar_ids)}
            final_df['_order'] = final_df['ID'].map(id_order)
            final_df = final_df.sort_values('_order')
            final_df = final_df.drop('_order', axis=1)
            sort_elapsed = time.time() - sort_start
            log(f"  ✓ Sorting completed in {sort_elapsed:.2f}s")
        else:
            log("\nSkipping pvar order check (sorting by chromosome and position)...")
            sort_start = time.time()
            final_df = final_df.sort_values(['#CHROM', 'POS'])
            sort_elapsed = time.time() - sort_start
            log(f"  ✓ Sorting completed in {sort_elapsed:.2f}s")
        
        # Write final output (uncompressed first, then bgzip for tabix compatibility)
        log(f"\nWriting final output to {args.output}...")
        write_start = time.time()
        
        # Write uncompressed TSV temporarily
        temp_tsv = args.output.replace('.gz', '') if args.output.endswith('.gz') else args.output + '.tmp'
        final_df.to_csv(temp_tsv, sep='\t', index=False)
        log(f"  ✓ Uncompressed file written: {temp_tsv}")
        
        # Compress with bgzip for tabix compatibility
        log("  Compressing with bgzip...")
        try:
            # Remove existing .gz file if present
            if os.path.exists(args.output):
                os.remove(args.output)
            
            # Use bgzip (part of htslib, same directory as tabix)
            bgzip_path = os.path.join(os.path.dirname(args.tabix), 'bgzip')
            subprocess.run([bgzip_path, '-f', temp_tsv], 
                          check=True, capture_output=True, text=True)
            
            # bgzip creates .gz file automatically
            bgzipped_file = temp_tsv + '.gz'
            if bgzipped_file != args.output:
                os.rename(bgzipped_file, args.output)
            
            write_elapsed = time.time() - write_start
            file_size = os.path.getsize(args.output) / (1024**2)
            log(f"  ✓ File written and compressed in {write_elapsed:.2f}s")
            log(f"  File size: {file_size:.2f} MB")
            
        except FileNotFoundError:
            log("  ⚠ bgzip not found, using gzip compression (tabix may not work)", level='WARNING')
            final_df.to_csv(args.output, sep='\t', index=False, compression='gzip')
            write_elapsed = time.time() - write_start
            file_size = os.path.getsize(args.output) / (1024**2)
            log(f"  ✓ File written with gzip in {write_elapsed:.2f}s")
            log(f"  File size: {file_size:.2f} MB")
        
        # Create tabix index
        log("\nCreating tabix index...")
        tabix_start = time.time()
        try:
            # tabix parameters: -s (sequence/chrom col), -b (begin col), -e (end col)
            # Our format: #CHROM(1) POS(2) ID REF ALT ...
            # Use POS for both begin and end since variants are single positions
            subprocess.run([args.tabix, '-s', '1', '-b', '2', '-e', '2', args.output],
                          check=True, capture_output=True, text=True)
            tabix_elapsed = time.time() - tabix_start
            log(f"  ✓ Tabix index created in {tabix_elapsed:.2f}s: {args.output}.tbi")
            log(f"  Query example: tabix {args.output} chr1:12345-12345")
        except subprocess.CalledProcessError as e:
            log(f"  ✗ Failed to create tabix index", level='ERROR')
            log(f"  Error output: {e.stderr}", level='ERROR')
            log(f"  This may happen if file is not bgzip-compressed or not properly sorted", level='WARNING')
            # Don't raise - tabix is optional
            log(f"  ⚠ Continuing without tabix index", level='WARNING')
        
    finally:
        # Clean up temporary directory
        log("Cleaning up temporary files...")
        try:
            shutil.rmtree(temp_dir)
            log(f"  ✓ Removed temporary directory: {temp_dir}")
        except Exception as e:
            log(f"  ⚠ Failed to remove temporary directory: {str(e)}", level='WARNING')
    
    overall_elapsed = time.time() - overall_start
    log("="*80)
    log(f"Pipeline completed successfully!")
    log(f"Total time: {overall_elapsed:.2f} seconds ({overall_elapsed/60:.2f} minutes)")
    log(f"Average time per chromosome: {overall_elapsed/len(chromosomes):.2f} seconds")
    log("="*80)


if __name__ == '__main__':
    main()
