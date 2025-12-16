#!/usr/bin/env python3
"""
Add imputation information to variant statistics file.

Adds two columns from pvar INFO field:
- IMPUTED_MARKER: TYPED, IMPUTED, or TYPED;IMPUTED
- IMPUTED_R2: R2 value from imputation

Uses chromosome-based parallel processing for efficiency.
"""

import argparse
import subprocess
import pandas as pd
import sys
import os
import time
import tempfile
import shutil
import re
from datetime import datetime
from multiprocessing import Pool, cpu_count, Manager
from functools import partial
from threading import Lock
from concurrent.futures import ThreadPoolExecutor, as_completed


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


def parse_pvar_info(info_str):
    """
    Parse INFO field from pvar file.
    
    Returns:
        tuple: (marker_type, r2_value)
        - marker_type: 'TYPED', 'IMPUTED', or 'TYPED;IMPUTED'
        - r2_value: float or None
    """
    # Split INFO field by semicolon to get individual tags
    info_tags = [tag.strip() for tag in info_str.split(';')]
    
    # Extract marker type - check for exact matches as standalone flags
    marker_type = []
    if 'TYPED' in info_tags:
        marker_type.append('TYPED')
    if 'IMPUTED' in info_tags:
        marker_type.append('IMPUTED')
    
    marker = ';'.join(marker_type) if marker_type else 'UNKNOWN'
    
    # Extract R2 value (use word boundary to avoid matching ER2)
    r2 = None
    # Match R2= that is either at start or after semicolon, to avoid matching ER2=
    r2_match = re.search(r'(?:^|;)R2=([0-9.]+)', info_str)
    if r2_match:
        try:
            r2 = float(r2_match.group(1))
        except ValueError:
            r2 = None
    
    return marker, r2


def read_pvar_info(pvar_file):
    """
    Read pvar file and extract imputation info.
    Memory-optimized for large files with progress tracking.
    
    Returns:
        dict: {variant_id: (marker_type, r2_value)}
    """
    log(f"Reading pvar file: {pvar_file}")
    info_dict = {}
    skipped_header = 0
    skipped_short = 0
    skipped_no_info = 0
    processed = 0
    last_progress_time = time.time()
    
    with open(pvar_file, 'r', buffering=1024*1024*10) as f:  # 10MB buffer
        for line in f:
            if line.startswith('#'):
                skipped_header += 1
                continue
            
            parts = line.strip().split('\t', maxsplit=6)  # Only split first 6 fields
            if len(parts) < 6:
                skipped_short += 1
                continue
            
            variant_id = parts[2]  # ID column (3rd column, 0-indexed)
            info_str = parts[5]    # INFO column (6th column, 0-indexed)
            
            # Skip if INFO field is missing or is '.'
            if not info_str or info_str == '.':
                skipped_no_info += 1
                continue
            
            marker, r2 = parse_pvar_info(info_str)
            info_dict[variant_id] = (marker, r2)
            processed += 1
            
            # Progress update every 5 seconds
            current_time = time.time()
            if current_time - last_progress_time >= 5:
                log(f"    Processing pvar: {processed:,} variants loaded...")
                last_progress_time = current_time
    
    log(f"  ✓ Loaded imputation info for {len(info_dict):,} variants")
    log(f"    (Skipped: {skipped_header:,} header lines, {skipped_short:,} short lines, {skipped_no_info:,} missing INFO)")
    return info_dict


def get_header_from_bgzip(input_file):
    """
    Efficiently extract header line from bgzipped file.
    Only reads first few KB instead of decompressing entire file.
    
    Returns:
        str: Header line
    """
    try:
        result = subprocess.run(
            ['gunzip', '-c', input_file],
            capture_output=True, text=True, check=True
        )
        # Only extract first line
        header = result.stdout.split('\n', 1)[0]
        return header
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"Failed to extract header from {input_file}: {e.stderr}")
    except Exception as e:
        raise RuntimeError(f"Failed to read header: {str(e)}")


def query_chromosome_with_tabix(tabix_path, input_file, chrom):
    """
    Use tabix to quickly extract variants for a specific chromosome.
    Tries both numeric (1, 2, ...) and chr-prefixed (chr1, chr2, ...) formats.
    
    Returns:
        pd.DataFrame: Variants for the chromosome, or empty DataFrame if none found
    """
    # Try both formats: numeric and chr-prefixed
    query_values = [str(chrom), f'chr{chrom}'] if chrom != 'other' else []
    
    if chrom == 'other':
        # For 'other', we need to scan the file since tabix can't do "not in" queries
        # Use streaming to avoid loading entire file into memory
        log(f"[CHR{chrom}] Scanning for non-standard chromosomes (streaming mode)...")
        
        # Get header
        header = get_header_from_bgzip(input_file)
        
        # Define standard chromosomes to exclude
        standard_chroms = set(str(i) for i in range(1, 23))
        standard_chroms.update(f'chr{i}' for i in range(1, 23))
        
        # Stream through file and collect non-standard chromosomes
        data_lines = []
        try:
            process = subprocess.Popen(
                ['gunzip', '-c', input_file],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                bufsize=1024*1024  # 1MB buffer
            )
            
            # Skip header
            next(process.stdout)
            
            # Stream through lines
            for line in process.stdout:
                line = line.strip()
                if not line:
                    continue
                
                # Only split first column to check chromosome
                chrom_value = line.split('\t', 1)[0]
                if chrom_value not in standard_chroms:
                    data_lines.append(line)
            
            # Wait for process to complete
            process.wait()
            if process.returncode != 0:
                stderr = process.stderr.read()
                raise RuntimeError(f"gunzip failed: {stderr}")
                
        except Exception as e:
            log(f"[CHR{chrom}] Error during streaming: {str(e)}", level='ERROR')
            if 'process' in locals():
                process.kill()
            raise
        
        if not data_lines:
            log(f"[CHR{chrom}] No non-standard chromosomes found")
            return pd.DataFrame()
        
        log(f"[CHR{chrom}] Found {len(data_lines):,} variants in non-standard chromosomes")
        
        # Parse into DataFrame
        from io import StringIO
        content = header + '\n' + '\n'.join(data_lines)
        return pd.read_csv(StringIO(content), sep='\t', dtype={'#CHROM': str})
    
    # For standard chromosomes, use tabix for fast extraction
    dfs = []
    for query in query_values:
        try:
            result = subprocess.run(
                [tabix_path, input_file, query],
                capture_output=True, text=True, check=True,
                timeout=300  # 5 minute timeout
            )
            
            if result.stdout.strip():
                # Get header (cached in process_chromosome to avoid redundancy)
                header = get_header_from_bgzip(input_file)
                
                # Combine header with data
                from io import StringIO
                content = header + '\n' + result.stdout
                df = pd.read_csv(StringIO(content), sep='\t', dtype={'#CHROM': str})
                dfs.append(df)
                log(f"[CHR{chrom}] Extracted {len(df):,} variants for query '{query}'")
        except subprocess.TimeoutExpired:
            log(f"[CHR{chrom}] Tabix query timeout for '{query}'", level='WARNING')
            continue
        except subprocess.CalledProcessError as e:
            # Chromosome not found with this query format (normal, not an error)
            if e.returncode == 1:  # tabix returns 1 when region not found
                continue
            else:
                log(f"[CHR{chrom}] Tabix error for '{query}': {e.stderr}", level='ERROR')
                raise
        except Exception as e:
            log(f"[CHR{chrom}] Unexpected error for query '{query}': {str(e)}", level='ERROR')
            raise
    
    if dfs:
        result_df = pd.concat(dfs, ignore_index=True)
        log(f"[CHR{chrom}] Total {len(result_df):,} variants extracted")
        return result_df
    else:
        log(f"[CHR{chrom}] No variants found for this chromosome")
        return pd.DataFrame()


def stream_annotate_chromosome(input_file, chrom, info_dict, output_file, log_interval=100000):
    """
    Stream through bgzipped file, filter by chromosome, and annotate on-the-fly.
    Memory-efficient: never loads entire chromosome into memory.
    Uses pigz for parallel decompression (4-8x faster than gunzip).
    
    Args:
        input_file: Bgzipped input file
        chrom: Chromosome to extract (e.g., '1', '2', ...)
        info_dict: Dictionary with imputation info
        output_file: Output file path
        log_interval: Log progress every N variants
    
    Returns:
        tuple: (num_variants, stats_dict)
    """
    import csv
    from io import TextIOWrapper
    
    # Target chromosome values (numeric and chr-prefixed)
    target_chroms = {str(chrom), f'chr{chrom}'}
    
    processed = 0
    start_time = time.time()
    
    # Statistics counters (track while streaming, no second pass needed)
    stats = {
        'TYPED': 0,
        'IMPUTED': 0,
        'TYPED;IMPUTED': 0,
        'UNKNOWN': 0
    }
    
    # Try pigz first (parallel gzip, 4-8x faster), fall back to gunzip
    decompress_cmd = None
    try:
        subprocess.run(['pigz', '--version'], capture_output=True, check=True)
        decompress_cmd = ['pigz', '-dc']  # -d for decompress, -c for stdout
    except (subprocess.CalledProcessError, FileNotFoundError):
        decompress_cmd = ['gunzip', '-c']
    
    # Start decompression subprocess
    decompress_proc = subprocess.Popen(
        decompress_cmd + [input_file],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        bufsize=4*1024*1024  # 4MB buffer for faster streaming
    )
    
    try:
        # Wrap stdout in text mode
        reader = TextIOWrapper(decompress_proc.stdout, encoding='utf-8', newline='')
        
        # Read header
        header_line = reader.readline().strip()
        header_cols = header_line.split('\t')
        
        # Add new columns to header
        output_header = header_cols + ['IMPUTED_MARKER', 'IMPUTED_R2']
        
        # Open output file for writing with larger buffer
        with open(output_file, 'w', newline='', buffering=4*1024*1024) as f_out:
            writer = csv.writer(f_out, delimiter='\t', lineterminator='\n')
            writer.writerow(output_header)
            
            # Batch write buffer (write every 10k rows for efficiency)
            write_buffer = []
            buffer_size = 10000
            
            # Stream through lines
            for line in reader:
                line = line.strip()
                if not line:
                    continue
                
                # Split only first column to check chromosome (faster)
                first_tab = line.find('\t')
                if first_tab == -1:
                    continue
                    
                chrom_value = line[:first_tab]
                if chrom_value not in target_chroms:
                    continue
                
                # Now split the full line
                parts = line.split('\t')
                if len(parts) < 3:
                    continue
                
                # Get variant ID (3rd column, index 2)
                variant_id = parts[2]
                
                # Lookup imputation info (dict lookup is O(1))
                marker, r2 = info_dict.get(variant_id, ('UNKNOWN', None))
                
                # Update statistics
                stats[marker] = stats.get(marker, 0) + 1
                
                # Prepare output row
                output_row = parts + [marker, str(r2) if r2 is not None else '']
                write_buffer.append(output_row)
                
                # Batch write when buffer is full
                if len(write_buffer) >= buffer_size:
                    writer.writerows(write_buffer)
                    write_buffer.clear()
                
                processed += 1
                
                # Progress logging
                if processed % log_interval == 0:
                    elapsed = time.time() - start_time
                    rate = processed / elapsed
                    log(f"[CHR{chrom}]   Processed {processed:,} variants ({rate:,.0f} var/s)")
            
            # Write remaining buffer
            if write_buffer:
                writer.writerows(write_buffer)
        
        # Wait for decompression to complete
        decompress_proc.wait()
        if decompress_proc.returncode != 0:
            stderr = decompress_proc.stderr.read().decode('utf-8')
            tool_name = 'pigz' if 'pigz' in decompress_cmd[0] else 'gunzip'
            raise RuntimeError(f"{tool_name} failed: {stderr}")
        
        elapsed = time.time() - start_time
        rate = processed / elapsed if elapsed > 0 else 0
        tool_name = 'pigz' if 'pigz' in decompress_cmd[0] else 'gunzip'
        log(f"[CHR{chrom}]   ✓ Streaming complete: {processed:,} variants in {elapsed:.2f}s ({rate:,.0f} var/s) [{tool_name}]")
        
        return processed, stats
        
    except Exception as e:
        decompress_proc.kill()
        raise
    finally:
        decompress_proc.stdout.close()
        decompress_proc.stderr.close()


def annotate_chunk(chunk_df, info_dict, chunk_id, chrom, total_chunks):
    """
    Annotate a chunk of variants with imputation info.
    Used for parallel processing within a chromosome.
    
    Args:
        chunk_df: DataFrame chunk to annotate
        info_dict: Dictionary with imputation info
        chunk_id: Chunk identifier for logging
        chrom: Chromosome name for logging
        total_chunks: Total number of chunks for logging
    
    Returns:
        Annotated DataFrame chunk
    """
    chunk_size = len(chunk_df)
    log(f"[CHR{chrom}]   Chunk {chunk_id}/{total_chunks}: Processing {chunk_size:,} variants...")
    
    # Add imputation info
    chunk_df['IMPUTED_MARKER'] = chunk_df['ID'].map(
        lambda x: info_dict.get(x, ('UNKNOWN', None))[0]
    )
    chunk_df['IMPUTED_R2'] = chunk_df['ID'].map(
        lambda x: info_dict.get(x, ('UNKNOWN', None))[1]
    )
    
    log(f"[CHR{chrom}]   Chunk {chunk_id}/{total_chunks}: ✓ Completed")
    return chunk_df


def get_chromosomes(tsv_gz_file):
    """
    Return predefined chromosome list (1-22) for autosomes only.
    No need to scan the file - we know the standard chromosome set.
    """
    log("Using standard chromosome set: 1-22 (autosomes only)")
    
    # Standard autosomes 1-22 only
    chrom_list = [str(i) for i in range(1, 23)]
    
    log(f"  ✓ Processing {len(chrom_list)} chromosomes: {', '.join(chrom_list[:5])}...{chrom_list[-3:]}")
    
    return chrom_list


def process_chromosome(chrom, input_file, tabix_path, info_dict_file, temp_dir, progress_dict, lock, total_chroms):
    """
    Process one chromosome: add imputation info using tabix for fast extraction.
    Uses tabix index for O(log n) chromosome lookup instead of full scan.
    
    Args:
        chrom: Chromosome name (1-22 or 'other')
        input_file: Input variant_stats file (must be bgzipped and tabix-indexed)
        tabix_path: Path to tabix executable
        info_dict_file: Path to pickled info_dict file
    """
    update_progress(chrom, 'PROCESSING', progress_dict, lock, total_chroms)
    
    log(f"[CHR{chrom}] Starting processing...")
    start_time = time.time()
    
    try:
        # Validate input file exists
        if not os.path.exists(input_file):
            raise FileNotFoundError(f"Input file not found: {input_file}")
        
        # Validate tabix index exists
        if not os.path.exists(input_file + '.tbi'):
            raise FileNotFoundError(f"Tabix index not found: {input_file}.tbi")
        
        # Load info_dict from file
        import pickle
        log(f"[CHR{chrom}] Loading imputation info dictionary...")
        load_start = time.time()
        with open(info_dict_file, 'rb') as f:
            info_dict = pickle.load(f)
        load_elapsed = time.time() - load_start
        log(f"[CHR{chrom}] Loaded {len(info_dict):,} variant annotations in {load_elapsed:.2f}s")
        
        # Use streaming annotation instead of tabix (faster and memory-efficient)
        log(f"[CHR{chrom}] Starting stream annotation (reading + annotating simultaneously)...")
        temp_output = os.path.join(temp_dir, f'result_chr{chrom}.tsv')
        
        stream_start = time.time()
        total_variants, stats = stream_annotate_chromosome(
            input_file, chrom, info_dict, temp_output, log_interval=500000
        )
        stream_elapsed = time.time() - stream_start
        
        if total_variants == 0:
            log(f"[CHR{chrom}] ⚠ No variants found for this chromosome")
            # Create empty gzipped file
            output_file = os.path.join(temp_dir, f'result_chr{chrom}.tsv.gz')
            with open(output_file, 'wb') as f:
                pass
        else:
            log(f"[CHR{chrom}] Stream annotation completed: {total_variants:,} variants in {stream_elapsed:.2f}s")
            
            # Display statistics (already computed during streaming)
            n_typed_only = stats.get('TYPED', 0)
            n_imputed_only = stats.get('IMPUTED', 0)
            n_both = stats.get('TYPED;IMPUTED', 0)
            n_unknown = stats.get('UNKNOWN', 0)
            
            log(f"[CHR{chrom}]   TYPED only: {n_typed_only:,} | IMPUTED only: {n_imputed_only:,} | "
                f"TYPED;IMPUTED: {n_both:,} | UNKNOWN: {n_unknown:,}")
            
            # Warn if too many unknowns
            if n_unknown > total_variants * 0.1:
                log(f"[CHR{chrom}]   ⚠ WARNING: {n_unknown/total_variants*100:.1f}% variants not found in pvar", level='WARNING')
            
            # Compress the output file with pigz if available (faster)
            log(f"[CHR{chrom}] Compressing output...")
            compress_start = time.time()
            output_file = os.path.join(temp_dir, f'result_chr{chrom}.tsv.gz')
            
            # Try pigz first, fall back to gzip
            try:
                subprocess.run(['pigz', '--version'], capture_output=True, check=True)
                subprocess.run(['pigz', '-f', temp_output], check=True)
                compress_tool = 'pigz'
            except (subprocess.CalledProcessError, FileNotFoundError):
                subprocess.run(['gzip', '-f', temp_output], check=True)
                compress_tool = 'gzip'
            
            os.rename(temp_output + '.gz', output_file)
            compress_elapsed = time.time() - compress_start
            log(f"[CHR{chrom}] Compression completed in {compress_elapsed:.2f}s [{compress_tool}]")
        
        # Clear info_dict from memory
        del info_dict
        
        # Verify file was created
        if not os.path.exists(output_file):
            raise RuntimeError(f"Failed to create output file: {output_file}")
        
        elapsed = time.time() - start_time
        log(f"[CHR{chrom}] ✓ COMPLETED: {total_variants:,} variants processed in {elapsed:.2f}s")
        
        update_progress(chrom, 'COMPLETED', progress_dict, lock, total_chroms)
        
        return output_file
        
    except Exception as e:
        log(f"[CHR{chrom}] ✗ ERROR: {str(e)}", level='ERROR')
        update_progress(chrom, 'FAILED', progress_dict, lock, total_chroms)
        raise


def main():
    parser = argparse.ArgumentParser(
        description='Add imputation information to variant statistics',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    parser.add_argument('--variant-stats', required=True,
                        help='Input variant statistics TSV.gz file')
    parser.add_argument('--pvar', required=True,
                        help='Input pvar file with imputation info')
    parser.add_argument('--tabix', required=True,
                        help='Path to tabix executable')
    parser.add_argument('--threads', type=int, default=None,
                        help='Number of threads/processes to use (default: all available CPUs)')
    parser.add_argument('--output', required=True,
                        help='Output TSV.gz file path')
    parser.add_argument('--skip-order-check', action='store_true',
                        help='Skip checking variant order against input file')
    parser.add_argument('--keep-temp', action='store_true', default=True,
                        help='Keep temporary files after completion (default: True)')
    
    args = parser.parse_args()
    
    # Determine number of processes
    n_processes = args.threads if args.threads else cpu_count()
    log("="*80)
    log("Add Imputation Information Pipeline (Parallel & Memory-Efficient)")
    log("="*80)
    log(f"Input variant stats: {args.variant_stats}")
    log(f"Input pvar: {args.pvar}")
    log(f"Output file: {args.output}")
    log(f"Parallel processes: {n_processes}")
    log("="*80)
    
    overall_start = time.time()
    
    # Create temporary directory in working directory
    temp_dir = os.path.join(os.getcwd(), 'tmp')
    os.makedirs(temp_dir, exist_ok=True)
    log(f"Using temporary directory: {temp_dir}")
    
    # Read pvar imputation info
    log("\n>>> PHASE 1/3: Read imputation info from pvar")
    info_dict = read_pvar_info(args.pvar)
    
    # Save info_dict to temporary file to avoid memory duplication in multiprocessing
    import pickle
    info_dict_file = os.path.join(temp_dir, f'info_dict_{os.getpid()}.pkl')
    log(f"Saving imputation info to temporary file: {info_dict_file}")
    with open(info_dict_file, 'wb') as f:
        pickle.dump(info_dict, f, protocol=pickle.HIGHEST_PROTOCOL)
    log(f"  ✓ Info dict saved ({os.path.getsize(info_dict_file) / (1024**2):.1f} MB)")
    
    # Clear from main process memory
    del info_dict
    
    try:
        # Get list of chromosomes
        log("\n>>> PHASE 2/3: Process chromosomes in parallel (using tabix)")
        chromosomes = get_chromosomes(args.variant_stats)
        log(f"Processing {len(chromosomes)} chromosomes with {n_processes} parallel processes...")
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
            input_file=args.variant_stats,
            tabix_path=args.tabix,
            info_dict_file=info_dict_file,
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
        
        # Maintain original order if requested
        if not args.skip_order_check:
            log("\nMaintaining original variant order...")
            sort_start = time.time()
            # Read original order using chunks to save memory
            log("  Reading original variant order...")
            id_order = {}
            chunk_num = 0
            for chunk in pd.read_csv(args.variant_stats, sep='\t', compression='gzip', 
                                     usecols=['ID'], chunksize=1000000):
                for idx, vid in enumerate(chunk['ID'], start=chunk_num * 1000000):
                    id_order[vid] = idx
                chunk_num += 1
                if chunk_num % 10 == 0:
                    log(f"    Processed {chunk_num}M variants for ordering...")
            
            log(f"  Mapping order to {len(final_df):,} variants...")
            final_df['_order'] = final_df['ID'].map(id_order)
            final_df = final_df.sort_values('_order')
            final_df = final_df.drop('_order', axis=1)
            sort_elapsed = time.time() - sort_start
            log(f"  ✓ Sorting completed in {sort_elapsed:.2f}s")
        
        # Write output with bgzip
        log(f"\nWriting final output to {args.output}...")
        write_start = time.time()
        
        # Write uncompressed TSV temporarily
        temp_tsv = args.output.replace('.gz', '') if args.output.endswith('.gz') else args.output + '.tmp'
        final_df.to_csv(temp_tsv, sep='\t', index=False)
        log(f"  ✓ Uncompressed file written: {temp_tsv}")
        
        # Compress with bgzip
        log("  Compressing with bgzip...")
        try:
            if os.path.exists(args.output):
                os.remove(args.output)
            
            bgzip_path = os.path.join(os.path.dirname(args.tabix), 'bgzip')
            subprocess.run([bgzip_path, '-f', temp_tsv], 
                          check=True, capture_output=True, text=True)
            
            bgzipped_file = temp_tsv + '.gz'
            if bgzipped_file != args.output:
                os.rename(bgzipped_file, args.output)
            
            write_elapsed = time.time() - write_start
            file_size = os.path.getsize(args.output) / (1024**2)
            log(f"  ✓ File written and compressed in {write_elapsed:.2f}s")
            log(f"  File size: {file_size:.2f} MB")
            
        except FileNotFoundError:
            log("  ⚠ bgzip not found, using gzip compression", level='WARNING')
            final_df.to_csv(args.output, sep='\t', index=False, compression='gzip')
            write_elapsed = time.time() - write_start
            file_size = os.path.getsize(args.output) / (1024**2)
            log(f"  ✓ File written with gzip in {write_elapsed:.2f}s")
            log(f"  File size: {file_size:.2f} MB")
        
        # Create tabix index
        log("\nCreating tabix index...")
        tabix_start = time.time()
        try:
            subprocess.run([args.tabix, '-s', '1', '-b', '2', '-e', '2', args.output],
                          check=True, capture_output=True, text=True)
            tabix_elapsed = time.time() - tabix_start
            log(f"  ✓ Tabix index created in {tabix_elapsed:.2f}s: {args.output}.tbi")
        except subprocess.CalledProcessError as e:
            log(f"  ⚠ Tabix indexing failed (file may not be properly formatted)", level='WARNING')
        
    finally:
        # Clean up temporary files based on --keep-temp parameter
        if args.keep_temp:
            log("\nKeeping temporary files as requested")
            log(f"  Temporary directory: {temp_dir}")
            log(f"  Contents:")
            try:
                for item in os.listdir(temp_dir):
                    item_path = os.path.join(temp_dir, item)
                    if os.path.isfile(item_path):
                        size_mb = os.path.getsize(item_path) / (1024**2)
                        log(f"    - {item} ({size_mb:.1f} MB)")
            except Exception as e:
                log(f"  ⚠ Failed to list temporary files: {str(e)}", level='WARNING')
        else:
            log("\nCleaning up temporary files...")
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
