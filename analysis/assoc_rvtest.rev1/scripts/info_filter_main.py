#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
VCF INFO Field Filtering Command Line Tool
==========================================

Command line interface for the info_filter_tools module, used to filter VCF files by INFO fields.

Typical Usage:
    python info_filter_main.py -i input.vcf.gz -k impact -v HIGH MODERATE -o filtered
    python info_filter_main.py -i input.vcf.gz -k impact -v HIGH -t 32 --check-chr-prefix

Author: ZHAO TIE
Date: 2025-10-14
"""

import sys
import os
import importlib
import argparse
import logging
import json
import time
from datetime import datetime

# Dynamically load custom module (maintain reload mode)
# Get current script directory
MODULE_DIR = os.path.dirname(os.path.abspath(__file__))
if MODULE_DIR not in sys.path:
    sys.path.insert(0, MODULE_DIR)

import info_filter_tools
importlib.reload(info_filter_tools)  # Maintain reload mode

from info_filter_tools import filter_vcf_by_info, get_vcf_stats


def setup_logging(verbose: bool = False) -> None:
    """Setup logging configuration"""
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )


def parse_arguments() -> argparse.Namespace:
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description='VCF INFO Field Filtering Tool',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s -i input.vcf.gz -k impact -v HIGH MODERATE -o filtered
  %(prog)s -i input.vcf.gz -k impact -v HIGH -t 32 --check-chr-prefix
  %(prog)s -i input.vcf.gz -k consequence -v "missense_variant" --no-keep-chr-prefix

Notes:
  - Input file must be bgzip compressed VCF (.vcf.gz)
  - values parameter supports multiple values separated by space
  - Automatically creates tabix index file
        """
    )
    
    # Required arguments
    parser.add_argument(
        '-i', '--input', '--anno-vcf',
        required=True,
        type=str,
        help='Input VCF file path (must be .vcf.gz)',
        metavar='FILE'
    )
    
    parser.add_argument(
        '-k', '--info-key',
        required=True,
        type=str,
        help='INFO field name to filter (e.g., impact, consequence)',
        metavar='KEY'
    )
    
    parser.add_argument(
        '-v', '--values',
        required=True,
        nargs='+',
        type=str,
        help='Target values for INFO field (multiple values supported, e.g., HIGH MODERATE)',
        metavar='VALUE'
    )
    
    # Optional arguments
    parser.add_argument(
        '-o', '--out-prefix',
        type=str,
        default='filtered',
        help='Output file prefix (default: filtered)',
        metavar='PREFIX'
    )
    
    parser.add_argument(
        '-t', '--threads',
        type=int,
        default=16,
        help='Number of parallel threads (default: 16)',
        metavar='N'
    )
    
    parser.add_argument(
        '--check-chr-prefix',
        action='store_true',
        help='Check chromosome prefix (default: False)'
    )
    
    parser.add_argument(
        '--keep-chr-prefix',
        action='store_true',
        help='Keep chromosome prefix (default: False)'
    )
    
    parser.add_argument(
        '--no-keep-chr-prefix',
        action='store_true',
        help='Do not keep chromosome prefix (opposite of --keep-chr-prefix)'
    )
    
    # Program options
    parser.add_argument(
        '--verbose',
        action='store_true',
        help='Verbose output (debug mode)'
    )
    
    parser.add_argument(
        '--dry-run',
        action='store_true',
        help='Show command to be executed only, do not run'
    )
    
    parser.add_argument(
        '--version',
        action='version',
        version='%(prog)s 1.1.0'
    )
    
    return parser.parse_args()


def validate_arguments(args: argparse.Namespace) -> None:
    """Validate command line arguments"""
    # Check input file
    if not os.path.exists(args.input):
        raise FileNotFoundError(f"Input file not found: {args.input}")
    
    if not args.input.endswith('.vcf.gz'):
        raise ValueError("Input file must be bgzip compressed VCF (.vcf.gz)")
    
    # Check threads
    if args.threads <= 0:
        raise ValueError(f"Threads must be greater than 0: {args.threads}")
    
    # Handle chromosome prefix conflict
    if args.keep_chr_prefix and args.no_keep_chr_prefix:
        raise ValueError("--keep-chr-prefix and --no-keep-chr-prefix cannot be used together")


def main() -> None:
    """Main function"""
    args = None
    start_time = time.time()
    
    try:
        # Parse arguments
        args = parse_arguments()
        
        # Setup logging
        setup_logging(args.verbose)
        
        # Validate arguments
        validate_arguments(args)
        
        # Determine keep_chr_prefix value
        if args.no_keep_chr_prefix:
            keep_chr_prefix = False
        else:
            keep_chr_prefix = args.keep_chr_prefix
        
        logging.info(f"Start processing VCF file: {args.input}")
        logging.info(f"INFO field: {args.info_key}")
        logging.info(f"Filter values: {args.values}")
        logging.info(f"Threads: {args.threads}")
        logging.info(f"Output prefix: {args.out_prefix}")
        
        if args.dry_run:
            logging.info("DRY RUN mode - Showing parameters only, no actual execution")
            filter_vcf_by_info(
                vcf_path=args.input,
                info_key=args.info_key,
                values=args.values,
                threads=args.threads,
                check_chr_prefix=args.check_chr_prefix,
                keep_chr_prefix=keep_chr_prefix,
                out_prefix=args.out_prefix,
                dry_run=True
            )
            return

        # 1. Collect Input Stats
        logging.info("Collecting input file statistics...")
        input_stats = get_vcf_stats(args.input)
        logging.info(f"Input stats: {input_stats}")
        
        # 2. Execute Filtering
        filtered_vcf = filter_vcf_by_info(
            vcf_path=args.input,
            info_key=args.info_key,
            values=args.values,
            threads=args.threads,
            check_chr_prefix=args.check_chr_prefix,
            keep_chr_prefix=keep_chr_prefix,
            out_prefix=args.out_prefix,
        )
        
        logging.info(f"Filtering complete! Output file: {filtered_vcf}")
        
        # 3. Collect Output Stats and Verify
        file_size_mb = 0
        output_stats = {}
        if os.path.exists(filtered_vcf):
            file_size_mb = os.path.getsize(filtered_vcf) / (1024**2)
            logging.info(f"Output file size: {file_size_mb:.2f} MB")
            
            logging.info("Collecting output file statistics...")
            output_stats = get_vcf_stats(filtered_vcf)
            logging.info(f"Output stats: {output_stats}")
        
        # 4. Generate Summary Log
        end_time = time.time()
        elapsed_time = end_time - start_time
        
        summary_data = {
            "date": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "input_file": os.path.abspath(args.input),
            "output_file": os.path.abspath(filtered_vcf),
            "filter_criteria": {
                "key": args.info_key,
                "values": args.values
            },
            "statistics": {
                "input": input_stats,
                "output": output_stats
            },
            "performance": {
                "elapsed_seconds": round(elapsed_time, 2),
                "output_size_mb": round(file_size_mb, 2)
            },
            "parameters": {
                "threads": args.threads,
                "check_chr_prefix": args.check_chr_prefix,
                "keep_chr_prefix": keep_chr_prefix
            }
        }
        
        # Construct summary filename
        # If output prefix is a path, place summary next to it
        if args.out_prefix.endswith('.vcf.gz'):
             # Should not happen typically based on usage, but handling just in case
             summary_path = args.out_prefix + ".summary.json"
        else:
             # If out_prefix is a directory or file prefix
             # We know filtered_vcf is the full path
             summary_path = filtered_vcf + ".summary.json"
             
        with open(summary_path, 'w', encoding='utf-8') as f:
            json.dump(summary_data, f, indent=4, ensure_ascii=False)
            
        logging.info(f"Summary log written to: {summary_path}")

    except KeyboardInterrupt:
        logging.error("Operation interrupted by user")
        sys.exit(130)
    except Exception as e:
        logging.error(f"Program execution error: {e}")
        if args and hasattr(args, 'verbose') and args.verbose:
            import traceback
            traceback.print_exc()
        sys.exit(1)


if __name__ == '__main__':
    main()
