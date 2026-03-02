#!/usr/bin/env python3
"""
snpEff VCF Annotation Tool - Command Line Interface

This tool performs snpEff functional annotation on VCF files, supporting parallel processing.
Key Features:
1. Automatically adds 'chr' prefix to VCF if needed.
2. Annotates VCF with functional information using pre-computed snpEff TSV files.
3. Supports both parallel and serial processing modes.
4. Provides detailed logging and progress tracking.
5. Generates comprehensive statistics on variants and their classification by chromosome.

Author: ZHAO TIE
Date: Jan 19, 2026
Version: 1.1.0
"""

import sys
import argparse
import logging
import subprocess
from pathlib import Path
from datetime import datetime
from collections import defaultdict, Counter
import shutil

# Dynamic module path configuration to ensure robustness
current_dir = Path(__file__).resolve().parent
if str(current_dir) not in sys.path:
    sys.path.insert(0, str(current_dir))

# Import custom modules
try:
    from snpeff_anno_tools import add_chr_prefix_to_vcf, annotate_vcf_with_snpeff_tsv
except ImportError as e:
    print(f"Error importing snpeff_anno_tools: {e}")
    sys.path.append("/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/analysis/assoc_rvtest.rev1/scripts")
    try:
        from snpeff_anno_tools import add_chr_prefix_to_vcf, annotate_vcf_with_snpeff_tsv
    except ImportError:
         print("Critical Error: Could not load necessary modules.")
         sys.exit(1)


def setup_argument_parser():
    """
    Set up command line argument parser
    
    Returns:
        argparse.ArgumentParser: Configured argument parser
    """
    parser = argparse.ArgumentParser(
        description="snpEff VCF Annotation Tool - Add functional annotations to VCF files",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Usage Examples:
  # Basic usage (using default parameters)
  python snpeff_anno_main.py
  
  # Custom input file
  python snpeff_anno_main.py --vcf-path /path/to/your/file.vcf.gz
  
  # Custom annotation directory
  python snpeff_anno_main.py --snpeff-dir /path/to/annotations
  
  # Adjust performance parameters
  python snpeff_anno_main.py --max-workers 8 --threads 4
  
  # Use serial mode
  python snpeff_anno_main.py --no-parallel
  
  # Keep cache files for debugging
  python snpeff_anno_main.py --keep-cache

Notes:
  - Defaults use project-preset paths suitable for the current environment.
  - Supports parallel processing to improve efficiency for large files.
  - Automatically generates detailed log files for debugging and monitoring.
        """
    )
    
    # Input options
    parser.add_argument(
        "--vcf-path", 
        type=str,
        default="/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/analysis/assoc_rvtest/results/01.rvtest_prepare/cteph_agp3k.rare.rand300x10000.all.nochr.norm.vcf.gz",
        help="Input VCF file path (Default: Project preset path)"
    )
    
    parser.add_argument(
        "--snpeff-dir", 
        type=str,
        default="/LARGE1/gr10478/platform/JHRPv4/workspace/pipeline/output/snpEff.v4.index",
        help="snpEff annotation directory path (Default: Project preset path)"
    )
    
    # Output options
    output_group = parser.add_argument_group("Output Options")
    output_group.add_argument(
        "-o", "--output", 
        type=str, 
        default=None,
        help="Output VCF file path (Default: *.snpeff.vcf.gz in working directory)"
    )
    output_group.add_argument(
        "--force", 
        action="store_true",
        default=True,
        help="Force overwrite existing output files (Default: Enabled)"
    )
    
    # Performance options
    performance_group = parser.add_argument_group("Performance Options")
    performance_group.add_argument(
        "--parallel", 
        action="store_true", 
        default=True,
        help="Enable parallel mode (Default: Enabled)"
    )
    performance_group.add_argument(
        "--no-parallel", 
        action="store_true",
        help="Disable parallel mode, use serial processing"
    )
    performance_group.add_argument(
        "--max-workers", 
        type=int, 
        default=16,
        help="Maximum parallel workers (Default: 16)"
    )
    performance_group.add_argument(
        "--threads", 
        type=int, 
        default=8,
        help="Number of threads for bcftools (Default: 8)"
    )
    performance_group.add_argument(
        "--allow-ref-alt-swap", 
        action="store_true",
        help="Allow annotation if REF/ALT are strictly swapped (Start Flip). Checks REF=ALT_ann and ALT=REF_ann."
    )
    
    # Log and Debug options
    debug_group = parser.add_argument_group("Log and Debug Options")
    debug_group.add_argument(
        "--log-file", 
        type=str, 
        default=None,
        help="Log file path (Default: snpeff_annotation.log)"
    )
    debug_group.add_argument(
        "--log-level", 
        choices=["DEBUG", "INFO", "WARNING", "ERROR"], 
        default="INFO",
        help="Log level (Default: INFO)"
    )
    debug_group.add_argument(
        "--keep-cache", 
        action="store_true",
        help="Keep intermediate cache files for debugging (Default: Delete)"
    )
    
    # Tool paths
    tool_group = parser.add_argument_group("Tool Configuration")
    tool_group.add_argument(
        "--bcftools-path",
        type=str,
        default="/home/b/b37974/bcftools/bcftools",
        help="Path to bcftools executable (Default: /home/b/b37974/bcftools/bcftools)"
    )
    
    # Version
    parser.add_argument(
        "--version", 
        action="version", 
        version="snpEff VCF Annotation Tool v1.0.0"
    )
    
    return parser


def validate_inputs(args):
    """
    Validate input files and directories.
    
    Args:
        args: Command line arguments
        
    Raises:
        SystemExit: If input is invalid
    """
    # Validate VCF file
    vcf_path = Path(args.vcf_path)
    if not vcf_path.exists():
        print(f"❌ Error: VCF file not found - {vcf_path}")
        sys.exit(1)
    
    if not str(vcf_path).endswith('.vcf.gz'):
        print(f"❌ Error: Input file must be .vcf.gz format - {vcf_path}")
        sys.exit(1)
    
    # Validate Annotation Directory
    snpeff_dir = Path(args.snpeff_dir)
    if not snpeff_dir.exists():
        print(f"❌ Error: snpEff annotation directory not found - {snpeff_dir}")
        sys.exit(1)
    
    # Check Annotation files
    tsv_files = list(snpeff_dir.glob("*.tsv.2.gz"))
    if not tsv_files:
        print(f"❌ Error: No *.tsv.2.gz files found in annotation directory - {snpeff_dir}")
        sys.exit(1)
        
    # Check bcftools
    if not Path(args.bcftools_path).exists() and shutil.which(args.bcftools_path) is None:
         # Try fallback to just 'bcftools'
         if shutil.which("bcftools"):
             print(f"⚠️ Warning: Specified bcftools not found at {args.bcftools_path}, using system 'bcftools'.")
             args.bcftools_path = "bcftools"
         else:
             print(f"❌ Error: bcftools not found at {args.bcftools_path} and not in PATH.")
             sys.exit(1)
    
    print(f"✓ Found {len(tsv_files)} annotation files")


def print_configuration(args):
    """
    Print task configuration
    
    Args:
        args: Command line arguments
    """
    print("=" * 70)
    print("snpEff VCF Annotation Tool - Configuration")
    print("=" * 70)
    print(f"Input VCF       : {args.vcf_path}")
    print(f"Annotation Dir  : {args.snpeff_dir}")
    print(f"Output File     : {args.output if args.output else 'Auto-generated'}")
    print(f"Mode            : {'Parallel' if args.parallel else 'Serial'}")
    if args.parallel:
        print(f"Max Workers     : {args.max_workers}")
    print(f"bcftools Path   : {args.bcftools_path}")
    print(f"bcftools Threads: {args.threads}")
    print(f"Force Overwrite : {args.force}")
    print(f"Keep Cache      : {args.keep_cache}")
    print(f"Log Level       : {args.log_level}")
    print("=" * 70)


def print_statistics(stats, logger, n_samples="Unknown"):
    """
    Print detailed statistics from the stats dictionary.
    """
    logger.info("="*60)
    logger.info("DETAILED VCF STATISTICS REPORT (Simultaneous Calculation)")
    logger.info("="*60)
    
    # Calculate global totals
    global_total = 0
    global_impact = Counter()
    global_effect = Counter()
    global_biotype = Counter()
    global_other = Counter()
    global_combinations = Counter()

    for chrom_key in stats:
        data = stats[chrom_key]
        # Robustly handle if data is just dict or simpler structure (backward compatibility safety)
        if "total" in data and isinstance(data["total"], int):
            global_total += data["total"]
        if "impact" in data: global_impact.update(data["impact"])
        if "effect" in data: global_effect.update(data["effect"])
        if "biotype" in data: global_biotype.update(data["biotype"])
        if "other_details" in data: global_other.update(data["other_details"])
        if "combinations" in data: global_combinations.update(data["combinations"])

    # --- Write Detailed Combination Log ---
    try:
        # Use a fixed name or timestamped? User said "generate a log file", 
        # timestamped is safer to avoid overwriting distinct runs if parallel
        combo_log_name = f"snpeff_stats_combinations_{datetime.now().strftime('%Y%m%d_%H%M%S')}.tsv"
        combo_log_path = Path.cwd() / combo_log_name
        
        with open(combo_log_path, "w") as f:
            f.write("Impact\tEffect\tBiotype\tCount\n")
            for (imp, eff, bio), count in global_combinations.most_common():
                 f.write(f"{imp}\t{eff}\t{bio}\t{count}\n")
        logger.info(f"Report: Detailed Combination Statistics saved to {combo_log_name}")
    except Exception as e:
        logger.error(f"Failed to generate combination statistics log: {e}")

    # --- Global Summary ---
    logger.info(f"Total Samples  : {n_samples}")
    logger.info(f"Total Variants : {global_total:,}")
    logger.info("-" * 40)
    
    def log_counter(title, counter, limit=None):
        logger.info(f"{title}:")
        items = counter.most_common()
        if not items:
            logger.info("  (No data)")
            return

        if limit and len(items) > limit:
            items_show = items[:limit]
            w_remaining = sum(c for _, c in items[limit:])
            items_show.append((f"Others ({len(items)-limit} types)", w_remaining))
            items = items_show
        
        # Determine width
        max_len = max([len(str(k)) for k, _ in items]) if items else 10
        for k, v in items:
            logger.info(f"  {str(k):<{max_len}} : {v:,}")
    
    log_counter("Global Impact Counts", global_impact)
    
    if global_other:
        logger.info("-" * 20)
        log_counter("Unclassified Impact (OTHER) Details", global_other, limit=10)
    
    logger.info("-" * 40)
    log_counter("Global Effect Counts (Top 20)", global_effect, limit=20)
    
    logger.info("-" * 40)
    log_counter("Global Biotype Counts (Top 20)", global_biotype, limit=20)
    
    logger.info("-" * 40)
    logger.info(f"Per-Chromosome Statistics:")
    
    # Sort chromosomes naturally
    def sort_key(k):
        try:
            k_clean = k.replace("chr", "")
            if k_clean.isdigit(): return int(k_clean)
            if k_clean == "X": return 23
            if k_clean == "Y": return 24
            if k_clean in ["M", "MT"]: return 25
            return k
        except: return k
        
    sorted_chrs = sorted(stats.keys(), key=sort_key)
    
    # Table Header
    header = f"{'CHROM':<8} {'TOTAL':<12} {'HIGH':<10} {'MODERATE':<10} {'LOW':<10} {'MODIFIER':<10}"
    logger.info(header)
    
    for chrom in sorted_chrs:
        d = stats[chrom]
        total = d.get('total', 0)
        imps = d.get('impact', {})
        
        h = imps.get('HIGH', 0)
        m = imps.get('MODERATE', 0)
        l = imps.get('LOW', 0)
        mod = imps.get('MODIFIER', 0)
        
        line = f"{chrom:<8} {total:<12,} {h:<10,} {m:<10,} {l:<10,} {mod:<10,}"
        logger.info(line)
        
    logger.info("="*60)


def main():
    """
    Main function - Execute full VCF annotation workflow
    """
    # Parse args
    parser = setup_argument_parser()
    args = parser.parse_args()
    
    # Handle conflict
    if args.no_parallel:
        args.parallel = False
    
    # Validate
    validate_inputs(args)
    
    # Print config
    print_configuration(args)
    
    try:
        total_start_time = datetime.now()
        
        # Step 1: Add chr prefix
        print("\n🔄 Step 1: Checking and adding 'chr' prefix...")
        step1_start = datetime.now()
        
        vcf_with_chr = add_chr_prefix_to_vcf(
            vcf_path=args.vcf_path,
            bcftools_path=args.bcftools_path,
            force=False 
        )
        
        step1_duration = datetime.now() - step1_start
        print(f"✓ Prefix Check Complete, Duration: {step1_duration}")
        print(f"  Processed VCF: {vcf_with_chr}")
        
        # Step 2: snpEff annotation (with simultaneous stats)
        print("\n🔄 Step 2: Executing snpEff functional annotation (and stats)...")
        step2_start = datetime.now()
        
        annotated_vcf, collected_stats = annotate_vcf_with_snpeff_tsv(
            vcf_path=vcf_with_chr,
            snpeff_tsv_dir=args.snpeff_dir,
            bcftools_path=args.bcftools_path,
            out_path=args.output,
            threads=args.threads,
            parallel=args.parallel,
            max_workers=args.max_workers,
            force=args.force,
            remove_cache=not args.keep_cache,
            log_file=args.log_file,
            log_level=args.log_level,
        )
        
        step2_duration = datetime.now() - step2_start
        total_duration = datetime.now() - total_start_time
        
        # Determine logger
        log_file_path = args.log_file if args.log_file else "snpeff_annotation.log"
        logger = logging.getLogger('snpeff_main_stats')
        logger.setLevel(logging.INFO)
        fh = logging.FileHandler(log_file_path)
        fh.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
        logger.addHandler(fh)
        sh = logging.StreamHandler(sys.stdout)
        sh.setFormatter(logging.Formatter('%(message)s'))
        logger.addHandler(sh)

        # Step 3: Print Statistics (using collected data)
        print("\n🔄 Step 3: Reporting Statistics...")
        
        # Optional: Get sample count (still need to query as stats collection didn't count samples)
        # Sample count is constant, so just quick query
        try:
            res = subprocess.run([args.bcftools_path, "query", "-l", str(annotated_vcf)], 
                               stdout=subprocess.PIPE, text=True)
            n_samples = len(res.stdout.splitlines())
        except:
            n_samples = "Unknown"
            
        print_statistics(collected_stats, logger, n_samples)

        # Summary
        print("\n" + "=" * 70)
        print("🎉 Task Completion Summary")
        print("=" * 70)
        print(f"✓ Total Duration : {total_duration}")
        print(f"  - Step 1 (Pre): {step1_duration}")
        print(f"  - Step 2 (Ann): {step2_duration}")
        print(f"✓ Final Output   : {annotated_vcf}")
        
        if annotated_vcf.exists():
            file_size = annotated_vcf.stat().st_size / (1024 * 1024)  # MB
            print(f"✓ Output Size    : {file_size:.2f} MB")
            
            tbi_file = Path(str(annotated_vcf) + ".tbi")
            if tbi_file.exists():
                print(f"✓ Index Generated: {tbi_file}")
        
        print(f"✓ Detailed Log   : {log_file_path}")
        print("=" * 70)
        print("🎉 All tasks completed successfully!")
        
    except KeyboardInterrupt:
        print("\n\n⚠️  User Interrupted")
        sys.exit(1)
        
    except Exception as e:
        print(f"\n\n❌ Execution Failed: {e}")
        log_path = args.log_file if args.log_file else "snpeff_annotation.log"
        print(f"Check log for details: {log_path}")
        sys.exit(1)


if __name__ == "__main__":
    main()


