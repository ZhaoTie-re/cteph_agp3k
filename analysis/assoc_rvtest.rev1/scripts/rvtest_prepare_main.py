#!/usr/bin/env python3
"""
rvtest_prepare_main.py - RVTest Data Preprocessing Main Program

This script is a command-line tool designed to automate all data preprocessing steps 
required for RVTest (Rare Variant Association Testing).

Key Features:
1. Genotype Data Conversion: Converts PLINK format (.bed/.bim/.fam) to standardized VCF format.
2. Phenotype/Covariate Reformatting: reformats phenotype and covariate files to match RVTest input requirements.
3. Gene Annotation File Processing: Removes chromosome prefixes from refFlat files to ensure format consistency.

Design Principles:
- Full command-line argument support.
- Detailed logging and error handling.
- Modular design for easy maintenance and extension.
- Support for parallel processing to improve efficiency.

Author: ZHAO TIE
Created: Oct 8, 2025
Updated: Jan 14, 2026

Usage Example:
    python rvtest_prepare_main.py \\
        --bed-prefix /path/to/plink_data \\
        --pheno-path /path/to/phenotype.csv \\
        --covar-path /path/to/covariate.csv \\
        --refflat-path /path/to/refFlat.txt.gz \\
        --threads 8 \\
        --no-norm
"""

import sys
import os
import argparse
import logging
from datetime import datetime
from typing import Tuple

# Add module path dynamically
current_dir = os.path.dirname(os.path.abspath(__file__))
if current_dir not in sys.path:
    sys.path.insert(0, current_dir)

import rvtest_prepare_tools


def setup_logging(verbose: bool = False, log_file: str = None) -> None:
    """
    Configure the logging system.
    
    Args:
        verbose (bool): Whether to enable verbose logging (DEBUG level).
        log_file (str): Path to the log file. If None, generates a timestamped name.
    """
    log_level = logging.DEBUG if verbose else logging.INFO
    log_format = '%(asctime)s - %(levelname)s - %(message)s'
    
    if not log_file:
         log_file = f'rvtest_prepare_{datetime.now().strftime("%Y%m%d_%H%M%S")}.log'

    logging.basicConfig(
        level=log_level,
        format=log_format,
        handlers=[
            logging.StreamHandler(sys.stdout),
            logging.FileHandler(log_file)
        ]
    )


def parse_arguments() -> argparse.Namespace:
    """
    Parse command-line arguments.
    
    Returns:
        argparse.Namespace: Object containing parsed command-line arguments.
    """
    parser = argparse.ArgumentParser(
        description='RVTest Data Preprocessing Tool - Automates genotype, phenotype, and annotation data processing',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Usage Examples:
  # Basic usage - Process all required files
  %(prog)s --bed-prefix /path/to/data \\
           --pheno-path /path/to/pheno.csv \\
           --covar-path /path/to/covar.csv \\
           --refflat-path /path/to/refFlat.txt.gz

  # Advanced usage - Custom threads, output format, and normalization control
  %(prog)s --bed-prefix /path/to/data \\
           --pheno-path /path/to/pheno.csv \\
           --covar-path /path/to/covar.csv \\
           --refflat-path /path/to/refFlat.txt.gz \\
           --threads 16 \\
           --keep-chr-prefix \\
           --snps-only \\
           --norm \\
           --verbose

Notes:
  - All input file paths must be absolute paths.
  - Ensure sufficient disk space for intermediate and final output files.
  - Using more threads is recommended for large datasets.
        """
    )
    
    # Required arguments group
    required_group = parser.add_argument_group('Required Arguments')
    required_group.add_argument(
        '--bed-prefix',
        type=str,
        required=True,
        help='Prefix for PLINK format data files (excluding .bed/.bim/.fam suffixes)'
    )
    required_group.add_argument(
        '--pheno-path',
        type=str,
        required=True,
        help='Path to phenotype data file (CSV format)'
    )
    required_group.add_argument(
        '--covar-path',
        type=str,
        required=True,
        help='Path to covariate data file (CSV format)'
    )
    required_group.add_argument(
        '--refflat-path',
        type=str,
        required=True,
        help='Path to refFlat gene annotation file (.txt.gz format)'
    )
    
    # Optional arguments group
    optional_group = parser.add_argument_group('Optional Arguments')
    optional_group.add_argument(
        '--threads',
        type=int,
        default=8,
        help='Number of threads for parallel processing (Default: 8)'
    )
    optional_group.add_argument(
        '--snps-only',
        action='store_true',
        help='Process only A/C/G/T standard SNP variants (Default: Process all variant types)'
    )
    optional_group.add_argument(
        '--keep-chr-prefix',
        action='store_true',
        help='Keep "chr" prefix in VCF file (Default: Remove prefix)'
    )
    optional_group.add_argument(
        '--output-dir',
        type=str,
        default=None,
        help='Path to output directory (Default: Current working directory)'
    )
    optional_group.add_argument(
        '--norm',
        action='store_true',
        help='Perform bcftools normalization (Default: Do NOT normalize)'
    )
    optional_group.add_argument(
        '--verbose',
        action='store_true',
        help='Enable verbose logging output'
    )
    optional_group.add_argument(
        '--log-file',
        type=str,
        default=None,
        help='Log file path (Default: rvtest_prepare_YYYYMMDD_HHMMSS.log)'
    )
    
    return parser.parse_args()


def validate_input_files(args: argparse.Namespace) -> None:
    """
    Validate the existence and readability of input files.
    
    Args:
        args (argparse.Namespace): Command-line arguments object.
        
    Raises:
        FileNotFoundError: If a required input file does not exist.
        PermissionError: If file permissions are insufficient.
    """
    logging.info("Validating input files...")
    
    # Check PLINK files
    plink_extensions = ['.bed', '.bim', '.fam']
    for ext in plink_extensions:
        file_path = args.bed_prefix + ext
        if not os.path.exists(file_path):
            raise FileNotFoundError(f"PLINK file not found: {file_path}")
        if not os.access(file_path, os.R_OK):
            raise PermissionError(f"Cannot read file: {file_path}")
    
    # Check phenotype and covariate files
    for file_path in [args.pheno_path, args.covar_path, args.refflat_path]:
        if not os.path.exists(file_path):
            raise FileNotFoundError(f"Input file not found: {file_path}")
        if not os.access(file_path, os.R_OK):
            raise PermissionError(f"Cannot read file: {file_path}")
    
    logging.info("✓ All input files validated successfully.")


def process_plink_to_vcf(args: argparse.Namespace) -> str:
    """
    Execute PLINK to VCF conversion.
    
    Args:
        args (argparse.Namespace): Command-line arguments object.
        
    Returns:
        str: Path to the generated VCF file.
    """
    logging.info("="*60)
    logging.info("Step 1/3: Starting PLINK to VCF conversion")
    logging.info("="*60)
    logging.info(f"Input file prefix: {args.bed_prefix}")
    logging.info(f"Threads: {args.threads}")
    logging.info(f"SNPs only: {'Yes' if args.snps_only else 'No'}")
    logging.info(f"Keep chr prefix: {'Yes' if args.keep_chr_prefix else 'No'}")
    logging.info(f"Perform normalization: {'Yes' if args.norm else 'No'}")
    
    try:
        vcf_path = rvtest_prepare_tools.plink_to_vcf_raw(
            bed_prefix=args.bed_prefix,
            threads=args.threads,
            snps_only_just_acgt=args.snps_only,
            keep_chr_prefix=args.keep_chr_prefix,
            perform_norm=args.norm
        )
        logging.info(f"✓ VCF conversion completed. Output file: {vcf_path}")
        return vcf_path
    except Exception as e:
        logging.error(f"✗ VCF conversion failed: {str(e)}")
        raise


def process_pheno_covar(args: argparse.Namespace) -> Tuple[str, str]:
    """
    Execute reformatting of phenotype and covariate files.
    
    Args:
        args (argparse.Namespace): Command-line arguments object.
        
    Returns:
        Tuple[str, str]: (New phenotype file path, New covariate file path)
    """
    logging.info("="*60)
    logging.info("Step 2/3: Starting Phenotype and Covariate Reformatting")
    logging.info("="*60)
    logging.info(f"Phenotype file: {args.pheno_path}")
    logging.info(f"Covariate file: {args.covar_path}")
    
    try:
        new_pheno, new_covar = rvtest_prepare_tools.reformat_pheno_covar(
            pheno_path=args.pheno_path,
            covar_path=args.covar_path,
            out_dir=args.output_dir,
            out_sep="\t",  # RVTest recommends tab separator
            force=True,
        )
        logging.info(f"✓ Phenotype file reformatting completed: {new_pheno}")
        logging.info(f"✓ Covariate file reformatting completed: {new_covar}")
        return new_pheno, new_covar
    except Exception as e:
        logging.error(f"✗ Phenotype/Covariate file processing failed: {str(e)}")
        raise


def process_refflat(args: argparse.Namespace) -> str:
    """
    Execute refFlat annotation file processing.
    
    Args:
        args (argparse.Namespace): Command-line arguments object.
        
    Returns:
        str: Path to the processed refFlat file.
    """
    logging.info("="*60)
    logging.info("Step 3/3: Starting refFlat Annotation File Processing")
    logging.info("="*60)
    logging.info(f"Input file: {args.refflat_path}")
    
    try:
        processed_refflat = rvtest_prepare_tools.reformat_refflat_remove_chr(
            refflat_path=args.refflat_path,
            out_dir=args.output_dir
        )
        logging.info(f"✓ refFlat file processing completed: {processed_refflat}")
        return processed_refflat
    except Exception as e:
        logging.error(f"✗ refFlat file processing failed: {str(e)}")
        raise


def print_summary(vcf_path: str, pheno_path: str, covar_path: str, refflat_path: str) -> None:
    """
    Print processing result summary.
    
    Args:
        vcf_path (str): VCF file path.
        pheno_path (str): Phenotype file path.
        covar_path (str): Covariate file path.
        refflat_path (str): refFlat file path.
    """
    logging.info("="*60)
    logging.info("RVTest Data Preprocessing Completed!")
    logging.info("="*60)
    logging.info("Processing Summary:")
    logging.info(f"  Genotype Data (VCF):     {vcf_path}")
    logging.info(f"  Phenotype Data:          {pheno_path}")
    logging.info(f"  Covariate Data:          {covar_path}")
    logging.info(f"  Gene Annotation (refFlat): {refflat_path}")
    logging.info("")
    logging.info("You can now use these files for RVTest Rare Variant Association Analysis.")
    logging.info("It is recommended to check the output files to confirm correct formatting.")


def main() -> None:
    """
    Main program entry point.
    
    Executes the complete RVTest data preprocessing workflow:
    1. Parse command-line arguments.
    2. Validate input files.
    3. Convert PLINK to VCF.
    4. Reformat phenotype and covariates.
    5. Process refFlat annotation file.
    6. Output processing summary.
    """
    # Parse command-line arguments
    args = parse_arguments()
    
    # Setup logging
    setup_logging(args.verbose, args.log_file)
    
    # Print program start info
    logging.info("="*60)
    logging.info("RVTest Data Preprocessing Tool Started")
    logging.info(f"Start Time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    logging.info("="*60)
    
    try:
        # Validate input files
        validate_input_files(args)
        
        # Step 1: PLINK to VCF conversion
        vcf_path = process_plink_to_vcf(args)
        
        # Step 2: Phenotype and Covariate reformatting
        pheno_path, covar_path = process_pheno_covar(args)
        
        # Step 3: refFlat file processing
        refflat_path = process_refflat(args)
        
        # Print result summary
        print_summary(vcf_path, pheno_path, covar_path, refflat_path)
        
        logging.info("="*60)
        logging.info("Program execution completed successfully!")
        logging.info(f"End Time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        logging.info("="*60)
        
    except Exception as e:
        logging.error("="*60)
        logging.error("An error occurred during program execution!")
        logging.error(f"Error Message: {str(e)}")
        logging.error("="*60)
        sys.exit(1)


if __name__ == "__main__":
    main()
