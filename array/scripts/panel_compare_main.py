#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script Name: panel_compare_main.py
==================================

Overview:
This script orchestrates a three-step pipeline:
1) Runs PLINK2 variant QC summary (`run_plink2_variant_qc`) to generate `*.variant_qc_summary.tsv`.
2) Appends ToMMo panel annotations (VCF) to the summary: `IN_TOMMO, TOMMO_AAF, TOMMO_FILTER` (`run_plink2_variant_qc_with_tommo`).
3) Generates a multi-page PDF comparison report (`plot_tommo_panel_compare_pdf`).

Use Case:
- Post-QC of WGS/WES/Array data, aligning with ToMMo (or other external VCF panels) to quickly produce statistical tables and visualization reports.

Inputs & Outputs:
- Input:
  - `--bed_prefix`: PLINK binary genotype prefix (.bed/.bim/.fam).
  - `--tommo_vcf_path`: ToMMo panel bgzip-compressed VCF path (must have .tbi index).
- Output:
  - `{output_prefix}.variant_qc_summary.tsv`
  - `{output_prefix}.variant_qc_summary.variant_qc_with_tommo.tsv`
  - `{output_prefix}.variant_qc_with_tommo.panel_compare.pdf`

Performance Parameters:
- `--threads`: Number of threads passed to `run_plink2_variant_qc` and `bcftools view --threads`.
- `--chunk_size`: Chunk size (rows) for streaming merge of ToMMo annotations.
- `--max_workers`: Max concurrency for chromosome-parallel processing.

Author: ZHAO TIE
"""

import argparse
import os
import sys
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend for PDF/PNG generation
import matplotlib.pyplot as plt  # noqa: F401

# Add current script directory to sys.path to ensure local modules can be imported
HERE = os.path.abspath(os.path.dirname(__file__))
if HERE not in sys.path:
    sys.path.append(HERE)

import importlib
import variant_qc_calculator
import panel_compare_tools

# Reload for development (can be removed in production)
importlib.reload(variant_qc_calculator)
importlib.reload(panel_compare_tools)

from variant_qc_calculator import run_plink2_variant_qc
from panel_compare_tools import (
    run_plink2_variant_qc_with_tommo,
    plot_tommo_panel_compare_pdf,
)


def parse_args() -> argparse.Namespace:
    """Parse command line arguments."""
    p = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=(
            "Integrated script for PLINK2 Variant QC -> ToMMo Annotation -> PDF Report.\n"
            "Author: ZHAO TIE"
        ),
    )
    p.add_argument(
        "--bed_prefix", required=True,
        help="PLINK binary genotype prefix (same prefix for .bed/.bim/.fam)",
    )
    p.add_argument(
        "--output_prefix", required=True,
        help="Output file prefix (used for filenames in all three steps)",
    )
    p.add_argument(
        "--threads", type=int, default=16,
        help="Number of threads: used for PLINK2 and bcftools view --threads",
    )
    p.add_argument(
        "--tommo_vcf_path", required=True,
        help="ToMMo panel VCF (.vcf.gz), must have accompanying .tbi index",
    )
    p.add_argument(
        "--chunk_size", type=int, default=500_000,
        help="Chunk size (rows) for merging ToMMo annotations",
    )
    p.add_argument(
        "--max_workers", type=int, default=4,
        help="Max concurrency for chromosome-parallel processing",
    )
    p.add_argument(
        "--skip_pdf", action="store_true",
        help="Only generate TSV with ToMMo annotations, skip PDF generation",
    )
    p.add_argument(
        "--regions_chunk_lines", type=int, default=50_000,
        help="Chunk size (lines) for splitting chromosome regions, used for progress control",
    )
    p.add_argument(
        "--tabix_path", default="/usr/bin/tabix",
        help="Path to tabix executable",
    )
    p.add_argument(
        "--grouping_metric", default="MAF_ALL",
        help="Metric column name used for grouping (default MAF_ALL)",
    )
    p.add_argument(
        "--grouping_threshold", type=float, default=0.05,
        help="Grouping threshold (default 0.05). If set to -1, no grouping (All Variants)",
    )
    return p.parse_args()


def _progress(msg: str):
    print(f"[panel_compare_main] {msg}", file=sys.stderr, flush=True)


def main():
    args = parse_args()

    bed_prefix     = args.bed_prefix
    output_prefix  = args.output_prefix
    threads        = int(args.threads)
    tommo_vcf_path = args.tommo_vcf_path
    chunk_size     = int(args.chunk_size)
    max_workers    = int(args.max_workers)
    regions_chunk_lines = int(args.regions_chunk_lines)
    tabix_path     = args.tabix_path
    grouping_metric = args.grouping_metric
    grouping_threshold = args.grouping_threshold if args.grouping_threshold >= 0 else None

    _progress("Step1: Running PLINK2 Variant QC Summary ...")
    variant_qc_summary = run_plink2_variant_qc(
        bed_prefix=bed_prefix,
        output_prefix=output_prefix,
        threads=threads,
        tabix_path=tabix_path,
    )
    _progress(f"Output: {variant_qc_summary}")

    _progress("Step2: Batch Annotating with ToMMo (IN_TOMMO/TOMMO_AAF/TOMMO_FILTER) ...")
    
    # Calculate threads per worker to avoid oversubscription
    # Total threads = max_workers * threads_per_worker
    # We want Total threads <= args.threads
    if max_workers > 0:
        bcftools_threads = max(1, threads // max_workers)
    else:
        bcftools_threads = 1
        
    _progress(f"Configuration: Total Threads={threads}, Max Workers={max_workers}, Threads per Worker={bcftools_threads}")

    out_tsv = run_plink2_variant_qc_with_tommo(
        variant_qc_summary=variant_qc_summary,
        tommo_vcf_path=tommo_vcf_path,
        threads=bcftools_threads, # Passed to bcftools view --threads (per worker)
        chunk_size=chunk_size,   # Streaming merge, memory friendly
        max_workers=max_workers, # Chromosome parallelism (I/O friendly 2~4)
        regions_chunk_lines=regions_chunk_lines,
        tabix_path=tabix_path,
    )
    _progress(f"Output: {out_tsv}")

    if args.skip_pdf:
        _progress("Step3: Skipping PDF generation (--skip_pdf is on)")
        print(out_tsv)
        return

    _progress("Step3: Generating Panel Comparison PDF Report ...")
    pdf_path = plot_tommo_panel_compare_pdf(
        variant_qc_with_tommo=out_tsv,
        split_metric=grouping_metric,
        split_threshold=grouping_threshold,
    )
    _progress(f"Output: {pdf_path}")

    # Terminal output for pipeline collection
    print(out_tsv)
    print("Saved:", pdf_path)


if __name__ == "__main__":
    main()
