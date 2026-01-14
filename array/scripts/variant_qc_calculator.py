"""
Module Name: variant_qc_calculator.py
=====================================

Overview:
This module provides tools for calculating variant-level Quality Control (QC) metrics on large-scale genotype data
(PLINK binary format: .bed/.bim/.fam). It automates the aggregation of QC statistics and supports filtering of
low-quality variants to generate cleaned PLINK datasets. The workflow uses plink2 as the core engine, supplemented
by pandas for chunked streaming and multiprocessing, ensuring both reproducibility and performance.

Key Features:
- Standard variant-level QC statistics (AAF/MAF, Missing Rate, HWE) for GWAS/Candidate gene analysis.
- Rapid calculation of statistics for Case/Control groups.
- Generation of traceable QC summary tables and exclusion lists.

Inputs:
- PLINK binary file prefix (.bed/.bim/.fam).
- Phenotype in .fam file must be encoded as 1=Control, 2=Case.

Outputs:
1. Variant QC Summary Table (`*.variant_qc_summary.tsv`):
   - VARIANT_ID, MAF, VMISS, CASE/CTRL AAF, CASE/CTRL HWE, etc.
2. Low Quality Variant Flags (`*.with_flags.tsv`) and Exclusion List (`*.variant_ids.tsv`).
3. Filtered Dataset (new PLINK binary files).

Dependencies:
- Python 3.9+ (pandas, numpy)
- plink2 executable

Author: ZHAO TIE
"""

import subprocess
import os
import tempfile
import csv
import math
import uuid
import gzip
import shutil
from typing import Optional, Tuple
import concurrent.futures
import pandas as pd
import numpy as np

_global_lookup = {}

def init_globals_for_chunk(vmiss_dict, case_vmiss_dict, ctrl_vmiss_dict, 
                           case_aaf_dict, ctrl_aaf_dict, 
                           case_hwe_dict, ctrl_hwe_dict, all_hwe_dict):
    """Initialize global lookup dictionaries for worker processes."""
    global _global_lookup
    _global_lookup['vmiss_dict'] = vmiss_dict
    _global_lookup['case_vmiss_dict'] = case_vmiss_dict
    _global_lookup['ctrl_vmiss_dict'] = ctrl_vmiss_dict
    _global_lookup['case_aaf_dict'] = case_aaf_dict
    _global_lookup['ctrl_aaf_dict'] = ctrl_aaf_dict
    _global_lookup['case_hwe_dict'] = case_hwe_dict
    _global_lookup['ctrl_hwe_dict'] = ctrl_hwe_dict
    _global_lookup['all_hwe_dict'] = all_hwe_dict

def process_chunk(chunk, idx, tmpdir):
    """
    Process a chunk of variant data to extract QC metrics and write to a temporary file.
    """
    vmiss_dict = _global_lookup.get('vmiss_dict', {})
    case_vmiss_dict = _global_lookup.get('case_vmiss_dict', {})
    ctrl_vmiss_dict = _global_lookup.get('ctrl_vmiss_dict', {})
    case_aaf_dict = _global_lookup.get('case_aaf_dict', {})
    ctrl_aaf_dict = _global_lookup.get('ctrl_aaf_dict', {})
    case_hwe_dict = _global_lookup.get('case_hwe_dict', {})
    ctrl_hwe_dict = _global_lookup.get('ctrl_hwe_dict', {})
    all_hwe_dict = _global_lookup.get('all_hwe_dict', {})

    tmp_output = os.path.join(tmpdir, f"chunk_{idx}_{uuid.uuid4().hex}.tsv")
    with open(tmp_output, "w", newline="") as fout:
        writer = csv.writer(fout, delimiter="\t")
        for _, row in chunk.iterrows():
            vid = row["ID"]  # Variant ID
            
            # Parse ID (chrom:pos:ref:alt)
            try:
                parts = vid.split(':')
                if len(parts) >= 4:
                    chrom = parts[0]
                    pos = parts[1]
                    ref = parts[2]
                    alt = parts[3]
                else:
                    # Fallback if ID format is unexpected
                    chrom, pos, ref, alt = ".", ".", ".", "."
            except Exception:
                chrom, pos, ref, alt = ".", ".", ".", "."

            # AAF / MAF
            aaf_all = row["ALT_FREQS"]
            maf_all = min(aaf_all, 1 - aaf_all) if pd.notnull(aaf_all) else float("nan")
            
            case_aaf = case_aaf_dict.get(vid, float("nan"))
            maf_case = min(case_aaf, 1 - case_aaf) if pd.notnull(case_aaf) else float("nan")
            
            ctrl_aaf = ctrl_aaf_dict.get(vid, float("nan"))
            maf_ctrl = min(ctrl_aaf, 1 - ctrl_aaf) if pd.notnull(ctrl_aaf) else float("nan")

            # VMISS
            vmiss_all = vmiss_dict.get(vid, float("nan"))
            vmiss_case = case_vmiss_dict.get(vid, float("nan"))
            vmiss_ctrl = ctrl_vmiss_dict.get(vid, float("nan"))

            # HWE
            hwe_all = all_hwe_dict.get(vid, float("nan"))
            hwe_case = case_hwe_dict.get(vid, float("nan"))
            hwe_ctrl = ctrl_hwe_dict.get(vid, float("nan"))

            # Write columns: 
            # #CHROM, POS, ID, REF, ALT, 
            # MAF_ALL, MAF_CASE, MAF_CTRL, 
            # VMISS_CASE, VMISS_CTRL, VMISS_ALL, 
            # AAF_CASE, AAF_CTRL, AAF_ALL, 
            # HWE_CASE, HWE_CTRL, HWE_ALL
            writer.writerow([
                chrom, pos, vid, ref, alt,
                maf_all, maf_case, maf_ctrl,
                vmiss_case, vmiss_ctrl, vmiss_all,
                case_aaf, ctrl_aaf, aaf_all,
                hwe_case, hwe_ctrl, hwe_all
            ])
    return tmp_output

def run_plink2_variant_qc(
    bed_prefix: str,
    tmpdir: str = "/tmp/variant_qc",
    plink2_path: str = "/home/b/b37974/plink2",
    threads: int = 8,
    output_prefix: str = "cteph_agp3k",
    verbose: bool = True,
    tabix_path: str = "/usr/bin/tabix"
) -> str:
    """
    Run plink2 to calculate variant-level QC metrics for a PLINK genotype dataset.

    Parameters
    ----------
    bed_prefix : str
        Prefix of the PLINK binary files (.bed/.bim/.fam).
    tmpdir : str
        Directory for temporary files.
    plink2_path : str
        Path to the plink2 executable.
    threads : int
        Number of threads for plink2.
    output_prefix : str
        Prefix for the output files.
    verbose : bool
        Whether to print progress messages.
    tabix_path : str
        Path to the tabix executable for indexing.

    Returns
    -------
    str
        Path to the generated variant QC summary file (*.variant_qc_summary.tsv.gz).
    """

    # Ensure tmpdir exists
    if tmpdir is None:
        tmpdir = tempfile.mkdtemp()
    else:
        os.makedirs(tmpdir, exist_ok=True)
    if verbose:
        print(f"[INFO] Using temporary directory: {tmpdir}")

    def run_cmd(cmd, desc: Optional[str] = None):
        if verbose and desc:
            print(f"[INFO] Running: {desc}")
        try:
            subprocess.run(cmd, check=True)
        except subprocess.CalledProcessError as e:
            print(f"[ERROR] Command failed: {' '.join(cmd)}")
            raise e

    # Step 1: freq + vmiss + hardy (ALL)
    out_all = os.path.join(tmpdir, "all_samples")
    run_cmd([
        plink2_path, "--bfile", bed_prefix, "--threads", str(threads),
        "--freq", "--missing", "--hardy", "--out", out_all
    ], "All Samples AAF + VMISS + HWE")

    # Step 2: Split IIDs by Phenotype
    fam_df = pd.read_csv(f"{bed_prefix}.fam", sep=r"\s+", header=None)
    fam_df.columns = ["FID", "IID", "PID", "MID", "SEX", "PHENO"]
    case_iids = fam_df[fam_df["PHENO"] == 2][["FID", "IID"]]
    ctrl_iids = fam_df[fam_df["PHENO"] == 1][["FID", "IID"]]

    case_iid_path = os.path.join(tmpdir, "case_iids.txt")
    ctrl_iid_path = os.path.join(tmpdir, "ctrl_iids.txt")
    case_iids.to_csv(case_iid_path, sep="\t", index=False, header=False)
    ctrl_iids.to_csv(ctrl_iid_path, sep="\t", index=False, header=False)

    # Step 3: Grouped freq + HWE + missing
    out_case = os.path.join(tmpdir, "case")
    out_ctrl = os.path.join(tmpdir, "ctrl")
    
    # Case
    run_cmd([
        plink2_path, "--bfile", bed_prefix, "--keep", case_iid_path,
        "--threads", str(threads), "--freq", "--hardy", "--missing", "--out", out_case
    ], "Case AAF + HWE + VMISS")

    # Control
    run_cmd([
        plink2_path, "--bfile", bed_prefix, "--keep", ctrl_iid_path,
        "--threads", str(threads), "--freq", "--hardy", "--missing", "--out", out_ctrl
    ], "Control AAF + HWE + VMISS")

    # Step 4: Load lookup tables
    if verbose:
        print("[INFO] Loading auxiliary statistics tables...")

    # VMISS
    vmiss_dict = dict(pd.read_csv(out_all + ".vmiss", sep=r"\s+")[["ID", "F_MISS"]].values)
    case_vmiss_dict = dict(pd.read_csv(out_case + ".vmiss", sep=r"\s+")[["ID", "F_MISS"]].values)
    ctrl_vmiss_dict = dict(pd.read_csv(out_ctrl + ".vmiss", sep=r"\s+")[["ID", "F_MISS"]].values)

    # AAF
    case_aaf_dict = dict(pd.read_csv(out_case + ".afreq", sep=r"\s+")[["ID", "ALT_FREQS"]].values)
    ctrl_aaf_dict = dict(pd.read_csv(out_ctrl + ".afreq", sep=r"\s+")[["ID", "ALT_FREQS"]].values)

    # HWE
    hwe_all_df = pd.read_csv(out_all + ".hardy", sep=r"\s+")
    hwe_case_df = pd.read_csv(out_case + ".hardy", sep=r"\s+")
    hwe_ctrl_df = pd.read_csv(out_ctrl + ".hardy", sep=r"\s+")
    
    all_hwe_dict = dict(hwe_all_df[["ID", "P"]].values)
    case_hwe_dict = dict(hwe_case_df[["ID", "P"]].values)
    ctrl_hwe_dict = dict(hwe_ctrl_df[["ID", "P"]].values)

    # Step 5: Parallel chunk processing
    raw_output_file = os.path.join(tmpdir, "raw_variant_qc_summary.tsv")
    
    chunk_files = []
    reader = pd.read_csv(out_all + ".afreq", sep=r"\s+", chunksize=100000)
    
    with concurrent.futures.ProcessPoolExecutor(
        max_workers=4, 
        initializer=init_globals_for_chunk,
        initargs=(vmiss_dict, case_vmiss_dict, ctrl_vmiss_dict, 
                  case_aaf_dict, ctrl_aaf_dict, 
                  case_hwe_dict, ctrl_hwe_dict, all_hwe_dict)
    ) as executor:
        futures = []
        for i, chunk in enumerate(reader):
            if verbose:
                print(f"[INFO] Submitting chunk {i + 1} ...")
            futures.append(
                executor.submit(
                    process_chunk,
                    chunk, i, tmpdir
                )
            )
        for i, future in enumerate(futures):
            chunk_file = future.result()
            chunk_files.append(chunk_file)

    # Merge chunk files
    with open(raw_output_file, "w", newline="") as fout:
        writer = csv.writer(fout, delimiter="\t")
        # Write header
        writer.writerow([
            "#CHROM", "POS", "ID", "REF", "ALT",
            "MAF_ALL", "MAF_CASE", "MAF_CTRL",
            "VMISS_CASE", "VMISS_CTRL", "VMISS_ALL",
            "AAF_CASE", "AAF_CTRL", "AAF_ALL",
            "HWE_CASE", "HWE_CTRL", "HWE_ALL"
        ])
        # Concatenate all chunk results
        for chunk_file in chunk_files:
            with open(chunk_file, "r") as fin:
                for line in fin:
                    fout.write(line)

    # Step 6: bgzip and tabix
    final_output_file = output_prefix + ".variant_qc_summary.tsv.gz"
    
    # Find bgzip (assume it's in the same dir as tabix or in PATH)
    bgzip_path = shutil.which("bgzip")
    if not bgzip_path and tabix_path:
        potential_bgzip = os.path.join(os.path.dirname(tabix_path), "bgzip")
        if os.path.exists(potential_bgzip):
            bgzip_path = potential_bgzip
    
    if not bgzip_path:
        print("[WARNING] bgzip not found, using gzip instead (no indexing possible)")
        # Fallback to gzip if bgzip is missing (though user requested tabix)
        with open(raw_output_file, 'rb') as f_in:
            with gzip.open(final_output_file, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
    else:
        # Use bgzip
        with open(raw_output_file, "rb") as f_in, open(final_output_file, "wb") as f_out:
             subprocess.run([bgzip_path, "-c"], stdin=f_in, stdout=f_out, check=True)
        
        # Index with tabix
        # -s 1 (#CHROM), -b 2 (POS), -e 2 (POS)
        if tabix_path and os.path.exists(tabix_path):
            run_cmd([tabix_path, "-s", "1", "-b", "2", "-e", "2", "-f", final_output_file], "Tabix Indexing")
        else:
            print(f"[WARNING] Tabix path invalid: {tabix_path}")

    if verbose:
        print(f"[INFO] Output completed: {final_output_file}")
    
    return final_output_file


def extract_maf0_or_vmiss1_variants_streaming(input_file: str,
                                               base_prefix: str = "maf0_or_vmiss1_variants",
                                               chunksize: int = 100000) -> Tuple[str, str]:
    """
    Extract low-quality variants (MAF=0, VMISS=1, or NA values) from the QC summary file using streaming.
    
    Parameters
    ----------
    input_file : str
        Path to the variant QC summary file (TSV).
    base_prefix : str
        Prefix for the output files (default: "maf0_or_vmiss1_variants").
    chunksize : int
        Number of rows per chunk for reading.

    Returns
    -------
    Tuple[str, str]
        - Path to the TSV file with flags (VARIANT_ID, MAF0_FLAG, VMISS1_FLAG, MAF_NA_FLAG, VMISS_NA_FLAG, CASE_AAF_NA_FLAG, CTRL_AAF_NA_FLAG).
        - Path to the variant ID list file (for plink2 --exclude).
    """
    input_dir = os.path.dirname(os.path.abspath(input_file))
    full_info_path = os.path.join(input_dir, f"{base_prefix}.with_flags.tsv")
    variant_only_path = os.path.join(input_dir, f"{base_prefix}.variant_ids.tsv")

    with open(full_info_path, 'w') as f_full, open(variant_only_path, 'w') as f_ids:
        # Write header, adding MAF_NA_FLAG, VMISS_NA_FLAG, CASE_AAF_NA_FLAG, CTRL_AAF_NA_FLAG
        f_full.write("VARIANT_ID\tMAF0_FLAG\tVMISS1_FLAG\tMAF_NA_FLAG\tVMISS_NA_FLAG\tCASE_AAF_NA_FLAG\tCTRL_AAF_NA_FLAG\n")

        reader = pd.read_csv(
            input_file,
            sep="\t",
            usecols=["VARIANT_ID", "MAF", "VMISS", "CASE_AAF", "CTRL_AAF"],
            dtype={"VARIANT_ID": str},
            chunksize=chunksize
        )

        for chunk in reader:
            if chunk.empty:
                continue

            for row in chunk.itertuples(index=False):
                try:
                    maf_val = float(row.MAF)
                except (ValueError, TypeError):
                    maf_val = float("nan")
                try:
                    vmiss_val = float(row.VMISS)
                except (ValueError, TypeError):
                    vmiss_val = float("nan")
                try:
                    case_aaf_val = float(row.CASE_AAF)
                except (ValueError, TypeError):
                    case_aaf_val = float("nan")
                try:
                    ctrl_aaf_val = float(row.CTRL_AAF)
                except (ValueError, TypeError):
                    ctrl_aaf_val = float("nan")
                is_maf0 = maf_val == 0.0 if pd.notnull(maf_val) else False
                is_vmiss1 = vmiss_val == 1.0 if pd.notnull(vmiss_val) else False
                is_maf_na = pd.isna(maf_val)
                is_vmiss_na = pd.isna(vmiss_val)
                is_case_aaf_na = pd.isna(case_aaf_val)
                is_ctrl_aaf_na = pd.isna(ctrl_aaf_val)
                if is_maf0 or is_vmiss1 or is_maf_na or is_vmiss_na or is_case_aaf_na or is_ctrl_aaf_na:
                    # Write file with flags
                    f_full.write(f"{row.VARIANT_ID}\t{is_maf0}\t{is_vmiss1}\t{is_maf_na}\t{is_vmiss_na}\t{is_case_aaf_na}\t{is_ctrl_aaf_na}\n")
                    # Write file with only variant IDs for plink2 exclusion
                    f_ids.write(f"{row.VARIANT_ID}\n")

    return full_info_path, variant_only_path


def run_plink2_exclude_variants(
    bed_prefix: str,
    output_prefix: str,
    exclude_variants_file: str,
    plink2_path: str = "/home/b/b37974/plink2",
    threads: int = 8
) -> str:
    """
    Use plink2 to remove specified variants and generate a new PLINK dataset.

    Parameters
    ----------
    bed_prefix : str
        Input PLINK binary file prefix.
    output_prefix : str
        Output PLINK binary file prefix.
    exclude_variants_file : str
        File containing variant IDs to exclude.
    plink2_path : str
        Path to the plink2 executable.
    threads : int
        Number of threads.

    Returns
    -------
    str
        Prefix of the new PLINK dataset.
    """
    cmd = [
        plink2_path,
        "--bfile", bed_prefix,
        "--exclude", exclude_variants_file,
        "--make-bed",
        "--out", output_prefix,
        "--threads", str(threads)
    ]

    subprocess.run(cmd, check=True)
    return output_prefix
