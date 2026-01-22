"""
rvtest_prepare_tools.py - RVTest Data Preprocessing Toolkit

This module provides a set of robust tool functions for preparing data files required for RVTest.
It features comprehensive logging of data statistics (sample counts, variant counts) at each processing step.

Functionalities:
1. plink_to_vcf_raw: Converts PLINK data to standardized VCF format with optional normalization and filtering.
2. reformat_pheno_covar: Reformats phenotype and covariate files, handling header cleaning and ID matching.
3. reformat_refflat_remove_chr: Processes refFlat gene annotation files to ensure chromosome naming consistency.

Author: ZHAO TIE
Last Updated: Jan 16, 2026
"""

import os
import shutil
import subprocess
import logging
import pandas as pd
import gzip
from typing import Tuple

def get_plink_stats(fam_path: str, bim_path: str) -> Tuple[int, int]:
    """Count lines in .fam and .bim files to get sample and variant counts."""
    n_samples, n_variants = 0, 0
    try:
        if os.path.exists(fam_path):
            with open(fam_path, 'r') as f:
                n_samples = sum(1 for _ in f)
        if os.path.exists(bim_path):
            with open(bim_path, 'r') as f:
                n_variants = sum(1 for _ in f)
    except Exception as e:
        logging.warning(f"Could not read PLINK stats from {fam_path}/{bim_path}: {e}")
    return n_samples, n_variants

def get_vcf_stats(vcf_path: str, bcftools_path: str = "bcftools", threads: int = 1) -> Tuple[int, int]:
    """Get sample and variant counts from a VCF file using bcftools."""
    n_samples, n_variants = 0, 0
    try:
        # Get variant count from index (fastest) -n
        cmd_idx = [bcftools_path, "index", "-n", "--threads", str(threads), vcf_path]
        res_idx = subprocess.run(cmd_idx, capture_output=True, text=True)
        if res_idx.returncode == 0:
            try:
                n_variants = int(res_idx.stdout.strip())
                logging.debug(f"Retrieved variant count from index for {vcf_path}")
            except ValueError:
                pass 

        # Get sample count (fast) - query -l
        cmd_n_samples = [bcftools_path, "query", "-l", vcf_path]
        res_ns = subprocess.run(cmd_n_samples, capture_output=True, text=True)
        if res_ns.returncode == 0:
            samples = res_ns.stdout.strip().splitlines()
            n_samples = len(samples)
            
    except Exception as e:
        logging.warning(f"Could not read VCF stats for {vcf_path}: {e}")
    return n_samples, n_variants

def plink_to_vcf_raw(
    bed_prefix: str,
    plink_path: str = "/home/b/b37974/plink2_alpha6/plink2",
    bcftools_path: str = "bcftools",
    ref_fa: str = "/LARGE0/gr10478/b37974/Pulmonary_Hypertension/nagasaki_pipeline/data/hs38DH.fa",
    threads: int = 6,
    rename_chr: str = "/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/info/rename_chr.txt",
    snps_only_just_acgt: bool = False,
    tabix_path: str = "tabix",
    keep_chr_prefix: bool = True,
    perform_norm: bool = False,
) -> str:
    """
    Converts PLINK format (.bed/.bim/.fam) to **bgzip compressed VCF (.vcf.gz)**.
    It optionally performs **bcftools norm** normalization using the reference genome
    and finally generates a `.tbi` index using **tabix**.
    
    Logs sample and variant counts at each major step.

    Changes (compared to older versions)
    ------------------------------------
    1) All output files are saved to the **current working directory** (instead of `./tmp`).
       Only the **final** `*.vcf.gz` and `*.vcf.gz.tbi` are kept; intermediate files from export/rename stages are deleted.
    2) New parameter `snps_only_just_acgt` (default False): Adds `--snps-only just-acgt` to plink export to extract only A/C/G/T SNPs (sample names remain IID via `--recode vcf-iid bgz`).
    3) Final indexing uses **tabix** (default path `"tabix"`); other stages use `bcftools index` if needed.
    4) Acceleration & Memory Overflow Prevention: All subprocesses **do not capture STDOUT** (redirected to `/dev/null`), capturing only `stderr`; adds `--threads` to `bcftools` commands where possible.
    5) **Naming Convention**: If only SNPs are exported (`snps_only_just_acgt=True`), the final filename appends `.snp`; otherwise, it appends `.all` (e.g., `{basename}.snp.norm.vcf.gz` vs `{basename}.all.norm.vcf.gz`).
    6) New parameter `keep_chr_prefix` (default True): When False, removes 'chr' prefix from the #CHROM column in the VCF file.
    7) New parameter `perform_norm` (default False): Control whether to perform bcftools normalization.
    8) **Stats Logging**: Logs input PLINK stats and output VCF stats.

    Args:
        bed_prefix (str): Prefix for PLINK input (without suffix), e.g., "/path/to/data".
        plink_path (str): Path to plink executable (default: plink2).
        bcftools_path (str): Path to bcftools executable.
        ref_fa (str): Reference genome FASTA file (must match current coordinates and should have .fai index).
        threads (int): Number of parallel threads (passed to bcftools; plink also supports --threads).
        rename_chr (str): Path to `bcftools annotate --rename-chrs` map file (optional; applied before `bcftools norm` if exists).
        snps_only_just_acgt (bool): If True, adds `--snps-only just-acgt` during plink export to keep only A/C/G/T SNPs.
        tabix_path (str): Path to tabix executable (for final `.tbi` index).
        keep_chr_prefix (bool): If True, keeps 'chr' prefix in chromosome names; if False, removes 'chr' prefix.
        perform_norm (bool): If True, performs bcftools normalization.

    Returns:
        str: Absolute path to the finalized `.vcf.gz` file.
    """
    import os
    import shutil
    import subprocess

    def _run(cmd, desc: str):
        """Run external command: do not capture STDOUT, capture and echo STDERR on failure."""
        proc = subprocess.run(
            cmd,
            stdout=subprocess.DEVNULL,  # Prevent large output from consuming memory
            stderr=subprocess.PIPE,
            text=True,
            check=False,
        )
        if proc.returncode != 0:
            raise RuntimeError(
                f"{desc} Failed\nCommand: {' '.join(cmd)}\nSTDERR:\n" + (proc.stderr or "")
            )
        return proc

    # --- 0) Initial Stats (PLINK Input) ---
    logging.info("--- Step 0: Initial PLINK Input Stats ---")
    try:
        n_fam, m_bim = get_plink_stats(f"{bed_prefix}.fam", f"{bed_prefix}.bim")
        logging.info(f"Input PLINK: {n_fam} samples, {m_bim} variants")
    except Exception as e:
        logging.warning(f"Could not get PLINK stats: {e}")

    # --- Verify Executables ---
    if not (os.path.isfile(plink_path) or shutil.which(plink_path)):
        # Try to find in path if just a name is given
        if not shutil.which(plink_path):
             # Fallback to absolute path check or raise error
             pass # Will fail later or assume user environment is correct if shutil.which fails but it's in path via shell expansion (unlikely with subprocess.run)
             # Better to trust shutil.which or os.path.isfile
             if not (plink_path == "plink" and shutil.which("plink")): # Common case
                 pass 
                 
    # simplified check: strictly rely on shutil.which if it's a command name, else isfile
    def check_executable(path):
        if os.sep in path:
            return os.path.isfile(path) and os.access(path, os.X_OK)
        return shutil.which(path) is not None

    if not check_executable(plink_path):
         logging.warning(f"PLINK executable not found at {plink_path}. Please check path or install plink2.")

    # --- Verify Input Files ---
    req = [f"{bed_prefix}.bed", f"{bed_prefix}.bim", f"{bed_prefix}.fam"]
    missing = [p for p in req if not os.path.exists(p) or os.path.getsize(p) == 0]
    if missing:
        raise FileNotFoundError("Missing PLINK input files or files are empty: " + ", ".join(missing))

    # --- Thread Parameters ---
    threads = int(threads) if isinstance(threads, (int, str)) else 6
    plink_threads = max(1, int(threads))

    # --- Output Prefix: CWD + Basename + Suffix based on SNP-only flag ---
    cwd = os.getcwd()
    base_name = os.path.basename(bed_prefix)
    suffix_tag = ".snp" if snps_only_just_acgt else ".all"
    out_prefix = os.path.join(cwd, f"{base_name}{suffix_tag}")

    # Path Definitions
    vcf_gz_path = f"{out_prefix}.vcf.gz"         # Direct export from plink
    tbi_path = f"{vcf_gz_path}.tbi"
    renamed_vcf = f"{out_prefix}.ren.vcf.gz"     # After chromosome renaming
    renamed_tbi = f"{renamed_vcf}.tbi"
    norm_vcf_gz_path = f"{out_prefix}.norm.vcf.gz"  # After normalization
    norm_tbi_path = f"{norm_vcf_gz_path}.tbi"
    nochr_vcf_gz_path = f"{out_prefix}.nochr.vcf.gz"  # After removing chr prefix (final if no norm and remove chr)
    
    # Determine the current VCF file to work on
    current_vcf = vcf_gz_path
    current_tbi = tbi_path
    
    cleanup_files = []

    # --- 1) Connect PLINK to export bgzip VCF (Sample names as IID; optional SNPs only) ---
    export_cmd = [
        plink_path,
        "--bfile", bed_prefix,
        "--export", "vcf", "id-paste=iid", "bgz",
        "--threads", str(plink_threads),
        "--out", out_prefix,
    ]
    if snps_only_just_acgt:
        export_cmd.extend(["--snps-only", "just-acgt"])  # Only A/C/G/T SNPs

    print("[Info] Exporting bgzip VCF using plink2 (Sample names as IID)...")
    _run(export_cmd, "plink export")
    if not os.path.exists(vcf_gz_path) or os.path.getsize(vcf_gz_path) == 0:
        raise RuntimeError(f"Exported .vcf.gz file not found or empty: {vcf_gz_path}")
    
    # No need to add vcf_gz_path to cleanup yet, it's the base.

    # Index the initial VCF
    index_cmd_1 = [bcftools_path, "index", "-t", "--threads", str(plink_threads), vcf_gz_path]
    print("[Info] Indexing initial .vcf.gz using bcftools...")
    _run(index_cmd_1, "bcftools index (Initial VCF)")
    if not os.path.exists(tbi_path) or os.path.getsize(tbi_path) == 0:
        raise RuntimeError(f"Initial tbi index not found or empty: {tbi_path}")

    # --- Step 1: Stats After PLINK Export ---
    logging.info("--- Step 1: Stats After PLINK Export (Genotype) ---")
    try:
        n_vcf, m_vcf = get_vcf_stats(vcf_gz_path, bcftools_path, threads=plink_threads)
        logging.info(f"After PLINK Export: {n_vcf} samples, {m_vcf} variants")
    except Exception as e:
        logging.warning(f"Could not get VCF stats after export: {e}")
    
    cleanup_files.extend([vcf_gz_path, tbi_path])

    # --- 2) Optional: Chromosome Renaming (bcftools annotate --rename-chrs) ---
    if rename_chr and isinstance(rename_chr, str) and rename_chr.strip() and os.path.exists(rename_chr):
        print(f"[Info] rename_chr map file detected, performing chromosome renaming: {rename_chr}")
        rename_cmd = [
            bcftools_path, "annotate",
            "--rename-chrs", rename_chr,
            "--threads", str(plink_threads),
            "-Oz",
            "-o", renamed_vcf,
            current_vcf,
        ]
        _run(rename_cmd, "bcftools annotate --rename-chrs")

        # Index after renaming
        index_cmd_rename = [bcftools_path, "index", "-t", "--threads", str(plink_threads), renamed_vcf]
        print("[Info] Indexing renamed VCF...")
        _run(index_cmd_rename, "bcftools index (Renamed VCF)")
        if not os.path.exists(renamed_tbi) or os.path.getsize(renamed_tbi) == 0:
            raise RuntimeError(f"Renamed tbi index not found or empty: {renamed_tbi}")
        
        current_vcf = renamed_vcf
        current_tbi = renamed_tbi
        cleanup_files.extend([renamed_vcf, renamed_tbi])
    else:
        if rename_chr and isinstance(rename_chr, str) and rename_chr.strip():
            print(f"[Warning] Specified rename_chr file does not exist, skipping chromosome renaming: {rename_chr}")

    # --- 3) Optional: bcftools norm Normalization ---
    if perform_norm:
        if not os.path.exists(ref_fa) or os.path.getsize(ref_fa) == 0:
            raise FileNotFoundError(f"Reference genome FASTA not found: {ref_fa}")

        norm_cmd = [
            bcftools_path, "norm",
            "--multiallelics", "-any",
            "--fasta-ref", ref_fa,
            "--check-ref", "s",
            "--threads", str(plink_threads),
            current_vcf,
            "-Oz",
            "-o", norm_vcf_gz_path,
        ]
        print("[Info] Performing normalization using bcftools norm (--multiallelics -any, --check-ref s)...")
        _run(norm_cmd, "bcftools norm")
        if not os.path.exists(norm_vcf_gz_path) or os.path.getsize(norm_vcf_gz_path) == 0:
            raise RuntimeError(f"Normalized .vcf.gz file not found or empty: {norm_vcf_gz_path}")
        
        # Index normalized VCF
        norm_tbi_cmd = [bcftools_path, "index", "-t", "--threads", str(plink_threads), norm_vcf_gz_path]
        _run(norm_tbi_cmd, "bcftools index (Normalized VCF)")

        current_vcf = norm_vcf_gz_path
        current_tbi = norm_tbi_path
        cleanup_files.extend([norm_vcf_gz_path, norm_tbi_path])

        # --- Step 2: Stats After Normalization ---
        logging.info("--- Step 2: Stats After Normalization ---")
        try:
            n_norm, m_norm = get_vcf_stats(norm_vcf_gz_path, bcftools_path, threads=plink_threads)
            logging.info(f"After Normalization: {n_norm} samples, {m_norm} variants")
        except Exception as e:
            logging.warning(f"Could not get VCF stats after normalization: {e}")

    else:
        logging.info("--- Step 2: Normalization Skipped ---")
    
    # --- 4) Optional: Remove 'chr' Prefix ---
    final_vcf_gz = current_vcf
    
    if not keep_chr_prefix:
        print("[Info] Removing 'chr' prefix from chromosome names...")
        # Create temporary chromosome rename map (chr -> no prefix)
        import tempfile
        prefix_removed_vcf = f"{out_prefix}.nochr.vcf.gz"
        
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.chr_remap.txt') as tmp_chr_file:
            # Write chromosome rename map: chr1 -> 1, chr2 -> 2, etc.
            for i in range(1, 23):  # 1-22
                tmp_chr_file.write(f"chr{i}\t{i}\n")
            tmp_chr_file.write("chrX\tX\n")
            tmp_chr_file.write("chrY\tY\n")
            tmp_chr_file.write("chrMT\tMT\n")
            tmp_chr_file.write("chrM\tM\n")
            tmp_chr_remap_path = tmp_chr_file.name
        
        try:
            # Use bcftools annotate --rename-chrs to remove chr prefix
            remove_chr_cmd = [
                bcftools_path, "annotate",
                "--rename-chrs", tmp_chr_remap_path,
                "--threads", str(plink_threads),
                "-Oz",
                "-o", prefix_removed_vcf,
                current_vcf,
            ]
            _run(remove_chr_cmd, "Remove chr prefix")
            if not os.path.exists(prefix_removed_vcf) or os.path.getsize(prefix_removed_vcf) == 0:
                raise RuntimeError(f"VCF file with chr prefix removed not found or empty: {prefix_removed_vcf}")
            
            final_vcf_gz = prefix_removed_vcf

            # --- Step 3: Stats After Remove Chr ---
            logging.info("--- Step 3: Stats After Removing 'chr' Prefix ---")
            try:
                # Indexing needed for stats
                tmp_index_cmd = [bcftools_path, "index", "-t", "--threads", str(plink_threads), prefix_removed_vcf]
                _run(tmp_index_cmd, "Index for stats (No Chr)")
                
                n_nochr, m_nochr = get_vcf_stats(prefix_removed_vcf, bcftools_path, threads=plink_threads)
                logging.info(f"After Remove Chr: {n_nochr} samples, {m_nochr} variants")
                
                # We created an index, so we should clean it up if we don't want it, 
                # but usually having an index is good. The caller expects a VCF.
                # The tbi file will be prefix_removed_vcf + ".tbi"
                if os.path.exists(prefix_removed_vcf + ".tbi"):
                    cleanup_files.append(prefix_removed_vcf + ".tbi")
            except Exception as e:
                logging.warning(f"Could not get VCF stats after remove chr: {e}")

            cleanup_files.append(prefix_removed_vcf) # Will be removed from cleanup list at the end
            
        finally:
            # Clean up temporary file
            try:
                os.remove(tmp_chr_remap_path)
            except Exception:
                pass
    
    # --- 5) Final Indexing using tabix ---
    final_tbi = f"{final_vcf_gz}.tbi"
    print(f"[Info] Creating tbi index for final VCF using tabix: {final_vcf_gz}")
    tabix_cmd = [tabix_path, "-f", "-p", "vcf", final_vcf_gz]
    _run(tabix_cmd, "tabix index (Final VCF)")
    
    if not os.path.exists(final_tbi) or os.path.getsize(final_tbi) == 0:
        raise RuntimeError(f"Final tbi index not found or empty: {final_tbi}")

    # --- 6) Clean up intermediate files (Keep only final files) ---
    # Convert paths to absolute for safe comparison
    abs_final_vcf = os.path.abspath(final_vcf_gz)
    abs_final_tbi = os.path.abspath(final_tbi)
    
    for f in cleanup_files:
        try:
            abs_f = os.path.abspath(f)
            if abs_f != abs_final_vcf and abs_f != abs_final_tbi:
                if os.path.exists(f):
                    os.remove(f)
        except Exception as e:
            print(f"[Warning] Failed to delete intermediate file: {f}, Error: {e}")

    print(f"[Complete] Final compressed VCF generated: {abs_final_vcf}")
    print(f"[Complete] Final index generated: {abs_final_tbi}")
    return abs_final_vcf


def reformat_pheno_covar(pheno_path: str, covar_path: str, out_dir: str = None, out_sep: str = " ", force: bool = True) -> tuple[str, str]: # type: ignore
    """
    Reformats phenotype and covariate files for RVTest:
    1. Converts all column names to lowercase.
    2. Removes '#' characters from column names.
    3. Outputs new pheno and covar files to specified directory (default: CWD).
    4. Configurable output separator (default: space " ").
    5. If force=True, inserts fatid, matid columns (filled with 0) after iid in phenotype file, then inserts sex column (from covariate file).
    
    Args:
        pheno_path (str): Path to raw phenotype file (tab-separated).
        covar_path (str): Path to raw covariate file (tab-separated).
        out_dir (str, optional): Output directory. Defaults to CWD.
        out_sep (str, optional): Output separator. Defaults to " ".
        force (bool, default True): Whether to force insertion of fatid, matid, sex columns into phenotype file.
    
    Returns:
        tuple[str, str]: (New phenotype file path, New covariate file path)
    """
    import os
    import pandas as pd

    if out_dir is None:
        out_dir = os.getcwd()

    def _process_file(path: str, suffix: str) -> str:
        df = pd.read_csv(path, sep="\t", dtype=str)
        logging.info(f"Loaded {os.path.basename(path)}: {len(df)} rows")
        # Column name processing
        new_cols = [c.lower().replace("#", "") for c in df.columns]
        df.columns = new_cols
        out_path = os.path.join(out_dir, os.path.basename(path).replace(".txt", f".{suffix}.csv"))
        df.to_csv(out_path, sep=out_sep, index=False)
        logging.info(f"Written {out_path}: {len(df)} rows")
        return out_path

    # Process Covariate File
    new_covar = _process_file(covar_path, "covar.reformat")
    
    # Process Phenotype File
    if force:
        # Read and process phenotype file
        pheno_df = pd.read_csv(pheno_path, sep="\t", dtype=str)
        logging.info(f"Loaded Phenotype (Force): {len(pheno_df)} samples")
        new_cols = [c.lower().replace("#", "") for c in pheno_df.columns]
        pheno_df.columns = new_cols
        
        # Read covariate file to get sex info
        covar_df = pd.read_csv(covar_path, sep="\t", dtype=str)
        logging.info(f"Loaded Covariate: {len(covar_df)} samples")
        covar_cols = [c.lower().replace("#", "") for c in covar_df.columns]
        covar_df.columns = covar_cols
        
        # Find iid column position
        if 'iid' not in pheno_df.columns:
            raise ValueError("Column 'iid' not found in phenotype file")
        
        iid_pos = list(pheno_df.columns).index('iid')
        
        # Insert fatid and matid columns after iid (filled with "0")
        pheno_df.insert(iid_pos + 1, 'fatid', '0')
        pheno_df.insert(iid_pos + 2, 'matid', '0')
        
        # Get sex info from covariate file and insert
        if 'sex' in covar_df.columns and 'fid' in covar_df.columns and 'iid' in covar_df.columns:
            # Create dictionary for matching
            sex_dict = {}
            for _, row in covar_df.iterrows():
                key = (str(row['fid']), str(row['iid']))
                sex_dict[key] = str(row['sex'])
            
            # Add sex column to phenotype DF
            sex_values = []
            for _, row in pheno_df.iterrows():
                key = (str(row['fid']), str(row['iid']))
                sex_values.append(sex_dict.get(key, '0'))  # Default to '0' if not found
            
            pheno_df.insert(iid_pos + 3, 'sex', sex_values)
        else:
            print("[Warning] Necessary columns (fid, iid, sex) missing in covariate file. 'sex' column will be filled with '0'")
            pheno_df.insert(iid_pos + 3, 'sex', '0')
        
        # Save processed phenotype file
        pheno_out_path = os.path.join(out_dir, os.path.basename(pheno_path).replace(".txt", ".pheno.reformat.csv"))
        pheno_df.to_csv(pheno_out_path, sep=out_sep, index=False)
        new_pheno = pheno_out_path
    else:
        # Use original processing method if not forced
        new_pheno = _process_file(pheno_path, "pheno.reformat")

    print(f"[Complete] New phenotype file output: {new_pheno}")
    print(f"[Complete] New covariate file output: {new_covar}")

    return new_pheno, new_covar


def reformat_refflat_remove_chr(refflat_path: str, out_dir: str = None) -> str:  # type: ignore
    """
    Reformats refFlat file, removing 'chr' prefix from the 3rd column (chromosome column).
    Optimized version: Uses stream processing (reading/writing line by line) to significantly reduce memory usage and improve speed.
    
    refFlat file is typically tab-separated:
    geneName, name, chrom, strand, txStart, txEnd, cdsStart, cdsEnd, exonCount, exonStarts, exonEnds
    
    Args:
        refflat_path (str): Input refFlat file path (.txt.gz format).
        out_dir (str, optional): Output directory. Defaults to CWD.
        
    Returns:
        str: absolute path to the processed refFlat file (.txt.gz format).
    """
    import os
    import gzip
    
    if out_dir is None:
        out_dir = os.getcwd()
    
    # Check if input file exists
    if not os.path.exists(refflat_path):
        raise FileNotFoundError(f"refFlat file not found: {refflat_path}")
    
    # Generate output filename
    base_name = os.path.basename(refflat_path)
    if base_name.endswith('.txt.gz'):
        out_name = base_name.replace('.txt.gz', '.nochr.txt.gz')
    elif base_name.endswith('.gz'):
        out_name = base_name.replace('.gz', '.nochr.gz')
    else:
        out_name = base_name + '.nochr.gz'
    
    out_path = os.path.join(out_dir, out_name)
    
    print(f"[Info] Processing refFlat file: {refflat_path}")
    print(f"[Info] Output file: {out_path}")
    
    try:
        # Statistics
        line_count = 0
        chroms_before = set()
        chroms_after = set()
        
        # Stream processing
        with gzip.open(refflat_path, 'rt', encoding='utf-8') as infile, \
             gzip.open(out_path, 'wt', encoding='utf-8', compresslevel=6) as outfile:
            
            for line in infile:
                line = line.rstrip('\n\r')
                if not line:  # Skip empty lines
                    continue
                
                fields = line.split('\t')
                
                # Check column count (on first line)
                if line_count == 0 and len(fields) < 3:
                    raise ValueError(f"refFlat file has insufficient columns. Need at least 3, found {len(fields)}")
                
                # Process 3rd column (index 2) - chromosome info
                if len(fields) > 2:
                    original_chrom = fields[2]
                    chroms_before.add(original_chrom)
                    
                    # Remove 'chr' prefix
                    if original_chrom.startswith('chr'):
                        fields[2] = original_chrom[3:]  # Slice off 'chr'
                    
                    chroms_after.add(fields[2])
                
                # Write processed line
                outfile.write('\t'.join(fields) + '\n')
                line_count += 1
                
                # Progress update every 1M lines
                if line_count % 1000000 == 0:
                    print(f"[Progress] Processed {line_count:,} lines")
        
        print(f"[Info] Total processed: {line_count:,} lines")
        
        # Display chromosome stats (limited)
        chroms_before_sorted = sorted(list(chroms_before))
        chroms_after_sorted = sorted(list(chroms_after))
        
        print(f"[Info] Chromosomes before ({len(chroms_before_sorted)} types): {chroms_before_sorted[:10]}{'...' if len(chroms_before_sorted) > 10 else ''}")
        print(f"[Info] Chromosomes after ({len(chroms_after_sorted)} types): {chroms_after_sorted[:10]}{'...' if len(chroms_after_sorted) > 10 else ''}")
        
        print(f"[Complete] Generated refFlat file without chr prefix: {out_path}")
        return os.path.abspath(out_path)
        
    except Exception as e:
        # Cleanup on failure
        try:
            if os.path.exists(out_path):
                os.remove(out_path)
        except Exception:
            pass
        raise RuntimeError(f"Error processing refFlat file: {str(e)}")
