"""
Module Name: panel_compare_tools.py
===================================

Overview:
This module provides a set of tools for cross-referencing variant sites with the ToMMo reference database
and generating statistical visualization reports. It is suitable for Quality Control (QC) of Whole Genome Sequencing (WGS)
or Array data and preparation for subsequent association analyses.

Key Functions:
1. run_plink2_variant_qc_with_tommo
   - Cross-references `variant_qc_summary` with ToMMo VCF to generate an extended table with ToMMo information.
   - Outputs include IN_TOMMO, TOMMO_AAF, and TOMMO_FILTER columns.

2. plot_tommo_panel_compare_pdf
   - Generates a multi-page PDF comparison report based on `variant_qc_with_tommo.tsv`.
   - Includes statistical tables, scatter plots (PASS/Non-PASS, SNP/InDel), and histograms.

Author: ZHAO TIE
"""

import os
import sys
import subprocess
import shlex
import tempfile
import uuid
import time
import shutil
import gzip
import concurrent.futures
from datetime import datetime
from typing import Optional, Dict, Tuple, Iterable, Set

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import FuncFormatter

# --- Matplotlib Configuration ---
plt.rcParams['pdf.compression'] = 0

# Publication-oriented defaults
plt.rcParams.update({
    'font.sans-serif': ['Arial', 'DejaVu Sans', 'Liberation Sans'],
    'font.size': 10,
    'axes.titlesize': 11,
    'axes.labelsize': 10,
    'xtick.labelsize': 9,
    'ytick.labelsize': 9,
    'legend.fontsize': 9,
    'figure.titlesize': 12,
    'axes.linewidth': 0.8,
    'grid.linewidth': 0.4,
    'grid.color': '#CCCCCC',
    'savefig.bbox': 'tight',
    'savefig.pad_inches': 0.05,
})

plt.style.use('default')


def _log_info(msg: str):
    """Helper to print timestamped log messages to stderr."""
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{timestamp}] [panel_compare_tools] {msg}", file=sys.stderr, flush=True)


def _run_cmd(cmd_list: Iterable[str], stdout_path: Optional[str] = None, shell: bool = False) -> int:
    """Helper to run shell commands."""
    try:
        if stdout_path is None:
            proc = subprocess.run(cmd_list, check=False, shell=shell)
        else:
            with open(stdout_path, "wb") as fo:
                proc = subprocess.run(cmd_list, check=False, stdout=fo, shell=shell)
        return proc.returncode
    except Exception as e:
        _log_info(f"Command execution failed: {e}")
        return 1


def run_plink2_variant_qc_with_tommo(
    variant_qc_summary: str,
    tommo_vcf_path: str,
    output_path: Optional[str] = None,
    bcftools_path: str = "bcftools",
    threads: int = 8,
    chunk_size: int = 500_000,
    max_workers: Optional[int] = None,
    regions_chunk_lines: int = 50_000,
    keep_tmp: bool = False,
    tabix_path: str = "/usr/bin/tabix"
) -> str:
    """
    Function: run_plink2_variant_qc_with_tommo
    ==========================================
    Functionality:
    - For very large `variant_qc_summary` files (containing VARIANT_ID=CHROM:POS:REF:ALT),
      it first generates bcftools-readable regions files based on CHROM:POS (split by chromosome and deduplicated).
    - Runs `bcftools query -R` in parallel to extract site information from ToMMo VCF,
      using per-allele expansion to ensure exact REF/ALT matching.
    - Reads the original table in a **streaming chunked** manner and merges 3 columns:
        * IN_TOMMO: bool (whether an exact match for CHROM:POS:REF:ALT exists)
        * TOMMO_AAF: float (ToMMo's INFO/AF, corresponding to the ALT allele)
        * TOMMO_FILTER: str (FILTER for that record)
    - Finally writes out a TSV (default suffix `.variant_qc_with_tommo.tsv.gz`).

    Key Implementation Details:
    - Regions file uses two columns 1-based `CHROM\tPOS` (do NOT mix with BED coordinates).
    - Uses `bcftools query` per-allele expansion:
        Format string: `%CHROM\t%POS[\t%REF\t%ALT\t%FILTER\t%INFO/AF]\n`
      The brackets `[]` expand ALT alleles one by one, ensuring REF/ALT correspondence.
    - The merge phase does not load the entire ToMMo or summary into memory:
        * Step 1: Generate deduplicated position lists per chromosome (disk-based sort -u).
        * Step 2: Run `bcftools query` independently for each chromosome and write intermediate mapping tables.
        * Step 3: Read summary in chunks, load only the necessary subset of keys from the mapping table
          for the current chromosome, and discard after mapping.
    - Chromosome names must match exactly with VCF (e.g., `chr20` vs `20`).

    Parameters
    ----------
    variant_qc_summary : str
        Path to `*.variant_qc_summary.tsv` produced by `run_plink2_variant_qc`.
    tommo_vcf_path : str
        Path to ToMMo bgzip-compressed VCF (must have .tbi index).
    output_path : Optional[str]
        Output path; defaults to input name with suffix replaced by `.variant_qc_with_tommo.tsv.gz`.
    bcftools_path : str
        Path to `bcftools` executable (default uses `bcftools` from environment).
    threads : int
        Number of threads for bcftools (read/write/decompression).
    chunk_size : int
        Pandas chunk size for reading `variant_qc_summary`.
    max_workers : Optional[int]
        Max concurrency for running `bcftools query`; defaults to `min(4, available CPUs)`.
    regions_chunk_lines : int
        Split region lists into smaller chunks (lines) for progress tracking and controlling single call size.
    keep_tmp : bool
        Whether to keep the temporary directory for debugging.
    tabix_path : str
        Path to `tabix` executable.

    Returns
    -------
    str
        Path to the generated `variant_qc_summary_with_tommo` file.
    """
    
    # 1) Setup temporary directory & output path
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    workdir = os.path.abspath(os.path.join(os.getcwd(), f"tommo_merge_{ts}_{uuid.uuid4().hex[:8]}"))
    os.makedirs(workdir, exist_ok=True)
    
    # Determine output path
    if output_path is None:
        if variant_qc_summary.endswith(".tsv.gz"):
            base = variant_qc_summary[:-7]
        elif variant_qc_summary.endswith(".tsv"):
            base = variant_qc_summary[:-4]
        else:
            base = variant_qc_summary
        output_path = base + ".variant_qc_with_tommo.tsv.gz"
    
    _log_info(f"Temporary directory: {workdir}")

    # Progress log file
    progress_log_path = os.path.join(workdir, "progress.log")
    def _log_file(msg: str):
        _log_info(msg)
        try:
            with open(progress_log_path, 'a') as pf:
                pf.write(msg + "\n")
        except Exception:
            pass

    # 2) Step 1: Scan summary, write CHROM\tPOS region files per chromosome
    region_tmp_files: Dict[str, str] = {}
    
    # Read only #CHROM and POS columns
    reader = pd.read_csv(
        variant_qc_summary, sep='\t', usecols=["#CHROM", "POS"], dtype={"#CHROM": "string", "POS": "int64"},
        chunksize=chunk_size, engine='c'
    )
    total_rows = 0
    for chunk in reader:
        total_rows += len(chunk)
        # Write per chromosome (allow duplicates, sort -u later)
        for chrom, sub in chunk.groupby("#CHROM"):
            path = region_tmp_files.get(chrom)
            if path is None:
                path = os.path.join(workdir, f"regions.{chrom}.tsv")
                region_tmp_files[chrom] = path
            sub[["#CHROM", "POS"]].to_csv(
                path, sep='\t', header=False, index=False, mode='a'
            )
    _log_info(f"Scanned {total_rows:,} rows, generated {len(region_tmp_files)} chromosome region files (raw).")

    # 3) Deduplicate regions for each chromosome using sort -u
    uniq_region_files: Dict[str, str] = {}
    for chrom, raw_path in region_tmp_files.items():
        uniq_path = os.path.join(workdir, f"regions.{chrom}.uniq.tsv")
        # Use system sort -u, sort by (CHROM, POS) numerically
        cmd = [
            "bash", "-lc",
            f"LC_ALL=C sort -u -t$'\t' -k1,1 -k2,2n {shlex.quote(raw_path)} > {shlex.quote(uniq_path)}"
        ]
        ret = _run_cmd(cmd)
        if ret != 0:
            raise RuntimeError(f"sort -u failed for: {raw_path}")
        uniq_region_files[chrom] = uniq_path
    _log_info("Completed region deduplication for all chromosomes.")

    # 4) Run bcftools query in parallel to generate allele-level mapping tables
    mapping_files: Dict[str, str] = {}
    # Note: bcftools `[]` iterates over FORMAT/sample fields, not ALT/INFO arrays;
    # So we print comma-separated lists for ALT and INFO/AF, and expand them in Python.
    fmt = "%CHROM\t%POS\t%REF\t%ALT\t%FILTER\t%INFO/AF\n"

    def _run_bcftools_for_chrom(chrom: str) -> Tuple[str, str]:
        region_file = uniq_region_files[chrom]
        out_map = os.path.join(workdir, f"tommo.map.{chrom}.tsv")

        # Split region file into smaller chunks
        parts_dir = os.path.join(workdir, f"regions.{chrom}.parts")
        os.makedirs(parts_dir, exist_ok=True)
        part_paths = []
        with open(region_file, 'r') as fin:
            part_idx = 0
            buf = []
            for ln, line in enumerate(fin, start=1):
                buf.append(line)
                if ln % regions_chunk_lines == 0:
                    part_idx += 1
                    p = os.path.join(parts_dir, f"part_{part_idx:05d}.tsv")
                    with open(p, 'w') as fout:
                        fout.writelines(buf)
                    part_paths.append(p)
                    buf.clear()
            if buf:
                part_idx += 1
                p = os.path.join(parts_dir, f"part_{part_idx:05d}.tsv")
                with open(p, 'w') as fout:
                    fout.writelines(buf)
                part_paths.append(p)
        _log_file(f"[{chrom}] Region split into {len(part_paths)} chunks (<= {regions_chunk_lines:,} lines each)")

        # Clear/Create output file
        open(out_map, 'wb').close()

        t0 = time.time()
        for i, p in enumerate(part_paths, start=1):
            if isinstance(threads, int) and threads > 1:
                fmt_q = shlex.quote(fmt)
                cmdline = (
                    f"{shlex.quote(bcftools_path)} view --threads {threads} -R {shlex.quote(p)} "
                    f"-Ou {shlex.quote(tommo_vcf_path)} | "
                    f"{shlex.quote(bcftools_path)} query -f {fmt_q} >> {shlex.quote(out_map)}"
                )
                cmd = ["bash", "-lc", cmdline]
                ret = _run_cmd(cmd)
            else:
                cmd = [
                    bcftools_path, "query",
                    "-R", p,
                    "-f", fmt,
                    tommo_vcf_path,
                ]
                with open(out_map, 'ab') as fout:
                    proc = subprocess.run(cmd, check=False, stdout=fout)
                    ret = proc.returncode
            if ret != 0:
                raise RuntimeError(f"bcftools query failed: Chromosome {chrom} Chunk {i}/{len(part_paths)}")
            
            # Progress logging
            done = i
            total = len(part_paths)
            pct = done / total if total else 1.0
            elapsed = time.time() - t0
            rate = done / max(elapsed, 1e-6)
            _log_file(f"[{chrom}] {done}/{total} ({pct*100:.1f}%) Rate {rate:.2f} chunk/s")

        return chrom, out_map

    if max_workers is None:
        try:
            import multiprocessing as _mp
            max_workers = max(1, min(4, _mp.cpu_count()))
        except Exception:
            max_workers = 2

    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as ex:
        futs = [ex.submit(_run_bcftools_for_chrom, chrom) for chrom in uniq_region_files.keys()]
        for fut in concurrent.futures.as_completed(futs):
            chrom, out_map = fut.result()
            mapping_files[chrom] = out_map
            _log_info(f"bcftools finished: {chrom}")
    _log_info(f"Generated {len(mapping_files)} chromosome mapping tables.")

    # Diagnostic: Check for empty mapping files
    empty_maps = 0
    for _c, _path in mapping_files.items():
        try:
            _size = os.path.getsize(_path)
        except OSError:
            _size = 0
        if _size == 0:
            empty_maps += 1
    if empty_maps == len(mapping_files) and len(mapping_files) > 0:
        _log_info("[WARN] All chromosome mapping tables are empty. Please check: 1) Chromosome naming consistency (e.g., chr1 vs 1); 2) Region file format; 3) VCF index.")

    # 5) Merge Phase: Load mapping subset on demand, write output in chunks
    
    # Get header
    if variant_qc_summary.endswith(".gz"):
        with gzip.open(variant_qc_summary, 'rt') as fi:
            header_line = fi.readline().rstrip('\n')
    else:
        with open(variant_qc_summary, 'r') as fi:
            header_line = fi.readline().rstrip('\n')
            
    base_cols = header_line.split('\t')
    out_cols = base_cols + ["IN_TOMMO", "TOMMO_AAF", "TOMMO_FILTER"]

    # Raw output path
    raw_output_path = os.path.join(workdir, "raw_variant_qc_with_tommo.tsv")
    with open(raw_output_path, 'w') as fo:
        fo.write('\t'.join(out_cols) + '\n')

    def _load_mapping_subset_for_chrom(chrom: str, needed_keys: Set[str]) -> Dict[str, Tuple[str, str]]:
        """Load only needed keys for the chromosome. Returns {VID: (FILTER, AF)}."""
        out: Dict[str, Tuple[str, str]] = {}
        map_path = mapping_files.get(chrom)
        if (map_path is None) or (not os.path.exists(map_path)):
            return out
        with open(map_path, 'r') as f:
            for line in f:
                line = line.rstrip('\n')
                if not line:
                    continue
                if "\\t" in line and "\t" not in line:
                    line = line.replace("\\t", "\t")
                cols = line.split('\t')
                if len(cols) < 6:
                    continue
                c, p, r, alts_s, flt, afs_s = cols[0], cols[1], cols[2], cols[3], cols[4], cols[5]
                
                alts = alts_s.split(",") if alts_s != "." else []
                afs = afs_s.split(",") if afs_s not in (".", "") else []
                
                n = min(len(alts), len(afs)) if afs else len(alts)
                for j in range(n):
                    a = alts[j]
                    af = afs[j] if j < len(afs) else "."
                    vid = f"{c}:{p}:{r}:{a}"
                    if (not needed_keys) or (vid in needed_keys):
                        out[vid] = (flt, af)
        return out

    # Read and merge
    reader2 = pd.read_csv(
        variant_qc_summary, sep='\t', dtype="string", chunksize=chunk_size, engine='c'
    )

    processed = 0
    for chunk in reader2:
        processed += len(chunk)
        chunk["IN_TOMMO"] = False
        chunk["TOMMO_AAF"] = pd.Series([pd.NA] * len(chunk), dtype="string")
        chunk["TOMMO_FILTER"] = pd.Series([pd.NA] * len(chunk), dtype="string")

        chunk["__CHROM__"] = chunk["#CHROM"].astype("string")

        for chrom, idx in chunk.groupby("__CHROM__").groups.items():
            sub = chunk.loc[idx]
            need_keys = set(sub["ID"].tolist())
            mapping = _load_mapping_subset_for_chrom(chrom, need_keys)
            if not mapping:
                continue
            
            hit_mask = sub["ID"].isin(mapping.keys())
            if not hit_mask.any():
                continue
            vids_hit = sub.loc[hit_mask, "ID"]
            
            to_filter = {k: v[0] for k, v in mapping.items()}
            to_af = {k: v[1] for k, v in mapping.items()}

            chunk.loc[idx[hit_mask], "IN_TOMMO"] = True
            chunk.loc[idx[hit_mask], "TOMMO_FILTER"] = vids_hit.map(to_filter).astype("string").values
            
            af_str = vids_hit.map(to_af).astype("string")
            af_num = pd.to_numeric(af_str, errors='coerce')
            chunk.loc[idx[hit_mask], "TOMMO_AAF"] = af_num.astype("Float32").astype("string")

        chunk = chunk.drop(columns=["__CHROM__"])

        chunk.to_csv(
            raw_output_path, sep='\t', header=False, index=False, mode='a', na_rep='nan'
        )
        _log_info(f"Processed {processed:,} rows")

    # 6) bgzip and tabix
    bgzip_path = shutil.which("bgzip")
    if not bgzip_path and tabix_path:
        potential_bgzip = os.path.join(os.path.dirname(tabix_path), "bgzip")
        if os.path.exists(potential_bgzip):
            bgzip_path = potential_bgzip
    
    if not bgzip_path:
        _log_info("[WARNING] bgzip not found, using gzip instead (no indexing possible)")
        with open(raw_output_path, 'rb') as f_in:
            with gzip.open(output_path, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
    else:
        with open(raw_output_path, "rb") as f_in, open(output_path, "wb") as f_out:
             subprocess.run([bgzip_path, "-c"], stdin=f_in, stdout=f_out, check=True)
        
        if tabix_path and os.path.exists(tabix_path):
            cmd = [tabix_path, "-s", "1", "-b", "2", "-e", "2", "-f", output_path]
            if _run_cmd(cmd) != 0:
                _log_info(f"[WARNING] Tabix indexing failed for {output_path}")
        else:
            _log_info(f"[WARNING] Tabix path invalid: {tabix_path}")

    # 7) Cleanup
    if keep_tmp:
        _log_info(f"Keeping temporary directory: {workdir}")
    else:
        try:
            shutil.rmtree(workdir)
            _log_info(f"Cleaned up temporary directory: {workdir}")
        except Exception as e:
            _log_info(f"[WARN] Failed to clean up temporary directory: {e}")

    _log_info(f"Done. Output: {output_path}")
    return output_path


def plot_tommo_panel_compare_pdf(
    variant_qc_with_tommo: str,
    output_pdf: Optional[str] = None,
    split_metric: str = 'MAF_ALL',
    split_threshold: Optional[float] = 0.05,
    max_points: Optional[int] = None,
    png_dpi: int = 600,
    page_figsize: tuple = (14, 6),
    base_fontsize: int = 10,
    theme: str = 'academic',
    snp_color: Optional[str] = None,
    indel_color: Optional[str] = None,
):
    """
    Function: plot_tommo_panel_compare_pdf
    ======================================
    Functionality:
    - Takes the output table from `run_plink2_variant_qc_with_tommo` (*.variant_qc_with_tommo.tsv)
      and generates a multi-page PDF report.
    - Supports two modes:
      1. All variants (split_threshold=None).
      2. Split by threshold (default): e.g., MAF_ALL < 0.05 and MAF_ALL >= 0.05.

    Page Content:
    1. **Page 1**: Statistics Table (Academic Style).
    2. **Page 2**: Scatter Plot (AAF (Study) vs AAF (ToMMo)), colored by PASS/Non-PASS.
    3. **Page 3**: Scatter Plot (AAF (Study) vs AAF (ToMMo)), colored by SNP/InDel (All).
    4. **Page 4**: Scatter Plot (AAF (Study) vs AAF (ToMMo)), colored by SNP/InDel (PASS only).
    5. **Page 5**: Histogram (DIFF), PASS variants only.
    6. **Page 6**: Histogram (DIFF), PASS variants only, with IQR-based outliers removed.

    Parameters
    ----------
    variant_qc_with_tommo : str
        Path to input table.
    output_pdf : Optional[str]
        Path to output PDF.
    split_metric : str
        Metric column used for grouping (default 'MAF_ALL').
    split_threshold : Optional[float]
        Grouping threshold (default 0.05). If None, no grouping (All Variants).
    max_points : Optional[int]
        Max points to sample per group for scatter plots.
    png_dpi : int
        DPI for scatter plot PNG rendering.
    page_figsize : tuple
        Page size (width, height).
    base_fontsize : int
        Base font size.
    theme : str
        Color theme ('academic', 'okabe_ito', etc.).
    snp_color : Optional[str]
        Custom color for SNPs.
    indel_color : Optional[str]
        Custom color for InDels.

    Returns
    -------
    str
        Path to the output PDF file.
    """
    
    # --- Styling knobs (publication-ready) ---
    plt.rcParams['font.size'] = base_fontsize
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']
    plt.rcParams['axes.linewidth'] = 1.0
    plt.rcParams['xtick.major.width'] = 1.0
    plt.rcParams['ytick.major.width'] = 1.0
    
    # Professional Academic Colors
    if theme == 'academic':
        # PASS: Deep Teal, Non-PASS: Muted Coral
        pass_color = '#006d77'     
        nonpass_color = '#e29578'  
        refline_color = '#d62828'  # Strong Red
        
        # SNP: Slate Blue (Not Black), InDel: Deep Red/Magenta
        if snp_color is None: snp_color = '#457b9d' # Slate Blue
        if indel_color is None: indel_color = '#d00000' # Deep Red
        
        hist_color = '#8ecae6'
        hist_edge = '#219ebc'
    elif theme == 'okabe_ito':
        pass_color = '#0072B2'     # blue
        nonpass_color = '#ff6f01'  # vermillion
        refline_color = '#CC0000'
        if snp_color is None: snp_color = '#56B4E9' # Sky Blue
        if indel_color is None: indel_color = '#CC79A7' # Reddish Purple
        hist_color = '#56B4E9'
        hist_edge = 'white'
    else:
        pass_color = '#1f77b4'
        nonpass_color = '#ff7f0e'
        refline_color = 'red'
        if snp_color is None: snp_color = '#2ca02c'
        if indel_color is None: indel_color = '#d62728'
        hist_color = 'gray'
        hist_edge = 'white'

    ms_pass = 12
    ms_nonpass = 14
    alpha_pass = 0.6
    alpha_nonpass = 0.7
    refline_ls = (0, (5, 5))
    refline_lw = 1.2

    # Determine output path
    if output_pdf is None:
        base = os.path.splitext(variant_qc_with_tommo)[0]
        if base.endswith('.variant_qc_with_tommo'):
            base = base[:-22]
        output_pdf = base + ".panel_compare.pdf"
    
    report_path = output_pdf.replace('.pdf', '.report.txt')
    log_lines = []
    def _log(msg):
        log_lines.append(msg)
        print(msg)

    _log("========================================================")
    _log("          PANEL COMPARISON REPORT                       ")
    _log("========================================================")
    _log(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    _log(f"Input File: {variant_qc_with_tommo}")
    _log(f"Output PDF: {output_pdf}")
    _log(f"Split Metric: {split_metric}")
    _log(f"Split Threshold: {split_threshold}")
    _log("--------------------------------------------------------")

    # Read data in chunks
    needed_cols = [
        'ID', 'AAF_CTRL', 'TOMMO_AAF', 'IN_TOMMO', 'TOMMO_FILTER', split_metric
    ]
    needed_cols = list(set(needed_cols))
    
    chunk_size = 500_000
    dfs = []
    total_rows_read = 0
    
    try:
        header = pd.read_csv(variant_qc_with_tommo, sep='\t', nrows=0)
        missing = [c for c in needed_cols if c not in header.columns]
        if missing:
            raise ValueError(f"Input file missing columns: {missing}")
    except Exception as e:
        raise ValueError(f"Failed to read header: {e}")

    for chunk in pd.read_csv(variant_qc_with_tommo, sep='\t', usecols=needed_cols,
                             dtype={'ID': 'string',
                                    'IN_TOMMO': 'object',
                                    'TOMMO_FILTER': 'string'},
                             chunksize=chunk_size, engine='c'):
        if 'IN_TOMMO' in chunk.columns:
            chunk['IN_TOMMO'] = chunk['IN_TOMMO'].map(
                lambda x: True if x in (True, 1, '1', 'True', 'TRUE') else (False if x in (False, 0, '0', 'False', 'FALSE') else pd.NA)
            )
        dfs.append(chunk)
        total_rows_read += len(chunk)
    
    df = pd.concat(dfs, ignore_index=True, copy=False)
    _log(f"Total rows loaded: {total_rows_read:,}")

    # Type conversion
    for col in ['AAF_CTRL', 'TOMMO_AAF', split_metric]:
        df[col] = pd.to_numeric(df[col], errors='coerce')
    df['IN_TOMMO'] = df['IN_TOMMO'].astype('boolean')
    df['TOMMO_FILTER'] = df['TOMMO_FILTER'].astype(str).str.upper()

    # Grouping logic
    groups = []
    if split_threshold is not None:
        # Mode 2: Split by threshold
        mask_lt = df[split_metric] < split_threshold
        mask_ge = df[split_metric] >= split_threshold
        g1 = df[mask_lt].copy()
        g2 = df[mask_ge].copy()
        groups.append((f"{split_metric} < {split_threshold}", g1))
        groups.append((f"{split_metric} >= {split_threshold}", g2))
        _log(f"Group 1 ({split_metric} < {split_threshold}): {len(g1):,} variants")
        _log(f"Group 2 ({split_metric} >= {split_threshold}): {len(g2):,} variants")
    else:
        # Mode 1: All variants
        groups.append(("All Variants", df))
        _log(f"Group All: {len(df):,} variants")

    # Helper functions
    def _subset_for_scatter(g):
        g2 = g.copy()
        g2 = g2[['AAF_CTRL', 'TOMMO_AAF', 'TOMMO_FILTER']].dropna()
        if max_points is not None and len(g2) > max_points:
            g2 = g2.sample(n=max_points, random_state=42)
        return g2

    def _subset_for_hist(g):
        g2 = g.copy()
        g2 = g2[['AAF_CTRL', 'TOMMO_AAF']].dropna()
        if max_points is not None and len(g2) > max_points:
            g2 = g2.sample(n=max_points, random_state=42)
        g2['DIFF'] = g2['AAF_CTRL'] - g2['TOMMO_AAF']
        return g2

    def _classify_variant_type_from_vid(series_vid):
        atcg = {"A", "T", "C", "G"}
        def _one(vid):
            if not isinstance(vid, str): return 'InDel'
            parts = vid.split(':', 3)
            if len(parts) != 4: return 'InDel'
            ref, alt = parts[2], parts[3]
            if ref in atcg and alt in atcg and len(ref) == 1 and len(alt) == 1:
                return 'SNP'
            return 'InDel'
        return series_vid.apply(_one)

    # Start plotting
    with PdfPages(output_pdf) as pdf, tempfile.TemporaryDirectory(prefix="panel_png_") as pngdir:
        
        # --- Page 1: Statistics Table (Academic Style) ---
        n_groups = len(groups)
        fig, axes = plt.subplots(1, n_groups, figsize=page_figsize, squeeze=False)
        axes = axes.flatten()
        
        for i, (title, g) in enumerate(groups):
            ax = axes[i]
            total = len(g)
            in_tommo = (g['IN_TOMMO'] == True).sum()
            pass_tommo = (g['TOMMO_FILTER'] == 'PASS').sum()
            
            def _pct(n, d): return (n/d*100.0) if d>0 else 0.0

            rows = ['Total Variants', 'In ToMMo', 'Pass ToMMo']
            counts = [total, in_tommo, pass_tommo]
            pcts = [100.0, _pct(in_tommo, total), _pct(pass_tommo, total)]
            
            cell_text = []
            for c, p in zip(counts, pcts):
                cell_text.append([f"{c:,}", f"{p:.1f}%"])
            
            ax.axis('off')
            ax.set_title(title, fontsize=14, weight='bold', pad=20)
            
            table = ax.table(cellText=cell_text,
                             rowLabels=rows,
                             colLabels=['Count', 'Percent'],
                             loc='center',
                             cellLoc='center',
                             bbox=[0.15, 0.3, 0.7, 0.5])
            
            table.auto_set_font_size(False)
            table.set_fontsize(11)
            
            for (row, col), cell in table.get_celld().items():
                cell.set_edgecolor('black')
                cell.set_linewidth(0)
                
                if row == 0:
                    cell.set_text_props(weight='bold')
                    cell.set_linewidth(1.5)
                    cell.set_edgecolor('black')
                    cell.visible_edges = "TB"
                elif row == len(rows):
                    cell.set_linewidth(1.5)
                    cell.set_edgecolor('black')
                    cell.visible_edges = "B"
                else:
                    cell.set_linewidth(0.5)
                    cell.set_edgecolor('#d3d3d3')
                    cell.visible_edges = "B"
                
                if col == -1:
                    cell.set_text_props(weight='bold')
                    cell.set_facecolor('white')
            
            _log(f"\n--- Stats for {title} ---")
            for r, c, p in zip(rows, counts, pcts):
                _log(f"{r}: {c:,} ({p:.1f}%)")

        fig.text(0.05, 0.05, f"Split Metric: {split_metric}\nThreshold: {split_threshold}", fontsize=9, color='gray')
        pdf.savefig(fig, bbox_inches='tight')
        plt.close(fig)

        # --- Helper: Compose PNGs to PDF ---
        def _compose_pngs(png_paths, titles, suffix=''):
            fig, axes = plt.subplots(1, n_groups, figsize=page_figsize, squeeze=False, gridspec_kw={'wspace': 0.05})
            axes = axes.flatten()
            
            for i, (ax, path, t) in enumerate(zip(axes, png_paths, titles)):
                full_t = t + suffix
                if path is None:
                    ax.text(0.5, 0.5, 'No Data', ha='center', va='center')
                    ax.axis('off')
                else:
                    img = plt.imread(path)
                    ax.imshow(img, interpolation='none', aspect='auto')
                    ax.axis('off')
                
                bbox = ax.get_position()
                cx = (bbox.x0 + bbox.x1) / 2.0
                fig.text(cx, bbox.y1 + 0.01, full_t, ha='center', va='bottom', fontsize=12, weight='bold')

            pdf.savefig(fig, bbox_inches='tight')
            plt.close(fig)

        # --- Page 2: Scatter (PASS/Non-PASS) ---
        pngs = []
        ts = []
        for title, g in groups:
            g2 = _subset_for_scatter(g)
            if g2.empty:
                pngs.append(None); ts.append(title); continue
            
            g2 = g2.sort_values('AAF_CTRL', ascending=True)
            
            n_total = len(g2)
            pearson_r = g2['AAF_CTRL'].corr(g2['TOMMO_AAF'])
            _log(f"Scatter Stats (PASS/Non-PASS) - {title}: N={n_total:,}, Pearson r={pearson_r:.4f}")

            f, ax = plt.subplots(figsize=(6, 6))
            is_pass = (g2['TOMMO_FILTER'] == 'PASS')
            
            n_pass = is_pass.sum()
            n_nonpass = (~is_pass).sum()
            
            if n_pass >= n_nonpass:
                ax.scatter(g2.loc[is_pass, 'TOMMO_AAF'], g2.loc[is_pass, 'AAF_CTRL'],
                           s=ms_pass, alpha=alpha_pass, c=pass_color, label='PASS', linewidths=0, zorder=2)
                ax.scatter(g2.loc[~is_pass, 'TOMMO_AAF'], g2.loc[~is_pass, 'AAF_CTRL'],
                           s=ms_nonpass, alpha=alpha_nonpass, c=nonpass_color, label='Non-PASS', linewidths=0, zorder=3)
            else:
                ax.scatter(g2.loc[~is_pass, 'TOMMO_AAF'], g2.loc[~is_pass, 'AAF_CTRL'],
                           s=ms_nonpass, alpha=alpha_nonpass, c=nonpass_color, label='Non-PASS', linewidths=0, zorder=2)
                ax.scatter(g2.loc[is_pass, 'TOMMO_AAF'], g2.loc[is_pass, 'AAF_CTRL'],
                           s=ms_pass, alpha=alpha_pass, c=pass_color, label='PASS', linewidths=0, zorder=3)
            
            ax.plot([0, 1], [0, 1], ls=refline_ls, lw=refline_lw, c=refline_color, zorder=10)
            ax.set_xlabel('AAF (ToMMo)'); ax.set_ylabel('AAF (Study)')
            ax.set_xlim(0, 1); ax.set_ylim(0, 1)
            ax.set_aspect('equal', adjustable='box')
            ax.legend(loc='upper left', frameon=True, fontsize=9, framealpha=0.9)
            ax.grid(True, ls=':', alpha=0.4)
            
            out = os.path.join(pngdir, f"p2_{i}.png")
            f.savefig(out, dpi=png_dpi, bbox_inches='tight')
            plt.close(f)
            pngs.append(out); ts.append(title)
        
        _compose_pngs(pngs, ts)

        # --- Page 3: Scatter (SNP/InDel - All) ---
        pngs = []
        ts = []
        for title, g in groups:
            cols = ['ID', 'AAF_CTRL', 'TOMMO_AAF']
            g2 = g[cols].dropna()
            if max_points and len(g2) > max_points: g2 = g2.sample(n=max_points, random_state=42)
            
            if g2.empty:
                pngs.append(None); ts.append(title); continue
            
            g2 = g2.sort_values('AAF_CTRL', ascending=True)
            
            n_total = len(g2)
            pearson_r = g2['AAF_CTRL'].corr(g2['TOMMO_AAF'])
            _log(f"Scatter Stats (SNP/InDel All) - {title}: N={n_total:,}, Pearson r={pearson_r:.4f}")

            g2['TYPE'] = _classify_variant_type_from_vid(g2['ID'])
            is_snp = g2['TYPE'] == 'SNP'
            
            n_snp = is_snp.sum()
            n_indel = (~is_snp).sum()
            
            f, ax = plt.subplots(figsize=(6, 6))
            
            if n_snp >= n_indel:
                ax.scatter(g2.loc[is_snp, 'TOMMO_AAF'], g2.loc[is_snp, 'AAF_CTRL'],
                           s=ms_pass, alpha=alpha_pass, c=snp_color, label='SNP', linewidths=0, zorder=2)
                ax.scatter(g2.loc[~is_snp, 'TOMMO_AAF'], g2.loc[~is_snp, 'AAF_CTRL'],
                           s=ms_nonpass, alpha=alpha_nonpass, c=indel_color, label='InDel', linewidths=0, zorder=3)
            else:
                ax.scatter(g2.loc[~is_snp, 'TOMMO_AAF'], g2.loc[~is_snp, 'AAF_CTRL'],
                           s=ms_nonpass, alpha=alpha_nonpass, c=indel_color, label='InDel', linewidths=0, zorder=2)
                ax.scatter(g2.loc[is_snp, 'TOMMO_AAF'], g2.loc[is_snp, 'AAF_CTRL'],
                           s=ms_pass, alpha=alpha_pass, c=snp_color, label='SNP', linewidths=0, zorder=3)
            
            ax.plot([0, 1], [0, 1], ls=refline_ls, lw=refline_lw, c=refline_color, zorder=10)
            ax.set_xlabel('AAF (ToMMo)'); ax.set_ylabel('AAF (Study)')
            ax.set_xlim(0, 1); ax.set_ylim(0, 1)
            ax.set_aspect('equal', adjustable='box')
            ax.legend(loc='upper left', frameon=True, fontsize=9, framealpha=0.9)
            ax.grid(True, ls=':', alpha=0.4)
            
            out = os.path.join(pngdir, f"p3_{i}.png")
            f.savefig(out, dpi=png_dpi, bbox_inches='tight')
            plt.close(f)
            pngs.append(out); ts.append(title + " (All)")
            
        _compose_pngs(pngs, ts)

        # --- Page 4: Scatter (SNP/InDel - PASS Only) ---
        pngs = []
        ts = []
        for title, g in groups:
            g_pass = g[g['TOMMO_FILTER'] == 'PASS'].copy()
            cols = ['ID', 'AAF_CTRL', 'TOMMO_AAF']
            g2 = g_pass[cols].dropna()
            if max_points and len(g2) > max_points: g2 = g2.sample(n=max_points, random_state=42)
            
            if g2.empty:
                pngs.append(None); ts.append(title); continue
            
            g2 = g2.sort_values('AAF_CTRL', ascending=True)
            
            n_total = len(g2)
            pearson_r = g2['AAF_CTRL'].corr(g2['TOMMO_AAF'])
            _log(f"Scatter Stats (SNP/InDel PASS) - {title}: N={n_total:,}, Pearson r={pearson_r:.4f}")

            g2['TYPE'] = _classify_variant_type_from_vid(g2['ID'])
            is_snp = g2['TYPE'] == 'SNP'
            
            n_snp = is_snp.sum()
            n_indel = (~is_snp).sum()
            
            f, ax = plt.subplots(figsize=(6, 6))
            
            if n_snp >= n_indel:
                ax.scatter(g2.loc[is_snp, 'TOMMO_AAF'], g2.loc[is_snp, 'AAF_CTRL'],
                           s=ms_pass, alpha=alpha_pass, c=snp_color, label='SNP', linewidths=0, zorder=2)
                ax.scatter(g2.loc[~is_snp, 'TOMMO_AAF'], g2.loc[~is_snp, 'AAF_CTRL'],
                           s=ms_nonpass, alpha=alpha_nonpass, c=indel_color, label='InDel', linewidths=0, zorder=3)
            else:
                ax.scatter(g2.loc[~is_snp, 'TOMMO_AAF'], g2.loc[~is_snp, 'AAF_CTRL'],
                           s=ms_nonpass, alpha=alpha_nonpass, c=indel_color, label='InDel', linewidths=0, zorder=2)
                ax.scatter(g2.loc[is_snp, 'TOMMO_AAF'], g2.loc[is_snp, 'AAF_CTRL'],
                           s=ms_pass, alpha=alpha_pass, c=snp_color, label='SNP', linewidths=0, zorder=3)
            
            ax.plot([0, 1], [0, 1], ls=refline_ls, lw=refline_lw, c=refline_color, zorder=10)
            ax.set_xlabel('AAF (ToMMo)'); ax.set_ylabel('AAF (Study)')
            ax.set_xlim(0, 1); ax.set_ylim(0, 1)
            ax.set_aspect('equal', adjustable='box')
            ax.legend(loc='upper left', frameon=True, fontsize=9, framealpha=0.9)
            ax.grid(True, ls=':', alpha=0.4)
            
            out = os.path.join(pngdir, f"p4_{i}.png")
            f.savefig(out, dpi=png_dpi, bbox_inches='tight')
            plt.close(f)
            pngs.append(out); ts.append(title + " (PASS)")
            
        _compose_pngs(pngs, ts)

        # --- Page 5: Histogram (PASS Only) ---
        fig, axes = plt.subplots(1, n_groups, figsize=page_figsize, squeeze=False)
        axes = axes.flatten()
        
        for i, (title, g) in enumerate(groups):
            ax = axes[i]
            g_pass = g[g['TOMMO_FILTER'] == 'PASS']
            g3 = _subset_for_hist(g_pass)
            
            if g3.empty:
                ax.text(0.5, 0.5, 'No Data', ha='center', va='center')
            else:
                data = g3['DIFF']
                mean_diff = np.nanmean(data)
                std_diff = np.nanstd(data)
                
                ax.hist(data, bins=100, alpha=0.85, color=hist_color, edgecolor=hist_edge, linewidth=0.5)
                
                ax.axvline(mean_diff, c='#e63946', lw=2, ls='-', label=f'Mean: {mean_diff:.4f}')
                ax.axvline(0, c='black', ls='--', lw=1, alpha=0.5)
                
                ax.set_title(f"{title} (PASS)\n$\\mu={mean_diff:.4f}, \\sigma={std_diff:.4f}$", fontsize=11)
                ax.legend(loc='upper right', fontsize=8, frameon=True, framealpha=0.9)
                ax.yaxis.set_major_formatter(FuncFormatter(lambda x, p: f'{int(x):,}'))
            
            ax.set_xlabel('Difference (Study AAF - ToMMo AAF)')
            ax.set_ylabel('Frequency')
            ax.grid(True, ls=':', alpha=0.4)

        plt.subplots_adjust(bottom=0.15)
        pdf.savefig(fig, bbox_inches='tight')
        plt.close(fig)

        # --- Page 6: Histogram (PASS Only, Outliers Removed) ---
        fig, axes = plt.subplots(1, n_groups, figsize=page_figsize, squeeze=False)
        axes = axes.flatten()
        
        outlier_method_text = "Outlier Removal: IQR Method (Q1 - 1.5*IQR, Q3 + 1.5*IQR)"

        for i, (title, g) in enumerate(groups):
            ax = axes[i]
            g_pass = g[g['TOMMO_FILTER'] == 'PASS']
            g3 = _subset_for_hist(g_pass)
            
            if g3.empty:
                ax.text(0.5, 0.5, 'No Data', ha='center', va='center')
            else:
                data = g3['DIFF']
                # Outlier removal using IQR (Boxplot method)
                q1 = np.nanpercentile(data, 25)
                q3 = np.nanpercentile(data, 75)
                iqr = q3 - q1
                lower = q1 - 1.5 * iqr
                upper = q3 + 1.5 * iqr
                
                data_filtered = data[(data >= lower) & (data <= upper)]
                
                n_total_pass = len(data)
                n_filtered = len(data_filtered)
                n_outliers = n_total_pass - n_filtered
                pct_outliers = (n_outliers / n_total_pass * 100.0) if n_total_pass > 0 else 0.0
                
                if len(data_filtered) == 0:
                     ax.text(0.5, 0.5, 'No Data after filtering', ha='center', va='center')
                else:
                    mean_diff = np.nanmean(data_filtered)
                    std_diff = np.nanstd(data_filtered)
                    
                    ax.hist(data_filtered, bins=100, alpha=0.85, color=hist_color, edgecolor=hist_edge, linewidth=0.5)
                    
                    ax.axvline(mean_diff, c='#e63946', lw=2, ls='-', label=f'Mean: {mean_diff:.4f}')
                    ax.axvline(0, c='black', ls='--', lw=1, alpha=0.5)
                    
                    ax.set_title(f"{title} (PASS, No Outliers)\n$\\mu={mean_diff:.4f}, \\sigma={std_diff:.4f}$", fontsize=11)
                    ax.legend(loc='upper right', fontsize=8, frameon=True, framealpha=0.9)
                    ax.yaxis.set_major_formatter(FuncFormatter(lambda x, p: f'{int(x):,}'))
                    
                    stats_text = f"Outliers removed: {n_outliers:,}\n({pct_outliers:.2f}% of PASS)"
                    ax.text(0.05, 0.95, stats_text, transform=ax.transAxes, fontsize=9,
                            verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
            
            ax.set_xlabel('Difference (Study AAF - ToMMo AAF)')
            ax.set_ylabel('Frequency')
            ax.grid(True, ls=':', alpha=0.4)

        plt.subplots_adjust(bottom=0.15)
        fig.text(0.05, 0.02, outlier_method_text, fontsize=9, color='gray', ha='left')
        pdf.savefig(fig, bbox_inches='tight')
        plt.close(fig)

    # Write report
    with open(report_path, 'w') as f:
        f.write('\n'.join(log_lines))
    
    return output_pdf

