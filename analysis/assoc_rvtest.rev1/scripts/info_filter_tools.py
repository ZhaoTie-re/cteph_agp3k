#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
VCF INFO Field Filtering Tool
=============================

This module provides efficient INFO field filtering for VCF files based on bcftools,
optimized for large genomic datasets.

Main Features
-------------
1. **INFO Field Filtering**: Filter VCF variants based on specified INFO field values.
2. **Chromosome Prefix Normalization**: Automatically detect and unify 'chr' prefix formats.
3. **Multiple Match Modes**: Supports exact match, substring match, and regex match.
4. **Parallel Processing**: Utilizes multi-threading to accelerate compression and processing.
5. **Automatic Indexing**: Automatically creates tabix index for output files.

Use Cases
---------
- Filter high or moderate impact variants (e.g., from SnpEff annotations).
- Select specific variant types based on functional annotations.
- Batch process large WGS/WES datasets.
- Prepare subset data for association analysis.

Performance Features
--------------------
- Processes bgzip-compressed VCF files directly without decompression.
- Uses native bcftools filtering for low memory footprint.
- Supports multi-threaded parallel compression.
- Atomic operations ensure file integrity.

Dependencies
------------
- bcftools: VCF file manipulation and filtering.
- tabix: VCF file indexing.
- bgzip: File compression.

Typical Usage Example
---------------------
```python
from info_filter_tools import filter_vcf_by_info

# Filter variants with HIGH and MODERATE impact
output_path = filter_vcf_by_info(
    vcf_path="input.vcf.gz",
    info_key="impact",
    values=["HIGH", "MODERATE"],
    out_prefix="filtered_variants"
)

# Use regex matching
output_path = filter_vcf_by_info(
    vcf_path="input.vcf.gz",
    info_key="consequence",
    values="missense.*",
    match_mode="regex"
)
```

Author: ZHAO TIE
Date: 2025-10-14
"""

import os
import shlex
import subprocess
import re
import difflib
import shutil
from typing import Iterable, List, Optional, Union, Dict, Any
import logging
from datetime import datetime
import tempfile


def get_vcf_stats(vcf_path: str, bcftools_path: str = "bcftools") -> Dict[str, Any]:
    """
    Get statistics for a VCF file: number of samples and number of records.
    
    Args:
        vcf_path: Path to the VCF file.
        bcftools_path: Path to the bcftools executable.
        
    Returns:
        dict: A dictionary containing 'num_samples' and 'num_records'.
    """
    stats = {
        'num_samples': 0,
        'num_records': 0
    }
    
    if not os.path.exists(vcf_path):
        return stats
        
    # Get number of samples
    try:
        # bcftools query -l lists samples
        cmd_samples = [bcftools_path, "query", "-l", vcf_path]
        result = subprocess.run(cmd_samples, capture_output=True, text=True, check=True)
        samples = result.stdout.strip().splitlines()
        stats['num_samples'] = len(samples)
    except Exception as e:
        logging.getLogger("info_filter_tools").warning(f"Failed to count samples: {e}")
        
    # Get number of records
    # Try using index first (fastest)
    try:
        # Check for .tbi or .csi
        if os.path.exists(vcf_path + ".tbi") or os.path.exists(vcf_path + ".csi"):
            cmd_index = [bcftools_path, "index", "--nrecords", vcf_path]
            result = subprocess.run(cmd_index, capture_output=True, text=True, check=True)
            try:
                stats['num_records'] = int(result.stdout.strip())
            except ValueError:
                 # Fallback if output is not just a number
                 stats['num_records'] = 0
        else:
            # Fallback to counting lines (slower)
            # Use -H to suppress header
            logging.getLogger("info_filter_tools").info("No index found, counting records by scanning file (this may be slow)...")
            cmd_count = [bcftools_path, "view", "-H", vcf_path]
            # Pipe to wc -l
            p1 = subprocess.Popen(cmd_count, stdout=subprocess.PIPE)
            p2 = subprocess.Popen(["wc", "-l"], stdin=p1.stdout, stdout=subprocess.PIPE, text=True)
            if p1.stdout:
                p1.stdout.close()
            output, _ = p2.communicate()
            try:
                stats['num_records'] = int(output.strip())
            except ValueError:
                stats['num_records'] = 0
    except Exception as e:
        logging.getLogger("info_filter_tools").warning(f"Failed to count records: {e}")
        
    return stats


# Helper: collect INFO IDs from VCF header using bcftools
def _collect_info_ids_from_header(vcf_path: str, bcftools_path: str) -> set:
    """Use bcftools to read VCF header and extract IDs of all INFO fields.
    Raises exception if bcftools call fails.
    """
    proc = subprocess.run(
        [bcftools_path, "view", "-h", vcf_path],
        check=True,
        capture_output=True,
        text=True,
    )
    info_ids = set()
    for line in proc.stdout.splitlines():
        # Format: ##INFO=<ID=impact,Number=.,Type=String,Description="...">
        if line.startswith("##INFO=<ID="):
            # Extract part between ID= and subsequent comma
            try:
                id_part = line.split("##INFO=<ID=", 1)[1]
                field_id = id_part.split(",", 1)[0].strip()
                if field_id:
                    info_ids.add(field_id)
            except Exception:
                pass
    return info_ids


def _detect_chr_prefix(vcf_path: str, bcftools_path: str, sample_lines: int = 200) -> Optional[bool]:
    """Detect if VCF CHROM column has 'chr' prefix.
    Returns: True if present; False if absent; None if undetermined (e.g., no variants).
    Note: Samples the first few variant lines to avoid reading the whole file.
    """
    try:
        # Output body (no header) first few lines
        # -H only output variant lines; no region limit, just sample first sample_lines
        cmd = [bcftools_path, "view", "-H", vcf_path]
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True)
        has_chr = None
        cnt = 0
        if proc.stdout is None:
            return None
        for line in proc.stdout:
            if not line.strip():
                continue
            chrom = line.split("\t", 1)[0]
            if chrom.startswith("#"):
                # Safety check: skip unexpected header lines
                continue
            has_chr = chrom.startswith("chr")
            cnt += 1
            if cnt >= sample_lines:
                break
        proc.stdout.close()
        proc.wait()
        return has_chr
    except Exception:
        return None


def filter_vcf_by_info(
    vcf_path: str = "cteph_agp3k.rare.all.nochr.norm.chrprefix.snpeff.vcf.gz",
    *,
    bcftools_path: str = "/home/b/b37974/bcftools/bcftools",
    tabix_path: str = "/home/b/b37974/htslib-1.9/tabix",
    info_key: str = "INFO/impact",
    values: Union[str, Iterable[str]] = ("HIGH", "MODERATE"),
    logic: str = "any",
    match_mode: str = "exact",
    check_chr_prefix: bool = False,
    keep_chr_prefix: bool = False,
    threads: int = 8,
    out_prefix: Optional[str] = None,
    out_dir: Optional[str] = os.getcwd(),
    index_output: bool = True,
    index_type: str = "tbi",
    overwrite: bool = True,
    dry_run: bool = False,
) -> str:
    """Filter large VCF.GZ by INFO field (using bcftools) and return the generated VCF path.

    Args
    ----
    vcf_path : str
        Input .vcf.gz file path (must be bgzip compressed).
    bcftools_path : str
        Path to bcftools executable.
    tabix_path : str
        Path to tabix executable (used for output indexing).
    info_key : str
        INFO key to filter by. Can be "impact" or "INFO/impact" (auto-completed).
    values : str | Iterable[str]
        Values to keep (single or multiple). E.g., "HIGH" or ["HIGH", "MODERATE"].
    logic : {"any", "all"}
        Logic for multiple values. "any" means OR, "all" means AND.
    match_mode : {"exact", "contains", "regex"}
        Matching mode. exact: string equality; contains: regex substring; regex: treated as regex pattern.
    check_chr_prefix : bool
        Whether to check and unify CHROM 'chr' prefix (default False).
        If True, detects style, then strictly follows keep_chr_prefix to rename if needed.
    keep_chr_prefix : bool
        Only effective if check_chr_prefix=True.
        If input has 'chr' and keep_chr_prefix=False -> Remove prefix (chr1->1).
        If input has no 'chr' and keep_chr_prefix=True -> Add prefix (1->chr1).
    threads : int
        Threads for compression (passed to bcftools --threads).
    out_prefix : Optional[str]
        Output file prefix (without extension). If None, auto-generated.
    out_dir : Optional[str]
        Output directory.
    index_output : bool
        Whether to tabix index the output .vcf.gz (default True).
    index_type : {"tbi", "csi"}
        Index type. Default tbi.
    overwrite : bool
        Whether to overwrite existing output.
    dry_run : bool
        If True, only return the command to be executed (for debugging).

    Returns
    -------
    str
        Path to the filtered .vcf.gz file.

    Notes and Performance Tips
    --------------------------
    1. Uses INFO expression filtering (-i), no decompression needed.
    2. Parallel compression via --threads.
    3. Output uses atomic write (write to tmp, then move).
    """

    # ---------- Argument and Path Validation ----------
    if not os.path.exists(vcf_path):
        raise FileNotFoundError(f"VCF file not found: {vcf_path}")
    
    # Check executables (basic existence check)
    if shutil.which(bcftools_path) is None and not os.path.exists(bcftools_path):
        raise FileNotFoundError(f"bcftools not found: {bcftools_path}")
    
    if index_output:
        if shutil.which(tabix_path) is None and not os.path.exists(tabix_path):
             # Try find tabix in path if not absolute
             if shutil.which("tabix"):
                 tabix_path = "tabix"
             else:
                 raise FileNotFoundError(f"tabix not found: {tabix_path}")

    # ---------- Logging Initialization ----------
    log_dir = os.getcwd()
    log_file = os.path.join(log_dir, f"filter_vcf_by_info_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log")

    # Use independent logger
    logger = logging.getLogger("info_filter_tools.filter_vcf_by_info")
    logger.setLevel(logging.DEBUG)
    # Avoid duplicate handlers
    for h in list(logger.handlers):
        logger.removeHandler(h)
    fh = logging.FileHandler(log_file, encoding='utf-8')
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(logging.Formatter('%(asctime)s [%(levelname)s] %(message)s'))
    logger.addHandler(fh)
    logger.propagate = False

    logger.info("=== Start filter_vcf_by_info ===")
    logger.info(f"Input file: {vcf_path}")
    logger.info(f"Output prefix: {out_prefix}")
    logger.info(f"INFO key: {info_key}")
    logger.info(f"Values: {values}")
    logger.info(f"Logic: {logic}, Match mode: {match_mode}")
    logger.info(f"Threads: {threads}")
    logger.info(f"check_chr_prefix: {check_chr_prefix}, keep_chr_prefix: {keep_chr_prefix}")

    # ---------- (Optional) Unify CHROM Prefix ----------
    normalized_vcf_path = vcf_path
    tmp_map_file = None
    tmp_norm_vcf = None

    if check_chr_prefix:
        detected = _detect_chr_prefix(vcf_path, bcftools_path)
        logger.info(f"CHROM prefix detection: {detected} (True='chr' present, False=absent, None=undetermined)")
        action = None
        if detected is True and not keep_chr_prefix:
            action = "remove_chr_prefix"
        elif detected is False and keep_chr_prefix:
            action = "add_chr_prefix"

        if action is None:
            logger.info("No CHROM prefix modification needed.")
        else:
            # Construct rename-chrs map file
            mapping_lines = []
            if action == "remove_chr_prefix":
                for i in range(1, 23):
                    mapping_lines.append(f"chr{i}\t{i}\n")
                mapping_lines.extend([
                    "chrX\tX\n",
                    "chrY\tY\n",
                    "chrM\tM\n",
                ])
            elif action == "add_chr_prefix":
                for i in range(1, 23):
                    mapping_lines.append(f"{i}\tchr{i}\n")
                mapping_lines.extend([
                    "X\tchrX\n",
                    "Y\tchrY\n",
                    "M\tchrM\n",
                ])

            # Write temp map file
            tf = tempfile.NamedTemporaryFile("w", delete=False, prefix="rename_chrs_", suffix=".tsv")
            tf.writelines(mapping_lines)
            tf.flush()
            tf.close()
            tmp_map_file = tf.name
            logger.info(f"Generated rename-chrs map: {tmp_map_file}, {len(mapping_lines)} lines; action={action}")

            # Generate temp normalized VCF (bgzip compressed)
            in_dir = os.path.dirname(os.path.abspath(vcf_path))
            tmp_norm_vcf = os.path.join(in_dir, f"chrfix.{os.path.basename(vcf_path)}")
            ann_cmd = [
                bcftools_path,
                "annotate",
                "--rename-chrs",
                tmp_map_file,
                "-Oz",
                "-o",
                tmp_norm_vcf,
                "--threads",
                str(int(threads) if threads and threads > 0 else 1),
                vcf_path,
            ]
            logger.info("bcftools (normalize CHROM) command: " + " ".join(ann_cmd))

            if dry_run:
                chrfix_cmd_str = " ".join(shlex.quote(c) for c in ann_cmd)
            else:
                try:
                    if os.path.exists(tmp_norm_vcf):
                        os.remove(tmp_norm_vcf)
                    subprocess.run(ann_cmd, check=True)
                    normalized_vcf_path = tmp_norm_vcf
                    logger.info(f"CHROM normalization complete: {normalized_vcf_path}")
                except subprocess.CalledProcessError as e:
                    logger.error(f"CHROM normalization failed: {e}. Command: {' '.join(ann_cmd)}")
                    # Cleanup
                    if os.path.exists(tmp_norm_vcf):
                        try:
                            os.remove(tmp_norm_vcf)
                        except Exception:
                            pass
                    raise RuntimeError(
                        f"bcftools annotate --rename-chrs failed (exit={e.returncode})."
                    ) from e

    # Normalize info_key
    if not info_key.startswith("INFO/"):
        info_key = f"INFO/{info_key}"

    # Verify INFO key in header
    key_name = info_key.split("/", 1)[1]
    try:
        info_ids = _collect_info_ids_from_header(normalized_vcf_path, bcftools_path)
    except subprocess.CalledProcessError as e:
        raise RuntimeError(
            f"Failed to read VCF header for INFO validation. exit={e.returncode}"
        ) from e

    if key_name not in info_ids:
        suggestions = difflib.get_close_matches(key_name, sorted(info_ids), n=5, cutoff=0.6)
        hint = ("; Did you mean: " + ", ".join(suggestions)) if suggestions else ""
        raise ValueError(
            f"INFO field '{key_name}' not found in VCF header{hint}. Available fields: "
            + ", ".join(list(sorted(info_ids))[:10])
            + (" ..." if len(info_ids) > 10 else "")
        )

    # Normalize values list
    if isinstance(values, (str, bytes)):
        values_list: List[str] = [str(values)]
    else:
        values_list = [str(v) for v in values]
    if len(values_list) == 0:
        raise ValueError("values must contain at least one value")

    # ---------- Build Filter Expression ----------
    # bcftools expression e.g.: INFO/impact=="HIGH" || INFO/impact=="MODERATE"
    ops = {
        "exact": lambda v: f'{info_key}=="{v}"',
        "contains": lambda v: f'{info_key} ~ "{re.escape(v)}"',
        "regex": lambda v: f'{info_key} ~ "{v}"',
    }
    if match_mode not in ops:
        raise ValueError("match_mode must be 'exact' | 'contains' | 'regex'")

    clauses = [ops[match_mode](v) for v in values_list]
    joiner = " || " if logic == "any" else " && "
    expr = joiner.join(clauses)

    logger.info(f"Constructed filter expression: {expr}")

    # ---------- Output Path ----------
    in_dir = os.path.dirname(os.path.abspath(vcf_path))
    auto_name = (
        f"infofilter_{info_key.split('/',1)[1]}_" + "-".join(values_list)
    )
    base_prefix = out_prefix if out_prefix else auto_name

    # Handle output directory
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
        if os.path.basename(base_prefix) == base_prefix:
            base_prefix = os.path.join(out_dir, base_prefix)
    else:
        # Default to input directory if only filename provided
        if os.path.basename(base_prefix) == base_prefix:
            base_prefix = os.path.join(in_dir, base_prefix)

    out_vcf = f"{base_prefix}.vcf.gz"
    tmp_vcf = f"{out_vcf}.tmp"

    if os.path.exists(out_vcf) and not overwrite:
        raise FileExistsError(f"Output exists and overwrite=False: {out_vcf}")

    # ---------- Build Command ----------
    cmd = [
        bcftools_path,
        "view",
        "-i",
        expr,
        "-Oz",
        "-o",
        tmp_vcf,
        "--threads",
        str(int(threads) if threads and threads > 0 else 1),
        normalized_vcf_path,
    ]

    logger.info(f"bcftools command: {' '.join(cmd)}")

    # ---------- dry-run ----------
    if dry_run:
        filter_cmd_str = " ".join(shlex.quote(c) for c in cmd)
        if check_chr_prefix and 'chrfix_cmd_str' in locals():
            return chrfix_cmd_str + " && " + filter_cmd_str # type: ignore[return-value]
        return filter_cmd_str

    # ---------- Execute Filtering ----------
    try:
        if os.path.exists(tmp_vcf):
            os.remove(tmp_vcf)
        subprocess.run(cmd, check=True)
        # Atomic replace
        os.replace(tmp_vcf, out_vcf)
        logger.info(f"Filter complete. Output: {out_vcf}")
    except subprocess.CalledProcessError as e:
        logger.error(f"Execution failed: {e}. Command: {' '.join(cmd)}")
        if os.path.exists(tmp_vcf):
            try:
                os.remove(tmp_vcf)
            except Exception:
                pass
        raise RuntimeError(
            f"bcftools view filtering failed (exit={e.returncode})."
        ) from e

    # ---------- Indexing ----------
    if index_output:
        idx_cmd = [tabix_path, "-f", "-p", "vcf"]
        if index_type.lower() == "csi":
            idx_cmd.append("-C")
        idx_cmd.append(out_vcf)
        try:
            subprocess.run(idx_cmd, check=True)
            logger.info(f"Indexing complete: {out_vcf}.{index_type}")
            logger.info("=== Finished successfully ===")
        except subprocess.CalledProcessError as e:
            raise RuntimeError(
                f"tabix indexing failed (exit={e.returncode})."
            ) from e

    # ---------- Temporary File Cleanup ----------
    try:
        if tmp_map_file and os.path.exists(tmp_map_file):
            os.remove(tmp_map_file)
        if tmp_norm_vcf and os.path.exists(tmp_norm_vcf):
            try:
                os.remove(tmp_norm_vcf)
                if os.path.exists(tmp_norm_vcf + ".tbi"):
                    os.remove(tmp_norm_vcf + ".tbi")
            except Exception:
                pass
    except Exception:
        logger.warning("Non-fatal error cleaning up temp files", exc_info=True)

    return out_vcf
