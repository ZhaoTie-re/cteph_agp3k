"""
extract_vcf_format.py
本脚本用于从两个VCF文件（如WGS与Array）中提取指定染色体上共享变异位点和样本的FORMAT字段矩阵，适用于基因型一致性评估等分析场景。
主要功能：
1. 并行提取两个VCF文件在指定染色体上的变异ID集合，并找出共享变异ID。
2. 提取两个VCF文件的样本ID集合，并找出共享样本ID。
3. 将共享变异ID和样本ID分别输出到文本文件。
4. 按块（chunk）方式，利用bcftools高效提取共享变异和样本的指定FORMAT字段（如GT、DP、GQ、AF等），并输出为压缩的矩阵文件（tsv.gz）。
参数说明：
- --call_prefix: 结果输出文件的前缀（如wgs）。
- --true_prefix: 结果输出文件的前缀（如array）。
- --call_vcf: 输入的call VCF文件路径（必需）。
- --true_vcf: 输入的true VCF文件路径（必需）。
- --chr: 需要处理的染色体（必需）。
- --max_variants: 最多处理的变异位点数（调试可选，默认全部）。
- --num_threads: bcftools并行线程数（默认4）。
- --call_format_fields: 需要提取的FORMAT字段列表（默认["GT", "DP", "GQ", "AF"]）。
输出文件：
- {chr}.shared_variant_ids.txt: 共享变异ID列表。
- {chr}.shared_samples.txt: 共享样本ID列表。
- {chr}.{call_prefix}.{field}.tsv.gz: call VCF中每个FORMAT字段的样本×变异矩阵。
- {chr}.{true_prefix}.GT.tsv.gz: true VCF中GT字段的样本×变异矩阵。
依赖：
- bcftools
- numpy
- Python标准库（concurrent.futures, subprocess, argparse, tempfile, os, gc, shutil）
用法示例：
python extract_vcf_format.py --call_vcf call.vcf.gz --true_vcf truth.vcf.gz --chr chr1 --call_prefix wgs --true_prefix array --num_threads 8 --call_format_fields GT DP GQ AF
"""
from concurrent.futures import ThreadPoolExecutor
import numpy as np
import subprocess
from typing import List
import subprocess
import argparse
import tempfile
import os
import gc
import shutil
import gzip

# %%
# 设置参数
parser = argparse.ArgumentParser(description="Extract VCF format fields for shared variants and samples.")
parser.add_argument("--call_prefix", type=str, default="wgs", help="Prefix for call VCF files (default: 'wgs').")
parser.add_argument("--true_prefix", type=str, default="array", help="Prefix for true VCF files (default: 'array').")
parser.add_argument("--call_vcf", type=str, required=True, help="Path to the call VCF file.")
parser.add_argument("--true_vcf", type=str, required=True, help="Path to the true VCF file.")
parser.add_argument("--chr", type=str, required=True, help="Chromosome to process.")
parser.add_argument("--max_variants", type=int, default=None, help="Maximum number of variants to process (default: None, process all).")
parser.add_argument(
    "--num_threads", 
    type=int, 
    default=4, 
    help="Number of threads to use for bcftools (default: 4)."
)
parser.add_argument(
    "--call_format_fields", 
    type=str, 
    nargs='+', 
    default=["GT", "DP", "GQ", "AF"],
    help="List of FORMAT fields to extract (e.g. --call_format_fields GT DP GQ AF)"
)
args = parser.parse_args()

# %%
call_prefix = args.call_prefix
true_prefix = args.true_prefix

call_vcf_path = args.call_vcf
true_vcf_path = args.true_vcf

chr = args.chr

# %%

from typing import Optional

def extract_variant_ids(vcf_path: str, chr: str, max_variants: Optional[int] = None) -> set:
    """提取 VCF 中 chr 上的所有 ID 字段为集合"""
    query_cmd = [
        "bcftools", "query",
        "-r", chr,
        "-f", "%ID\n",
        vcf_path
    ]
    proc = subprocess.Popen(query_cmd, stdout=subprocess.PIPE, text=True, bufsize=8192)

    variant_ids = set()
    if proc.stdout is not None:
        for i, line in enumerate(proc.stdout):
            if max_variants is not None and i >= max_variants:
                break
            line = line.rstrip('\n')
            if line != ".":
                variant_ids.add(line)
        proc.stdout.close()
    proc.wait()
    return variant_ids


def extract_sample_ids(vcf_path: str) -> list:
    """提取 VCF 中的样本 ID 列表"""
    header = subprocess.check_output(["bcftools", "view", "-h", vcf_path], text=True)
    for line in header.strip().splitlines():
        if line.startswith("#CHROM"):
            return line.strip().split("\t")[9:]
    return []


# === 设置参数 ===

max_variants = args.max_variants

# === 并行提取 ID ===
from concurrent.futures import wait

with ThreadPoolExecutor(max_workers=2) as executor:
    futures = {
        "true": executor.submit(extract_variant_ids, true_vcf_path, chr, max_variants),
        "call": executor.submit(extract_variant_ids, call_vcf_path, chr, max_variants)
    }
    done, not_done = wait(futures.values())
    for fut in done:
        exc = fut.exception()
        if exc is not None:
            raise exc
    true_ids = futures["true"].result()
    call_ids = futures["call"].result()

# === 提取 shared IDs 和 shared samples ===
shared_variant_ids = true_ids & call_ids
shared_variant_ids = sorted(shared_variant_ids, key=lambda x: int(x.split(":")[1]))
true_samples = extract_sample_ids(true_vcf_path)
call_samples = extract_sample_ids(call_vcf_path)
shared_samples = sorted(set(true_samples) & set(call_samples))  # 保证顺序


# === 可选输出检查 ===
print(f"Shared samples: {len(shared_samples)}")
print(f"Shared variant IDs on {chr}: {len(shared_variant_ids)}")

# === 显式释放内存 ===
del true_ids, call_ids, true_samples, call_samples
gc.collect()

# %%
# === 输出共享变异ID和样本ID到文件 ===

with open(f"{chr}.shared_variant_ids.txt", "w") as f:
    for vid in shared_variant_ids:
        f.write(vid + "\n")

with open(f"{chr}.shared_samples.txt", "w") as f:
    for sid in shared_samples:
        f.write(sid + "\n")

# %%
def extract_format_matrices_to_disk(
    vcf_path: str,
    sample_ids: List[str],
    focus_loci: List[str],
    format_fields_to_extract: List[str] = ["GT"],
    use_gt_map: bool = True,
    num_threads: int = 4
) -> dict:
    """
    从VCF文件中提取指定样本和变异位点的FORMAT字段矩阵，并返回字典结果。

    参数:
        vcf_path (str): 输入VCF文件路径。
        sample_ids (List[str]): 需要提取的样本ID列表。
        focus_loci (List[str]): 需要提取的变异位点ID列表，格式为 'chr:pos:ref:alt'。
        format_fields_to_extract (List[str], 可选): 需要提取的FORMAT字段名列表，默认["GT"]。
        use_gt_map (bool, 可选): 是否将GT字段转换为数字编码，默认True。
        num_threads (int, 可选): bcftools使用的线程数，默认4。

    返回:
        dict: {field_name: {sample_id: [values按locus顺序]}}
    """
    tmp_dir = tempfile.mkdtemp(prefix="matrix_tmp_")

    # VCF字段索引
    VCF_COLUMNS = {
        "CHROM": 0,
        "POS": 1,
        "ID": 2,
        "REF": 3,
        "ALT": 4,
        "QUAL": 5,
        "FILTER": 6,
        "INFO": 7,
        "FORMAT": 8,
        "SAMPLES_START": 9
    }

    # GT字段转换为数字
    gt_map = {
        "0/0": 0, "0|0": 0,
        "0/1": 1, "1/0": 1, "0|1": 1, "1|0": 1,
        "1/1": 2, "1|1": 2,
        "./.": np.nan, ".|.": np.nan
    }

    # 构建查询区域列表与定位字典
    loci_dict = {}
    for locus in focus_loci:
        chrom, pos, ref, alt = locus.split(":")
        key = (chrom, int(pos), ref, alt)
        loci_dict[key] = locus
    
    region_file_path = os.path.join(tmp_dir, "regions.txt")
    with open(region_file_path, "w") as f:
        for locus in focus_loci:
            chrom, pos, ref, alt = locus.split(":")
            f.write(f"{chrom}\t{pos}\t{pos}\n")
    
    # 获取VCF中样本名
    header = subprocess.check_output(["bcftools", "view", "-h", vcf_path], text=True)
    vcf_samples = None
    for line in header.strip().splitlines():
        if line.startswith("#CHROM"):
            vcf_samples = line.strip().split('\t')[VCF_COLUMNS["SAMPLES_START"]:]
            break
    if vcf_samples is None:
        raise RuntimeError("Could not find #CHROM header line in VCF to extract sample IDs.")

    # 初始化结果字典
    result = {field: {s: [] for s in sample_ids} for field in format_fields_to_extract}

    # bcftools 命令
    cmd = [
        "bcftools", "view",
        "--threads", str(num_threads),
        "-R", region_file_path,
        vcf_path
    ]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True)

    if proc.stdout is not None:
        for line in proc.stdout:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            chrom = fields[VCF_COLUMNS["CHROM"]]
            pos = int(fields[VCF_COLUMNS["POS"]])
            ref = fields[VCF_COLUMNS["REF"]]
            alt = fields[VCF_COLUMNS["ALT"]].split(',')[0]
            key = (chrom, pos, ref, alt)
            if key not in loci_dict:
                continue
            locus_id = loci_dict[key]

            format_keys = fields[VCF_COLUMNS["FORMAT"]].split(":")
            format_index = {k: i for i, k in enumerate(format_keys)}
            sample_fields = fields[VCF_COLUMNS["SAMPLES_START"]:]

            for field in format_fields_to_extract:
                values = []
                for sample_str in sample_fields:
                    parts = sample_str.split(":")
                    raw_value = parts[format_index[field]] if field in format_index and len(parts) > format_index[field] else None
                    if field == "GT":
                        key = raw_value if raw_value is not None else "."
                        val = gt_map.get(key, np.nan) if use_gt_map else raw_value
                    elif field == "AD":
                        try:
                            val = tuple(map(int, raw_value.split(","))) if raw_value and "," in raw_value else np.nan
                        except:
                            val = np.nan
                    else:
                        try:
                            val = float(raw_value) if (raw_value is not None and raw_value != ".") else np.nan
                        except:
                            val = np.nan
                    values.append(val)
                per_sample = dict(zip(vcf_samples, values))
                for s in sample_ids:
                    result[field][s].append(per_sample.get(s, np.nan))

        proc.stdout.close()
    proc.wait()

    shutil.rmtree(tmp_dir)

    return result

def write_matrix_chunk_v2(output_path, sample_ids, loci_chunk, data_dict, append=False):
    mode = "at" if append else "wt"
    with gzip.open(output_path, mode) as out:
        if not append:
            out.write("sample_id\t" + "\t".join(loci_chunk) + "\n")
        for sid in sample_ids:
            row = [sid]
            values = data_dict.get(sid, [])
            for v in values:
                if isinstance(v, tuple):
                    val_str = ",".join(map(str, v))
                elif v is np.nan or (isinstance(v, float) and np.isnan(v)):
                    val_str = ""
                else:
                    val_str = str(v)
                row.append(val_str)
            out.write("\t".join(row) + "\n")

def merge_matrix_chunks(output_path, sample_ids, all_loci, chunk_paths):
    with gzip.open(output_path, "wt") as out:
        out.write("sample_id\t" + "\t".join(all_loci) + "\n") #写入所有列名
        for sid in sample_ids:
            row = [sid]
            for tmp_path in chunk_paths:
                with gzip.open(tmp_path, "rt") as f:
                    header = f.readline()
                    found = False
                    for line in f:
                        if line.startswith(sid + "\t"): # 找到对应样本的行
                            row.extend(line.strip().split("\t")[1:]) # 添加对应样本的值
                            found = True
                            break
                    if not found:
                        row.extend([""] * (len(header.strip().split("\t")) - 1))
            out.write("\t".join(row) + "\n") # 写入每个样本的行数据


# === 并行块处理与合并优化实现 ===
from concurrent.futures import ThreadPoolExecutor, as_completed

def process_vcf_in_chunks_parallel(vcf_path, sample_ids, variant_ids, output_prefix, format_fields_to_extract, num_threads):
    chunk_size = 10000
    n = len(variant_ids)
    tmp_dir = tempfile.mkdtemp(prefix="matrix_chunk_")

    chunk_files_by_field = {field: {} for field in format_fields_to_extract}
    loci_chunks = {}

    def process_chunk(i, loci_chunk):
        data = extract_format_matrices_to_disk(
            vcf_path=vcf_path,
            sample_ids=sample_ids,
            focus_loci=loci_chunk,
            format_fields_to_extract=format_fields_to_extract,
            use_gt_map=True,
            num_threads=num_threads
        )
        paths = {}
        for field in format_fields_to_extract:
            tmp_path = os.path.join(tmp_dir, f"{field}.chunk{i}.tmp.gz")
            write_matrix_chunk_v2(tmp_path, sample_ids, loci_chunk, data[field], append=False)
            paths[field] = tmp_path
        return i, loci_chunk, paths

    with ThreadPoolExecutor(max_workers=num_threads) as executor:
        futures = []
        for i in range(0, n, chunk_size):
            loci_chunk = variant_ids[i:i+chunk_size]
            loci_chunks[i] = loci_chunk
            futures.append(executor.submit(process_chunk, i, loci_chunk))
        for future in as_completed(futures):
            i, loci_chunk, paths = future.result()
            for field in format_fields_to_extract:
                chunk_files_by_field[field][i] = paths[field]

    def merge_chunks_for_field(field):
        output_path = f"{chr}.{output_prefix}.{field}.tsv.gz"
        sorted_indices = sorted(chunk_files_by_field[field].keys())
        chunk_paths = [chunk_files_by_field[field][i] for i in sorted_indices]
        for tmp_path in chunk_paths:
            if not os.path.exists(tmp_path):
                raise FileNotFoundError(f"Missing chunk file: {tmp_path}")
        all_loci = []
        for i in sorted_indices:
            all_loci.extend(loci_chunks[i])
        merge_matrix_chunks(output_path, sample_ids, all_loci, chunk_paths)

    merge_threads = min(num_threads, len(format_fields_to_extract))
    with ThreadPoolExecutor(max_workers=merge_threads) as merge_executor:
        merge_futures = [merge_executor.submit(merge_chunks_for_field, field) for field in format_fields_to_extract]
        for f in as_completed(merge_futures):
            f.result()

    shutil.rmtree(tmp_dir)

# 并发处理call_vcf和true_vcf
with ThreadPoolExecutor(max_workers=2) as outer_executor:
    f_call = outer_executor.submit(
        process_vcf_in_chunks_parallel,
        call_vcf_path,
        shared_samples,
        shared_variant_ids,
        call_prefix,
        args.call_format_fields,
        args.num_threads
    )
    f_true = outer_executor.submit(
        process_vcf_in_chunks_parallel,
        true_vcf_path,
        shared_samples,
        shared_variant_ids,
        true_prefix,
        ["GT"],
        args.num_threads
    )
    f_call.result()
    f_true.result()