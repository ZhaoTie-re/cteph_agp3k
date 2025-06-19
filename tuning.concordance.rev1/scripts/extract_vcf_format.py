# %%
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

def extract_variant_ids(vcf_path: str, chr: str, max_variants: int = None) -> set:
    """提取 VCF 中 chr 上的所有 ID 字段为集合"""
    query_cmd = [
        "bcftools", "query",
        "-r", chr,
        "-f", "%ID\n",
        vcf_path
    ]
    proc = subprocess.Popen(query_cmd, stdout=subprocess.PIPE, text=True, bufsize=8192)

    variant_ids = set()
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
with ThreadPoolExecutor(max_workers=2) as executor:
    future_true = executor.submit(extract_variant_ids, true_vcf_path, chr, max_variants)
    future_call = executor.submit(extract_variant_ids, call_vcf_path, chr, max_variants)
    true_ids = future_true.result()
    call_ids = future_call.result()

# === 提取 shared IDs 和 shared samples ===
shared_variant_ids = true_ids & call_ids
true_samples = extract_sample_ids(true_vcf_path)
call_samples = extract_sample_ids(call_vcf_path)
shared_samples = set(true_samples) & set(call_samples)

# === 可选输出检查 ===
print(f"Shared samples: {len(shared_samples)}")
print(f"Shared variant IDs on {chr}: {len(shared_variant_ids)}")

# === 显式释放内存 ===
del true_ids, call_ids, true_samples, call_samples
gc.collect()

# %%
# === 排序并转换为列表 ===
shared_variant_ids = sorted(shared_variant_ids)
shared_samples = sorted(shared_samples)

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
    for line in header.strip().splitlines():
        if line.startswith("#CHROM"):
            vcf_samples = line.strip().split('\t')[VCF_COLUMNS["SAMPLES_START"]:]
            break

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
                    val = gt_map.get(raw_value, np.nan) if use_gt_map else raw_value
                elif field == "AD":
                    try:
                        val = tuple(map(int, raw_value.split(","))) if raw_value and "," in raw_value else np.nan
                    except:
                        val = np.nan
                else:
                    try:
                        val = float(raw_value) if raw_value not in [None, "."] else np.nan
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
    import gzip
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

def process_vcf_in_chunks(vcf_path, sample_ids, variant_ids, output_prefix, format_fields_to_extract, num_threads):
    chunk_size = 10000
    n = len(variant_ids)
    for field in format_fields_to_extract:
        output_path = f"{chr}.{output_prefix}.{field}.tsv.gz"
        # Remove file if exists to start fresh
        if os.path.exists(output_path):
            os.remove(output_path)

    for i in range(0, n, chunk_size):
        loci_chunk = variant_ids[i:i+chunk_size]
        data = extract_format_matrices_to_disk(
            vcf_path=vcf_path,
            sample_ids=sample_ids,
            focus_loci=loci_chunk,
            format_fields_to_extract=format_fields_to_extract,
            use_gt_map=True,
            num_threads=num_threads
        )
        # write chunk to output files
        for idx_field, field in enumerate(format_fields_to_extract):
            output_path = f"{chr}.{output_prefix}.{field}.tsv.gz"
            append = (i != 0)
            write_matrix_chunk_v2(output_path, sample_ids, loci_chunk, data[field], append=append)
        print(f"[✓] Processed chunk {i//chunk_size + 1} / {(n + chunk_size - 1)//chunk_size} for {output_prefix}")

# %%
process_vcf_in_chunks(
    vcf_path=call_vcf_path,
    sample_ids=shared_samples,
    variant_ids=shared_variant_ids,
    output_prefix=call_prefix,
    format_fields_to_extract=args.call_format_fields,
    num_threads=args.num_threads
)

# %%
process_vcf_in_chunks(
    vcf_path=true_vcf_path,
    sample_ids=shared_samples,
    variant_ids=shared_variant_ids,
    output_prefix=true_prefix,
    format_fields_to_extract=["GT"],
    num_threads= args.num_threads
)