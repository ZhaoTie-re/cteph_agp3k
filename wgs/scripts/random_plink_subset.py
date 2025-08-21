#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
工具名称：random_plink_subset
================================
【功能】
- 从给定的 PLINK 二进制数据（.bed/.bim/.fam）中，随机抽取指定数量的样本与变体，
  并用 plink2 生成一个新的子集数据集（.bed/.bim/.fam）。

- 本模块暴露了一个可调用函数 `random_plink_subset`，方便从其他 Python 代码中调用。
- 已移除命令行接口部分，改为函数调用。

【函数接口示例】
```python
from random_plink_subset import random_plink_subset

result = random_plink_subset(
    bed_prefix="data/mydata",
    n_samples=200,
    n_variants=5000,
    output_prefix="data/mydata_subset",
    seed=12345,
    plink2_path="/usr/local/bin/plink2",
    threads=4,
    tmpdir="/tmp"
)
print("生成文件路径:", result)
```

【输入参数】
- bed_prefix:  原始 PLINK 前缀（不带扩展名）
- n_samples:   随机抽取的样本数（默认 300）
- n_variants:  随机抽取的变体数（默认 10000）
- output_prefix: 输出前缀（默认在原前缀后加 `.rand300x10000`）
- seed:        随机种子（默认 20250820，确保可复现）
- plink2_path: plink2 可执行文件路径（默认 "/home/b/b37974/plink2"）
- threads:     线程数（默认 8）
- tmpdir:      临时目录，默认 None 表示自动创建 /tmp/plink_subset_user_<uuid>

【输出】
- {output_prefix}.bed/.bim/.fam/.log
- 返回包含输出文件路径的字典

【注意】
- 若数据中实际样本/变体数量少于请求数量，则自动取“能取到的最大值”（不报错）。
- 临时文件会创建于 tmpdir 并在结束时清理。
"""

import os
import sys
import random
import subprocess
from pathlib import Path
import shutil
import uuid
from typing import Optional, Dict

def read_fam_iids(fam_path: str):
    """读取 .fam（6 列，空白分隔），返回 [(FID, IID), ...] 列表"""
    iids = []
    with open(fam_path, "r") as f:
        for line in f:
            if not line.strip():
                continue
            parts = line.strip().split()
            # FID, IID 位于前两列
            iids.append((parts[0], parts[1]))
    return iids

def read_bim_variant_ids(bim_path: str):
    """读取 .bim（6 列，空白分隔），返回变体ID列表（第2列SNP ID）"""
    vids = []
    with open(bim_path, "r") as f:
        for line in f:
            if not line.strip():
                continue
            parts = line.strip().split()
            vids.append(parts[1])  # 第2列：variant/SNP ID
    return vids

def write_keep(path: str, pairs):
    """写 FID IID 两列（无表头）"""
    with open(path, "w") as f:
        for fid, iid in pairs:
            f.write(f"{fid}\t{iid}\n")

def write_extract(path: str, vids):
    """写每行一个变体ID（无表头）"""
    with open(path, "w") as f:
        for vid in vids:
            f.write(f"{vid}\n")

def run_plink_subset(plink2: str, bed_prefix: str, keep_path: str, extract_path: str,
                     out_prefix: str, threads: int, tmp_dir: str):
    # Use the provided tmp_dir for intermediate files and plink output
    # Output files will be generated under tmp_dir, but with the same out_prefix basename
    out_basename = os.path.basename(out_prefix)
    tmp_out_prefix = os.path.join(tmp_dir, out_basename)

    cmd = [
        plink2,
        "--bfile", bed_prefix,
        "--keep", keep_path,
        "--extract", extract_path,
        "--make-bed",
        "--out", tmp_out_prefix,
        "--threads", str(threads),
        "--allow-no-sex"
    ]
    print("[INFO] Running:", " ".join(cmd))
    subprocess.run(cmd, check=True)

    # Move only .bed, .bim, .fam, .log to current directory with requested out_prefix
    for ext in [".bed", ".bim", ".fam", ".log"]:
        src = tmp_out_prefix + ext
        dst = out_prefix + ext
        if os.path.exists(src):
            shutil.move(src, dst)
        else:
            print(f"[WARN] Expected file not found: {src}")

def random_plink_subset(
    bed_prefix: str,
    n_samples: int = 300,
    n_variants: int = 10000,
    output_prefix: Optional[str] = None,
    seed: int = 20250820,
    plink2_path: str = "/home/b/b37974/plink2",
    threads: int = 8,
    tmpdir: Optional[str] = None,
) -> Dict[str, str]:
    """
    从给定 PLINK 前缀中随机抽取 n_samples 个样本与 n_variants 个变体，
    仅在工作目录生成 {output_prefix}.bed/.bim/.fam/.log，其它文件写入 tmpdir 并在结束时清理。

    返回：包含输出前缀、四个关键文件路径的字典。
    """
    random.seed(seed)

    fam_path = bed_prefix + ".fam"
    bim_path = bed_prefix + ".bim"
    bed_path = bed_prefix + ".bed"

    # 基础检查
    for p in [fam_path, bim_path, bed_path]:
        if not Path(p).exists():
            raise FileNotFoundError(f"缺少文件：{p}")

    iids = read_fam_iids(fam_path)
    vids = read_bim_variant_ids(bim_path)
    n_total_samples = len(iids)
    n_total_variants = len(vids)

    if n_total_samples == 0 or n_total_variants == 0:
        raise ValueError("样本或变体数量为0，无法抽样。")

    n_take_samples = min(n_samples, n_total_samples)
    n_take_variants = min(n_variants, n_total_variants)

    print(f"[INFO] 总样本数: {n_total_samples}，请求抽取: {n_samples} -> 实际抽取: {n_take_samples}")
    print(f"[INFO] 总变体数: {n_total_variants}，请求抽取: {n_variants} -> 实际抽取: {n_take_variants}")
    print(f"[INFO] 随机种子: {seed}")

    sampled_iids = random.sample(iids, n_take_samples)
    sampled_vids = random.sample(vids, n_take_variants)

    out_prefix = output_prefix or f"{bed_prefix}.rand{n_take_samples}x{n_take_variants}"

    # Prepare tmp_dir
    if tmpdir is None:
        # Generate unique tmpdir under /tmp with user and uuid
        user = os.getenv("USER") or "user"
        unique_id = str(uuid.uuid4())
        tmp_dir = os.path.join("/tmp", f"plink_subset_{user}_{unique_id}")
    else:
        tmp_dir = os.path.abspath(tmpdir)
    os.makedirs(tmp_dir, exist_ok=True)

    out_basename = os.path.basename(out_prefix)
    keep_path = os.path.join(tmp_dir, f"{out_basename}.keep.txt")
    extract_path = os.path.join(tmp_dir, f"{out_basename}.extract.txt")

    write_keep(keep_path, sampled_iids)
    write_extract(extract_path, sampled_vids)

    run_plink_subset(
        plink2=plink2_path,
        bed_prefix=bed_prefix,
        keep_path=keep_path,
        extract_path=extract_path,
        out_prefix=out_prefix,
        threads=threads,
        tmp_dir=tmp_dir,
    )

    # Cleanup tmp_dir and its contents except the moved files
    try:
        shutil.rmtree(tmp_dir)
    except Exception as e:
        print(f"[WARN] Failed to remove temp dir {tmp_dir}: {e}", file=sys.stderr)

    return {
        "output_prefix": out_prefix,
        "bed": out_prefix + ".bed",
        "bim": out_prefix + ".bim",
        "fam": out_prefix + ".fam",
        "log": out_prefix + ".log",
    }