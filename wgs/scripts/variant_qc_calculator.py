"""
变异质控计算工具集（variant_qc_calculator.py）

【工具概述】
本工具用于在大规模基因型数据（PLINK 二进制格式：.bed/.bim/.fam）上计算变异层面的质量控制（QC）指标，自动汇总并支持后续一键排除低质量变异，生成新的 PLINK 数据集。流程以 plink2 为核心，辅以 pandas 的分块（streaming）与多进程处理，兼顾可复现性与性能。

【适用场景】
- GWAS/候选位点分析前的标准变异级 QC 统计与筛除；
- 快速获取群体与分组（病例/对照）的 AAF/MAF、缺失率（VMISS）、HWE p 值；
- 生成可追溯的 QC 汇总表与可直接用于 --exclude 的变异列表。

【输入与前提】
- 输入：PLINK 二进制文件前缀（bed_prefix），要求 .bed/.bim/.fam 三件套齐全；
- fam 文件的 PHENO 列使用 1=对照、2=病例 的常规编码（缺失或其他编码将导致分组统计为空或异常）；
- 需要可执行的 plink2（二进制路径可通过 plink2_path 指定）。

【核心输出】
1) 变异 QC 汇总表：`{output_prefix}.variant_qc_summary.tsv`
   - 列说明：
     - VARIANT_ID：PLINK 变异 ID（与 .bim 中 ID 一致）
     - MAF：基于全体样本 ALT_FREQS 计算的 minor allele frequency（min(p, 1-p)）
     - VMISS：全体样本缺失率（F_MISS）
     - CASE_AAF：病例组 ALT 等位基因频率（ALT_FREQS）
     - CTRL_AAF：对照组 ALT 等位基因频率（ALT_FREQS）
     - CTRL_MAF：对照组 MAF（min(CTRL_AAF, 1-CTRL_AAF)）
     - CASE_HWE：病例组 Hardy–Weinberg 平衡检验 p 值（plink2 --hardy）
     - CTRL_HWE：对照组 Hardy–Weinberg 平衡检验 p 值（plink2 --hardy）

2) 低质量变异标记与 ID 列表（由 `extract_maf0_or_vmiss1_variants_streaming` 生成）
   - `{base_prefix}.with_flags.tsv`：列含 VARIANT_ID、MAF0_FLAG、VMISS1_FLAG、MAF_NA_FLAG、VMISS_NA_FLAG
     （新增 MAF 和 VMISS 的 NA 值检测，方便识别缺失数据）
   - `{base_prefix}.variant_ids.tsv`：仅含 VARIANT_ID，可直接用于 plink2 `--exclude`

3) 过滤后的新数据集（由 `run_plink2_exclude_variants` 生成）
   - 输出前缀：调用时指定的 `output_prefix`，包含新 .bed/.bim/.fam

【中间文件（tmpdir）】
- 全体样本频率/缺失：`all_samples.afreq`、`all_samples.vmiss`
- 分组统计：`case.afreq`、`case.hardy`、`ctrl.afreq`、`ctrl.hardy`
- 并行分块临时结果：`chunk_*.tsv`
> 这些文件用于加速与溯源，流程结束后可按需保留或清理。

【阈值与解读建议（非强制，仅供参考）】
- MAF/VMISS：常见预筛包括 MAF==0 或 VMISS==1 的变异直接剔除；
- HWE（对照组）：常用阈值在 1e-6 到 1e-4 之间，需结合样本量与研究设计；
- 低频位点（0.01–0.05）在 HWE 上更敏感，可视研究策略单独设限；
- 本脚本仅负责“计算与汇总”，阈值落地由下游脚本/分析Notebook控制。

【使用示例】
1. 计算 QC 汇总：
   ```python
   out_tsv = run_plink2_variant_qc(
       bed_prefix="/path/to/dataset",
       tmpdir="/tmp/variant_qc",
       plink2_path="/home/b/b37974/plink2",
       threads=16,
       output_prefix="cteph_agp3k",
       verbose=True
   )
   ```
2. 提取 MAF=0 或 VMISS=1 及 MAF/VMISS NA 并生成 --exclude 列表：
   ```python
   full_info_path, variant_ids_path = extract_maf0_or_vmiss1_variants_streaming(out_tsv)
   ```
3. 生成剔除后的新数据：
   ```python
   new_prefix = run_plink2_exclude_variants(
       bed_prefix="/path/to/dataset",
       output_prefix="/path/to/dataset.rm_maf0_vmiss1",
       exclude_variants_file=variant_ids_path,
       plink2_path="/home/b/b37974/plink2",
       threads=16
   )
   ```

【性能与并行】
- 频率/缺失/HWE 由 plink2 多线程计算（`--threads`）；
- pandas 按 100k 行分块读取 `.afreq`，并用 `ProcessPoolExecutor` 并行写入临时块；
- 大队列建议将 `threads` 调大，并确保 `tmpdir` 位于高速磁盘（NVMe 本地盘优于网络盘）。

【错误处理与常见问题】
- fam 的 PHENO 若非 1/2 编码，将导致 case/ctrl 为空；
- 变异 ID 需与 .bim 一致；不同软件生成的数据集请先统一 ID 命名规则；
- 若 plink2 返回非零退出码，请优先检查输入路径/权限/磁盘可写与线程数设置。

【依赖】
- 运行环境：Python 3.9+（pandas、numpy）；系统可执行：plink2
- 主要库：pandas, numpy, csv, subprocess, tempfile, concurrent.futures, uuid, typing

【可复现性与记录】
- 建议将 `plink2` 版本与调用参数记录在项目 README 或运行日志中；
- 输出文件名包含 `output_prefix`，便于与批量任务的样本队列对齐。

【版权与维护】
- 作者：ZHAO TIE
- 文件：variant_qc_calculator.py

【更新记录】
- 2025-08-12：补充用途说明、字段定义、分块与并行策略、常见阈值参考与使用示例；明确中间文件与错误处理注意事项。
- 2025-08-13：在 `extract_maf0_or_vmiss1_variants_streaming` 中新增对 MAF 和 VMISS NA 值的检测与标记，方便识别缺失数据。
"""

import subprocess
import os
import tempfile
import csv
import math
import uuid
from typing import Optional
import concurrent.futures
import pandas as pd
import numpy as np

_global_lookup = {}

def init_globals_for_chunk(vmiss_dict, case_aaf_dict, ctrl_aaf_dict, case_hwe_dict, ctrl_hwe_dict):
    global _global_lookup
    _global_lookup['vmiss_dict'] = vmiss_dict
    _global_lookup['case_aaf_dict'] = case_aaf_dict
    _global_lookup['ctrl_aaf_dict'] = ctrl_aaf_dict
    _global_lookup['case_hwe_dict'] = case_hwe_dict
    _global_lookup['ctrl_hwe_dict'] = ctrl_hwe_dict

def process_chunk(chunk, idx, tmpdir):
    """
    对变异数据块进行处理，将每个变异的QC指标提取并写入临时文件。

    参数:
        chunk (DataFrame): 当前数据块，包含多个变异的统计信息
        idx (int): 数据块编号，用于命名输出文件
        tmpdir (str): 临时目录路径，用于存放中间结果

    返回:
        str: 写入结果的临时文件路径
    """
    vmiss_dict = _global_lookup.get('vmiss_dict', {})
    case_aaf_dict = _global_lookup.get('case_aaf_dict', {})
    ctrl_aaf_dict = _global_lookup.get('ctrl_aaf_dict', {})
    case_hwe_dict = _global_lookup.get('case_hwe_dict', {})
    ctrl_hwe_dict = _global_lookup.get('ctrl_hwe_dict', {})

    tmp_output = os.path.join(tmpdir, f"chunk_{idx}_{uuid.uuid4().hex}.tsv")
    with open(tmp_output, "w", newline="") as fout:
        writer = csv.writer(fout, delimiter="\t")
        for _, row in chunk.iterrows():
            vid = row["ID"]  # 变异ID
            aaf = row["ALT_FREQS"]  # 等位基因频率
            maf = min(aaf, 1 - aaf) if pd.notnull(aaf) else float("nan")  # 计算MAF
            vmiss = vmiss_dict.get(vid, float("nan"))  # 缺失率VMISS
            case_aaf = case_aaf_dict.get(vid, float("nan"))  # 病例组等位基因频率
            ctrl_aaf = ctrl_aaf_dict.get(vid, float("nan"))  # 对照组等位基因频率
            ctrl_maf = min(ctrl_aaf, 1 - ctrl_aaf) if pd.notnull(ctrl_aaf) else float("nan")  # 对照组MAF
            case_hwe = case_hwe_dict.get(vid, float("nan"))  # 病例组HWE p值
            ctrl_hwe = ctrl_hwe_dict.get(vid, float("nan"))  # 对照组HWE p值
            # 写入列顺序: 变异ID, MAF, VMISS, 病例组AAF, 对照组AAF, 对照组MAF, 病例组HWE p值, 对照组HWE p值
            writer.writerow([vid, maf, vmiss, case_aaf, ctrl_aaf, ctrl_maf, case_hwe, ctrl_hwe])
    return tmp_output

def run_plink2_variant_qc(
    bed_prefix: str,
    tmpdir: str = "/tmp/variant_qc",
    plink2_path: str = "/home/b/b37974/plink2",
    threads: int = 8,
    output_prefix: str = "cteph_agp3k",
    verbose: bool = True
) -> str:
    """
    使用 plink2 对 Plink 格式基因型文件进行变异层面的QC计算。

    输出字段包括：
        - MAF（全部样本）
        - VMISS（全部样本）
        - CASE/CONTROL AAF
        - CONTROL MAF
        - CASE/CONTROL HWE P值

    参数:
        bed_prefix (str): 输入文件的 Plink 数据前缀（.bed/.bim/.fam）
        tmpdir (str): 临时目录路径，用于存放中间结果
        plink2_path (str): plink2 执行路径
        threads (int): 并行使用的线程数
        output_prefix (str): 输出 QC 结果文件的前缀
        verbose (bool): 是否打印进度信息

    返回:
        str: 输出 QC 汇总结果文件的路径（.variant_qc_summary.tsv）
    """

    # 确保 tmpdir 目录存在
    if tmpdir is None:
        tmpdir = tempfile.mkdtemp()
    else:
        os.makedirs(tmpdir, exist_ok=True)
    if verbose:
        print(f"[INFO] 使用临时目录: {tmpdir}")

    def run_cmd(cmd, desc: Optional[str] = None):
        if verbose and desc:
            print(f"[INFO] 运行: {desc}")
        try:
            subprocess.run(cmd, check=True)
        except subprocess.CalledProcessError as e:
            print(f"[ERROR] 命令运行失败: {' '.join(cmd)}")
            raise e

    # Step 1: freq + vmiss
    out_all = os.path.join(tmpdir, "all_samples")
    # 计算全体样本的等位基因频率和缺失率
    run_cmd([
        plink2_path, "--bfile", bed_prefix, "--threads", str(threads),
        "--freq", "--missing", "--out", out_all
    ], "全体样本 AAF + VMISS")

    # Step 2: 分组 IID
    fam_df = pd.read_csv(f"{bed_prefix}.fam", sep=r"\s+", header=None)
    fam_df.columns = ["FID", "IID", "PID", "MID", "SEX", "PHENO"]
    # 根据PHENO字段区分病例和对照样本ID
    case_iids = fam_df[fam_df["PHENO"] == 2][["FID", "IID"]]
    ctrl_iids = fam_df[fam_df["PHENO"] == 1][["FID", "IID"]]

    case_iid_path = os.path.join(tmpdir, "case_iids.txt")
    ctrl_iid_path = os.path.join(tmpdir, "ctrl_iids.txt")
    # 保存病例和对照样本ID，用于plink2的--keep参数
    case_iids.to_csv(case_iid_path, sep="\t", index=False, header=False)
    ctrl_iids.to_csv(ctrl_iid_path, sep="\t", index=False, header=False)

    # Step 3: 分组 freq + HWE
    out_case = os.path.join(tmpdir, "case")
    out_ctrl = os.path.join(tmpdir, "ctrl")
    # 计算病例组的等位基因频率和HWE检验
    run_cmd([
        plink2_path, "--bfile", bed_prefix, "--keep", case_iid_path,
        "--threads", str(threads), "--freq", "--hardy", "--out", out_case
    ], "Case AAF + HWE")

    # 计算对照组的等位基因频率和HWE检验
    run_cmd([
        plink2_path, "--bfile", bed_prefix, "--keep", ctrl_iid_path,
        "--threads", str(threads), "--freq", "--hardy", "--out", out_ctrl
    ], "Control AAF + HWE")

    # Step 4: load lookup tables
    if verbose:
        print("[INFO] 读取辅助统计表...")

    # 读取全体样本缺失率字典
    vmiss_dict = dict(pd.read_csv(out_all + ".vmiss", sep=r"\s+")[["ID", "F_MISS"]].values)
    # 读取病例组等位基因频率字典
    case_aaf_dict = dict(pd.read_csv(out_case + ".afreq", sep=r"\s+")[["ID", "ALT_FREQS"]].values)
    # 读取对照组等位基因频率字典
    ctrl_aaf_dict = dict(pd.read_csv(out_ctrl + ".afreq", sep=r"\s+")[["ID", "ALT_FREQS"]].values)

    # 读取病例组和对照组HWE检验p值字典
    hwe_case_df = pd.read_csv(out_case + ".hardy", sep=r"\s+")
    hwe_ctrl_df = pd.read_csv(out_ctrl + ".hardy", sep=r"\s+")
    case_hwe_dict = dict(hwe_case_df[["ID", "P"]].values)
    ctrl_hwe_dict = dict(hwe_ctrl_df[["ID", "P"]].values)

    # Step 5: parallel chunk processing
    output_file = output_prefix + ".variant_qc_summary.tsv"
    if verbose:
        test_chunk = pd.read_csv(out_all + ".afreq", sep=r"\s+", nrows=5)
        print("[DEBUG] .afreq 字段名:", list(test_chunk.columns))


    chunk_files = []
    reader = pd.read_csv(out_all + ".afreq", sep=r"\s+", chunksize=100000)
    with concurrent.futures.ProcessPoolExecutor(max_workers=4, initializer=init_globals_for_chunk,
                                                initargs=(vmiss_dict, case_aaf_dict, ctrl_aaf_dict, case_hwe_dict, ctrl_hwe_dict)) as executor:
        futures = []
        for i, chunk in enumerate(reader):
            if verbose:
                print(f"[INFO] 正在提交 chunk {i + 1} 任务...")
            futures.append(
                executor.submit(
                    process_chunk,
                    chunk, i, tmpdir
                )
            )
        for i, future in enumerate(futures):
            chunk_file = future.result()
            if verbose:
                print(f"[INFO] chunk {i + 1} 处理完成，结果文件: {chunk_file}")
            chunk_files.append(chunk_file)

    # merge chunk files
    with open(output_file, "w", newline="") as fout:
        writer = csv.writer(fout, delimiter="\t")
        # 写入表头
        writer.writerow([
            "VARIANT_ID", "MAF", "VMISS", "CASE_AAF",
            "CTRL_AAF", "CTRL_MAF", "CASE_HWE", "CTRL_HWE"
        ])
        # 合并所有chunk的结果文件
        for chunk_file in chunk_files:
            with open(chunk_file, "r") as fin:
                for line in fin:
                    fout.write(line)

    if verbose:
        print(f"[INFO] 输出完成: {output_file}")
        print(f"[INFO] 结果文件路径（可用于后续加载）: {output_file}")
    return output_file


def extract_maf0_or_vmiss1_variants_streaming(input_file: str,
                                               base_prefix: str = "maf0_or_vmiss1_variants",
                                               chunksize: int = 100000) -> tuple[str, str]:
    """
    从 QC 结果文件中提取 MAF=0 或 VMISS=1 及 MAF/VMISS/CASE_AAF/CTRL_AAF NA 的低质量变异，采用分块流式处理，避免内存溢出。

    参数:
        input_file (str): 变异 QC 汇总结果文件路径（TSV 格式）
        base_prefix (str): 输出文件前缀名（默认: maf0_or_vmiss1_variants）
        chunksize (int): 每块读取的行数

    返回:
        Tuple[str, str]:
            - 含有标志位的TSV路径（包含VARIANT_ID, MAF0_FLAG, VMISS1_FLAG, MAF_NA_FLAG, VMISS_NA_FLAG, CASE_AAF_NA_FLAG, CTRL_AAF_NA_FLAG）
            - 仅含VARIANT_ID列的TSV路径（用于plink2 --exclude）
    """
    input_dir = os.path.dirname(os.path.abspath(input_file))
    full_info_path = os.path.join(input_dir, f"{base_prefix}.with_flags.tsv")
    variant_only_path = os.path.join(input_dir, f"{base_prefix}.variant_ids.tsv")

    with open(full_info_path, 'w') as f_full, open(variant_only_path, 'w') as f_ids:
        # 写入表头，新增 MAF_NA_FLAG、VMISS_NA_FLAG、CASE_AAF_NA_FLAG、CTRL_AAF_NA_FLAG
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
                    # 写入含标志的文件
                    f_full.write(f"{row.VARIANT_ID}\t{is_maf0}\t{is_vmiss1}\t{is_maf_na}\t{is_vmiss_na}\t{is_case_aaf_na}\t{is_ctrl_aaf_na}\n")
                    # 写入仅含变异ID的文件，方便plink2排除
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
    使用 plink2 删除指定的低质量变异（来自 ID 列表），输出新的 Plink 二进制文件。

    参数:
        bed_prefix (str): 输入 .bed/.bim/.fam 文件前缀
        output_prefix (str): 输出文件前缀
        exclude_variants_file (str): 要排除的变异列表（每行一个变异ID，无表头）
        plink2_path (str): plink2 执行路径
        threads (int): 并行线程数

    返回:
        str: 新的 Plink 数据前缀（无扩展名）
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