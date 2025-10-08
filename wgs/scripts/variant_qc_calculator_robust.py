"""
变异质控计算工具集（variant_qc_calculator_robust.py）

【工具概述】
本工具是在 variant_qc_flags 和 variant_qc_calculator 基础上的增强版补充脚本，专门用于在大规模基因型数据（PLINK 二进制格式：.bed/.bim/.fam）上计算变异层面的质量控制（QC）指标，自动汇总并支持后续一键排除低质量变异，生成新的 PLINK 数据集。相比基础版本，本脚本增加了分组深度分析（30X/15X）、可视化功能和更robust的错误处理机制。流程以 plink2 为核心，辅以 pandas 的分块（streaming）与多进程处理，兼顾可复现性与性能。

【与前版本的关系】
- variant_qc_flags: 提供基础的变异标记和过滤功能
- variant_qc_calculator: 提供核心的QC指标计算功能  
- variant_qc_calculator_robust: 本脚本，在前两者基础上增强，提供更全面的分组分析和可视化功能

【适用场景】
- GWAS/候选位点分析前的标准变异级 QC 统计与筛除；
- 快速获取群体与分组（病例/对照、30X/15X深度组）的 AAF/MAF、缺失率（VMISS）、HWE p 值；
- 生成可追溯的 QC 汇总表与可直接用于 --exclude 的变异列表；
- 提供基于MAF分层和深度分组的缺失率可视化分析；
- 作为 variant_qc_flags 和 variant_qc_calculator 的增强版本，支持更复杂的分组分析需求。

【输入与前提】
- 输入：PLINK 二进制文件前缀（bed_prefix），要求 .bed/.bim/.fam 三件套齐全；
- fam 文件的 PHENO 列使用 1=对照、2=病例 的常规编码（缺失或其他编码将导致分组统计为空或异常）；
- 需要可执行的 plink2（二进制路径可通过 plink2_path 指定）。

【核心功能增强】
相比基础版本（variant_qc_flags 和 variant_qc_calculator），本robust版本新增：
1. 深度分组分析：支持基于样本信息表中的目标深度（30X/15X）进行分组缺失率统计
2. 可视化功能：提供基于MAF分层的缺失率散点图和KDE密度图
3. 阈值配置：支持从JSON配置文件读取不同MAF类别的缺失率阈值
4. 更robust的错误处理：增强的分块处理和内存管理
5. 四象限统计：在可视化中提供详细的通过/失败变异计数统计

【核心输出】
1) 变异 QC 汇总表：`{output_prefix}.variant_qc_summary.tsv`
   - 列说明：
     - VARIANT_ID：PLINK 变异 ID（与 .bim 中 ID 一致）
     - MAF：基于全体样本 ALT_FREQS 计算的 minor allele frequency（min(p, 1-p)）
     - VMISS：全体样本缺失率（F_MISS）
     - CASE_VMISS：病例组缺失率（F_MISS）
     - CTRL_VMISS：对照组缺失率（F_MISS）
     - 30X_VMISS：30× 目标深度组缺失率（F_MISS）
     - 15X_VMISS：15× 目标深度组缺失率（F_MISS）
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

4) 可视化输出（由 `plot_vmiss_scatter_by_maf_category` 生成，robust版本新增功能）
   - `{output_prefix}.vmiss.{mode}.png`：基于MAF分层的缺失率散点图，包含KDE密度图和四象限统计
   - `vmiss_pass_variants.tsv`：通过所有阈值的变异ID列表

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
1. 基础QC计算（与基础版本兼容）：
   ```python
   out_tsv = run_plink2_variant_qc(
       bed_prefix="/path/to/dataset",
       tmpdir="/tmp/variant_qc",
       plink2_path="/home/b/b37974/plink2",
       threads=16,
       output_prefix="cteph_agp3k",
       verbose=True,
       info_path="/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/info/cteph_agp3k_jhrpv4.rev1.xlsx",
       sample_col="ID",
       target_dp_col="Target DP (JHRPv4)"
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
4. 可视化分析（robust版本新增功能）：
   ```python
   pass_variants_file = plot_vmiss_scatter_by_maf_category(
       variant_qc_summary=out_tsv,
       vmiss_json_path="vmiss_thresholds.json",
       mode="dp",  # 或 "case_ctrl"
       output_prefix="cteph_agp3k"
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
- 运行环境：Python 3.9+（pandas、numpy、matplotlib、seaborn）；系统可执行：plink2
- 主要库：pandas, numpy, csv, subprocess, tempfile, concurrent.futures, uuid, typing, matplotlib, seaborn
- 相比基础版本增加：matplotlib（可视化）、seaborn（KDE密度图）

【可复现性与记录】
- 建议将 `plink2` 版本与调用参数记录在项目 README 或运行日志中；
- 输出文件名包含 `output_prefix`，便于与批量任务的样本队列对齐。

【版权与维护】
- 作者：ZHAO TIE

【更新记录】
- 2025-08-12：补充用途说明、字段定义、分块与并行策略、常见阈值参考与使用示例；明确中间文件与错误处理注意事项。
- 2025-08-13：在 `extract_maf0_or_vmiss1_variants_streaming` 中新增对 MAF 和 VMISS NA 值的检测与标记，方便识别缺失数据。
- 2025-10-06：创建robust版本，作为variant_qc_flags和variant_qc_calculator的补充脚本，增加深度分组分析、可视化功能和更强的错误处理机制。
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
import matplotlib.pyplot as plt

plt.style.use('default')

_global_lookup = {}

def init_globals_for_chunk(vmiss_dict, case_aaf_dict, ctrl_aaf_dict, case_hwe_dict, ctrl_hwe_dict,
                           case_vmiss_dict, ctrl_vmiss_dict, vmiss_30x_dict, vmiss_15x_dict):
    global _global_lookup
    _global_lookup['vmiss_dict'] = vmiss_dict
    _global_lookup['case_aaf_dict'] = case_aaf_dict
    _global_lookup['ctrl_aaf_dict'] = ctrl_aaf_dict
    _global_lookup['case_hwe_dict'] = case_hwe_dict
    _global_lookup['ctrl_hwe_dict'] = ctrl_hwe_dict
    _global_lookup['case_vmiss_dict'] = case_vmiss_dict
    _global_lookup['ctrl_vmiss_dict'] = ctrl_vmiss_dict
    _global_lookup['vmiss_30x_dict'] = vmiss_30x_dict
    _global_lookup['vmiss_15x_dict'] = vmiss_15x_dict

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
    case_vmiss_dict = _global_lookup.get('case_vmiss_dict', {})
    ctrl_vmiss_dict = _global_lookup.get('ctrl_vmiss_dict', {})
    vmiss_30x_dict = _global_lookup.get('vmiss_30x_dict', {})
    vmiss_15x_dict = _global_lookup.get('vmiss_15x_dict', {})

    tmp_output = os.path.join(tmpdir, f"chunk_{idx}_{uuid.uuid4().hex}.tsv")
    with open(tmp_output, "w", newline="") as fout:
        writer = csv.writer(fout, delimiter="\t")
        for _, row in chunk.iterrows():
            vid = row["ID"]  # 变异ID
            aaf = row["ALT_FREQS"]  # 等位基因频率
            maf = min(aaf, 1 - aaf) if pd.notnull(aaf) else float("nan")  # 计算MAF
            vmiss = vmiss_dict.get(vid, float("nan"))  # 缺失率VMISS（全体）
            case_vmiss = case_vmiss_dict.get(vid, float("nan"))
            ctrl_vmiss = ctrl_vmiss_dict.get(vid, float("nan"))
            vmiss_30x = vmiss_30x_dict.get(vid, float("nan"))
            vmiss_15x = vmiss_15x_dict.get(vid, float("nan"))
            case_aaf = case_aaf_dict.get(vid, float("nan"))  # 病例组等位基因频率
            ctrl_aaf = ctrl_aaf_dict.get(vid, float("nan"))  # 对照组等位基因频率
            ctrl_maf = min(ctrl_aaf, 1 - ctrl_aaf) if pd.notnull(ctrl_aaf) else float("nan")  # 对照组MAF
            case_hwe = case_hwe_dict.get(vid, float("nan"))  # 病例组HWE p值
            ctrl_hwe = ctrl_hwe_dict.get(vid, float("nan"))  # 对照组HWE p值
            # 写入列顺序: 变异ID, MAF, VMISS(全体), CASE_VMISS, CTRL_VMISS, 30X_VMISS, 15X_VMISS,
            # CASE_AAF, CTRL_AAF, CTRL_MAF, CASE_HWE, CTRL_HWE
            writer.writerow([vid, maf, vmiss, case_vmiss, ctrl_vmiss, vmiss_30x, vmiss_15x,
                             case_aaf, ctrl_aaf, ctrl_maf, case_hwe, ctrl_hwe])
    return tmp_output

def run_plink2_variant_qc(
    bed_prefix: str,
    tmpdir: str = "/tmp/variant_qc",
    plink2_path: str = "/home/b/b37974/plink2",
    threads: int = 8,
    output_prefix: str = "cteph_agp3k",
    verbose: bool = True,
    info_path: str = "/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/info/cteph_agp3k_jhrpv4.rev1.xlsx",
    sample_col: str = "ID",
    target_dp_col: str = "Target DP (JHRPv4)"
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
        info_path (str): Excel 格式的样本信息文件路径，含目标深度等信息
        sample_col (str): Excel 文件中样本ID列名
        target_dp_col (str): Excel 文件中目标深度列名

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

    # Step 2.5: 通过 meta 表（info_path）构建 30x/15x 分组 IID 列表，并计算各自的缺失率
    try:
        info_df = pd.read_excel(info_path)
    except Exception as e:
        raise RuntimeError(f"无法读取 info_path: {info_path}. 原因: {e}")

    if sample_col not in info_df.columns or target_dp_col not in info_df.columns:
        raise ValueError(f"info_path 中缺少必要列: {sample_col} 或 {target_dp_col}")

    # 规范化目标深度列，以兼容 '30x'/'30X'/' 30x ' 等写法
    norm_dp = info_df[target_dp_col].astype(str).str.strip().str.lower()
    info_df = info_df.assign(_norm_dp=norm_dp)

    # 将 Excel 的样本 ID 与 fam 的 IID 对齐，获取 FID/IID 成对列表
    fam_cols = ["FID", "IID", "PHENO"]
    fam_sub = fam_df[fam_cols]
    merged = fam_sub.merge(info_df[[sample_col, "_norm_dp"]], left_on="IID", right_on=sample_col, how="inner")

    dp30 = merged[merged["_norm_dp"].isin(["30x", "30"])][["FID", "IID"]]
    dp15 = merged[merged["_norm_dp"].isin(["15x", "15"])][["FID", "IID"]]

    dp30_iid_path = os.path.join(tmpdir, "dp30x_iids.txt")
    dp15_iid_path = os.path.join(tmpdir, "dp15x_iids.txt")
    dp30.to_csv(dp30_iid_path, sep="\t", index=False, header=False)
    dp15.to_csv(dp15_iid_path, sep="\t", index=False, header=False)

    out_dp30 = os.path.join(tmpdir, "dp30x")
    out_dp15 = os.path.join(tmpdir, "dp15x")
    # 只需缺失率即可
    run_cmd([
        plink2_path, "--bfile", bed_prefix, "--keep", dp30_iid_path,
        "--threads", str(threads), "--missing", "--out", out_dp30
    ], "30X VMISS")
    run_cmd([
        plink2_path, "--bfile", bed_prefix, "--keep", dp15_iid_path,
        "--threads", str(threads), "--missing", "--out", out_dp15
    ], "15X VMISS")

    # Step 3: 分组 freq + HWE + missingness
    out_case = os.path.join(tmpdir, "case")
    out_ctrl = os.path.join(tmpdir, "ctrl")
    # 计算病例组的等位基因频率、HWE检验与缺失率
    run_cmd([
        plink2_path, "--bfile", bed_prefix, "--keep", case_iid_path,
        "--threads", str(threads), "--freq", "--hardy", "--missing", "--out", out_case
    ], "Case AAF + HWE + VMISS")

    # 计算对照组的等位基因频率、HWE检验与缺失率
    run_cmd([
        plink2_path, "--bfile", bed_prefix, "--keep", ctrl_iid_path,
        "--threads", str(threads), "--freq", "--hardy", "--missing", "--out", out_ctrl
    ], "Control AAF + HWE + VMISS")

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

    # 读取分组缺失率（case/ctrl/30x/15x）
    case_vmiss_dict = dict(pd.read_csv(out_case + ".vmiss", sep=r"\s+")[ ["ID", "F_MISS"] ].values)
    ctrl_vmiss_dict = dict(pd.read_csv(out_ctrl + ".vmiss", sep=r"\s+")[ ["ID", "F_MISS"] ].values)
    vmiss_30x_dict = dict(pd.read_csv(out_dp30 + ".vmiss", sep=r"\s+")[ ["ID", "F_MISS"] ].values) if os.path.exists(out_dp30 + ".vmiss") else {}
    vmiss_15x_dict = dict(pd.read_csv(out_dp15 + ".vmiss", sep=r"\s+")[ ["ID", "F_MISS"] ].values) if os.path.exists(out_dp15 + ".vmiss") else {}

    # Step 5: parallel chunk processing
    output_file = output_prefix + ".variant_qc_summary.tsv"
    if verbose:
        test_chunk = pd.read_csv(out_all + ".afreq", sep=r"\s+", nrows=5)
        print("[DEBUG] .afreq 字段名:", list(test_chunk.columns))


    chunk_files = []
    reader = pd.read_csv(out_all + ".afreq", sep=r"\s+", chunksize=100000) # 每10万行一个chunk
    with concurrent.futures.ProcessPoolExecutor(
        max_workers=4,
        initializer=init_globals_for_chunk,
        initargs=(vmiss_dict, case_aaf_dict, ctrl_aaf_dict, case_hwe_dict, ctrl_hwe_dict,
                  case_vmiss_dict, ctrl_vmiss_dict, vmiss_30x_dict, vmiss_15x_dict)
    ) as executor:
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
            "VARIANT_ID", "MAF", "VMISS", "CASE_VMISS", "CTRL_VMISS", "30X_VMISS", "15X_VMISS",
            "CASE_AAF", "CTRL_AAF", "CTRL_MAF", "CASE_HWE", "CTRL_HWE"
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



def plot_vmiss_scatter_by_maf_category(
    variant_qc_summary: str,
    vmiss_json_path: str,
    mode: str = "dp",  # "dp" or "case_ctrl"
    variant_id_col: str = "VARIANT_ID",
    output_tsv: str = "vmiss_pass_variants.tsv",
    output_prefix: str = "cteph_agp3k",
    plot_style: str = "hex",          # "hex" | "hist2d" | "scatter" | "kde2d" | "density"
    gridsize: int = 75,               # hexbin 网格密度（越大越细）
    hist_bins: int = 75,              # hist2d 的分箱数
    density_norm: str = "log",        # "log" or "linear" 颜色归一化
    bw_adjust: float = 1.0,           # KDE 平滑程度（越大越平滑）
    density_thresh: float = 0.02,     # density风格下的最小显示阈值（0~1，相对密度）
    overlay_points: bool = False      # 是否在密度图上叠加半透明小散点
) -> str:
    """
    【已重载为 VMISS 版本】
    根据 CTRL_MAF 分三类（Rare/Low Frequency/Common），并基于 VMISS 指标绘制二维散点：
      - 当 mode='dp' 时：x=30X_VMISS，y=15X_VMISS
      - 当 mode='case_ctrl' 时：x=CTRL_VMISS，y=CASE_VMISS
    阈值来源于 vmiss_json_path（见 vmiss.json 结构说明），用于可视化阈值线并进行通过/未通过统计与导出。

    参数
    ----
    variant_qc_summary : str
        变体质控汇总 TSV 路径（需包含列：VARIANT_ID, CTRL_MAF, CASE_VMISS, CTRL_VMISS, 30X_VMISS, 15X_VMISS）
    vmiss_json_path : str
        vmiss.json 文件路径，包含不同模式下各 MAF 类别的 VMISS 阈值。
    mode : {"dp","case_ctrl"}
        选择阈值与坐标轴来源；"dp" 使用 30X_VMISS vs 15X_VMISS，"case_ctrl" 使用 CTRL_VMISS vs CASE_VMISS。
    variant_id_col : str
        变体 ID 列名，默认 "VARIANT_ID"。
    output_tsv : str
        输出通过阈值（双向都通过）的变体 ID 列表文件路径。
    output_prefix : str
        输出图文件前缀；将保存为 `{output_prefix}.vmiss.{mode}.png`。
    plot_style : {"hex","hist2d","scatter","kde2d","density"}
        控制主面板点的呈现方式，默认 "hex"（六边形密度），可有效缓解过密覆盖问题。
    gridsize : int
        当 plot_style="hex" 时的六边形网格密度。
    hist_bins : int
        当 plot_style="hist2d" 时的二维直方图分箱数。
    density_norm : {"log","linear"}
        密度着色的归一化方式；数目跨度大时推荐 "log"。
    bw_adjust : float
        仅在 "kde2d"/"density" 风格下生效，控制核密度平滑程度（>1 更平滑）。
    density_thresh : float
        仅在 "density" 风格下生效，控制最小等高线显示阈值（避免过稀区域完全不可见）。
    overlay_points : bool
        为 True 时，在密度图上叠加少量半透明散点以增强边缘感知。

    返回
    ----
    str
        通过阈值变体 ID 列表文件路径。
    """
    import json
    import gc
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes
    import seaborn as sns  # 仅用于KDE，可选；若未安装，请注释相应KDE代码
    from matplotlib.colors import LogNorm, LinearSegmentedColormap, to_rgba

    # ---------- 基本检查 ----------
    if mode not in {"dp", "case_ctrl"}:
        raise ValueError("mode 必须为 'dp' 或 'case_ctrl'")

    # 读取阈值 JSON
    with open(vmiss_json_path, "r") as f:
        cfg = json.load(f)
    if mode not in cfg:
        raise KeyError(f"在 {vmiss_json_path} 中未找到模式 '{mode}' 的配置")

    # 需要的列
    base_cols = [variant_id_col, "CTRL_MAF"]
    if mode == "dp":
        x_col, y_col = "30X_VMISS", "15X_VMISS"
    else:
        x_col, y_col = "CTRL_VMISS", "CASE_VMISS"
    usecols = base_cols + [x_col, y_col]

    # --------- 将 CTRL_MAF 分类 ---------
    def classify_maf(maf):
        if pd.isna(maf) or maf < 0:
            return None
        if maf < 0.01:
            return "Rare Variant (<0.01)"
        elif 0.01 <= maf <= 0.05:
            return "Low Frequency Variant (0.01~0.05)"
        elif maf > 0.05:
            return "Common Variant (>0.05)"
        return None

    # 分块读取
    data_chunks = []
    reader = pd.read_csv(variant_qc_summary, sep="\t", usecols=usecols, chunksize=500000)
    for chunk in reader:
        chunk["Category"] = chunk["CTRL_MAF"].apply(classify_maf)
        chunk = chunk[chunk["Category"].notna()]
        data_chunks.append(chunk)
        del chunk
        gc.collect()
    if len(data_chunks) == 0:
        raise RuntimeError("未读取到任何可用数据，请检查输入列是否存在/非空")

    df = pd.concat(data_chunks, ignore_index=True)

    # --------- 准备绘图 ---------
    categories = [
        "Rare Variant (<0.01)",
        "Low Frequency Variant (0.01~0.05)",
        "Common Variant (>0.05)"
    ]
    colors = ["#1f77b4", "#2ca02c", "#ff7f0e"]

    fig = plt.figure(figsize=(14, 14))
    outer_gs = gridspec.GridSpec(2, 2, width_ratios=[1, 1], height_ratios=[1, 1])
    subplots = [outer_gs[0, 0], outer_gs[0, 1], outer_gs[1, 0]]

    summary_rows = []
    pass_mask_global = pd.Series(False, index=df.index)

    def _make_density_cmap(base_hex: str) -> LinearSegmentedColormap:
        """
        基于分类颜色生成密度图渐变色（由浅至深、带透明度），增强对比但保留低密度可见。
        """
        # 提高整体起始不透明度，同时拉大分段差距
        c0 = to_rgba(base_hex, 0.15)   # 最浅端：明确可见
        c1 = to_rgba(base_hex, 0.40)   # 浅
        c2 = to_rgba(base_hex, 0.80)   # 中深
        c3 = to_rgba(base_hex, 1.00)   # 最深
        return LinearSegmentedColormap.from_list("density_" + base_hex, [c0, c1, c2, c3])


    def plot_one(ax_main, ax_top, ax_right, sub_df, label, color, thr_x, thr_y):
        # x/y 即 VMISS（0~1），阈值解释为 "≤ 阈值 通过"
        x_raw = pd.to_numeric(sub_df[x_col], errors="coerce")
        y_raw = pd.to_numeric(sub_df[y_col], errors="coerce")

        # 通过判定（基于原始值）
        pass_mask = (x_raw <= thr_x) & (y_raw <= thr_y)
        # 在全局索引上记录通过
        nonlocal pass_mask_global
        pass_mask_global.loc[sub_df.index] = pass_mask | pass_mask_global.loc[sub_df.index]

        # 四象限计数（以阈值为界，将阈值位置画线）
        q_pass = pass_mask.sum()
        q_fail_x = ((x_raw > thr_x) & (y_raw <= thr_y)).sum()
        q_fail_y = ((x_raw <= thr_x) & (y_raw > thr_y)).sum()
        q_fail_both = ((x_raw > thr_x) & (y_raw > thr_y)).sum()

        # ---- 主面板可视化（根据 plot_style 切换） ----
        x = pd.to_numeric(x_raw, errors="coerce")
        y = pd.to_numeric(y_raw, errors="coerce")
        valid = x.notna() & y.notna()
        xv = x[valid]
        yv = y[valid]

        # 针对当前类别颜色生成专属密度 colormap
        cmap_local = _make_density_cmap(color)
        # 颜色归一化
        norm_obj = LogNorm(vmin=1) if density_norm == "log" else None

        if plot_style == "hex":
            hb = ax_main.hexbin(
                xv, yv,
                extent=[0, 1, 0, 1],
                gridsize=gridsize,
                mincnt=1,
                linewidths=0,
                norm=norm_obj,
                cmap=cmap_local
            )
        elif plot_style == "hist2d":
            h = ax_main.hist2d(
                xv, yv,
                bins=hist_bins,
                range=[[0, 1], [0, 1]],
                norm=norm_obj,
                cmap=cmap_local
            )
        elif plot_style == "density":
            # Tableau-like 密度风格：连续填色、无等高线、较高 levels、阈值裁剪
            sns.kdeplot(
                x=xv, y=yv,
                ax=ax_main,
                fill=True,
                thresh=density_thresh,    # 避免极稀薄区域完全不可见
                levels=100,               # 更平滑的连续密度
                bw_adjust=bw_adjust,
                cmap=cmap_local,
                linewidths=0
            )
            if overlay_points:
                ax_main.scatter(xv, yv, alpha=0.05, c=color, s=2, linewidths=0)
        elif plot_style == "kde2d":
            # 传统二维 KDE：带等高线
            sns.kdeplot(
                x=xv, y=yv,
                ax=ax_main,
                fill=True,
                thresh=0,
                levels=30,
                bw_adjust=bw_adjust,
                cmap=cmap_local,
                linewidths=0.8
            )
        else:  # "scatter"
            ax_main.scatter(xv, yv, alpha=0.1, c=color, s=5)

        # 阈值线与坐标
        ax_main.axvline(x=thr_x, color='red', linestyle='--')
        ax_main.axhline(y=thr_y, color='blue', linestyle='--')
        ax_main.set_xlim(0, 1)
        ax_main.set_ylim(0, 1)
        xlabel = f"{x_col}"
        ylabel = f"{y_col}"
        ax_main.set_xlabel(xlabel)
        ax_main.set_ylabel(ylabel)
        ax_main.grid(True, linestyle=":", linewidth=0.5)

        # ---- 顶部 KDE（x） ----
        x_plot = x[valid].dropna()
        if x_plot.nunique() >= 2:
            sns.kdeplot(x=x_plot, ax=ax_top, fill=True, color='gray', linewidth=1.5, cut=0)
            ax_top.set_xlim(0, 1)
            ax_top.axvline(x=thr_x, color='red', linestyle='--')
        else:
            ax_top.text(0.5, 0.5, 'KDE skipped', ha='center', va='center', transform=ax_top.transAxes, fontsize=8)
        ax_top.set_xlabel('')
        ax_top.set_ylabel('')
        ax_top.tick_params(axis='both', which='both', bottom=False, top=False, left=False, right=False,
                           labelbottom=False, labelleft=False)
        ax_top.grid(False)

        # ---- 右侧 KDE（y） ----
        y_plot = y[valid].dropna()
        if y_plot.nunique() >= 2:
            sns.kdeplot(y=y_plot, ax=ax_right, fill=True, color='gray', linewidth=1.5, cut=0)
            ax_right.set_ylim(0, 1)
            ax_right.axhline(y=thr_y, color='blue', linestyle='--')
        else:
            ax_right.text(0.5, 0.5, 'KDE skipped', ha='center', va='center', transform=ax_right.transAxes, fontsize=8)
        ax_right.set_xlabel('')
        ax_right.set_ylabel('')
        ax_right.tick_params(axis='both', which='both', bottom=False, top=False, left=False, right=False,
                             labelbottom=False, labelleft=False)
        ax_right.grid(False)

        # ---- 角落插图：计数摘要 ----
        top_ylim = ax_top.get_ylim()
        right_xlim = ax_right.get_xlim()
        inset_ax = inset_axes(ax_main, width="25%", height="25%", loc='upper right',
                              bbox_to_anchor=(0.28, 0.28, 1, 1), bbox_transform=ax_main.transAxes, borderpad=0)
        inset_ax.set_xlim(right_xlim)
        inset_ax.set_ylim(top_ylim)
        for spine in inset_ax.spines.values():
            spine.set_visible(True)
        inset_ax.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)
        inset_ax.grid(False)
        xmid = (right_xlim[0] + right_xlim[1]) / 2
        ymid = (top_ylim[0] + top_ylim[1]) / 2
        # 计算四块区域的中心点
        x_left = (right_xlim[0] + xmid) / 2
        x_right = (xmid + right_xlim[1]) / 2
        y_bottom = (top_ylim[0] + ymid) / 2
        y_top = (ymid + top_ylim[1]) / 2
        # 以竖线/横线分四块：通过（左下），仅X不通过（右下），仅Y不通过（左上），双不通过（右上）
        inset_ax.axvline(x=xmid, linestyle='--', color='red', linewidth=1.5)
        inset_ax.axhline(y=ymid, linestyle='--', color='blue', linewidth=1.5)
        # 左上：仅Y不通过
        inset_ax.text(x_left, y_top, f'{q_fail_y:,}', ha='center', va='center', fontsize=7, fontweight='bold')
        # 右上：X和Y都不通过
        inset_ax.text(x_right, y_top, f'{q_fail_both:,}', ha='center', va='center', fontsize=7, fontweight='bold')
        # 左下：双向通过
        inset_ax.text(x_left, y_bottom, f'{q_pass:,}', ha='center', va='center', fontsize=7, fontweight='bold')
        # 右下：仅X不通过
        inset_ax.text(x_right, y_bottom, f'{q_fail_x:,}', ha='center', va='center', fontsize=7, fontweight='bold')

    # 遍历三类分别绘图与统计
    for subplot_spec, category, color in zip(subplots, categories, colors):
        sub_df = df[df["Category"] == category]
        # 读取对应类别的阈值
        cat_cfg = cfg[mode].get(category, {})
        if mode == "dp":
            thr_x = float(cat_cfg.get("30X_VMISS"))
            thr_y = float(cat_cfg.get("15X_VMISS"))
        else:  # case_ctrl
            thr_x = float(cat_cfg.get("CTRL_VMISS"))
            thr_y = float(cat_cfg.get("CASE_VMISS"))

        # 内嵌 2x2（主图 + 顶部KDE + 右侧KDE）
        inner_gs = gridspec.GridSpecFromSubplotSpec(2, 2, subplot_spec=subplot_spec,
                                                    width_ratios=[6, 1.5], height_ratios=[1.5, 6],
                                                    wspace=0.05, hspace=0.05)
        ax_main = plt.Subplot(fig, inner_gs[1, 0])  # type: ignore
        ax_top = plt.Subplot(fig, inner_gs[0, 0], sharex=ax_main)  # type: ignore
        ax_right = plt.Subplot(fig, inner_gs[1, 1], sharey=ax_main)  # type: ignore
        fig.add_subplot(ax_main)
        fig.add_subplot(ax_top)
        fig.add_subplot(ax_right)

        # 汇总表（右下角表格使用）
        summary_rows.append({
            "": category.replace(" (", "\n("),
            "raw_category": category,
            f"{x_col} Threshold": f"{thr_x:.3f}",
            f"{y_col} Threshold": f"{thr_y:.3f}",
        })

        plot_one(ax_main, ax_top, ax_right, sub_df, category, color, thr_x, thr_y)

    # 右下角信息表
    ax_legend = fig.add_subplot(outer_gs[1, 1])
    ax_legend.axis('off')
    summary_df = pd.DataFrame(summary_rows)
    summary_df_display = summary_df.drop(columns=["raw_category"])
    table = ax_legend.table(cellText=summary_df_display.values,  # type: ignore
                            colLabels=summary_df_display.columns,  # type: ignore
                            cellLoc='center',
                            colWidths=[0.35, 0.325, 0.325],
                            loc='center')
    table.scale(1.2, 1.6)
    table.auto_set_font_size(False)
    table.set_fontsize(9)
    # 表头样式
    for col_idx in range(len(summary_df_display.columns)):
        cell = table[(0, col_idx)]
        cell.set_facecolor('#f0f0f0')
        cell.set_text_props(color='black', weight='bold')
    # 行名着色
    for row_idx, raw_category in enumerate(summary_df["raw_category"]):
        color = colors[categories.index(raw_category)]
        facecolor_rgba = plt.matplotlib.colors.to_rgba(color, alpha=0.35)  # type: ignore
        cell = table[(row_idx + 1, 0)]
        cell.set_facecolor(facecolor_rgba)
        cell.set_text_props(color='black', weight='bold')

    title_suffix = "30X_VMISS vs 15X_VMISS" if mode == "dp" else "CTRL_VMISS vs CASE_VMISS"
    fig.suptitle(f"{title_suffix} by CTRL_MAF Category", fontsize=16)
    plt.tight_layout(rect=[0, 0, 1, 0.96])  # type: ignore
    out_png = f"{output_prefix}.vmiss.{mode}.png"
    plt.savefig(out_png, dpi=600)
    # plt.show()

    # --------- 导出通过阈值的变体 ID ---------
    if variant_id_col not in df.columns:
        print(f"[Warning] {variant_id_col} column not found; skipping export of passed variants.")
        return output_tsv

    pass_variants = df.loc[pass_mask_global, variant_id_col].dropna().unique()
    pd.Series(pass_variants).to_csv(output_tsv, sep="\t", index=False, header=False)
    return output_tsv