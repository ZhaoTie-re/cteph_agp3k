"""
模块名称：bbj_projection_tools
=================================

作者：ZHAO TIE
日期：2025-08-20

【概述】
本模块围绕“将研究数据投影到 BBJ（Biobank Japan）PCA 空间，并据此进行样本选择与子集导出”的完整流程，
提供从基因型准备、PCA 输入生成、投影打分、可视化到样本 keep 导出的实用函数集合。

【功能总览】
1) prepare_bbj_genotype(...)：
   - 从已通过 QC 的 BBJ PLINK 前缀出发，导出 VCF（按 MAF 过滤）→ 染色体重命名 → 规范化/拆分多等位 → 设置变体 ID →
     最终生成标准化的 BBJ PLINK 三件套（bed/bim/fam）。
   - 所有中间产物与每步日志写入本地持久化临时目录（默认保留，便于复查）。

2) prepare_bbj_pca_inputs(...)：
   - 基于 BBJ PLINK 前缀，统一前缀化 FID/IID，执行只保留 A/C/G/T、去除高 LD 区域与 LD 修剪，
     输出 `no_high_ld` 前缀与 `prune.in/out`；随后以修剪后的变体计算 PCA（approx + allele-wts），
     产出 `*.eigenvec/*.eigenval` 供下游使用。

3) plot_pca_pairwise_pdf(...)：
   - 使用 PCA 结果绘制 PC 成对散点图（PC1vsPC2、PC3vsPC4…），在轴标中标注方差解释比例；
     采用“保存 PNG 再嵌入 PDF”的策略以保证渲染一致性。

4) run_bbj_projection(...)：
   - 将“我的数据”与 BBJ 在相同修剪集合下对齐并合并，随后用 BBJ 的频率与等位信息执行投影打分（plink2 --score），
     返回 `*.sscore` 文件。

5) plot_projection_pairwise_pdf(...)：
   - 基于投影结果（PC*_AVG 列）绘制成对散点图；按 IID 规则分组为 BBJ / 病例 / 对照，分层绘制以便对比。

6) plot_projection_pc1_pc2_and_select(...)：
   - 绘制 PC1_AVG vs PC2_AVG 的“主图 + 边际 KDE”，并在主图上按矩形范围选择病例/对照样本；
     生成高清 PNG 与 `keep.txt`（两列：#FID、IID）。

7) export_subset_bed(...)：
   - 基于 `keep.txt` 使用 plink2 --keep 从任意 PLINK 前缀导出子集（bed/bim/fam）。

【依赖与环境】
- 外部工具：plink2、plink、bcftools、samtools（建议用于生成 fasta .fai）。
- Python 库：pandas、matplotlib、seaborn。
- 建议在 Conda/venv 环境中运行，确保外部可执行程序在 $PATH 或通过函数参数提供。

【日志与可追溯性】
- `_run_cmd` 会将每个子命令的 stdout/stderr 追加写入对应 step*.log；
- `_progress` 会输出带时间戳的进度信息，并可选同步到 progress 日志；
- 关键产物缺失会抛出带有日志路径提示的异常，便于定位故障步骤。

【输入/输出规范（关键约定）】
- 所有以 “*_prefix” 结尾的参数均为**不含扩展名**的前缀路径；
- `keep.txt` 至少包含两列（#FID、IID），允许带表头，函数会在内部规范化为两列无表头再传给 plink2；
- 产出文件的默认写入目录：
    * 标准化产物与 PDF/PNG：当前工作目录；
    * 中间文件与分步日志：模块内部创建的本地持久化 tmp 目录（可通过参数选择保留/清理）。

【典型用法（摘要）】
- 生成 BBJ 标准化基因型与 PCA：
    1) prepare_bbj_genotype(...)
    2) prepare_bbj_pca_inputs(...)
- 将研究数据投影到 BBJ 空间：
    3) run_bbj_projection(...)
    4) plot_projection_pc1_pc2_and_select(...)
    5) export_subset_bed(...)

本模块力求：命令与参数与用户既有流程保持一致；出现错误时，提供最短路径的定位信息；
同时以“可复现、可追溯”为原则保留中间产物（可选自动清理）。
"""
from __future__ import annotations

import os
import shutil
import subprocess
import tempfile
from pathlib import Path

from typing import Tuple
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages


def _run_cmd(cmd: list[str], log_path: Path) -> None:
    """
    以追加方式将 stdout/stderr 写入 log_path，并在命令失败时抛出异常。
    """
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("a") as lf:
        lf.write("\n$ " + " ".join(cmd) + "\n")
        lf.flush()
        proc = subprocess.run(cmd, stdout=lf, stderr=lf, text=True)
        if proc.returncode != 0:
            raise RuntimeError(
                f"Command failed (exit {proc.returncode}). See log: {log_path}"
            )


# 轻量级进度记录器
def _progress(message: str, log_file: Path | None = None) -> None:
    """
    轻量级进度记录器：打印到控制台，并可选将同样信息附加写入一个 progress 日志文件。
    便于在失败时快速定位到出错步骤。
    """
    from datetime import datetime
    ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    line = f"[BBJ-PREP {ts}] {message}"
    print(line)
    if log_file is not None:
        log_file.parent.mkdir(parents=True, exist_ok=True)
        with log_file.open("a") as f:
            f.write(line + "\n")

def prepare_bbj_genotype(
    bbj_bed_prefix: str,
    threads: int,
    rename_chr: str,
    fasta_ref: str,
    maf: float = 0.05,
    plink2_path: str = "/home/b/b37974/plink2",
    plink1_path: str = "/home/b/b37974/plink",
    output_prefix: str = "bbj.b38.auto.sqc.vqc.rechr.norm.setid",
    keep_tmp: bool = True,
) -> Path:
    """
    准备 BBJ 基因型：从已通过样本/位点 QC 的 PLINK 二进制前缀出发，
    导出 VCF（按 MAF 过滤），重命名染色体，规范化/拆分，设置变体 ID，
    最终仅输出目标 PLINK 三件套与 .log：
        bbj.b38.auto.sqc.vqc.rechr.norm.setid.{bed,bim,fam,log}

    其余中间文件全部写入本地持久化临时目录（`./bbj_prep_tmp/...`），默认**不**自动清理，便于故障复查。

    参数
    ----
    bbj_bed_prefix : str
        输入 BBJ 的 PLINK 前缀（路径不带扩展名），需存在 .bed/.bim/.fam。
    threads : int
        线程数（用于 plink/plink2 与 bcftools）。
    rename_chr : str
        染色体重命名映射文件（bcftools annotate --rename-chrs 的文本映射表）。
    fasta_ref : str
        参考基因组 fasta（需包含 .fai 索引；用于 bcftools norm --fasta-ref）。
    maf : float, default 0.05
        导出 VCF 前的 MAF 过滤阈值。
    plink2_path : str, default "/home/b/b37974/plink2"
        plink2 可执行程序路径（用于从 bfile 导出 VCF）。
    plink1_path : str, default "/home/b/b37974/plink"
        plink1 可执行程序路径（用于 VCF → PLINK bed 转换）。
    output_prefix : str, default "bbj.b38.auto.sqc.vqc.rechr.norm.setid"
        最终输出前缀（无扩展名），会在工作目录下生成对应 bed/bim/fam/log。
    keep_tmp : bool, default True
        是否保留中间文件的临时目录。True 时保留，False 时在流程结束后自动清理。

    调试提示
    -------
    - 每个步骤都会通过 `_progress` 输出带时间戳的进度信息（同时写入 progress 日志）。
    - 每个子命令的完整 stdout/stderr 会写入临时目录中的 step*.log 文件（由 `_run_cmd` 负责）。
    - 若抛出异常，异常消息中会包含失败子步骤的 log 文件路径，便于快速定位。

    返回
    ----
    Path
        最终输出前缀（无扩展名）的绝对路径（output_prefix）。
    """
    progress_log = Path("bbj_prep.progress.log").resolve()
    _progress("启动 BBJ 基因型准备流程", progress_log)

    # 参数校验
    if threads <= 0:
        raise ValueError(f"threads 必须为正整数，当前: {threads}")
    if not (0 < maf < 0.5):
        raise ValueError(f"maf 应在 (0, 0.5) 区间，当前: {maf}")

    # 依赖软件预检
    bcftools_path = shutil.which("bcftools")
    if bcftools_path is None:
        raise FileNotFoundError("未找到 bcftools，请确认其已安装并在 $PATH 中。")
    _progress(f"检测到 bcftools: {bcftools_path}", progress_log)

    # 前置校验
    plink2 = Path(plink2_path).expanduser().resolve()
    plink1 = Path(plink1_path).expanduser().resolve()
    if not plink2.exists():
        raise FileNotFoundError(f"plink2 不存在: {plink2}")
    if not plink1.exists():
        raise FileNotFoundError(f"plink 不存在: {plink1}")
    _progress(f"检测到 plink2: {plink2}", progress_log)
    _progress(f"检测到 plink : {plink1}", progress_log)

    bbj_prefix = Path(bbj_bed_prefix).expanduser().resolve()
    for ext in (".bed", ".bim", ".fam"):
        if not (bbj_prefix.parent / f"{bbj_prefix.name}{ext}").exists():
            raise FileNotFoundError(f"缺少输入文件: {bbj_prefix}{ext}")
    _progress(f"输入前缀就绪: {bbj_prefix}", progress_log)

    rename_chr_path = Path(rename_chr).expanduser().resolve()
    if not rename_chr_path.exists():
        raise FileNotFoundError(f"rename_chr 文件不存在: {rename_chr_path}")
    fasta_ref_path = Path(fasta_ref).expanduser().resolve()
    if not fasta_ref_path.exists():
        raise FileNotFoundError(f"参考序列不存在: {fasta_ref_path}")
    _progress(f"重命名表: {rename_chr_path}", progress_log)
    _progress(f"参考序列: {fasta_ref_path}", progress_log)
    # .fai 索引建议存在
    if not (fasta_ref_path.parent / (fasta_ref_path.name + ".fai")).exists():
        _progress(
            "警告：未找到参考序列的 .fai 索引，后续 `bcftools norm` 可能失败；请先运行 `samtools faidx`。",
            progress_log,
        )

    # 目标输出基名（写到当前工作目录）
    bbj_prepared_prefix = Path(output_prefix).resolve()
    _progress(f"最终输出前缀: {bbj_prepared_prefix}", progress_log)
    bed_out = Path(str(bbj_prepared_prefix) + ".bed")
    bim_out = Path(str(bbj_prepared_prefix) + ".bim")
    fam_out = Path(str(bbj_prepared_prefix) + ".fam")
    log_out = Path(str(bbj_prepared_prefix) + ".log")

    # 使用本地持久化的临时目录保存所有中间文件（便于排错与复查，不自动清理）
    from datetime import datetime
    tmp_root = Path("./bbj_prep_tmp").resolve()
    tmp_root.mkdir(parents=True, exist_ok=True)
    tmp = tmp_root / datetime.now().strftime("run_%Y%m%d_%H%M%S_%f")
    tmp.mkdir(parents=True, exist_ok=False)
    _progress(f"创建本地临时目录: {tmp}", progress_log)

    # 1) 导出 VCF（bgzip），按 MAF 过滤
    vcf_base = tmp / f"bbj.b38.auto.qc.maf{maf:g}"
    _progress(f"[Step1] 导出 VCF（MAF={maf:g}）→ {vcf_base}.vcf.gz", progress_log)
    cmd_export = [
        str(plink2),
        "--bfile", str(bbj_prefix),
        "--export", "vcf", "bgz",
        "--maf", str(maf),
        "--out", str(vcf_base),
        "--threads", str(threads),
    ]
    _run_cmd(cmd_export, tmp / "step1_export_vcf.log")

    vcf_gz = Path(str(vcf_base) + ".vcf.gz")
    if not vcf_gz.exists():
        raise FileNotFoundError(f"导出 VCF 失败，未找到文件: {vcf_gz}")
    _progress(f"[Step1] 导出完成: {vcf_gz}", progress_log)

    # 2) 建立 tabix 索引
    _run_cmd(["bcftools", "index", "-t", str(vcf_gz), "--threads", str(threads)], tmp / "step1_index_vcf.log")
    _progress(f"[Step1] 已建立索引: {vcf_gz}.tbi", progress_log)

    # 3) 重命名染色体
    rechr_vcf = tmp / "bbj.b38.auto.sqc.vqc.rechr.vcf.gz"
    _progress(f"[Step2] 染色体重命名 → {rechr_vcf}", progress_log)
    _run_cmd([
        "bcftools", "annotate",
        str(vcf_gz),
        "--rename-chrs", str(rename_chr_path),
        "-Oz", "-o", str(rechr_vcf),
        "--threads", str(threads),
    ], tmp / "step2_rename_chr.log")
    _run_cmd(["bcftools", "index", "-t", str(rechr_vcf), "--threads", str(threads)], tmp / "step2_index_rechr.log")

    # 4) 标准化与拆分多等位（-m-），并检查参考等位基因
    norm_vcf = tmp / "bbj.b38.auto.sqc.vqc.rechr.norm.vcf.gz"
    _progress(f"[Step3] 规范化与拆分（-m-，--check-ref s）→ {norm_vcf}", progress_log)
    _run_cmd([
        "bcftools", "norm",
        "--fasta-ref", str(fasta_ref_path),
        "-m-", "--check-ref", "s",
        str(rechr_vcf),
        "-Oz", "-o", str(norm_vcf),
        "--threads", str(threads),
    ], tmp / "step3_norm.log")
    _run_cmd(["bcftools", "index", "-t", str(norm_vcf), "--threads", str(threads)], tmp / "step3_index_norm.log")

    # 5) 设定变体 ID（CHROM:POS:REF:ALT）
    setid_vcf = tmp / "bbj.b38.auto.sqc.vqc.rechr.norm.setid.vcf.gz"
    _progress(f"[Step4] 设置变体 ID（%CHROM:%POS:%REF:%ALT）→ {setid_vcf}", progress_log)
    _run_cmd([
        "bcftools", "annotate",
        "--set-id", "%CHROM:%POS:%REF:%ALT",
        str(norm_vcf),
        "-Oz", "-o", str(setid_vcf),
        "--threads", str(threads),
    ], tmp / "step4_setid.log")
    _run_cmd(["bcftools", "index", "-t", str(setid_vcf), "--threads", str(threads)], tmp / "step4_index_setid.log")

    # 6) VCF → PLINK bed，输出到工作目录，仅保留最终四件套
    _progress(f"[Step5] VCF→PLINK bed/fam/bim → {bbj_prepared_prefix}.bed|.bim|.fam", progress_log)
    _run_cmd([
        str(plink1),
        "--vcf", str(setid_vcf),
        "--make-bed",
        "--keep-allele-order",
        "--double-id",
        "--out", str(bbj_prepared_prefix),
    ], tmp / "step5_vcf2bed.log")
    _progress("[Step5] 转换完成", progress_log)
    _progress(f"中间文件与分步日志已保存在: {tmp}", progress_log)

    _progress("校验最终产物是否存在", progress_log)
    # 产物校验
    for p in (bed_out, bim_out, fam_out, log_out):
        if not p.exists():
            raise FileNotFoundError(f"缺少目标输出文件: {p}")

    # 自动清理临时目录
    if not keep_tmp:
        shutil.rmtree(tmp, ignore_errors=True)
        _progress(f"已删除临时目录: {tmp}", progress_log)

    _progress("流程完成 ✓", progress_log)
    return bbj_prepared_prefix



def prepare_bbj_pca_inputs(
    bbj_bed_prefix: str,
    output_prefix: str = "bbj.b38.auto.sqc.vqc.norm",
    threads: int = 16,
    high_ld: str = "/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/info/high-LD-regions-hg38-GRCh38.txt",
    plink2_path: str = "/home/b/b37974/plink2",
    keep_tmp: bool = True,
) -> tuple[Path, Path, Path, Path, Path]:
    """
    基于 BBJ 的 PLINK 前缀执行以下流程：
      1) 将 FID/IID 前缀化为 "bbj_" 并生成新的中间前缀 `output_prefix`；
      2) 进行只保留 A/C/G/T、去除高 LD 区域、LD 修剪（indep-pairwise 50 5 0.2），
         输出 `f"{output_prefix}.no_high_ld"` 的 PLINK 三件套及 `prune.in/out`；
      3) 使用修剪后的变体计算 PCA（approx + allele-wts），输出到
         `f"{output_prefix}.no_high_ld.prune.pca"`（与调用示例一致）。
      4) 可选是否保留中间缓存目录（keep_tmp 参数控制）。

    返回
    ----
    tuple[Path, Path, Path, Path, Path]
        分别为：
        - `no_high_ld` 的 PLINK 前缀 Path（不带扩展名）；
        - `prune.in` 文件 Path；
        - `prune.out` 文件 Path。
        - PCA 的 eigenvec 文件 Path（*.prune.pca.eigenvec）；
        - PCA 的 eigenval 文件 Path（*.prune.pca.eigenval）。
    """
    progress_log = Path("bbj_prep.progress.log").resolve()
    _progress("启动 BBJ PCA 输入准备流程", progress_log)

    # 参数与依赖检查
    if threads <= 0:
        raise ValueError(f"threads 必须为正整数，当前: {threads}")
    plink2 = Path(plink2_path).expanduser().resolve()
    if not plink2.exists():
        raise FileNotFoundError(f"plink2 不存在: {plink2}")
    bbj_prefix = Path(bbj_bed_prefix).expanduser().resolve()
    for ext in (".bed", ".bim", ".fam"):
        if not (bbj_prefix.parent / f"{bbj_prefix.name}{ext}").exists():
            raise FileNotFoundError(f"缺少输入文件: {bbj_prefix}{ext}")
    high_ld_path = Path(high_ld).expanduser().resolve()
    if not high_ld_path.exists():
        raise FileNotFoundError(f"高 LD 区域文件不存在: {high_ld_path}")

    # 本地临时目录保存更新 ID 的映射表
    from datetime import datetime
    tmp_root = Path("./bbj_prep_tmp").resolve()
    tmp_root.mkdir(parents=True, exist_ok=True)
    tmp = tmp_root / datetime.now().strftime("run_%Y%m%d_%H%M%S_%f")
    tmp.mkdir(parents=True, exist_ok=False)
    _progress(f"创建本地临时目录: {tmp}", progress_log)

    # 1) 读取 .fam 并构建 update-ids 表
    fam_file = Path(f"{bbj_prefix}.fam")
    _progress(f"[Step1] 读取 FAM 并生成 update-ids 表 → {fam_file}", progress_log)
    fam_df = pd.read_csv(fam_file, delim_whitespace=True, header=None, dtype=str)
    update_df = pd.DataFrame({
        "#OLD_FID": fam_df[0].astype(str),
        "OLD_IID": fam_df[1].astype(str),
        "NEW_FID": "bbj_" + fam_df[0].astype(str),
        "NEW_IID": "bbj_" + fam_df[1].astype(str),
    })
    update_path = tmp / "update_fam.txt"
    update_df.to_csv(update_path, sep="\t", header=True, index=False)

    # 2) 更新 ID，输出中间前缀：output_prefix
    _progress(f"[Step2] 更新 FID/IID 并生成中间前缀 → {output_prefix}", progress_log)
    _run_cmd([
        str(plink2),
        "--bfile", str(bbj_prefix),
        "--update-ids", str(update_path),
        "--make-bed",
        "--out", str(output_prefix),
        "--threads", str(threads),
    ], tmp / "step2_update_ids.log")

    # 3) 只保留 A/C/G/T，排除高 LD，LD 修剪并生成 no_high_ld 及 prune.in/out
    no_high_ld_prefix = Path(f"{output_prefix}.no_high_ld").resolve()
    _progress(f"[Step3] 高 LD 排除 + LD 修剪 → {no_high_ld_prefix}", progress_log)
    _run_cmd([
        str(plink2),
        "--bfile", str(output_prefix),
        "--snps-only", "just-acgt",
        "--exclude", "range", str(high_ld_path),
        "--indep-pairwise", "50", "5", "0.2",
        "--make-bed",
        "--out", str(no_high_ld_prefix),
        "--threads", str(threads),
    ], tmp / "step3_ld_prune.log")

    prune_in = Path(f"{no_high_ld_prefix}.prune.in").resolve()
    prune_out = Path(f"{no_high_ld_prefix}.prune.out").resolve()
    if not prune_in.exists() or not prune_out.exists():
        raise FileNotFoundError("未生成 prune.in/prune.out，请检查 plink2 日志。")

    # 4) PCA（与用户示例保持一致：freq counts + pca allele-wts approx）
    _progress("[Step4] 计算 PCA（approx + allele-wts）", progress_log)
    _run_cmd([
        str(plink2),
        "--bfile", str(no_high_ld_prefix),
        "--extract", str(prune_in),
        "--freq", "counts",
        "--pca", "allele-wts", "approx",
        "--out", f"{no_high_ld_prefix}.prune.pca",
        "--threads", str(threads),
    ], tmp / "step4_pca.log")

    # PCA 输出文件路径
    eigen_base = Path(f"{no_high_ld_prefix}.prune.pca").resolve()
    eigenvec_path = Path(str(eigen_base) + ".eigenvec").resolve()
    eigenval_path = Path(str(eigen_base) + ".eigenval").resolve()
    if not eigenvec_path.exists() or not eigenval_path.exists():
        raise FileNotFoundError("未找到 PCA 输出（eigenvec/eigenval），请检查 plink2 日志。")

    # 自动清理临时目录
    if not keep_tmp:
        shutil.rmtree(tmp, ignore_errors=True)
        _progress(f"已删除临时目录: {tmp}", progress_log)

    _progress("BBJ PCA 输入准备完成 ✓", progress_log)
    return no_high_ld_prefix, prune_in, prune_out, eigenvec_path, eigenval_path



def plot_pca_pairwise_pdf(
    eigenvec_path: str | Path,
    eigenval_path: str | Path,
    output_pdf: str = "bbj_pca_pairwise_plots.pdf",
) -> Path:
    """
    基于 prepare_bbj_pca_inputs() 产出的 PCA 结果，绘制“PC 成对散点图”多页 PDF。

    参数
    ----
    eigenvec_path : str | Path
        *.eigenvec 文件路径（包含列 FID、IID、PC1、PC2...）。
    eigenval_path : str | Path
        *.eigenval 文件路径（每行一个特征值）。
    output_pdf : str, default "bbj_pca_pairwise_plots.pdf"
        输出 PDF 文件名（写到当前工作目录）。

    说明
    ----
    - 默认生成 PC1 vs PC2、PC3 vs PC4、...、PC9 vs PC10 共 5 页（若 PC 数不足则按可用对数生成）。
    - 散点统一为灰色底（与用户示例一致），右侧添加小圆点图例。
    - X/Y 轴标题会附上对应 PC 的方差解释比例（%），基于 eigenval 的相对贡献计算。
    - 采用“先存 PNG，再嵌入 PDF”的方式以确保渲染一致性（与给定示例保持一致）。

    返回
    ----
    Path
        生成的 PDF 文件的绝对路径。
    """
    from matplotlib.lines import Line2D

    eigenvec_path = Path(eigenvec_path).expanduser().resolve()
    eigenval_path = Path(eigenval_path).expanduser().resolve()
    out_path = Path(output_pdf).resolve()

    if not eigenvec_path.exists():
        raise FileNotFoundError(f"未找到 eigenvec 文件: {eigenvec_path}")
    if not eigenval_path.exists():
        raise FileNotFoundError(f"未找到 eigenval 文件: {eigenval_path}")

    # 读取 PCA 数据
    pca_df = pd.read_csv(eigenvec_path, delim_whitespace=True)
    eigenval_df = pd.read_csv(eigenval_path, delim_whitespace=True, header=None)

    # 风格设置（与示例一致）
    plt.style.use("default")
    sns.set_theme(style="whitegrid", context="talk", font_scale=1.2)

    # 可用 PC 的最大编号
    # eigenvec 通常具有列 ['FID','IID','PC1','PC2',...]
    pc_cols = [c for c in pca_df.columns if str(c).startswith("PC")]
    if len(pc_cols) < 2:
        raise ValueError("eigenvec 文件中可用的 PC 列少于 2 列，无法绘制成对散点图。")

    # 能绘制的最大 PC 编号（偶数），例如 10 表示绘制到 PC10
    max_pc_index = min(10, len(pc_cols))
    if max_pc_index % 2 == 1:
        max_pc_index -= 1  # 确保成对

    # 预计算方差解释比例
    # eigenval_df 的第 1 列为特征值，取前 max_pc_index 项
    eigvals = eigenval_df.iloc[:max_pc_index, 0].astype(float)
    total = float(eigenval_df[0].sum())
    if total == 0:
        # 极端情况下避免除零
        var_ratios = [0.0 for _ in range(max_pc_index)]
    else:
        var_ratios = [(v / total) * 100.0 for v in eigvals]

    # 逐页绘制：PC1 vs PC2, PC3 vs PC4, ...
    with PdfPages(out_path) as pdf:
        for i in range(0, max_pc_index, 2):
            pc_x = f"PC{i+1}"
            pc_y = f"PC{i+2}"
            pc_x_var = var_ratios[i] if i < len(var_ratios) else 0.0
            pc_y_var = var_ratios[i+1] if (i+1) < len(var_ratios) else 0.0

            # 散点图
            plt.figure(figsize=(10, 8))
            ax = sns.scatterplot(
                data=pca_df,
                x=pc_x, y=pc_y,
                color="#C0C0C0",
                s=40, alpha=0.85,
                edgecolor="black", linewidth=0.4,
            )

            # 右侧点状图例（与示例一致）
            bbj_dot = Line2D(
                [0], [0],
                marker="o",
                color="w",
                label="BBJ",
                markerfacecolor="#C0C0C0",
                markeredgecolor="black",
                markersize=8,
                linewidth=0,
            )
            ax.legend(
                handles=[bbj_dot],
                loc="upper left",
                bbox_to_anchor=(1.02, 1),
                frameon=False,
            )

            plt.xlabel(f"{pc_x} ({pc_x_var:.2f}%)", fontsize=16)
            plt.ylabel(f"{pc_y} ({pc_y_var:.2f}%)", fontsize=16)
            plt.title(f"PCA: {pc_x} vs {pc_y}", fontsize=18, weight="bold")

            plt.grid(True, linestyle="--", linewidth=0.5, color="gray", alpha=0.4)
            plt.axhline(0, color="gray", linewidth=0.6, linestyle="--", alpha=0.5)
            plt.axvline(0, color="gray", linewidth=0.6, linestyle="--", alpha=0.5)

            plt.tight_layout()

            # 保存位图 PNG，再嵌入 PDF（使用 plt.imread 避免额外全局导入）
            tmp_png = f"tmp_pca_plot_{i+1}_{i+2}.png"
            plt.savefig(tmp_png, dpi=600)
            plt.close()

            fig = plt.figure(figsize=(10, 8))
            img = plt.imread(tmp_png)
            plt.imshow(img)
            plt.axis("off")
            pdf.savefig(fig)
            plt.close()

            # 清理临时 PNG
            try:
                os.remove(tmp_png)
            except OSError:
                pass

    return out_path


def run_bbj_projection(
    my_bed_prefix: str,
    bbj_bed_prefix: str,
    bbj_prune_in: str = "bbj.maf.no_high_ld.prune.in",
    my_prefix_out: str = "cteph_agp3k",
    bbj_prefix_out: str = "bbj",
    bbj_pca_acount: str = "bbj.maf.no_high_ld.prune.pca.acount",
    bbj_pca_eigenvec_allele: str = "bbj.maf.no_high_ld.prune.pca.eigenvec.allele",
    threads: int = 32,
    plink2_path: str = "/home/b/b37974/plink2",
    plink1_path: str = "/home/b/b37974/plink",
    keep_tmp: bool = True,
) -> Path:
    """
    运行 BBJ 投影（Projection）流程，将“我的数据”按 BBJ 的修剪集合进行对齐并投影到 BBJ PCA 空间。
    除最终返回的投影结果（*.sscore）外，**所有中间产物与日志均写入本地 tmp 目录**。
    **注意：内部所用命令（cmd）与用户给定示例一致，不做修改，仅封装为工具函数。**

    参数
    ----
    my_bed_prefix : str
        你的样本数据（已完成 QC）的 PLINK 前缀。
    bbj_bed_prefix : str
        BBJ 的 PLINK 前缀（已通过 bbj_prepare/bbj_pca 步骤）。
    bbj_prune_in : str
        BBJ 的 prune.in 文件路径。
    my_prefix_out : str
        你的数据在“BBJ 修剪集合”下的输出前缀基名（仅用于命名，实际中间文件写入 tmp）。
    bbj_prefix_out : str
        BBJ 数据在“共享SNP集合过滤”后的输出前缀基名（仅用于命名，实际中间文件写入 tmp）。
    bbj_pca_acount : str
        BBJ PCA 的 acount 文件路径，用于 --read-freq。
    bbj_pca_eigenvec_allele : str
        BBJ PCA 的 eigenvec.allele 文件路径，用于 --score。
    threads : int
        线程数量，传递给 plink2。
    plink2_path : str
        plink2 可执行文件路径。
    plink1_path : str
        plink1 可执行文件路径。
    keep_tmp : bool, default True
        是否保留 tmp 目录；False 时流程结束后自动删除。

    返回
    ----
    Path
        最终生成的投影得分文件（*.sscore）的绝对路径（位于当前工作目录）。
    """
    from pathlib import Path
    import pandas as pd
    import shutil
    from datetime import datetime

    # 创建本地持久化 tmp 目录（时间戳子目录）
    tmp_root = Path("./bbj_projection_tmp").resolve()
    tmp_root.mkdir(parents=True, exist_ok=True)
    tmp = tmp_root / datetime.now().strftime("run_%Y%m%d_%H%M%S_%f")
    tmp.mkdir(parents=True, exist_ok=False)

    progress_log = (tmp / "bbj_projection.progress.log").resolve()
    _progress(f"启动 BBJ 投影流程（tmp: {tmp}）", progress_log)

    # 路径与依赖检查
    plink2 = Path(plink2_path).expanduser().resolve()
    plink1 = Path(plink1_path).expanduser().resolve()
    if not plink2.exists():
        raise FileNotFoundError(f"plink2 不存在: {plink2}")
    if not plink1.exists():
        raise FileNotFoundError(f"plink 不存在: {plink1}")

    my_prefix = Path(my_bed_prefix).expanduser().resolve()
    bbj_prefix = Path(bbj_bed_prefix).expanduser().resolve()
    for ext in (".bed", ".bim", ".fam"):
        if not (my_prefix.parent / f"{my_prefix.name}{ext}").exists():
            raise FileNotFoundError(f"缺少你的输入文件: {my_prefix}{ext}")
        if not (bbj_prefix.parent / f"{bbj_prefix.name}{ext}").exists():
            raise FileNotFoundError(f"缺少BBJ输入文件: {bbj_prefix}{ext}")

    prune_in_path = Path(bbj_prune_in).expanduser().resolve()
    if not prune_in_path.exists():
        raise FileNotFoundError(f"未找到 BBJ prune.in: {prune_in_path}")

    acount_path = Path(bbj_pca_acount).expanduser().resolve()
    eigenvec_allele_path = Path(bbj_pca_eigenvec_allele).expanduser().resolve()
    if not acount_path.exists():
        raise FileNotFoundError(f"未找到 BBJ PCA acount 文件: {acount_path}")
    if not eigenvec_allele_path.exists():
        raise FileNotFoundError(f"未找到 BBJ PCA eigenvec.allele 文件: {eigenvec_allele_path}")

    # 1) 用 BBJ 的 prune.in 抽取“我的数据” → tmp/{my_prefix_out}.bbj_pruned
    _progress("[Step1] 对我的数据按 BBJ prune.in 进行抽取", progress_log)
    my_pruned_prefix = tmp / f"{my_prefix_out}.bbj_pruned"
    cmd_bed_prune = [
        str(plink2),
        "--bfile", str(my_prefix),
        "--extract", str(prune_in_path),
        "--make-bed",
        "--out", str(my_pruned_prefix),
        "--threads", str(threads),
    ]
    _run_cmd(cmd_bed_prune, tmp / "step1_my_prune.log")

    # 2) 生成共享 SNP 列表（基于我的 .bim）
    _progress("[Step2] 生成共享 SNP 列表（基于我的 .bim）", progress_log)
    bim_path = Path(str(my_pruned_prefix) + ".bim")
    if not bim_path.exists():
        raise FileNotFoundError(f"未找到 BIM 文件: {bim_path}")
    bim_df = pd.read_csv(bim_path, sep="\t", header=None)
    snp_list_path = tmp / f"{my_prefix_out}.bbj_pruned.snps.txt"
    bim_df[1].to_csv(snp_list_path, index=False, header=False)

    # 3) 用共享 SNP 列表抽取 BBJ 数据 → tmp/{bbj_prefix_out}.bbj_pruned
    _progress("[Step3] 对 BBJ 数据按共享 SNP 列表进行抽取", progress_log)
    bbj_pruned_prefix = tmp / f"{bbj_prefix_out}.bbj_pruned"
    cmd_share_snps = [
        str(plink2),
        "--bfile", str(bbj_prefix),
        "--extract", str(snp_list_path),
        "--make-bed",
        "--out", str(bbj_pruned_prefix),
        "--threads", str(threads),
    ]
    _run_cmd(cmd_share_snps, tmp / "step3_bbj_extract.log")

    # 4) 合并我的数据与 BBJ（使用 plink1）→ tmp/{my}.{bbj}.bbj_pruned.merged
    _progress("[Step4] 合并我的数据与 BBJ（plink --bmerge）", progress_log)
    merged_prefix = tmp / f"{my_prefix_out}.{bbj_prefix_out}.bbj_pruned.merged"
    cmd_merge = [
        str(plink1),
        "--bfile", str(my_pruned_prefix),
        "--bmerge", str(bbj_pruned_prefix),
        "--keep-allele-order",
        "--make-bed",
        "--out", str(merged_prefix),
    ]
    _run_cmd(cmd_merge, tmp / "step4_merge.log")

    # 5) 投影到 BBJ PCA 空间（plink2 --score）
    _progress("[Step5] 投影到 BBJ PCA 空间（plink2 --score）", progress_log)
    # 最终产物放在工作目录（非 tmp）
    final_out_prefix = Path(f"{my_prefix_out}.{bbj_prefix_out}.projection").resolve()
    cmd_projection = [
        str(plink2),
        "--bfile", str(merged_prefix),
        "--read-freq", str(acount_path),
        "--extract", str(snp_list_path),
        "--score", str(eigenvec_allele_path), "2", "6", "header-read", "no-mean-imputation", "variance-standardize", "list-variants",
        "--score-col-nums", "7-16",
        "--out", str(final_out_prefix),
        "--threads", str(threads),
    ]
    _run_cmd(cmd_projection, tmp / "step5_projection.log")

    # 返回 .sscore 的路径
    sscore_path = Path(str(final_out_prefix) + ".sscore").resolve()
    if not sscore_path.exists():
        raise FileNotFoundError(f"未生成投影得分文件: {sscore_path}")

    _progress("BBJ 投影流程完成 ✓", progress_log)

    # 自动清理 tmp 目录
    if not keep_tmp:
        shutil.rmtree(tmp, ignore_errors=True)
        _progress(f"已删除临时目录: {tmp}", None)

    return sscore_path



def plot_projection_pairwise_pdf(
    sscore_path: str | Path,
    case_prefix: str = "PHOM",
    case_name: str = "CTEPH",
    control_name: str = "AGP3K",
    output_pdf: str = "bbj_projection_pairwise_plots.pdf",
) -> Path:
    """
    基于投影打分文件（*.sscore）的 PC*_AVG 列，绘制“PC 成对散点图”多页 PDF。

    与 plot_pca_pairwise_pdf 的不同：
      1) 不显示方差解释比例；
      2) 需要根据 IID 分组三类：case(以 case_prefix 开头)、bbj(以 "bbj" 开头)、ctrl(其余)；
      3) 绘制顺序（由下到上）：bbj → ctrl → case（并设置一定透明度）。

    参数
    ----
    sscore_path : str | Path
        投影得分文件路径（列名包含：#FID, IID, PC1_AVG..PC10_AVG 等）。
    case_prefix : str, default "PHOM"
        用于识别病例样本的 IID 前缀。
    case_name : str, default "CTEPH"
        病例在图例中的显示名称。
    control_name : str, default "AGP3K"
        对照在图例中的显示名称。
    output_pdf : str, default "bbj_projection_pairwise_plots.pdf"
        输出 PDF 文件名（写到当前工作目录）。

    返回
    ----
    Path
        生成的 PDF 文件绝对路径。
    """
    from matplotlib.lines import Line2D
    import re

    sscore_path = Path(sscore_path).expanduser().resolve()
    out_path = Path(output_pdf).resolve()
    if not sscore_path.exists():
        raise FileNotFoundError(f"未找到 sscore 文件: {sscore_path}")

    # 读取 sscore（空白分隔，首行为表头）
    df = pd.read_csv(sscore_path, delim_whitespace=True, header=0, dtype={ "IID": str })
    if "IID" not in df.columns:
        raise ValueError("sscore 缺少 'IID' 列，请检查文件格式。")

    # 识别可用的 PC*_AVG 列
    pc_cols = [c for c in df.columns if re.fullmatch(r"PC\d+_AVG", str(c)) is not None]
    if len(pc_cols) < 2:
        raise ValueError("sscore 中可用的 PC*_AVG 列少于 2 列，无法绘制成对散点图。")

    # 按编号排序（PC1_AVG, PC2_AVG, ...）
    def pc_index(col: str) -> int:
        m = re.match(r"PC(\d+)_AVG", col)
        return int(m.group(1)) if m else 1_000_000
    pc_cols.sort(key=pc_index)

    # 生成分组标签
    def classify(iid: str) -> str:
        s = str(iid)
        if s.startswith("bbj"):
            return "BBJ"
        elif s.startswith(case_prefix):
            return case_name
        else:
            return control_name

    df["Group"] = df["IID"].astype(str).map(classify)

    # 主题样式
    plt.style.use("default")
    sns.set_theme(style="whitegrid", context="talk", font_scale=1.2)

    # 最多绘制到 PC10_AVG（若列更少则按可用对数）
    max_pairs = min(10, len(pc_cols))
    if max_pairs % 2 == 1:
        max_pairs -= 1

    # 绘制顺序与配色（由下到上）
    order = ["BBJ", control_name, case_name]
    palette = {
        "BBJ": "#BDBDBD",          # 浅灰
        control_name: "#6BAED6",   # 蓝
        case_name: "#FB6A4A",      # 红
    }
    alpha_map = {"BBJ": 0.55, control_name: 0.70, case_name: 0.90}

    # 自定义图例
    legend_handles = []
    for key in order:
        legend_handles.append(
            Line2D(
                [0], [0], marker="o", linestyle="",
                label=key,
                markerfacecolor=palette.get(key, "#C0C0C0"),
                markeredgecolor="black", markersize=8
            )
        )

    # 逐页绘制
    with PdfPages(out_path) as pdf:
        for i in range(0, max_pairs, 2):
            pc_x = pc_cols[i]
            pc_y = pc_cols[i+1]

            plt.figure(figsize=(10, 8))
            ax = plt.gca()

            # 由下到上分层绘制
            for grp in order:
                sub = df[df["Group"] == grp]
                if sub.empty:
                    continue
                sns.scatterplot(
                    data=sub, x=pc_x, y=pc_y, ax=ax,
                    s=40, alpha=alpha_map.get(grp, 0.8),
                    color=palette.get(grp, "#C0C0C0"),
                    edgecolor="black", linewidth=0.4
                )

            # 图例放右侧
            ax.legend(
                handles=legend_handles,
                loc="upper left",
                bbox_to_anchor=(1.02, 1),
                frameon=False
            )

            # 轴与标题（不显示方差解释）
            # PC*_AVG → PC*
            def labelize(col: str) -> str:
                idx = pc_index(col)
                return f"PC{idx}"
            plt.xlabel(labelize(pc_x), fontsize=16)
            plt.ylabel(labelize(pc_y), fontsize=16)
            plt.title(f"Projection: {labelize(pc_x)} vs {labelize(pc_y)}", fontsize=18, weight="bold")

            plt.grid(True, linestyle="--", linewidth=0.5, color="gray", alpha=0.4)
            plt.axhline(0, color="gray", linewidth=0.6, linestyle="--", alpha=0.5)
            plt.axvline(0, color="gray", linewidth=0.6, linestyle="--", alpha=0.5)

            plt.tight_layout()

            # 与现有函数一致：保存 PNG 再嵌入 PDF
            tmp_png = f"tmp_proj_plot_{i+1}_{i+2}.png"
            plt.savefig(tmp_png, dpi=600)
            plt.close()

            fig = plt.figure(figsize=(10, 8))
            img = plt.imread(tmp_png)
            plt.imshow(img)
            plt.axis("off")
            pdf.savefig(fig)
            plt.close()

            try:
                os.remove(tmp_png)
            except OSError:
                pass

    return out_path



def plot_projection_pc1_pc2_and_select(
    sscore_path: str | Path,
    case_prefix: str = "PHOM",
    case_name: str = "CTEPH",
    control_name: str = "AGP3K",
    prefix_out: str = "cteph_agp3k",
    rect_xlim: tuple[float, float] = (-0.025, 0.016),
    rect_ylim: tuple[float, float] = (-0.025, 0.025),
    output_png: str | None = None,
) -> tuple[Path, Path]:
    """
    基于投影打分文件（*.sscore）的 PC1_AVG/PC2_AVG，绘制“主图 + 边际KDE”的图像，并在主图上以矩形框筛选样本。
    - 分组三类：BBJ（IID 以 "bbj" 开头）、病例（IID 以 case_prefix 开头）、对照（其余）。
    - 绘制顺序（由下到上）：BBJ → 对照（control_name） → 病例（case_name）。
    - 仅对病例/对照绘制边际 KDE（顶部/右侧）。
    - 将矩形内的病例/对照样本（FID/IID）输出到 `{prefix_out}.bbj.projection.sample.keep.txt`。

    参数
    ----
    sscore_path : str | Path
        投影得分文件路径，需包含列：#FID、IID、PC1_AVG、PC2_AVG（以及可选 PC3_AVG...）。
    case_prefix : str, default "PHOM"
        病例样本的 IID 前缀。
    case_name : str, default "CTEPH"
        病例在图例中的标签名。
    control_name : str, default "AGP3K"
        对照在图例中的标签名。
    prefix_out : str, default "cteph_agp3k"
        用于输出文件名前缀（PNG 及 keep.txt）。
    rect_xlim : (float, float), default (-0.025, 0.016)
        主图中矩形框的 X 轴范围。
    rect_ylim : (float, float), default (-0.025, 0.025)
        主图中矩形框的 Y 轴范围。
    output_png : str | None, default None
        可选自定义 PNG 文件名；若为 None，则使用 `{prefix_out}.bbj_projection.pc1_pc2.png`。

    返回
    ----
    (Path, Path)
        (png_path, keep_txt_path) 的绝对路径。
    """
    import re
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    import seaborn as sns
    from matplotlib.lines import Line2D
    from matplotlib.patches import Rectangle
    import matplotlib.gridspec as gridspec

    sscore_path = Path(sscore_path).expanduser().resolve()
    if not sscore_path.exists():
        raise FileNotFoundError(f"未找到 sscore 文件: {sscore_path}")

    # 读取数据（空白分隔，保留表头）
    df = pd.read_csv(sscore_path, delim_whitespace=True, header=0, dtype={"IID": str})
    # 兼容 #FID / FID 两种写法
    fid_col = "#FID" if "#FID" in df.columns else ("FID" if "FID" in df.columns else None)
    if fid_col is None or "IID" not in df.columns:
        raise ValueError("sscore 缺少 '#FID'/'FID' 或 'IID' 列，请检查文件格式。")

    # 必需坐标列
    for col in ("PC1_AVG", "PC2_AVG"):
        if col not in df.columns:
            raise ValueError(f"sscore 缺少 '{col}' 列，请确认 --score 输出包含 PC*_AVG 列。")

    # 分组与配色
    def classify(iid: str) -> str:
        iid = str(iid)
        if iid.startswith("bbj"):
            return "BBJ"
        elif iid.startswith(case_prefix):
            return case_name
        else:
            return control_name

    df["Group"] = df["IID"].astype(str).map(classify)

    palette = {
        "BBJ": "#C0C0C0",         # 浅灰
        control_name: "#4DBBD5",  # 蓝绿
        case_name: "#E64B35",     # 红
    }
    # 绘制顺序（由下到上）
    draw_order = ["BBJ", control_name, case_name]
    alpha_map = {"BBJ": 0.55, control_name: 0.75, case_name: 0.90}

    # 画布与主题
    plt.style.use("default")
    sns.set_theme(style="whitegrid", context="talk", font_scale=1.2)

    fig = plt.figure(figsize=(8, 8))
    gs = gridspec.GridSpec(
        2, 2,
        width_ratios=[4, 1],
        height_ratios=[1, 4],
        wspace=0.05, hspace=0.05
    )

    ax_main = plt.subplot(gs[1, 0])
    ax_top  = plt.subplot(gs[0, 0], sharex=ax_main)
    ax_right= plt.subplot(gs[1, 1], sharey=ax_main)

    # 主散点：按顺序分层绘制
    for grp in draw_order:
        sub = df[df["Group"] == grp]
        if sub.empty:
            continue
        ax_main.scatter(
            sub["PC1_AVG"], sub["PC2_AVG"],
            c=palette.get(grp, "#C0C0C0"),
            s=10, alpha=alpha_map.get(grp, 0.8), linewidth=0
        )

    ax_main.set_xlabel("PC1_AVG", fontsize=14)
    ax_main.set_ylabel("PC2_AVG", fontsize=14)
    ax_main.grid(True, linestyle="--", linewidth=0.5, alpha=0.3)

    # 图例（右上角）
    legend_elements = [
        Line2D([0], [0], marker="o", color="w", label=case_name,   markerfacecolor=palette[case_name],   markersize=7),
        Line2D([0], [0], marker="o", color="w", label=control_name,markerfacecolor=palette[control_name],markersize=7),
        Line2D([0], [0], marker="o", color="w", label="BBJ",       markerfacecolor=palette["BBJ"],       markersize=7),
    ]
    ax_main.legend(handles=legend_elements, loc="upper left", frameon=False, fontsize=11)

    # 边际 KDE：仅病例/对照
    ax_top.tick_params(bottom=False, labelbottom=False)
    ax_right.tick_params(left=False, labelleft=False)

    for color, label in [(palette[case_name], case_name), (palette[control_name], control_name)]:
        subset = df[df["Group"] == label]
        if not subset.empty:
            try:
                sns.kdeplot(x=subset["PC1_AVG"], ax=ax_top, color=color, fill=True, linewidth=1, alpha=0.5)
                sns.kdeplot(y=subset["PC2_AVG"], ax=ax_right, color=color, fill=True, linewidth=1, alpha=0.5)
            except Exception:
                # KDE 可能因样本过少或全部相同值而失败，忽略并继续
                pass

    ax_top.set_ylabel("")
    ax_top.set_yticks([])
    ax_right.set_xlabel("")
    ax_right.set_xticks([])

    # 矩形框（筛选区域）
    rect = Rectangle(
        (rect_xlim[0], rect_ylim[0]),
        rect_xlim[1] - rect_xlim[0],
        rect_ylim[1] - rect_ylim[0],
        linewidth=1.5, edgecolor="black", facecolor="none", linestyle="--"
    )
    ax_main.add_patch(rect)

    # 选择矩形内的病例/对照样本
    in_rect = (
        (df["PC1_AVG"] >= rect_xlim[0]) & (df["PC1_AVG"] <= rect_xlim[1]) &
        (df["PC2_AVG"] >= rect_ylim[0]) & (df["PC2_AVG"] <= rect_ylim[1])
    )
    selected = df.loc[in_rect & df["Group"].isin([case_name, control_name]), [fid_col, "IID", "Group"]].copy()

    # 数量注释
    counts = selected["Group"].value_counts().to_dict()
    n_case = counts.get(case_name, 0)
    n_ctrl = counts.get(control_name, 0)
    text_str = f"{case_name}: n={n_case}\n{control_name}: n={n_ctrl}"
    ax_main.text(
        0.70, 0.02, text_str,
        transform=ax_main.transAxes,
        fontsize=12,
        va="bottom", ha="left",
        bbox=dict(boxstyle="round,pad=0.3", edgecolor="gray", facecolor="white", alpha=0.7)
    )

    plt.suptitle(f"{case_name} + {control_name} based on BBJ Projection", fontsize=16, weight="bold")
    plt.tight_layout(rect=[0, 0, 1, 0.98])

    # 输出路径
    png_path = Path(output_png) if output_png else Path(f"{prefix_out}.bbj_projection.pc1_pc2.png")
    keep_txt_path = Path(f"{prefix_out}.bbj.projection.sample.keep.txt")
    png_path = png_path.resolve()
    keep_txt_path = keep_txt_path.resolve()

    # 保存高清 PNG
    plt.savefig(png_path, dpi=600, bbox_inches="tight")
    plt.close(fig)

    # 输出 keep.txt（仅病例/对照）
    selected.rename(columns={fid_col: "#FID"}, inplace=True)
    selected.loc[:, ["#FID", "IID"]].to_csv(keep_txt_path, sep="\t", index=False, header=True)

    return png_path, keep_txt_path


def export_subset_bed(
    bed_prefix: str | Path,
    keep_txt_path: str | Path,
    out_prefix: str | Path | None = None,
    plink2_path: str = "/home/b/b37974/plink2",
    threads: int = 32,
) -> Path:
    """
    使用 plink2 --keep 基于样本清单导出子集基因型（bed/bim/fam）。

    参数
    ----
    bed_prefix : str | Path
        输入 PLINK 前缀（不带扩展名），需存在 .bed/.bim/.fam。
    keep_txt_path : str | Path
        样本保留清单，至少包含两列 (#FID, IID)。允许带表头。
        （若存在表头，将自动转换为无表头两列再传给 plink2）
    out_prefix : str | Path | None, default None
        子集输出前缀（不带扩展名）。若为 None，则默认使用
        `./{basename(bed_prefix)}.keep`。
    plink2_path : str, default "/home/b/b37974/plink2"
        plink2 可执行文件路径。
    threads : int, default 32
        plink2 线程数。

    返回
    ----
    Path
        子集数据的前缀绝对路径（不带扩展名）。
    """
    from datetime import datetime

    bed_prefix = Path(bed_prefix).expanduser().resolve()
    keep_txt_path = Path(keep_txt_path).expanduser().resolve()
    plink2 = Path(plink2_path).expanduser().resolve()

    # 基本校验
    if not plink2.exists():
        raise FileNotFoundError(f"plink2 不存在: {plink2}")
    for ext in (".bed", ".bim", ".fam"):
        if not (bed_prefix.parent / f"{bed_prefix.name}{ext}").exists():
            raise FileNotFoundError(f"缺少输入文件: {bed_prefix}{ext}")
    if not keep_txt_path.exists():
        raise FileNotFoundError(f"未找到 keep 文件: {keep_txt_path}")

    # 输出前缀
    out_prefix = Path(out_prefix) if out_prefix is not None else Path(f"{bed_prefix.name}.keep")
    out_prefix = out_prefix.resolve()

    # 构建临时目录与日志
    tmp_root = Path("./bbj_projection_tmp").resolve()
    tmp_root.mkdir(parents=True, exist_ok=True)
    tmp = tmp_root / datetime.now().strftime("run_%Y%m%d_%H%M%S_%f")
    tmp.mkdir(parents=True, exist_ok=False)
    progress_log = tmp / "subset.progress.log"
    _progress(f"启动子集导出：bed_prefix={bed_prefix}, keep={keep_txt_path}", progress_log)

    # 将 keep 文件规范化为两列无表头文本（以防表头或额外列导致报错）
    # 兼容制表/空白分隔
    try:
        _progress("规范化 keep 样本清单（两列无表头）", progress_log)
        df_keep = pd.read_csv(keep_txt_path, sep=None, engine="python", dtype=str)
        # 尝试自动识别列名
        cand_cols = [c for c in df_keep.columns if str(c).strip() in ("#FID", "FID", "fid", "Fid")]
        if cand_cols:
            fid_col = cand_cols[0]
        else:
            fid_col = df_keep.columns[0]

        cand_cols2 = [c for c in df_keep.columns if str(c).strip() in ("IID", "iid", "Iid")]
        if cand_cols2:
            iid_col = cand_cols2[0]
        else:
            iid_col = df_keep.columns[1] if df_keep.shape[1] > 1 else df_keep.columns[0]

        keep_nhdr = tmp / "keep.noheader.txt"
        df_keep[[fid_col, iid_col]].to_csv(keep_nhdr, sep="\t", index=False, header=False)
    except Exception as e:
        _progress(f"规范化 keep 文件失败（将直接传入原文件）：{e}", progress_log)
        keep_nhdr = keep_txt_path  # 回退：直接使用原文件

    # 运行 plink2 --keep
    _progress("运行 plink2 --keep 导出子集", progress_log)
    cmd = [
        str(plink2),
        "--bfile", str(bed_prefix),
        "--keep", str(keep_nhdr),
        "--make-bed",
        "--out", str(out_prefix),
        "--threads", str(threads),
    ]
    _run_cmd(cmd, tmp / "subset_keep.log")

    # 校验产物
    for ext in (".bed", ".bim", ".fam"):
        p = Path(str(out_prefix) + ext)
        if not p.exists():
            raise FileNotFoundError(f"子集导出失败，缺少文件: {p}")

    _progress(f"子集导出完成，输出前缀：{out_prefix}", progress_log)
    return out_prefix
