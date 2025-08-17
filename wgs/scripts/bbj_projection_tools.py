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