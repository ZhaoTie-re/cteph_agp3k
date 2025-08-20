"""
模块名称：pca_qc_tools
=================================

【概述】
本模块提供两个与群体结构相关的常用分析工具：
1) `run_prune_and_pca`：基于 PLINK2 的标准流程（MAF 过滤 → 去除高 LD 区域 → LD-Prune → PCA）生成 PCA 所需的文件。
2) `plot_pca_pairwise_pdf`：将 PCA 结果（*.eigenvec/*.eigenval）绘制为“PC 成对散点图”的多页 PDF（每页两个 PC：PC1 vs PC2、PC3 vs PC4…）。

【适用场景】
- WGS/WES/芯片数据在完成基础变异与个体质控后，进行群体结构探索与可视化；
- 与外部参考（如 BBJ/ToMMo 等）进行 PCA 位置比较（本模块仅产出 PCA 与散点图，跨面板投影需在外部完成）。

【输入与输出（核心函数）】
- run_prune_and_pca(
    bed_prefix: str,
    maf_threshold: float = 0.05,
    high_ld: str,
    threads: int = 16,
    output_prefix: str = "cteph_agp3k",
    plink2_path: str = "/home/b/b37974/plink2",
  ) -> Tuple[Path, Path, Path, Path, Path, Path]
  - 依赖：存在 `{bed_prefix}.bed/.bim/.fam` 与高 LD 区域 BED 文件。
  - 产物（按顺序返回）：
    1) `{output_prefix}.no_high_ld`（中间 bed/bim/fam 前缀）
    2) `{output_prefix}.no_high_ld.prune.in`
    3) `{output_prefix}.no_high_ld.prune.out`
    4) `{output_prefix}.no_high_ld.prune.pca.eigenvec`
    5) `{output_prefix}.no_high_ld.prune.pca.eigenval`
    6) `{output_prefix}.no_high_ld.prune.pca.eigenvec.allele`
  - 日志：写入 `{output_prefix}.logs/` 目录，分步保存。

- plot_pca_pairwise_pdf(
    eigenvec_file: str,
    eigenval_file: str,
    case_prefix: str,
    case_name: str = "CTEPH",
    control_name: str = "AGP3K",
    output_pdf: str = "pca_pairwise_plots.pdf",
  ) -> Path
  - 行为：读取含表头的 eigenvec（自动把首列 `#FID` 规范为 `FID`），识别 PC 列并与 eigenval 对齐；
    根据 IID 前缀将样本分为“病例/对照”，生成多页 PDF。

【依赖/环境】
- 运行时：Python ≥ 3.9（已使用 `from __future__ import annotations`）。
- 第三方包：pandas、matplotlib、seaborn。
- 外部可执行：PLINK2（提供完整路径或确保在 PATH 中）。

【约定与注意事项】
- MAF 过滤与去除高 LD 区域只保留 A/C/G/T 的单核苷酸位点（`--snps-only just-acgt`）。
- 高 LD 区域文件使用 `--exclude range` 语义，需为三列（chr, start, end）的 BED/范围格式。
- PCA 默认计算前 10 个主成分并输出 allele-wts；若需要更多 PC，请在函数内调整参数。
- 所有输出路径默认采用绝对路径生成，以避免工作目录变化导致的相对路径问题。
- 失败时抛出异常，并在日志中记录原始命令与 stderr/stdout，便于溯源。

【快速示例】
1) 运行 PCA 前处理：
   >>> run_prune_and_pca(
   ...     bed_prefix="/path/to/data/cohort",
   ...     maf_threshold=0.05,
   ...     high_ld="/path/to/high-LD-regions-hg38-GRCh38.txt",
   ...     threads=32,
   ...     output_prefix="cteph_agp3k",
   ...     plink2_path="/home/b/b37974/plink2",
   ... )

2) 生成配对 PC 图：
   >>> plot_pca_pairwise_pdf(
   ...     eigenvec_file="cteph_agp3k.no_high_ld.prune.pca.eigenvec",
   ...     eigenval_file="cteph_agp3k.no_high_ld.prune.pca.eigenval",
   ...     case_prefix="PHOM",
   ...     case_name="CTEPH",
   ...     control_name="AGP3K",
   ...     output_pdf="cteph_agp3k.pca_pairwise.pdf",
   ... )

作者: ZHAO TIE

【变更记录】
- 2025-08-17：新增模块级中文文档；`plot_pca_pairwise_pdf` 明确保留 eigenvec 表头并规范化 `#FID`→`FID`；
               强化输入校验与异常提示；统一日志与绝对路径策略。

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


def run_prune_and_pca(
    bed_prefix: str,
    maf_threshold: float = 0.05,
    high_ld: str | os.PathLike[str] = "",
    threads: int = 16,
    output_prefix: str = "cteph_agp3k",
    plink2_path: str = "/home/b/b37974/plink2",
) -> Tuple[Path, Path, Path, Path, Path, Path]:
    """
    使用 PLINK2 完成 MAF 过滤 → 去除高 LD → LD-Prune → PCA 的标准流程。

    参数
    ----
    bed_prefix : str
        原始 PLINK 二进制前缀（不带扩展名，需存在 .bed/.bim/.fam）。
    maf_threshold : float, default 0.05
        过滤阈值，仅保留 MAF > 阈值的变体。
    high_ld : str | PathLike, required
        高 LD 区域 BED 文件路径（用于 --exclude range）。
    threads : int, default 16
        PLINK2 线程数。
    output_prefix : str, default "cteph_agp3k"
        输出前缀（会在不同阶段追加后缀）。
    plink2_path : str, default "/home/b/b37974/plink2"
        PLINK2 可执行文件路径。

    返回
    ----
    Tuple[Path, Path, Path, Path, Path, Path]
        依次返回：
        - no_high_ld_prefix: 去除高 LD 后的中间数据前缀（无扩展名）
        - prune_in:  LD-Prune 的保留列表文件路径 (*.prune.in)
        - prune_out: LD-Prune 的剔除列表文件路径 (*.prune.out)
        - eigenvec_file: PCA 结果的个体坐标文件 (*.eigenvec)
        - eigenval_file: PCA 结果的特征值文件 (*.eigenval)
        - eigenvec_allele_file: PCA 载荷/等位基因权重文件 (*.eigenvec.allele，使用 --pca allele-wts 生成)

    流程
    ----
    1) --bfile {bed_prefix} --maf {maf_threshold} --snps-only just-acgt --make-bed
       --out {output_prefix}.maf{maf_threshold}
    2) --bfile {output_prefix}.maf{maf_threshold} --exclude range {high_ld}
       --indep-pairwise 50 5 0.2 --make-bed --out {output_prefix}.no_high_ld
       （会生成 *.prune.in 和 *.prune.out）
    3) --bfile {output_prefix}.no_high_ld --extract {output_prefix}.no_high_ld.prune.in
       --pca allele-wts 10 --out {output_prefix}.no_high_ld.prune.pca
    """
    # ---- 前置检查 ----
    plink2 = Path(plink2_path).expanduser().resolve()
    if not plink2.exists():
        raise FileNotFoundError(f"plink2 不存在: {plink2}")

    bed_prefix_path = Path(bed_prefix).expanduser().resolve()
    for ext in (".bed", ".bim", ".fam"):
        if not (bed_prefix_path.parent / f"{bed_prefix_path.name}{ext}").exists():
            raise FileNotFoundError(f"缺少输入文件: {bed_prefix_path}{ext}")

    high_ld_path = Path(high_ld).expanduser().resolve()
    if not high_ld_path.exists():
        raise FileNotFoundError(f"高 LD BED 文件不存在: {high_ld_path}")

    # 统一使用绝对路径，避免工作目录变化导致的问题
    out_prefix_maf = Path(f"{output_prefix}.maf{maf_threshold}").resolve()
    out_prefix_no_ld = Path(f"{output_prefix}.no_high_ld").resolve()
    out_prefix_pca = Path(f"{output_prefix}.no_high_ld.prune.pca").resolve()

    # log 目录与文件
    log_dir = Path(f"{output_prefix}.logs").resolve()
    log_dir.mkdir(parents=True, exist_ok=True)

    # ---- Step 1: MAF 过滤，保留 A/C/G/T 单核苷酸位点 ----
    cmd1 = [
        str(plink2),
        "--bfile", str(bed_prefix_path),
        "--maf", str(maf_threshold),
        "--snps-only", "just-acgt",
        "--make-bed",
        "--out", str(out_prefix_maf),
        "--threads", str(threads),
    ]
    _run_cmd(cmd1, log_dir / "step1_maf_filter.log")

    # ---- Step 2: 去除高 LD + 计算 LD-Prune（生成 prune.in/out），并输出新的中间 bed ----
    cmd2 = [
        str(plink2),
        "--bfile", str(out_prefix_maf),
        "--exclude", "range", str(high_ld_path),
        "--indep-pairwise", "50", "5", "0.2",
        "--make-bed",
        "--out", str(out_prefix_no_ld),
        "--threads", str(threads),
    ]
    _run_cmd(cmd2, log_dir / "step2_high_ld_and_prune.log")

    prune_in = Path(f"{out_prefix_no_ld}.prune.in")
    prune_out = Path(f"{out_prefix_no_ld}.prune.out")
    if not prune_in.exists() or not prune_out.exists():
        raise FileNotFoundError(
            "未找到 LD-Prune 输出文件（*.prune.in / *.prune.out）。请检查 Step 2 日志。"
        )

    # ---- Step 3: 基于 prune.in 提取并做 PCA ----
    cmd3 = [
        str(plink2),
        "--bfile", str(out_prefix_no_ld),
        "--extract", str(prune_in),
        "--pca", "allele-wts", "10",
        "--out", str(out_prefix_pca),
        "--threads", str(threads),
    ]
    _run_cmd(cmd3, log_dir / "step3_pca.log")

    eigenvec_file = Path(f"{out_prefix_pca}.eigenvec")
    eigenval_file = Path(f"{out_prefix_pca}.eigenval")
    eigenvec_allele_file = Path(f"{out_prefix_pca}.eigenvec.allele")
    if (not eigenvec_file.exists()) or (not eigenval_file.exists()) or (not eigenvec_allele_file.exists()):
        raise FileNotFoundError(
            "未找到 PCA 输出文件（*.eigenvec / *.eigenval / *.eigenvec.allele）。请检查 Step 3 日志。"
        )

    # 返回：no_high_ld 前缀、prune.in/out、eigenvec/eigenval/eigenvec_allele
    return out_prefix_no_ld, prune_in, prune_out, eigenvec_file, eigenval_file, eigenvec_allele_file


def plot_pca_pairwise_pdf(
    eigenvec_file: str,
    eigenval_file: str,
    case_prefix: str,
    case_name: str = "CTEPH",
    control_name: str = "AGP3K",
    output_pdf: str | os.PathLike[str] = "pca_pairwise_plots.pdf",
) -> Path:
    """
    基于 PLINK2 的 PCA 输出（*.eigenvec / *.eigenval）绘制成对主成分散点图，并写入一个多页 PDF。

    参数
    ----
    eigenvec_file : str
        PLINK2 生成的个体坐标文件（通常含表头 #FID IID PC1..PCk；本函数将保留该表头并把列名“#FID”规范化为“FID”）。
    eigenval_file : str
        PLINK2 生成的特征值文件（每行一个特征值）。
    case_prefix : str
        用于判定病例个体（IID 以该前缀开头）。
    case_name : str, default "CTEPH"
        病例组在图例中的显示名称。
    control_name : str, default "AGP3K"
        对照组在图例中的显示名称。
    output_pdf : str | PathLike, default "pca_pairwise_plots.pdf"
        输出 PDF 路径。

    返回
    ----
    Path
        生成的 PDF 文件路径。
    """
    eigenvec_path = Path(eigenvec_file).expanduser().resolve()
    eigenval_path = Path(eigenval_file).expanduser().resolve()
    if not eigenvec_path.exists():
        raise FileNotFoundError(f"eigenvec 文件不存在: {eigenvec_path}")
    if not eigenval_path.exists():
        raise FileNotFoundError(f"eigenval 文件不存在: {eigenval_path}")

    # 读取 eigenvec：文件包含表头（首列常为“#FID”）；保留表头并规范化列名
    pca_df = pd.read_csv(eigenvec_path, delim_whitespace=True, header=0)
    # 将“#FID”规范化为“FID”
    if pca_df.columns[0].lstrip("#") == "FID":
        pca_df = pca_df.rename(columns={pca_df.columns[0]: "FID"})
    # 基本校验
    if "FID" not in pca_df.columns or "IID" not in pca_df.columns:
        raise ValueError("eigenvec 表头缺少 FID/IID 列。请检查文件格式。")

    # 识别 PC 列并确保为数值类型
    pc_label_cols = [c for c in pca_df.columns if isinstance(c, str) and c.startswith("PC")]
    if len(pc_label_cols) < 2:
        raise ValueError("eigenvec 中的主成分列不足（需要至少 PC1 与 PC2）。")
    # 转为数值，非数值强制为 NaN
    pca_df[pc_label_cols] = pca_df[pc_label_cols].apply(pd.to_numeric, errors="coerce")
    # 丢弃含有缺失 PC 值的行，避免绘图时出现异常
    pca_df = pca_df.dropna(subset=pc_label_cols, how="any")

    # 读取 eigenval：单列无表头
    eigenval_df = pd.read_csv(eigenval_path, delim_whitespace=True, header=None)
    if eigenval_df.shape[1] < 1:
        raise ValueError("eigenval 文件格式异常，期望单列表。")

    # 推断可用 PC 的数量：受 eigenval 行数与 eigenvec 的 PC 列数共同限制
    pc_indices = sorted(
        [int(c[2:]) for c in pca_df.columns if isinstance(c, str) and c.startswith("PC") and c[2:].isdigit()]
    )
    max_pcs_from_vec = len(pc_indices)
    max_pcs_from_val = eigenval_df.shape[0]
    n_pcs = max(0, min(max_pcs_from_vec, max_pcs_from_val))
    if n_pcs < 2:
        raise ValueError("可用的主成分数量不足（<2），无法绘制成对 PC 图。")

    # 分组：根据 IID 前缀
    pca_df["Group"] = pca_df["IID"].astype(str).apply(
        lambda x: case_name if str(x).startswith(str(case_prefix)) else control_name
    )

    # 解释度（百分比）
    total_var = eigenval_df.iloc[:n_pcs, 0].sum()
    var_pct = (eigenval_df.iloc[:n_pcs, 0] / total_var) * 100.0

    # 主题风格
    plt.style.use("default")
    sns.set_theme(style="whitegrid", context="talk", font_scale=1.2)

    output_pdf_path = Path(output_pdf).resolve()
    with PdfPages(output_pdf_path) as pdf:
        # 两两成对：PC1 vs PC2, PC3 vs PC4, ...
        for i in range(0, n_pcs - 1, 2):
            pc_x = f"PC{i+1}"
            pc_y = f"PC{i+2}"
            pc_x_var = float(var_pct.iloc[i]) if i < len(var_pct) else 0.0
            pc_y_var = float(var_pct.iloc[i + 1]) if (i + 1) < len(var_pct) else 0.0

            plt.figure(figsize=(10, 8))
            ax = sns.scatterplot(
                data=pca_df,
                x=pc_x,
                y=pc_y,
                hue="Group",
                hue_order=[control_name, case_name],
                palette={f"{control_name}": "#4DBBD5", f"{case_name}": "#E64B35"},
                s=40,
                alpha=0.85,
                edgecolor="black",
                linewidth=0.4,
            )

            plt.xlabel(f"{pc_x} ({pc_x_var:.2f}%)", fontsize=16)
            plt.ylabel(f"{pc_y} ({pc_y_var:.2f}%)", fontsize=16)
            plt.title(f"PCA: {pc_x} vs {pc_y}", fontsize=18, weight="bold")

            plt.legend(
                title=None,
                loc="upper left",
                bbox_to_anchor=(1.02, 1),
                borderaxespad=0,
                frameon=False,
            )

            plt.grid(True, linestyle="--", linewidth=0.5, color="gray", alpha=0.4)
            plt.axhline(0, color="gray", linewidth=0.6, linestyle="--", alpha=0.5)
            plt.axvline(0, color="gray", linewidth=0.6, linestyle="--", alpha=0.5)

            plt.tight_layout()
            pdf.savefig()
            plt.close()

    return output_pdf_path
