# -*- coding: utf-8 -*-
"""
snpeff_anno_tools.py
====================
专业文档（Chinese Documentation）

【脚本简介】
本模块提供两类与注释相关的高性能工具函数：
  1) `add_chr_prefix_to_vcf`：使用 `bcftools annotate --rename-chrs` 为 VCF 的 `#CHROM` 添加 "chr" 前缀，并仅使用 `tabix` 建立索引。
  2) `annotate_vcf_with_snpeff_tsv`：利用按染色体切分且带索引（`.tbi`）的 snpEff TSV（`*.tsv.2.gz`）对压缩 VCF（`.vcf.gz`）进行批量 INFO 注释，支持并行/串行两种模式，并最终输出带索引的注释结果。

【适用场景】
- 需要将无前缀的染色体名称（如 `1..22,X,Y,M/MT`）统一重命名为带 `chr` 前缀（如 `chr1..chr22,chrX,chrY,chrM`）。
- 已有按染色体导出的 snpEff 注释结果（TSV 格式，包含 `#CHROM, POS, REF, ALT` 键列），希望将 effect/impact/gene/hgvs 等字段合并写入 VCF 的 INFO 中。
- 处理超大规模 VCF：模块通过“按染色体切分→并行注释→合并”的策略降低内存占用、提高吞吐。

【输入与输出】
- 输入：
  - BGZF 压缩的 VCF（.vcf.gz），要求内容规范且可被 `bcftools` 正确解析。
  - snpEff 注释 TSV（`*.tsv.2.gz`）与对应的索引（`.tbi`），按染色体拆分存放于同一目录。
- 输出：
  - 注释后的 BGZF 压缩 VCF（`.vcf.gz`）与其索引（`.tbi`），默认输出到当前工作目录，文件名在原始基础上添加 `.chrprefix` 或 `.snpeff` 后缀。

【依赖与环境】
- 外部软件：`bcftools`、`tabix`（htslib），需与输入文件格式兼容。
- Python 版本：≥3.8（已使用 `pathlib`、`typing` 等特性）。
- 运行环境建议：本地或 HPC 环境均可。并行模式下推荐根据 CPU/Core 数量设置 `max_workers/threads`。

【关键设计与实现要点】
- **只用 tabix 建索引**：避免 `bcftools index` 与 `tabix -p vcf` 潜在差异，统一采用后者。
- **稳定的列映射**：将 `-c` 的列映射外部化为 `-C` 文件（`snpeff_columns.txt`），确保重复运行结果稳定、可追踪。
- **INFO 头自动管理**：若未提供 header 文件，会在 `./tmp/snpEff_info.hdr` 自动生成包含 10 个 INFO tag 的头信息。
- **并行友好**：VCF 先按染色体切分，再逐个调用 `bcftools annotate`，有效避开单次大文件操作瓶颈。
- **日志与可观测性**：统一的 `_setup_logging()` 负责到文件与控制台的双通道输出；关键步骤（切分/注释/合并/索引）均有耗时统计与错误输出。

【性能与资源建议】
- `threads` 控制 `bcftools` 内部多线程；`max_workers` 控制 Python 侧的并行任务数。二者配合可在 I/O 与 CPU 之间取得平衡。
- I/O 绑定严重时宜提高 `threads`；CPU 计算占比高时适当提高 `max_workers`（但不要超过染色体文件数）。
- 极大 VCF 建议先确保输入已建立 CSI/TBI 索引，以便 `index -s` 与 `view -r` 高效工作。

【错误处理与常见问题（FAQ）】
- **输出为空**：多由列映射不匹配或注释文件与 VCF 的 `#CHROM/REF/ALT` 不一致导致。请检查 `-C` 映射与 TSV 头是否一致，确认 `chr` 前缀与坐标体系（hg19/hg38）一致。
- **缺少 .tbi**：模块在运行前会校验每个 `*.tsv.2.gz` 的 `.tbi` 索引，缺失会直接报错。
- **PAR 处理**：默认将 `PAR` 注释映射到 `chrX` 片段；如需更细分（PAR1/PAR2），可在 `_get_vcf_chr_key_for_annotation()` 中扩展映射逻辑。

作者: ZHAO TIE
日期: 2025-10-12
版本: 1.0.0
"""

from pathlib import Path
from typing import Optional, Union
import logging


def _setup_logging(log_file: Optional[Union[str, Path]] = None, log_level: str = "INFO") -> logging.Logger:
    """设置日志系统
    
    Args:
        log_file: 日志文件路径，如果为None则使用默认路径
        log_level: 日志级别（DEBUG, INFO, WARNING, ERROR）
        
    Returns:
        配置好的logger实例
    """
    import sys
    import datetime
    
    # 获取日志文件路径
    if log_file is None:
        log_file = Path.cwd() / "snpeff_annotation.log"
    else:
        log_file = Path(log_file)
    
    # 确保日志目录存在
    log_file.parent.mkdir(parents=True, exist_ok=True)
    
    # 设置日志级别
    numeric_level = getattr(logging, log_level.upper(), logging.INFO)
    
    # 创建logger
    logger = logging.getLogger('snpeff_annotation')
    logger.setLevel(numeric_level)
    
    # 清除已有的handlers
    for handler in logger.handlers[:]:
        logger.removeHandler(handler)
    
    # 创建格式器
    formatter = logging.Formatter(
        '%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    
    # 文件处理器
    file_handler = logging.FileHandler(log_file, mode='a', encoding='utf-8')
    file_handler.setLevel(numeric_level)
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)
    
    # 控制台处理器（只显示INFO及以上级别）
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(logging.INFO)
    console_handler.setFormatter(logging.Formatter('%(levelname)s: %(message)s'))
    logger.addHandler(console_handler)
    
    return logger


def _get_log_file_path(log_file: Optional[Union[str, Path]] = None) -> Path:
    """获取日志文件路径"""
    if log_file is None:
        return Path.cwd() / "snpeff_annotation.log"
    return Path(log_file)

def add_chr_prefix_to_vcf(
    vcf_path: Union[str, Path],
    bcftools_path: Union[str, Path] = "/home/b/b37974/bcftools/bcftools",
    tabix_path: Union[str, Path] = "/home/b/b37974/htslib-1.9/tabix",
    out_path: Optional[Union[str, Path]] = None,
    force: bool = False,
) -> Path:
    """
    使用 **bcftools annotate --rename-chrs** 为 VCF 的 `#CHROM` 添加 "chr" 前缀，并用 **tabix** 建立索引。

    ➤ 设计要点
    - **仅使用 tabix** 建索引（不使用 `bcftools index`）。
    - 自动检测输入是否已带 `chr` 前缀；如已带且未指定 `out_path` 且 `force=False`，仅保证索引并返回原文件。
    - 全程使用 `pathlib.Path` 处理路径；函数返回 **`Path`**（--> Path）。

    Args:
        vcf_path: 输入的 BGZF 压缩 VCF 文件路径（.vcf.gz）。
        bcftools_path: `bcftools` 可执行文件的路径。
        tabix_path: `tabix` 可执行文件的路径。
        out_path: 可选的输出路径（默认将输入文件名替换为 `*.chrprefix.vcf.gz`）。
        force: 若为 True，即便输入已带 `chr` 也会生成 `out_path`（或默认输出）。

    Returns:
        Path: 输出的 `.vcf.gz` 文件路径（已建立 `.tbi` 索引）。

    Raises:
        FileNotFoundError: 当输入文件/工具路径不存在时。
        ValueError: 当输入文件不是以 `.vcf.gz` 结尾时。
        RuntimeError: 当外部命令执行失败时。

    Example:
        >>> add_chr_prefix_to_vcf(
        ...     "/data/sample.vcf.gz",
        ...     bcftools_path="/home/b/b37974/bcftools/bcftools",
        ...     tabix_path="/home/b/b37974/htslib-1.9/tabix",
        ... )
        PosixPath('/data/sample.chrprefix.vcf.gz')
    """
    import os
    import gzip
    import tempfile
    import subprocess

    logger = logging.getLogger(__name__)

    vcf_path = Path(vcf_path)
    bcftools_path = Path(bcftools_path)
    tabix_path = Path(tabix_path)

    if not vcf_path.exists():
        raise FileNotFoundError(f"VCF 不存在: {vcf_path}")
    if not str(vcf_path).endswith(".vcf.gz"):
        raise ValueError("输入文件必须是以 .vcf.gz 结尾的 BGZF VCF")
    if not bcftools_path.exists():
        raise FileNotFoundError(f"bcftools 不存在: {bcftools_path}")
    if not tabix_path.exists():
        raise FileNotFoundError(f"tabix 不存在: {tabix_path}")

    # 读取第一条变异行，判断是否已经带有 'chr' 前缀
    already_has_chr = False
    with gzip.open(vcf_path, "rt") as fin:
        for line in fin:
            if line.startswith("#"):
                continue
            chrom = line.split("\t", 1)[0]
            already_has_chr = chrom.startswith("chr")
            break

    # 默认输出路径：当前工作目录
    if out_path is None:
        out_path = Path.cwd() / vcf_path.name.replace(".vcf.gz", ".chrprefix.vcf.gz")
    out_path = Path(out_path)

    # 如果输入已带 chr 且不强制输出，则仅保证索引并返回输入路径
    if already_has_chr and not force and out_path == vcf_path:
        logger.info("输入已包含 'chr' 前缀，跳过重命名，仅确保索引存在: %s", vcf_path)
        if not vcf_path.with_suffix(vcf_path.suffix + ".tbi").exists() and not Path(str(vcf_path) + ".tbi").exists():
            try:
                subprocess.run([str(tabix_path), "-f", "-p", "vcf", str(vcf_path)], check=True)
            except subprocess.CalledProcessError as e:
                raise RuntimeError(f"tabix 索引失败: {e}")
        return vcf_path

    # 构建静态 rename-chrs 映射
    mapping_lines = [
        *(f"{i}\tchr{i}" for i in range(1, 23)),
        "X\tchrX",
        "Y\tchrY",
        "MT\tchrM",
        "M\tchrM",
    ]

    with tempfile.NamedTemporaryFile("w", delete=False, prefix="rename_chr_", suffix=".txt") as tf:
        tf.write("\n".join(mapping_lines) + "\n")
        rename_map_path = Path(tf.name)

    try:
        # bcftools annotate
        cmd_annotate = [
            str(bcftools_path),
            "annotate",
            "--rename-chrs",
            str(rename_map_path),
            "-Oz",
            "-o",
            str(out_path),
            str(vcf_path),
        ]
        logger.info("运行: %s", " ".join(cmd_annotate))
        subprocess.run(cmd_annotate, check=True)

        # 建立索引（仅 tabix）
        cmd_tabix = [str(tabix_path), "-f", "-p", "vcf", str(out_path)]
        logger.info("运行: %s", " ".join(cmd_tabix))
        subprocess.run(cmd_tabix, check=True)

    except subprocess.CalledProcessError as e:
        # 清理部分产物
        try:
            if out_path.exists():
                out_path.unlink()
            tbi = Path(str(out_path) + ".tbi")
            if tbi.exists():
                tbi.unlink()
        except Exception:
            pass
        raise RuntimeError(f"外部命令执行失败: {e}")

    finally:
        # 清理临时映射文件
        try:
            if rename_map_path.exists():
                rename_map_path.unlink()
        except Exception:
            pass

    return out_path



# Function to use tsv-formatted snpEff output(*.tsv.2.gz) & corresponding index(*.tbi) to annotate VCF file(vcf.gz)
# Note: snpEff output columns(#CHROM,POS,REF,ALT,id,qual,filter,allele,effect,impact,gene,geneid,feature,featureid,biotype,rank,hgvs_c,hgvs_p,cdna_pos,cdna_len,cds_pos,cds_len,aa_pos,aa_len,distance,errors,lof_gene,lof_geneid,lof_numtr,lof_perc,nmd_gene,nmd_geneid,nmd_numtr,nmd_perc,ID)
# Note: snpEff output 存在 #CHROM,POS,REF,ALT 完全相同的多条记录, 这些记录的effect,impact,gene等字段可能不同的情况。针对这一情况, 需要将这些记录的effect,impact,gene等字段进行合并, 用逗号分隔

def annotate_vcf_with_snpeff_tsv(
    vcf_path: Union[str, Path],
    snpeff_tsv_dir: Union[str, Path],
    *,
    bcftools_path: Union[str, Path] = "/home/b/b37974/bcftools/bcftools",
    tabix_path: Union[str, Path] = "/home/b/b37974/htslib-1.9/tabix",
    out_path: Optional[Union[str, Path]] = None,
    threads: int = 16,
    header_path: Optional[Union[str, Path]] = None,
    force: bool = False,
    remove_cache: bool = False,
    parallel: bool = True,
    max_workers: Optional[int] = None,
    log_file: Optional[Union[str, Path]] = None,
    log_level: str = "INFO",
) -> Path:
    """
    使用按染色体切分并带 .tbi 的 snpEff TSV 注释文件( *.tsv.2.gz )为输入 VCF(.vcf.gz) 填充 INFO。

    关键点：
      - 支持并行注释：按染色体分割VCF，并行注释各染色体，最后合并结果
      - 串行模式：逐个染色体依次注释（保持原有逻辑作为备选）
      - 将输入 VCF 放在 `bcftools annotate` 的命令最前（与你已验证成功的形式一致）
      - 使用一个统一的列映射文件 (-C) 来声明列到 INFO 的映射；不使用 --merge-logic
      - 如未提供 header_path，则在工作目录的 ./tmp 目录生成包含 10 个 INFO 标签的 header
      - 仅使用 tabix 建索引
      - 支持详细日志记录到文件

    Args:
        parallel: 是否使用并行模式（默认True）
        max_workers: 最大并行worker数量（默认为CPU核心数）
        log_file: 日志文件路径（默认为工作目录下的 snpeff_annotation.log）
        log_level: 日志级别（DEBUG, INFO, WARNING, ERROR，默认INFO）
    """
    import os
    import re
    import shutil
    import subprocess
    import tempfile
    import datetime
    from concurrent.futures import ThreadPoolExecutor, as_completed
    from multiprocessing import cpu_count

    # 设置日志系统
    logger = _setup_logging(log_file, log_level)
    
    # 记录开始时间和参数
    start_time = datetime.datetime.now()
    logger.info("=" * 80)
    logger.info("snpEff VCF 注释任务开始")
    logger.info(f"开始时间: {start_time.strftime('%Y-%m-%d %H:%M:%S')}")
    logger.info(f"输入VCF: {vcf_path}")
    logger.info(f"注释目录: {snpeff_tsv_dir}")
    logger.info(f"模式: {'并行' if parallel else '串行'}")
    if parallel:
        logger.info(f"最大worker数: {max_workers if max_workers else 'auto'}")
    logger.info(f"线程数: {threads}")
    logger.info(f"强制覆盖: {force}")
    logger.info(f"删除缓存: {remove_cache}")
    logger.info("=" * 80)

    vcf_path = Path(vcf_path)
    snpeff_tsv_dir = Path(snpeff_tsv_dir)
    bcftools_path = Path(bcftools_path)
    tabix_path = Path(tabix_path)

    # 基本校验
    try:
        if not vcf_path.exists():
            raise FileNotFoundError(f"VCF 不存在: {vcf_path}")
        if not str(vcf_path).endswith(".vcf.gz"):
            raise ValueError("输入文件必须是 .vcf.gz")
        if not snpeff_tsv_dir.is_dir():
            raise FileNotFoundError(f"snpEff 目录不存在: {snpeff_tsv_dir}")
        if not bcftools_path.exists():
            raise FileNotFoundError(f"bcftools 不存在: {bcftools_path}")
        if not tabix_path.exists():
            raise FileNotFoundError(f"tabix 不存在: {tabix_path}")
    except Exception as e:
        logger.error(f"输入验证失败: {e}")
        raise

    # 缓存目录 ./tmp
    cache_dir = Path.cwd() / "tmp"
    cache_dir.mkdir(parents=True, exist_ok=True)

    # 收集注释文件（*.tsv.2.gz）并排序（1..22, X, Y, M/MT, PAR）
    gz_files = sorted(snpeff_tsv_dir.glob("*.tsv.2.gz"))
    if not gz_files:
        raise FileNotFoundError(f"未找到 *.tsv.2.gz: {snpeff_tsv_dir}")

    def _chr_key(p: Path):
        # 内联实现染色体提取逻辑（与 _extract_chr_from_filename 保持一致）
        filename = p.name
        
        # 优先匹配 PAR（没有chr前缀）
        if "PAR" in filename:
            chr_name = "PAR"
        else:
            # 匹配 chr 开头的染色体编号
            import re
            match = re.search(r'chr(\w+)', filename)
            chr_name = match.group(1) if match else None
        
        if not chr_name:
            return (999, filename)
        
        # 定义排序优先级
        rank = {str(i): i for i in range(1, 23)}
        rank.update({"X": 23, "Y": 24, "M": 25, "MT": 25, "PAR": 26})
        return (rank.get(chr_name, 998), filename)

    gz_files = sorted(gz_files, key=_chr_key)
    logger.info(f"找到 {len(gz_files)} 个snpEff注释文件")
    for gf in gz_files:
        chr_name = _extract_chr_from_filename(gf)
        logger.info(f"  - {gf.name} (染色体: {chr_name})")

    # 确认每个注释文件都存在 .tbi
    missing_tbi = [p for p in gz_files if not (p.parent / (p.name + ".tbi")).exists()]
    if missing_tbi:
        raise FileNotFoundError("以下注释文件缺少 .tbi 索引:\n" + "\n".join(map(str, missing_tbi)))

    # 输出路径（默认 CWD）
    if out_path is None:
        out_path = Path.cwd() / vcf_path.name.replace(".vcf.gz", ".snpeff.vcf.gz")
    out_path = Path(out_path)
    if out_path.exists() and not force:
        raise FileExistsError(f"输出已存在: {out_path}（设置 force=True 以覆盖）")

    # 构建/确认 header（10 个 INFO 标签）
    info_tags = [
        ("effect",    "String", ".", "snpEff effect (merged)"),
        ("impact",    "String", ".", "snpEff impact (merged)"),
        ("gene",      "String", ".", "snpEff gene symbol (merged)"),
        ("geneid",    "String", ".", "snpEff gene id (merged)"),
        ("feature",   "String", ".", "snpEff feature (merged)"),
        ("featureid", "String", ".", "snpEff feature id (merged)"),
        ("biotype",   "String", ".", "snpEff biotype (merged)"),
        ("rank",      "String", ".", "snpEff rank (merged)"),
        ("hgvs_c",    "String", ".", "snpEff HGVS c. (merged)"),
        ("hgvs_p",    "String", ".", "snpEff HGVS p. (merged)"),
    ]

    if header_path is None:
        header_path = cache_dir / "snpeff_info.hdr"
        hdr_lines = [f"##INFO=<ID={t[0]},Number={t[2]},Type={t[1]},Description=\"{t[3]}\">" for t in info_tags]
        header_path.write_text("\n".join(hdr_lines) + "\n")
    else:
        header_path = Path(header_path)
        if not header_path.exists():
            raise FileNotFoundError(f"头文件不存在: {header_path}")

    # 统一的列映射 -C 文件（等价于 -c 的外部化）
    columns_file = cache_dir / "snpeff_columns.txt"
    columns_file.write_text(
        "CHROM,POS,REF,ALT,-,-,-,-,.INFO/effect,.INFO/impact,.INFO/gene,.INFO/geneid,.INFO/feature,.INFO/featureid,.INFO/biotype,.INFO/rank,.INFO/hgvs_c,.INFO/hgvs_p\n"
    )

    # 选择并行或串行模式
    try:
        if parallel:
            logger.info("使用并行模式注释VCF")
            result_path = _annotate_parallel(
                vcf_path, gz_files, cache_dir, header_path, columns_file,
                bcftools_path, tabix_path, out_path, threads, max_workers, logger
            )
        else:
            logger.info("使用串行模式注释VCF")
            result_path = _annotate_sequential(
                vcf_path, gz_files, cache_dir, header_path, columns_file,
                bcftools_path, tabix_path, out_path, threads, logger
            )
    except Exception as e:
        logger.error(f"注释过程失败: {e}")
        logger.error(f"错误类型: {type(e).__name__}")
        raise

    # Header 自检
    chk = subprocess.run([str(bcftools_path), "view", "-h", str(result_path)], capture_output=True, text=True)
    if chk.returncode != 0:
        raise RuntimeError("无法读取输出 VCF header：" + chk.stderr)
    hdr_text = chk.stdout
    missing = [t[0] for t in info_tags if f"##INFO=<ID={t[0]}," not in hdr_text]
    if missing:
        raise RuntimeError("输出 VCF 缺少以下 INFO 标签: " + ",".join(missing))

    logger.info("总计注释步数: %d", len(gz_files))
    for ann in gz_files:
        logger.info("使用注释: %s", ann)
    logger.info("使用 header: %s", header_path)
    logger.info("使用 columns(-C): %s", columns_file)
    logger.info("缓存目录: %s", cache_dir)

    # 根据参数决定是否清理缓存
    if remove_cache:
        try:
            shutil.rmtree(cache_dir)
            logger.info("已删除缓存目录: %s", cache_dir)
        except Exception as _e:
            logger.warning("删除缓存目录失败(%s): %s", type(_e).__name__, _e)
    else:
        logger.info("保留缓存目录以便调试: %s", cache_dir)

    # 记录完成信息
    end_time = datetime.datetime.now()
    duration = end_time - start_time
    logger.info("=" * 80)
    logger.info("snpEff VCF 注释任务完成")
    logger.info(f"结束时间: {end_time.strftime('%Y-%m-%d %H:%M:%S')}")
    logger.info(f"总耗时: {duration}")
    logger.info(f"输出文件: {result_path}")
    logger.info(f"输出文件大小: {result_path.stat().st_size / (1024*1024):.2f} MB")
    logger.info("=" * 80)

    print(f"当前工作目录: {os.getcwd()}")
    print(f"注释完成，输出文件: {result_path}")
    print(f"详细日志已写入: {_get_log_file_path(log_file)}")
    return result_path


def _annotate_parallel(
    vcf_path, gz_files, cache_dir, header_path, columns_file,
    bcftools_path, tabix_path, out_path, threads, max_workers, logger
):
    """并行注释模式：按染色体分割VCF，并行注释，最后合并"""
    from concurrent.futures import ThreadPoolExecutor, as_completed
    from multiprocessing import cpu_count
    import subprocess
    import shutil
    import datetime
    
    parallel_start = datetime.datetime.now()
    logger.info("开始并行注释模式")
    
    # 设置worker数量
    if max_workers is None:
        max_workers = min(len(gz_files), cpu_count())
    
    logger.info(f"使用 {max_workers} 个并行worker处理 {len(gz_files)} 个染色体")
    
    # 1. 按染色体分割输入VCF
    logger.info("步骤1: 按染色体分割VCF文件")
    split_start = datetime.datetime.now()
    chr_vcf_map = _split_vcf_by_chromosome(vcf_path, cache_dir, bcftools_path, logger)
    split_duration = datetime.datetime.now() - split_start
    logger.info(f"VCF分割完成，耗时: {split_duration}")
    
    # 2. 并行注释各染色体
    logger.info("步骤2: 并行注释各染色体")
    annotate_start = datetime.datetime.now()
    annotated_files = []
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        # 提交所有注释任务
        future_to_chr = {}
        for ann_file in gz_files:
            # 从注释文件名提取染色体
            chr_name = _extract_chr_from_filename(ann_file)
            
            # 获取对应的VCF片段key（处理PAR等特殊情况）
            vcf_chr_key = _get_vcf_chr_key_for_annotation(chr_name)
            
            if vcf_chr_key in chr_vcf_map:
                chr_vcf = chr_vcf_map[vcf_chr_key]
                future = executor.submit(
                    _annotate_single_chromosome,
                    chr_vcf, ann_file, cache_dir, header_path, columns_file,
                    bcftools_path, threads, chr_name, logger
                )
                future_to_chr[future] = chr_name
                logger.info(f"已提交任务: 注释文件 {ann_file.name} (chr={chr_name}) -> VCF片段 {vcf_chr_key}")
            else:
                logger.info(f"跳过注释文件 {ann_file.name}：VCF中不存在对应染色体 {vcf_chr_key}")
        
        # 收集结果
        completed_count = 0
        total_tasks = len(future_to_chr)
        for future in as_completed(future_to_chr):
            chr_name = future_to_chr[future]
            try:
                annotated_file = future.result()
                if annotated_file.exists():
                    annotated_files.append(annotated_file)
                    completed_count += 1
                    logger.info(f"染色体 {chr_name} 注释完成 ({completed_count}/{total_tasks}): {annotated_file}")
            except Exception as e:
                logger.error(f"染色体 {chr_name} 注释失败: {e}")
                raise
    
    annotate_duration = datetime.datetime.now() - annotate_start
    logger.info(f"并行注释完成，耗时: {annotate_duration}")
    
    # 3. 合并所有注释文件
    logger.info("步骤3: 合并注释文件")
    merge_start = datetime.datetime.now()
    merged_file = _merge_annotated_vcfs(annotated_files, cache_dir, bcftools_path, logger)
    merge_duration = datetime.datetime.now() - merge_start
    logger.info(f"文件合并完成，耗时: {merge_duration}")
    
    # 4. 复制到最终输出位置并建索引
    logger.info("步骤4: 生成最终输出文件和索引")
    final_start = datetime.datetime.now()
    if merged_file != out_path:
        if out_path.exists():
            out_path.unlink()
        shutil.copy2(merged_file, out_path)
    
    # 建立索引
    proc = subprocess.run([str(tabix_path), "-f", "-p", "vcf", str(out_path)], capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(f"tabix 索引失败\nSTDERR:\n{proc.stderr.strip()}")
    
    final_duration = datetime.datetime.now() - final_start
    parallel_total_duration = datetime.datetime.now() - parallel_start
    logger.info(f"最终文件生成完成，耗时: {final_duration}")
    logger.info(f"并行模式总耗时: {parallel_total_duration}")
    
    return out_path


def _annotate_sequential(
    vcf_path, gz_files, cache_dir, header_path, columns_file,
    bcftools_path, tabix_path, out_path, threads, logger
):
    """串行注释模式：原有的逐个染色体依次注释逻辑"""
    import subprocess
    import shutil
    import datetime
    
    sequential_start = datetime.datetime.now()
    logger.info("开始串行注释模式")
    logger.info(f"将依次处理 {len(gz_files)} 个染色体注释文件")
    
    # 工作链：逐个注释，前一次输出作为下一次输入
    current_in = vcf_path
    temp_products = []
    try:
        for i, ann in enumerate(gz_files, start=1):
            step_start = datetime.datetime.now()
            chr_name = _extract_chr_from_filename(ann)
            logger.info(f"开始步骤 {i}/{len(gz_files)}: 注释染色体 {chr_name}")
            
            # 每个步骤输出到临时文件
            tmp_out = cache_dir / f"step_{i:02d}.vcf.gz"
            cmd = [str(bcftools_path), "annotate"]
            # 仅第一次注入 header，避免重复追加；将 -h 紧跟在子命令后，避免打断选项与其参数
            if i == 1:
                cmd += ["-h", str(header_path)]
            cmd += [
                str(current_in),
                "--threads", str(threads),
                "-a", str(ann),
                "-C", str(columns_file),
                "-Oz", "-o", str(tmp_out),
            ]

            logger.info("运行 annotate (%d/%d): %s", i, len(gz_files), " ".join(cmd))
            proc = subprocess.run(cmd, capture_output=True, text=True)
            if proc.returncode != 0:
                raise RuntimeError(
                    "bcftools annotate 失败\n"
                    f"STEP: {i}/{len(gz_files)}\n"
                    f"ANN: {ann}\n"
                    f"CMD: {' '.join(cmd)}\n"
                    f"STDERR:\n{proc.stderr.strip()}\nSTDOUT:\n{proc.stdout.strip()}"
                )

            step_duration = datetime.datetime.now() - step_start
            logger.info(f"步骤 {i}/{len(gz_files)} 完成，染色体 {chr_name}，耗时: {step_duration}")

            # 下游输入切换为本步输出
            current_in = tmp_out
            temp_products.append(tmp_out)

        # 全部步骤完成后，把最后产物移动/覆盖到 out_path
        logger.info("所有串行注释步骤完成，生成最终输出文件")
        if current_in != out_path:
            if out_path.exists():
                out_path.unlink()
            shutil.copy2(current_in, out_path)

        # 仅用 tabix 建索引
        logger.info("建立最终文件索引")
        proc2 = subprocess.run([str(tabix_path), "-f", "-p", "vcf", str(out_path)], capture_output=True, text=True)
        if proc2.returncode != 0:
            raise RuntimeError(
                "tabix 索引失败\n" +
                f"CMD: {tabix_path} -f -p vcf {out_path}\n" +
                f"STDERR:\n{proc2.stderr.strip()}\nSTDOUT:\n{proc2.stdout.strip()}"
            )

        sequential_total_duration = datetime.datetime.now() - sequential_start
        logger.info(f"串行模式总耗时: {sequential_total_duration}")

    except Exception:
        # 失败时尝试清理目标产物（保留中间文件用于调试）
        logger.error("串行注释过程中发生错误，正在清理...")
        try:
            if out_path.exists():
                out_path.unlink()
            tbi = Path(str(out_path) + ".tbi")
            if tbi.exists():
                tbi.unlink()
        except Exception:
            pass
        raise
    
    return out_path


def _split_vcf_by_chromosome(vcf_path, cache_dir, bcftools_path, logger):
    """按染色体分割VCF文件"""
    import subprocess
    
    logger.info(f"开始按染色体分割VCF文件...")
    
    # 获取VCF中存在的染色体列表
    cmd_chroms = [str(bcftools_path), "index", "-s", str(vcf_path)]
    proc = subprocess.run(cmd_chroms, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(f"获取染色体列表失败: {proc.stderr}")
    
    chr_vcf_map = {}
    available_chroms = [line.split('\t')[0] for line in proc.stdout.strip().split('\n') if line.strip()]
    logger.info(f"VCF中检测到的染色体: {available_chroms}")
    
    for chrom in available_chroms:
        # 标准化染色体名称用作key
        chr_key = chrom.replace('chr', '')
        chr_vcf = cache_dir / f"input_chr{chr_key}.vcf.gz"
        
        # 分割染色体
        cmd_split = [
            str(bcftools_path), "view",
            "-r", chrom,
            "-Oz", "-o", str(chr_vcf),
            str(vcf_path)
        ]
        
        proc_split = subprocess.run(cmd_split, capture_output=True, text=True)
        if proc_split.returncode != 0:
            logger.warning(f"分割染色体 {chrom} 失败: {proc_split.stderr}")
            continue
        
        # 检查文件是否为空
        if chr_vcf.stat().st_size > 100:  # 至少有header
            chr_vcf_map[chr_key] = chr_vcf
            logger.info(f"分割染色体 {chrom} 完成: {chr_vcf} (key: {chr_key})")
        else:
            chr_vcf.unlink()  # 删除空文件
            logger.info(f"删除空文件: {chr_vcf}")
    
    logger.info(f"成功分割 {len(chr_vcf_map)} 个染色体，可用keys: {list(chr_vcf_map.keys())}")
    return chr_vcf_map


def _extract_chr_from_filename(filepath):
    """从文件名提取染色体编号
    
    支持多种命名格式：
    - all.VQSR3.chr1.vcf_out.tsv.2.gz -> "1"
    - all.VQSR3.chrX.vcf_out.tsv.2.gz -> "X"
    - all.VQSR3.chrM.vcf_out.tsv.2.gz -> "M"
    - all.VQSR3.PAR.vcf_out.tsv.2.gz -> "PAR"
    - chr1.tsv.2.gz -> "1" (向后兼容)
    """
    import re
    
    # 优先匹配 PAR（没有chr前缀）
    if "PAR" in filepath.name:
        return "PAR"
    
    # 匹配 chr 开头的染色体编号
    match = re.search(r'chr(\w+)', filepath.name)
    if match:
        return match.group(1)
    
    return None


def _get_vcf_chr_key_for_annotation(annotation_chr):
    """获取注释文件对应的VCF染色体key
    
    处理特殊情况：
    - PAR注释文件对应chrX的VCF片段
    - 其他情况直接对应
    
    Args:
        annotation_chr: 从注释文件名提取的染色体（如"1", "X", "PAR", "MT"）
        
    Returns:
        str: VCF片段的key（如"1", "X", "MT"）
    """
    if annotation_chr == "PAR":
        return "X"  # PAR注释使用chrX的VCF片段
    return annotation_chr


def _annotate_single_chromosome(
    chr_vcf, ann_file, cache_dir, header_path, columns_file,
    bcftools_path, threads, chr_name, logger
):
    """注释单个染色体"""
    import subprocess
    import datetime
    
    start_time = datetime.datetime.now()
    output_file = cache_dir / f"annotated_chr{chr_name}.vcf.gz"
    
    # 记录输入文件信息
    input_size = chr_vcf.stat().st_size / (1024*1024)  # MB
    ann_size = ann_file.stat().st_size / (1024*1024)   # MB
    logger.debug(f"染色体 {chr_name}: 输入VCF大小 {input_size:.2f}MB, 注释文件大小 {ann_size:.2f}MB")
    
    cmd = [str(bcftools_path), "annotate"]
    cmd += ["-h", str(header_path)]  # 每个文件都需要header
    cmd += [
        str(chr_vcf),
        "--threads", str(threads),
        "-a", str(ann_file),
        "-C", str(columns_file),
        "-Oz", "-o", str(output_file),
    ]
    
    logger.debug(f"注释染色体 {chr_name}: {' '.join(cmd)}")
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(
            f"染色体 {chr_name} 注释失败\n"
            f"CMD: {' '.join(cmd)}\n"
            f"STDERR:\n{proc.stderr.strip()}"
        )
    
    # 记录输出信息
    duration = datetime.datetime.now() - start_time
    if output_file.exists():
        output_size = output_file.stat().st_size / (1024*1024)  # MB
        logger.debug(f"染色体 {chr_name}: 注释完成，输出大小 {output_size:.2f}MB，耗时 {duration}")
    
    return output_file


def _merge_annotated_vcfs(annotated_files, cache_dir, bcftools_path, logger):
    """合并注释后的VCF文件"""
    import subprocess
    import re
    
    if not annotated_files:
        raise RuntimeError("没有注释文件可合并")
    
    if len(annotated_files) == 1:
        return annotated_files[0]
    
    logger.info(f"开始合并 {len(annotated_files)} 个注释文件...")
    
    # 按染色体顺序排序
    def _get_chr_order(filepath):
        match = re.search(r'chr(\w+)', filepath.name)
        if not match:
            return (999, filepath.name)
        chr_token = match.group(1)
        
        # 定义排序优先级
        order_map = {str(i): i for i in range(1, 23)}
        order_map.update({"X": 23, "Y": 24, "M": 25, "MT": 25, "PAR": 26})
        return (order_map.get(chr_token, 998), filepath.name)
    
    sorted_files = sorted(annotated_files, key=_get_chr_order)
    
    # 使用bcftools concat合并
    merged_file = cache_dir / "merged_annotated.vcf.gz"
    cmd_concat = [
        str(bcftools_path), "concat",
        "-Oz", "-o", str(merged_file)
    ] + [str(f) for f in sorted_files]
    
    logger.info(f"合并命令: {' '.join(cmd_concat)}")
    proc = subprocess.run(cmd_concat, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(f"合并VCF失败: {proc.stderr}")
    
    logger.info(f"合并完成: {merged_file}")
    return merged_file