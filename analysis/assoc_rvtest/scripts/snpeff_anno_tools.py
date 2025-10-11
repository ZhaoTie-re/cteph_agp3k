from pathlib import Path
from typing import Optional, Union
import logging

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



# New function to use tsv-formatted snpEff output(*.tsv.2.gz) & corresponding index(*.tbi) to annotate VCF file(vcf.gz)
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
) -> Path:
    """
    使用按染色体切分并带 .tbi 的 snpEff TSV 注释文件( *.tsv.2.gz )为输入 VCF(.vcf.gz) 填充 INFO。

    仅修改本函数，其他代码保持不变。

    关键点：
      - 逐个染色体依次注释：每次以上一步产物为输入，避免一次性堆叠导致的 header 同步问题。
      - 将输入 VCF 放在 `bcftools annotate` 的命令最前（与你已验证成功的形式一致）。
      - 使用一个统一的列映射文件 (-C) 来声明列到 INFO 的映射；不使用 --merge-logic。
      - 如未提供 header_path，则在工作目录的 ./tmp 目录生成包含 10 个 INFO 标签的 header。
      - 仅使用 tabix 建索引。
    """
    import os
    import re
    import shutil
    import subprocess
    import tempfile

    logger = logging.getLogger(__name__)

    vcf_path = Path(vcf_path)
    snpeff_tsv_dir = Path(snpeff_tsv_dir)
    bcftools_path = Path(bcftools_path)
    tabix_path = Path(tabix_path)

    # 基本校验
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

    # 缓存目录 ./tmp
    cache_dir = Path.cwd() / "tmp"
    cache_dir.mkdir(parents=True, exist_ok=True)

    # 收集注释文件（*.tsv.2.gz）并排序（1..22, X, Y, M/MT, PAR）
    gz_files = sorted(snpeff_tsv_dir.glob("*.tsv.2.gz"))
    if not gz_files:
        raise FileNotFoundError(f"未找到 *.tsv.2.gz: {snpeff_tsv_dir}")

    def _chr_key(p: Path):
        m = re.search(r"chr(\w+)", p.name)
        if not m:
            return (999, p.name)
        tok = m.group(1)
        rank = {str(i): i for i in range(1, 23)}
        rank.update({"X": 23, "Y": 24, "M": 25, "MT": 25, "PAR": 26})
        return (rank.get(tok, 998), p.name)

    gz_files = sorted(gz_files, key=_chr_key)

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

    # 工作链：逐个注释，前一次输出作为下一次输入
    current_in = vcf_path
    temp_products = []
    try:
        for i, ann in enumerate(gz_files, start=1):
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

            # 下游输入切换为本步输出
            current_in = tmp_out
            temp_products.append(tmp_out)

        # 全部步骤完成后，把最后产物移动/覆盖到 out_path
        if current_in != out_path:
            if out_path.exists():
                out_path.unlink()
            shutil.copy2(current_in, out_path)

        # 仅用 tabix 建索引
        proc2 = subprocess.run([str(tabix_path), "-f", "-p", "vcf", str(out_path)], capture_output=True, text=True)
        if proc2.returncode != 0:
            raise RuntimeError(
                "tabix 索引失败\n" +
                f"CMD: {tabix_path} -f -p vcf {out_path}\n" +
                f"STDERR:\n{proc2.stderr.strip()}\nSTDOUT:\n{proc2.stdout.strip()}"
            )

        # Header 自检
        chk = subprocess.run([str(bcftools_path), "view", "-h", str(out_path)], capture_output=True, text=True)
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

    except Exception:
        # 失败时尝试清理目标产物（保留中间文件用于调试）
        try:
            if out_path.exists():
                out_path.unlink()
            tbi = Path(str(out_path) + ".tbi")
            if tbi.exists():
                tbi.unlink()
        except Exception:
            pass
        raise
    finally:
        # 根据参数决定是否清理缓存
        if remove_cache:
            try:
                shutil.rmtree(cache_dir)
                logger.info("已删除缓存目录: %s", cache_dir)
            except Exception as _e:
                logger.warning("删除缓存目录失败(%s): %s", type(_e).__name__, _e)
        else:
            logger.info("保留缓存目录以便调试: %s", cache_dir)

    print(f"当前工作目录: {os.getcwd()}")
    return out_path