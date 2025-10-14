#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
VCF INFO 字段过滤工具
===================

本模块提供了基于 bcftools 的高效 VCF 文件 INFO 字段过滤功能，专门针对大型基因组数据设计。

主要功能
--------
1. **INFO 字段过滤**: 根据指定的 INFO 字段值过滤 VCF 变异位点
2. **染色体前缀标准化**: 自动检测并统一 'chr' 前缀格式
3. **多种匹配模式**: 支持精确匹配、包含匹配和正则匹配
4. **并行处理**: 利用多线程加速压缩和处理
5. **自动索引**: 输出文件自动创建 tabix 索引

使用场景
--------
- 过滤高影响或中等影响的变异 (SnpEff 注释结果)
- 根据功能注释筛选特定类型变异
- 批量处理大型 WGS/WES 数据集
- 为关联性分析准备特定子集数据

性能特点
--------
- 直接处理 bgzip 压缩的 VCF 文件，无需解压
- 使用 bcftools 原生过滤，内存占用小
- 支持多线程并行压缩
- 原子操作保证文件完整性

依赖工具
--------
- bcftools: VCF 文件操作和过滤
- tabix: VCF 文件索引
- bgzip: 文件压缩

典型用法示例
-----------
```python
from info_filter_tools import filter_vcf_by_info

# 过滤高影响和中等影响的变异
output_path = filter_vcf_by_info(
    vcf_path="input.vcf.gz",
    info_key="impact",
    values=["HIGH", "MODERATE"],
    out_prefix="filtered_variants"
)

# 使用正则表达式匹配
output_path = filter_vcf_by_info(
    vcf_path="input.vcf.gz",
    info_key="consequence",
    values="missense.*",
    match_mode="regex"
)
```

作者: ZHAO TIE
日期: 2025-10-14
"""

import os
import shlex
import subprocess
import re
import difflib
from typing import Iterable, List, Optional, Union
import logging
from datetime import datetime
import tempfile


# Helper: collect INFO IDs from VCF header using bcftools
def _collect_info_ids_from_header(vcf_path: str, bcftools_path: str) -> set:
    """使用 bcftools 读取 VCF 头部，提取所有 INFO 的 ID 集合。
    若 bcftools 调用失败，将抛出异常。
    """
    proc = subprocess.run(
        [bcftools_path, "view", "-h", vcf_path],
        check=True,
        capture_output=True,
        text=True,
    )
    info_ids = set()
    for line in proc.stdout.splitlines():
        # 形如：##INFO=<ID=impact,Number=.,Type=String,Description="...">
        if line.startswith("##INFO=<ID="):
            # 提取 ID= 和 后续逗号 之间的部分
            try:
                id_part = line.split("##INFO=<ID=", 1)[1]
                field_id = id_part.split(",", 1)[0].strip()
                if field_id:
                    info_ids.add(field_id)
            except Exception:
                pass
    return info_ids


def _detect_chr_prefix(vcf_path: str, bcftools_path: str, sample_lines: int = 200) -> Optional[bool]:
    """检测 VCF 的 CHROM 是否带有 'chr' 前缀。
    返回值: True 表示带有前缀；False 表示不带前缀；None 表示无法判断（例如没有变体行）。
    说明：为了避免读取全文件，这里仅采样前若干变体行进行判断。
    """
    try:
        # 输出体（非头部）若干行进行判断
        # -H 仅输出变体行；这里不加区域限制，直接采样前 sample_lines 行
        cmd = [bcftools_path, "view", "-H", vcf_path]
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True)
        has_chr = None
        cnt = 0
        if proc.stdout is None:
            return None
        for line in proc.stdout:
            if not line.strip():
                continue
            chrom = line.split("\t", 1)[0]
            if chrom.startswith("#"):
                # 保险判断：如果有异常行
                continue
            has_chr = chrom.startswith("chr")
            cnt += 1
            if cnt >= sample_lines:
                break
        proc.stdout.close()
        proc.wait()
        return has_chr
    except Exception:
        return None


def filter_vcf_by_info(
    vcf_path: str = "cteph_agp3k.rare.all.nochr.norm.chrprefix.snpeff.vcf.gz",
    *,
    bcftools_path: str = "/home/b/b37974/bcftools/bcftools",
    tabix_path: str = "/home/b/b37974/htslib-1.9/tabix",
    info_key: str = "INFO/impact",
    values: Union[str, Iterable[str]] = ("HIGH", "MODERATE"),
    logic: str = "any",
    match_mode: str = "exact",
    check_chr_prefix: bool = False,
    keep_chr_prefix: bool = False,
    threads: int = 8,
    out_prefix: Optional[str] = None,
    out_dir: Optional[str] = os.getcwd(),
    index_output: bool = True,
    index_type: str = "tbi",
    overwrite: bool = True,
    dry_run: bool = False,
) -> str:
    """根据 INFO 字段过滤超大 VCF.GZ（基于 bcftools），并返回生成的 VCF.GZ 路径。

    参数
    ----
    vcf_path : str
        输入的 .vcf.gz 文件路径（必须已 bgzip 压缩），默认为项目罕见变异注释文件。
    bcftools_path : str
        bcftools 可执行程序路径。
    tabix_path : str
        tabix 可执行程序路径，用于输出索引。
    info_key : str
        需要筛选的 INFO 键。可写成 "impact" 或 "INFO/impact"（函数会自动补全）。
    values : str | Iterable[str]
        需要保留的取值（单值或多值）。例如 "HIGH" 或 ["HIGH", "MODERATE"].
    logic : {"any", "all"}
        多值时的逻辑，"any" 表示 OR，"all" 表示 AND（极少用）。
    match_mode : {"exact", "contains", "regex"}
        匹配方式。exact 使用等号匹配；contains 使用正则子串匹配；regex 完全按传入值当作正则。
    check_chr_prefix : bool
        是否检查并根据需要统一 CHROM 的 'chr' 前缀（默认 False）。为 True 时先探测输入的染色体命名风格，然后视 keep_chr_prefix 决定是否用 bcftools annotate --rename-chrs 进行标准化。
    keep_chr_prefix : bool
        仅当 check_chr_prefix=True 生效。若输入带有 'chr' 前缀且 keep_chr_prefix=False，则去掉前缀（chr1→1, chrX→X, chrY→Y, chrM→M）。若输入不带 'chr' 前缀且 keep_chr_prefix=True，则添加前缀（1→chr1, X→chrX, Y→chrY, M→chrM）。其他情况不做修改。
    threads : int
        压缩线程数（传给 bcftools 的 --threads，用于并行 bgzip 压缩）。
    out_prefix : Optional[str]
        输出文件前缀名（不含扩展名）。可以是绝对路径或仅文件名；若为空，则自动基于输入和条件生成。
    out_dir : Optional[str]
        输出目录；若 out_prefix 为相对/仅文件名且指定了 out_dir，则输出到该目录。
    index_output : bool
        是否对输出 .vcf.gz 执行 tabix 索引（默认 True）。
    index_type : {"tbi", "csi"}
        索引类型。默认 tbi；若染色体坐标非常大可选 csi。
    overwrite : bool
        若目标已存在，是否覆盖。
    dry_run : bool
        若 True，仅返回将要执行的命令并不实际运行（用于调试）。

    返回
    ----
    str
        过滤后的 .vcf.gz 文件路径。

    说明与加速要点
    --------------
    1) 仅使用 INFO 表达式过滤（-i），无需解压；--threads 用于并行压缩加速。
    2) 若上游已按染色体拆分，可在调用者层面并行多个染色体以进一步加速。
    3) 输出使用原子写入（先到 tmp，再替换），随后默认 tabix 建索引，方便下游随机访问。
    """

    # ---------- 参数与路径校验 ----------
    if not os.path.exists(vcf_path):
        raise FileNotFoundError(f"VCF 文件不存在: {vcf_path}")
    if not os.path.exists(bcftools_path):
        raise FileNotFoundError(f"bcftools 不存在: {bcftools_path}")
    if index_output and not os.path.exists(tabix_path):
        raise FileNotFoundError(f"tabix 不存在: {tabix_path}")

    # ---------- 日志初始化 ----------
    log_dir = os.getcwd()
    log_file = os.path.join(log_dir, f"filter_vcf_by_info_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log")

    # 使用独立 logger，避免被外部 basicConfig 覆盖
    logger = logging.getLogger("info_filter_tools.filter_vcf_by_info")
    logger.setLevel(logging.DEBUG)
    # 避免重复添加 handler
    for h in list(logger.handlers):
        logger.removeHandler(h)
    fh = logging.FileHandler(log_file, encoding='utf-8')
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(logging.Formatter('%(asctime)s [%(levelname)s] %(message)s'))
    logger.addHandler(fh)
    logger.propagate = False

    logger.info("=== 开始执行 filter_vcf_by_info ===")
    logger.info(f"输入文件: {vcf_path}")
    logger.info(f"输出前缀: {out_prefix}")
    logger.info(f"INFO 键: {info_key}")
    logger.info(f"取值: {values}")
    logger.info(f"逻辑: {logic}, 匹配模式: {match_mode}")
    logger.info(f"线程数: {threads}")
    logger.info(f"check_chr_prefix: {check_chr_prefix}, keep_chr_prefix: {keep_chr_prefix}")

    # ---------- （可选）统一 CHROM 前缀 ----------
    normalized_vcf_path = vcf_path
    tmp_map_file = None
    tmp_norm_vcf = None

    if check_chr_prefix:
        detected = _detect_chr_prefix(vcf_path, bcftools_path)
        logger.info(f"CHROM 前缀探测结果: {detected} (True=带有 'chr', False=不带, None=无法判断)")
        action = None
        if detected is True and not keep_chr_prefix:
            action = "remove_chr_prefix"
        elif detected is False and keep_chr_prefix:
            action = "add_chr_prefix"

        if action is None:
            logger.info("不需要修改 CHROM 前缀，继续后续过滤。")
        else:
            # 构造 rename-chrs 映射文件
            mapping_lines = []
            if action == "remove_chr_prefix":
                for i in range(1, 23):
                    mapping_lines.append(f"chr{i}\t{i}\n")
                mapping_lines.extend([
                    "chrX\tX\n",
                    "chrY\tY\n",
                    "chrM\tM\n",
                ])
            elif action == "add_chr_prefix":
                for i in range(1, 23):
                    mapping_lines.append(f"{i}\tchr{i}\n")
                mapping_lines.extend([
                    "X\tchrX\n",
                    "Y\tchrY\n",
                    "M\tchrM\n",
                ])

            # 写入临时映射文件
            tf = tempfile.NamedTemporaryFile("w", delete=False, prefix="rename_chrs_", suffix=".tsv")
            tf.writelines(mapping_lines)
            tf.flush()
            tf.close()
            tmp_map_file = tf.name
            logger.info(f"生成 rename-chrs 映射文件: {tmp_map_file}，共 {len(mapping_lines)} 行；action={action}")

            # 生成临时标准化 VCF（bgzip 压缩）
            in_dir = os.path.dirname(os.path.abspath(vcf_path))
            tmp_norm_vcf = os.path.join(in_dir, f"chrfix.{os.path.basename(vcf_path)}")
            ann_cmd = [
                bcftools_path,
                "annotate",
                "--rename-chrs",
                tmp_map_file,
                "-Oz",
                "-o",
                tmp_norm_vcf,
                "--threads",
                str(int(threads) if threads and threads > 0 else 1),
                vcf_path,
            ]
            logger.info("bcftools（标准化 CHROM）命令: " + " ".join(ann_cmd))

            if dry_run:
                # 仅在 dry-run 时把命令记下来并返回给调用者（稍后会与过滤命令一起返回）
                chrfix_cmd_str = " ".join(shlex.quote(c) for c in ann_cmd)
            else:
                try:
                    if os.path.exists(tmp_norm_vcf):
                        os.remove(tmp_norm_vcf)
                    subprocess.run(ann_cmd, check=True)
                    normalized_vcf_path = tmp_norm_vcf
                    logger.info(f"已完成 CHROM 标准化: {normalized_vcf_path}")
                except subprocess.CalledProcessError as e:
                    logger.error(f"CHROM 标准化失败: {e}. 命令: {' '.join(ann_cmd)}")
                    # 清理半成品
                    if os.path.exists(tmp_norm_vcf):
                        try:
                            os.remove(tmp_norm_vcf)
                        except Exception:
                            pass
                    raise RuntimeError(
                        f"bcftools annotate --rename-chrs 失败 (退出码={e.returncode}). 命令: {' '.join(shlex.quote(x) for x in ann_cmd)}"
                    ) from e

    # 规范化 info_key
    if not info_key.startswith("INFO/"):
        info_key = f"INFO/{info_key}"

    # 校验该 INFO key 是否存在于头部；若不存在给出相近建议
    key_name = info_key.split("/", 1)[1]
    try:
        info_ids = _collect_info_ids_from_header(normalized_vcf_path, bcftools_path)
    except subprocess.CalledProcessError as e:
        raise RuntimeError(
            f"读取 VCF 头部失败，无法校验 INFO 字段。命令退出码={e.returncode}"
        ) from e

    if key_name not in info_ids:
        suggestions = difflib.get_close_matches(key_name, sorted(info_ids), n=5, cutoff=0.6)
        hint = ("；你是否想用: " + ", ".join(suggestions)) if suggestions else ""
        raise ValueError(
            f"INFO 字段 '{key_name}' 不存在于 VCF 头部{hint}。可用字段示例: "
            + ", ".join(list(sorted(info_ids))[:10])
            + (" …" if len(info_ids) > 10 else "")
        )

    # 规范化 values 为列表
    if isinstance(values, (str, bytes)):
        values_list: List[str] = [str(values)]
    else:
        values_list = [str(v) for v in values]
    if len(values_list) == 0:
        raise ValueError("values 至少应包含一个取值")

    # ---------- 构建过滤表达式 ----------
    # bcftools 表达式示例： INFO/impact=="HIGH" || INFO/impact=="MODERATE"
    ops = {
        "exact": lambda v: f'{info_key}=="{v}"',
        "contains": lambda v: f'{info_key} ~ "{re.escape(v)}"',
        "regex": lambda v: f'{info_key} ~ "{v}"',
    }
    if match_mode not in ops:
        raise ValueError("match_mode 只能为 'exact' | 'contains' | 'regex'")

    clauses = [ops[match_mode](v) for v in values_list]
    joiner = " || " if logic == "any" else " && "
    expr = joiner.join(clauses)

    logger.info(f"构建的过滤表达式: {expr}")

    # ---------- 输出路径 ----------
    in_dir = os.path.dirname(os.path.abspath(vcf_path))
    auto_name = (
        f"infofilter_{info_key.split('/',1)[1]}_" + "-".join(values_list)
    )
    base_prefix = out_prefix if out_prefix else auto_name

    # 如果 out_prefix 没有目录，并且给了 out_dir，则输出到 out_dir
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
        if os.path.basename(base_prefix) == base_prefix:
            base_prefix = os.path.join(out_dir, base_prefix)
    else:
        # 相对/仅文件名：默认写到输入文件所在目录
        if os.path.basename(base_prefix) == base_prefix:
            base_prefix = os.path.join(in_dir, base_prefix)

    out_vcf = f"{base_prefix}.vcf.gz"
    tmp_vcf = f"{out_vcf}.tmp"

    if os.path.exists(out_vcf) and not overwrite:
        raise FileExistsError(f"输出已存在且不允许覆盖: {out_vcf}")

    # ---------- 组装命令 ----------
    # 注：--threads 仅对压缩阶段加速；表达式作为一个参数传入。
    cmd = [
        bcftools_path,
        "view",
        "-i",
        expr,
        "-Oz",
        "-o",
        tmp_vcf,
        "--threads",
        str(int(threads) if threads and threads > 0 else 1),
        normalized_vcf_path,
    ]

    logger.info(f"bcftools 命令: {' '.join(cmd)}")

    # ---------- dry-run ----------
    if dry_run:
        filter_cmd_str = " ".join(shlex.quote(c) for c in cmd)
        if check_chr_prefix and 'chrfix_cmd_str' in locals():
            return chrfix_cmd_str + " && " + filter_cmd_str # type: ignore[return-value]
        return filter_cmd_str

    # ---------- 执行过滤 ----------
    try:
        # 先删除旧的 tmp
        if os.path.exists(tmp_vcf):
            os.remove(tmp_vcf)
        subprocess.run(cmd, check=True)
        # 原子替换
        os.replace(tmp_vcf, out_vcf)
        logger.info(f"过滤完成，生成文件: {out_vcf}")
    except subprocess.CalledProcessError as e:
        logger.error(f"执行失败: {e}. 命令: {' '.join(cmd)}")
        # 若失败，清理 tmp
        if os.path.exists(tmp_vcf):
            try:
                os.remove(tmp_vcf)
            except Exception:
                pass
        raise RuntimeError(
            f"bcftools 过滤失败 (退出码={e.returncode}). 命令: {' '.join(shlex.quote(x) for x in cmd)}"
        ) from e

    # ---------- 索引输出 ----------
    if index_output:
        idx_cmd = [tabix_path, "-f", "-p", "vcf"]
        if index_type.lower() == "csi":
            idx_cmd.append("-C")  # 生成 CSI 索引
        idx_cmd.append(out_vcf)
        try:
            subprocess.run(idx_cmd, check=True)
            logger.info(f"索引完成: {out_vcf}.{index_type}")
            logger.info("=== 执行完成 ===")
        except subprocess.CalledProcessError as e:
            raise RuntimeError(
                f"tabix 索引失败 (退出码={e.returncode}). 命令: {' '.join(shlex.quote(x) for x in idx_cmd)}"
            ) from e

    # ---------- 清理临时文件 ----------
    try:
        if tmp_map_file and os.path.exists(tmp_map_file):
            os.remove(tmp_map_file)
        # 若生成了标准化的临时 VCF，使用后删除
        if tmp_norm_vcf and os.path.exists(tmp_norm_vcf):
            try:
                os.remove(tmp_norm_vcf)
                if os.path.exists(tmp_norm_vcf + ".tbi"):
                    os.remove(tmp_norm_vcf + ".tbi")
            except Exception:
                pass
    except Exception:
        # 清理失败不应影响主流程
        logger.warning("清理临时文件时发生非致命错误", exc_info=True)

    return out_vcf
