# -*- coding: utf-8 -*-
"""
功能：
  - 使用 plink2 将 .vcf.gz 转换为 PLINK 二进制（.bed/.bim/.fam），并确保 FID 与 IID 都等于样本 ID（--double-id）。
  - （可选）提供一个 FAM 文件，统计交并集数量，并仅保留在 FAM 中的样本生成新的 PLINK 输出。

进程提示：
  - 函数运行时会打印关键步骤（开始/完成）和统计信息，便于在终端追踪进度。

输出：
  - 返回 dict，其中包含输入/输出路径、样本数量统计以及最终汇总字符串（summary_text）。
"""

from __future__ import annotations

import os
import subprocess
import re
import pandas as pd
from pathlib import Path
from typing import Optional, Dict, Any, Tuple, Set


def _echo(msg: str) -> None:
    """统一的进程提示输出。"""
    print(f"[array_check_prepare_tools] {msg}")


# === 新增：根据 IID 前缀统计病例/对照数量 ===
def _count_case_ctrl(iids: Set[str], case_prefix: str) -> Dict[str, int]:
    """根据 IID 前缀统计病例/对照数量（病例 IID 以 case_prefix 开头）。"""
    n_case = sum(1 for x in iids if str(x).startswith(case_prefix))
    n_total = len(iids)
    n_ctrl = n_total - n_case
    return {"n_case": n_case, "n_ctrl": n_ctrl, "n_total": n_total}


# === 子进程执行工具 ===
def _run_cmd(cmd: list[str]) -> None:
    """以检查错误的方式运行命令。"""
    proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    if proc.returncode != 0:
        raise RuntimeError(f"Command failed (exit={proc.returncode}):\n{' '.join(cmd)}\n--- OUTPUT ---\n{proc.stdout}")
    _echo("命令执行成功。")


# === FAM 读取与统计工具 ===
def _read_fam_iids(fam_path: Path) -> Tuple[Set[str], int]:
    """读取 .fam 文件，返回 IID 集合与总行数（样本数）。"""
    iids: Set[str] = set()
    total = 0
    with fam_path.open("r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            fid, iid = parts[0], parts[1]
            # 我们只基于 IID 统计重叠（因为转换产物使用 --double-id，FID==IID）
            iids.add(iid)
            total += 1
    return iids, total


# === 汇总统计文件处理工具 ===
def _read_sum_stat_valid_variants(sum_stat_path: Path, id_col: str = "ID", p_col: str = "P") -> Tuple[Set[str], int, int]:
    """
    读取 GWAS 汇总统计文件，返回 P 值有效的变体 ID 集合。
    
    参数
    ------
    sum_stat_path : Path
        汇总统计文件路径
    id_col : str, default "ID"
        变体 ID 列名
    p_col : str, default "P"
        P 值列名
        
    返回
    ------
    Tuple[Set[str], int, int]
        (有效变体ID集合, 总变体数, 无效P值变体数)
    """
    import pandas as pd
    import numpy as np
    
    # 读取汇总统计文件
    try:
        # 尝试不同的分隔符
        if sum_stat_path.suffix.lower() == '.csv':
            df = pd.read_csv(sum_stat_path)
        else:
            # 尝试制表符分隔
            df = pd.read_csv(sum_stat_path, sep='\t')
            # 如果只有一列，尝试空格分隔
            if df.shape[1] == 1:
                df = pd.read_csv(sum_stat_path, sep=r'\s+')
    except Exception as e:
        raise ValueError(f"无法读取汇总统计文件 {sum_stat_path}: {e}")
    
    # 检查必需的列是否存在
    if id_col not in df.columns:
        raise ValueError(f"汇总统计文件中未找到 ID 列: {id_col}。可用列: {list(df.columns)}")
    if p_col not in df.columns:
        raise ValueError(f"汇总统计文件中未找到 P 值列: {p_col}。可用列: {list(df.columns)}")
    
    total_variants = len(df)
    
    # 识别有效的 P 值
    # 转换 P 值列为数值，无效值会变成 NaN
    p_values = pd.to_numeric(df[p_col], errors='coerce')
    
    # 找出有效的 P 值（非 NaN、非负、<=1）
    valid_p_mask = (
        ~p_values.isna() &  # 不是 NaN
        (p_values >= 0) &   # 非负
        (p_values <= 1)     # 小于等于 1
    )
    
    valid_variants = set(df.loc[valid_p_mask, id_col].astype(str))
    invalid_p_count = total_variants - len(valid_variants)
    
    return valid_variants, total_variants, invalid_p_count


# === BIM 读取与统计工具 ===
def _read_bim_variants(bim_path: Path) -> Tuple[Set[str], int]:
    """读取 .bim 文件，返回变体 ID 集合与总行数（变体数）。"""
    variant_ids: Set[str] = set()
    total = 0
    with bim_path.open("r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            # BIM 文件格式：CHR SNP_ID GENETIC_DIST BP_COORD ALLELE1 ALLELE2
            # 我们使用第2列的 SNP_ID 作为变体标识符
            variant_id = parts[1]
            variant_ids.add(variant_id)
            total += 1
    return variant_ids, total


def _parse_and_reorder_covar_names(covar_file: str, covar_name: str) -> str:
    """
    读取协变量文件头部，按文件中的顺序重新排列协变量名称。
    支持 PC 范围表示法（如 PC1_AVG-PC10_AVG）。
    
    参数
    ------
    covar_file : str
        协变量文件路径
    covar_name : str
        用户提供的协变量名称字符串，用逗号分隔
        
    返回
    ------
    str
        按文件头部顺序重新排列的协变量名称字符串
    """
    import re
    
    # 解析用户提供的协变量名称
    requested_covars = []
    for name in covar_name.split(','):
        name = name.strip()
        # 检查是否是 PC 范围表示法（如 PC1_AVG-PC10_AVG）
        pc_range_match = re.match(r'PC(\d+)_AVG-PC(\d+)_AVG', name)
        if pc_range_match:
            start_pc = int(pc_range_match.group(1))
            end_pc = int(pc_range_match.group(2))
            for i in range(start_pc, end_pc + 1):
                requested_covars.append(f'PC{i}_AVG')
        else:
            requested_covars.append(name)
    
    # 读取文件头部
    with open(covar_file, 'r') as f:
        header_line = f.readline().strip()
        # 处理可能的制表符分隔
        if '\t' in header_line:
            file_columns = [col.strip() for col in header_line.split('\t')]
        else:
            file_columns = [col.strip() for col in header_line.split(',')]
    
    # 移除可能的 # 前缀
    file_columns = [col.lstrip('#') for col in file_columns]
    
    # 按文件中的顺序重新排列协变量
    reordered_covars = []
    for col in file_columns:
        if col in requested_covars:
            reordered_covars.append(col)
    
    # 检查是否有请求的协变量不在文件中
    missing_covars = set(requested_covars) - set(file_columns)
    if missing_covars:
        _echo(f"警告：以下协变量在文件中未找到：{missing_covars}")
    
    reordered_str = ', '.join(reordered_covars)
    _echo(f"原始协变量表示：{covar_name}")
    _echo(f"展开后的协变量：{', '.join(requested_covars)}")
    _echo(f"按文件顺序重排：{reordered_str}")
    
    return reordered_str


def _format_summary(result: Dict[str, Any], case_prefix: str = "PHOM") -> str:
    lines = [
        "=== 运行统计汇总 ===",
        f"输入 VCF：{result.get('vcf_gz','')}",
        f"初次转换输出前缀：{result.get('plink_out','')}",
        f"初次转换样本数：{result.get('n_samples_plink','NA')}",
        f"初次转换变体数：{result.get('n_variants_plink','NA')}",
        f"初次转换病例数({case_prefix}*)：{result.get('n_case_plink','NA')}",
        f"初次转换对照数：{result.get('n_ctrl_plink','NA')}",
    ]
    
    # 样本筛选信息
    if 'keep_fam_path' in result:
        lines.extend([
            f"--- 样本筛选 ---",
            f"keep FAM：{result.get('keep_fam_path','')}",
            f"与 FAM 重叠样本数：{result.get('n_in_keep','NA')}",
            f"不在 FAM 的样本数：{result.get('n_not_in_keep','NA')}",
            f"保留集病例数({case_prefix}*)：{result.get('n_case_keep','NA')}",
            f"保留集对照数：{result.get('n_ctrl_keep','NA')}",
            f"未保留病例数({case_prefix}*)：{result.get('n_case_not_in_keep','NA')}",
            f"未保留对照数：{result.get('n_ctrl_not_in_keep','NA')}",
            f"样本筛选后输出前缀：{result.get('filtered_out','(未启用)')}",
        ])
    
    # 变体筛选信息
    if 'extract_bim_path' in result:
        lines.extend([
            f"--- 变体筛选 ---",
            f"extract BIM：{result.get('extract_bim_path','')}",
            f"与 BIM 重叠变体数：{result.get('n_variants_extract','NA')}",
            f"不在 BIM 的变体数：{result.get('n_variants_not_extract','NA')}",
            f"变体筛选后输出前缀：{result.get('extracted_out','(未启用)')}",
        ])
    
    # 汇总统计筛选信息
    if 'sum_stat_path' in result:
        lines.extend([
            f"--- 汇总统计 P 值筛选 ---",
            f"汇总统计文件：{result.get('sum_stat_path','')}",
            f"汇总统计总变体数：{result.get('n_variants_sumstat_total','NA')}",
            f"汇总统计文件中P值异常变体数：{result.get('n_variants_sumstat_invalid_file','NA')}",
            f"当前数据中与P值异常变体交集：{result.get('n_variants_sumstat_invalid','NA')} 个（已过滤）",
            f"过滤后最终保留变体数：{result.get('n_variants_sumstat_valid','NA')} 个",
            f"汇总统计筛选后输出前缀：{result.get('sumstat_out','(未启用)')}",
        ])
    
    # 最终输出
    lines.extend([
        f"--- 最终结果 ---",
        f"最终输出前缀：{result.get('final_out',result.get('plink_out',''))}",
    ])
    
    return "\n".join(lines)


def vcf_to_plink_with_optional_keep(
    vcf_gz: str,
    out_prefix: str,
    plink2_path: str = "/home/b/b37974/plink2",
    enable_keep: bool = False,
    keep_fam_path: Optional[str] = None,
    filtered_out_prefix: Optional[str] = None,
    enable_extract: bool = False,
    extract_bim_path: Optional[str] = None,
    extract_out_prefix: Optional[str] = None,
    sum_stat: Optional[str] = None,
    sum_stat_id_col: str = "ID",
    sum_stat_p_col: str = "P",
    sum_stat_out_prefix: Optional[str] = None,
    threads: Optional[int] = None,
    out_dir: Optional[str] = None,
    case_prefix: str = "PHOM",
) -> Dict[str, Any]:
    """
    使用 plink2 将 VCF（.vcf.gz）转换为 PLINK 二进制格式（.bed/.bim/.fam），并可选基于给定 FAM 进行样本筛选、基于 BIM 进行变体筛选，以及基于汇总统计文件过滤 P 值异常的变体。

    参数
    ------
    vcf_gz : str
        输入的压缩 VCF 路径（.vcf.gz）。
    out_prefix : str
        首次转换后的输出前缀（将生成 .bed/.bim/.fam）。
    plink2_path : str, default "/home/b/b37974/plink2"
        plink2 可执行文件路径。
    enable_keep : bool, default False
        是否启用基于 FAM 的样本筛选。
    keep_fam_path : Optional[str], default None
        当 enable_keep=True 时，需要提供的 FAM 路径（两列 FID、IID 也可）。
    filtered_out_prefix : Optional[str], default None
        启用 keep 时筛选后输出的前缀；若为 None，则使用 f"{out_prefix}.keep"。
    enable_extract : bool, default False
        是否启用基于 BIM 的变体筛选。
    extract_bim_path : Optional[str], default None
        当 enable_extract=True 时，需要提供的 BIM 路径，用于提取共有变体。
    extract_out_prefix : Optional[str], default None
        启用 extract 时筛选后输出的前缀；若为 None，则使用 f"{out_prefix}.extract"。
    sum_stat : Optional[str], default None
        GWAS 汇总统计文件路径，用于过滤 P 值异常的变体。
    sum_stat_id_col : str, default "ID"
        汇总统计文件中变体 ID 列名。
    sum_stat_p_col : str, default "P"
        汇总统计文件中 P 值列名。
    sum_stat_out_prefix : Optional[str], default None
        基于汇总统计筛选后输出的前缀；若为 None，则使用 f"{current_prefix}.sumstat"。
    threads : Optional[int], default None
        若提供，则在 plink2 调用中添加 "--threads {threads}"。
    out_dir : Optional[str], default None
        输出目录，若为 None 则使用 tmp/ 目录。
    case_prefix : str, default "PHOM"
        用于区分病例的 IID 前缀，默认为 "PHOM"。

    返回
    ------
    Dict[str, Any]
        关键字段包含：
        - "vcf_gz": 输入 VCF 路径
        - "plink_out": 初次转换的 bed/bim/fam 前缀
        - "filtered_out": 二次筛选输出前缀（若启用 keep）
        - "extracted_out": 变体筛选输出前缀（若启用 extract）
        - "sumstat_out": 汇总统计筛选输出前缀（若提供 sum_stat）
        - "final_out": 最终输出前缀（经过所有筛选步骤后）
        - "n_samples_plink": 初次转换得到的样本数
        - "n_variants_plink": 初次转换得到的变体数
        - "n_case_plink": 初次转换病例数
        - "n_ctrl_plink": 初次转换对照数
        - "n_in_keep": 与 keep FAM 交集的样本数（若启用 keep）
        - "n_not_in_keep": 不在 keep FAM 中的样本数（若启用 keep）
        - "n_variants_extract": 与 extract BIM 交集的变体数（若启用 extract）
        - "n_variants_not_extract": 不在 extract BIM 中的变体数（若启用 extract）
        - "n_variants_sumstat_valid": P 值有效的变体数（若提供 sum_stat）
        - "n_variants_sumstat_invalid": P 值异常的变体数（若提供 sum_stat）
        - "case_prefix": 用于病例识别的前缀
    """

    _echo("开始任务：VCF 转换为 PLINK 二进制格式（--double-id）")

    vcf_path = Path(vcf_gz)
    if not vcf_path.exists():
        raise FileNotFoundError(f"VCF not found: {vcf_gz}")

    # === 默认输出目录逻辑：若未指定 out_dir，则输出到当前目录下 tmp/，并自动创建该目录 ===
    # 若 out_dir 未指定，则默认输出到 ./tmp/ 目录
    # 中文注释：若用户未指定输出目录，默认所有输出文件存放在当前工作目录下的 tmp/ 子目录
    if out_dir is None:
        default_tmp_dir = Path.cwd() / "tmp"
        if not default_tmp_dir.exists():
            default_tmp_dir.mkdir(parents=True, exist_ok=True)
        out_prefix_path = default_tmp_dir / out_prefix
        _echo(f"未指定 out_dir，默认输出到 tmp/ 目录：{default_tmp_dir}")
    else:
        out_dir_path = Path(out_dir)
        if not out_dir_path.exists():
            out_dir_path.mkdir(parents=True, exist_ok=True)
        out_prefix_path = out_dir_path / out_prefix
        _echo(f"指定输出目录：{out_dir_path}")

    # 保留原有输出路径逻辑：out_prefix_path 已拼接在上面
    out_dir_final = out_prefix_path.parent
    if out_dir_final and not out_dir_final.exists():
        out_dir_final.mkdir(parents=True, exist_ok=True)
    _echo(f"输出目录已准备：{out_dir_final if str(out_dir_final) else '.'}")

    # 组装 plink2 导入命令
    cmd = [
        plink2_path,
        "--vcf", str(vcf_path),
        "--double-id",  # 使 FID 与 IID 都等于样本 ID
        "--make-bed",
        "--out", str(out_prefix_path),
    ]
    # 可选线程数（HPC 环境建议显式设置）
    if threads is not None:
        cmd.extend(["--threads", str(threads)])

    _echo("步骤 1/2：执行 plink2 VCF 导入与 --make-bed（确保 FID=IID）")
    _echo("命令：" + ' '.join(cmd))
    _run_cmd(cmd)

    fam_path = out_prefix_path.with_suffix(".fam")
    if not fam_path.exists():
        # 有些 plink2 会生成 .fam 而非 .psam；若没有 .fam 则报错
        raise FileNotFoundError(f"Expected FAM not found: {fam_path}")
    _echo(f".fam 文件就绪：{fam_path}")

    # 统计初次转换样本数和变体数
    plink_iids, n_samples_plink = _read_fam_iids(fam_path)
    bim_path = out_prefix_path.with_suffix(".bim")
    plink_variants, n_variants_plink = _read_bim_variants(bim_path)
    _echo(f"转换完成。样本数 = {n_samples_plink}，变体数 = {n_variants_plink}")

    # 新增：统计初次转换的病例/对照数量
    cc_plink = _count_case_ctrl(plink_iids, case_prefix)
    _echo(f"初次转换：病例({case_prefix}*) = {cc_plink['n_case']}，对照 = {cc_plink['n_ctrl']}，合计 = {cc_plink['n_total']}")

    result: Dict[str, Any] = {
        "vcf_gz": str(vcf_path),
        "plink_out": str(out_prefix_path),
        "n_samples_plink": n_samples_plink,
        "n_variants_plink": n_variants_plink,
        "n_case_plink": cc_plink["n_case"],
        "n_ctrl_plink": cc_plink["n_ctrl"],
        "case_prefix": case_prefix,
    }

    # 记录当前处理阶段的输出前缀，用于后续步骤的输入
    current_prefix = out_prefix_path

    # 2) 可选：基于 FAM 筛选样本
    if enable_keep:
        _echo("步骤 2/N：启用 --keep 基于 FAM 的样本筛选")
        if not keep_fam_path:
            raise ValueError("enable_keep=True 时必须提供 keep_fam_path")
        keep_fam = Path(keep_fam_path)
        if not keep_fam.exists():
            raise FileNotFoundError(f"Keep FAM not found: {keep_fam_path}")

        keep_iids, _ = _read_fam_iids(keep_fam)
        _echo(f"读取 keep FAM：{keep_fam}")

        # 计算与 keep FAM 的交集/差集样本数量
        n_in = len(plink_iids & keep_iids)
        n_not_in = len(plink_iids - keep_iids)
        _echo(f"与转换结果重叠样本数 = {n_in}；不在 FAM 中 = {n_not_in}")

        # 新增：统计保留集和未保留集的病例/对照数量
        kept_iids = plink_iids & keep_iids
        notin_iids = plink_iids - keep_iids
        cc_keep = _count_case_ctrl(kept_iids, case_prefix)
        cc_notin = _count_case_ctrl(notin_iids, case_prefix)
        _echo(f"筛选后（保留集）：病例({case_prefix}*) = {cc_keep['n_case']}，对照 = {cc_keep['n_ctrl']}，合计 = {cc_keep['n_total']}")
        _echo(f"未保留（差集）：病例({case_prefix}*) = {cc_notin['n_case']}，对照 = {cc_notin['n_ctrl']}，合计 = {cc_notin['n_total']}")

        result.update({
            "keep_fam_path": str(keep_fam),
            "n_in_keep": n_in,
            "n_not_in_keep": n_not_in,
            "n_case_keep": cc_keep["n_case"],
            "n_ctrl_keep": cc_keep["n_ctrl"],
            "n_case_not_in_keep": cc_notin["n_case"],
            "n_ctrl_not_in_keep": cc_notin["n_ctrl"],
        })

        out2_prefix = Path(filtered_out_prefix) if filtered_out_prefix else Path(str(out_prefix_path) + ".keep")
        _echo(f"输出样本筛选后前缀：{out2_prefix}")

        cmd2 = [
            plink2_path,
            "--bfile", str(current_prefix),
            "--keep", str(keep_fam),  # 直接使用用户给的 FAM/两列文件
            "--make-bed",
            "--out", str(out2_prefix),
        ]
        if threads is not None:
            cmd2.extend(["--threads", str(threads)])

        _echo("命令：" + ' '.join(cmd2))
        _run_cmd(cmd2)
        _echo("样本筛选完成并已输出新的 PLINK 文件。")
        result["filtered_out"] = str(out2_prefix)
        current_prefix = out2_prefix  # 更新当前前缀

    # 3) 可选：基于 BIM 筛选变体
    if enable_extract:
        _echo("步骤 3/N：启用 --extract 基于 BIM 的变体筛选")
        if not extract_bim_path:
            raise ValueError("enable_extract=True 时必须提供 extract_bim_path")
        extract_bim = Path(extract_bim_path)
        if not extract_bim.exists():
            raise FileNotFoundError(f"Extract BIM not found: {extract_bim_path}")

        extract_variants, _ = _read_bim_variants(extract_bim)
        _echo(f"读取 extract BIM：{extract_bim}")

        # 获取当前数据的变体集合
        current_bim_path = current_prefix.with_suffix(".bim")
        current_variants, _ = _read_bim_variants(current_bim_path)

        # 计算与 extract BIM 的交集/差集变体数量
        n_variants_in = len(current_variants & extract_variants)
        n_variants_not_in = len(current_variants - extract_variants)
        _echo(f"与当前数据重叠变体数 = {n_variants_in}；不在 extract BIM 中 = {n_variants_not_in}")

        result.update({
            "extract_bim_path": str(extract_bim),
            "n_variants_extract": n_variants_in,
            "n_variants_not_extract": n_variants_not_in,
        })

        # 创建变体 ID 列表文件用于 --extract
        extract_list_file = current_prefix.parent / f"{current_prefix.name}.extract_list.txt"
        with extract_list_file.open("w") as f:
            for variant_id in sorted(current_variants & extract_variants):
                f.write(f"{variant_id}\n")
        _echo(f"创建变体列表文件：{extract_list_file}")

        out3_prefix = Path(extract_out_prefix) if extract_out_prefix else Path(str(current_prefix) + ".extract")
        _echo(f"输出变体筛选后前缀：{out3_prefix}")

        cmd3 = [
            plink2_path,
            "--bfile", str(current_prefix),
            "--extract", str(extract_list_file),
            "--make-bed",
            "--out", str(out3_prefix),
        ]
        if threads is not None:
            cmd3.extend(["--threads", str(threads)])

        _echo("命令：" + ' '.join(cmd3))
        _run_cmd(cmd3)
        _echo("变体筛选完成并已输出新的 PLINK 文件。")
        result["extracted_out"] = str(out3_prefix)
        current_prefix = out3_prefix  # 更新当前前缀

    # 4) 可选：基于汇总统计文件筛选 P 值有效的变体
    if sum_stat:
        _echo("步骤 4/N：启用基于汇总统计文件的 P 值筛选")
        sum_stat_path = Path(sum_stat)
        if not sum_stat_path.exists():
            raise FileNotFoundError(f"汇总统计文件未找到: {sum_stat}")

        try:
            valid_variants, total_sumstat, invalid_p_count = _read_sum_stat_valid_variants(
                sum_stat_path, sum_stat_id_col, sum_stat_p_col
            )
            _echo(f"读取汇总统计文件：{sum_stat_path}")
            _echo(f"汇总统计文件总变体数：{total_sumstat}")
            _echo(f"汇总统计文件中 P 值异常变体数：{invalid_p_count}")
        except Exception as e:
            raise ValueError(f"处理汇总统计文件时出错: {e}")

        # 获取当前数据的变体集合
        current_bim_path = current_prefix.with_suffix(".bim")
        current_variants, _ = _read_bim_variants(current_bim_path)

        # 计算汇总统计文件中的无效变体集合
        import pandas as pd
        import numpy as np
        
        # 重新读取汇总统计文件获取无效变体
        if sum_stat_path.suffix.lower() == '.csv':
            df = pd.read_csv(sum_stat_path)
        else:
            df = pd.read_csv(sum_stat_path, sep='\t')
            if df.shape[1] == 1:
                df = pd.read_csv(sum_stat_path, sep=r'\s+')
        
        p_values = pd.to_numeric(df[sum_stat_p_col], errors='coerce')
        invalid_p_mask = (
            p_values.isna() |       # 是 NaN
            (p_values < 0) |        # 负值
            (p_values > 1)          # 大于 1
        )
        invalid_variants = set(df.loc[invalid_p_mask, sum_stat_id_col].astype(str))
        
        # 计算当前数据中与汇总统计无效变体的交集（需要过滤的变体）
        invalid_overlap = current_variants & invalid_variants
        n_invalid_overlap = len(invalid_overlap)
        
        # 计算与汇总统计有效变体的交集（保留的变体）
        valid_overlap = current_variants & valid_variants
        n_valid_overlap = len(valid_overlap)
        n_current_total = len(current_variants)
        
        _echo(f"当前数据变体数：{n_current_total}")
        _echo(f"当前数据中与汇总统计P值异常变体交集：{n_invalid_overlap} 个（需要过滤）")
        _echo(f"当前数据中与汇总统计P值有效变体交集：{n_valid_overlap} 个（保留）")
        _echo(f"因此过滤掉 {n_invalid_overlap} 个变体，最终保留 {n_valid_overlap} 个变体")

        result.update({
            "sum_stat_path": str(sum_stat_path),
            "n_variants_sumstat_total": total_sumstat,
            "n_variants_sumstat_invalid_file": invalid_p_count,  # 汇总统计文件中P值异常的变体数
            "n_variants_sumstat_valid": n_valid_overlap,
            "n_variants_sumstat_invalid": n_invalid_overlap,  # 当前数据中与P值异常变体的交集数
        })

        # 只有当存在需要过滤的变体时才进行筛选
        if n_invalid_overlap > 0:
            # 创建有效变体 ID 列表文件用于 --extract
            sumstat_list_file = current_prefix.parent / f"{current_prefix.name}.sumstat_valid_list.txt"
            with sumstat_list_file.open("w") as f:
                for variant_id in sorted(valid_overlap):
                    f.write(f"{variant_id}\n")
            _echo(f"创建 P 值有效变体列表文件：{sumstat_list_file}")

            out4_prefix = Path(sum_stat_out_prefix) if sum_stat_out_prefix else Path(str(current_prefix) + ".sumstat")
            _echo(f"输出汇总统计筛选后前缀：{out4_prefix}")

            cmd4 = [
                plink2_path,
                "--bfile", str(current_prefix),
                "--extract", str(sumstat_list_file),
                "--make-bed",
                "--out", str(out4_prefix),
            ]
            if threads is not None:
                cmd4.extend(["--threads", str(threads)])

            _echo("命令：" + ' '.join(cmd4))
            _run_cmd(cmd4)
            _echo("汇总统计 P 值筛选完成并已输出新的 PLINK 文件。")
            result["sumstat_out"] = str(out4_prefix)
            current_prefix = out4_prefix  # 更新当前前缀
        else:
            _echo("所有变体的 P 值都有效，无需进一步筛选。")
            result["sumstat_out"] = str(current_prefix)

    # 记录最终输出前缀
    result["final_out"] = str(current_prefix)

    # === 最终结果汇总输出 ===
    summary_text = _format_summary(result, case_prefix=case_prefix)
    _echo("\n" + summary_text)
    result["summary_text"] = summary_text

    return result


def run_plink2_gwas_association(
    bed_prefix: str,
    pheno_file: str,
    pheno_name: str = "PHENO1",
    covar_file: Optional[str] = None,
    covar_name: Optional[str] = None,
    model: str = "additive",
    out_prefix: str = "gwas_result",
    plink2_path: str = "/home/b/b37974/plink2",
    threads: int = 16,
    ci_level: float = 0.95,
    out_dir: Optional[str] = None,
    extra_opts: Optional[list] = None,
) -> Dict[str, Any]:
    """
    使用 plink2 进行 GWAS 关联分析。

    参数
    ------
    bed_prefix : str
        输入的 PLINK 二进制文件前缀（.bed/.bim/.fam）。
    pheno_file : str
        表型文件路径（CSV 格式）。
    pheno_name : str, default "PHENO1"
        表型列名。
    covar_file : Optional[str], default None
        协变量文件路径（CSV 格式）。
    covar_name : Optional[str], default None
        协变量列名，多个用逗号分隔，如 "SEX, PC1_AVG-PC10_AVG"。
    model : str, default "additive"
        遗传模型，可选：'additive', 'dominant', 'recessive', 'genotypic'。
    out_prefix : str, default "gwas_result"
        输出文件前缀。
    plink2_path : str, default "/home/b/b37974/plink2"
        plink2 可执行文件路径。
    threads : int, default 16
        线程数。
    ci_level : float, default 0.95
        置信区间水平。
    out_dir : Optional[str], default None
        输出目录，若为 None 则使用 tmp/ 目录。
    extra_opts : Optional[list], default None
        额外的 plink2 选项列表，如 ["--maf", "0.01", "--geno", "0.1"]。

    返回
    ------
    Dict[str, Any]
        包含输入参数、输出路径、命令等信息的结果字典。
    """

    _echo("开始任务：使用 plink2 进行 GWAS 关联分析")

    # 检查输入文件是否存在
    bed_path = Path(bed_prefix + ".bed")
    bim_path = Path(bed_prefix + ".bim")
    fam_path = Path(bed_prefix + ".fam")
    pheno_path = Path(pheno_file)

    if not bed_path.exists():
        raise FileNotFoundError(f"BED file not found: {bed_path}")
    if not bim_path.exists():
        raise FileNotFoundError(f"BIM file not found: {bim_path}")
    if not fam_path.exists():
        raise FileNotFoundError(f"FAM file not found: {fam_path}")
    if not pheno_path.exists():
        raise FileNotFoundError(f"Phenotype file not found: {pheno_path}")

    if covar_file:
        covar_path = Path(covar_file)
        if not covar_path.exists():
            raise FileNotFoundError(f"Covariate file not found: {covar_path}")

    # 设置输出目录
    if out_dir is None:
        default_tmp_dir = Path.cwd() / "tmp"
        if not default_tmp_dir.exists():
            default_tmp_dir.mkdir(parents=True, exist_ok=True)
        out_full_prefix = default_tmp_dir / out_prefix
        _echo(f"未指定 out_dir，默认输出到 tmp/ 目录：{default_tmp_dir}")
    else:
        out_dir_path = Path(out_dir)
        if not out_dir_path.exists():
            out_dir_path.mkdir(parents=True, exist_ok=True)
        out_full_prefix = out_dir_path / out_prefix
        _echo(f"指定输出目录：{out_dir_path}")

    # 构建 plink2 命令
    cmd = [
        plink2_path,
        "--bfile", bed_prefix,
        "--pheno", str(pheno_path),
        "--pheno-name", pheno_name,
    ]

    # 添加协变量（如果提供）
    if covar_file and covar_name:
        # 清理协变量名称：去除逗号后的空格，plink2 对空格敏感
        cleaned_covar_name = ','.join([name.strip() for name in covar_name.split(',')])
        cmd.extend(["--covar", str(covar_path)])
        cmd.extend(["--covar-name", cleaned_covar_name])
        _echo(f"原始协变量名称：{covar_name}")
        _echo(f"清理后协变量名称：{cleaned_covar_name}")
        # 更新 covar_name 为清理后的值，用于后续记录
        covar_name = cleaned_covar_name

    # 设置遗传模型
    if model == "additive":
        model_opt = "--glm"
    else:
        model_opt = f"--glm {model}"
    
    cmd.extend(model_opt.split())
    
    # 添加标准选项
    cmd.extend(["omit-ref", "no-firth", "hide-covar"])
    
    # 输出和置信区间
    cmd.extend(["--out", str(out_full_prefix)])
    cmd.extend(["--ci", str(ci_level)])
    
    # 线程数
    cmd.extend(["--threads", str(threads)])

    # 添加额外选项
    if extra_opts:
        cmd.extend(extra_opts)

    _echo(f"执行 GWAS 关联分析")
    _echo("命令：" + ' '.join(cmd))
    
    # 执行命令
    _run_cmd(cmd)

    # 检查输出文件
    expected_output = out_full_prefix.with_suffix(".PHENO1.glm.logistic")
    if not expected_output.exists():
        # 尝试其他可能的输出文件扩展名
        possible_exts = [".PHENO1.glm.linear", f".{pheno_name}.glm.logistic", f".{pheno_name}.glm.linear"]
        found_output = None
        for ext in possible_exts:
            test_file = Path(str(out_full_prefix) + ext)
            if test_file.exists():
                found_output = test_file
                break
        
        if found_output:
            expected_output = found_output
        else:
            _echo(f"警告：未找到预期的输出文件，请检查输出目录：{out_full_prefix.parent}")

    result = {
        "bed_prefix": bed_prefix,
        "pheno_file": str(pheno_path),
        "pheno_name": pheno_name,
        "covar_file": str(covar_path) if covar_file else None,
        "covar_name": covar_name,
        "model": model,
        "out_prefix": str(out_full_prefix),
        "output_file": str(expected_output) if expected_output.exists() else None,
        "command": ' '.join(cmd),
        "threads": threads,
        "ci_level": ci_level,
    }

    # 汇总输出
    summary_lines = [
        "=== GWAS 关联分析完成 ===",
        f"输入 BED 前缀：{bed_prefix}",
        f"表型文件：{pheno_path}",
        f"表型名称：{pheno_name}",
        f"协变量文件：{covar_path if covar_file else '未使用'}",
        f"协变量名称：{covar_name if covar_name else '未使用'}",
        f"遗传模型：{model}",
        f"输出前缀：{out_full_prefix}",
        f"结果文件：{expected_output if expected_output.exists() else '请检查输出目录'}",
        f"线程数：{threads}",
        f"置信区间：{ci_level}",
    ]
    
    summary_text = "\n".join(summary_lines)
    _echo("\n" + summary_text)
    result["summary_text"] = summary_text

    return result


def harmonize_and_convert_sumstat(
    plink_stats_file: str,
    ref_seq_path: str = "/LARGE0/gr10478/b37974/Pulmonary_Hypertension/nagasaki_pipeline/data/hs38DH.fa",
    output_prefix: Optional[str] = None,
    log_path: Optional[str] = None,
    verbose: Optional[bool] = False,
    enable_harmonize: bool = False,
) -> Dict[str, Any]:
    """
    功能：读取 PLINK2 glm.logistic 汇总统计文件，对整个数据集进行 OR 到 beta 转换和可选的等位基因方向协调。
    
    参数说明：
        plink_stats_file (str): PLINK2 结果文件路径（*.glm.logistic，gwaslab 支持 fmt="plink2"）。
        ref_seq_path (str): 参考基因组 fasta 路径（用于 harmonize），默认与 fine_map_susie_tools.py 中一致。
        output_prefix (Optional[str]): 输出文件前缀。若为 None，则自动用输入文件名生成。
        log_path (Optional[str]): 日志文件路径；若为 None，则自动使用 `{output_prefix}.log`。
        verbose (Optional[bool]): 是否在调用 gwaslab 的相关函数时开启详细输出；默认 False。
        enable_harmonize (bool): 是否启用等位基因方向协调，默认 False（不启用）。
    
    输出制品：
        - {output_prefix}.harmonized_sumstat.tsv    ：协调后的汇总统计文件（Tab 分隔）
        - {output_prefix}.summary.json              ：结构化记录输入参数与输出文件位置
        - {output_prefix}.log                       ：同步中文日志
    
    返回：
        Dict[str, Any]：包含输入/输出路径、处理统计信息等的结果字典。
    
    依赖：
        - gwaslab（函数内延迟导入）
        - pandas / numpy（已在脚本层导入）
    
    备注：
        - 本函数对 gwaslab 的 `verbose` 参数采用"尽力而为"的兼容策略。
        - 使用与 fine_map_susie_tools.py 相同的逻辑进行 OR->beta 转换。
        - 等位基因协调功能默认关闭，可通过 enable_harmonize=True 启用。
    """
    import os
    import json
    import time
    import logging
    from datetime import datetime
    from pathlib import Path
    
    def _call_with_verbose(fn, *args, **kwargs):
        """兼容 gwaslab verbose 参数的调用封装"""
        if verbose is not None:
            kwargs_with_verbose = dict(kwargs)
            kwargs_with_verbose["verbose"] = verbose
            try:
                return fn(*args, **kwargs_with_verbose)
            except TypeError:
                # 兼容旧版本 gwaslab：该方法不接受 verbose 参数时退回不传
                return fn(*args, **kwargs)
        return fn(*args, **kwargs)

    # 延迟导入：仅在函数内部导入 gwaslab
    try:
        import gwaslab as gl # type: ignore
    except ImportError:
        raise ImportError("需要安装 gwaslab 库：pip install gwaslab")

    # ---------- 准备路径 ----------
    if output_prefix is None:
        base = os.path.splitext(os.path.basename(plink_stats_file))[0]
        output_prefix = os.path.join(os.path.dirname(plink_stats_file), f"{base}.harmonized")

    harmonized_tsv = f"{output_prefix}.harmonized_sumstat.tsv"
    summary_json = f"{output_prefix}.summary.json"
    if log_path is None:
        log_path = f"{output_prefix}.log"

    output_dir = os.path.dirname(output_prefix)
    if output_dir == "":
        output_dir = os.getcwd()
        output_prefix = os.path.join(output_dir, output_prefix)
    os.makedirs(output_dir, exist_ok=True)

    # ---------- 设置日志（中文） ----------
    logger = logging.getLogger("harmonize_sumstat")
    logger.setLevel(logging.INFO)
    # 清除重复 handler（避免在交互式环境重复添加）
    if logger.handlers:
        for h in list(logger.handlers):
            logger.removeHandler(h)
    fh = logging.FileHandler(log_path, mode="w", encoding="utf-8")
    fmt = logging.Formatter("[%(asctime)s] %(levelname)s - %(message)s")
    fh.setFormatter(fmt)
    logger.addHandler(fh)

    logger.info("启动流程：汇总统计文件协调与 beta 转换")
    logger.info(f"输入文件：{plink_stats_file}")
    logger.info(f"参考序列（用于harmonize）：{ref_seq_path}")
    logger.info(f"是否启用等位基因协调：{enable_harmonize}")

    t0 = time.time()

    # ---------- 读取与基础检查 ----------
    try:
        sumstats = gl.Sumstats(
            plink_stats_file,
            fmt="plink2",
            build="38",
            ea="A1",
            nea="OMITTED",
            OR_95L="L95",
            OR_95U="U95",
            verbose=verbose, # type: ignore
        )
    except TypeError:
        # 兼容不支持 verbose 参数的 gwaslab 版本
        sumstats = gl.Sumstats(
            plink_stats_file,
            fmt="plink2",
            build="38",
            ea="A1",
            nea="OMITTED",
            OR_95L="L95",
            OR_95U="U95",
        )
    
    logger.info("已创建 gwaslab.Sumstats 对象，开始 basic_check()")
    _call_with_verbose(sumstats.basic_check)
    logger.info("basic_check() 完成")

    # 记录原始数据统计
    n_variants_original = len(sumstats.data)
    logger.info(f"原始变体数：{n_variants_original}")

    # ---------- OR 到 beta 转换 ----------
    logger.info("开始 OR 到 beta 转换")
    _call_with_verbose(sumstats.fill_data, to_fill=["BETA"]) # type: ignore
    logger.info("OR 到 beta 转换完成")

    # ---------- 等位基因方向协调（可选） ----------
    if enable_harmonize:
        logger.info("开始等位基因方向协调")
        _call_with_verbose(sumstats.harmonize, basic_check=False, ref_seq=ref_seq_path) # type: ignore
        logger.info("等位基因方向协调完成")
    else:
        logger.info("跳过等位基因方向协调（enable_harmonize=False）")

    # 记录处理后数据统计
    n_variants_final = len(sumstats.data)
    logger.info(f"最终变体数：{n_variants_final}")

    # ---------- 保存协调后的数据 ----------
    sumstats.data.to_csv(harmonized_tsv, sep="\t", index=False) # type: ignore
    logger.info(f"协调后汇总统计已写出：{harmonized_tsv}")

    # ---------- 写出与返回 JSON 摘要 ----------
    meta = {
        "created_at": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        "inputs": {
            "plink_stats_file": plink_stats_file,
            "ref_seq_path": ref_seq_path,
            "enable_harmonize": enable_harmonize,
        },
        "outputs": {
            "harmonized_tsv": harmonized_tsv,
            "log_file": log_path,
        },
        "statistics": {
            "n_variants_original": int(n_variants_original),
            "n_variants_final": int(n_variants_final),
            "n_variants_filtered": int(n_variants_original - n_variants_final),
        },
    }
    
    with open(summary_json, "w", encoding="utf-8") as f:
        json.dump(meta, f, ensure_ascii=False, indent=2)

    dt = time.time() - t0
    logger.info(f"JSON 摘要已写出：{summary_json}")
    logger.info(f"流程结束，总耗时：{dt:.2f} 秒")

    # 准备返回结果
    result = {
        "success": True,
        "plink_stats_file": plink_stats_file,
        "harmonized_tsv": harmonized_tsv,
        "summary_json": summary_json,
        "log_file": log_path,
        "enable_harmonize": enable_harmonize,
        "n_variants_original": n_variants_original,
        "n_variants_final": n_variants_final,
        "n_variants_filtered": n_variants_original - n_variants_final,
        "processing_time": dt,
    }

    # 格式化汇总文本
    summary_lines = [
        "=== 汇总统计文件协调与转换完成 ===",
        f"输入文件：{plink_stats_file}",
        f"输出文件：{harmonized_tsv}",
        f"参考序列：{ref_seq_path}",
        f"启用等位基因协调：{'是' if enable_harmonize else '否'}",
        f"原始变体数：{n_variants_original:,}",
        f"最终变体数：{n_variants_final:,}",
        f"过滤变体数：{n_variants_original - n_variants_final:,}",
        f"处理耗时：{dt:.2f} 秒",
        f"日志文件：{log_path}",
        f"摘要文件：{summary_json}",
    ]
    
    summary_text = "\n".join(summary_lines)
    _echo("\n" + summary_text)
    result["summary_text"] = summary_text

    return result



def merge_harmonized_sumstats(
    wgs_path: str = "cteph_agp3k.wgs.harmonized_sumstat.tsv",
    array_path: str = "cteph_agp3k.array.harmonized_sumstat.tsv",
    output_prefix: Optional[str] = None,
    out_dir: Optional[str] = None,
) -> Dict[str, Any]:
    """
    功能：合并两个 harmonized_sumstat.tsv 文件，按照 SNPID 一致进行内连接合并。
    
    参数说明：
        wgs_path (str): WGS harmonized sumstat 文件路径，默认为 "cteph_agp3k.wgs.harmonized_sumstat.tsv"。
        array_path (str): Array harmonized sumstat 文件路径，默认为 "cteph_agp3k.array.harmonized_sumstat.tsv"。
        output_prefix (Optional[str]): 输出文件前缀。若为 None，则自动生成为 "merged_harmonized_sumstat"。
        out_dir (Optional[str]): 输出目录，若为 None 则使用当前工作目录。
    
    提取字段：
        - SNPID: 变体ID（用于合并）
        - EAF: 效应等位基因频率
        - BETA: beta系数
        - SE: 标准误
        - OR: 比值比
        - OR_95L: OR 95%置信区间下限
        - OR_95U: OR 95%置信区间上限
        - P: P值
        - N: 样本量
    
    输出制品：
        - {output_prefix}.merged_sumstat.tsv    ：合并后的汇总统计文件（Tab 分隔）
        - {output_prefix}.summary.json          ：结构化记录输入参数与输出文件位置
    
    返回：
        Dict[str, Any]：包含输入/输出路径、处理统计信息等的结果字典。
    
    备注：
        - 只保留两个文件中 SNPID 一致的变体（内连接）
        - WGS 数据的列名添加 "_wgs" 后缀，Array 数据的列名添加 "_array" 后缀
        - SNPID 列保持原名，不添加后缀
    """
    import os
    import json
    from datetime import datetime
    from pathlib import Path
    
    _echo("开始任务：合并两个 harmonized sumstat 文件")
    
    # 检查输入文件是否存在
    wgs_file = Path(wgs_path)
    array_file = Path(array_path)
    
    if not wgs_file.exists():
        raise FileNotFoundError(f"WGS sumstat 文件未找到: {wgs_path}")
    if not array_file.exists():
        raise FileNotFoundError(f"Array sumstat 文件未找到: {array_path}")
    
    # 设置输出目录和前缀
    if out_dir is None:
        out_dir_path = Path.cwd()
        _echo(f"未指定 out_dir，默认输出到当前工作目录：{out_dir_path}")
    else:
        out_dir_path = Path(out_dir)
        if not out_dir_path.exists():
            out_dir_path.mkdir(parents=True, exist_ok=True)
        _echo(f"指定输出目录：{out_dir_path}")
    
    if output_prefix is None:
        output_prefix = "merged_harmonized_sumstat"
    
    # 输出文件路径
    merged_tsv = out_dir_path / f"{output_prefix}.merged_sumstat.tsv"
    summary_json = out_dir_path / f"{output_prefix}.summary.json"
    
    # 定义需要提取的列
    required_cols = ["SNPID", "EAF", "BETA", "SE", "OR", "OR_95L", "OR_95U", "P", "N"]
    
    _echo(f"读取 WGS 文件：{wgs_file}")
    try:
        wgs_df = pd.read_csv(wgs_file, sep="\t")
        _echo(f"WGS 文件读取成功，行数：{len(wgs_df)}")
    except Exception as e:
        raise ValueError(f"读取 WGS 文件失败: {e}")
    
    _echo(f"读取 Array 文件：{array_file}")
    try:
        array_df = pd.read_csv(array_file, sep="\t")
        _echo(f"Array 文件读取成功，行数：{len(array_df)}")
    except Exception as e:
        raise ValueError(f"读取 Array 文件失败: {e}")
    
    # 检查必需列是否存在
    wgs_missing = [col for col in required_cols if col not in wgs_df.columns]
    array_missing = [col for col in required_cols if col not in array_df.columns]
    
    if wgs_missing:
        raise ValueError(f"WGS 文件缺少必需列: {wgs_missing}。可用列: {list(wgs_df.columns)}")
    if array_missing:
        raise ValueError(f"Array 文件缺少必需列: {array_missing}。可用列: {list(array_df.columns)}")
    
    # 提取所需列
    wgs_subset = wgs_df[required_cols].copy()
    array_subset = array_df[required_cols].copy()
    
    _echo(f"WGS 提取字段后行数：{len(wgs_subset)}")
    _echo(f"Array 提取字段后行数：{len(array_subset)}")
    
    # 重命名列（除了 SNPID）
    cols_to_rename = [col for col in required_cols if col != "SNPID"]
    
    wgs_rename_dict = {col: f"{col}_wgs" for col in cols_to_rename}
    array_rename_dict = {col: f"{col}_array" for col in cols_to_rename}
    
    wgs_subset = wgs_subset.rename(columns=wgs_rename_dict)
    array_subset = array_subset.rename(columns=array_rename_dict)
    
    _echo("重命名列完成")
    _echo(f"WGS 列名: {list(wgs_subset.columns)}")
    _echo(f"Array 列名: {list(array_subset.columns)}")
    
    # 按 SNPID 进行内连接合并
    _echo("开始按 SNPID 合并数据...")
    merged_df = pd.merge(wgs_subset, array_subset, on="SNPID", how="inner")
    
    n_merged = len(merged_df)
    n_wgs_unique = len(wgs_subset["SNPID"].unique())
    n_array_unique = len(array_subset["SNPID"].unique())
    n_wgs_only = n_wgs_unique - n_merged
    n_array_only = n_array_unique - n_merged
    
    _echo(f"合并完成，结果统计：")
    _echo(f"  WGS 独有变体数：{n_wgs_only}")
    _echo(f"  Array 独有变体数：{n_array_only}")
    _echo(f"  共同变体数（合并后）：{n_merged}")
    
    if n_merged == 0:
        _echo("警告：没有共同的 SNPID，合并结果为空")
    
    # 保存合并结果
    merged_df.to_csv(merged_tsv, sep="\t", index=False)
    _echo(f"合并结果已保存：{merged_tsv}")
    
    # 创建摘要信息
    meta = {
        "created_at": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        "inputs": {
            "wgs_path": str(wgs_file.absolute()),
            "array_path": str(array_file.absolute()),
            "required_columns": required_cols,
        },
        "outputs": {
            "merged_tsv": str(merged_tsv.absolute()),
        },
        "statistics": {
            "n_wgs_variants": int(len(wgs_subset)),
            "n_array_variants": int(len(array_subset)),
            "n_wgs_unique_snpid": int(n_wgs_unique),
            "n_array_unique_snpid": int(n_array_unique),
            "n_merged_variants": int(n_merged),
            "n_wgs_only": int(n_wgs_only),
            "n_array_only": int(n_array_only),
        },
        "column_mapping": {
            "wgs_columns": list(wgs_subset.columns),
            "array_columns": list(array_subset.columns),
            "merged_columns": list(merged_df.columns),
        },
    }
    
    with open(summary_json, "w", encoding="utf-8") as f:
        json.dump(meta, f, ensure_ascii=False, indent=2)
    
    _echo(f"摘要信息已保存：{summary_json}")
    
    # 准备返回结果
    result = {
        "success": True,
        "wgs_path": str(wgs_file.absolute()),
        "array_path": str(array_file.absolute()),
        "merged_tsv": str(merged_tsv.absolute()),
        "summary_json": str(summary_json.absolute()),
        "n_wgs_variants": len(wgs_subset),
        "n_array_variants": len(array_subset),
        "n_merged_variants": n_merged,
        "n_wgs_only": n_wgs_only,
        "n_array_only": n_array_only,
    }
    
    # 格式化汇总文本
    summary_lines = [
        "=== 合并 harmonized sumstat 文件完成 ===",
        f"WGS 输入文件：{wgs_file}",
        f"Array 输入文件：{array_file}",
        f"合并输出文件：{merged_tsv}",
        f"WGS 变体数：{len(wgs_subset):,}",
        f"Array 变体数：{len(array_subset):,}",
        f"WGS 独有变体：{n_wgs_only:,}",
        f"Array 独有变体：{n_array_only:,}",
        f"共同变体数：{n_merged:,}",
        f"摘要文件：{summary_json}",
    ]
    
    summary_text = "\n".join(summary_lines)
    _echo("\n" + summary_text)
    result["summary_text"] = summary_text
    
    return result


def plot_wgs_array_beta_comparison(
    merged_sumstat_file: str = "cteph_agp3k.array_wgs_compration.merged_sumstat.tsv",
    output_prefix: Optional[str] = None,
    out_dir: Optional[str] = None,
    figsize: tuple = (10, 8),
    alpha: float = 0.6,
    point_size: float = 20,
    error_alpha: float = 0.3,
    add_diagonal: bool = True,
    add_correlation: bool = True,
    p_threshold: Optional[float] = None,
    dpi: int = 300,
    save_pdf: bool = True,
    chunk_size: int = 100000,
) -> Dict[str, Any]:
    """
    功能：基于合并的 sumstat 文件绘制 WGS vs Array beta 值比较的学术散点图。
    
    参数说明：
        merged_sumstat_file (str): 合并后的 sumstat 文件路径，默认为 "cteph_agp3k.array_wgs_compration.merged_sumstat.tsv"。
        output_prefix (Optional[str]): 输出图片文件前缀。若为 None，则自动生成为 "wgs_array_beta_comparison"。
        out_dir (Optional[str]): 输出目录，若为 None 则使用当前工作目录。
        figsize (tuple): 图片尺寸 (width, height)，默认 (10, 8)。
        alpha (float): 散点透明度，默认 0.6。
        point_size (float): 散点大小，默认 20。
        error_alpha (float): 误差线透明度，默认 0.3。
        add_diagonal (bool): 是否添加对角线（y=x），默认 True。
        add_correlation (bool): 是否添加相关系数信息，默认 True。
        p_threshold (Optional[float]): P值阈值，仅绘制小于此阈值的点。若为 None 则绘制所有点。
        dpi (int): 图片分辨率，默认 300。
        save_pdf (bool): 是否保存PDF格式文件，默认 True。
        chunk_size (int): 分块读取数据的大小，优化内存使用，默认 100000。
    
    输出制品：
        - {output_prefix}.png                    ：高分辨率 PNG 图片
        - {output_prefix}.pdf                    ：矢量 PDF 图片（当 save_pdf=True 时）
        - {output_prefix}.plot_summary.json      ：绘图统计信息
    
    返回：
        Dict[str, Any]：包含输入/输出路径、绘图统计信息等的结果字典。
    
    依赖：
        - matplotlib
        - seaborn
        - numpy
        - scipy.stats (用于计算相关系数)
    
    备注：
        - 自动处理缺失值和异常值
        - 误差线基于各自的标准误 (SE)
        - 支持 P 值过滤以突出显著性结果
    """
    from datetime import datetime
    from pathlib import Path
    import numpy as np
    import scipy.stats as stats
    import json
    
    # 延迟导入绘图库
    try:
        import matplotlib.pyplot as plt
        import seaborn as sns
    except ImportError:
        raise ImportError("需要安装绘图库：pip install matplotlib seaborn")
    
    _echo("开始任务：绘制 WGS vs Array beta 值比较散点图")
    
    # 检查输入文件是否存在
    input_file = Path(merged_sumstat_file)
    if not input_file.exists():
        raise FileNotFoundError(f"合并的 sumstat 文件未找到: {merged_sumstat_file}")
    
    # 设置输出目录和前缀
    if out_dir is None:
        out_dir_path = Path.cwd()
        _echo(f"未指定 out_dir，默认输出到当前工作目录：{out_dir_path}")
    else:
        out_dir_path = Path(out_dir)
        if not out_dir_path.exists():
            out_dir_path.mkdir(parents=True, exist_ok=True)
        _echo(f"指定输出目录：{out_dir_path}")
    
    if output_prefix is None:
        output_prefix = "wgs_array_beta_comparison"
    
    # 输出文件路径
    png_file = out_dir_path / f"{output_prefix}.png"
    pdf_file = out_dir_path / f"{output_prefix}.pdf" if save_pdf else None
    summary_json = out_dir_path / f"{output_prefix}.plot_summary.json"
    
    # 读取数据（优化内存使用）
    _echo(f"读取合并的 sumstat 文件：{input_file}")
    try:
        # 分块读取以优化内存使用
        _echo("使用分块读取优化内存使用")
        chunks = []
        for chunk in pd.read_csv(input_file, sep="\t", chunksize=chunk_size):
            chunks.append(chunk)
        df = pd.concat(chunks, ignore_index=True)
        del chunks  # 释放内存
        _echo(f"数据读取成功，行数：{len(df)}")
    except Exception as e:
        raise ValueError(f"读取文件失败: {e}")
    
    # 检查必需列是否存在
    required_cols = ["BETA_wgs", "BETA_array", "SE_wgs", "SE_array"]
    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        raise ValueError(f"文件缺少必需列: {missing_cols}。可用列: {list(df.columns)}")
    
    # 立即只保留必需的列以减少内存占用
    available_cols = required_cols.copy()
    if p_threshold is not None:
        p_cols = ["P_wgs", "P_array"]
        available_p_cols = [col for col in p_cols if col in df.columns]
        available_cols.extend(available_p_cols)
    
    df = df[available_cols].copy()
    _echo(f"只保留必需列以优化内存，当前列数：{len(df.columns)}")
    
    # 可选：P值过滤
    if p_threshold is not None:
        # 检查P值列是否存在
        p_cols = ["P_wgs", "P_array"]
        available_p_cols = [col for col in p_cols if col in df.columns]
        
        if available_p_cols:
            # 使用任意一个P值列进行过滤，或者两个都要满足条件
            if len(available_p_cols) == 2:
                # 两个P值列都存在，使用更严格的条件（两个都要小于阈值）
                mask = (df["P_wgs"] < p_threshold) & (df["P_array"] < p_threshold)
                filter_desc = f"P_wgs < {p_threshold} AND P_array < {p_threshold}"
            else:
                # 只有一个P值列存在
                mask = df[available_p_cols[0]] < p_threshold
                filter_desc = f"{available_p_cols[0]} < {p_threshold}"
            
            n_before = len(df)
            df = df[mask].copy()
            n_after = len(df)
            _echo(f"P值过滤 ({filter_desc})：{n_before} -> {n_after} 个变体")
        else:
            _echo(f"警告：未找到P值列 {p_cols}，跳过P值过滤")
    
    # 移除缺失值
    n_before_dropna = len(df)
    df = df.dropna(subset=required_cols)
    n_after_dropna = len(df)
    if n_before_dropna != n_after_dropna:
        _echo(f"移除缺失值：{n_before_dropna} -> {n_after_dropna} 个变体")
    
    if len(df) == 0:
        raise ValueError("经过过滤和缺失值处理后，没有有效数据可以绘图")
    
    # 提取数据并进行类型优化
    _echo("优化数据类型以减少内存占用")
    n_variants = len(df)  # 保存变体数量
    beta_wgs = df["BETA_wgs"].astype(np.float32)
    beta_array = df["BETA_array"].astype(np.float32)
    se_wgs = df["SE_wgs"].astype(np.float32)
    se_array = df["SE_array"].astype(np.float32)
    
    # 释放原始DataFrame内存
    del df
    
    # 计算统计量
    _echo("计算相关系数")
    corr_coef, corr_p = stats.pearsonr(beta_wgs, beta_array)
    
    # 计算数据范围
    wgs_min, wgs_max = float(min(beta_wgs)), float(max(beta_wgs))
    array_min, array_max = float(min(beta_array)), float(max(beta_array))
    
    # 计算数据范围用于设置坐标轴
    beta_range = max(wgs_max, array_max) - min(wgs_min, array_min)
    axis_margin = beta_range * 0.1
    axis_min = min(wgs_min, array_min) - axis_margin
    axis_max = max(wgs_max, array_max) + axis_margin
    
    _echo(f"绘图数据统计：")
    _echo(f"  有效变体数：{n_variants}")
    _echo(f"  WGS Beta 范围：[{wgs_min:.4f}, {wgs_max:.4f}]")
    _echo(f"  Array Beta 范围：[{array_min:.4f}, {array_max:.4f}]")
    _echo(f"  Pearson 相关系数：{corr_coef:.4f} (P = {corr_p:.2e})")
    
    # 设置绘图样式
    plt.style.use('default')
    sns.set_palette("husl")
    
    # 创建图形
    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
    
    # 优化绘图性能：减少误差线密度
    _echo("开始绘制散点图（优化性能）")
    
    # 绘制散点图和误差线（优化版本）
    ax.errorbar(
        beta_wgs, beta_array,
        xerr=se_wgs, yerr=se_array,
        fmt='o',
        markersize=np.sqrt(point_size),
        alpha=alpha,
        ecolor=(0.5, 0.5, 0.5, error_alpha),  # 使用 RGBA 格式设置误差线颜色和透明度
        elinewidth=0.5,
        capsize=0,
        capthick=0.5,
        errorevery=1,
        markeredgewidth=0.5,
        markeredgecolor='white',
        label=f'Variants (n={n_variants:,})',
        rasterized=True  # 栅格化以提高大数据集性能
    )
    
    # 添加对角线
    if add_diagonal:
        ax.plot([axis_min, axis_max], [axis_min, axis_max], 
                'k--', alpha=0.7, linewidth=1.5, label='y = x')
    
    # 设置坐标轴
    ax.set_xlim(axis_min, axis_max)
    ax.set_ylim(axis_min, axis_max)
    ax.set_xlabel('WGS Beta', fontsize=14, fontweight='bold')
    ax.set_ylabel('Array Beta', fontsize=14, fontweight='bold')
    
    # 设置标题
    title = 'WGS vs Array Beta Comparison'
    if p_threshold is not None:
        title += f'\n(P < {p_threshold})'
    ax.set_title(title, fontsize=16, fontweight='bold', pad=20)
    
    # 添加相关系数信息
    if add_correlation:
        corr_text = f'Pearson r = {corr_coef:.4f}\nP = {corr_p:.2e}\nn = {n_variants:,}'
        ax.text(0.05, 0.95, corr_text, transform=ax.transAxes, 
                fontsize=12, verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    # 美化坐标轴
    ax.tick_params(axis='both', which='major', labelsize=12)
    ax.grid(True, alpha=0.3)
    ax.set_aspect('equal', adjustable='box')
    
    # 添加图例
    ax.legend(loc='lower right', fontsize=11)
    
    # 调整布局
    plt.tight_layout()
    
    # 保存图片
    _echo(f"保存 PNG 图片：{png_file}")
    fig.savefig(png_file, dpi=dpi, bbox_inches='tight', facecolor='white')
    
    if save_pdf:
        _echo(f"保存 PDF 图片：{pdf_file}")
        fig.savefig(pdf_file, bbox_inches='tight', facecolor='white')
    
    plt.close(fig)
    
    # 创建摘要信息
    meta = {
        "created_at": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        "inputs": {
            "merged_sumstat_file": str(input_file.absolute()),
            "p_threshold": p_threshold,
            "figsize": figsize,
            "alpha": alpha,
            "point_size": point_size,
            "save_pdf": save_pdf,
        },
        "outputs": {
            "png_file": str(png_file.absolute()),
            "pdf_file": str(pdf_file.absolute()) if save_pdf else None,
        },
        "statistics": {
            "n_variants_plotted": int(n_variants),
            "beta_range_wgs": [wgs_min, wgs_max],
            "beta_range_array": [array_min, array_max],
            "pearson_correlation": float(corr_coef),
            "correlation_p_value": float(corr_p),
        },
        "plot_settings": {
            "add_diagonal": add_diagonal,
            "add_correlation": add_correlation,
            "dpi": dpi,
        },
    }
    
    with open(summary_json, "w", encoding="utf-8") as f:
        json.dump(meta, f, ensure_ascii=False, indent=2)
    
    _echo(f"绘图摘要已保存：{summary_json}")
    
    # 准备返回结果
    result = {
        "success": True,
        "input_file": str(input_file.absolute()),
        "png_file": str(png_file.absolute()),
        "pdf_file": str(pdf_file.absolute()) if save_pdf else None,
        "summary_json": str(summary_json.absolute()),
        "n_variants_plotted": n_variants,
        "pearson_correlation": corr_coef,
        "correlation_p_value": corr_p,
        "beta_range_wgs": [wgs_min, wgs_max],
        "beta_range_array": [array_min, array_max],
        "save_pdf": save_pdf,
    }
    
    # 格式化汇总文本
    summary_lines = [
        "=== WGS vs Array Beta 比较绘图完成 ===",
        f"输入文件：{input_file}",
        f"PNG 输出：{png_file}",
    ]
    
    if save_pdf:
        summary_lines.append(f"PDF 输出：{pdf_file}")
    
    summary_lines.extend([
        f"绘制变体数：{n_variants:,}",
        f"Pearson 相关系数：{corr_coef:.4f} (P = {corr_p:.2e})",
        f"WGS Beta 范围：[{wgs_min:.4f}, {wgs_max:.4f}]",
        f"Array Beta 范围：[{array_min:.4f}, {array_max:.4f}]",
        f"图片分辨率：{dpi} DPI",
        f"摘要文件：{summary_json}",
    ])
    
    if p_threshold is not None:
        summary_lines.insert(4, f"P值过滤阈值：{p_threshold}")
    
    summary_text = "\n".join(summary_lines)
    _echo("\n" + summary_text)
    result["summary_text"] = summary_text
    
    return result
