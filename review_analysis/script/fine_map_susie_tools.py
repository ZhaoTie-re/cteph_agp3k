#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import json
import time
import logging
import subprocess
from datetime import datetime
from typing import Optional, Dict

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt


def build_lead_and_locus_tables(
    plink_stats_file: str,
    range_bp: int = 1_000_000,
    windowsizekb: int = 500,
    sig_level: float = 5e-8,
    ref_seq_path: Optional[str] = "/LARGE0/gr10478/b37974/Pulmonary_Hypertension/nagasaki_pipeline/data/hs38DH.fa",
    enable_harmonize: bool = True,
    output_prefix: Optional[str] = None,
    log_path: Optional[str] = None,
    verbose: Optional[bool] = False,
) -> str:
    """
    功能：基于 GWAS 汇总结果（PLINK2 glm.logistic）自动提取 genome-wide 显著的 lead 变体，并以每个 lead 为中心生成指定碱基范围（range_bp）的 locus 数据表；输出 lead 位点表（TSV）、locus 数据集合（PKL）、运行摘要（JSON）与中文日志。
    
    参数说明：
        plink_stats_file (str): PLINK2 结果文件路径（*.glm.logistic，gwaslab 支持 fmt="plink2"）。
        range_bp (int): 每个 locus 的总碱基范围（对称窗口），例如 1_000_000 表示 ±500kb。
        windowsizekb (int): 判定 lead 变体时的窗口大小（单位 KB），传入 gwaslab.get_lead。
        sig_level (float): genome-wide 显著性阈值，例如 5e-8。
        ref_seq_path (Optional[str]): 参考基因组 fasta 路径（用于 harmonize）。若为 None 或 enable_harmonize=False，则跳过 harmonize。
        enable_harmonize (bool): 是否对 locus 数据进行 harmonize 处理，默认 True。设为 False 时不进行 harmonize，无论 ref_seq_path 是否提供。
        output_prefix (Optional[str]): 输出文件前缀。若仅为文件名不含路径，则默认写入当前工作目录；若为 None，则自动用输入文件名生成 `{basename}.susie_prep` 前缀。
        log_path (Optional[str]): 日志文件路径；若为 None，则自动使用 `{output_prefix}.log`。
        verbose (Optional[bool]): 是否在调用 gwaslab 的相关函数时开启详细输出；默认 False。内部对不支持 verbose 参数的 gwaslab 版本会自动降级为不传该参数。
    
    输出制品：
        - {output_prefix}.lead_variants.tsv     ：lead 位点表（Tab 分隔）
        - {output_prefix}.locus_summaries.pkl   ：dict[str, pandas.DataFrame]，每个 lead 的 locus 数据
        - {output_prefix}.summary.json          ：结构化记录输入参数与输出文件位置
        - {output_prefix}.log                   ：同步中文日志
    
    返回：
        str：`summary.json` 的绝对路径。
    
    依赖：
        - gwaslab（函数内延迟导入）
        - pandas / numpy / matplotlib（已在脚本层导入）
    
    备注：
        - 若未检测到显著 lead，会写出带表头的空 TSV，并写出空的 PKL 与 JSON，便于下游流程兼容。
        - 本函数对 gwaslab 的 `verbose` 参数采用“尽力而为”的兼容策略：若当前 gwaslab 版本的方法不支持该参数，会自动重试为不带该参数调用。
    """
    def _call_with_verbose(fn, *args, **kwargs):
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
    import gwaslab as gl

    # ---------- 准备路径 ----------
    if output_prefix is None:
        base = os.path.splitext(os.path.basename(plink_stats_file))[0]
        output_prefix = os.path.join(os.path.dirname(plink_stats_file), f"{base}.susie_prep")

    lead_tsv = f"{output_prefix}.lead_variants.tsv"
    locus_pkl = f"{output_prefix}.locus_summaries.pkl"
    summary_json = f"{output_prefix}.summary.json"
    if log_path is None:
        log_path = f"{output_prefix}.log"

    output_dir = os.path.dirname(output_prefix)
    if output_dir == "":
        output_dir = os.getcwd()
        output_prefix = os.path.join(output_dir, output_prefix)
    os.makedirs(output_dir, exist_ok=True)

    # ---------- 设置日志（中文） ----------
    logger = logging.getLogger("susie_prep")
    logger.setLevel(logging.INFO)
    # 清除重复 handler（避免在交互式环境重复添加）
    if logger.handlers:
        for h in list(logger.handlers):
            logger.removeHandler(h)
    fh = logging.FileHandler(log_path, mode="w", encoding="utf-8")
    fmt = logging.Formatter("[%(asctime)s] %(levelname)s - %(message)s")
    fh.setFormatter(fmt)
    logger.addHandler(fh)

    logger.info("启动流程：SuSiE 精细定位前置（lead与locus生成）")
    logger.info(f"输入文件：{plink_stats_file}")
    logger.info(f"参数：range_bp={range_bp}；windowsizekb={windowsizekb}；sig_level={sig_level}")
    logger.info(f"参考序列（用于harmonize）：{ref_seq_path}")
    logger.info(f"是否启用harmonize：{enable_harmonize}")

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

    # ---------- 获取 lead 变体 ----------
    logger.info("开始识别 lead 变体")
    lead_df = _call_with_verbose(sumstats.get_lead, windowsizekb=windowsizekb, sig_level=sig_level)
    n_lead = 0 if lead_df is None else len(lead_df) # type: ignore
    logger.info(f"识别到 {n_lead} 个 lead 变体")
    # 写出 lead 表
    if lead_df is None or lead_df.empty: # type: ignore
        # 仍写一个空表头，避免下游出错
        empty_cols = ["CHR", "POS", "SNPID"]  # 常见字段，具体以 gwaslab 输出为准
        pd.DataFrame(columns=empty_cols).to_csv(lead_tsv, sep="\t", index=False)
        logger.warning("未检测到显著 lead 变体，已写出空的 TSV 表头")
        locus_summaries = {}
    else:
        lead_df.to_csv(lead_tsv, sep="\t", index=False) # type: ignore
        logger.info(f"lead 变体已写出：{lead_tsv}")

        # ---------- 为每个 lead 构建 locus ----------
        locus_summaries: Dict[str, pd.DataFrame] = {}
        half = range_bp // 2

        for idx, row in lead_df.iterrows(): # type: ignore
            chr_val = row["CHR"]
            pos_val = row["POS"]
            snpid = row.get("SNPID", f"CHR{chr_val}:{pos_val}")

            start_pos = int(pos_val) - half
            end_pos = int(pos_val) + half
            # gwaslab 的 filter_value 语法：使用 pandas query 表达式
            filter_condition = f"CHR=={chr_val} & POS>{start_pos} & POS<{end_pos}"

            logger.info(f"构建 locus：{snpid}；区间：{chr_val}:{start_pos}-{end_pos}")
            locus = _call_with_verbose(sumstats.filter_value, filter_condition)

            # 填充与协调等位方向
            _call_with_verbose(locus.fill_data, to_fill=["BETA"]) # type: ignore
            
            # 根据参数决定是否进行 harmonize
            if enable_harmonize and ref_seq_path is not None:
                _call_with_verbose(locus.harmonize, basic_check=False, ref_seq=ref_seq_path) # type: ignore
                logger.info(f"已对 locus {snpid} 完成 harmonize 处理")
            else:
                logger.info(f"跳过 locus {snpid} 的 harmonize 处理 (enable_harmonize={enable_harmonize}, ref_seq_path={ref_seq_path})")

            # 保存 DataFrame
            locus_summaries[str(snpid)] = locus.data # type: ignore
            logger.info(f"完成处理 locus：{snpid}（记录数：{len(locus.data)})") # type: ignore

    # ---------- 序列化 locus 字典 ----------
    pd.to_pickle(locus_summaries, locus_pkl) # type: ignore
    logger.info(f"locus_summaries 已写出（pickle）：{locus_pkl}")

    # ---------- 写出与返回 JSON 摘要 ----------
    meta = {
        "created_at": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        "inputs": {
            "plink_stats_file": plink_stats_file,
            "range_bp": int(range_bp),
            "windowsizekb": int(windowsizekb),
            "sig_level": float(sig_level),
            "ref_seq_path": ref_seq_path,
            "enable_harmonize": bool(enable_harmonize),
        },
        "outputs": {
            "lead_tsv": lead_tsv,
            "locus_pkl": locus_pkl,
            "log_file": log_path,
        },
    }
    with open(summary_json, "w", encoding="utf-8") as f:
        json.dump(meta, f, ensure_ascii=False, indent=2)

    dt = time.time() - t0
    logger.info(f"JSON 摘要已写出：{summary_json}")
    logger.info(f"流程结束，总耗时：{dt:.2f} 秒")

    return summary_json



def build_ld_matrices_from_summary(
    summary_json_path: str,
    bed_prefix: str = "/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/wgs/19.tommo_panel_filter/cteph_agp3k.lowfreq_common",
    plink_path: str = "/home/b/b37974/plink",
    out_dir: Optional[str] = None,
    strict: bool = False,
    id_column: str = "SNPID",
    logger: Optional[logging.Logger] = None,
    output_prefix: Optional[str] = None,
    case_prefix: str = "PHOM",
    sample_mode: str = "all",
) -> str:
    """
    功能：
        基于由 `build_lead_and_locus_tables` 生成的 summary.json，逐个 lead（SNPID）
        提取其对应的 sum_stat DataFrame（来自 locus_pkl），并调用 PLINK 对该 sum_stat 中的
        变体集合计算：
          1) LD 相关系数矩阵（R；--r square）
          2) LD 决定系数矩阵（R^2；--r2 square）
        最终为每个 lead 写出三份 TSV（sum_stat / R / R2），并汇总生成新的 JSON 后返回其路径。

    矩阵标签说明：
        默认将每个矩阵（R 与 R²）的行名/列名标注为该 lead 对应的 `extract` 列表中的 SNPID（按行写出的顺序）。
        若提取列表数量与实际矩阵阶数不一致，则尝试按参考面板 `.bim` 的顺序与 `extract` 求交来推断最终顺序；
        若仍不一致，则降级为使用 0-based 整数索引。

    参数：
        summary_json_path (str): `build_lead_and_locus_tables` 产生的 JSON 路径。
        bed_prefix (str): 参考面板的 PLINK 二进制前缀（不含扩展名），用于计算 LD。
        plink_path (str): plink 可执行文件路径。
        out_dir (Optional[str]): 输出目录。默认为当前工作目录（os.getcwd()）。
        strict (bool): 若 True，遇到 lead 在 pkl 中缺失时立即抛错；默认 False 则跳过并记录。
        id_column (str): 在 sum_stat DataFrame 中作为变体 ID 的列名，默认 "SNPID"。
        logger (Optional[logging.Logger]): 可选日志对象；若不提供，使用 print 输出关键信息。
        output_prefix (Optional[str]): 控制输出 JSON 文件名的前缀；若为空则使用默认文件名。
        case_prefix (str): 病例样本 ID 的前缀（用于 .fam 的 IID 前缀匹配），默认 "PHOM"。
        sample_mode (str): 样本选择模式，"all"（默认，使用所有样本）|"case"（仅病例，IID 以 case_prefix 开头）|"ctrl"（仅对照，IID 不以 case_prefix 开头）。

    返回：
        str：新生成的汇总 JSON 文件路径。

    说明与注意事项：
        - 路径解析：若 `lead_tsv` / `locus_pkl` 为相对路径，则以 summary_json 所在目录为基准。
        - 变体匹配：`--extract` 使用的 ID 必须与参考面板的 .bim 第二列一致；若不完全一致，
          PLINK 会自动忽略缺失 ID。本函数会统计保留数量。
        - 输出文件：
            {out_dir}/{lead_sanitized}.sum_stat.tsv
            {out_dir}/{lead_sanitized}.ld_r.tsv
            {out_dir}/{lead_sanitized}.ld_r2.tsv
          其中 lead_sanitized 会将不安全字符替换为下划线。
        - PLINK 输出解析：使用 `--r square` / `--r2 square` 时，PLINK 会生成 `.ld` 矩阵文件。
          本函数直接以整数索引顺序读取，不再依赖 .ld.id 文件。
        - 汇总 JSON 文件名：当提供 output_prefix 时，写为
          {out_dir}/{output_prefix}.ld_matrices_by_lead.summary.json；
          否则为 {out_dir}/ld_matrices_by_lead.summary.json
        - 样本选择：依据 {bed_prefix}.fam 中的 IID 与 case_prefix 前缀匹配实现；当 sample_mode != "all" 时将通过 `--keep` 传递给 PLINK。会在汇总 JSON 中记录模式与纳入样本数。
    """
    # ---------- 帮助函数 ----------
    def _log(msg: str):
        if logger is not None:
            logger.info(msg)
        else:
            print("[build_ld_matrices_from_summary]", msg)
    
    def _resolve(base: str, p: str) -> str:
        return p if os.path.isabs(p) else os.path.join(base, p)
    
    def _sanitize(name: str) -> str:
        # 仅保留常见安全字符
        return "".join(ch if ch.isalnum() or ch in ("_", "-", ".", ":") else "_" for ch in str(name))
    
    def _read_plink_square_matrix(matrix_path: str) -> pd.DataFrame:
        """
        读取 PLINK `--r square` 或 `--r2 square` 生成的矩阵文件，并返回以整数索引为标签的 DataFrame。
        当前版本不再读取 .ld.id 文件，所有输出均以 0-based 整数索引编号。
        """
        if not os.path.exists(matrix_path):
            raise FileNotFoundError(f"缺少 PLINK 矩阵文件：{matrix_path}")
        mat = np.loadtxt(matrix_path)
        if mat.ndim == 1:
            mat = np.array(mat, ndmin=2)
        n = mat.shape[0]
        idx = list(range(n))
        df = pd.DataFrame(mat, index=idx, columns=idx)
        return df
    
    # ---------- 读入 summary.json ----------
    base_dir = os.path.dirname(os.path.abspath(summary_json_path))
    with open(summary_json_path, "r", encoding="utf-8") as f:
        meta = json.load(f)
    outputs = meta.get("outputs", {})
    lead_tsv_path = outputs.get("lead_tsv")
    locus_pkl_path = outputs.get("locus_pkl")
    if lead_tsv_path is None or locus_pkl_path is None:
        raise ValueError("summary.json 缺少必须的 outputs.lead_tsv 或 outputs.locus_pkl 字段。")
    
    lead_tsv_path = _resolve(base_dir, lead_tsv_path)
    locus_pkl_path = _resolve(base_dir, locus_pkl_path)
    if out_dir is None or out_dir == "":
        out_dir = os.getcwd()
    os.makedirs(out_dir, exist_ok=True)
    # 统一目录策略：
    # - 汇总 JSON 写在工作目录 out_dir
    # - 其他中间/明细输出写在 out_dir/tmp
    tmp_dir = os.path.join(out_dir, "tmp")
    os.makedirs(tmp_dir, exist_ok=True)

    # ---------- 样本选择（case/ctrl/all） ----------
    fam_path = f"{bed_prefix}.fam"
    if not os.path.exists(fam_path):
        raise FileNotFoundError(f"缺少 FAM 文件：{fam_path}")
    fam_pairs = []
    with open(fam_path, "r") as ff:
        for line in ff:
            parts = line.rstrip("\n").split()
            if len(parts) >= 2:
                fam_pairs.append((parts[0], parts[1]))  # (FID, IID)
    n_total_samples = len(fam_pairs)
    kept_pairs = fam_pairs
    keep_path = None

    smode = str(sample_mode).lower().strip()
    if smode not in ("all", "case", "ctrl"):
        raise ValueError(f"sample_mode 必须是 'all'|'case'|'ctrl' 之一，当前：{sample_mode}")

    if smode == "case":
        kept_pairs = [p for p in fam_pairs if p[1].startswith(case_prefix)]
    elif smode == "ctrl":
        kept_pairs = [p for p in fam_pairs if not p[1].startswith(case_prefix)]

    n_kept = len(kept_pairs)
    if smode != "all":
        if n_kept == 0:
            raise ValueError(f"样本选择结果为空：sample_mode={sample_mode}, case_prefix={case_prefix}")
        keep_path = os.path.join(tmp_dir, f"samples.keep.{smode}.txt")
        with open(keep_path, "w") as kf:
            for fid, iid in kept_pairs:
                kf.write(f"{fid}\t{iid}\n")
    _log(f"样本选择：mode={smode}，case_prefix={case_prefix}，纳入 {n_kept}/{n_total_samples} 个样本")
    
    # ---------- 载入 lead 表与 locus 字典 ----------
    lead_df = pd.read_csv(lead_tsv_path, sep="\t", dtype={id_column: str})
    if id_column not in lead_df.columns:
        raise ValueError(f"lead_tsv 缺少 {id_column} 列：{lead_tsv_path}")
    locus_dict = pd.read_pickle(locus_pkl_path)
    if not isinstance(locus_dict, dict):
        raise TypeError(f"locus_pkl 应为 dict[str, DataFrame]，实际类型：{type(locus_dict)}")
    
    # ---------- 遍历每个 lead ----------
    out_index = {
        "created_at": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        "source_summary_json": os.path.abspath(summary_json_path),
        "inputs": {
            "bed_prefix": bed_prefix,
            "plink_path": plink_path,
            "out_dir": os.path.abspath(out_dir),
            "id_column": id_column,
            "output_prefix": output_prefix,
            "case_prefix": case_prefix,
            "sample_mode": smode,
        },
        "sample_selection": {
            "mode": smode,
            "case_prefix": case_prefix,
            "n_total": n_total_samples,
            "n_kept": n_kept,
            "keep_file": keep_path,
        },
        "per_lead_outputs": {},
    }
    _log(f"输出目录策略：汇总 JSON -> {out_dir}；明细与PLINK产物 -> {tmp_dir}")
    
    def _load_extract_ids(path: str) -> list:
        ids = []
        with open(path, "r") as f:
            for line in f:
                s = line.strip()
                if s:
                    ids.append(s)
        return ids

    def _infer_labels_from_bim(bed_prefix: str, extract_ids: list, target_n: int) -> Optional[list]:
        """
        当提取列表数量与矩阵阶数不一致时，按 .bim 的变体顺序与 extract 列表求交，
        以期获得与 PLINK 实际方阵相同的保留顺序；仅当数量恰好等于 target_n 时返回该列表，否则返回 None。
        """
        bim_path = f"{bed_prefix}.bim"
        if not os.path.exists(bim_path):
            return None
        extract_set = set(extract_ids)
        kept = []
        with open(bim_path, "r") as bf:
            for line in bf:
                parts = line.rstrip("\n").split()
                if len(parts) >= 2 and parts[1] in extract_set:
                    kept.append(parts[1])
        if len(kept) == target_n:
            return kept
        return None

    for _, row in lead_df.iterrows():
        lead_id_original = str(row[id_column])
        lead_id = lead_id_original
        sum_stat = locus_dict.get(lead_id)
        # 尝试兼容 key 形式
        if sum_stat is None:
            chr_val = row.get("CHR")
            pos_val = row.get("POS")
            if pd.notna(chr_val) and pd.notna(pos_val):
                for k in (f"CHR{chr_val}:{int(pos_val)}", f"{chr_val}:{int(pos_val)}"):
                    if k in locus_dict:
                        sum_stat = locus_dict[k]
                        lead_id = k
                        break
        if sum_stat is None:
            msg = f"跳过：在 pkl 中找不到该 lead 的 sum_stat -> {lead_id_original}"
            if strict:
                raise KeyError(msg)
            _log(msg)
            continue
        if not isinstance(sum_stat, pd.DataFrame):
            raise TypeError(f"sum_stat 类型异常（非 DataFrame）：{type(sum_stat)} for lead {lead_id}")
    
        # 变体列表
        if id_column not in sum_stat.columns:
            raise ValueError(f"sum_stat 缺少 {id_column} 列 for lead {lead_id}")
        variant_ids = sum_stat[id_column].astype(str).dropna().unique().tolist()
        if len(variant_ids) == 0:
            _log(f"lead {lead_id} 的变体集合为空，跳过。")
            continue
    
        lead_tag = _sanitize(lead_id)
        prefix = os.path.join(tmp_dir, lead_tag)
    
        # 写出 sum_stat
        sum_stat_tsv = f"{prefix}.sum_stat.tsv"
        sum_stat.to_csv(sum_stat_tsv, sep="\t", index=False)
    
        # 写出提取列表
        extract_path = f"{prefix}.extract.txt"
        with open(extract_path, "w") as ef:
            ef.write("\n".join(variant_ids) + "\n")
    
        # 调用 PLINK 计算 R
        r_prefix = f"{prefix}.r"
        cmd_r = [
            plink_path,
            "--bfile", bed_prefix,
            "--keep-allele-order",
            "--r", "square",
            "--extract", extract_path,
            "--out", r_prefix,
        ]
        if keep_path is not None:
            cmd_r.extend(["--keep", keep_path])
        _log(f"运行：{' '.join(cmd_r)}")
        r_log_path = f"{r_prefix}.plink.log"
        with open(r_log_path, "w") as rlog:
            subprocess.run(cmd_r, check=True, stdout=rlog, stderr=subprocess.STDOUT)
        r_ld_path = f"{r_prefix}.ld"
        if not os.path.exists(r_ld_path):
            _log(f"警告：未生成 {r_ld_path}，可能是提取后无有效变体；跳过该 lead。")
            continue
        try:
            r_df = _read_plink_square_matrix(r_ld_path)
        except FileNotFoundError as e:
            _log(f"警告：{e}；该 lead 将跳过 R/R2 输出。")
            continue
        # 为 R 矩阵设置标签：优先使用 extract 列表；必要时尝试 .bim 对齐
        extract_ids = _load_extract_ids(extract_path)
        if len(extract_ids) == r_df.shape[0]:
            r_df.index = extract_ids #type: ignore
            r_df.columns = extract_ids
        else:
            _log(f"警告：R 阶数({r_df.shape[0]})与 extract 数量({len(extract_ids)})不一致，尝试基于 .bim 对齐。")
            inferred = _infer_labels_from_bim(bed_prefix, extract_ids, r_df.shape[0])
            if inferred is not None:
                r_df.index = inferred #type: ignore
                r_df.columns = inferred
            else:
                _log("警告：.bim 对齐未获得一致数量，保留整数索引。")

        r_tsv_path = f"{prefix}.ld_r.tsv"
        r_df.to_csv(r_tsv_path, sep="\t", index=True)

        # 调用 PLINK 计算 R2
        r2_prefix = f"{prefix}.r2"
        cmd_r2 = [
            plink_path,
            "--bfile", bed_prefix,
            "--keep-allele-order",
            "--r2", "square",
            "--extract", extract_path,
            "--out", r2_prefix,
        ]
        if keep_path is not None:
            cmd_r2.extend(["--keep", keep_path])
        _log(f"运行：{' '.join(cmd_r2)}")
        r2_log_path = f"{r2_prefix}.plink.log"
        with open(r2_log_path, "w") as r2log:
            subprocess.run(cmd_r2, check=True, stdout=r2log, stderr=subprocess.STDOUT)
        r2_ld_path = f"{r2_prefix}.ld"
        if not os.path.exists(r2_ld_path):
            _log(f"警告：未生成 {r2_ld_path}（通常由于有效变体数量不足或被过滤）；本次仅输出 R。")
        try:
            r2_df = _read_plink_square_matrix(r2_ld_path)
        except FileNotFoundError as e:
            _log(f"警告：{e}；仅输出 R。")
            r2_df = None

        # 保存 R2（若成功）
        r2_tsv_path = None
        if r2_df is not None:
            # 为 R² 矩阵设置标签：与 R 相同策略
            extract_ids = _load_extract_ids(extract_path)
            if len(extract_ids) == r2_df.shape[0]:
                r2_df.index = extract_ids #type: ignore
                r2_df.columns = extract_ids
            else:
                _log(f"警告：R² 阶数({r2_df.shape[0]})与 extract 数量({len(extract_ids)})不一致，尝试基于 .bim 对齐。")
                inferred = _infer_labels_from_bim(bed_prefix, extract_ids, r2_df.shape[0])
                if inferred is not None:
                    r2_df.index = inferred
                    r2_df.columns = inferred
                else:
                    _log("警告：.bim 对齐未获得一致数量，保留整数索引。")

            r2_tsv_path = f"{prefix}.ld_r2.tsv"
            r2_df.to_csv(r2_tsv_path, sep="\t", index=True)

        # 汇总记录
        out_index["per_lead_outputs"][lead_id] = {
            "sum_stat_tsv": sum_stat_tsv,
            "ld_r_tsv": r_tsv_path,
            "ld_r2_tsv": r2_tsv_path,
            "extract_list": extract_path,
            "plink_outputs": {
                "r_matrix_ld": r_ld_path,
                "r2_matrix_ld": r2_ld_path
            },
            "n_variants_requested": len(variant_ids),
            "n_variants_matrix": r_df.shape[0],
        }
    
    # 写出汇总 JSON
    if output_prefix is None or str(output_prefix).strip() == "":
        final_json = os.path.join(out_dir, "ld_matrices_by_lead.summary.json")
    else:
        final_json = os.path.join(out_dir, f"{output_prefix}.ld_matrices_by_lead.summary.json")
    with open(final_json, "w", encoding="utf-8") as jf:
        json.dump(out_index, jf, ensure_ascii=False, indent=2)
    _log(f"完成：输出索引 JSON -> {final_json}")
    return final_json


def plot_ld_heatmaps_from_index(
    index_json_path: str,
    out_pdf_path: Optional[str] = None,
    threads: int = 8,
    dpi: int = 150,
    show_lead_cross: bool = True,
    output_prefix: Optional[str] = None,
) -> str:
    """
    基于 `build_ld_matrices_from_summary` 产出的汇总 JSON（index_json_path），
    为每个 lead 生成一页 PDF：左侧为 R（相关系数，范围 [-1, 1]）热图，右侧为 R²（决定系数，范围 [0, 1]）热图。
    若某 lead 缺少 R²，则该页仅绘制 R 并在标题中标注“无 R²”。

    参数：
        index_json_path (str): `ld_matrices_by_lead.summary.json` 的路径（或带前缀版本，如 `{prefix}.ld_matrices_by_lead.summary.json`）。
        out_pdf_path (Optional[str]): 输出 PDF 路径；若为空则与 JSON 同目录同前缀，名为 `{basename}.ld_heatmaps.pdf`。
        output_prefix (Optional[str]): 输出 PDF 文件名前缀。若提供，则默认输出为 `{output_prefix}.ld_heatmaps.pdf`。
        threads (int): 并行加载与预处理的最大线程数，上限 8（绘图在主进程顺序执行以避免 matplotlib 并发问题）。
        dpi (int): 每页图像的渲染分辨率 DPI。
        show_lead_cross (bool): 是否绘制 lead 变体的交叉引导线（默认 True 为绘制；False 不绘制）。

    返回：
        str：输出 PDF 的绝对路径。
    """
    # 读取索引 JSON
    with open(index_json_path, "r", encoding="utf-8") as f:
        index_meta = json.load(f)
    per_lead = index_meta.get("per_lead_outputs", {})
    if not isinstance(per_lead, dict) or len(per_lead) == 0:
        raise ValueError("索引 JSON 中缺少 per_lead_outputs 或为空。")

    # 使用默认样式
    plt.style.use("default")

    # 解析输出 PDF 路径
    if output_prefix is not None and (out_pdf_path is None or str(out_pdf_path).strip() == ""):
        base_dir = os.path.dirname(os.path.abspath(index_json_path))
        out_pdf_path = os.path.join(base_dir, f"{output_prefix}.ld_heatmaps.pdf")
    elif out_pdf_path is None or str(out_pdf_path).strip() == "":
        base_dir = os.path.dirname(os.path.abspath(index_json_path))
        base_name = os.path.splitext(os.path.basename(index_json_path))[0]
        out_pdf_path = os.path.join(base_dir, f"{base_name}.ld_heatmaps.pdf")
    out_pdf_path = os.path.abspath(out_pdf_path)

    # 并行加载矩阵（仅做 I/O 与数组预处理；绘图在主线程）
    from concurrent.futures import ThreadPoolExecutor, as_completed
    max_workers = max(1, min(int(threads), 8))

    # ---------- helpers: 染色体/位置排序 ----------
    def _chr_rank(ch: str) -> int:
        """将 chr 字符串映射为自然排序数值：1..22->1..22, X->23, Y->24, M/MT->25, 其他->9999"""
        s = str(ch).strip()
        if s.upper().startswith("CHR"):
            s = s[3:]
        s_up = s.upper()
        if s_up in ("X",):
            return 23
        if s_up in ("Y",):
            return 24
        if s_up in ("M", "MT"):
            return 25
        try:
            val = int(s)
            if 1 <= val <= 22:
                return val
            return 9999
        except Exception:
            return 9999

    def _parse_lead_sig(sig: str):
        """解析 'CHROM:POS:REF:ALT' 或相近格式；返回 (chr_rank, pos)；失败返回 (9999, inf)"""
        if sig is None:
            return (9999, float("inf"))
        parts = str(sig).split(":")
        if len(parts) < 2:
            return (9999, float("inf"))
        chr_key = _chr_rank(parts[0])
        try:
            pos = int(parts[1])
        except Exception:
            pos = float("inf")
        return (chr_key, pos)

    # ---------- helpers: 定位 lead 在标签列表中的索引，并在热图上标记 ----------
    def _normalize_chr(s: str) -> str:
        s = s.strip()
        if s.upper().startswith("CHR"):
            return s[3:]
        return s

    def _locate_lead_index(labels: list, lead: str) -> Optional[int]:
        """在标签列表中查找 lead 的位置。支持以下匹配：
        1) 完整 ID 精确匹配（如 CHR3:154069965:A:G）
        2) 仅按 CHR:POS 匹配（忽略 chr/CHR 前缀差异）
        """
        if labels is None or len(labels) == 0:
            return None
        # 1) 精确匹配
        try:
            return labels.index(lead)
        except ValueError:
            pass
        # 2) 位置匹配
        #   解析 lead 的 chr 和 pos
        parts = str(lead).split(":")
        if len(parts) >= 2:
            lead_chr = _normalize_chr(parts[0])
            try:
                lead_pos = int(parts[1])
            except Exception:
                lead_pos = None
            if lead_pos is not None:
                for i, lab in enumerate(labels):
                    p = str(lab).split(":")
                    if len(p) >= 2:
                        if _normalize_chr(p[0]) == lead_chr:
                            try:
                                if int(p[1]) == lead_pos:
                                    return i
                            except Exception:
                                continue
        return None

    def _mark_lead(ax: plt.Axes, idx: Optional[int], n: int): # type: ignore
        """在热图上画出 lead 的位置十字标记（行/列）"""
        if idx is None or idx < 0 or idx >= n:
            return
        ax.axhline(y=idx, color="k", linewidth=1.2, alpha=0.9)
        ax.axvline(x=idx, color="k", linewidth=1.2, alpha=0.9)
        # 叠加一层浅色线增强可见性（兼容深色背景配色）
        ax.axhline(y=idx, color="w", linewidth=0.6, alpha=0.7)
        ax.axvline(x=idx, color="w", linewidth=0.6, alpha=0.7)

    def _read_matrix(tsv_path: str):
        if tsv_path is None:
            return None, None, None
        if not os.path.exists(tsv_path):
            return None, None, f"文件缺失：{tsv_path}"
        try:
            # 以第一列为索引读取；若文件为方阵且带有 SNPID 行列标签，可直接读取为 DataFrame
            df = pd.read_csv(tsv_path, sep="\t", index_col=0)
            # 确保为数值类型（有时 read_csv 可能把列名解析为字符串）
            df = df.apply(pd.to_numeric, errors="coerce")
            arr = df.to_numpy(dtype=float)
            labels = df.index.tolist()
            # 简单一致性检查（方阵）
            if arr.shape[0] != arr.shape[1]:
                return None, None, f"非方阵：{tsv_path} -> 形状 {arr.shape}"
            return (arr, labels, None)
        except Exception as e:
            return None, None, f"读取失败：{tsv_path}；原因：{e}"

    def _job(lead_id: str, lead_rec: dict):
        r_path = lead_rec.get("ld_r_tsv")
        r2_path = lead_rec.get("ld_r2_tsv")
        r_arr, r_labels, r_err = _read_matrix(r_path) # type: ignore
        r2_arr, r2_labels, r2_err = _read_matrix(r2_path) if r2_path else (None, None, None)
        return {
            "lead_id": lead_id,
            "r": r_arr,
            "r_labels": r_labels,
            "r_err": r_err,
            "r2": r2_arr,
            "r2_labels": r2_labels,
            "r2_err": r2_err,
        }

    results = []
    with ThreadPoolExecutor(max_workers=max_workers) as ex:
        future_map = {ex.submit(_job, lid, rec): (lid, rec) for lid, rec in per_lead.items()}
        for fut in as_completed(future_map):
            results.append(fut.result())

    # 按基因组顺序排序：chr1..22, X, Y, M；POS 升序；兼容不带 'chr' 的形式
    sort_items = []
    for res in results:
        lead_id = res["lead_id"]
        lead_rec = per_lead.get(lead_id, {})
        # 优先使用 per_lead 的 lead_sig（已知格式 CHROM:POS:REF:ALT）
        lead_sig = lead_rec.get("lead_sig", None)
        # 若缺失则退回用 lead_id 自身尝试解析
        sig_for_sort = lead_sig if lead_sig else lead_id
        chr_key, pos = _parse_lead_sig(sig_for_sort)
        sort_items.append((chr_key, pos, str(lead_id), res))
    sort_items.sort(key=lambda t: (t[0], t[1], t[2]))
    results = [t[3] for t in sort_items]

    # 逐页绘图（主进程，避免 matplotlib 并发问题）
    from matplotlib.backends.backend_pdf import PdfPages
    with PdfPages(out_pdf_path) as pdf:
        for res in results:
            lead_id = res["lead_id"]
            r_arr, r2_arr = res["r"], res["r2"]
            r_err, r2_err = res["r_err"], res["r2_err"]

            if r_arr is None and r2_arr is None:
                # 两者都缺失，跳过该 lead
                continue

            # 决定本页 subplot 列数
            ncols = 2 if r_arr is not None and r2_arr is not None else 1
            fig_w = 10 if ncols == 2 else 6
            fig_h = 5
            fig, axes = plt.subplots(1, ncols, figsize=(fig_w, fig_h), dpi=dpi)
            # 两个分图之间增加空隙
            fig.subplots_adjust(wspace=0.35, hspace=0.1)
            if ncols == 1:
                axes = [axes]

            # 绘制 R
            if r_arr is not None:
                ax = axes[0]
                im = ax.imshow(r_arr, vmin=-1.0, vmax=1.0, interpolation="nearest", aspect="equal")
                ax.set_aspect('equal', adjustable='box')
                try:
                    ax.set_box_aspect(1)
                except Exception:
                    pass
                # 标记 lead 在热图中的位置
                lead_idx = _locate_lead_index(res.get("r_labels"), lead_id)
                if show_lead_cross:
                    _mark_lead(ax, lead_idx, r_arr.shape[0])
                ax.set_title(f"{lead_id} — R (n={r_arr.shape[0]})" + (f"\n注意：{r_err}" if r_err else ""))
                ax.set_xlabel("Variants")
                ax.set_ylabel("Variants")
                # 大矩阵时隐藏刻度标签以提升绘制速度
                ax.set_xticks([])
                ax.set_yticks([])
                cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
                cbar.set_label("LD R")
            else:
                axes[0].set_visible(False)

            # 绘制 R²
            if ncols == 2 and r2_arr is not None:
                ax = axes[1]
                im2 = ax.imshow(r2_arr, vmin=0.0, vmax=1.0, interpolation="nearest", aspect="equal")
                ax.set_aspect('equal', adjustable='box')
                try:
                    ax.set_box_aspect(1)
                except Exception:
                    pass
                # 标记 lead 在热图中的位置
                lead_idx2 = _locate_lead_index(res.get("r2_labels"), lead_id)
                if show_lead_cross:
                    _mark_lead(ax, lead_idx2, r2_arr.shape[0])
                ax.set_title(f"{lead_id} — R² (n={r2_arr.shape[0]})" + (f"\n注意：{r2_err}" if r2_err else ""))
                ax.set_xlabel("Variants")
                ax.set_ylabel("Variants")
                ax.set_xticks([])
                ax.set_yticks([])
                cbar2 = fig.colorbar(im2, ax=ax, fraction=0.046, pad=0.04)
                cbar2.set_label("LD R²")
            elif ncols == 1 and r2_arr is not None and r_arr is None:
                # 只有 R² 的情况
                ax = axes[0]
                im2 = ax.imshow(r2_arr, vmin=0.0, vmax=1.0, interpolation="nearest", aspect="equal")
                ax.set_aspect('equal', adjustable='box')
                try:
                    ax.set_box_aspect(1)
                except Exception:
                    pass
                # 标记 lead 在热图中的位置
                lead_idx2 = _locate_lead_index(res.get("r2_labels"), lead_id)
                if show_lead_cross:
                    _mark_lead(ax, lead_idx2, r2_arr.shape[0])
                ax.set_title(f"{lead_id} — R² (n={r2_arr.shape[0]})" + (f"\n注意：{r2_err}" if r2_err else ""))
                ax.set_xlabel("Variants")
                ax.set_ylabel("Variants")
                ax.set_xticks([])
                ax.set_yticks([])
                cbar2 = fig.colorbar(im2, ax=ax, fraction=0.046, pad=0.04)
                cbar2.set_label("LD R²")
            elif r2_arr is None and r_arr is not None:
                # 缺少 R²，给出提示
                axes[0].set_title(f"{lead_id} — 仅有 R；无 R² 输出" + (f"\n注意：{r_err}" if r_err else ""))

            fig.tight_layout()
            pdf.savefig(fig)
            plt.close(fig)

    return out_pdf_path


def integrate_susie_results_with_sumstat(
    susie_summary_json_path: str,
    ld_matrices_summary_json_path: str,
    output_dir: Optional[str] = None,
    output_prefix: str = "integrated_susie_results",
    logger: Optional[logging.Logger] = None,
) -> dict:
    """
    整合 SuSiE 结果与汇总统计量数据，并提取可信集信息。

    参数：
        susie_summary_json_path (str): susie_summary.json 文件路径
        ld_matrices_summary_json_path (str): cteph_agp3k.lowfreq_common.ld_matrices_by_lead.summary.json 文件路径
        output_dir (Optional[str]): 输出目录，默认为 ld_matrices_summary_json 同目录
        output_prefix (str): 输出文件前缀，默认为 "integrated_susie_results"
        logger (Optional[logging.Logger]): 日志记录器

    返回：
        dict: 包含输出文件路径和处理状态的字典
    """
    # ---------- 辅助函数 ----------
    def _log(msg: str):
        if logger:
            logger.info(msg)
        else:
            print(f"[INFO] {msg}")

    def _sanitize_filename(name: str) -> str:
        """将文件名中的特殊字符替换为下划线"""
        import re
        return re.sub(r'[^\w\-_\.]', '_', name)

    # ---------- 读取输入文件 ----------
    _log("开始读取输入文件...")
    
    # 读取 susie_summary.json
    with open(susie_summary_json_path, "r", encoding="utf-8") as f:
        susie_summary = json.load(f)
    
    # 读取 ld_matrices_summary.json
    with open(ld_matrices_summary_json_path, "r", encoding="utf-8") as f:
        ld_matrices_summary = json.load(f)
    
    # 确定输出目录
    if output_dir is None:
        output_dir = os.path.dirname(os.path.abspath(ld_matrices_summary_json_path))
    os.makedirs(output_dir, exist_ok=True)
    
    _log(f"输出目录: {output_dir}")
    
    # ---------- 处理每个 lead variant ----------
    per_lead_outputs = ld_matrices_summary.get("per_lead_outputs", {})
    integrated_results = {
        "created_at": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        "inputs": {
            "susie_summary_json": os.path.abspath(susie_summary_json_path),
            "ld_matrices_summary_json": os.path.abspath(ld_matrices_summary_json_path),
            "output_dir": os.path.abspath(output_dir),
            "output_prefix": output_prefix,
        },
        "per_lead_integrated": {},
        "credible_sets_summary": {},
    }
    
    all_credible_sets = []
    
    for lead_variant, lead_data in per_lead_outputs.items():
        _log(f"处理 lead variant: {lead_variant}")
        
        # 检查是否有对应的 SuSiE 结果
        if lead_variant not in susie_summary:
            _log(f"警告: 在 susie_summary 中未找到 {lead_variant}，跳过")
            continue
        
        susie_info = susie_summary[lead_variant]
        if not susie_info.get("json_exists", False):
            _log(f"警告: {lead_variant} 的 SuSiE JSON 文件不存在，跳过")
            continue
        
        # 读取 sum_stat_tsv
        sum_stat_tsv_path = lead_data.get("sum_stat_tsv")
        if not sum_stat_tsv_path or not os.path.exists(sum_stat_tsv_path):
            _log(f"警告: {lead_variant} 的 sum_stat_tsv 文件不存在: {sum_stat_tsv_path}")
            continue
        
        try:
            sum_stat_df = pd.read_csv(sum_stat_tsv_path, sep="\t", dtype={"SNPID": str})
        except Exception as e:
            _log(f"错误: 读取 {sum_stat_tsv_path} 失败: {e}")
            continue
        
        # 读取 SuSiE JSON 结果
        susie_json_path = susie_info.get("json_file")
        if not susie_json_path or not os.path.exists(susie_json_path):
            _log(f"警告: {lead_variant} 的 SuSiE JSON 文件不存在: {susie_json_path}")
            continue
        
        try:
            with open(susie_json_path, "r", encoding="utf-8") as f:
                susie_results = json.load(f)
        except Exception as e:
            _log(f"错误: 读取 {susie_json_path} 失败: {e}")
            continue
        
        # ---------- 提取 SuSiE 结果中的 PIP 和 LD 信息 ----------
        # SuSiE 结果的 pip 是一个列表，每个元素包含 variant_id, pip, ld_r_with_lead, ld_r2_with_lead
        pip_list = susie_results.get("pip", [])
        
        pip_data = {}
        ld_r_data = {}
        ld_r2_data = {}
        
        if isinstance(pip_list, list):
            for item in pip_list:
                if isinstance(item, dict) and 'variant_id' in item:
                    variant_id = item['variant_id']
                    pip_data[variant_id] = item.get('pip', 0.0)
                    ld_r_data[variant_id] = item.get('ld_r_with_lead', 0.0)
                    ld_r2_data[variant_id] = item.get('ld_r2_with_lead', 0.0)
        else:
            _log(f"警告: {lead_variant} 的 pip 数据不是列表格式: {type(pip_list)}")
        
        _log(f"提取到 {len(pip_data)} 个变体的 PIP 信息")
        
        # 合并 PIP 和 LD 信息到 sum_stat_df
        sum_stat_df["PIP"] = sum_stat_df["SNPID"].map(pip_data).fillna(0.0)
        sum_stat_df["LD_R_WITH_LEAD"] = sum_stat_df["SNPID"].map(ld_r_data).fillna(0.0)
        sum_stat_df["LD_R2_WITH_LEAD"] = sum_stat_df["SNPID"].map(ld_r2_data).fillna(0.0)
        
        # 保存增强后的 sum_stat 文件
        sanitized_lead = _sanitize_filename(lead_variant)
        enhanced_sumstat_path = os.path.join(
            output_dir, f"{output_prefix}.{sanitized_lead}.enhanced_sumstat.tsv"
        )
        sum_stat_df.to_csv(enhanced_sumstat_path, sep="\t", index=False)
        
        # ---------- 提取可信集信息 ----------
        # SuSiE 结果的 credible_sets 是一个列表，每个元素包含一个可信集的信息
        credible_sets_list = susie_results.get("credible_sets", [])
        lead_credible_sets = []
        
        if isinstance(credible_sets_list, list):
            for cs_data in credible_sets_list:
                if not isinstance(cs_data, dict):
                    continue
                
                cs_index = cs_data.get("cs_index", "unknown")
                cs_name = f"CS{cs_index}"
                cs_coverage = cs_data.get("coverage", None)
                cs_size = cs_data.get("size", None)
                cs_variants = cs_data.get("variants", [])
                
                # 提取 purity 信息中的相关系数统计
                purity_info = cs_data.get("purity", {})
                cs_min_abs_corr = None
                cs_mean_abs_corr = None
                cs_median_abs_corr = None
                
                if isinstance(purity_info, dict):
                    # 提取 min_abs_corr
                    min_abs_corr_list = purity_info.get("min_abs_corr", [])
                    if isinstance(min_abs_corr_list, list) and len(min_abs_corr_list) > 0:
                        first_item = min_abs_corr_list[0]
                        if isinstance(first_item, dict):
                            cs_min_abs_corr = first_item.get("min.abs.corr", None)
                    
                    # 提取 mean_abs_corr
                    mean_abs_corr_list = purity_info.get("mean_abs_corr", [])
                    if isinstance(mean_abs_corr_list, list) and len(mean_abs_corr_list) > 0:
                        first_item = mean_abs_corr_list[0]
                        if isinstance(first_item, dict):
                            cs_mean_abs_corr = first_item.get("mean.abs.corr", None)
                    
                    # 提取 median_abs_corr
                    median_abs_corr_list = purity_info.get("median_abs_corr", [])
                    if isinstance(median_abs_corr_list, list) and len(median_abs_corr_list) > 0:
                        first_item = median_abs_corr_list[0]
                        if isinstance(first_item, dict):
                            cs_median_abs_corr = first_item.get("median.abs.corr", None)
                
                if isinstance(cs_variants, list):
                    for i, variant_item in enumerate(cs_variants):
                        # 从可信集变体中提取变体ID
                        if isinstance(variant_item, dict):
                            variant_id = variant_item.get('variant_id', str(variant_item))
                        else:
                            variant_id = str(variant_item)
                        
                        # 从 pip_data 中获取该变体的 PIP 值
                        pip_val = pip_data.get(variant_id, 0.0)
                        
                        cs_record = {
                            "locus_lead_variant": lead_variant,  # 更清晰的命名：该位点的主导变体
                            "credible_set": cs_name,
                            "variant_id": variant_id,
                            "is_lead_variant": variant_id == lead_variant,  # 标识是否为主导变体
                            "pip": float(pip_val),
                            "rank_in_cs": i + 1,
                            "cs_coverage": float(cs_coverage) if cs_coverage is not None else None,
                            "cs_min_abs_corr": float(cs_min_abs_corr) if cs_min_abs_corr is not None else None,
                            "cs_mean_abs_corr": float(cs_mean_abs_corr) if cs_mean_abs_corr is not None else None,
                            "cs_median_abs_corr": float(cs_median_abs_corr) if cs_median_abs_corr is not None else None,
                            "cs_size": int(cs_size) if cs_size is not None else len(cs_variants),
                        }
                        
                        # 从 sum_stat_df 中补充变体的统计信息
                        variant_stats = sum_stat_df[sum_stat_df["SNPID"] == variant_id]
                        if not variant_stats.empty:
                            variant_row = variant_stats.iloc[0]
                            cs_record.update({
                                "chr": int(val) if pd.notna(val := variant_row.get("CHR")) else None,  # type: ignore
                                "pos": int(val) if pd.notna(val := variant_row.get("POS")) else None,  # type: ignore
                                "ea": str(val) if pd.notna(val := variant_row.get("EA")) else None,  # type: ignore
                                "nea": str(val) if pd.notna(val := variant_row.get("NEA")) else None,  # type: ignore
                                "eaf": float(val) if pd.notna(val := variant_row.get("EAF")) else None,  # type: ignore
                                "beta": float(val) if pd.notna(val := variant_row.get("BETA")) else None,  # type: ignore
                                "se": float(val) if pd.notna(val := variant_row.get("SE")) else None,  # type: ignore
                                "z": float(val) if pd.notna(val := variant_row.get("Z")) else None,  # type: ignore
                                "p": float(val) if pd.notna(val := variant_row.get("P")) else None,  # type: ignore
                                "or": float(val) if pd.notna(val := variant_row.get("OR")) else None,  # type: ignore
                                "or_95l": float(val) if pd.notna(val := variant_row.get("OR_95L")) else None,  # type: ignore
                                "or_95u": float(val) if pd.notna(val := variant_row.get("OR_95U")) else None,  # type: ignore
                                "ld_r_with_lead": float(val) if pd.notna(val := variant_row.get("LD_R_WITH_LEAD")) else None,  # type: ignore
                                "ld_r2_with_lead": float(val) if pd.notna(val := variant_row.get("LD_R2_WITH_LEAD")) else None,  # type: ignore
                            })
                        
                        lead_credible_sets.append(cs_record)
                        all_credible_sets.append(cs_record)
        else:
            _log(f"警告: {lead_variant} 的 credible_sets 不是列表格式: {type(credible_sets_list)}")
        
        _log(f"提取到 {len(lead_credible_sets)} 个可信集变体")
        
        # 为每个 lead variant 单独保存可信集文件
        lead_credible_sets_path = None
        if lead_credible_sets:
            lead_credible_sets_df = pd.DataFrame(lead_credible_sets)
            lead_credible_sets_path = os.path.join(
                output_dir, f"{output_prefix}.{sanitized_lead}.credible_sets.tsv"
            )
            lead_credible_sets_df.to_csv(lead_credible_sets_path, sep="\t", index=False)
            _log(f"保存 {lead_variant} 可信集文件: {lead_credible_sets_path} ({len(lead_credible_sets_df)} 行)")

        # 记录结果
        integrated_results["per_lead_integrated"][lead_variant] = {
            "enhanced_sumstat_path": enhanced_sumstat_path,
            "credible_sets_path": lead_credible_sets_path,
            "original_sumstat_path": sum_stat_tsv_path,
            "susie_json_path": susie_json_path,
            "n_variants_total": len(sum_stat_df),
            "n_variants_with_pip": len(sum_stat_df[sum_stat_df["PIP"] > 0]),
            "n_credible_sets": len(set(cs["credible_set"] for cs in lead_credible_sets)),
            "n_credible_variants": len(lead_credible_sets),
        }
        
        integrated_results["credible_sets_summary"][lead_variant] = lead_credible_sets
        
        _log(f"完成 {lead_variant}: {len(sum_stat_df)} 个变体，{len(lead_credible_sets)} 个可信集变体")
    
    # ---------- 保存可信集汇总表（合并所有 lead 的可信集） ----------
    merged_credible_sets_tsv_path = None
    if all_credible_sets:
        credible_sets_df = pd.DataFrame(all_credible_sets)
        merged_credible_sets_tsv_path = os.path.join(output_dir, f"{output_prefix}.merged_credible_sets.tsv")
        credible_sets_df.to_csv(merged_credible_sets_tsv_path, sep="\t", index=False)
        _log(f"保存合并可信集汇总表: {merged_credible_sets_tsv_path} ({len(credible_sets_df)} 行)")
    else:
        _log("未找到任何可信集数据")
    
    # ---------- 保存整合结果的 JSON ----------
    integrated_json_path = os.path.join(output_dir, f"{output_prefix}.integrated_results.json")
    
    # 收集所有单独的可信集文件路径
    individual_credible_sets_files = {}
    for lead_variant in integrated_results["per_lead_integrated"]:
        cs_path = integrated_results["per_lead_integrated"][lead_variant].get("credible_sets_path")
        if cs_path:
            individual_credible_sets_files[lead_variant] = cs_path
    
    integrated_results["outputs"] = {
        "merged_credible_sets_tsv": merged_credible_sets_tsv_path,  # 合并所有 lead 的可信集文件
        "individual_credible_sets": individual_credible_sets_files,  # 每个 lead 的单独可信集文件
        "integrated_json": integrated_json_path,
    }
    
    with open(integrated_json_path, "w", encoding="utf-8") as f:
        json.dump(integrated_results, f, ensure_ascii=False, indent=2)
    
    _log(f"保存整合结果 JSON: {integrated_json_path}")
    _log(f"处理完成: {len(integrated_results['per_lead_integrated'])} 个 lead variants")
    
    return {
        "success": True,
        "integrated_json_path": integrated_json_path,
        "merged_credible_sets_tsv_path": merged_credible_sets_tsv_path,  # 合并所有 lead 的可信集文件
        "individual_credible_sets_files": individual_credible_sets_files,  # 每个 lead 的单独文件字典
        "n_leads_processed": len(integrated_results["per_lead_integrated"]),
        "n_credible_variants": len(all_credible_sets),
    }


def plot_manhattan_plots_from_integrated_results(
    integrated_json_path: str,
    out_pdf_path: Optional[str] = None,
    output_prefix: Optional[str] = None,
    highlight_credible_sets: bool = True,
    dpi: int = 150,
    figsize: tuple = (12, 8),
    logger: Optional[logging.Logger] = None,
) -> str:
    """
    基于 integrate_susie_results_with_sumstat 产生的整合结果 JSON，
    为每个 lead variant 生成曼哈顿图：上方显示 -log10(P)，下方显示 PIP。
    每个 lead variant 占用 PDF 的一页。

    参数：
        integrated_json_path (str): integrated_results.json 文件路径
        out_pdf_path (Optional[str]): 输出 PDF 路径；若为空则自动生成
        output_prefix (Optional[str]): 输出文件前缀
        highlight_credible_sets (bool): 是否高亮显示可信集变体
        dpi (int): 图像分辨率
        figsize (tuple): 图像尺寸 (width, height)
        logger (Optional[logging.Logger]): 日志记录器

    返回：
        str: 输出 PDF 的绝对路径
    """
    # ---------- 辅助函数 ----------
    def _log(msg: str):
        if logger:
            logger.info(msg)
        else:
            print(f"[plot_manhattan_plots] {msg}")

    def _parse_chr_pos(variant_id: str) -> tuple:
        """解析变体ID获取染色体和位置信息"""
        try:
            parts = str(variant_id).split(":")
            if len(parts) >= 2:
                chr_str = parts[0].replace("CHR", "").replace("chr", "")
                pos = int(parts[1])
                # 处理染色体编号
                if chr_str.upper() == "X":
                    chr_num = 23
                elif chr_str.upper() == "Y":
                    chr_num = 24
                elif chr_str.upper() in ("M", "MT"):
                    chr_num = 25
                else:
                    chr_num = int(chr_str)
                return chr_num, pos
        except:
            pass
        return 0, 0

    # ---------- 读取整合结果 ----------
    _log(f"读取整合结果: {integrated_json_path}")
    with open(integrated_json_path, "r", encoding="utf-8") as f:
        integrated_results = json.load(f)

    per_lead_integrated = integrated_results.get("per_lead_integrated", {})
    if not per_lead_integrated:
        raise ValueError("整合结果中没有找到 per_lead_integrated 数据")

    # 确定输出路径
    if output_prefix is not None and (out_pdf_path is None or str(out_pdf_path).strip() == ""):
        base_dir = os.path.dirname(os.path.abspath(integrated_json_path))
        out_pdf_path = os.path.join(base_dir, f"{output_prefix}.manhattan_plots.pdf")
    elif out_pdf_path is None or str(out_pdf_path).strip() == "":
        base_dir = os.path.dirname(os.path.abspath(integrated_json_path))
        base_name = os.path.splitext(os.path.basename(integrated_json_path))[0]
        out_pdf_path = os.path.join(base_dir, f"{base_name}.manhattan_plots.pdf")
    
    out_pdf_path = os.path.abspath(out_pdf_path)
    _log(f"输出PDF路径: {out_pdf_path}")

    # ---------- 为每个 lead variant 生成图表 ----------
    from matplotlib.backends.backend_pdf import PdfPages
    import matplotlib.pyplot as plt
    
    # 设置学术发表标准的图表样式
    # 优先使用 Arial 或 Helvetica，如果不可用则回退到其他无衬线字体
    plt.rcParams.update({
        'font.family': 'sans-serif',   # 使用无衬线字体
        'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans', 'Bitstream Vera Sans', 'sans-serif'],
        'font.size': 12,               # 增大基础字体
        'axes.titlesize': 14,          # 增大标题字体
        'axes.labelsize': 14,          # 增大轴标签字体
        'xtick.labelsize': 10,         # 增大X轴刻度字体
        'ytick.labelsize': 10,         # 增大Y轴刻度字体
        'legend.fontsize': 10,         # 增大图例字体
        'figure.titlesize': 16,        # 增大图像标题字体
        'axes.linewidth': 1.5,         # 加粗轴线
        'grid.linewidth': 1.0,         # 网格线宽度
        'lines.linewidth': 2.0,        # 线条宽度
        'patch.linewidth': 1.5,        # 补丁线宽度
        'axes.edgecolor': 'black',     # 轴边框颜色
        'axes.facecolor': 'white',     # 图表背景色
        'figure.facecolor': 'white',   # 图像背景色
        'grid.alpha': 0.3,             # 网格透明度
        'axes.axisbelow': True,        # 网格在数据下方
        'figure.autolayout': False,    # 禁用自动布局，手动控制
        'savefig.bbox': 'tight',       # 保存时紧凑布局
        # 强制设置字体颜色为黑色，防止暗色主题导致文字变白（不可见）
        'text.color': 'black',
        'axes.labelcolor': 'black',
        'xtick.color': 'black',
        'ytick.color': 'black',
    })
    
    with PdfPages(out_pdf_path) as pdf:
        for lead_variant, lead_data in per_lead_integrated.items():
            _log(f"处理 lead variant: {lead_variant}")
            
            # 读取增强的汇总统计数据
            enhanced_sumstat_path = lead_data.get("enhanced_sumstat_path")
            if not enhanced_sumstat_path or not os.path.exists(enhanced_sumstat_path):
                _log(f"警告: 跳过 {lead_variant}，增强汇总统计文件不存在")
                continue
            
            try:
                df = pd.read_csv(enhanced_sumstat_path, sep="\t", dtype={"SNPID": str})
            except Exception as e:
                _log(f"错误: 读取 {enhanced_sumstat_path} 失败: {e}")
                continue
            
            # 检查必需的列
            required_cols = ["SNPID", "P", "PIP"]
            missing_cols = [col for col in required_cols if col not in df.columns]
            if missing_cols:
                _log(f"警告: 跳过 {lead_variant}，缺少必需列: {missing_cols}")
                continue
            
            # 过滤有效数据
            df = df.dropna(subset=["P", "PIP"])
            df = df[df["P"] > 0]  # 确保P值为正数以计算-log10
            
            if df.empty:
                _log(f"警告: 跳过 {lead_variant}，没有有效数据")
                continue
            
            # 计算-log10(P)
            df["-log10P"] = -np.log10(df["P"])
            
            # 解析染色体和位置信息
            chr_pos_info = df["SNPID"].apply(_parse_chr_pos)
            df["CHR_NUM"] = [x[0] for x in chr_pos_info]
            df["POS_PARSED"] = [x[1] for x in chr_pos_info]
            
            # 按染色体和位置排序
            df = df.sort_values(["CHR_NUM", "POS_PARSED"])
            
            # 获取可信集信息
            credible_sets_df = None
            if highlight_credible_sets:
                credible_sets_path = lead_data.get("credible_sets_path")
                if credible_sets_path and os.path.exists(credible_sets_path):
                    try:
                        credible_sets_df = pd.read_csv(credible_sets_path, sep="\t", dtype={"variant_id": str})
                    except Exception as e:
                        _log(f"警告: 读取可信集文件失败: {e}")
            
            # 创建图表 - 学术发表标准布局
            # 调整为纵向布局 (width=10, height=12) 以获得更好的可读性
            fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 12), dpi=dpi, 
                                         gridspec_kw={'height_ratios': [1, 1], 'hspace': 0.3})
            
            # 设置图表背景和边框
            for ax in [ax1, ax2]:
                ax.spines['top'].set_visible(False)
                ax.spines['right'].set_visible(False)
                ax.spines['left'].set_linewidth(1.2)
                ax.spines['bottom'].set_linewidth(1.2)
                ax.tick_params(direction='out', length=4, width=1.2)
            
            # 准备基于 LD_R2_WITH_LEAD 的颜色映射
            import matplotlib.ticker as ticker
            
            # 使用物理位置作为 X 轴
            x_values = df["POS_PARSED"]
            
            # 确定染色体名称用于 X 轴标签
            unique_chrs = df["CHR_NUM"].unique()
            valid_chrs = [c for c in unique_chrs if c > 0]
            if len(valid_chrs) == 1:
                chr_num = int(valid_chrs[0])
                if chr_num == 23:
                    chr_name = "X"
                elif chr_num == 24:
                    chr_name = "Y"
                elif chr_num == 25:
                    chr_name = "MT"
                else:
                    chr_name = str(chr_num)
                xlabel_text = f'Chr{chr_name}'
            else:
                xlabel_text = 'Chr'
            
            # 检查是否有 LD_R2_WITH_LEAD 列，如果没有则使用默认值 0
            if "LD_R2_WITH_LEAD" in df.columns:
                ld_r2_values = df["LD_R2_WITH_LEAD"].fillna(0.0)
            else:
                _log(f"警告: {lead_variant} 缺少 LD_R2_WITH_LEAD 列，使用默认值 0")
                ld_r2_values = pd.Series([0.0] * len(df), index=df.index)
            
            # 使用 viridis 颜色映射，范围 [0, 1]
            # viridis: 深紫色(低LD) -> 蓝色 -> 绿色 -> 黄绿色(高LD)
            # 优点：感知均匀、色盲友好、全范围高对比度、科学标准
            cmap = plt.cm.viridis
            vmin, vmax = 0.0, 1.0
            
            # 上方子图：-log10(P) 曼哈顿图 - 学术标准样式
            scatter1 = ax1.scatter(x_values, df["-log10P"], c=ld_r2_values, s=25, alpha=0.8, 
                                  cmap=cmap, vmin=vmin, vmax=vmax, edgecolors='none', rasterized=True)
            
            # 标记lead variant - 更显眼的样式
            lead_mask = df["SNPID"] == lead_variant
            if lead_mask.any():
                lead_idx = df[lead_mask].index[0]
                lead_pos = df.loc[lead_idx, "POS_PARSED"]
                ax1.scatter(lead_pos, df.loc[lead_idx, "-log10P"], 
                           c='darkred', s=120, marker='D', edgecolors='white', linewidth=2, 
                           label='Lead Variant', zorder=10, alpha=0.9)
            
            # 高亮可信集变体 - 支持多个可信集，用不同颜色和形状区分
            if credible_sets_df is not None and not credible_sets_df.empty:
                # 获取所有可信集的名称
                unique_cs = credible_sets_df["credible_set"].unique()
                # 选择与 viridis 配色协调的边框颜色 - 高对比度，易于在紫-蓝-绿-黄背景上识别
                cs_colors = ['red', 'orange', 'white', 'black', 'magenta', 'cyan']  # 高对比度边框，与viridis形成清晰区分
                cs_markers = ['o', 's', '^', 'v', '<', '>']  # 多个可信集的形状
                
                for i, cs_name in enumerate(unique_cs):
                    cs_variants_in_set = set(credible_sets_df[credible_sets_df["credible_set"] == cs_name]["variant_id"].unique())
                    cs_mask_in_set = df["SNPID"].isin(cs_variants_in_set)
                    
                    if cs_mask_in_set.any():
                        cs_indices_in_set = df.loc[df[cs_mask_in_set].index, "POS_PARSED"]
                        # 获取这些变体的 LD R² 值用于颜色映射
                        cs_ld_r2_values = df.loc[df[cs_mask_in_set].index, "LD_R2_WITH_LEAD"].fillna(0.0) if "LD_R2_WITH_LEAD" in df.columns else [0.0] * len(cs_indices_in_set)
                        
                        color = cs_colors[i % len(cs_colors)]
                        marker = cs_markers[i % len(cs_markers)]
                        
                        # 使用 LD R² 值映射颜色，但用固定颜色的边框来区分可信集 - 学术样式
                        ax1.scatter(cs_indices_in_set, df.loc[df[cs_mask_in_set].index, "-log10P"],
                                   c=cs_ld_r2_values, s=60, marker=marker, cmap=cmap, vmin=vmin, vmax=vmax,
                                   edgecolors=color, linewidth=2.0, alpha=0.9, zorder=5, rasterized=True)
            
            # 添加显著性阈值线 - 学术标准样式
            ax1.axhline(y=-np.log10(5e-8), color='#d62728', linestyle='--', linewidth=2, alpha=0.8, 
                       label='Genome-wide Sig (5×10⁻⁸)', zorder=1)
            ax1.axhline(y=-np.log10(1e-5), color='#1f77b4', linestyle='--', linewidth=2, alpha=0.8, 
                       label='Suggestive (1×10⁻⁵)', zorder=1)
            
            # 为图例添加清晰的可信集标识（仅显示边框颜色和形状）
            if credible_sets_df is not None and not credible_sets_df.empty:
                unique_cs = credible_sets_df["credible_set"].unique()
                cs_colors = ['red', 'orange', 'white', 'black', 'magenta', 'cyan']
                cs_markers = ['o', 's', '^', 'v', '<', '>']
                for i, cs_name in enumerate(unique_cs):
                    color = cs_colors[i % len(cs_colors)]
                    marker = cs_markers[i % len(cs_markers)]
                    # 添加专门的图例项，使用更好的学术样式
                    ax1.scatter([], [], c='lightgray', s=60, marker=marker, 
                               edgecolors=color, linewidth=2.0, alpha=0.9,
                               label=f'Credible Set {cs_name[2:]}')
            
            ax1.set_xlabel(xlabel_text, fontweight='bold')
            ax1.set_ylabel('-log$_{10}$(P)', fontweight='bold')
            ax1.set_title(f'{lead_variant} — Association P-values', fontweight='bold', pad=15)
            
            # 格式化 X 轴标签，添加千分位分隔符
            ax1.xaxis.set_major_formatter(ticker.StrMethodFormatter('{x:,.0f}'))
            
            # 优化图例样式
            legend1 = ax1.legend(loc='upper right', fontsize=10, frameon=True, fancybox=True, 
                               shadow=True, framealpha=0.9, edgecolor='gray',
                               ncol=2 if credible_sets_df is not None and len(credible_sets_df["credible_set"].unique()) > 3 else 1)
            legend1.get_frame().set_linewidth(1.2)
            ax1.grid(True, alpha=0.4, linestyle='-', linewidth=0.8)
            
            # 下方子图：PIP 曼哈顿图 - 学术标准样式
            scatter2 = ax2.scatter(x_values, df["PIP"], c=ld_r2_values, s=25, alpha=0.8, 
                                  cmap=cmap, vmin=vmin, vmax=vmax, edgecolors='none', rasterized=True)
            
            # 标记lead variant - 更显眼的样式
            if lead_mask.any():
                ax2.scatter(lead_pos, df.loc[lead_idx, "PIP"], 
                           c='darkred', s=120, marker='D', edgecolors='white', linewidth=2, 
                           label='Lead Variant', zorder=10, alpha=0.9)
            
            # 高亮可信集变体 - 支持多个可信集，用不同颜色和形状区分
            if credible_sets_df is not None and not credible_sets_df.empty:
                # 获取所有可信集的名称
                unique_cs = credible_sets_df["credible_set"].unique()
                # 选择与 viridis 配色协调的边框颜色 - 高对比度，易于在紫-蓝-绿-黄背景上识别
                cs_colors = ['red', 'orange', 'white', 'black', 'magenta', 'cyan']  # 高对比度边框，与viridis形成清晰区分
                cs_markers = ['o', 's', '^', 'v', '<', '>']  # 多个可信集的形状
                
                for i, cs_name in enumerate(unique_cs):
                    cs_variants_in_set = set(credible_sets_df[credible_sets_df["credible_set"] == cs_name]["variant_id"].unique())
                    cs_mask_in_set = df["SNPID"].isin(cs_variants_in_set)
                    
                    if cs_mask_in_set.any():
                        cs_indices_in_set = df.loc[df[cs_mask_in_set].index, "POS_PARSED"]
                        # 获取这些变体的 LD R² 值用于颜色映射
                        cs_ld_r2_values = df.loc[df[cs_mask_in_set].index, "LD_R2_WITH_LEAD"].fillna(0.0) if "LD_R2_WITH_LEAD" in df.columns else [0.0] * len(cs_indices_in_set)
                        
                        color = cs_colors[i % len(cs_colors)]
                        marker = cs_markers[i % len(cs_markers)]
                        
                        # 使用 LD R² 值映射颜色，但用固定颜色的边框来区分可信集 - 学术样式
                        ax2.scatter(cs_indices_in_set, df.loc[df[cs_mask_in_set].index, "PIP"],
                                   c=cs_ld_r2_values, s=60, marker=marker, cmap=cmap, vmin=vmin, vmax=vmax,
                                   edgecolors=color, linewidth=2.0, alpha=0.9, zorder=5, rasterized=True)
            
            # 添加PIP阈值线 - 学术标准样式
            ax2.axhline(y=0.1, color='#2ca02c', linestyle='--', linewidth=2, alpha=0.8, 
                       label='PIP = 0.1', zorder=1)
            ax2.axhline(y=0.5, color='#ff7f0e', linestyle='--', linewidth=2, alpha=0.8, 
                       label='PIP = 0.5', zorder=1)
            ax2.axhline(y=0.9, color='#d62728', linestyle='--', linewidth=2, alpha=0.8, 
                       label='PIP = 0.9', zorder=1)
            
            # 为图例添加清晰的可信集标识（仅显示边框颜色和形状）
            if credible_sets_df is not None and not credible_sets_df.empty:
                unique_cs = credible_sets_df["credible_set"].unique()
                cs_colors = ['red', 'orange', 'white', 'black', 'magenta', 'cyan']
                cs_markers = ['o', 's', '^', 'v', '<', '>']
                for i, cs_name in enumerate(unique_cs):
                    color = cs_colors[i % len(cs_colors)]
                    marker = cs_markers[i % len(cs_markers)]
                    # 添加专门的图例项，使用更好的学术样式
                    ax2.scatter([], [], c='lightgray', s=60, marker=marker, 
                               edgecolors=color, linewidth=2.0, alpha=0.9,
                               label=f'Credible Set {cs_name[2:]}')
            
            ax2.set_xlabel(xlabel_text, fontweight='bold')
            ax2.set_ylabel('Posterior Inclusion Probability (PIP)', fontweight='bold')
            ax2.set_title(f'{lead_variant} — SuSiE Posterior Inclusion Probabilities', fontweight='bold', pad=15)
            
            # 格式化 X 轴标签，添加千分位分隔符
            ax2.xaxis.set_major_formatter(ticker.StrMethodFormatter('{x:,.0f}'))
            
            # 优化图例样式
            legend2 = ax2.legend(loc='upper right', fontsize=10, frameon=True, fancybox=True, 
                               shadow=True, framealpha=0.9, edgecolor='gray',
                               ncol=2 if credible_sets_df is not None and len(credible_sets_df["credible_set"].unique()) > 3 else 1)
            legend2.get_frame().set_linewidth(1.2)
            ax2.grid(True, alpha=0.4, linestyle='-', linewidth=0.8)
            ax2.set_ylim(-0.05, 1.05)  # 稍微扩展Y轴范围以获得更好的视觉效果
            
            # 添加统一的颜色条（基于 LD R² with Lead） - 学术样式
            # 在右侧添加颜色条，跨越两个子图的高度
            cbar = fig.colorbar(scatter2, ax=[ax1, ax2], fraction=0.015, pad=0.02, aspect=25, shrink=0.8)
            cbar.set_label('LD R² with Lead Variant', rotation=270, labelpad=20, fontweight='bold', fontsize=12)
            cbar.ax.tick_params(labelsize=10, width=1.2)
            cbar.outline.set_linewidth(1.2)
            
            # 添加整体标题，包含可信集信息 - 学术样式
            n_credible_sets = len(credible_sets_df["credible_set"].unique()) if credible_sets_df is not None and not credible_sets_df.empty else 0
            n_credible_variants = len(credible_sets_df) if credible_sets_df is not None and not credible_sets_df.empty else 0
            
            fig.suptitle(f'Fine-mapping Results for {lead_variant}\n'
                        f'Total variants: {len(df):,} | '
                        f'Max -log$_{{10}}$(P): {df["-log10P"].max():.2f} | '
                        f'Max PIP: {df["PIP"].max():.3f}\n'
                        f'Credible Sets: {n_credible_sets} | CS variants: {n_credible_variants}', 
                        fontsize=16, fontweight='bold', y=0.98)
            
            # 优化布局和保存 - 学术发表标准
            plt.tight_layout(rect=[0, 0, 0.98, 0.95])  # 为suptitle和colorbar留出空间
            pdf.savefig(fig, bbox_inches='tight', dpi=dpi, facecolor='white', edgecolor='none')
            plt.close(fig)
            
            _log(f"完成 {lead_variant}: {len(df)} 个变体")
    
    # 恢复默认matplotlib设置
    plt.rcdefaults()
    
    _log(f"曼哈顿图已保存到: {out_pdf_path}")
    return out_pdf_path


