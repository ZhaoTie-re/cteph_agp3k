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
    ref_seq_path: str = "/LARGE0/gr10478/b37974/Pulmonary_Hypertension/nagasaki_pipeline/data/hs38DH.fa",
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
        ref_seq_path (str): 参考基因组 fasta 路径（用于 harmonize）。
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
            _call_with_verbose(locus.harmonize, basic_check=False, ref_seq=ref_seq_path) # type: ignore

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
            r_df.index = extract_ids
            r_df.columns = extract_ids
        else:
            _log(f"警告：R 阶数({r_df.shape[0]})与 extract 数量({len(extract_ids)})不一致，尝试基于 .bim 对齐。")
            inferred = _infer_labels_from_bim(bed_prefix, extract_ids, r_df.shape[0])
            if inferred is not None:
                r_df.index = inferred
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
                r2_df.index = extract_ids
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
    final_json = os.path.join(out_dir, "ld_matrices_by_lead.summary.json")
    with open(final_json, "w", encoding="utf-8") as jf:
        json.dump(out_index, jf, ensure_ascii=False, indent=2)
    _log(f"完成：输出索引 JSON -> {final_json}")
    return final_json

