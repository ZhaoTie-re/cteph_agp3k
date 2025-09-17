"""
蛋白质关联分析工具集

本模块提供用于CTEPH研究中基因型-蛋白质表达关联分析的工具函数，包括：

主要功能模块：
1. PLINK基因型数据处理
   - export_plink_geno_matrix: 从PLINK二进制文件导出基因型矩阵
   - summarize_gt_counts: 统计基因型计数信息

2. 变体-基因-蛋白质数据整合
   - assemble_variant_gene_protein_table: 整合变体、基因和蛋白质数据

3. 蛋白质表达可视化分析
   - plot_protein_boxplot_per_variant: 生成按变体分组的蛋白质表达箱线图

4. 辅助工具函数
   - 数据清洗和标准化函数
   - 统计分析函数
   - 图表渲染函数

依赖软件：
- PLINK2: 用于基因型数据处理
- bcftools: 用于VCF文件操作
- pandas, numpy, matplotlib, seaborn: 数据分析和可视化

作者：ZHAO TIE
最后更新：2025年9月
"""

import os
import sys
import shutil
import subprocess
import tempfile
import time
import shlex
from typing import Tuple, Optional, List, Dict, Any

import pandas as pd

# PLINK2和bcftools的默认路径配置
DEFAULT_PLINK2 = "/home/b/b37974/plink2"
DEFAULT_BCFTOOLS = "/home/b/b37974/bcftools/bcftools"

# ---- 外部工具执行的共享辅助函数 ----

def _log_info(msg: str):
    """记录信息级别的日志消息。优先使用全局log_info函数，否则输出到stderr。"""
    if "log_info" in globals() and callable(globals()["log_info"]):
        globals()["log_info"](msg)
    else:
        import sys
        print(f"[INFO] {msg}", file=sys.stderr, flush=True)

def _log_warn(msg: str):
    """记录警告级别的日志消息。优先使用全局log_warn函数，否则输出到stderr。"""
    if "log_warn" in globals() and callable(globals()["log_warn"]):
        globals()["log_warn"](msg)
    else:
        import sys
        print(f"[WARN] {msg}", file=sys.stderr, flush=True)

def _run(cmd):
    """
    执行外部命令的通用函数。
    
    优先使用全局_run_cmd函数，如果不可用则使用subprocess作为备用方案。
    如果命令执行失败，会抛出RuntimeError异常并包含错误信息。
    
    参数：
        cmd: 要执行的命令，可以是字符串或列表格式
        
    返回：
        subprocess.CompletedProcess对象（当使用备用方案时）
    """
    if "_run_cmd" in globals() and callable(globals()["_run_cmd"]):
        return globals()["_run_cmd"](cmd)
    if isinstance(cmd, (list, tuple)):
        cmd_disp = " ".join(cmd)
    else:
        cmd_disp = cmd
    _log_info(f"执行命令（fallback）: {cmd_disp}")
    import subprocess
    ret = subprocess.run(cmd, shell=isinstance(cmd, str), stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if ret.returncode != 0:
        tail_err = "\n".join((ret.stderr or "").strip().splitlines()[-20:])
        if tail_err:
            _log_warn(tail_err)
        raise RuntimeError(f"命令执行失败（退出码 {ret.returncode}）：{cmd_disp}")
    return ret

def _ensure(path: str):
    """
    确保指定目录存在。
    
    优先使用全局_ensure_dir函数，如果不可用则使用os.makedirs作为备用方案。
    
    参数：
        path: 要创建的目录路径
    """
    if "_ensure_dir" in globals() and callable(globals()["_ensure_dir"]):
        return globals()["_ensure_dir"](path)
    import os
    os.makedirs(path, exist_ok=True)


def export_plink_geno_matrix(
    plink_prefix,
    output_tsv,
    variant_ids=None,
    sample_ids=None,
    work_dir=None,
    plink_threads=8,
    keep_temp=False,
):
    """
    从 PLINK 二进制数据（.bed/.bim/.fam）导出“按 ALT 计数”的基因型矩阵（TSV），
    行为变体（ID=CHR:POS:REF:ALT），列为样本 IID，元素为 {0,1,2,'.'}（'.' 表示缺失）。

    参数
    ----
    plink_prefix : str
        PLINK 前缀（.bed/.bim/.fam）
    output_tsv : str
        输出 TSV 路径（例如 /path/to/chr22.post.gt.tsv）
    variant_ids : list[str] | None
        需要导出的变体 ID 列表（格式 CHR:POS:REF:ALT）。若为 None 或空，则导出全部变体。
        支持输入 '1'/'chr1' 前缀差异，会自动与 .bim 对齐。
    sample_ids : list[str] | None
        需要导出的样本 IID 列表。若为 None 或空，则导出全部样本。
    work_dir : str | None
        中间文件目录（默认：与 output_tsv 同目录）；中间文件写入 <work_dir>/tmp_export/
    plink_threads : int
        plink2 线程数
    keep_temp : bool
        是否保留中间文件目录

    返回
    ----
    str
        output_tsv 路径（便于链式调用）

    说明
    ----
    - 使用 plink2 `--export A-transpose` + `--export-allele`（从 .bim 的 ID 解析出 ALT）确保 0/1/2 为 ALT 计数；
    - 仅输出 `ID + 所选样本(IID)` 两类列，缺失为 '.'；
    - 自动处理 `.traw` 表头中的 FID_IID → IID 的映射（含启发式兜底）。
    """
    # ---- 兼容外部工具函数/常量（若未定义则提供兜底） ----
    plink_bin = globals().get("DEFAULT_PLINK2", DEFAULT_PLINK2)

    # warnings and summary
    warnings_log: List[str] = []
    summary_json = output_tsv + ".summary.json"
    summary: Dict[str, Any] = {
        "plink_prefix": plink_prefix,
        "output_tsv": output_tsv,
        "work_dir": work_dir if work_dir else os.path.dirname(os.path.abspath(output_tsv)),
        "threads": {"plink": plink_threads},
        "totals": {"variants_bim": 0, "samples_fam": 0},
        "request": {"variants": None, "samples": None},
        "hit": {"variants": None, "samples": None},
        "miss": {"variants": None, "samples": None},
        "pct_hit": {"variants": None, "samples": None},
        "warnings": warnings_log,
        "result": {"gt_tsv_path": None, "variants": None, "samples": None}
    }

    # defaults for summary
    total_samples_fam = 0
    total_variants_bim = 0
    hit_variants = None
    hit_samples = None

    # thresholds for small requests
    SMALL_N_VAR = 20
    SMALL_N_SAMPLE = 20

    # ---------- 路径检查 ----------
    fam_path = plink_prefix + ".fam"
    bim_path = plink_prefix + ".bim"
    if not os.path.exists(fam_path):
        raise FileNotFoundError(f"FAM not found: {fam_path}")
    if not os.path.exists(bim_path):
        raise FileNotFoundError(f"BIM not found: {bim_path}")

    # ---------- 工作/缓存目录 ----------
    out_dir = os.path.dirname(os.path.abspath(output_tsv)) or os.getcwd()
    base_tmp = work_dir if work_dir else out_dir
    _ensure(base_tmp)
    cache_dir = os.path.join(base_tmp, "tmp_export")
    _ensure(cache_dir)

    _log_info(f"输出：{output_tsv}")
    _log_info(f"工作目录：{base_tmp}；缓存目录：{cache_dir}")

    # ---------- 读取 FAM，建立 FID_IID → IID 映射 ----------
    fid_iid_to_iid = {}
    iid_to_fid = {}
    iid_set = set()
    with open(fam_path, "r") as f:
        for line in f:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split()
            if len(parts) < 2:
                continue
            fid, iid = parts[0], parts[1]
            fid_iid_to_iid[f"{fid}_{iid}"] = iid
            iid_to_fid[iid] = fid
            iid_set.add(iid)
    _log_info(f"读取 FAM：样本总数={len(iid_set)}")
    total_samples_fam = len(iid_set)
    summary["totals"]["samples_fam"] = total_samples_fam

    # ---- 若用户仅请求少量位点/样本，按需构建索引，避免全量加载占内存 ----
    def _strip_chr(c: str) -> str:
        return c[3:] if c.lower().startswith("chr") else c

    requested_pos_keys = None  # set of (chrom_wo_chr, pos_int)
    if variant_ids:
        requested_pos_keys = set()
        for _vid in variant_ids:
            _parts = _vid.split(":")
            if len(_parts) >= 2:
                _chrom = _strip_chr(_parts[0])
                try:
                    _pos = int(_parts[1])
                except ValueError:
                    continue
                requested_pos_keys.add((_chrom, _pos))

    # ---------- 读取 BIM（仅统计总数 + 按需为请求位点建立索引） ----------
    from collections import defaultdict
    bim_has_chr_prefix = None
    total_variants_bim = 0
    bim_pos_to_ids = defaultdict(list)  # 仅当位点在请求集合中时才填充
    with open(bim_path, "r") as f:
        for i, line in enumerate(f):
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            bid = parts[1]
            # 仅用首条有效记录判断 chr 风格
            if bim_has_chr_prefix is None:
                chrom_part = bid.split(":")[0]
                bim_has_chr_prefix = chrom_part.lower().startswith("chr")
            total_variants_bim += 1
            # 仅当用户给了变体列表时，才按需建立 CHR:POS -> BIM IDs 索引
            if requested_pos_keys:
                _p = bid.split(":")
                if len(_p) >= 2:
                    key_chrom = _strip_chr(_p[0])
                    try:
                        key_pos = int(_p[1])
                    except ValueError:
                        continue
                    key = (key_chrom, key_pos)
                    if key in requested_pos_keys:
                        bim_pos_to_ids[key].append(bid)

    if bim_has_chr_prefix is None:
        bim_has_chr_prefix = True
    _log_info(f"BIM ID 风格：{'chr*' if bim_has_chr_prefix else 'numeric'}")
    _log_info(f"BIM 变体总数={total_variants_bim}")
    summary["totals"]["variants_bim"] = total_variants_bim

    # ---------- 生成 export-allele 映射（ID \t ALT） ----------
    allele_map_path = os.path.join(cache_dir, "export_allele.from_bim_id_ALT.txt")
    _run(
        "awk -F'\\t' '{split($2,a,\":\"); if (length(a)>=4) print $2\"\\t\"a[4]; else {print \"[ERR] Bad BIM ID: \"$2 > \"/dev/stderr\"; exit 1}}' "
        + shlex.quote(bim_path) + " > " + shlex.quote(allele_map_path)
    )

    # ---------- 可选：--extract 变体列表（以 BIM 第二列为准：按 CHR:POS 匹配，REF/ALT 完全取自 BIM） ----------
    extract_path = None
    if variant_ids:
        requested_variants = len(variant_ids)
        selected_bim_ids = []     # 去重后的 BIM ID 列表（用于 --extract）
        seen = set()
        missing_variant_items = []  # 记录未命中的请求（按位点）
        used_order = []            # 输出行顺序（按请求顺序，将该请求命中的 BIM ID 依次追加）

        def parse_vid_to_pos(vid: str):
            parts = vid.split(":")
            if len(parts) < 2:
                return None
            chrom = parts[0]
            pos = parts[1]
            try:
                pos_int = int(pos)
            except ValueError:
                return None
            chrom_key = _strip_chr(chrom)
            return (chrom_key, pos_int)

        for vid in variant_ids:
            key = parse_vid_to_pos(vid)
            if key is None:
                missing_variant_items.append({"requested": vid, "reason": "bad_format_or_pos"})
                continue
            cands = bim_pos_to_ids.get(key, [])
            if not cands:
                missing_variant_items.append({"requested": vid, "reason": "pos_not_found"})
                continue
            # 将该请求命中的 BIM ID 依次追加，且全局去重
            for bid in cands:
                if bid not in seen:
                    selected_bim_ids.append(bid)
                    seen.add(bid)
            used_order.extend([bid for bid in cands])

        hit_variants = len(selected_bim_ids)
        # 未命中按“请求的位点中没有对应 POS”的数量统计
        miss_variants = sum(1 for x in missing_variant_items if x.get("reason") == "pos_not_found" or x.get("reason") == "bad_format_or_pos")
        pct_hit_var = ( (requested_variants - miss_variants) / requested_variants * 100.0 ) if requested_variants > 0 else 0.0
        _log_info(f"[变体选择] 请求={requested_variants}；按位点命中请求={requested_variants - miss_variants}（{pct_hit_var:.2f}%）；BIM提取ID数={hit_variants}")
        if requested_variants <= SMALL_N_VAR and miss_variants > 0:
            sneak = ", ".join([m['requested'] for m in missing_variant_items[:10]])
            msg = f"小规模变体请求：未命中 {miss_variants}/{requested_variants}（示例：{sneak}）"
            _log_warn(msg); warnings_log.append(msg)
        if (requested_variants - miss_variants) == 0:
            raise RuntimeError("提供的 variant_ids 在 BIM 中均无法按 CHR:POS 匹配。")

        # 写入 --extract（以 BIM 原始 ID 写入；REF/ALT 由 BIM 决定）
        extract_path = os.path.join(cache_dir, "extract.ids")
        with open(extract_path, "w") as f:
            for v in selected_bim_ids:
                f.write(v + "\n")
        _log_info(f"--extract 使用 {len(selected_bim_ids)} 个 BIM ID（由 {requested_variants} 个请求位点映射而来）")
        summary["request"]["variants"] = requested_variants
        summary["hit"]["variants"] = requested_variants - miss_variants
        summary["miss"]["variants"] = miss_variants
        summary["pct_hit"]["variants"] = round(pct_hit_var, 4)
        summary["missing_detail_variants"] = missing_variant_items[:100]
        summary["used_variant_order"] = used_order  # 用于输出行按请求位点顺序展开
        summary["paths_extract_ids"] = extract_path
    else:
        _log_info(f"未提供 variant_ids：将导出全部变体（BIM总数={total_variants_bim}）")
        summary["request"]["variants"] = total_variants_bim
        summary["hit"]["variants"] = total_variants_bim
        summary["miss"]["variants"] = 0
        summary["pct_hit"]["variants"] = 100.0

    # ---------- 可选：--keep 样本列表（FID IID） ----------
    keep_path = None
    if sample_ids:
        keep_pairs = []
        miss_cnt = 0
        for iid in sample_ids:
            fid = iid_to_fid.get(iid)
            if fid is None:
                miss_cnt += 1
                continue
            keep_pairs.append((fid, iid))
        missing_sample_ids = [iid for iid in sample_ids if iid_to_fid.get(iid) is None]
        requested_samples = len(sample_ids)
        hit_samples = len(keep_pairs)
        pct_hit_s = (hit_samples / requested_samples * 100.0) if requested_samples > 0 else 0.0
        _log_info(f"[样本选择] 请求={requested_samples}；命中={hit_samples}（{pct_hit_s:.2f}%）；未命中={miss_cnt}；FAM总数={total_samples_fam}")
        if hit_samples == 0:
            msg = "请求的样本一个都没命中 FAM（请检查 IID 列表）。"
            _log_warn(msg); warnings_log.append(msg)
        else:
            if requested_samples <= SMALL_N_SAMPLE:
                if miss_cnt > 0:
                    preview = ", ".join(missing_sample_ids[:10])
                    msg = f"小规模样本请求：未命中 {miss_cnt}/{requested_samples}（示例：{preview}）"
                    _log_warn(msg); warnings_log.append(msg)
            elif pct_hit_s < 80.0:
                msg = f"请求样本命中率偏低（{pct_hit_s:.2f}%），可能存在 IID 不一致或空白。"
                _log_warn(msg); warnings_log.append(msg)
        if len(keep_pairs) == 0:
            raise RuntimeError("提供的 sample_ids 在 FAM 中均不存在。请检查 IID。")
        summary["request"]["samples"] = requested_samples
        summary["hit"]["samples"] = hit_samples
        summary["miss"]["samples"] = miss_cnt
        summary["pct_hit"]["samples"] = round(pct_hit_s, 4)
        summary["missing_detail_samples"] = missing_sample_ids[:100]
        keep_path = os.path.join(cache_dir, "keep.samples.txt")
        with open(keep_path, "w") as f:
            for fid, iid in keep_pairs:
                f.write(f"{fid}\t{iid}\n")
        _log_info(f"--keep 使用 {len(keep_pairs)} 个样本")
        summary["paths_keep_samples"] = keep_path
    else:
        _log_info(f"未提供 sample_ids：将导出全部样本（FAM总数={total_samples_fam})")
        summary["request"]["samples"] = total_samples_fam
        summary["hit"]["samples"] = total_samples_fam
        summary["miss"]["samples"] = 0
        summary["pct_hit"]["samples"] = 100.0

    # ---------- 调用 plink2 导出 .traw ----------
    traw_prefix = os.path.join(cache_dir, "export.genotype")
    cmd = [
        plink_bin,
        "--bfile", plink_prefix,
        "--export", "A-transpose",
        "--export-allele", allele_map_path,
        "--threads", str(plink_threads),
        "--out", traw_prefix
    ]
    if extract_path:
        cmd.extend(["--extract", extract_path])
    if keep_path:
        cmd.extend(["--keep", keep_path])
    # 导出计划（与请求规模对比，而非与全库对比）
    var_plan = "全部" if not variant_ids else f"{hit_variants}/{summary['request']['variants']}"
    samp_plan = "全部" if not sample_ids else f"{hit_samples}/{summary['request']['samples']}"
    _log_info(f"[导出计划] 变体：{var_plan}；样本：{samp_plan}；线程：plink={plink_threads}")
    _run(cmd)

    traw_file = traw_prefix + ".traw"
    if not os.path.exists(traw_file):
        raise RuntimeError(f"未生成预期的 .traw 文件：{traw_file}")

    # ---------- 清洗 .traw → 只留 ID + IID，值为 0/1/2/'.' ----------
    hdr = pd.read_csv(traw_file, sep="\t", nrows=0)
    meta_cols = ["CHR","SNP","(C)M","POS","COUNTED","ALT"]
    sample_cols = [c for c in hdr.columns if c not in meta_cols]
    usecols = ["SNP"] + sample_cols

    # 小规模请求通常行列很少：减少 dtype 推断开销
    # （保持现有逻辑以确保缺失写为'.'，无需大改）
    df = pd.read_csv(traw_file, sep="\t", usecols=usecols, dtype=str)
    df = df.rename(columns={"SNP":"ID"})

    # 将样本列名标准化为 IID
    iid_cols = []
    heuristic_hits = 0
    for c in sample_cols:
        if c in fid_iid_to_iid:
            iid_cols.append(fid_iid_to_iid[c])
        else:
            if c in iid_set:
                iid_cols.append(c)
            else:
                if "_" in c:
                    after_first = c.split("_", 1)[1]
                    if after_first in iid_set:
                        iid_cols.append(after_first)
                        heuristic_hits += 1
                        continue
                    last_seg = c.rsplit("_", 1)[-1]
                    iid_cols.append(last_seg)
                    heuristic_hits += 1
                else:
                    iid_cols.append(c)
    if heuristic_hits > 0:
        _log_warn(f"有 {heuristic_hits} 个样本列名通过启发式推断 IID，请确认 FID/IID 映射是否完整")

    df.columns = ["ID"] + iid_cols

    # 统一缺失为 '.'；转换为可空整数以便写出 '.' 作为缺失
    for c in iid_cols:
        s = pd.to_numeric(df[c], errors="coerce")
        try:
            df[c] = s.astype("Int8")
        except TypeError:
            df[c] = s.astype("Int64")

    # 若提供了 variant_ids，则按请求顺序（归一化到 BIM 风格后）重排行顺序
    if variant_ids:
        order_ids = summary.get("used_variant_order")
        if order_ids:
            df = df.set_index("ID").reindex(order_ids).reset_index()

    # 若用户指定了 sample_ids，则按用户次序重排列（仅保留存在的）
    if sample_ids:
        ordered_cols = ["ID"] + [iid for iid in sample_ids if iid in iid_cols]
        df = df[ordered_cols]

    # 写出 TSV
    os.makedirs(os.path.dirname(os.path.abspath(output_tsv)), exist_ok=True)
    df.to_csv(output_tsv, sep="\t", index=False, na_rep=".")
    _log_info(f"已写出基因型矩阵：{output_tsv}（行={df.shape[0]}，列={df.shape[1]-1} 样本）")
    _log_info(f"[导出完成] 实际导出：变体={df.shape[0]}；样本={df.shape[1]-1}")

    # 更新结果尺寸并写出 JSON 汇总
    summary["result"]["gt_tsv_path"] = os.path.abspath(output_tsv)
    summary["result"]["variants"] = int(df.shape[0])
    summary["result"]["samples"] = int(df.shape[1] - 1)
    try:
        import json
        with open(summary_json, "w", encoding="utf-8") as jf:
            json.dump(summary, jf, ensure_ascii=False, indent=2)
        _log_info(f"[SUMMARY] 已写出 JSON 汇总：{summary_json}")
    except Exception as e:
        msg = f"写出 JSON 汇总失败：{e}"
        _log_warn(msg); warnings_log.append(msg)

    # 清理
    if not keep_temp:
        try:
            if os.path.isdir(cache_dir):
                shutil.rmtree(cache_dir)
        except Exception as e:
            _log_warn(f"临时目录清理失败：{e}")

    return output_tsv


def summarize_gt_counts(
    gt_tsv_path: str,
    out_path: Optional[str] = None,
    chunksize: int = 5000,
    preview_n: int = 10,
    print_preview: bool = True,
) -> str:
    """
    统计基因型矩阵中每个变体的基因型分布。
    
    从基因型矩阵（TSV格式）中统计每个变体的GT=0/1/2/缺失的样本数量，
    生成按行对齐的汇总表格。支持大文件的分块处理以节省内存。
    
    参数
    ----
    gt_tsv_path : str
        输入的基因型矩阵文件路径（由export_plink_geno_matrix生成）。
        格式：第一列为变体ID，其余列为样本，值为'0'/'1'/'2'/'.'
        
    out_path : str, optional
        输出汇总文件路径。如果未指定，将在输入文件同目录下生成
        `<basename>.gt_summary.tsv`文件。
        
    chunksize : int, default=5000
        分块读取的行数，用于处理大文件时控制内存使用。
        
    preview_n : int, default=10
        预览显示的行数，便于在notebook中快速查看结果。
        
    print_preview : bool, default=True
        是否将前preview_n行打印到stderr（在notebook中也会显示）。
        
    返回
    ----
    str
        输出汇总文件的路径。
        
    输出格式
    --------
    生成的TSV文件包含以下列：
    - ID: 变体标识符
    - n0: 基因型为0的样本数量
    - n1: 基因型为1的样本数量  
    - n2: 基因型为2的样本数量
    - nmiss: 缺失基因型('.')的样本数量
    - ncalled: 有效基因型的样本数量（nsamples - nmiss）
    - nsamples: 总样本数量
    
    使用场景
    --------
    - 质量控制：检查各变体的基因型分布
    - 预筛选：识别适合特定遗传模型分析的变体
    - 统计报告：生成基因型分布摘要
    
    性能优化
    --------
    - 使用向量化计算提高处理速度
    - 分块读取避免大文件内存溢出
    - 追加写入模式处理超大数据集
    """
    import os
    import pandas as pd
    import numpy as np

    if not os.path.exists(gt_tsv_path):
        raise FileNotFoundError(f"GT 矩阵不存在：{gt_tsv_path}")

    if out_path is None:
        base = os.path.basename(gt_tsv_path)
        if base.endswith(".tsv"):
            base = base[:-4]
        out_path = os.path.join(os.path.dirname(os.path.abspath(gt_tsv_path)), f"{base}.gt_summary.tsv")

    _log_info(f"[GT-SUM] 输入矩阵：{gt_tsv_path}")
    _log_info(f"[GT-SUM] 输出汇总：{out_path}")

    # 写表头
    with open(out_path, "w") as fout:
        fout.write("ID\tn0\tn1\tn2\tnmiss\tncalled\tnsamples\n")

    total_rows = 0
    # 分块读取：首列为 ID，其余列为样本；全部按字符串读入，便于识别 '.'
    reader = pd.read_csv(gt_tsv_path, sep="\t", dtype=str, chunksize=chunksize)
    preview_buf = None

    for ichunk, df in enumerate(reader):
        if df.shape[1] < 2:
            raise ValueError("输入矩阵列数不足：必须至少包含 'ID' 与一个样本列")

        ids = df.iloc[:, 0].astype(str).values
        arr = df.iloc[:, 1:].astype(str).to_numpy()

        # 计数（向量化）
        c0 = (arr == '0').sum(axis=1)
        c1 = (arr == '1').sum(axis=1)
        c2 = (arr == '2').sum(axis=1)
        cm = (arr == '.').sum(axis=1)
        nsamples = arr.shape[1]
        ccalled = nsamples - cm

        out_df = pd.DataFrame({
            "ID": ids,
            "n0": c0,
            "n1": c1,
            "n2": c2,
            "nmiss": cm,
            "ncalled": ccalled,
            "nsamples": nsamples,
        })

        # 预览缓冲前 preview_n 行
        if ichunk == 0 and print_preview:
            preview_buf = out_df.head(preview_n)

        # 追加写出
        out_df.to_csv(out_path, sep="\t", index=False, header=False, mode="a")
        total_rows += out_df.shape[0]

    _log_info(f"[GT-SUM] 完成，累计变体数={total_rows} 行；每行含 n0/n1/n2/nmiss/ncalled/nsamples。")

    if print_preview and preview_buf is not None:
        _log_info(f"[GT-SUM] 预览前 {min(preview_n, len(preview_buf))} 行：")
        # 打印为制表符分隔，便于在 notebook 输出对齐
        preview_text = preview_buf.to_csv(sep="\t", index=False)
        # 打印到 stderr（notebook 也会显示）
        import sys as _sys
        print(preview_text.strip(), file=_sys.stderr)

    return out_path

def assemble_variant_gene_protein_table(
    variant_ids: List[str],
    gwas_path: str,
    pro_apt_path: str,
    out_path: Optional[str] = None,
    preview_n: int = 10,
    print_preview: bool = True,
) -> str:
    """
    整合变体、基因和蛋白质注释信息，生成综合注释表格。
    
    根据变体ID列表，从GWAS结果文件和蛋白质适配子(apt)文件中提取并整合
    相关注释信息，包括rsID、基因名称、蛋白质SeqID和UniProt信息。
    
    参数
    ----
    variant_ids : List[str]
        要查询的变体ID列表，格式通常为 "CHROM:POS:REF:ALT"。
        
    gwas_path : str  
        GWAS结果文件路径，必须包含列：
        - 'ID': 变体标识符（用于与variant_ids匹配）
        - 'rsID': dbSNP参考编号
        - 'Gene': 基因名称（可能包含格式如"{'GENE_NAME'}"的字符串）
        
    pro_apt_path : str
        蛋白质适配子注释文件路径，必须包含列：
        - 'EntrezGeneSymbol': 基因符号（用于与GWAS中的Gene匹配）
        - 'SeqId': 蛋白质序列标识符
        - 'UniProt': UniProt数据库标识符
        
    out_path : str, optional
        输出文件路径。如果未指定，将在gwas_path同目录下生成
        "variant_gene_protein.tsv"文件。
        
    preview_n : int, default=10
        预览显示的行数。
        
    print_preview : bool, default=True
        是否打印预览结果到stderr。
        
    返回
    ----
    str
        输出文件路径。
        
    输出格式  
    --------
    生成的TSV文件包含以下列：
    - ID: 变体标识符（来自variant_ids）
    - rsID: dbSNP编号（来自GWAS文件，未匹配时为'na'）
    - Gene: 基因名称（清理后的格式，未匹配时为'na'）
    - SeqID: 蛋白质序列ID（多个ID用';'分隔，未匹配时为'na'）
    - UniProt: UniProt标识符（多个ID用';'分隔，未匹配时为'na'）
    - Warn: 警告信息（记录匹配过程中的问题，如缺失或未找到）
    
    数据处理逻辑
    -----------
    1. 基因名称清理：从格式如"{'GENE_NAME'}"中提取实际基因名
    2. 模糊匹配：支持基因名的大小写不敏感匹配
    3. 多重匹配：一个基因可能对应多个蛋白质，结果用';'分隔
    4. 错误记录：所有匹配问题都记录在Warn列中
    
    使用场景
    --------
    - GWAS后分析：为显著变体添加功能注释
    - 蛋白质关联研究：连接遗传变异与蛋白质表达
    - 数据整合：合并来自不同数据库的注释信息
    """
    import os
    import re
    import pandas as pd

    if not variant_ids:
        raise ValueError("variant_ids 不能为空")

    if not os.path.exists(gwas_path):
        raise FileNotFoundError(f"GWAS 文件不存在：{gwas_path}")
    if not os.path.exists(pro_apt_path):
        raise FileNotFoundError(f"protein apt 文件不存在：{pro_apt_path}")

    # 输出路径
    if out_path is None:
        out_dir = os.path.dirname(os.path.abspath(gwas_path))
        out_path = os.path.join(out_dir, "variant_gene_protein.tsv")

    _log_info(f"[AGG] 读取 GWAS：{gwas_path}")
    # 读取 gwas（优先 \t，其次 ,，最后自动推断）
    def _read_any(p):
        try:
            df = pd.read_csv(p, sep="\t", dtype=str)
            if df.shape[1] == 1:
                df = pd.read_csv(p, sep=",", dtype=str)
        except Exception:
            df = pd.read_csv(p, sep=None, dtype=str, engine="python")
        return df

    gwas = _read_any(gwas_path)
    gwas_cols = {c.lower(): c for c in gwas.columns}
    for need in ["id", "rsid", "gene"]:
        if need not in gwas_cols:
            raise KeyError(f"GWAS 缺少必要列：{need}")

    col_ID = gwas_cols["id"]
    col_rs = gwas_cols["rsid"]
    col_gene = gwas_cols["gene"]

    # 建 ID -> (rsID, Gene)
    gwas[col_ID] = gwas[col_ID].astype(str)
    gwas[col_rs] = gwas[col_rs].astype(str)
    gwas[col_gene] = gwas[col_gene].astype(str)
    gwas_map = {row[col_ID]: (row[col_rs], row[col_gene]) for _, row in gwas[[col_ID, col_rs, col_gene]].iterrows()}

    _log_info(f"[AGG] 读取 pro_apt：{pro_apt_path}")
    apt = _read_any(pro_apt_path)
    apt_cols = {c.lower(): c for c in apt.columns}
    for need in ["entrezgenesymbol", "seqid", "uniprot"]:
        if need not in apt_cols:
            raise KeyError(f"pro_apt 缺少必要列：{need}")

    col_gene_sym = apt_cols["entrezgenesymbol"]
    col_seqid = apt_cols["seqid"]
    col_uniprot = apt_cols["uniprot"]

    # GeneSymbol（不区分大小写） -> {SeqId集合, UniProt集合}
    apt[col_gene_sym] = apt[col_gene_sym].astype(str)
    apt[col_seqid] = apt[col_seqid].astype(str)
    apt[col_uniprot] = apt[col_uniprot].astype(str)

    gene_map: Dict[str, Dict[str, set]] = {}
    for _, r in apt[[col_gene_sym, col_seqid, col_uniprot]].iterrows():
        g = (r[col_gene_sym] or "").strip()
        if not g or g.lower() == "nan":
            continue
        k = g.upper()
        if k not in gene_map:
            gene_map[k] = {"SeqId": set(), "UniProt": set()}
        if r[col_seqid] and r[col_seqid].lower() != "nan":
            gene_map[k]["SeqId"].add(r[col_seqid])
        if r[col_uniprot] and r[col_uniprot].lower() != "nan":
            gene_map[k]["UniProt"].add(r[col_uniprot])

    def _clean_gene_field(val: str) -> List[str]:
        """将 Gene 列值规范为基因列表。
        支持形式："A"、"{'A'}"、"['A','B']"、"{\"A\", \"B\"}" 等；
        先尝试 ast.literal_eval，再回退到宽松正则；保持顺序并去重。
        """
        if val is None:
            return []
        s = str(val).strip()
        if s == "" or s.lower() == "nan":
            return []
        # 1) 优先尝试安全解析（能处理大多数形如集合/列表的字符串）
        try:
            import ast
            parsed = ast.literal_eval(s)
            if isinstance(parsed, (list, tuple, set)):
                items = [str(x).strip() for x in parsed if str(x).strip()]
                # 去重保持顺序
                seen = set(); out = []
                for x in items:
                    if x not in seen:
                        out.append(x); seen.add(x)
                return out
            elif isinstance(parsed, str):
                v = parsed.strip()
                return [v] if v else []
        except Exception:
            pass
        # 2) 回退：删除各种括号和引号（包含排版引号）后按逗号切分
        # 包含 ASCII quote '" 以及常见排版引号 ’ ‘ “ ”，以及花括号/方括号/空格
        s2 = re.sub(r"[\{\}\[\]\s'\"’‘“”]", "", s)
        if s2 == "":
            return []
        parts = [p for p in s2.split(",") if p]
        seen = set(); out = []
        for p in parts:
            if p not in seen:
                out.append(p); seen.add(p)
        return out

    def _uniq_keep(seq: List[str]) -> List[str]:
        seen = set()
        out = []
        for x in seq:
            if x not in seen:
                out.append(x)
                seen.add(x)
        return out

    rows = []
    for vid in variant_ids:
        warn_msgs = []
        rs, gene_raw = ("na", None)
        tup = gwas_map.get(vid)
        if tup is None:
            warn_msgs.append("ID_not_in_gwas")
        else:
            rs, gene_raw = tup
            if not isinstance(rs, str) or rs.strip() == "" or rs.lower() == "nan":
                rs = "na"

        gene_list = _clean_gene_field(gene_raw)
        gene_display = ";".join(gene_list) if gene_list else "na"
        if not gene_list:
            warn_msgs.append("no_gene_in_gwas")

        seqids_accum: List[str] = []
        uniprot_accum: List[str] = []
        for g in gene_list:
            info = gene_map.get(g.upper())
            if not info:
                warn_msgs.append(f"no_pro_apt_for_gene:{g.strip()}")
                continue
            if info["SeqId"]:
                seqids_accum.extend(list(info["SeqId"]))
            if info["UniProt"]:
                uniprot_accum.extend(list(info["UniProt"]))

        seqids_out = ";".join(_uniq_keep([s for s in seqids_accum if s and s.lower() != "nan"])) if seqids_accum else "na"
        uniprot_out = ";".join(_uniq_keep([u for u in uniprot_accum if u and u.lower() != "nan"])) if uniprot_accum else "na"
        warn_out = ";".join(warn_msgs) if warn_msgs else ""

        rows.append({
            "ID": vid,
            "rsID": rs,
            "Gene": gene_display,
            "SeqID": seqids_out,
            "UniProt": uniprot_out,
            "Warn": warn_out
        })

    out_df = pd.DataFrame(rows, columns=["ID", "rsID", "Gene", "SeqID", "UniProt", "Warn"])
    os.makedirs(os.path.dirname(os.path.abspath(out_path)), exist_ok=True)
    out_df.to_csv(out_path, sep="\t", index=False, na_rep="na")
    _log_info(f"[AGG] 已写出：{out_path}（{out_df.shape[0]} 行）")

    if print_preview:
        import sys as _sys
        _log_info(f"[AGG] 预览前 {min(preview_n, len(out_df))} 行：")
        print(out_df.head(preview_n).to_csv(sep="\t", index=False).strip(), file=_sys.stderr)

    return out_path


# ------------------------------
# 每个变体的蛋白质表达 vs 基因型箱线图
# ------------------------------

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
plt.style.use('default')

from textwrap import fill as _tw_fill

import re

# --- 安全导入统计和格式化模块 ---
try:
    from scipy import stats as _sp_stats
except Exception:
    _sp_stats = None

_DEF_FP = "{:.3g}"

# --- 标准化可选字符串字段显示的辅助函数（Gene、UniProt等） ---
def _norm_opt_str(x: Any) -> str:
    """标准化可选字符串字段用于显示。将None、空字符串或'nan'转换为'na'。"""
    s = str(x).strip() if x is not None else ""
    if s == "" or s.lower() == "nan":
        return "na"
    return s


def _clean_sample_id(sample_id: Any) -> Optional[str]:
    """
    标准化样本ID以改善表达量和基因型数据之间的对齐。
    
    处理步骤：
    - 去除空格
    - 移除尾部的 `_day<数字><任意内容>`
    - 移除尾部的 `.数字`（如合并后缀中的".1"）
    - 如果为空则返回None
    
    参数：
        sample_id: 原始样本ID
        
    返回：
        清理后的样本ID，如果为空则返回None
    """
    if sample_id is None:
        return None
    s = str(sample_id).strip()
    if s == "" or s.lower() == "nan":
        return None
    s = s.replace(" ", "")
    # 移除 _day* 后缀
    s = re.sub(r"_day\d+.*$", "", s, flags=re.IGNORECASE)
    # 移除 .<数字> 后缀，如 .1
    s = re.sub(r"\.(\d+)$", "", s)
    return s


def _ensure_df(obj, sep="\t", dtype=str):
    """辅助函数：接受DataFrame或文件路径，返回DataFrame。"""
    if isinstance(obj, pd.DataFrame):
        return obj.copy()
    if isinstance(obj, str):
        # 先尝试TSV，然后CSV，最后使用python引擎自动检测
        try:
            df = pd.read_csv(obj, sep=sep, dtype=dtype)
            if df.shape[1] == 1:
                df = pd.read_csv(obj, sep=",", dtype=dtype)
        except Exception:
            df = pd.read_csv(obj, sep=None, dtype=dtype, engine="python")
        return df
    raise TypeError("期望DataFrame或文件路径")


def _strip_day_suffix(sample_id: str) -> str:
    """如果存在，移除样本ID尾部的 `_day*` 后缀。"""
    if sample_id is None:
        return sample_id
    s = str(sample_id)
    key = "_day"
    pos = s.find(key)
    return s[:pos] if pos >= 0 else s


def _normalize_model(model: str) -> str:
    """标准化遗传模型字符串，必须是ADD、DOM或REC之一。"""
    if not isinstance(model, str):
        raise ValueError("model必须是字符串，可选值：{ADD, DOM, REC}")
    m = model.strip().upper()
    if m not in {"ADD", "DOM", "REC"}:
        raise ValueError("model必须是以下之一：ADD, DOM, REC")
    return m


def _select_seqid(seqid_field: str) -> str:
    """从SeqID字段中选择第一个可用的ID。字段可能是'na'或用';'分隔的多个ID。"""
    if seqid_field is None:
        return "na"
    s = str(seqid_field).strip()
    if s == "" or s.lower() == "na" or s.lower() == "nan":
        return "na"
    # 按';'分割并取第一个非空部分
    parts = [p.strip() for p in s.split(";") if p.strip()]
    return parts[0] if parts else "na"


def _get_group_labels_and_values(expr_s: pd.Series, gt_s: pd.Series, model: str):
    """
    根据遗传模型构建用于绘图的分组。
    
    参数：
        expr_s: 蛋白质表达数据Series
        gt_s: 基因型数据Series
        model: 遗传模型（'ADD'/'DOM'/'REC'）
        
    返回：
        Tuple[List[str], List[np.ndarray], Dict[str,int]]: 
        (标签列表, 数值数组列表, 计数字典)
        
    注意：
        缺失基因型'.'会被忽略。
    """
    model = _normalize_model(model)
    # 只保留基因型和表达量都非缺失的样本
    df = pd.DataFrame({"expr": expr_s, "gt": gt_s}).dropna()
    df = df[df["gt"].astype(str) != "."]
    # 将基因型强制转换为整数
    try:
        df["gt"] = df["gt"].astype(int)
    except Exception:
        df = df[pd.to_numeric(df["gt"], errors="coerce").notna()]
        df["gt"] = df["gt"].astype(int)

    if model == "ADD":  # 加性模型
        groups = {
            "GT=0": df.loc[df["gt"] == 0, "expr"].to_numpy(),
            "GT=1": df.loc[df["gt"] == 1, "expr"].to_numpy(),
            "GT=2": df.loc[df["gt"] == 2, "expr"].to_numpy(),
        }
        labels = ["GT=0", "GT=1", "GT=2"]
    elif model == "DOM":  # 显性模型
        groups = {
            "GT=0": df.loc[df["gt"] == 0, "expr"].to_numpy(),
            "GT=1|2": df.loc[df["gt"].isin([1, 2]), "expr"].to_numpy(),
        }
        labels = ["GT=0", "GT=1|2"]
    else:  # REC - 隐性模型
        groups = {
            "GT=0|1": df.loc[df["gt"].isin([0, 1]), "expr"].to_numpy(),
            "GT=2": df.loc[df["gt"] == 2, "expr"].to_numpy(),
        }
        labels = ["GT=0|1", "GT=2"]

    values = [groups[lbl] for lbl in labels]
    counts = {lbl: int(len(groups[lbl])) for lbl in labels}
    return labels, values, counts


# --- 成对统计分析计算辅助函数 ---
def _compute_pairwise_stats(labels: List[str], values: List[np.ndarray]) -> pd.DataFrame:
    """
    计算非空组之间的成对统计分析。
    
    参数：
        labels: 组标签列表
        values: 对应的数值数组列表
        
    返回：
        包含以下列的DataFrame：
        ['group_a','group_b','n_a','n_b','mean_a','mean_b','mean_diff',
         'median_a','median_b','median_diff','t_stat','t_p','mw_u','mw_p']
        
    注意：
        如果SciPy不可用，t检验和Mann-Whitney U检验字段将填充'na'。
    """
    rows = []
    k = len(labels)
    for i in range(k):
        ai = values[i]
        if ai is None or len(ai) == 0:
            continue
        for j in range(i+1, k):
            bj = values[j]
            if bj is None or len(bj) == 0:
                continue
            n_a = int(np.isfinite(ai).sum())
            n_b = int(np.isfinite(bj).sum())
            if n_a == 0 or n_b == 0:
                continue
            mean_a = float(np.nanmean(ai))
            mean_b = float(np.nanmean(bj))
            med_a = float(np.nanmedian(ai))
            med_b = float(np.nanmedian(bj))
            mean_diff = mean_a - mean_b
            med_diff = med_a - med_b
            t_stat = t_p = mw_u = mw_p = None
            if _sp_stats is not None:
                try:
                    t_stat, t_p = _sp_stats.ttest_ind(ai, bj, equal_var=False, nan_policy='omit')
                except Exception:
                    t_stat, t_p = None, None
                try:
                    # Use two-sided Mann–Whitney U; require finite values only
                    a_fin = np.asarray(ai, float)
                    b_fin = np.asarray(bj, float)
                    a_fin = a_fin[np.isfinite(a_fin)]
                    b_fin = b_fin[np.isfinite(b_fin)]
                    if len(a_fin) > 0 and len(b_fin) > 0:
                        mw_u, mw_p = _sp_stats.mannwhitneyu(a_fin, b_fin, alternative='two-sided')
                    else:
                        mw_u, mw_p = None, None
                except Exception:
                    mw_u, mw_p = None, None
            rows.append({
                'group_a': labels[i], 'group_b': labels[j],
                'n_a': n_a, 'n_b': n_b,
                'mean_a': mean_a, 'mean_b': mean_b, 'mean_diff': mean_diff,
                'median_a': med_a, 'median_b': med_b, 'median_diff': med_diff,
                't_stat': t_stat, 't_p': t_p,
                'mw_u': mw_u, 'mw_p': mw_p,
            })
    return pd.DataFrame(rows)


def _check_group_feasibility(n0: int, n1: int, n2: int, model: str) -> (bool, str, str):
    """
    检查基因型分组的可行性。
    
    参数：
        n0, n1, n2: 基因型0/1/2的样本数量（来自汇总统计）
        model: 遗传模型（'ADD'/'DOM'/'REC'）
        
    返回：
        Tuple[bool, str, str]: (是否可行, 中文原因, 英文原因)
    """
    model = _normalize_model(model)
    if model == "ADD":
        # 至少两个基因型非零
        nonzero = sum([n0 > 0, n1 > 0, n2 > 0])
        if nonzero >= 2:
            return True, "", ""
        return (
            False,
            f"ADD条件不满足：三种基因型中至少两组需要非空（n0={n0}, n1={n1}, n2={n2})",
            f"ADD not satisfied: at least two genotype groups must be non-empty (n0={n0}, n1={n1}, n2={n2}).",
        )
    if model == "DOM":
        if (n0 > 0) and ((n1 + n2) > 0):
            return True, "", ""
        return (
            False,
            f"DOM条件不满足：要求GT=0与(GT=1+GT=2)均非空（n0={n0}, n1+n2={n1+n2})",
            f"DOM not satisfied: both GT=0 and (GT=1|2) must be non-empty (n0={n0}, n1+n2={n1+n2}).",
        )
    # REC - 隐性模型
    if (n0 + n1) > 0 and (n2 > 0):
        return True, "", ""
    return (
        False,
        f"REC条件不满足：要求(GT=0+GT=1)与GT=2均非空（n0+n1={n0+n1}, n2={n2})",
        f"REC not satisfied: both (GT=0|1) and GT=2 must be non-empty (n0+n1={n0+n1}, n2={n2}).",
    )



# --- 绘图辅助函数：箱线图和均值 ---
def _draw_box_with_means(ax, labels: List[str], values: List[np.ndarray], counts: Dict[str,int], title: str, y_label: str):
    """
    绘制箱线图并叠加均值标记。
    
    参数：
        ax: matplotlib轴对象
        labels: 组标签列表
        values: 对应的数值数组列表
        counts: 各组的样本计数
        title: 图标题
        y_label: Y轴标签
        
    注意：
        此函数不进行数据转换。
    """
    bp = ax.boxplot(values, labels=[f"{lbl} (n={counts[lbl]})" for lbl in labels], showfliers=True)
    ax.set_title(title, fontsize=10, loc='center', weight='bold')
    ax.set_ylabel(y_label)
    ax.grid(True, axis="y", linestyle=":", alpha=0.4)
    # 用蓝色三角形叠加均值
    mean_plotted = False
    for i, arr in enumerate(values, start=1):
        if arr is None or len(arr) == 0:
            continue
        try:
            m = float(np.nanmean(arr))
        except Exception:
            continue
        if np.isfinite(m):
            ax.plot(i, m, marker='^', color='blue', markersize=7, linestyle='None', label=('Mean' if not mean_plotted else None))
            mean_plotted = True
    if mean_plotted:
        ax.legend(loc='best', frameon=False, fontsize=9)


# --- 统计表格渲染辅助函数 ---
def _render_stats_table(ax, stats_df: pd.DataFrame):
    """
    在指定的轴上渲染统计分析结果表格。
    
    参数：
        ax: matplotlib轴对象
        stats_df: 包含成对统计分析结果的DataFrame
    """
    ax.set_axis_off()
    if stats_df is None or stats_df.empty:
        ax.text(0.5, 0.5, "No pairwise stats (need \u2265 two non-empty groups)", ha='center', va='center')
        return

    # Display-friendly subset & formatting
    disp = stats_df.copy()
    disp_cols = [
        ('group_a','A'), ('group_b','B'),
        ('mean_diff','\u0394mean'), ('median_diff','\u0394median'),
        ('t_p','t p'), ('mw_p','MWU p'), ('n_a','nA'), ('n_b','nB')
    ]
    cols = [c for c,_ in disp_cols]
    heads = [h for _,h in disp_cols]

    def _fmt(x):
        if x is None:
            return 'na'
        try:
            if isinstance(x,str):
                return x
            if np.isnan(x):
                return 'na'
        except Exception:
            pass
        if isinstance(x,(int,np.integer)):
            return str(int(x))
        try:
            return _DEF_FP.format(float(x))
        except Exception:
            return str(x)

    cell_text = [[_fmt(v) for v in row] for row in disp[cols].itertuples(index=False, name=None)]

    # Auto sizing: make table occupy more width, slightly shorter height
    bbox = [0.0, 0.20, 1.0, 0.50]
    n_cols = len(heads)
    base_width = 1.0 / n_cols
    # Widen delta & p-value columns, shrink counts
    col_widths = []
    for h in heads:
        if 'Δ' in h or 'p' in h:
            col_widths.append(base_width * 1.5)
        elif h in ('nA','nB'):
            col_widths.append(base_width * 0.7)
        else:
            col_widths.append(base_width * 1.25)
    # Normalize widths
    total = sum(col_widths)
    col_widths = [w/total for w in col_widths]

    table = ax.table(
        cellText=cell_text,
        colLabels=heads,
        loc='center',
        colWidths=col_widths,
        cellLoc='center',
        bbox=bbox,
    )
    table.auto_set_font_size(False)
    table.set_fontsize(8)  # slightly larger font
    table.scale(1.2, 0.4)  # wider, a bit less tall

    # Header styling
    for j, _h in enumerate(heads):
        cell = table[0, j]
        cell.set_text_props(weight='bold')
        cell.set_facecolor('#f0f0f0')
        cell.set_edgecolor('0.55')
        cell.set_linewidth(0.8)

    # Body cell styling & alignment
    n_rows = len(cell_text)
    for i in range(1, n_rows+1):
        for j in range(len(heads)):
            cell = table[i, j]
            cell.set_edgecolor('0.7')
            cell.set_linewidth(0.6)
            # align numbers to right for deltas/p-values; names stay centered
            if j <= 1:
                cell._loc = 'center'
            else:
                cell._loc = 'right'

    # Footnote clarifying delta direction
    ax.text(0.5, 0.06, "\u0394 = A \u2212 B  (both mean and median)", ha='center', va='center', fontsize=9)


# --- 生成简洁"未绘制"页面的辅助函数 ---
def _render_not_plotted_page(pdf, page_title: str, message: str):
    """
    渲染一个简洁的"未绘制"页面，具有一致的版式。
    
    使用图形级别的文本和约束布局，确保标题和正文在PDF查看器中
    视觉居中且不被裁剪。
    
    参数：
        pdf: PdfPages对象
        page_title: 页面标题
        message: 消息内容
    """
    # 为保持一致性，匹配正常绘图页面的大小
    fig, ax = plt.subplots(figsize=(6, 4), constrained_layout=True)
    ax.axis('off')

    # 标题（粗体，居中，顶部留有充足边距）
    fig.text(0.5, 0.88, str(page_title), ha='center', va='center', fontsize=10, weight='bold')

    # 正文消息（自动换行并垂直居中）
    wrapped = _tw_fill(str(message), width=60)
    fig.text(0.5, 0.55, wrapped, ha='center', va='center', fontsize=9)

    # 可选：底部的淡色引导文字以保持视觉平衡
    # fig.text(0.5, 0.12, "This page is intentionally left blank.", ha='center', fontsize=8, alpha=0.45)

    # 按原样保存（不使用bbox_inches='tight'以避免意外裁剪）
    pdf.savefig(fig)
    plt.close(fig)


def plot_protein_boxplot_per_variant(
    variant_meta,
    summary_path,
    pro_ex_path,
    gt_mx,
    model: str,
    out_pdf: Optional[str] = None,
    log2_transform: bool = True,
) -> str:
    """
    为每个变体生成蛋白质表达箱线图PDF报告。
    
    根据变体注释信息、基因型矩阵和蛋白质表达数据，按照指定的遗传模型
    （ADD/DOM/REC）为每个变体生成一页箱线图，展示不同基因型组的蛋白质表达分布。
    
    参数
    ----
    variant_meta : DataFrame 或 str
        变体元数据，必须包含列 ['ID','SeqID']。
        仅当 SeqID != 'na' 时才会生成图表，否则跳过并记录到日志。
        可选列：'Gene', 'UniProt', 'rsID'（用于图表标题）。
        
    summary_path : str
        基因型计数汇总文件路径，由 summarize_gt_counts() 生成。
        必须包含列 ['ID','n0','n1','n2']，分别表示基因型0/1/2的样本数。
        
    pro_ex_path : DataFrame 或 str  
        蛋白质表达矩阵文件。
        - 行：样本ID（末尾的 `_day*` 后缀会被自动移除以匹配基因型数据）
        - 列：蛋白质SeqID（需与variant_meta.SeqID匹配）
        - 值：蛋白质表达水平
        
    gt_mx : DataFrame 或 str
        基因型矩阵文件。
        - 行：变体ID（需与variant_meta.ID对齐）
        - 列：样本ID
        - 值：基因型编码 {0,1,2,'.'} ，其中'.'表示缺失
        
    model : str
        遗传模型，可选值：
        - 'ADD': 加性模型（GT=0 vs GT=1 vs GT=2）
        - 'DOM': 显性模型（GT=0 vs GT=1|2）  
        - 'REC': 隐性模型（GT=0|1 vs GT=2）
        
    out_pdf : str, optional
        输出PDF文件路径。如果未指定，将根据model自动命名。
        
    log2_transform : bool, default=True
        是否对蛋白质表达数据进行log2转换。
        
    返回
    ----
    str
        生成的PDF文件路径。同时会在相同目录生成`.log.txt`日志文件。
        
    输出文件
    --------
    - PDF报告：每个变体一页，包含箱线图和统计分析表格
    - 日志文件：记录处理过程、跳过的变体和错误信息
    
    图表内容
    --------
    每页包含：
    1. 箱线图：显示不同基因型组的蛋白质表达分布
    2. 均值标记：用蓝色三角形标注各组均值
    3. 统计表格：组间比较的t检验和Mann-Whitney U检验结果
    4. 变体信息：ID、基因、UniProt、rsID等注释
    
    注意事项
    --------
    - 如果变体的SeqID为'na'，将跳过绘图
    - 如果基因型分组不满足指定模型的要求（如某些组为空），将跳过绘图
    - 缺失的表达数据和基因型'.'会被自动排除
    - 样本ID匹配时会自动处理常见的后缀差异
    """
    # Load inputs
    vmeta = _ensure_df(variant_meta)
    if not set(["ID", "SeqID"]).issubset(set(vmeta.columns)):
        raise KeyError("variant_meta 需要包含列: ID, SeqID")

    vmeta_cols = {c.lower(): c for c in vmeta.columns}
    col_gene_opt = vmeta_cols.get("gene")
    col_uniprot_opt = vmeta_cols.get("uniprot")
    col_rsid_opt = vmeta_cols.get("rsid")

    summ = _ensure_df(summary_path)
    summ_cols = {c.lower(): c for c in summ.columns}
    for need in ["id", "n0", "n1", "n2"]:
        if need not in summ_cols:
            raise KeyError("summary_path 缺少必要列：ID, n0, n1, n2")
    sid_col = summ_cols["id"]; sc0 = summ_cols["n0"]; sc1 = summ_cols["n1"]; sc2 = summ_cols["n2"]
    summ = summ[[sid_col, sc0, sc1, sc2]].copy()

    # Expression matrix: rows = samples (with _day* suffix), cols = SeqID
    ex_df = _ensure_df(pro_ex_path)
    # normalize index (samples)
    # Case A: If index already meaningful, keep; else try to use first column as index if it looks like sample IDs
    if ex_df.index.name is None or isinstance(ex_df.index, pd.RangeIndex) or (not ex_df.index.is_unique):
        first_col = ex_df.columns[0]
        # Heuristic: if first column is unique and not too many NaNs, use it as index
        if ex_df[first_col].notna().sum() >= len(ex_df) * 0.95 and ex_df[first_col].nunique(dropna=True) >= len(ex_df) * 0.95:
            ex_df = ex_df.set_index(first_col)
    # Clean sample IDs in index and strip `_day*` etc.
    ex_df.index = [ _clean_sample_id(x) for x in ex_df.index ]
    # Drop any rows that became None after cleaning
    ex_df = ex_df[~pd.isna(ex_df.index)]

    # Genotype matrix: rows = variant IDs, cols = samples
    gt_df = _ensure_df(gt_mx)
    # try to set index/columns if not labeled
    if "ID" in gt_df.columns:
        gt_df = gt_df.set_index("ID")
    gt_df.index = gt_df.index.astype(str)
    # Clean genotype sample column names as well
    gt_df.columns = [ _clean_sample_id(c) for c in gt_df.columns.astype(str) ]

    model = _normalize_model(model)

    # build summary map: ID -> (n0,n1,n2)
    cnt_map = {row[sid_col]: (int(row[sc0]), int(row[sc1]), int(row[sc2])) for _, row in summ.iterrows()}

    # Prepare output paths
    if out_pdf is None:
        out_dir = os.getcwd()
        out_pdf = os.path.join(out_dir, f"protein_genotype_boxplots.{model}.pdf")
    out_log = out_pdf + ".log"
    logs: List[str] = []
    logs.append(f"BEGIN plotting: model={model}, log2_transform={log2_transform}")

    # Start PDF
    with PdfPages(out_pdf) as pdf:
        for _, r in vmeta.iterrows():
            vid = str(r["ID"]).strip()
            seqid = _select_seqid(r["SeqID"])  # choose first SeqID if multiple

            gene_disp = _norm_opt_str(r[col_gene_opt]) if col_gene_opt else "na"
            uniprot_disp = _norm_opt_str(r[col_uniprot_opt]) if col_uniprot_opt else "na"
            rs_disp = _norm_opt_str(r[col_rsid_opt]) if col_rsid_opt else "na"
            _title_line1 = f"Variant: {vid} | rsID: {rs_disp} | SeqID: {seqid}"
            _title_line2 = f"Gene: {gene_disp} | UniProt: {uniprot_disp} | Model: {model}"
            page_title = _title_line1 + "\n" + _title_line2

            # 1) SeqID must not be 'na'
            if seqid == "na":
                logs.append(f"SKIP ID={vid}: SeqID 为 'na'")
                _render_not_plotted_page(pdf, page_title, "Not plotted: SeqID is 'na'")
                continue

            # 2) SeqID existence in expression matrix
            if seqid not in ex_df.columns:
                logs.append(f"SKIP ID={vid}: 表达矩阵中不存在 SeqID={seqid}")
                _render_not_plotted_page(pdf, page_title, f"Not plotted: SeqID not found in expression matrix (SeqID={seqid})")
                continue

            # 3) Counts feasibility by model
            if vid not in cnt_map:
                logs.append(f"SKIP ID={vid}: 未在 summary 中找到计数信息")
                _render_not_plotted_page(pdf, page_title, "Not plotted: n0/n1/n2 counts missing in summary")
                continue

            n0, n1, n2 = cnt_map[vid]
            ok, reason_cn, reason_en = _check_group_feasibility(n0, n1, n2, model)
            if not ok:
                logs.append(f"SKIP ID={vid}: {reason_cn}")
                _render_not_plotted_page(pdf, page_title, f"Not plotted: {reason_en}")
                continue

            # 4) Assemble aligned series (intersection of samples)
            # Ensure expression is numeric; coerce non-numeric to NaN
            raw_expr = ex_df[seqid]
            expr_s = pd.to_numeric(raw_expr, errors="coerce")
            # Log and handle negative values for log2
            neg_cnt = int((expr_s < 0).sum()) if expr_s.notna().any() else 0
            nan_cnt = int(expr_s.isna().sum())
            if log2_transform:
                # add a small epsilon if needed to avoid -inf; assume values >= 0 for log2
                pos_mask = expr_s > 0
                min_positive = float(expr_s[pos_mask].min()) if pos_mask.any() else 0.0
                eps = min(1e-9, min_positive * 0.5) if min_positive > 0 else 1e-9
                if neg_cnt > 0:
                    logs.append(f"ID={vid}: 表达含 {neg_cnt} 个负值，log2 前按 0 下限裁剪")
                # clip negatives to 0 to avoid NaN in log2; keep NaN as missing
                expr_s = np.log2(expr_s.clip(lower=0).astype(float) + eps)
            else:
                # no transform: keep numeric; negatives allowed
                expr_s = expr_s.astype(float)
            expr_s.index = expr_s.index.astype(str)
            if log2_transform:
                logs.append(f"ID={vid}: applied log2 transform to expression (with epsilon and non-negative clipping)")
            else:
                logs.append(f"ID={vid}: no log transform applied to expression")

            if vid not in gt_df.index:
                logs.append(f"SKIP ID={vid}: 在基因型矩阵中不存在该变体")
                _render_not_plotted_page(pdf, page_title, "Not plotted: variant not found in genotype matrix")
                continue

            gt_s = gt_df.loc[vid]
            # Align by sample intersection (after cleaning)
            expr_idx = expr_s.index.astype(str)
            gt_idx = gt_s.index.astype(str)
            common = expr_idx.intersection(gt_idx)
            if len(common) == 0:
                # write diagnostics to logs for easier debugging
                logs.append(
                    f"ID={vid}: 无交集诊断 → expr_n={len(expr_idx)}, gt_n={len(gt_idx)}, expr_example={list(expr_idx[:5])}, gt_example={list(gt_idx[:5])}"
                )
                _render_not_plotted_page(pdf, page_title, "Not plotted: no overlapping samples between expression and genotype")
                continue

            expr_s = expr_s.loc[common]
            gt_s = gt_s.loc[common]

            # 5) Build groups and compute pairwise stats
            labels, values, counts = _get_group_labels_and_values(expr_s, gt_s, model)
            stats_df = _compute_pairwise_stats(labels, values)

            # 6) Plot: boxplot + stats table side-by-side
            from matplotlib import gridspec as _gridspec
            fig = plt.figure(figsize=(10.4, 4.2))
            gs = _gridspec.GridSpec(1, 2, width_ratios=[2.8, 3.2], wspace=0.26)
            ax_box = fig.add_subplot(gs[0, 0])
            ax_tbl = fig.add_subplot(gs[0, 1])
            _draw_box_with_means(
                ax_box, labels, values, counts, title=page_title,
                y_label=("Protein expression (log2)" if log2_transform else "Protein expression")
            )
            _render_stats_table(ax_tbl, stats_df)
            pdf.savefig(fig)
            plt.close(fig)

            # collect stats for optional TSV export
            if 'all_stats_rows' not in locals():
                all_stats_rows = []
            if stats_df is not None and not stats_df.empty:
                sdf = stats_df.copy()
                sdf.insert(0, 'ID', vid)
                sdf.insert(1, 'SeqID', seqid)
                try:
                    sdf.insert(2, 'rsID', rs_disp)
                except Exception:
                    pass
                all_stats_rows.append(sdf)

    # Optional: write aggregated pairwise stats to TSV alongside the PDF
    try:
        if 'all_stats_rows' in locals() and len(all_stats_rows) > 0:
            stats_out = out_pdf + '.pairwise_stats.tsv'
            all_stats_df = pd.concat(all_stats_rows, ignore_index=True)
            all_stats_df.to_csv(stats_out, sep='\t', index=False)
            logs.append(f"Wrote pairwise stats: {stats_out} (rows={all_stats_df.shape[0]})")
    except Exception as _e:
        logs.append(f"Pairwise stats export failed: {_e}")

    # append END marker before writing logs
    logs.append(f"END plotting: output={out_pdf}, pages={len([1 for _ in vmeta.iterrows()])}")
    # write log file
    try:
        with open(out_log, "w", encoding="utf-8") as f:
            for line in logs:
                f.write(line + "\n")
        _log_info(f"[PLOT] 日志写出：{out_log}，共 {len(logs)} 条记录")
    except Exception as e:
        _log_warn(f"写出日志失败：{e}")

    _log_info(f"[PLOT] PDF 生成完成：{out_pdf}")
    return out_pdf
