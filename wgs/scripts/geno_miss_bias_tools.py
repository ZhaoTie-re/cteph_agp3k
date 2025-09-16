#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import sys
import shutil
import subprocess
import tempfile
import time
import shlex
from typing import Tuple, Optional, List, Dict, Any

import pandas as pd

# Default paths for plink2 and bcftools
DEFAULT_PLINK2 = "/home/b/b37974/plink2"
DEFAULT_BCFTOOLS = "/home/b/b37974/bcftools/bcftools"

# ---- Lightweight logging (stderr + optional file in work_dir) ----
from datetime import datetime
LOG_FILE: Optional[str] = None  # will be set in build_aligned_geno_mats()

def _now() -> str:
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")

def _emit(level: str, msg: str):
    line = f"[{_now()}][{level}] {msg}"
    # to stderr
    print(line, file=sys.stderr, flush=True)
    # to file (best-effort)
    if LOG_FILE:
        try:
            with open(LOG_FILE, "a") as lf:
                lf.write(line + "\n")
        except Exception:
            pass

def log_info(msg: str): _emit("INFO", msg)

def log_warn(msg: str): _emit("WARN", msg)

def log_err(msg: str):  _emit("ERROR", msg)

def _run_cmd(cmd, cwd=None):
    """Run external command, raise on failure, log tails of stdout/stderr.
    Accepts list or str.
    """
    if isinstance(cmd, str):
        shell = True
        cmd_disp = cmd
    else:
        shell = False
        cmd_disp = " ".join(cmd)
    log_info(f"执行命令：{cmd_disp}")
    try:
        ret = subprocess.run(cmd, cwd=cwd, shell=shell, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    except Exception as e:
        log_err(f"命令启动失败：{e}")
        raise
    tail_out = "\n".join((ret.stdout or "").strip().splitlines()[-10:])
    tail_err = "\n".join((ret.stderr or "").strip().splitlines()[-20:])
    if tail_out:
        log_info(f"[stdout-tail]\n{tail_out}")
    if ret.returncode != 0:
        if tail_err:
            log_err(f"[stderr-tail]\n{tail_err}")
        raise RuntimeError(f"命令执行失败（退出码 {ret.returncode}）：{cmd_disp}")
    if tail_err:
        # some tools write warnings to stderr even on success
        log_warn(f"[stderr-tail]\n{tail_err}")
    return ret


def _ensure_dir(path: str):
    os.makedirs(path, exist_ok=True)


def _vcf_gt_to012(gt: str) -> Optional[int]:
    """将 VCF GT 字符串映射为 0/1/2 或 '.'（缺失）"""
    if gt in ('.', './.', '.|.'):
        return '.'
    if gt in ('0/0', '0|0'):
        return 0
    if gt in ('0/1', '1/0', '0|1', '1|0'):
        return 1
    if gt in ('1/1', '1|1'):
        return 2
    # 其他非常见编码：保守按缺失处理
    return '.'


def _run_pool(cmds: List[List[str]], max_parallel: int):
    """并行运行多个命令，最多同时运行 max_parallel 个进程，等待所有完成。"""
    procs = []
    cmds_iter = iter(cmds)
    active = 0
    try:
        # 启动初始进程
        for _ in range(min(max_parallel, len(cmds))):
            cmd = next(cmds_iter)
            log_info(f"[POOL] 启动：{' '.join(cmd)}")
            p = subprocess.Popen(cmd)
            procs.append(p)
            active += 1
        # 循环监控，启动新进程
        while active > 0:
            for i, p in enumerate(procs):
                if p is not None and p.poll() is not None:
                    if p.returncode != 0:
                        raise RuntimeError(f"Parallel command failed (exit {p.returncode}): {' '.join(cmds[i])}")
                    procs[i] = None
                    active -= 1
                    try:
                        cmd = next(cmds_iter)
                        log_info(f"[POOL] 启动：{' '.join(cmd)}")
                        p_new = subprocess.Popen(cmd)
                        procs[i] = p_new
                        active += 1
                    except StopIteration:
                        pass
            time.sleep(0.1)
    except Exception as e:
        # 尝试杀死所有子进程
        for p in procs:
            if p is not None and p.poll() is None:
                try:
                    p.kill()
                except Exception:
                    pass
        raise e


# 核心函数

def build_aligned_geno_mats(
    plink_prefix: str,
    vcf_path: str,
    work_dir: Optional[str] = None,
    keep_temp: bool = False,
    output_prefix: Optional[str] = None,
    plink_threads: int = 8,
    bcftools_threads: int = 8,
    chrom: Optional[str] = None,
    n_chunks: Optional[int] = None,
    chunk_size: Optional[int] = None,
    max_parallel: int = 4,
    return_paths_only: bool = True,
) -> Tuple[str, str]:
    """
    生成对齐的基因型矩阵（PLINK vs VCF），数值统一为 0/1/2 或 '.'（缺失）。
    假设两侧变体 ID 均为 CHR:POS:REF:ALT（本函数不在此修改 ID）。

    处理流程（概览）：
      1) 从 `plink_prefix`.fam 导出样本顺序（仅 IID）与表型；
      2) 从 `.bim` 第二列解析 ALT，构造 `--export-allele` 映射（确保 PLINK 导出按 ALT 计数）；
      3) 构建 PLINK 与 VCF 变体 ID 的交集；
         - VCF 侧 ID 提取策略：先尝试 `-r <chrom>`，不行试 `-r chr<chrom>`，仍不行则不加 `-r`；
      4) 将交集按 `n_chunks`/`chunk_size` 切成多个分块；
      5) 并行导出：
         - PLINK：`--export A-transpose`，随后清洗为 `ID + IID`、整数 0/1/2，缺失 '.'；
         - VCF：`bcftools query` 导出 GT 并映射为 0/1/2 或 '.'，列为 IID；
      6) 写出分块清单（manifest），并按顺序拼接为最终矩阵：
         - `<work_dir>/chr{chrom_or_ALL}.pre.gt.tsv`
         - `<work_dir>/chr{chrom_or_ALL}.post.gt.tsv`
      7) 日志：所有进度/告警/错误信息同步写入 `<work_dir>/geno_miss_bias.prepare.log`，同时打印到 stderr。

    参数（常用）：
      plink_prefix      : str  PLINK 二进制文件前缀（.bed/.bim/.fam）
      vcf_path          : str  输入 VCF.gz（需配套 .tbi 或 .csi 索引）
      work_dir          : str  工作目录（默认=当前工作目录）；缓存写入 `<work_dir>/tmp/`
      keep_temp         : bool 是否在结束时删除 `<work_dir>/tmp/`
      plink_threads     : int  plink2 线程数
      bcftools_threads  : int  bcftools 线程数
      chrom             : str  仅处理该染色体（如 '22'、'chr22'、'X'；命名不一致会自动回退策略）
      n_chunks          : int  分块数量（优先于 `chunk_size`）
      chunk_size        : int  每块行数（当未指定 `n_chunks` 时生效）
      max_parallel      : int  最大并行任务数

    返回：
      (pre_path, post_path) 两个字符串
        - pre_path  : `<work_dir>/chr{chrom_or_ALL}.pre.gt.tsv`
        - post_path : `<work_dir>/chr{chrom_or_ALL}.post.gt.tsv`

    重要说明：
      - 计数等位基因为 ALT（由 `.bim` 第二列解析），PLINK 与 VCF 数值语义对齐；
      - 样本列名仅为 IID（去掉 FID_ 前缀）；缺失一律输出为 '.'；
      - 交集为 0 时会导出少量样例 ID 用于人工比对；
      - 运行日志保存在 `<work_dir>/geno_miss_bias.prepare.log`，建议先查看日志定位问题（如 chr 前缀、REF/ALT 顺序、VCF 空 ID 等）。
    """
    # ---------- 路径与工作目录 ----------
    plink_bin = DEFAULT_PLINK2
    bcftools_bin = DEFAULT_BCFTOOLS
    # 工作目录：默认使用当前工作目录；缓存统一写入 <work_dir>/tmp
    base_tmp = work_dir if work_dir else os.getcwd()
    _ensure_dir(base_tmp)

    cache_dir = os.path.join(base_tmp, "tmp")
    _ensure_dir(cache_dir)
    log_info(f"工作目录：{base_tmp}；缓存目录：{cache_dir}")
    # initialize log file in work_dir
    global LOG_FILE
    LOG_FILE = os.path.join(base_tmp, "geno_miss_bias.prepare.log")
    try:
        with open(LOG_FILE, "a") as lf:
            lf.write("="*72 + "\n")
            lf.write(f"[{_now()}][INFO] 开始执行 build_aligned_geno_mats\n")
            lf.write(f"[{_now()}][INFO] 工作目录={base_tmp}；缓存目录={cache_dir}\n")
            lf.write(f"[{_now()}][INFO] plink2路径={DEFAULT_PLINK2}；bcftools路径={DEFAULT_BCFTOOLS}\n")
    except Exception:
        pass
    log_info(f"日志文件：{LOG_FILE}")
    # preflight checks
    for exe, name in [(DEFAULT_PLINK2, "plink2"), (DEFAULT_BCFTOOLS, "bcftools")]:
        if not os.path.exists(exe):
            log_err(f"未找到可执行文件 {name} ：'{exe}'")
            raise FileNotFoundError(f"{name} not found at '{exe}'")
        if not os.access(exe, os.X_OK):
            log_err(f"{name} 不可执行：'{exe}'")
            raise PermissionError(f"{name} is not executable: '{exe}'")
    if not (os.path.exists(vcf_path) and (os.path.exists(vcf_path + '.tbi') or os.path.exists(vcf_path + '.csi'))):
        log_warn("未发现 VCF 索引(.tbi/.csi)，bcftools 可能会很慢或失败。")

    def _p(*names: str) -> str:
        """在缓存子目录下拼出路径"""
        return os.path.join(cache_dir, *names)

    # 输出/中间文件前缀
    postvars_prefix = _p("postvars")
    sample_order_path = _p("sample.order.txt")
    sample_pheno_path = _p("sample.pheno.tsv")

    try:
        # ---------- 1) 样本顺序 & 表型 ----------
        fam_path = plink_prefix + ".fam"
        if not os.path.exists(fam_path):
            raise FileNotFoundError(f"FAM not found: {fam_path}")

        _run_cmd(f"awk '{{print $2}}' {fam_path} > {sample_order_path}")
        _run_cmd(f"awk '{{print $2\"\\t\"$6}}' {fam_path} > {sample_pheno_path}")
        log_info(f"样本顺序已写出：{sample_order_path}；样本数：{sum(1 for _ in open(sample_order_path))}")
        try:
            vcfl = subprocess.check_output(["bash", "-lc", f"{bcftools_bin} query -l {vcf_path} | wc -l"], text=True).strip()
            log_info(f"VCF 表头样本数：{vcfl}")
        except Exception as e:
            log_warn(f"统计 VCF 样本数失败：{e}")

        # 构建 FID_IID → IID 的映射，以正确处理 IID 本身含有下划线的情况
        fid_iid_to_iid: Dict[str, str] = {}
        iid_set: set = set()
        with open(fam_path, "r") as f:
            for line in f:
                if not line.strip():
                    continue
                parts = line.rstrip("\n").split()
                if len(parts) < 2:
                    continue
                fid, iid = parts[0], parts[1]
                fid_iid_to_iid[f"{fid}_{iid}"] = iid
                iid_set.add(iid)
        log_info(f"已建立 FID_IID→IID 映射：{len(fid_iid_to_iid)} 条；IID 数量：{len(iid_set)}")

        # ---------- 2) 解析 ALT ----------
        bim_path = plink_prefix + ".bim"
        if not os.path.exists(bim_path):
            raise FileNotFoundError(f"BIM not found: {bim_path}")
        allele_map_path = _p("export_allele.from_bim_id_ALT.txt")
        # 解析 VAR_ID=CHR:POS:REF:ALT，取 ALT (第4段)。若格式异常则报错。
        _run_cmd(
            "awk -F'\\t' '{split($2,a,\":\"); if (length(a)>=4) print $2\"\\t\"a[4]; else {print \"[ERR] Bad BIM ID: \"$2 > \"/dev/stderr\"; exit 1}}' "
            + bim_path + " > " + allele_map_path
        )
        log_info(f"已从 BIM 解析 ALT（基于第二列 CHR:POS:REF:ALT）：{allele_map_path}")

        # ---------- 3) 生成交集变体列表 ----------
        post_snplist_path = postvars_prefix + ".snplist"
        snplist_cmd = [
            plink_bin,
            "--bfile", plink_prefix,
            "--write-snplist",
            "--threads", str(plink_threads),
            "--out", postvars_prefix
        ]
        if chrom is not None:
            snplist_cmd.extend(["--chr", chrom])
        _run_cmd(snplist_cmd)

        # --- Build pre.id list with region (if given), then normalize both sides and compute intersection ---
        post_ids_norm = _p("post.snplist.norm.txt")
        pre_ids_txt = _p("pre.ids.txt")
        pre_ids_norm = _p("pre.ids.norm.txt")

        # Try region filter strategies: raw chrom, then 'chr'+chrom, then no -r
        def try_region(region_arg: Optional[str], label: str) -> bool:
            try:
                if region_arg:
                    cmd = f"{bcftools_bin} view --threads {bcftools_threads} -Ou -r {region_arg} {vcf_path} | {bcftools_bin} query -f '%ID\\n' > {pre_ids_txt}"
                else:
                    cmd = f"{bcftools_bin} view --threads {bcftools_threads} -Ou {vcf_path} | {bcftools_bin} query -f '%ID\\n' > {pre_ids_txt}"
                _run_cmd(f"bash -lc \"{cmd}\"")
                _run_cmd(f"bash -lc \"sed -e 's/\\r$//' {pre_ids_txt} | grep -v '^$' | LC_ALL=C sort -u > {pre_ids_norm}\"")
                # quick non-empty check
                out = subprocess.check_output(["wc", "-l", pre_ids_norm], text=True)
                n = int(out.strip().split()[0])
                if n > 0:
                    log_info(f"使用 {label} 提取到 VCF 变体ID 数量：{n}")
                    return True
                else:
                    log_warn(f"使用 {label} 提取的 VCF 变体ID 为空")
                    return False
            except Exception as e:
                log_warn(f"使用 {label} 提取 VCF 变体ID 失败：{e}")
                return False

        ok = False
        if chrom is not None:
            # try plain chrom
            ok = try_region(chrom, f"-r {chrom}")
            if not ok and not chrom.startswith("chr"):
                ok = try_region("chr"+chrom, f"-r chr{chrom}")
        if not ok:
            log_warn("回退为全局提取 VCF 变体ID（不加 -r）")
            ok = try_region(None, "no -r")
        if not ok:
            raise RuntimeError("即使不加 -r 也未能提取到任何 VCF 变体ID")

        # Normalize line endings/whitespace and sort-unique post IDs
        _run_cmd(f"bash -lc \"sed -e 's/\\r$//' {post_snplist_path} | grep -v '^$' | LC_ALL=C sort -u > {post_ids_norm}\"")

        # Initial intersection
        intersect_ids = _p("intersect.ids")
        comm_cmd = f"LC_ALL=C comm -12 {post_ids_norm} {pre_ids_norm} > {intersect_ids}"
        try:
            _run_cmd(f"bash -lc \"{comm_cmd}\"")
        except RuntimeError:
            # Fallback: handle potential locale/ordering issues without requiring sorted inputs
            log_warn("comm -12 失败（可能与区域设置/排序有关），回退为 grep -Fxf + sort -u。")
            fallback_cmd = f"grep -Fxf {post_ids_norm} {pre_ids_norm} | LC_ALL=C sort -u > {intersect_ids}"
            _run_cmd(f"bash -lc \"{fallback_cmd}\"")

        # Count lines for diagnostics
        def _wc_l(path: str) -> int:
            try:
                out = subprocess.check_output(["wc", "-l", path], text=True)
                return int(out.strip().split()[0])
            except Exception:
                return -1

        post_n = _wc_l(post_ids_norm)
        pre_n  = _wc_l(pre_ids_norm)
        inter_n = _wc_l(intersect_ids)
        log_info(f"PLINK 变体ID数：{post_n}；VCF 变体ID数（可能含区域过滤）：{pre_n}；交集：{inter_n}")
        if post_n == 0:
            log_err("PLINK 侧 snplist 为空。请检查 --bfile 是否正确、--chr 是否匹配 BIM 染色体编码。")
        if pre_n == 0:
            log_err("VCF 侧 ID 列表为空。常见原因：-r 区域不匹配（22 vs chr22）、VCF 缺少 ID 字段或索引损坏。")
        if inter_n == 0:
            log_warn("Empty intersection: likely REF/ALT order or chr prefix mismatch; VCF may also contain '.' IDs.")
            # Dump a few example IDs from both sides to help debugging
            diag_post = _p("post.ids.head.txt")
            diag_pre  = _p("pre.ids.head.txt")
            _run_cmd(f"bash -lc \"head -n 10 {post_ids_norm} > {diag_post}\"")
            _run_cmd(f"bash -lc \"head -n 10 {pre_ids_norm} > {diag_pre}\"")
            raise RuntimeError("No intersecting variants found between PLINK and VCF. Check naming (chr prefix), REF/ALT order in IDs, and whether VCF has missing '.' IDs. See head files in work_dir for examples.")

        # ---------- 4) 分块切割 ----------
        # 计算总变体数
        try:
            wc_out = subprocess.check_output(["wc", "-l", intersect_ids], text=True)
            total_variants = int(wc_out.strip().split()[0])
        except Exception as e:
            raise RuntimeError(f"Failed to count lines in intersect.ids: {e}")

        if total_variants == 0:
            raise RuntimeError("No intersecting variants found between PLINK and VCF.")

        # 计算分块数量和每块大小
        if n_chunks is None and chunk_size is None:
            n_chunks = 8
        if n_chunks is not None:
            n_chunks = max(1, n_chunks)
            lines_per_chunk = (total_variants + n_chunks - 1) // n_chunks
        else:
            chunk_size = max(1, chunk_size)
            lines_per_chunk = chunk_size
            n_chunks = (total_variants + chunk_size - 1) // chunk_size
        log_info(f"交集位点总数：{total_variants}；计划分块数：{n_chunks}；每块行数：{lines_per_chunk}")

        # split 命令分割文件
        split_prefix = _p("chunk_")
        split_cmd = [
            "split",
            "-d",
            "-l", str(lines_per_chunk),
            "--additional-suffix=.ids",
            intersect_ids,
            split_prefix
        ]
        _run_cmd(split_cmd)

        # 获取所有 chunk 文件名，按编号排序
        chunk_files = []
        for i in range(n_chunks):
            chunk_id = f"{i:02d}"
            path = f"{split_prefix}{chunk_id}.ids"
            if not os.path.exists(path):
                # 可能最后一个chunk少于n_chunks，跳过不存在文件
                continue
            chunk_files.append(path)

        if len(chunk_files) == 0:
            raise RuntimeError("No chunk files generated from intersect.ids.")

        # ---------- 5) 并行运行各块命令 ----------
        plink_cmds = []
        bcftools_cmds = []
        chunk_infos = []

        for idx, chunk_path in enumerate(chunk_files):
            chunk_id = f"{idx:02d}"
            post_gt_tsv = _p(f"post.gt.chunk_{chunk_id}.tsv")
            pre_gt_tsv = _p(f"pre.gt.chunk_{chunk_id}.tsv")

            # plink2 命令 (A-transpose (variant-major 0/1/2); dosage counts controlled by --export-allele (ALT from BIM ID))
            cmd_plink = [
                plink_bin,
                "--bfile", plink_prefix,
                "--extract", chunk_path,
                "--export", "A-transpose",
                "--export-allele", allele_map_path,
                "--threads", str(plink_threads),
                "--out", post_gt_tsv
            ]
            if chrom is not None:
                cmd_plink.extend(["--chr", chrom])
            plink_cmds.append((cmd_plink, post_gt_tsv))

            # 不再使用 -r 以免染色体命名差异（22 vs chr22）导致整块被过滤掉；
            # 只用 -i 'ID=@file' 通过交集ID做精确筛选。
            cmd_bcftools = [
                "bash", "-lc",
                f"{bcftools_bin} view --threads {bcftools_threads} -i 'ID=@{chunk_path}' "
                f"-S {sample_order_path} -Ou {vcf_path} | {bcftools_bin} query -f '%ID[\\t%GT]\\n' > {pre_gt_tsv}"
            ]
            bcftools_cmds.append(cmd_bcftools)

        log_info(f"开始并行运行 PLINK（每块导出 ALT 剂量 0/1/2）...")
        plink_cmds_only = [c[0] for c in plink_cmds]
        try:
            _run_pool(plink_cmds_only, max_parallel=max_parallel)
        except Exception as e:
            raise RuntimeError(f"PLINK chunk processing failed: {e}")

        # 清理 .traw → .tsv（两阶段：先读表头决定列，再按需读列；节省内存）
        import numpy as np
        for _, out_tsv in plink_cmds:
            traw_file = out_tsv + ".traw"
            if not os.path.exists(traw_file):
                raise RuntimeError(f"Expected .traw file not found: {traw_file}")
            # 先只读表头，确定样本列
            hdr = pd.read_csv(traw_file, sep="\t", nrows=0)
            meta_cols = ["CHR","SNP","(C)M","POS","COUNTED","ALT"]
            sample_cols = [c for c in hdr.columns if c not in meta_cols]
            usecols = ["SNP"] + sample_cols
            # 只读需要的列
            df = pd.read_csv(traw_file, sep="\t", usecols=usecols, dtype=str)
            df = df.rename(columns={"SNP":"ID"})
            # 新的 IID 解析逻辑，兼容 FID_IID、IID含下划线等情况
            iid_cols = []
            heuristic_hits = 0
            for c in sample_cols:
                if c in fid_iid_to_iid:
                    iid_cols.append(fid_iid_to_iid[c])
                else:
                    # 兼容异常：若表头直接就是 IID
                    if c in iid_set:
                        iid_cols.append(c)
                    else:
                        # 尝试用第一个下划线分割后的部分，如果在 IID 集合中，则采用
                        if "_" in c:
                            after_first = c.split("_", 1)[1]
                            if after_first in iid_set:
                                iid_cols.append(after_first)
                                heuristic_hits += 1
                                continue
                            # 最后兜底：保留原逻辑取最后一段，但记录为启发式
                            last_seg = c.rsplit("_", 1)[-1]
                            iid_cols.append(last_seg)
                            heuristic_hits += 1
                        else:
                            # 无下划线，无法分割，保留原样
                            iid_cols.append(c)
            if heuristic_hits > 0:
                log_warn(f"样本列名有 {heuristic_hits} 项通过启发式规则推断 IID（请确认 FID/IID 映射是否完整）")
            df.columns = ["ID"] + iid_cols
            # 将 0/1/2/NA 规范成 可空整数；写出缺失为 '.'
            for c in iid_cols:
                s = pd.to_numeric(df[c], errors="coerce")
                try:
                    df[c] = s.astype("Int8")
                except TypeError:
                    df[c] = s.astype("Int64")
            df.to_csv(out_tsv, sep="\t", index=False, na_rep=".")
            os.remove(traw_file)
        log_info("PLINK 各块 .traw 已清洗为 .tsv（仅 ID+IID；整数 0/1/2；缺失 '.'），并删除 .traw")

        log_info(f"开始并行运行 bcftools（每块导出 GT→0/1/2/'.'）...")
        try:
            _run_pool(bcftools_cmds, max_parallel=max_parallel)
        except Exception as e:
            raise RuntimeError(f"bcftools chunk processing failed: {e}")

        # 诊断：如出现空 pre.gt.chunk_*.tsv，提示用户检查 ID/样本名
        empty_chunks = []
        for idx, (_, post_gt_tsv) in enumerate(plink_cmds):
            pre_gt_tsv = _p(f"pre.gt.chunk_{idx:02d}.tsv")
            try:
                sz = os.path.getsize(pre_gt_tsv)
                if sz == 0:
                    empty_chunks.append(pre_gt_tsv)
            except FileNotFoundError:
                empty_chunks.append(pre_gt_tsv)
        if empty_chunks:
            log_warn("检测到部分 VCF 分块输出为空：")
            for pth in empty_chunks[:5]:
                log_warn(f"  空文件示例：{pth}")
            log_warn("建议排查：1) chunk.ids 是否存在于 VCF；2) sample.order 是否与 VCF HEADER 样本名一致；3) VCF 是否存在 '.' 空ID；4) 如使用了 chrom 过滤，尝试去掉 -r 再取交集。")

        # 将 VCF 侧 pre.gt.chunk_*.tsv 的 GT 字符串直接映射为 0/1/2/'.'，并补上表头(ID+IID)
        with open(sample_order_path, "r") as f:
            sample_ids = [x.strip() for x in f if x.strip()]
        for idx, (_, _) in enumerate(plink_cmds):
            chunk_id = f"{idx:02d}"
            pre_path = _p(f"pre.gt.chunk_{chunk_id}.tsv")
            if not os.path.exists(pre_path):
                continue
            tmp_out = pre_path + ".tmp"
            with open(pre_path, "r") as fin, open(tmp_out, "w") as fout:
                # 写表头：ID + IID（来自 sample.order）
                fout.write("ID\t" + "\t".join(sample_ids) + "\n")
                for line in fin:
                    parts = line.rstrip("\n").split("\t")
                    if not parts or len(parts) < 2:
                        continue
                    vid = parts[0]
                    gts = parts[1:]
                    mapped = [_vcf_gt_to012(gt) for gt in gts]
                    # 将整数映射为字符串，其它（'.'）保持
                    mapped_str = [str(x) for x in mapped]
                    fout.write(vid + "\t" + "\t".join(mapped_str) + "\n")
            os.replace(tmp_out, pre_path)
        log_info("VCF 各块 GT 已映射为 0/1/2/'.' 并补上表头")

        # ---------- 6) 生成 manifest 文件 ----------
        manifest_path = _p("aligned_geno_manifest.tsv")
        with open(manifest_path, "w") as mf:
            mf.write("chunk_id\tnum_variants\tpre_tsv\tpost_tsv\n")
            for idx, (_, post_gt_tsv) in enumerate(plink_cmds):
                chunk_id = f"{idx:02d}"
                # 统计变体数（行数减1）
                try:
                    out = subprocess.check_output(["wc", "-l", post_gt_tsv], text=True)
                    n_lines = int(out.strip().split()[0])
                    n_vars = n_lines - 1 if n_lines > 0 else 0
                except Exception:
                    n_vars = 0
                pre_gt_tsv = _p(f"pre.gt.chunk_{chunk_id}.tsv")
                mf.write(f"{chunk_id}\t{n_vars}\t{pre_gt_tsv}\t{post_gt_tsv}\n")

        # ---------- 7) 合并所有分块至工作目录，并返回路径 ----------
        # 以 manifest 顺序合并，首块保留表头，其后跳过表头
        post_files = []
        pre_files = []
        with open(manifest_path, "r") as mf:
            next(mf)
            for line in mf:
                chunk_id, nvars, pre_tsv, post_tsv = line.strip().split('\t')
                pre_files.append(pre_tsv)
                post_files.append(post_tsv)
        for fpath in pre_files + post_files:
            if not os.path.exists(fpath):
                log_warn(f"missing chunk file: {fpath}")
            else:
                try:
                    if os.path.getsize(fpath) == 0:
                        log_warn(f"empty chunk file: {fpath}")
                except Exception:
                    pass

        def concat_tsv(tsv_files: List[str], out_path: str):
            with open(out_path, "w") as fout:
                for i, fpath in enumerate(tsv_files):
                    with open(fpath, "r") as fin:
                        for j, line in enumerate(fin):
                            if i > 0 and j == 0:
                                continue  # skip header except first file
                            fout.write(line)

        chrom_tag = f"chr{chrom}" if chrom is not None else "chrALL"
        pre_merged = os.path.join(base_tmp, f"{chrom_tag}.pre.gt.tsv")
        post_merged = os.path.join(base_tmp, f"{chrom_tag}.post.gt.tsv")

        log_info("开始合并所有分块 TSV...")
        concat_tsv(pre_files, pre_merged)
        concat_tsv(post_files, post_merged)
        log_info(f"合并完成：pre → {pre_merged}；post → {post_merged}")
        # 统计合并后的行列数（仅读表头 + wc）
        try:
            out = subprocess.check_output(["wc", "-l", pre_merged], text=True)
            pre_lines = int(out.strip().split()[0])
            out = subprocess.check_output(["wc", "-l", post_merged], text=True)
            post_lines = int(out.strip().split()[0])
            # 读一个文件表头统计样本数
            with open(post_merged, "r") as f:
                header_cols = f.readline().rstrip("\n").split("\t")
            nsamples = max(0, len(header_cols) - 1)
            log_info(f"合并后行数：pre={max(0, pre_lines-1)}，post={max(0, post_lines-1)}；样本数={nsamples}")
        except Exception as e:
            log_warn(f"统计合并规模失败：{e}")
        try:
            with open(LOG_FILE, "a") as lf:
                lf.write(f"[{_now()}][INFO] 完成 build_aligned_geno_mats\n")
                lf.write("="*72 + "\n")
        except Exception:
            pass
        # 函数最终返回两个文件路径（pre, post）
        return pre_merged, post_merged

    finally:
        if not keep_temp:
            try:
                if os.path.isdir(cache_dir):
                    shutil.rmtree(cache_dir)
            except Exception as e:
                log_warn(f"清理临时目录失败: {e}")
                


def harmonize_gt_matrices(
    pre_mt: str,
    post_mt: str,
    out_pre: Optional[str] = None,
    out_post: Optional[str] = None,
    max_check_rows: int = 200000,
) -> Tuple[str, str]:
    """
    对齐/修复两个大型基因型矩阵（pre / post）：
      1) 检查并确保 ID 列（第一列）为 CHROM:POS:REF:ALT 格式；
      2) 检查 ID 是否按 POS 升序；若不是则按 POS 进行稳定排序（保留表头在首行）；
      3) 检查两侧列名（样本）内容/数量/顺序是否一致；若不一致，按照 post 的顺序修复：
         - 若集合一致，仅重排 pre 的列顺序；
         - 若集合不同，取两侧交集并警告，分别重写 pre/post 为交集顺序；
      4) 过程产生的日志会写入 `<work_dir>/geno_miss_bias.log`（保持现有日志格式标签）。

    为避免内存溢出，所有检查与修复均采用**流式处理**：不将全表读入内存。

    参数：
      pre_mt   : VCF → 矩阵路径（形如 `chr*.pre.gt.tsv`），第一列为 `ID`，其余列为 IID；
      post_mt  : PLINK → 矩阵路径（形如 `chr*.post.gt.tsv`），第一列为 `ID`，其余列为 IID；
      out_pre  : 修复后的 pre 输出路径（默认=原路径加后缀 `.harm.tsv`）；
      out_post : 修复后的 post 输出路径（默认=原路径加后缀 `.harm.tsv`）；
      max_check_rows : 验证排序时最多检查的行数（全量检查仍为流式；该参数仅控制提前退出阈值）。

    返回：
      (out_pre_path, out_post_path)
    """
    assert os.path.exists(pre_mt), f"pre_mt not found: {pre_mt}"
    assert os.path.exists(post_mt), f"post_mt not found: {post_mt}"

    out_pre = out_pre or (pre_mt + ".harm.tsv")
    out_post = out_post or (post_mt + ".harm.tsv")

    # 初始化日志文件（若未由其他函数设定）。日志路径默认放在输出目录公共前缀下。
    global LOG_FILE
    if LOG_FILE is None:
        pre_dir = os.path.dirname(os.path.abspath(out_pre)) or os.getcwd()
        post_dir = os.path.dirname(os.path.abspath(out_post)) or os.getcwd()
        try:
            base_dir = os.path.commonpath([pre_dir, post_dir])
        except Exception:
            base_dir = pre_dir or post_dir or os.getcwd()
        LOG_FILE = os.path.join(base_dir, "geno_miss_bias.harmonize.log")
        try:
            with open(LOG_FILE, "a") as lf:
                lf.write("="*72 + "\n")
                lf.write(f"[{_now()}][INFO] 开始执行 harmonize_gt_matrices\n")
                lf.write(f"[{_now()}][INFO] pre={pre_mt}; post={post_mt}\n")
                lf.write(f"[{_now()}][INFO] out_pre={out_pre}; out_post={out_post}\n")
        except Exception:
            pass
    log_info(f"日志文件：{LOG_FILE}")

    # --- Logging: banner and context ---
    log_info(f"开始对齐/修复矩阵：pre={pre_mt}；post={post_mt}")
    log_info(f"输出路径：out_pre={out_pre}；out_post={out_post}")
    log_info(f"排序抽检最大行数（max_check_rows）={max_check_rows}")

    # --- 小工具 ---
    import re
    id_pat = re.compile(r"^[^:\t]+:\d+:[ACGTN]+:[ACGTN]+$")

    def _read_header(path: str) -> List[str]:
        with open(path, "r") as f:
            hdr = f.readline().rstrip("\n").split("\t")
        if len(hdr) < 2 or hdr[0] != "ID":
            log_err(f"{os.path.basename(path)} 表头首列不是 'ID'，实际为：{hdr[0] if hdr else '<空>'}")
            raise ValueError(f"Bad header in {path}: first column must be 'ID'")
        return hdr

    def _check_id_format(path: str, sample_n: int = 10000):
        bad = 0
        with open(path, "r") as f:
            next(f)  # skip header
            for i, line in enumerate(f, start=1):
                if not line.strip():
                    continue
                vid = line.split("\t", 1)[0]
                if not id_pat.match(vid):
                    bad += 1
                    if bad <= 5:
                        log_warn(f"{os.path.basename(path)} 非法ID示例：{vid}")
                if i >= sample_n and bad == 0:
                    break
        if bad > 0:
            log_warn(f"{os.path.basename(path)} 中检测到 {bad} 个不符合 CHROM:POS:REF:ALT 的 ID（仅抽样统计）。")

    def _is_sorted_by_pos(path: str, max_check: int = max_check_rows) -> bool:
        last_pos = -1
        checked = 0
        with open(path, "r") as f:
            next(f)  # header
            for line in f:
                if not line:
                    continue
                parts = line.split("\t", 1)
                if not parts:
                    continue
                vid = parts[0]
                a = vid.split(":")
                if len(a) < 2 or not a[1].isdigit():
                    return False
                pos = int(a[1])
                if pos < last_pos:
                    return False
                last_pos = pos
                checked += 1
                if checked >= max_check:
                    break
        return True

    def _sort_by_pos(path: str, out_path: str):
        tmp_body_pref = out_path + ".bodypref.tmp"
        tmp_body_sorted = out_path + ".bodysorted.tmp"
        tmp_out = out_path + ".tmp"
        # 1) Python 流式：为每行前置 POS（无法解析则置 999999999），避免 awk 在超长行/引号上的不确定性
        total = 0
        with open(path, "r") as fin, open(tmp_body_pref, "w") as fpref:
            header_line = fin.readline()  # 仅读取一次表头
            for line in fin:
                if not line:
                    continue
                s = line.rstrip("\r\n")
                if not s:
                    continue
                # 提取首列 ID 的 POS
                try:
                    first_tab = s.find("\t")
                    idv = s if first_tab == -1 else s[:first_tab]
                    pos_str = idv.split(":", 2)[1]
                    pos = int(pos_str)
                except Exception:
                    pos = 999_999_999
                fpref.write(f"{pos}\t{s}\n")
                total += 1
        log_info(f"[POS-SORT] 需要排序的主体行数（不含表头）={total}")
        if total == 0:
            # 主体为空：直接复制并返回
            log_warn("[POS-SORT] 主体行为 0，排序无意义：保留原文件")
            if path != out_path:
                shutil.copyfile(path, out_path)
            return

        # 2) 用 GNU sort 对前置 POS 做数值排序
        qpref = shlex.quote(tmp_body_pref)
        qsorted = shlex.quote(tmp_body_sorted)
        _run_cmd(f"bash -lc \"LC_ALL=C sort -t $'\\t' -k1,1n {qpref} > {qsorted}\"")

        # 3) 合并表头 + 排序后的主体，并去掉临时 POS 列
        with open(tmp_out, "w") as fout, open(path, "r") as fin_h, open(tmp_body_sorted, "r") as fsorted:
            # 输出表头原样
            fout.write(fin_h.readline())
            # 输出主体：去掉前缀 POS 列
            out_rows = 0
            for row in fsorted:
                tab = row.find("\t")
                if tab != -1:
                    fout.write(row[tab+1:])
                else:
                    # 理论兜底：保留整行
                    fout.write(row)
                out_rows += 1
        log_info(f"[POS-SORT] 排序后主体行数（不含表头）={out_rows}")

        # 4) 行数一致性校验：排序前（含表头） vs 排序后（含表头）
        try:
            c1 = int(subprocess.check_output(['bash','-lc', f"wc -l < {shlex.quote(path)}"], text=True).strip())
            c2 = int(subprocess.check_output(['bash','-lc', f"wc -l < {shlex.quote(tmp_out)}"], text=True).strip())
            if c1 != c2:
                log_warn(f"POS 排序前后行数不一致：原始={c1}，排序后={c2}。保留未排序版本作为输出；排序结果保留：{tmp_out}")
                return
        except Exception as e:
            log_warn(f"排序后行数校验失败：{e}")

        os.replace(tmp_out, out_path)
        # 清理临时文件
        for _p in (tmp_body_pref, tmp_body_sorted):
            try:
                os.remove(_p)
            except Exception:
                pass
        log_info(f"已按 POS 升序排序：{out_path}")

    def _reorder_stream(in_path: str, out_path: str, desired_cols: List[str]):
        """按 desired_cols（IID 顺序）重排/裁剪列：流式处理，仅持有当前行。"""
        with open(in_path, "r") as fin:
            header = fin.readline().rstrip("\n").split("\t")
            cur_cols = header[1:]
            idx_map = {name: i for i, name in enumerate(cur_cols)}
            missing = [c for c in desired_cols if c not in idx_map]
            extra = [c for c in cur_cols if c not in set(desired_cols)]
            if missing:
                log_warn(f"列缺失（将被跳过）：{len(missing)} 个，例如：{missing[:5]}")
            if extra:
                log_warn(f"多余列（将被丢弃）：{len(extra)} 个，例如：{extra[:5]}")
            # 仅输出交集，顺序按 desired_cols
            keep_idx = [idx_map[c] for c in desired_cols if c in idx_map]
            with open(out_path, "w") as fout:
                fout.write("ID\t" + "\t".join([cur_cols[i] for i in keep_idx]) + "\n")
                for line in fin:
                    if not line:
                        continue
                    parts = line.rstrip("\n").split("\t")
                    idv = parts[0]
                    row = [parts[1 + i] if 1 + i < len(parts) else '.' for i in keep_idx]
                    fout.write(idv + "\t" + "\t".join(row) + "\n")

    # 读取表头
    post_hdr = _read_header(post_mt)
    pre_hdr = _read_header(pre_mt)
    post_samples = post_hdr[1:]
    pre_samples = pre_hdr[1:]
    log_info(f"表头校验完成：pre 列数（含ID）={len(pre_hdr)}，post 列数（含ID）={len(post_hdr)}")
    log_info(f"样本数：pre={len(pre_samples)}；post={len(post_samples)}")

    # 1) ID 格式检查（抽样）
    _check_id_format(pre_mt)
    _check_id_format(post_mt)
    log_info("ID 格式抽检完成（若存在异常会在上方 WARN 提示示例）")

    # 2) POS 升序检查；不满足则稳定排序
    # 为避免覆写原文件，这里生成排序后的临时文件，并作为后续输入
    log_info("开始检查 POS 升序……")
    pre_sorted = pre_mt if _is_sorted_by_pos(pre_mt) else (pre_mt + ".possorted.tsv")
    post_sorted = post_mt if _is_sorted_by_pos(post_mt) else (post_mt + ".possorted.tsv")
    log_info(f"POS 升序检查结果：pre_sorted={'是' if pre_sorted==pre_mt else '否（将排序）'}；post_sorted={'是' if post_sorted==post_mt else '否（将排序）'}")
    pre_need_sort = (pre_sorted != pre_mt)
    post_need_sort = (post_sorted != post_mt)
    if pre_need_sort:
        log_warn("检测到 pre_mt 非 POS 升序，将进行排序……")
        _sort_by_pos(pre_mt, pre_sorted)
    if post_need_sort:
        log_warn("检测到 post_mt 非 POS 升序，将进行排序……")
        _sort_by_pos(post_mt, post_sorted)

    # 3) 列名一致性与顺序修复
    log_info("开始检查样本列名集合与顺序……")
    same_set = set(pre_samples) == set(post_samples)
    same_order = pre_samples == post_samples
    log_info(f"样本集合是否一致：{same_set}；样本顺序是否一致：{same_order}")

    if same_set and same_order:
        if not pre_need_sort and not post_need_sort:
            # 完全无需修复：直接返回原路径（不做任何复制）
            log_info("无需修复：直接返回原路径（ID 已按 POS 升序，样本集合与顺序一致）")
            log_info(f"对齐完成：out_pre={pre_mt}；out_post={post_mt}")
            try:
                with open(LOG_FILE, "a") as lf:
                    lf.write(f"[{_now()}][INFO] 完成 harmonize_gt_matrices\n")
                    lf.write("="*72 + "\n")
            except Exception:
                pass
            return pre_mt, post_mt
        else:
            # 仅涉及排序：复制排序后的文件到输出
            log_info("无需更改样本集合/顺序：仅输出按 POS 排序后的文件")
            if pre_sorted != out_pre:
                shutil.copyfile(pre_sorted, out_pre)
            else:
                out_pre = pre_sorted
            if post_sorted != out_post:
                shutil.copyfile(post_sorted, out_post)
            else:
                out_post = post_sorted
            log_info(f"对齐完成：out_pre={out_pre}；out_post={out_post}")
            try:
                with open(LOG_FILE, "a") as lf:
                    lf.write(f"[{_now()}][INFO] 完成 harmonize_gt_matrices\n")
                    lf.write("="*72 + "\n")
            except Exception:
                pass
            return out_pre, out_post

    if same_set and not same_order:
        log_info("仅重排 pre 的样本顺序以匹配 post（列集合一致）")
        # 仅重排 pre 为 post 顺序；post 直接拷贝
        _reorder_stream(pre_sorted, out_pre, post_samples)
        if post_sorted != out_post:
            shutil.copyfile(post_sorted, out_post)
        else:
            out_post = post_sorted
        log_info(f"对齐完成：out_pre={out_pre}；out_post={out_post}")
        try:
            with open(LOG_FILE, "a") as lf:
                lf.write(f"[{_now()}][INFO] 完成 harmonize_gt_matrices\n")
                lf.write("="*72 + "\n")
        except Exception:
            pass
        return out_pre, out_post

    # 集合不同：取交集并分别重写两侧
    inter = [s for s in post_samples if s in set(pre_samples)]
    if not inter:
        log_err("两侧样本集合没有交集，无法对齐矩阵。")
        raise ValueError("No overlapping samples between pre_mt and post_mt")
    dropped_post = [s for s in post_samples if s not in set(inter)]
    dropped_pre = [s for s in pre_samples if s not in set(inter)]
    if dropped_post:
        log_warn(f"post_mt 中将被丢弃的样本：{len(dropped_post)} 个，例如：{dropped_post[:5]}")
    if dropped_pre:
        log_warn(f"pre_mt 中将被丢弃的样本：{len(dropped_pre)} 个，例如：{dropped_pre[:5]}")

    log_info(f"按交集重写两侧：交集样本数={len(inter)}；pre丢弃={len(dropped_pre)}；post丢弃={len(dropped_post)}")
    _reorder_stream(pre_sorted, out_pre, inter)
    _reorder_stream(post_sorted, out_post, inter)
    log_info(f"对齐完成：out_pre={out_pre}；out_post={out_post}")
    try:
        with open(LOG_FILE, "a") as lf:
            lf.write(f"[{_now()}][INFO] 完成 harmonize_gt_matrices\n")
            lf.write("="*72 + "\n")
    except Exception:
        pass
    return out_pre, out_post


import sqlite3

def reorder_pre_to_post(
    pre_mt: str,
    post_mt: str,
    out_pre: Optional[str] = None,
    out_post: Optional[str] = None,
) -> Tuple[str, str]:
    """
    使用 `post_mt` 的**行顺序与列顺序**，对超大的 `pre_mt` 进行严格对齐重排，输出两个新文件路径：
      - `out_pre`：按 `post_mt` 的列（样本 IID）与行（变体 ID）顺序重排后的 pre 矩阵；
      - `out_post`：`post_mt` 的副本（作为配对输出，便于下游统一使用）。

    设计目标：
      - **不把整表读入内存**。通过磁盘 SQLite 索引（ID -> 文件偏移）实现 O(1) 随机读取 pre 的行；
      - 列重排按 post 的 IID 顺序进行；pre 中缺失的样本列用 '.' 补齐，pre 中多余列被丢弃并记录日志；
      - 若某个 post 的 ID 在 pre 中不存在，整行以 '.' 填充并 WARN 计数；
      - 表头首列要求为 'ID'，分隔符为 TAB。

    返回：(out_pre_path, out_post_path)
    """
    assert os.path.exists(pre_mt), f"pre_mt not found: {pre_mt}"
    assert os.path.exists(post_mt), f"post_mt not found: {post_mt}"

    out_pre = out_pre or (pre_mt + ".reordered.tsv")
    out_post = out_post or (post_mt + ".reordered.tsv")

    # 准备日志（如未设置全局 LOG_FILE，则在 CWD 建一个专用日志）
    global LOG_FILE
    if LOG_FILE is None:
        try:
            base_dir = os.getcwd()
        except Exception:
            base_dir = "."
        LOG_FILE = os.path.join(base_dir, "geno_miss_bias.reorder.log")
        try:
            with open(LOG_FILE, "a") as lf:
                lf.write("="*72 + "\n")
                lf.write(f"[{_now()}][INFO] 开始执行 reorder_pre_to_post\n")
                lf.write(f"[{_now()}][INFO] pre={pre_mt}; post={post_mt}\n")
                lf.write(f"[{_now()}][INFO] out_pre={out_pre}; out_post={out_post}\n")
        except Exception:
            pass
    log_info(f"日志文件：{LOG_FILE}")

    # 读取表头，建立列顺序映射
    with open(post_mt, "r") as fpost:
        post_header = fpost.readline().rstrip("\n").split("\t")
    if len(post_header) < 2 or post_header[0] != "ID":
        log_err("post_mt 表头首列必须为 'ID'")
        raise ValueError("Bad post header")
    post_iids = post_header[1:]

    with open(pre_mt, "r") as fpre:
        pre_header = fpre.readline().rstrip("\n").split("\t")
    if len(pre_header) < 2 or pre_header[0] != "ID":
        log_err("pre_mt 表头首列必须为 'ID'")
        raise ValueError("Bad pre header")
    pre_iids = pre_header[1:]

    # 列重排映射：按 post 的顺序提取 pre 的列；缺失用 '.'
    pre_idx_map = {name: i for i, name in enumerate(pre_iids)}
    keep_idx = []
    missing_cols = []
    for name in post_iids:
        if name in pre_idx_map:
            keep_idx.append(pre_idx_map[name])
        else:
            keep_idx.append(None)
            missing_cols.append(name)
    extra_cols = [c for c in pre_iids if c not in set(post_iids)]
    if missing_cols:
        log_warn(f"pre 中缺失 {len(missing_cols)} 个样本列，将以 '.' 补齐，例如：{missing_cols[:5]}")
    if extra_cols:
        log_warn(f"pre 中存在 {len(extra_cols)} 个多余样本列，将被丢弃，例如：{extra_cols[:5]}")

    # 为 post 建立行顺序（按 post 行序输出）
    # 同时为 pre 建立磁盘索引：ID -> 文件偏移
    tmp_db = None
    conn = None
    try:
        tmp_db = tempfile.NamedTemporaryFile(prefix="pre_index_", suffix=".sqlite", delete=False)
        tmp_db.close()
        conn = sqlite3.connect(tmp_db.name)
        conn.execute("PRAGMA synchronous=OFF")
        conn.execute("PRAGMA journal_mode=OFF")
        conn.execute("CREATE TABLE idx (id TEXT PRIMARY KEY, off INTEGER)")
        cur = conn.cursor()

        # 扫描 pre 文件，建立偏移索引
        log_info("开始建立 pre 行索引（ID -> 文件偏移）……")
        with open(pre_mt, "rb") as fpreb:
            header_line = fpreb.readline()  # skip header
            off = fpreb.tell()
            batch = []
            batch_size = 100000
            n_rows = 0
            while True:
                pos = fpreb.tell()
                line = fpreb.readline()
                if not line:
                    break
                if line == b"\n":
                    off = fpreb.tell()
                    continue
                # 仅提取第一列 ID（直到第一个 tab）
                try:
                    tab = line.find(b"\t")
                    if tab == -1:
                        off = fpreb.tell()
                        continue
                    vid = line[:tab].decode("utf-8", errors="ignore")
                except Exception:
                    off = fpreb.tell()
                    continue
                batch.append((vid, pos))
                n_rows += 1
                if len(batch) >= batch_size:
                    cur.executemany("INSERT OR REPLACE INTO idx(id, off) VALUES(?, ?)", batch)
                    conn.commit()
                    batch.clear()
                off = fpreb.tell()
            if batch:
                cur.executemany("INSERT OR REPLACE INTO idx(id, off) VALUES(?, ?)", batch)
                conn.commit()
        log_info("pre 行索引建立完成")

        # 写出对齐后的 pre，并生成 post 的副本
        log_info("开始根据 post 行/列顺序输出对齐后的 pre……")
        with open(pre_mt, "rb") as fpreb, \
             open(post_mt, "r") as fpost, \
             open(out_pre, "w") as fpre_out:
            # 写表头：ID + post_iids
            fpre_out.write("ID\t" + "\t".join(post_iids) + "\n")
            # 跳过 pre 头
            _ = fpreb.readline()
            # 遍历 post 的每一行，按顺序输出对应的 pre 行
            missing_ids = 0
            total_ids = 0
            for j, line in enumerate(fpost):
                if j == 0:
                    continue  # 已读过表头
                parts = line.rstrip("\n").split("\t", 1)
                if not parts or parts[0] == "":
                    continue
                vid = parts[0]
                total_ids += 1
                row = None
                # 在索引中查找偏移
                r = conn.execute("SELECT off FROM idx WHERE id=?", (vid,)).fetchone()
                if r is not None:
                    off = r[0]
                    fpreb.seek(off)
                    raw = fpreb.readline().decode("utf-8", errors="ignore").rstrip("\n")
                    cols = raw.split("\t")
                    # cols[0] 应该是 ID
                    if cols and cols[0] == vid:
                        # 重排列
                        vals = cols[1:]
                        out_vals = []
                        for k in keep_idx:
                            if k is None:
                                out_vals.append('.')
                            else:
                                v = vals[k] if k < len(vals) else ''
                                out_vals.append(v if v != '' else '.')
                        row = vid + "\t" + "\t".join(out_vals)
                if row is None:
                    # pre 中不存在该 ID：整行 '.'
                    missing_ids += 1
                    row = vid + "\t" + "\t".join(['.'] * len(post_iids))
                fpre_out.write(row + "\n")
        if missing_ids > 0:
            log_warn(f"在 post 的 {total_ids} 个 ID 中，pre 缺失 {missing_ids} 个；对应行已用 '.' 填充")
        else:
            log_info("pre 与 post 的变体 ID 全部匹配，无缺失")

        # 复制 post 为配对输出文件
        if out_post == post_mt:
            log_info("目标 out_post 等于原始 post 路径，跳过复制")
        else:
            shutil.copyfile(post_mt, out_post)
            log_info(f"已复制 post 到：{out_post}")

        try:
            with open(LOG_FILE, "a") as lf:
                lf.write(f"[{_now()}][INFO] 完成 reorder_pre_to_post\n")
                lf.write("="*72 + "\n")
        except Exception:
            pass
        return out_pre, out_post

    finally:
        try:
            if conn is not None:
                conn.close()
        except Exception:
            pass
        if tmp_db and os.path.exists(tmp_db.name):
            try:
                os.remove(tmp_db.name)
            except Exception:
                pass



from typing import List

def build_sample_group_info(
    info_xls: str,
    sample_list: List[str],   # 改成 List
    id_col: str = "ID",
    use_cols: Optional[List[str]] = None,
) -> pd.DataFrame:
    """
    读取样本信息表，并按给定样本列表筛选与清洗列：
      - 仅保留 `id_col` 在 `sample_list` 内的行（`id_col` 按字符串读取）。
      - 默认选择列：['WGS', 'Target DP (JHRPv4)', 'DP (JHRPv4)']。
      - 对含有多项且以 '|' 分割的单元格：
          * `Target DP (JHRPv4)`: 若为 "15x | 30x" 等，取第一个非空项并去空格（例如 "15x"）。
          * `WGS`: 取第一个非空项并去空格（例如 "HiSeqX 15x | NovaSeq 30x" → "HiSeqX 15x"）。
          * `DP (JHRPv4)`: 取第一个非空项并去空格（例如 "16.5307 | 16.5307" → "16.5307"），并转换为 float。
    返回：包含 `id_col` 与所选列的 DataFrame（`id_col` 为列，不设为索引）。
    """

    # 初始化独立日志文件（仅当未由其它流程设定 LOG_FILE 时）
    global LOG_FILE
    if LOG_FILE is None:
        # 默认将日志写到当前工作目录
        try:
            base_dir = os.getcwd()
        except Exception:
            base_dir = "."
        LOG_FILE = os.path.join(base_dir, "group_info.log")
        try:
            with open(LOG_FILE, "a") as lf:
                lf.write("="*72 + "\n")
                lf.write(f"[{_now()}][INFO] 开始执行 build_sample_group_info\n")
                lf.write(f"[{_now()}][INFO] info_xls={info_xls}\n")
                lf.write(f"[{_now()}][INFO] cwd={os.path.abspath(base_dir)}\n")
        except Exception:
            pass
    log_info(f"日志文件：{LOG_FILE}")

    assert os.path.exists(info_xls), f"info_xls not found: {info_xls}"

    default_cols = ["WGS", "Target DP (JHRPv4)", "DP (JHRPv4)"]
    if use_cols is None:
        use_cols = default_cols
    else:
        use_cols = list(use_cols)

    read_cols = [id_col] + use_cols
    try:
        df = pd.read_excel(info_xls, dtype={id_col: str}, engine=None, usecols=lambda c: c in set(read_cols))
    except Exception as e:
        log_warn(f"read_excel 失败（可能不是 Excel？）：{e}，尝试按 TSV 读取……")
        try:
            df = pd.read_csv(info_xls, sep="\t", dtype={id_col: str}, usecols=lambda c: c in set(read_cols))
        except Exception as e2:
            log_err(f"既非可读 Excel 也非 TSV：{e2}")
            raise

    df = df.loc[df[id_col].notna()].copy()
    df[id_col] = df[id_col].astype(str)

    # 仅保留在 sample_list 内的样本
    before = len(df)
    df = df[df[id_col].isin(sample_list)].copy()
    after = len(df)
    log_info(f"样本筛选：输入 {before} 行，匹配 sample_list 后保留 {after} 行")

    for col in use_cols:
        if col not in df.columns:
            df[col] = pd.NA
            log_warn(f"缺少列：{col}，将以缺失填充")

    def _first_token(val: Any) -> Any:
        if pd.isna(val):
            return val
        s = str(val)
        parts = [p.strip() for p in s.split('|')]
        for p in parts:
            if p != "":
                return p
        return ""

    if "Target DP (JHRPv4)" in df.columns:
        df["Target DP (JHRPv4)"] = df["Target DP (JHRPv4)"].map(_first_token)
        df["Target DP (JHRPv4)"] = df["Target DP (JHRPv4)"].astype(str).str.strip()

    if "WGS" in df.columns:
        df["WGS"] = df["WGS"].map(_first_token)

    if "DP (JHRPv4)" in df.columns:
        df["DP (JHRPv4)"] = df["DP (JHRPv4)"].map(_first_token)
        def _to_float(x):
            try:
                return float(str(x).strip()) if pd.notna(x) and str(x).strip() != "" else pd.NA
            except Exception:
                return pd.NA
        df["DP (JHRPv4)"] = df["DP (JHRPv4)"].map(_to_float)

    df = df[[id_col] + use_cols]

    log_info(f"已生成分组信息表：行数={len(df)}；列={ [id_col] + use_cols }")
    try:
        with pd.option_context('display.max_columns', None, 'display.width', 200):
            log_info("示例预览（前5行）：\n" + df.head(5).to_string(index=False))
    except Exception:
        pass

    try:
        with open(LOG_FILE, "a") as lf:
            lf.write(f"[{_now()}][INFO] 完成 build_sample_group_info\n")
            lf.write("="*72 + "\n")
    except Exception:
        pass
    return df


def compute_coverage_transition_counts(
    out_pre_order: str,
    out_post_order: str,
    df_info: pd.DataFrame,
    id_col: str = "ID",
    coverage_col: str = "Target DP (JHRPv4)",
    coverage_labels: Tuple[str, str] = ("15x", "30x"),
    out_path: Optional[str] = None,
    chrom: Optional[str] = None,
) -> str:
    """
    基于**已对齐行列顺序**的超大矩阵（`out_pre_order`, `out_post_order`），在
    **15x/30x 覆盖分组**内逐位点统计以下转换计数：
      - 0 -> .
      - 1 -> .
      - 2 -> .
      - . -> .

    约束/假设：
      - 两个矩阵均为制表符分隔，**第一列为 `ID`（变体ID）**，首行是表头：`ID + sample_ids`；
      - 两矩阵的样本列顺序**一致**（若不一致请先用 `reorder_pre_to_post()` 处理）。
    设计：
      - **流式逐行**读取，常数内存；
      - 覆盖分组基于 `df_info[id_col, coverage_col]`，若某样本缺少覆盖度信息则跳过；
      - 覆盖文本先取 `|` 的第一个 token 并归一化为 `"15x"/"30x"`。

    输出：写入 TSV 文件（默认与 `out_post_order` 同目录；若提供 `chrom`，命名为 `chr{chrom}.coverage_transitions.tsv`，否则命名为 `chrALL.coverage_transitions.tsv`），
      列为：
        ID,
        15x_0_to_missing, 15x_1_to_missing, 15x_2_to_missing, 15x_missing_to_missing,
        30x_0_to_missing, 30x_1_to_missing, 30x_2_to_missing, 30x_missing_to_missing
    返回该 TSV 路径。
    """
    assert os.path.exists(out_pre_order), f"pre matrix not found: {out_pre_order}"
    assert os.path.exists(out_post_order), f"post matrix not found: {out_post_order}"

    # 初始化日志（若未设置）
    global LOG_FILE
    if LOG_FILE is None:
        try:
            base_dir = os.getcwd()
        except Exception:
            base_dir = "."
        LOG_FILE = os.path.join(base_dir, "geno_miss_bias.coverage_transitions.log")
        try:
            with open(LOG_FILE, "a") as lf:
                lf.write("="*72 + "\n")
                lf.write(f"[{_now()}][INFO] 开始执行 compute_coverage_transition_counts\n")
                lf.write(f"[{_now()}][INFO] pre={out_pre_order}; post={out_post_order}\n")
        except Exception:
            pass
    log_info(f"日志文件：{LOG_FILE}")

    # 读取表头并校验样本顺序
    with open(out_pre_order, "r") as fpre, open(out_post_order, "r") as fpost:
        hdr_pre = fpre.readline().rstrip("\n").split("\t")
        hdr_post = fpost.readline().rstrip("\n").split("\t")
    if not hdr_pre or hdr_pre[0] != "ID":
        log_err("pre 矩阵首列必须为 'ID'")
        raise ValueError("Bad header in pre matrix")
    if not hdr_post or hdr_post[0] != "ID":
        log_err("post 矩阵首列必须为 'ID'")
        raise ValueError("Bad header in post matrix")
    if hdr_pre[1:] != hdr_post[1:]:
        log_err("pre/post 样本列顺序不一致，请先调用 reorder_pre_to_post()")
        raise ValueError("Sample order mismatch between matrices")

    iids = hdr_post[1:]
    iid_to_idx = {iid: i for i, iid in enumerate(iids)}

    # 从 df_info 构造覆盖分组索引（列下标列表）
    if id_col not in df_info.columns:
        log_err(f"df_info 缺少列：{id_col}")
        raise KeyError(f"df_info missing column {id_col}")
    if coverage_col not in df_info.columns:
        log_err(f"df_info 缺少列：{coverage_col}")
        raise KeyError(f"df_info missing column {coverage_col}")

    def _first_token(s: Any) -> str:
        if pd.isna(s):
            return ""
        tok = str(s).split('|')[0].strip()
        return tok

    def _norm_cov(tok: str) -> str:
        # 归一化为 "15x" / "30x"（若无法解析则返回原样小写去空格）
        import re
        if not tok:
            return ""
        m = re.search(r"(\d+)\s*[xX]", tok)
        if m:
            return f"{m.group(1)}x"
        return tok.lower().strip()

    cov_a, cov_b = coverage_labels
    meta = df_info[[id_col, coverage_col]].copy()
    meta[id_col] = meta[id_col].astype(str)
    meta["__cov"] = meta[coverage_col].map(_first_token).map(_norm_cov)

    cov_idx = {cov_a: [], cov_b: []}
    missing_iids = 0
    for _, r in meta.iterrows():
        iid = r[id_col]
        cov = r["__cov"]
        if iid in iid_to_idx and cov in cov_idx:
            cov_idx[cov].append(iid_to_idx[iid])
        else:
            if iid in iid_to_idx:
                missing_iids += 1  # 样本存在但没有识别到 15x/30x
    if missing_iids:
        log_warn(f"有 {missing_iids} 个矩阵样本在 df_info 中未映射到 15x/30x，将不计入任何分组")
    log_info(f"覆盖分组样本数：{cov_a}={len(cov_idx[cov_a])}；{cov_b}={len(cov_idx[cov_b])}")

    # 组内样本总数（用于计算未 drop 到缺失的数量）
    n_cov_a = len(cov_idx[cov_a])
    n_cov_b = len(cov_idx[cov_b])
    # 输出路径与表头
    base_dir = os.path.dirname(os.path.abspath(out_post_order)) or os.getcwd()
    chrom_tag = f"chr{chrom}" if chrom else "chrALL"
    out_path = out_path or os.path.join(base_dir, f"{chrom_tag}.coverage_transitions.tsv")
    log_info(f"coverage 转换计数输出：{out_path}（chrom={chrom if chrom else 'ALL'}）")
    header_cols = [
        f"{cov_a}_0_to_missing", f"{cov_a}_1_to_missing", f"{cov_a}_2_to_missing", f"{cov_a}_missing_to_missing", f"{cov_a}_not_missing",
        f"{cov_b}_0_to_missing", f"{cov_b}_1_to_missing", f"{cov_b}_2_to_missing", f"{cov_b}_missing_to_missing", f"{cov_b}_not_missing",
    ]
    with open(out_path, "w") as fout:
        fout.write("ID\t" + "\t".join(header_cols) + "\n")

    # 流式逐行统计
    def _count_transitions(vals_pre: List[str], vals_post: List[str], idxs: List[int]) -> Tuple[int,int,int,int]:
        c0 = c1 = c2 = cm = 0
        for j in idxs:
            if j >= len(vals_pre) or j >= len(vals_post):
                continue
            a = vals_pre[j]
            b = vals_post[j]
            if a == '.' and b == '.':
                cm += 1
            elif a == '0' and b == '.':
                c0 += 1
            elif a == '1' and b == '.':
                c1 += 1
            elif a == '2' and b == '.':
                c2 += 1
        return c0, c1, c2, cm

    with open(out_pre_order, "r") as fpre, open(out_post_order, "r") as fpost, open(out_path, "a") as fout:
        _ = fpre.readline(); _ = fpost.readline()  # skip headers
        line_no = 0
        while True:
            lpre = fpre.readline()
            lpost = fpost.readline()
            if not lpre and not lpost:
                break
            if not lpre or not lpost:
                log_err("pre/post 行数不一致")
                raise RuntimeError("Row count mismatch between pre and post")
            line_no += 1
            sp = lpre.rstrip("\n").split("\t")
            sq = lpost.rstrip("\n").split("\t")
            vid_p = sp[0]; vid_q = sq[0]
            if vid_p != vid_q:
                log_err(f"第 {line_no} 行 ID 不一致：pre={vid_p}；post={vid_q}")
                raise RuntimeError("Row ID mismatch")
            vals_pre = sp[1:]
            vals_post = sq[1:]
            a0, a1, a2, am = _count_transitions(vals_pre, vals_post, cov_idx[cov_a])
            b0, b1, b2, bm = _count_transitions(vals_pre, vals_post, cov_idx[cov_b])
            a_keep = n_cov_a - (a0 + a1 + a2 + am)
            b_keep = n_cov_b - (b0 + b1 + b2 + bm)
            if a_keep < 0: a_keep = 0
            if b_keep < 0: b_keep = 0
            fout.write(vid_p + "\t" + "\t".join(map(str, (a0,a1,a2,am,a_keep,b0,b1,b2,bm,b_keep))) + "\n")

    try:
        with open(LOG_FILE, "a") as lf:
            lf.write(f"[{_now()}][INFO] 完成 compute_coverage_transition_counts；输出={out_path}\n")
            lf.write("="*72 + "\n")
    except Exception:
        pass

    return out_path


from contextlib import nullcontext

def summarize_coverage_transition_significance(
    transitions_tsv: str,
    coverage_labels: Tuple[str, str] = ("15x", "30x"),
    out_summary: Optional[str] = None,
    n_resamples: int = 9999,
    rng: int = 42,
    chrom: Optional[str] = None,
) -> str:
    """
    基于 `compute_coverage_transition_counts()` 产出的 `coverage_transitions.tsv`，为每个变体计算三层平台敏感性检验：
      L1) 2×2 精确检验（Fisher exact, two-sided）：
          行 = {15x, 30x}；列 = {missing_to_missing, other}
          其中 other = 0_to_missing + 1_to_missing + 2_to_missing + not_missing。
          结果命名为 `p_cov_L1`；若成功计算则 `warn_L1` 为空，否则写入原因。
      L2) 2×2 精确检验（Fisher exact, two-sided）：
          行 = {15x, 30x}；列 = {not_missing, drop_to_missing}
          其中 drop_to_missing = 0_to_missing + 1_to_missing + 2_to_missing。
          结果命名为 `p_cov_L2`；若成功计算则 `warn_L2` 为空，否则写入原因。
      L3) 2×3 Fisher Monte Carlo（仅使用精确蒙特卡罗，不再使用卡方近似）：
          行 = {15x, 30x}；列 = {0_to_missing, 1_to_missing, 2_to_missing}。
          **若至少两列在两行合计均为 0（即该列总和为 0），则不进行检验**，记录 `warn_L3` 并将 `p_cov_L3` 置为 `NA`；
          其余情况使用 `scipy.stats.fisher_exact` 的 `MonteCarloMethod` 进行检验，成功则 `warn_L3` 为空。

    规模可能很大：本函数**流式逐行**读取并写出 `stat_summary.tsv`，避免将整表载入内存。

    参数：
      transitions_tsv : `compute_coverage_transition_counts` 输出路径
      coverage_labels : 覆盖标签（默认 ("15x","30x")），需与 transitions 表头匹配
      out_summary     : 统计汇总输出路径（默认与输入同目录 `stat_summary.tsv`）
      n_resamples     : Monte Carlo 重抽样次数（若支持 RxC Fisher Monte Carlo）
      rng             : 随机数种子

    返回： out_summary 路径
    """
    assert os.path.exists(transitions_tsv), f"transitions file not found: {transitions_tsv}"

    # 初始化日志（若未设置）
    global LOG_FILE
    if LOG_FILE is None:
        try:
            base_dir = os.getcwd()
        except Exception:
            base_dir = "."
        LOG_FILE = os.path.join(base_dir, "geno_miss_bias.summary.log")
        try:
            with open(LOG_FILE, "a") as lf:
                lf.write("="*72 + "\n")
                lf.write(f"[{_now()}][INFO] 开始执行 summarize_coverage_transition_significance\n")
                lf.write(f"[{_now()}][INFO] transitions_tsv={transitions_tsv}\n")
        except Exception:
            pass
    log_info(f"日志文件：{LOG_FILE}")

    cov_a, cov_b = coverage_labels

    # 解析表头，确定所需列索引
    with open(transitions_tsv, "r") as fin:
        header = fin.readline().rstrip("\n").split("\t")
    required_cols = [
        f"{cov_a}_0_to_missing", f"{cov_a}_1_to_missing", f"{cov_a}_2_to_missing", f"{cov_a}_missing_to_missing", f"{cov_a}_not_missing",
        f"{cov_b}_0_to_missing", f"{cov_b}_1_to_missing", f"{cov_b}_2_to_missing", f"{cov_b}_missing_to_missing", f"{cov_b}_not_missing",
    ]
    missing = [c for c in required_cols if c not in header]
    if missing:
        log_err(f"coverage_transitions 缺少必要列：{missing}")
        raise KeyError("Missing required columns in transitions file")

    # 建立列索引
    idx = {name: header.index(name) for name in required_cols}
    id_idx = 0  # 首列 ID

    # 输出路径
    base_dir = os.path.dirname(os.path.abspath(transitions_tsv)) or os.getcwd()
    chrom_tag = f"chr{chrom}" if chrom else "chrALL"
    out_summary = out_summary or os.path.join(base_dir, f"{chrom_tag}.stat_summary.tsv")
    with open(out_summary, "w") as fout:
        fout.write(
            "ID\t"
            "p_missing_to_missing_vs_other\t"
            "p_not_missing_vs_drop_to_missing\t"
            "p_genotype_drop_composition\t"
            "warn_missing_to_missing_vs_other\t"
            "warn_not_missing_vs_drop_to_missing\t"
            "warn_genotype_drop_composition\n"
        )
    log_info(f"统计汇总输出：{out_summary}（chrom={chrom if chrom else 'ALL'}）")

    # 函数内工具：安全转换为整数
    def to_int(x: str) -> int:
        try:
            return int(x)
        except Exception:
            return 0

    # SciPy 检验函数准备（2x2 使用 fisher_exact；2x3 使用 fisher_exact + MonteCarloMethod）
    fisher_exact_fn = None
    MonteCarloMethodCls = None
    try:
        from scipy.stats import fisher_exact as _fisher_exact
        fisher_exact_fn = _fisher_exact
    except Exception as e:
        log_warn(f"无法导入 scipy.stats.fisher_exact：{e}")
    try:
        from scipy.stats import MonteCarloMethod as _MonteCarloMethod
        MonteCarloMethodCls = _MonteCarloMethod
    except Exception as e:
        log_warn(f"无法导入 scipy.stats.MonteCarloMethod（RxC Monte Carlo 将不可用）：{e}")

    # 逐行流式处理
    n_rows = 0
    n_warn = 0
    with open(transitions_tsv, "r") as fin, open(out_summary, "a") as fout:
        _ = fin.readline()  # skip header
        for line in fin:
            if not line:
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) <= max(idx.values()):
                continue
            vid = parts[id_idx]
            a0 = to_int(parts[idx[f"{cov_a}_0_to_missing"]])
            a1 = to_int(parts[idx[f"{cov_a}_1_to_missing"]])
            a2 = to_int(parts[idx[f"{cov_a}_2_to_missing"]])
            a_mm = to_int(parts[idx[f"{cov_a}_missing_to_missing"]])
            a_keep = to_int(parts[idx[f"{cov_a}_not_missing"]])
            b0 = to_int(parts[idx[f"{cov_b}_0_to_missing"]])
            b1 = to_int(parts[idx[f"{cov_b}_1_to_missing"]])
            b2 = to_int(parts[idx[f"{cov_b}_2_to_missing"]])
            b_mm = to_int(parts[idx[f"{cov_b}_missing_to_missing"]])
            b_keep = to_int(parts[idx[f"{cov_b}_not_missing"]])

            # --- (L1) 2x2 Fisher exact: missing_to_missing vs others ---
            # other = 0_to_missing + 1_to_missing + 2_to_missing + not_missing
            a_other = (a0 + a1 + a2 + a_keep)
            b_other = (b0 + b1 + b2 + b_keep)
            pL1 = "NA"; warn_L1 = ""
            try:
                if fisher_exact_fn is not None:
                    _, p = fisher_exact_fn([[a_mm, a_other], [b_mm, b_other]], alternative='two-sided')
                    pL1 = f"{p:.6g}"
                else:
                    warn_L1 = "scipy.stats.fisher_exact 不可用"
            except Exception as e:
                warn_L1 = f"L1 检验失败：{e}"

            # --- (L2) 2x2 Fisher exact: not_missing vs drop_to_missing ---
            drop_a = a0 + a1 + a2
            drop_b = b0 + b1 + b2
            pL2 = "NA"; warn_L2 = ""
            try:
                if fisher_exact_fn is not None:
                    _, p = fisher_exact_fn([[a_keep, drop_a], [b_keep, drop_b]], alternative='two-sided')
                    pL2 = f"{p:.6g}"
                else:
                    warn_L2 = "scipy.stats.fisher_exact 不可用"
            except Exception as e:
                warn_L2 = f"L2 检验失败：{e}"

            # --- (L3) 2x3 Fisher Monte Carlo（仅蒙特卡罗；无卡方近似） ---
            warn_L3 = ""
            cols_sum = [a0 + b0, a1 + b1, a2 + b2]
            zero_cols = sum(1 for s in cols_sum if s == 0)
            pL3 = "NA"
            if zero_cols >= 2:
                warn_L3 = (
                    "2x3 跳过：至少两列总和为0，信息不足；counts="
                    f"{{0:{(a0,b0)}, 1:{(a1,b1)}, 2:{(a2,b2)}}}"
                )
                n_warn += 1
            else:
                table = [[a0, a1, a2], [b0, b1, b2]]
                if fisher_exact_fn is not None and MonteCarloMethodCls is not None:
                    try:
                        import numpy as _np
                        rng_obj = _np.random.default_rng(rng)
                        method = MonteCarloMethodCls(n_resamples=n_resamples, rng=rng_obj)
                        res = fisher_exact_fn(table, method=method)
                        pL3 = f"{res.pvalue:.6g}" if hasattr(res, 'pvalue') else f"{res:.6g}"
                    except Exception as e:
                        warn_L3 = f"L3 Monte Carlo 失败：{e}"
                else:
                    warn_L3 = "RxC MonteCarloMethod 不可用（缺少新版本 SciPy），pL3=NA"

            # 写出一行
            fout.write(f"{vid}\t{pL1}\t{pL2}\t{pL3}\t{warn_L1}\t{warn_L2}\t{warn_L3}\n")
            n_rows += 1

    log_info(f"统计完成：总计 {n_rows} 个变体；警告 {n_warn} 条（2x3 未检验或蒙特卡罗失败）")
    try:
        with open(LOG_FILE, "a") as lf:
            lf.write(f"[{_now()}][INFO] 完成 summarize_coverage_transition_significance；输出={out_summary}\n")
            lf.write("="*72 + "\n")
    except Exception:
        pass

    return out_summary



def merge_bias_results_json(
    json_path: str,
    out_dir: Optional[str] = None,
    out_prefix: str = "all",
) -> Dict[str, str]:
    """
    合并 bias_results.json 中各染色体的结果文件，并产出新的汇总 JSON。

    输入：
      json_path : 由流程产出的 bias_results.json 路径。支持两种结构：
                  1) {"bias_results": [ { "chr": "chr1", "coverage_transitions": "...", "stat_summary": "..." }, ... ]}
                  2) [ { "chr": "chr1", "coverage_transitions": "...", "stat_summary": "..." }, ... ]
      out_dir   : （已忽略）输出目录强制为 **当前 Python 工作目录** `os.getcwd()`，便于在工作目录直接收集产物
      out_prefix: 合并结果文件名前缀（默认 "all"），将生成：
                  - <cwd>/<out_prefix>.coverage_transitions.tsv
                  - <cwd>/<out_prefix>.stat_summary.tsv
                  - <cwd>/bias_results.merged.json

    行为：
      - 按 chr1→chr22 的**自然顺序**进行拼接；
      - 仅第一个文件保留表头，后续文件跳过表头；
      - 流式拼接，避免大内存；
      - 生成新的 JSON，记录合并文件路径、包含的染色体列表和源 JSON 路径。

    返回：
      dict，包含：
        {
          "merged_coverage_transitions": "<path>",
          "merged_stat_summary": "<path>",
          "merged_json": "<path>"
        }
    """
    assert os.path.exists(json_path), f"json file not found: {json_path}"

    # 解析 JSON
    import json as _json
    with open(json_path, "r") as f:
        try:
            data = _json.load(f)
        except Exception as e:
            log_err(f"读取 JSON 失败：{e}")
            raise

    # 兼容两种顶层结构
    if isinstance(data, dict) and "bias_results" in data and isinstance(data["bias_results"], list):
        items = data["bias_results"]
    elif isinstance(data, list):
        items = data
    else:
        log_err("JSON 结构不符合预期：应为对象含 bias_results 列表或顶层为列表")
        raise ValueError("Unexpected JSON structure for bias_results")

    # 规范输出目录：强制使用当前 Python 工作目录
    out_dir = os.getcwd()
    _ensure_dir(out_dir)
    log_info(f"输出目录已固定为当前工作目录：{out_dir}")

    # 目标输出路径
    out_cov = os.path.join(out_dir, f"{out_prefix}.coverage_transitions.tsv")
    out_sum = os.path.join(out_dir, f"{out_prefix}.stat_summary.tsv")
    out_json = os.path.join(out_dir, "bias_results.merged.json")

    # 记录日志
    log_info(f"开始合并 bias_results：源={json_path}")
    log_info(f"输出：coverage_transitions={out_cov}；stat_summary={out_sum}；json={out_json}")

    # 去重并按 chr1..chr22 排序
    def _chr_key(c: str) -> int:
        # 提取数字部分，无法解析则置为 1e9 以排在末尾
        import re
        s = str(c)
        m = re.search(r'(\d+)$', s)
        return int(m.group(1)) if m else 10**9

    seen = set()
    by_chr: Dict[str, Dict[str, str]] = {}
    for it in items:
        try:
            c = it["chr"]
            cov = it["coverage_transitions"]
            summ = it["stat_summary"]
        except Exception:
            log_warn(f"条目缺失关键键，将跳过：{it}")
            continue
        if c in seen:
            log_warn(f"检测到重复染色体 {c}，仅保留首次出现的路径")
            continue
        seen.add(c)
        by_chr[c] = {"coverage": cov, "summary": summ}

    # 仅保留 chr1..chr22
    wanted = [f"chr{i}" for i in range(1, 23)]
    present = [c for c in wanted if c in by_chr]
    missing = [c for c in wanted if c not in by_chr]
    if missing:
        log_warn(f"以下染色体在 JSON 中缺失：{', '.join(missing)}")
    log_info(f"将按以下顺序拼接：{', '.join(present)}")

    # 工具：按顺序拼接 TSV（仅首个保留表头）
    def _concat_tsv_ordered(paths: List[str], out_path: str):
        n_written = 0
        with open(out_path, "w") as fout:
            for i, pth in enumerate(paths):
                if not os.path.exists(pth):
                    log_warn(f"文件不存在，跳过：{pth}")
                    continue
                if os.path.getsize(pth) == 0:
                    log_warn(f"文件为空，跳过：{pth}")
                    continue
                with open(pth, "r") as fin:
                    for j, line in enumerate(fin):
                        if i > 0 and j == 0:
                            continue  # 跳过表头
                        fout.write(line)
                        n_written += 1
        return n_written

    # 按序列出两个类型的文件清单
    cov_files = [by_chr[c]["coverage"] for c in present]
    sum_files = [by_chr[c]["summary"] for c in present]

    # 拼接
    nw_cov = _concat_tsv_ordered(cov_files, out_cov)
    nw_sum = _concat_tsv_ordered(sum_files, out_sum)
    log_info(f"拼接完成：coverage_transitions 写入 {nw_cov} 行；stat_summary 写入 {nw_sum} 行")

    # 生成新的 JSON
    import json as _json2
    merged_obj = {
        "source_json": os.path.abspath(json_path),
        "chromosomes": present,
        "merged": {
            "coverage_transitions": os.path.abspath(out_cov),
            "stat_summary": os.path.abspath(out_sum),
        }
    }
    try:
        with open(out_json, "w") as jf:
            _json2.dump(merged_obj, jf, ensure_ascii=False, indent=2)
        log_info(f"已写出汇总 JSON：{out_json}")
    except Exception as e:
        log_err(f"写出汇总 JSON 失败：{e}")
        raise

    # 确保返回绝对路径
    out_cov = os.path.abspath(out_cov)
    out_sum = os.path.abspath(out_sum)
    out_json = os.path.abspath(out_json)
    return {
        "merged_coverage_transitions": out_cov,
        "merged_stat_summary": out_sum,
        "merged_json": out_json,
    }


def adjust_stat_summary_fdr(
    merged_json_path: str,
    out_path: Optional[str] = None,
    pcols: Tuple[str, str, str] = (
        "p_missing_to_missing_vs_other",
        "p_not_missing_vs_drop_to_missing",
        "p_genotype_drop_composition",
    ),
    float_fmt: str = ".6g",
) -> str:
    """
    基于合并 JSON（由 `merge_bias_results_json()` 产出）的 `stat_summary` 路径，
    对三个 p 值列分别执行 Benjamini–Hochberg (BH/FDR) 矫正，并将校正后的 q 值
    作为新列写回一个新的 TSV。

    参数
    ----
    merged_json_path : str
        `bias_results.merged.json` 的路径。其 JSON 结构需包含：
        {
          "merged": { "stat_summary": "/path/to/all.stat_summary.tsv", ... },
          ...
        }
    out_path : Optional[str]
        输出 TSV 路径；默认在 `stat_summary` 同目录、文件名追加后缀 `.fdr.tsv`。
    pcols : Tuple[str, str, str]
        需校正的三列名称，默认：
          - "p_missing_to_missing_vs_other"
          - "p_not_missing_vs_drop_to_missing"
          - "p_genotype_drop_composition"
    float_fmt : str
        输出数值格式（传给 f-string，如 ".6g"）。

    行为
    ----
    * 对每一列 **独立** 进行 BH 校正；
    * 无法解析为浮点的值（如 "NA"/"na"/空）**不参与**校正，同时在输出中原样回填；
    * 为避免内存溢出，采用“两遍扫描”并**流式写出**结果：
        - 第1遍：仅收集可用 p 值，完成每列的 BH 校正，构建“原始 token → 按升序出现的校正值队列”的映射；
        - 第2遍：逐行读取原始文件，遇到可用 p 值则从对应队列**弹出**一个校正值写出，否则原样写出无效 token。

    返回
    ----
    str : 新的带 FDR 列的 stat_summary 路径。
    """
    assert os.path.exists(merged_json_path), f"merged json not found: {merged_json_path}"

    # 初始化日志文件（若尚未设置）
    global LOG_FILE
    if LOG_FILE is None:
        try:
            base_dir = os.getcwd()
        except Exception:
            base_dir = "."
        LOG_FILE = os.path.join(base_dir, "geno_miss_bias.fdr.log")
        try:
            with open(LOG_FILE, "a") as lf:
                lf.write("="*72 + "\n")
                lf.write(f"[{_now()}][INFO] 开始执行 adjust_stat_summary_fdr\n")
                lf.write(f"[{_now()}][INFO] merged_json_path={merged_json_path}\n")
        except Exception:
            pass
    log_info(f"日志文件：{LOG_FILE}")

    # 读取 JSON，获取 stat_summary 路径
    import json as _json
    with open(merged_json_path, "r") as jf:
        data = _json.load(jf)
    try:
        stat_summary_path = data["merged"]["stat_summary"]
    except Exception:
        log_err("合并 JSON 中缺少 merged.stat_summary 路径")
        raise KeyError("merged.stat_summary missing")

    assert os.path.exists(stat_summary_path), f"stat_summary not found: {stat_summary_path}"

    # 输出路径：默认写在当前工作目录
    if out_path is None:
        base_dir = os.getcwd()
        root, ext = os.path.splitext(os.path.basename(stat_summary_path))
        out_path = os.path.join(base_dir, root + ".fdr.tsv")
    log_info(f"输入 stat_summary：{stat_summary_path}")
    log_info(f"输出（带 FDR）：{out_path}")

    # 解析表头，定位列索引
    with open(stat_summary_path, "r") as fin:
        header = fin.readline().rstrip("\n").split("\t")
    missing_cols = [c for c in pcols if c not in header]
    if missing_cols:
        log_err(f"stat_summary 缺少必要 p 列：{missing_cols}")
        raise KeyError(f"Missing p columns: {missing_cols}")
    pidx = tuple(header.index(c) for c in pcols)

    # 工具：数字解析（严格 0-1 也可放宽，这里仅作 float 解析）
    def _parse_p(tok: str):
        if tok is None:
            return None, False
        s = tok.strip()
        if s == "" or s.lower() in ("na", "nan", "null", "."):
            return None, False
        try:
            v = float(s)
            if not (v >= 0.0 and v <= 1.0):
                # 超界 p 值视作无效
                return None, False
            return v, True
        except Exception:
            return None, False

    # 第一遍：收集各列的可用 p 值（按**出现顺序**记录），为 FDR 做准备
    import numpy as _np
    col_vals = [[], [], []]  # 每列一个列表，按出现顺序追加 float p

    with open(stat_summary_path, "r") as fin:
        _ = fin.readline()  # 跳过表头
        for line in fin:
            if not line:
                continue
            parts = line.rstrip("\n").split("\t")
            for ci, pi in enumerate(pidx):
                tok = parts[pi] if pi < len(parts) else ""
                v, ok = _parse_p(tok)
                if ok:
                    col_vals[ci].append(float(v))

    # 使用 SciPy 的 BH/FDR（Benjamini–Hochberg）实现
    try:
        from scipy import stats as _stats
    except Exception as e:
        log_err(f"无法导入 scipy.stats.false_discovery_control：{e}")
        raise

    # 对三列分别独立校正；返回值顺序与输入 p 顺序一致
    qvals_list = []
    for ci in range(3):
        if len(col_vals[ci]) == 0:
            qvals_list.append(_np.array([], dtype=float))
        else:
            try:
                qv = _stats.false_discovery_control(_np.asarray(col_vals[ci], dtype=float), method='bh')
            except TypeError:
                # 兼容较老的 SciPy（某些版本需要 keyword `method` 或不支持）；如失败则报错
                log_err("当前 SciPy 版本不支持 stats.false_discovery_control(method='bh')，请升级 SciPy ≥ 1.11")
                raise
            qvals_list.append(_np.asarray(qv, dtype=float))

    # 为第二遍输出准备按出现顺序的队列（deque），遇到一个可用 p 就弹出一个对应的 q
    from collections import deque
    q_deques = [deque(qvals_list[0].tolist()), deque(qvals_list[1].tolist()), deque(qvals_list[2].tolist())]

    # 第二遍：逐行写出，追加 3 列 FDR（无效 token 原样回填）
    fdr_headers = [c.replace("p_", "fdr_") for c in pcols]
    with open(out_path, "w") as fout, open(stat_summary_path, "r") as fin:
        # 写表头
        hdr = fin.readline().rstrip("\n")
        fout.write(hdr + "\t" + "\t".join(fdr_headers) + "\n")
        # 逐行
        for line in fin:
            if not line:
                continue
            parts = line.rstrip("\n").split("\t")
            fdr_out = []
            for ci, pi in enumerate(pidx):
                tok = parts[pi] if pi < len(parts) else ""
                v, ok = _parse_p(tok)
                if ok:
                    if q_deques[ci]:
                        q = q_deques[ci].popleft()
                        fdr_out.append(format(float(q), float_fmt))
                    else:
                        # 意外：队列已空，保底写回原始 p
                        fdr_out.append(tok)
                else:
                    fdr_out.append(tok)
            # 写行
            fout.write("\t".join(parts + fdr_out) + "\n")

    log_info("BH/FDR 校正完成并写出新文件")
    return out_path

