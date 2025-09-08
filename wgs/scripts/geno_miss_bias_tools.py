#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import sys
import shutil
import subprocess
import tempfile
import time
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
) -> Any:
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
      7) 日志：所有进度/告警/错误信息同步写入 `<work_dir>/geno_miss_bias.log`，同时打印到 stderr。

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
      - 运行日志保存在 `<work_dir>/geno_miss_bias.log`，建议先查看日志定位问题（如 chr 前缀、REF/ALT 顺序、VCF 空 ID 等）。
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
    LOG_FILE = os.path.join(base_tmp, "geno_miss_bias.log")
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
            iid_cols = [c.split("_")[-1] for c in sample_cols]
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
