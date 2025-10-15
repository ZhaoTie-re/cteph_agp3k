
"""
summary_tools.py — 关联结果汇总与导出工具
================================================

概述
----
本脚本用于在 RVTESTS/SAIGE 等基因集合关联分析流程结束后，对多个条件（如 no_macmin、macmin_5、macmin_10）的结果进行**并行读取与统计**，并将关键产出物与统计信息**结构化写回 JSON**。同时，它还能基于给定 VCF 路径以 **bcftools→awk→bgzip** 的**流式**方式快速导出等位计数表，避免高内存占用。

主要功能
--------
1. 读取 `association_results` 中的 `.assoc` 文件，统计：
   - `n_genes_total`：不含表头的基因条目数；
   - `n_genes_pass`：满足 `NumVar ≥ num_var_thr` 的条目数；
   - `sig_level`：`0.05 / n_genes_pass`；
   - 生成各条件的 `*.significant.tsv`（含 `Gene,RANGE,NumVar,Pvalue,SigLevel`，按 `Pvalue` 升序；`Gene/RANGE` 会按逗号切分去重后再合并）。
2. 生成整合表 `summary.significant.csv`：
   - 列：`Gene,RANGE,NumVar,Pvalue,SigLevel,Reason_no_macmin,Reason_macmin_5,Reason_macmin_10`；
   - 其中 `NumVar/Pvalue/SigLevel` 为按条件组织的 JSON；
   - `Reason_*` 原因编码：`1`=原始结果缺失（如 MAC 过滤导致）、`2`=存在记录但 `NumVar` 未达阈值、`3`=存在记录且达阈值但不显著。
3. 基于 `input_files.vcf_file` 导出 `allele_counts.tsv.gz`：
   - 字段：`CHROM,POS,ID,REF,ALT,AN,AC,MAC`；
   - 多等位位点按 `ALT` 展开为多行；
   - 使用 `tabix -S 1` 建立索引，自动记录到 JSON 的 `exported_files`。
4. 统计 `allele_counts.tsv.gz` 的观测量：
   - `records_excl_header`：总行数（排除表头）；
   - `per_condition_records`：`no_macmin`（等于总行数）、`macmin_5`（`MAC≥5` 行数）、`macmin_10`（`MAC≥10` 行数）。

输出 JSON 的关键字段
--------------------
- `association_summary.num_var_thr`：本次阈值；
- `association_summary.non_sig_reason_legend`：原因编码图例；
- `association_summary.{no_macmin|macmin_5|macmin_10}`：各条件统计（见上）；
- `association_summary.summary_significant_csv`：整合 CSV 绝对路径；
- `association_summary.allele_counts_summary`：导出表的行数统计与条件行数；
- `exported_files.allele_counts_tsv_gz` 及其索引路径。

依赖与环境
----------
- 需要可执行：`bcftools`、`awk`、`bgzip`、`tabix`；
- 本脚本默认 `bcftools_path="/home/b/b37974/bcftools/bcftools"`；
- `bgzip -@ <threads>` 负责压缩时的多线程；`bcftools query` 不依赖 `--threads`；
- 日志 Logger 名称为 `summary_tools`，输出中文提示便于定位问题。

使用示例
--------
>>> from summary_tools import update_json_manifest
>>> update_json_manifest("meta.updated.json", num_var_thr=3)
"/当前工作目录/meta.updated.json"

注意事项
--------
- `INFO/AN` 与 `INFO/AC` 缺失的记录会在 `awk` 阶段被跳过；
- 头行以 `tabix -S 1` 跳过建索引；
- `Gene`/`RANGE` 会被按逗号分隔去重后再写出；
- 仅在存在显著基因时才会生成对应的 `*.significant.tsv` 与整合 CSV。

作者
----
- ZHAO TIE
"""

def update_json_manifest(
    json_path: str, 
    num_var_thr: int = 2, 
    out_path: str = None, # type: ignore
    bcftools_path: str = "/home/b/b37974/bcftools/bcftools", 
    threads: int = 8
    ) -> str: 
    """
    基于清单 JSON 统计关联结果并导出汇总文件（简要说明）。

    功能：
    - 并行读取 `no_macmin`、`macmin_5`、`macmin_10` 的 `.assoc` 文件，统计 `n_genes_total`、`n_genes_pass`、`sig_level`；
    - 生成各条件的 `*.significant.tsv`（`Gene,RANGE,NumVar,Pvalue,SigLevel`，`Pvalue` 升序；`Gene/RANGE` 去重合并）；
    - 合并生成 `summary.significant.csv`，并给出每个条件未显著的原因编码（1/2/3，详见模块文档）；
    - 基于 `input_files.vcf_file` 流式导出 `allele_counts.tsv.gz` 并索引；统计总行数及 `MAC≥5/10` 行数；
    - 将以上结果结构化写回新的 JSON。

    参数：
    - json_path (str)：输入 JSON 路径；
    - num_var_thr (int, 默认 2)：`NumVar` 判定阈值；
    - out_path (str|None)：输出 JSON 路径（默认写到当前工作目录）；
    - bcftools_path (str, 默认 "/home/b/b37974/bcftools/bcftools")：bcftools 可执行路径；
    - threads (int, 默认 8)：压缩与相关步骤的线程数。

    返回：
    - (str) 新 JSON 文件的绝对路径。
    """
    import json, os, re, math, datetime, concurrent.futures, subprocess, shlex, logging, gzip

    def _read_assoc_counts(path: str, thr: int, sig_out_path: str = None) -> dict: # type: ignore
        """读取 .assoc 文件并返回统计字典。
        允许以空白分隔（\t/空格），自动识别表头，大小写不敏感地寻找 'NumVar' 列。
        返回: {"n_genes_total": int, "n_genes_pass": int, "sig_level": float or None, "n_genes_sig": int, "sig_list_path": str or None, "raw_records": dict}
        若文件不存在或为空，返回 0/0/None。
        """
        if not path or not os.path.exists(path):
            return {"n_genes_total": 0, "n_genes_pass": 0, "sig_level": None, "n_genes_sig": 0, "sig_list_path": None, "raw_records": {}}

        n_total = 0
        n_pass = 0
        numvar_idx = None
        gene_idx = None
        range_idx = None
        pval_idx = None
        header_seen = False
        candidates = []  # (gene, range_, pval, nv)
        raw_records = {}  # 记录所有观测，用于后续不显著原因分析

        # 为了兼容 RVTESTS 的输出，这里使用正则按任意空白分割
        splitter = re.compile(r"\s+")

        def _unique_join(val: str) -> str:
            parts = [v.strip() for v in val.split(',') if v.strip()]
            uniq = list(dict.fromkeys(parts))
            return ','.join(uniq)

        with open(path, "r", encoding="utf-8", errors="ignore") as fh:
            for line in fh:
                line = line.rstrip("\n")
                if not line:
                    continue
                cols = splitter.split(line)
                # 尝试识别表头：包含 NumVar/NUMVAR 字段
                if not header_seen:
                    lower = [c.lower() for c in cols]
                    if "numvar" in lower:
                        numvar_idx = lower.index("numvar")
                    elif "nvar" in lower:
                        numvar_idx = lower.index("nvar")
                    else:
                        # 如果第一行不是表头，则继续读取直到遇到含 NumVar 的行（极端情况）
                        continue
                    # 识别其他列索引
                    gene_idx = lower.index("gene") if "gene" in lower else None
                    range_idx = lower.index("range") if "range" in lower else None
                    if "pvalue" in lower:
                        pval_idx = lower.index("pvalue")
                    elif "pval" in lower:
                        pval_idx = lower.index("pval")
                    else:
                        pval_idx = None
                    header_seen = True
                    continue  # 表头不计入
                # 到这里说明已识别表头
                if numvar_idx is None or numvar_idx >= len(cols):
                    # 异常行，跳过
                    continue
                n_total += 1
                try:
                    nv = int(float(cols[numvar_idx]))
                except ValueError:
                    continue
                # 记录所有观测，用于后续不显著原因分析
                rec_gene = _unique_join(cols[gene_idx]) if gene_idx is not None and gene_idx < len(cols) else ''
                rec_range = _unique_join(cols[range_idx]) if range_idx is not None and range_idx < len(cols) else ''
                rec_pval = None
                try:
                    if pval_idx is not None and pval_idx < len(cols):
                        rec_pval = float(cols[pval_idx])
                except ValueError:
                    rec_pval = None
                raw_records[(rec_gene, rec_range)] = {"NumVar": nv, "Pvalue": rec_pval}
                if nv >= thr:
                    n_pass += 1
                    pval = None
                    try:
                        if pval_idx is not None and pval_idx < len(cols):
                            pval = float(cols[pval_idx])
                    except ValueError:
                        pval = None
                    gene = rec_gene
                    range_ = rec_range
                    candidates.append((gene, range_, pval, nv))

        sig = (0.05 / n_pass) if n_pass > 0 else None

        n_sig = 0
        written_path = None
        if sig is not None and candidates:
            # Collect rows meeting significance, including NumVar and Pvalue
            sig_rows = [ (g, r, nv, p) for (g, r, p, nv) in candidates if isinstance(p, float) and p < sig ]
            sig_rows.sort(key=lambda x: (float('inf') if x[3] is None else x[3]))
            n_sig = len(sig_rows)
            if sig_out_path and n_sig > 0:
                abs_path = os.path.abspath(sig_out_path)
                with open(abs_path, 'w', encoding='utf-8') as wf:
                    wf.write('Gene\tRANGE\tNumVar\tPvalue\tSigLevel\n')
                    for g, r, nv, p in sig_rows:
                        # Use compact formatting for p-values
                        p_str = (f"{p:.6g}" if isinstance(p, float) else str(p))
                        wf.write(f"{g}\t{r}\t{nv}\t{p_str}\t{sig}\n")
                written_path = abs_path

        return {"n_genes_total": n_total, "n_genes_pass": n_pass, "sig_level": sig, "n_genes_sig": n_sig, "sig_list_path": written_path, "raw_records": raw_records}

    # 读取原 JSON
    with open(json_path, "r", encoding="utf-8") as f:
        meta = json.load(f)

    assoc = (meta.get("association_results") or {})
    paths = {
        "no_macmin": assoc.get("no_macmin"),
        "macmin_5": assoc.get("macmin_5"),
        "macmin_10": assoc.get("macmin_10"),
    }

    cwd = os.getcwd()
    base_name = os.path.basename(json_path)
    stem, _ = os.path.splitext(base_name)
    sig_paths = {
        "no_macmin": os.path.join(cwd, f"{stem}.no_macmin.significant.tsv"),
        "macmin_5": os.path.join(cwd, f"{stem}.macmin_5.significant.tsv"),
        "macmin_10": os.path.join(cwd, f"{stem}.macmin_10.significant.tsv"),
    }

    summary = {}
    with concurrent.futures.ThreadPoolExecutor(max_workers=min(4, len(paths))) as ex:
        future_map = {ex.submit(_read_assoc_counts, paths[k], num_var_thr, sig_paths[k]): k for k in paths} # type: ignore
        for fut in concurrent.futures.as_completed(future_map):
            k = future_map[fut]
            try:
                summary[k] = fut.result()
            except Exception as e:
                summary[k] = {"n_genes_total": 0, "n_genes_pass": 0, "sig_level": None, "n_genes_sig": 0, "sig_list_path": None, "error": str(e)}

    # ---------- 生成整合的 summary.significant.csv ----------
    # 优先使用实际写入路径（summary[k]['sig_list_path']），以防某些条件下无显著结果
    cond_order = ["no_macmin", "macmin_5", "macmin_10"]
    cond_to_path = {k: (summary.get(k, {}) or {}).get("sig_list_path") for k in cond_order}

    merged = {}  # key=(Gene, RANGE) -> {"NumVar": {cond: int}, "Pvalue": {cond: float}, "SigLevel": {cond: float}}

    def _safe_float(x):
        try:
            return float(x)
        except Exception:
            return None

    def _unique_join(val: str) -> str:
        parts = [v.strip() for v in val.split(',') if v.strip()]
        uniq = list(dict.fromkeys(parts))
        return ','.join(uniq)

    # Build sig_key_map: cond -> set of (gene, rng) that are significant in that condition
    sig_key_map = {}
    for cond in cond_order:
        tsv = cond_to_path.get(cond)
        sig_keys = set()
        if not tsv or not os.path.exists(tsv):
            sig_key_map[cond] = sig_keys
            continue
        with open(tsv, "r", encoding="utf-8") as fh:
            header = fh.readline().rstrip("\n").split("\t")
            # Expect: Gene, RANGE, NumVar, Pvalue, SigLevel
            try:
                gene_i = header.index("Gene")
                range_i = header.index("RANGE")
                nv_i = header.index("NumVar")
                pv_i = header.index("Pvalue")
                sl_i = header.index("SigLevel")
            except ValueError:
                # 表头异常则跳过该文件
                sig_key_map[cond] = sig_keys
                continue
            for line in fh:
                if not line.strip():
                    continue
                cols = line.rstrip("\n").split("\t")
                if max(gene_i, range_i, nv_i, pv_i, sl_i) >= len(cols):
                    continue
                gene = cols[gene_i]
                rng = cols[range_i]
                gene = _unique_join(gene)
                rng = _unique_join(rng)
                key = (gene, rng)
                sig_keys.add(key)
                try:
                    nv = int(float(cols[nv_i]))
                except Exception:
                    nv = None
                pv = _safe_float(cols[pv_i])
                sl = _safe_float(cols[sl_i])
                if key not in merged:
                    merged[key] = {"NumVar": {}, "Pvalue": {}, "SigLevel": {}}
                if nv is not None:
                    merged[key]["NumVar"][cond] = nv
                if pv is not None:
                    merged[key]["Pvalue"][cond] = pv
                if sl is not None:
                    merged[key]["SigLevel"][cond] = sl
        sig_key_map[cond] = sig_keys

    # 写出 CSV（仅当 merged 非空）
    summary_sig_csv = None
    if merged:
        summary_sig_csv = os.path.abspath(os.path.join(cwd, f"{stem}.summary.significant.csv"))
        import csv as _csv

        def _reason_for(cond, key):
            # 若该条件已显著，返回空字符串
            if key in sig_key_map.get(cond, set()):
                return ""
            # 取该条件的原始记录与阈值
            raw = (summary.get(cond) or {}).get("raw_records") or {}
            sig_lv = (summary.get(cond) or {}).get("sig_level")
            info = raw.get(key)
            if info is None:
                return "1"  # 1. 原始结果中不存在（例如 MAC 策略导致）
            nv = info.get("NumVar")
            pv = info.get("Pvalue")
            if nv is None or nv < num_var_thr:
                return "2"  # 2. 存在记录但 NumVar < 阈值
            # 到这里说明 NumVar>=thr，但不显著
            return "3"     # 3. 有记录且满足阈值，但 Pvalue>=SigLevel

        with open(summary_sig_csv, "w", encoding="utf-8", newline="") as wf:
            writer = _csv.writer(wf)
            writer.writerow(["Gene", "RANGE", "NumVar", "Pvalue", "SigLevel", "Reason_no_macmin", "Reason_macmin_5", "Reason_macmin_10"])
            for (gene, rng), vals in merged.items():
                # 为保证键顺序一致，按 cond_order 排序
                def _ordered(d):
                    return {k: d[k] for k in cond_order if k in d}
                numvar_json = json.dumps(_ordered(vals["NumVar"]), ensure_ascii=False, separators=(",", ":"))
                pvalue_json = json.dumps(_ordered(vals["Pvalue"]), ensure_ascii=False, separators=(",", ":"))
                siglvl_json = json.dumps(_ordered(vals["SigLevel"]), ensure_ascii=False, separators=(",", ":"))
                r1 = _reason_for("no_macmin", (gene, rng))
                r2 = _reason_for("macmin_5", (gene, rng))
                r3 = _reason_for("macmin_10", (gene, rng))
                writer.writerow([gene, rng, numvar_json, pvalue_json, siglvl_json, r1, r2, r3])

    # ---------- 依据 VCF 导出等位计数表 (tsv.gz) ----------
    logger = logging.getLogger("summary_tools")
    if not logger.handlers:
        _h = logging.StreamHandler()
        _fmt = logging.Formatter("[%(asctime)s] %(levelname)s - %(message)s")
        _h.setFormatter(_fmt)
        logger.addHandler(_h)
        logger.setLevel(logging.INFO)

    vcf_file = (meta.get("input_files") or {}).get("vcf_file")
    allele_counts_tsv = None
    allele_counts_tbi = None
    if vcf_file and os.path.exists(vcf_file):
        allele_counts_tsv = os.path.abspath(os.path.join(cwd, f"{stem}.allele_counts.tsv.gz"))
        allele_counts_tbi = allele_counts_tsv + ".tbi"
        logger.info("开始导出等位计数表（流式）：CHROM POS ID REF ALT AN AC MAC ..")
        # 使用 bcftools query 流式输出，然后用 awk 展开多等位并计算 MAC，再 bgzip 压缩
        # 字段说明：
        # $1=CHROM $2=POS $3=ID $4=REF $5=ALT(逗号分隔) $6=AN $7=AC(逗号分隔)
        awk_code = r'BEGIN{OFS="\t"; print "CHROM","POS","ID","REF","ALT","AN","AC","MAC"} {split($5, alts, ","); split($7, acs, ","); an=$6; n=(length(alts)>length(acs)?length(alts):length(acs)); for(i=1;i<=n;i++){ alt=alts[i]; ac=acs[i]; if(alt==""||ac==""||ac==".") continue; mac=(an-ac); if(ac<mac) mac=ac; print $1,$2,$3,$4,alt,an,ac,mac; }}'
        fmt = "%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/AN\t%INFO/AC\n"
        query_cmd = f"{shlex.quote(bcftools_path)} query -u -f {shlex.quote(fmt)} {shlex.quote(vcf_file)}"
        shell_cmd = (
            f"set -euo pipefail; \n" \
            f"{query_cmd} | awk '{awk_code}' | bgzip -@ {int(threads)} -c > {shlex.quote(allele_counts_tsv)} && "
            f"tabix -f -s 1 -b 2 -e 2 -S 1 {shlex.quote(allele_counts_tsv)}"
        )
        try:
            logger.info("执行命令：" + shell_cmd)
            subprocess.run(["/bin/bash", "-lc", shell_cmd], check=True)
            logger.info("等位计数表导出完成：" + allele_counts_tsv)
        except subprocess.CalledProcessError as e:
            logger.error(f"等位计数表导出失败，退出码={e.returncode}")
            allele_counts_tsv = None
            allele_counts_tbi = None
    else:
        logger.warning("未找到可用的 VCF 文件路径，跳过等位计数表导出。")

    # ---------- 统计 allele_counts.tsv.gz 的行数（排除表头）及按 MAC 阈值的行数 ----------
    allele_counts_src = allele_counts_tsv
    if not allele_counts_src:
        allele_counts_src = ((meta.get("exported_files") or {}).get("allele_counts_tsv_gz"))
    allele_counts_summary = None
    if allele_counts_src and os.path.exists(allele_counts_src):
        total_excl_header = 0
        ge5 = 0
        ge10 = 0
        mac_idx = None
        try:
            with gzip.open(allele_counts_src, "rt", encoding="utf-8", errors="ignore") as fh:
                header = fh.readline()
                if header:
                    cols = header.rstrip("\n").split("\t")
                    # 期望表头: CHROM POS ID REF ALT AN AC MAC
                    lower = [c.lower() for c in cols]
                    if "mac" in lower:
                        mac_idx = lower.index("mac")
                for line in fh:
                    if not line.strip():
                        continue
                    total_excl_header += 1
                    if mac_idx is None:
                        continue
                    parts = line.rstrip("\n").split("\t")
                    if mac_idx >= len(parts):
                        continue
                    try:
                        mac_val = float(parts[mac_idx])
                    except Exception:
                        continue
                    if mac_val >= 5:
                        ge5 += 1
                    if mac_val >= 10:
                        ge10 += 1
        except Exception as e:
            logger.warning(f"读取 {allele_counts_src} 统计失败: {e}")
        allele_counts_summary = {
            "source": os.path.abspath(allele_counts_src),
            "records_excl_header": total_excl_header,
            "per_condition_records": {
                "no_macmin": total_excl_header,  # 不过滤
                "macmin_5": ge5,
                "macmin_10": ge10,
            },
        }
    else:
        logger.warning("未找到 allele_counts.tsv.gz，无法统计观测数量。")

    # 把统计信息写回到新的字段中
    meta.setdefault("association_summary", {})
    meta["association_summary"]["num_var_thr"] = num_var_thr
    meta["association_summary"]["non_sig_reason_legend"] = {
        "1": "not present in association results (e.g., MAC setting)",
        "2": "present but filtered by NumVar < threshold",
        "3": "present and NumVar>=threshold but not significant (Pvalue >= SigLevel)"
    }
    if summary_sig_csv:
        meta["association_summary"]["summary_significant_csv"] = summary_sig_csv
    if allele_counts_tsv:
        meta.setdefault("exported_files", {})
        meta["exported_files"]["allele_counts_tsv_gz"] = allele_counts_tsv
        if allele_counts_tbi is not None and os.path.exists(allele_counts_tbi):
            meta["exported_files"]["allele_counts_tsv_gz_tbi"] = allele_counts_tbi
    for k, stat in summary.items():
        cleaned = dict(stat) if isinstance(stat, dict) else {}
        # 移除内部使用的原始记录，避免 tuple 作为键导致 JSON 序列化失败
        if "raw_records" in cleaned:
            cleaned.pop("raw_records", None)
        meta["association_summary"][k] = cleaned
    if allele_counts_summary is not None:
        meta["association_summary"]["allele_counts_summary"] = allele_counts_summary

    # 生成输出文件路径
    if out_path is None:
        cwd = os.getcwd()
        base_name = os.path.basename(json_path)
        stem, _ = os.path.splitext(base_name)
        out_path = os.path.join(cwd, f"{stem}.updated.json")
    # 附加更新时间戳
    meta["updated_time"] = datetime.datetime.now().isoformat(timespec="seconds")

    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(meta, f, ensure_ascii=False, indent=2)

    return out_path