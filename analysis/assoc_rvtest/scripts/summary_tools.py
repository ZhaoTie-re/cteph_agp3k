
def update_json_manifest(json_path: str, num_var_thr: int = 2, out_path: str = None) -> str: # type: ignore
    """
    基于输入的 JSON 清单文件，读取 association 结果（no_macmin / macmin_5 / macmin_10）
    的 .assoc 文件，统计：
      1) 各文件包含的基因数（不含表头行）
      2) 各文件中满足 NumVar >= num_var_thr 的基因数（不含表头行）
      3) 计算显著性阈值 Sig_level = 0.05 / 通过数量
    并把统计结果回写到新的 JSON 文件中。

    补充：生成 `<stem>.summary.significant.csv`，在整合显著基因的基础上，新增 3 列
    `Reason_no_macmin`、`Reason_macmin_5`、`Reason_macmin_10` 标记各条件下未显著的原因：
      1 = 原始 association 结果中不存在该基因（如 MAC 过滤导致）
      2 = 存在记录，但 NumVar < num_var_thr 因而被过滤
      3 = 存在记录且满足 NumVar 阈值，但 Pvalue >= 对应 SigLevel
    若该条件已显著，则对应的 Reason 留空。

    参数
    ------
    json_path : str
        现有清单 JSON 文件路径，包含 association_results 字段。
    num_var_thr : int, default=2
        NumVar 的阈值（大于等于该阈值计数）。
    out_path : Optional[str]
        输出 JSON 文件路径；若为 None，则在同目录创建 `<原名>.updated.json`。

    返回
    ------
    str
        新 JSON 文件路径。
    """
    import json, os, re, math, datetime, concurrent.futures

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
    for k, stat in summary.items():
        cleaned = dict(stat) if isinstance(stat, dict) else {}
        # 移除内部使用的原始记录，避免 tuple 作为键导致 JSON 序列化失败
        if "raw_records" in cleaned:
            cleaned.pop("raw_records", None)
        meta["association_summary"][k] = cleaned

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
