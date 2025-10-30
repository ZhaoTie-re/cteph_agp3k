"""
rvtest_prepare_tools.py - RVTest 数据预处理工具集

本模块提供了一系列用于准备 RVTest（Rare Variant Tests）分析所需数据文件的工具函数。
主要包括以下功能：

1. plink_to_vcf_raw: 将 PLINK 格式数据转换为标准化的 VCF 格式
2. reformat_pheno_covar: 重新格式化表型和协变量文件以适配 RVTest
3. reformat_refflat_remove_chr: 处理 refFlat 基因注释文件，去除染色体前缀

作者: ZHAO TIE
最后更新: 2025年10月8日
"""

def plink_to_vcf_raw(
    bed_prefix: str,
    plink_path: str = "/home/b/b37974/plink",
    bcftools_path: str = "/home/b/b37974/bcftools/bcftools",
    ref_fa: str = "/LARGE0/gr10478/b37974/Pulmonary_Hypertension/nagasaki_pipeline/data/hs38DH.fa",
    threads: int = 6,
    rename_chr: str = "/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/info/rename_chr.txt",
    snps_only_just_acgt: bool = False,
    tabix_path: str = "/home/b/b37974/htslib-1.9/tabix",
    keep_chr_prefix: bool = True,
) -> str:
    """
    将 PLINK 格式（.bed/.bim/.fam）转换为 **bgzip 压缩的 VCF（.vcf.gz）**，随后用参考基因组进行 **bcftools norm** 规范化，并最终使用 **tabix** 生成 `.tbi` 索引。

    变更点（与旧版相比）
    --------------------
    1) 所有输出文件保存到**当前工作目录**（不再使用 `./tmp`），并且**仅保留最终的** `*.vcf.gz` 与 `*.vcf.gz.tbi`；导出/重命名阶段的中间文件全部删除。
    2) 新增参数 `snps_only_just_acgt`（默认 False）：为 plink 导出添加 `--snps-only just-acgt` 仅提取 A/C/G/T 的 SNP（样本名仍为 IID via `--recode vcf-iid bgz`）。
    3) 最终索引改用 **tabix**（默认路径 `tabix_path="/home/b/b37974/htslib-1.9/tabix"`）；其余阶段若需索引仍使用 `bcftools index`。
    4) 加速&防止内存溢出：所有子进程**不捕获 STDOUT**（重定向到 `/dev/null`），仅保留 `stderr`；并尽可能为 `bcftools` 命令添加 `--threads`。
    5) **命名规则**：若仅导出 SNP（`snps_only_just_acgt=True`），最终文件名添加后缀 `.snp`；否则添加 `.all`（例如：`{basename}.snp.norm.vcf.gz` vs `{basename}.all.norm.vcf.gz`）。
    6) 新增参数 `keep_chr_prefix`（默认 True）：当为 False 时，删除 VCF 文件中 #CHROM 列的 'chr' 前缀。

    参数
    ----
    bed_prefix : str
        PLINK 输入前缀（不含后缀），例如 "/path/to/data"。
    plink_path : str
        plink 1.9 可执行文件路径；若在 $PATH 中可直接写 "plink"。
    bcftools_path : str
        bcftools 可执行文件路径。
    ref_fa : str
        参考基因组 FASTA 文件（需与当前坐标系一致，且建议已建 .fai 索引）。
    threads : int
        并行线程数（传递给 bcftools；plink 也支持 --threads）。
    rename_chr : str
        `bcftools annotate --rename-chrs` 映射文件路径（可选；存在则会在 `bcftools norm` 之前应用）。
    snps_only_just_acgt : bool
        True 时在 plink 导出时添加 `--snps-only just-acgt` 仅保留 A/C/G/T 的 SNP。
    tabix_path : str
        tabix 可执行文件路径（用于最终 `.tbi` 索引）。
    keep_chr_prefix : bool
        True 时保留染色体名称中的 'chr' 前缀；False 时删除 'chr' 前缀。

    返回
    ----
    str
        规范化后的 `.vcf.gz` 文件的绝对路径。
    """
    import os
    import shutil
    import subprocess

    def _run(cmd, desc: str):
        """运行外部命令：不捕获 STDOUT，捕获并在失败时回显 STDERR。"""
        proc = subprocess.run(
            cmd,
            stdout=subprocess.DEVNULL,  # 防止大输出占用内存
            stderr=subprocess.PIPE,
            text=True,
            check=False,
        )
        if proc.returncode != 0:
            raise RuntimeError(
                f"{desc} 失败\n命令: {' '.join(cmd)}\nSTDERR:\n" + (proc.stderr or "")
            )
        return proc

    # --- 校验可执行文件 ---
    if not (os.path.isfile(plink_path) or shutil.which(plink_path)):
        raise FileNotFoundError(f"找不到 plink：{plink_path}")
    if not (os.path.isfile(bcftools_path) or shutil.which(bcftools_path)):
        raise FileNotFoundError(f"找不到 bcftools：{bcftools_path}")
    if not (os.path.isfile(tabix_path) or shutil.which(tabix_path)):
        raise FileNotFoundError(f"找不到 tabix：{tabix_path}")
    if not os.path.exists(ref_fa) or os.path.getsize(ref_fa) == 0:
        raise FileNotFoundError(f"找不到参考基因组 FASTA：{ref_fa}")

    # --- 校验输入文件 ---
    req = [f"{bed_prefix}.bed", f"{bed_prefix}.bim", f"{bed_prefix}.fam"]
    missing = [p for p in req if not os.path.exists(p) or os.path.getsize(p) == 0]
    if missing:
        raise FileNotFoundError("缺少 PLINK 输入文件或文件为空：" + ", ".join(missing))

    # --- 线程参数 ---
    threads = int(threads) if isinstance(threads, (int, str)) else 6
    plink_threads = max(1, int(threads))

    # --- 输出前缀：工作目录 + 源文件名（basename） + 按是否仅 SNP 加后缀 ---
    cwd = os.getcwd()
    base_name = os.path.basename(bed_prefix)
    suffix_tag = ".snp" if snps_only_just_acgt else ".all"
    out_prefix = os.path.join(cwd, f"{base_name}{suffix_tag}")

    # 路径定义
    vcf_gz_path = f"{out_prefix}.vcf.gz"         # plink 直接导出
    tbi_path = f"{vcf_gz_path}.tbi"
    renamed_vcf = f"{out_prefix}.ren.vcf.gz"     # 染色体重命名后
    renamed_tbi = f"{renamed_vcf}.tbi"
    norm_vcf_gz_path = f"{out_prefix}.norm.vcf.gz"  # 规范化后
    norm_tbi_path = f"{norm_vcf_gz_path}.tbi"
    nochr_vcf_gz_path = f"{out_prefix}.nochr.norm.vcf.gz"  # 删除chr前缀后（最终保留）
    nochr_tbi_path = f"{nochr_vcf_gz_path}.tbi"

    # --- 1) 使用 plink 导出 bgzip VCF（样本名仅 IID；可选仅 SNP） ---
    export_cmd = [
        plink_path,
        "--bfile", bed_prefix,
        "--recode", "vcf-iid", "bgz",
        "--threads", str(plink_threads),
        "--out", out_prefix,
    ]
    if snps_only_just_acgt:
        export_cmd.extend(["--snps-only", "just-acgt"])  # 仅 A/C/G/T SNP

    print("[信息] 正在用 plink 导出 bgzip VCF（样本名为 IID）…")
    _run(export_cmd, "plink 导出")
    if not os.path.exists(vcf_gz_path) or os.path.getsize(vcf_gz_path) == 0:
        raise RuntimeError(f"未找到导出的 .vcf.gz 文件或文件大小为 0：{vcf_gz_path}")

    # 为导出的 .vcf.gz 建立临时索引（bcftools）
    index_cmd_1 = [bcftools_path, "index", "-t", "--threads", str(plink_threads), vcf_gz_path]
    print("[信息] 正在用 bcftools 为初始 .vcf.gz 建立 tbi 索引…")
    _run(index_cmd_1, "bcftools index（初始 VCF）")
    if not os.path.exists(tbi_path) or os.path.getsize(tbi_path) == 0:
        raise RuntimeError(f"未找到生成的初始 tbi 索引或文件大小为 0：{tbi_path}")

    # --- 2) 可选：染色体重命名（bcftools annotate --rename-chrs） ---
    use_renamed = False
    if rename_chr and isinstance(rename_chr, str) and rename_chr.strip() and os.path.exists(rename_chr):
        print(f"[信息] 检测到 rename_chr 映射文件，执行染色体重命名：{rename_chr}")
        rename_cmd = [
            bcftools_path, "annotate",
            "--rename-chrs", rename_chr,
            "--threads", str(plink_threads),
            "-Oz",
            "-o", renamed_vcf,
            vcf_gz_path,
        ]
        _run(rename_cmd, "bcftools annotate --rename-chrs")

        # 重命名后建立临时索引（bcftools）
        index_cmd_rename = [bcftools_path, "index", "-t", "--threads", str(plink_threads), renamed_vcf]
        print("[信息] 正在为重命名后的 VCF 建立 tbi 索引…")
        _run(index_cmd_rename, "bcftools index（重命名 VCF）")
        if not os.path.exists(renamed_tbi) or os.path.getsize(renamed_tbi) == 0:
            raise RuntimeError(f"未找到重命名后的 tbi 索引或文件大小为 0：{renamed_tbi}")
        use_renamed = True
    else:
        if rename_chr and isinstance(rename_chr, str) and rename_chr.strip():
            print(f"[警告] 指定的 rename_chr 文件不存在，跳过染色体重命名：{rename_chr}")

    # --- 3) bcftools norm 规范化（多等位拆分 + 参考等位基因校验/左对齐） ---
    input_for_norm = renamed_vcf if use_renamed else vcf_gz_path
    norm_cmd = [
        bcftools_path, "norm",
        "--multiallelics", "-any",
        "--fasta-ref", ref_fa,
        "--check-ref", "s",
        "--threads", str(plink_threads),
        input_for_norm,
        "-Oz",
        "-o", norm_vcf_gz_path,
    ]
    print("[信息] 正在用 bcftools norm 进行规范化（--multiallelics -any, --check-ref s）…")
    _run(norm_cmd, "bcftools norm")
    if not os.path.exists(norm_vcf_gz_path) or os.path.getsize(norm_vcf_gz_path) == 0:
        raise RuntimeError(f"未找到规范化后的 .vcf.gz 文件或文件大小为 0：{norm_vcf_gz_path}")

    # --- 4) 最终索引：使用 tabix 生成 .tbi ---
    print("[信息] 正在用 tabix 为规范化后的 .vcf.gz 建立 tbi 索引…")
    tabix_cmd = [tabix_path, "-f", "-p", "vcf", norm_vcf_gz_path]
    _run(tabix_cmd, "tabix index（规范化 VCF）")
    if not os.path.exists(norm_tbi_path) or os.path.getsize(norm_tbi_path) == 0:
        raise RuntimeError(f"未找到规范化后的 tbi 索引或文件大小为 0：{norm_tbi_path}")

    # --- 5) 可选：删除 chr 前缀 ---
    final_vcf_gz = norm_vcf_gz_path
    final_tbi = norm_tbi_path
    cleanup_files = [vcf_gz_path, tbi_path, renamed_vcf, renamed_tbi]
    
    if not keep_chr_prefix:
        print("[信息] 正在删除染色体名称中的 'chr' 前缀…")
        # 创建临时的染色体重命名文件（chr -> 无前缀）
        import tempfile
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.chr_remap.txt') as tmp_chr_file:
            # 写入染色体重命名映射：chr1 -> 1, chr2 -> 2, etc.
            for i in range(1, 23):  # 1-22
                tmp_chr_file.write(f"chr{i}\t{i}\n")
            tmp_chr_file.write("chrX\tX\n")
            tmp_chr_file.write("chrY\tY\n")
            tmp_chr_file.write("chrMT\tMT\n")
            tmp_chr_file.write("chrM\tM\n")
            tmp_chr_remap_path = tmp_chr_file.name
        
        try:
            # 使用 bcftools annotate --rename-chrs 删除 chr 前缀
            remove_chr_cmd = [
                bcftools_path, "annotate",
                "--rename-chrs", tmp_chr_remap_path,
                "--threads", str(plink_threads),
                "-Oz",
                "-o", nochr_vcf_gz_path,
                norm_vcf_gz_path,
            ]
            _run(remove_chr_cmd, "删除 chr 前缀")
            if not os.path.exists(nochr_vcf_gz_path) or os.path.getsize(nochr_vcf_gz_path) == 0:
                raise RuntimeError(f"未找到删除 chr 前缀后的 .vcf.gz 文件或文件大小为 0：{nochr_vcf_gz_path}")
        finally:
            # 清理临时文件
            try:
                os.remove(tmp_chr_remap_path)
            except Exception:
                pass
        
        # 为删除 chr 前缀后的文件建立索引
        print("[信息] 正在用 tabix 为删除 chr 前缀后的 .vcf.gz 建立 tbi 索引…")
        tabix_nochr_cmd = [tabix_path, "-f", "-p", "vcf", nochr_vcf_gz_path]
        _run(tabix_nochr_cmd, "tabix index（删除 chr 前缀 VCF）")
        if not os.path.exists(nochr_tbi_path) or os.path.getsize(nochr_tbi_path) == 0:
            raise RuntimeError(f"未找到删除 chr 前缀后的 tbi 索引或文件大小为 0：{nochr_tbi_path}")
        
        # 更新最终文件路径和清理列表
        final_vcf_gz = nochr_vcf_gz_path
        final_tbi = nochr_tbi_path
        cleanup_files.extend([norm_vcf_gz_path, norm_tbi_path])

    # --- 6) 清理中间文件（仅保留最终文件） ---
    for f in cleanup_files:
        try:
            if f and os.path.exists(f):
                os.remove(f)
        except Exception as e:
            print(f"[警告] 删除中间文件失败：{f}，错误：{e}")

    abs_final_vcf_gz = os.path.abspath(final_vcf_gz)
    print(f"[完成] 生成最终的压缩 VCF：{abs_final_vcf_gz}")
    print(f"[完成] 生成最终的索引：{os.path.abspath(final_tbi)}")
    return abs_final_vcf_gz



def reformat_pheno_covar(pheno_path: str, covar_path: str, out_dir: str = None, out_sep: str = " ", force: bool = True) -> tuple[str, str]: # type: ignore
    """
    重新格式化 phenotype 和 covariate 文件以适配 rvtest:
    1. 将所有列名转为小写
    2. 删除列名中的 '#' 符号
    3. 输出新的 pheno 文件和 covar 文件到指定目录（默认当前工作目录）
    4. 输出分隔符可配置，默认空格(" ")
    5. 当 force=True 时，在 phenotype 文件中的 iid 列后插入 fatid、matid 列（填充为0），然后插入 sex 列（从 covariate 文件获取）
    
    参数
    ----
    pheno_path : str
        原始 phenotype 文件路径（制表符分隔）
    covar_path : str
        原始 covariate 文件路径（制表符分隔）
    out_dir : str, 默认 None
        输出目录，若为 None 则为当前工作目录
    out_sep : str, 默认 " "
        输出文件的分隔符，默认使用空格。可根据需要设置为"\t"等其它分隔符。
    force : bool, 默认 True
        是否在 phenotype 文件中强制插入 fatid、matid、sex 列
    
    返回
    ----
    (str, str)
        新的 phenotype 文件路径, 新的 covariate 文件路径
    """
    import os
    import pandas as pd

    if out_dir is None:
        out_dir = os.getcwd()

    def _process_file(path: str, suffix: str) -> str:
        df = pd.read_csv(path, sep="\t", dtype=str)
        # 列名处理
        new_cols = [c.lower().replace("#", "") for c in df.columns]
        df.columns = new_cols
        out_path = os.path.join(out_dir, os.path.basename(path).replace(".txt", f".{suffix}.csv"))
        df.to_csv(out_path, sep=out_sep, index=False)
        return out_path

    # 处理 covariate 文件
    new_covar = _process_file(covar_path, "covar.reformat")
    
    # 处理 phenotype 文件
    if force:
        # 读取并处理 phenotype 文件
        pheno_df = pd.read_csv(pheno_path, sep="\t", dtype=str)
        new_cols = [c.lower().replace("#", "") for c in pheno_df.columns]
        pheno_df.columns = new_cols
        
        # 读取 covariate 文件以获取 sex 信息
        covar_df = pd.read_csv(covar_path, sep="\t", dtype=str)
        covar_cols = [c.lower().replace("#", "") for c in covar_df.columns]
        covar_df.columns = covar_cols
        
        # 找到 iid 列的位置
        if 'iid' not in pheno_df.columns:
            raise ValueError("phenotype 文件中未找到 'iid' 列")
        
        iid_pos = list(pheno_df.columns).index('iid')
        
        # 在 iid 列后插入 fatid 和 matid 列（填充为 "0"）
        pheno_df.insert(iid_pos + 1, 'fatid', '0')
        pheno_df.insert(iid_pos + 2, 'matid', '0')
        
        # 从 covariate 文件中获取 sex 信息并插入
        if 'sex' in covar_df.columns and 'fid' in covar_df.columns and 'iid' in covar_df.columns:
            # 创建用于匹配的字典
            sex_dict = {}
            for _, row in covar_df.iterrows():
                key = (str(row['fid']), str(row['iid']))
                sex_dict[key] = str(row['sex'])
            
            # 为 phenotype 文件添加 sex 列
            sex_values = []
            for _, row in pheno_df.iterrows():
                key = (str(row['fid']), str(row['iid']))
                sex_values.append(sex_dict.get(key, '0'))  # 如果找不到匹配，默认为 '0'
            
            pheno_df.insert(iid_pos + 3, 'sex', sex_values)
        else:
            print("[警告] covariate 文件中缺少必要的列（fid、iid、sex），sex 列将填充为 '0'")
            pheno_df.insert(iid_pos + 3, 'sex', '0')
        
        # 保存处理后的 phenotype 文件
        pheno_out_path = os.path.join(out_dir, os.path.basename(pheno_path).replace(".txt", ".pheno.reformat.csv"))
        pheno_df.to_csv(pheno_out_path, sep=out_sep, index=False)
        new_pheno = pheno_out_path
    else:
        # 使用原始处理方法
        new_pheno = _process_file(pheno_path, "pheno.reformat")

    print(f"[完成] 输出新的 phenotype 文件: {new_pheno}")
    print(f"[完成] 输出新的 covariate 文件: {new_covar}")

    return new_pheno, new_covar


def reformat_refflat_remove_chr(refflat_path: str, out_dir: str = None) -> str:  # type: ignore
    """
    整理 refFlat 文件，去掉第3列（染色体列）的 'chr' 前缀。
    优化版本：使用流式处理，逐行读取和写入，显著减少内存占用并提高处理速度。
    
    refFlat 文件通常为制表符分隔的文件，格式为：
    geneName, name, chrom, strand, txStart, txEnd, cdsStart, cdsEnd, exonCount, exonStarts, exonEnds
    
    参数
    ----
    refflat_path : str
        输入的 refFlat 文件路径（.txt.gz 格式）
    out_dir : str, 默认 None
        输出目录，若为 None 则为当前工作目录
        
    返回
    ----
    str
        处理后的 refFlat 文件路径（.txt.gz 格式）
    """
    import os
    import gzip
    
    if out_dir is None:
        out_dir = os.getcwd()
    
    # 检查输入文件是否存在
    if not os.path.exists(refflat_path):
        raise FileNotFoundError(f"refFlat 文件不存在: {refflat_path}")
    
    # 生成输出文件名
    base_name = os.path.basename(refflat_path)
    if base_name.endswith('.txt.gz'):
        out_name = base_name.replace('.txt.gz', '.nochr.txt.gz')
    elif base_name.endswith('.gz'):
        out_name = base_name.replace('.gz', '.nochr.gz')
    else:
        out_name = base_name + '.nochr.gz'
    
    out_path = os.path.join(out_dir, out_name)
    
    print(f"[信息] 正在处理 refFlat 文件: {refflat_path}")
    print(f"[信息] 输出文件: {out_path}")
    
    try:
        # 用于统计处理信息
        line_count = 0
        chroms_before = set()
        chroms_after = set()
        
        # 流式处理：逐行读取和写入
        with gzip.open(refflat_path, 'rt', encoding='utf-8') as infile, \
             gzip.open(out_path, 'wt', encoding='utf-8', compresslevel=6) as outfile:
            
            for line in infile:
                line = line.rstrip('\n\r')
                if not line:  # 跳过空行
                    continue
                
                fields = line.split('\t')
                
                # 检查列数（第一行时检查）
                if line_count == 0 and len(fields) < 3:
                    raise ValueError(f"refFlat 文件列数不足，需要至少3列，但只有 {len(fields)} 列")
                
                # 处理第3列（索引为2）的染色体信息
                if len(fields) > 2:
                    original_chrom = fields[2]
                    chroms_before.add(original_chrom)
                    
                    # 去掉 'chr' 前缀（使用简单字符串操作，比正则表达式快）
                    if original_chrom.startswith('chr'):
                        fields[2] = original_chrom[3:]  # 去掉前3个字符 'chr'
                    
                    chroms_after.add(fields[2])
                
                # 写入处理后的行
                outfile.write('\t'.join(fields) + '\n')
                line_count += 1
                
                # 每处理 100万行显示一次进度
                if line_count % 1000000 == 0:
                    print(f"[进度] 已处理 {line_count:,} 行")
        
        print(f"[信息] 总共处理 {line_count:,} 行数据")
        
        # 显示染色体统计信息（限制显示数量）
        chroms_before_sorted = sorted(list(chroms_before))
        chroms_after_sorted = sorted(list(chroms_after))
        
        print(f"[信息] 处理前的染色体 ({len(chroms_before_sorted)} 种): {chroms_before_sorted[:10]}{'...' if len(chroms_before_sorted) > 10 else ''}")
        print(f"[信息] 处理后的染色体 ({len(chroms_after_sorted)} 种): {chroms_after_sorted[:10]}{'...' if len(chroms_after_sorted) > 10 else ''}")
        
        print(f"[完成] 已生成去除 chr 前缀的 refFlat 文件: {out_path}")
        return os.path.abspath(out_path)
        
    except Exception as e:
        # 如果处理失败，清理可能的部分输出文件
        try:
            if os.path.exists(out_path):
                os.remove(out_path)
        except Exception:
            pass
        raise RuntimeError(f"处理 refFlat 文件时出错: {str(e)}")



