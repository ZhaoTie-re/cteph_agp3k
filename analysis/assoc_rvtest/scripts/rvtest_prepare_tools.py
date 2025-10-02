def plink_to_vcf_raw(
    bed_prefix: str,
    plink_path: str = "/home/b/b37974/plink",
    bcftools_path: str = "/home/b/b37974/bcftools/bcftools",
    ref_fa: str = "/LARGE0/gr10478/b37974/Pulmonary_Hypertension/nagasaki_pipeline/data/hs38DH.fa",
    threads: int = 6,
    rename_chr: str = "/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/info/rename_chr.txt",
    snps_only_just_acgt: bool = False,
    tabix_path: str = "/home/b/b37974/htslib-1.9/tabix",
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
    norm_vcf_gz_path = f"{out_prefix}.norm.vcf.gz"  # 规范化后（最终保留）
    norm_tbi_path = f"{norm_vcf_gz_path}.tbi"

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
    _run(tabix_cmd, "tabix index（最终 VCF）")
    if not os.path.exists(norm_tbi_path) or os.path.getsize(norm_tbi_path) == 0:
        raise RuntimeError(f"未找到规范化后的 tbi 索引或文件大小为 0：{norm_tbi_path}")

    # --- 5) 清理中间文件（仅保留最终 *.norm.vcf.gz 及其 .tbi） ---
    for f in [vcf_gz_path, tbi_path, renamed_vcf, renamed_tbi]:
        try:
            if f and os.path.exists(f):
                os.remove(f)
        except Exception as e:
            print(f"[警告] 删除中间文件失败：{f}，错误：{e}")

    abs_norm_vcf_gz = os.path.abspath(norm_vcf_gz_path)
    print(f"[完成] 生成规范化后的压缩 VCF：{abs_norm_vcf_gz}")
    print(f"[完成] 生成规范化后的索引：{os.path.abspath(norm_tbi_path)}")
    return abs_norm_vcf_gz


# new funtion to reformat pheno & covar file for rvtest
