def plink_to_vcf_raw(
    bed_prefix: str,
    plink_path: str = "/home/b/b37974/plink",
    bcftools_path: str = "/home/b/b37974/bcftools/bcftools",
    ref_fa: str = "/LARGE0/gr10478/b37974/Pulmonary_Hypertension/nagasaki_pipeline/data/hs38DH.fa",
    threads: int = 6,
    delete_tmp: bool = False,
    rename_chr: str = "/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/info/rename_chr.txt",
) -> str:
    """
    将 PLINK 格式（.bed/.bim/.fam）转换为 **bgzip 压缩的 VCF（.vcf.gz）**，随后用参考基因组进行 **bcftools norm** 规范化，并建立 **tbi** 索引。
    
    变更点（与旧版相比）
    --------------------
    1) 使用 plink 1.9：`--recode vcf-iid bgz` 导出 `.vcf.gz`（样本名只使用 IID）。
    2) 导出后使用 `bcftools index -t` 为初始 `.vcf.gz` 建立索引。
    3) 使用 `bcftools norm`（`--multiallelics -any --fasta-ref {ref_fa} --check-ref s`）对 VCF 进行规范化，输出新的 `.vcf.gz`。
    4) 对规范化后的 `.vcf.gz` 再次建立索引（`.tbi`），并 **删除导出阶段的旧 `.vcf.gz` 与 `.tbi`**。
    5) 输出文件默认写入当前工作目录下的 `./tmp` 文件夹（若不存在则自动创建）；`delete_tmp=True` 仍会在最后删除该文件夹（包含本次输出）。
    
    参数
    ----
    bed_prefix : str
        PLINK 输入前缀（不含后缀），例如 "/path/to/data"。
    plink_path : str, 默认 "/home/b/b37974/plink"
        plink 1.9 可执行文件路径；若在 $PATH 中可直接写 "plink"。
    bcftools_path : str, 默认 "/home/b/b37974/bcftools/bcftools"
        bcftools 可执行文件路径。
    ref_fa : str, 默认 "/LARGE0/gr10478/b37974/Pulmonary_Hypertension/nagasaki_pipeline/data/hs38DH.fa"
        参考基因组 FASTA 文件（需与当前坐标系一致，且建议已建 .fai 索引）。
    threads : int, 默认 6
        并行线程数（传递给 bcftools；plink 也支持 --threads）。
    delete_tmp : bool, 默认 False
        是否在完成后删除 `./tmp` 文件夹（**会删除本次输出文件**）。
    rename_chr : str, 默认 "/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/info/rename_chr.txt"
        bcftools `--rename-chrs` 映射文件路径（可选；如果存在则在 `bcftools norm` 之前应用）。
    
    返回
    ----
    str
        规范化后的 `.vcf.gz` 文件的绝对路径（若 `delete_tmp=True` 则返回的路径随后将失效）。
    """
    import os
    import shutil
    import subprocess

    # --- 校验可执行文件 ---
    if not (os.path.isfile(plink_path) or shutil.which(plink_path)):
        raise FileNotFoundError(f"找不到 plink：{plink_path}")
    if not (os.path.isfile(bcftools_path) or shutil.which(bcftools_path)):
        raise FileNotFoundError(f"找不到 bcftools：{bcftools_path}")
    if not os.path.exists(ref_fa) or os.path.getsize(ref_fa) == 0:
        raise FileNotFoundError(f"找不到参考基因组 FASTA：{ref_fa}")

    # --- 校验输入文件 ---
    req = [f"{bed_prefix}.bed", f"{bed_prefix}.bim", f"{bed_prefix}.fam"]
    missing = [p for p in req if not os.path.exists(p) or os.path.getsize(p) == 0]
    if missing:
        raise FileNotFoundError("缺少 PLINK 输入文件或文件为空：" + ", ".join(missing))

    # --- 规范化线程参数 ---
    threads = int(threads) if isinstance(threads, (int, str)) else 6
    plink_threads = max(1, int(threads))

    # --- 输出目录：当前工作目录下的 ./tmp ---
    cwd = os.getcwd()
    out_dir = os.path.join(cwd, "tmp")
    os.makedirs(out_dir, exist_ok=True)

    # 以 bed_prefix 的 **basename** 作为输出前缀名，避免把文件写回到数据源目录
    out_prefix = os.path.join(out_dir, os.path.basename(bed_prefix))

    # 初始导出 VCF 路径（由 plink 直接生成）
    vcf_gz_path = f"{out_prefix}.vcf.gz"
    tbi_path = f"{vcf_gz_path}.tbi"

    # 规范化后的最终输出路径
    norm_vcf_gz_path = f"{out_prefix}.norm.vcf.gz"
    norm_tbi_path = f"{norm_vcf_gz_path}.tbi"

    # 用于记录所有中间文件，方便最后删除
    tmp_files = [vcf_gz_path, tbi_path]

    # --- 1) 使用 plink 导出 bgzip VCF（样本名仅 IID） ---
    export_cmd = [
        plink_path,
        "--bfile", bed_prefix,
        "--recode", "vcf-iid", "bgz",
        "--threads", str(plink_threads),
        "--out", out_prefix,
    ]

    print("[信息] 正在用 plink (1.9) 导出 bgzip 压缩 VCF（.vcf.gz，样本名为 IID）…")
    proc = subprocess.run(export_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if proc.returncode != 0:
        raise RuntimeError(
            "plink 导出失败\n命令: " + " ".join(export_cmd) +
            "\nSTDOUT:\n" + (proc.stdout or "") +
            "\nSTDERR:\n" + (proc.stderr or "")
        )
    if not os.path.exists(vcf_gz_path) or os.path.getsize(vcf_gz_path) == 0:
        raise RuntimeError(f"未找到导出的 .vcf.gz 文件或文件大小为 0：{vcf_gz_path}")

    # --- 2) 为导出的 .vcf.gz 建立 tbi 索引 ---
    index_cmd_1 = [bcftools_path, "index", "-t", "--threads", str(plink_threads), vcf_gz_path]
    print("[信息] 正在用 bcftools 为初始 .vcf.gz 建立 tbi 索引…")
    proc_index1 = subprocess.run(index_cmd_1, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if proc_index1.returncode != 0:
        raise RuntimeError(
            "bcftools index（初始 VCF）失败\n命令: " + " ".join(index_cmd_1) +
            "\nSTDOUT:\n" + (proc_index1.stdout or "") +
            "\nSTDERR:\n" + (proc_index1.stderr or "")
        )
    if not os.path.exists(tbi_path) or os.path.getsize(tbi_path) == 0:
        raise RuntimeError(f"未找到生成的初始 tbi 索引或文件大小为 0：{tbi_path}")

    # --- 新增步骤：根据 rename_chr 文件进行染色体重命名 ---
    if rename_chr and isinstance(rename_chr, str) and rename_chr.strip():
        if os.path.exists(rename_chr):
            renamed_vcf = f"{out_prefix}.ren.vcf.gz"
            renamed_tbi = f"{renamed_vcf}.tbi"
            rename_cmd = [
                bcftools_path, "annotate",
                "--rename-chrs", rename_chr,
                "--threads", str(plink_threads),
                "-Oz",
                "-o", renamed_vcf,
                vcf_gz_path,
            ]
            print(f"[信息] 检测到 rename_chr 文件，正在用 bcftools annotate --rename-chrs 进行染色体重命名，映射文件：{rename_chr}")
            proc_rename = subprocess.run(rename_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            if proc_rename.returncode != 0:
                raise RuntimeError(
                    "bcftools annotate --rename-chrs 失败\n命令: " + " ".join(rename_cmd) +
                    "\nSTDOUT:\n" + (proc_rename.stdout or "") +
                    "\nSTDERR:\n" + (proc_rename.stderr or "")
                )
            # 建立重命名后 VCF 索引
            index_cmd_rename = [bcftools_path, "index", "-t", "--threads", str(plink_threads), renamed_vcf]
            print("[信息] 正在为重命名后的 VCF 建立 tbi 索引…")
            proc_index_rename = subprocess.run(index_cmd_rename, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            if proc_index_rename.returncode != 0:
                raise RuntimeError(
                    "bcftools index（重命名 VCF）失败\n命令: " + " ".join(index_cmd_rename) +
                    "\nSTDOUT:\n" + (proc_index_rename.stdout or "") +
                    "\nSTDERR:\n" + (proc_index_rename.stderr or "")
                )
            if not os.path.exists(renamed_tbi) or os.path.getsize(renamed_tbi) == 0:
                raise RuntimeError(f"未找到重命名后的 tbi 索引或文件大小为 0：{renamed_tbi}")

            # 更新路径变量，后续使用重命名后的 VCF 进行规范化
            vcf_gz_path = renamed_vcf
            tbi_path = renamed_tbi

            # 记录重命名中间文件，方便最后删除
            tmp_files.extend([renamed_vcf, renamed_tbi])
            print(f"[信息] 已应用染色体重命名，使用映射文件：{rename_chr}")
        else:
            print(f"[警告] 指定的 rename_chr 文件不存在，跳过染色体重命名步骤，继续使用自动协调逻辑：{rename_chr}")

    # --- 3) bcftools norm 规范化（多等位拆分 + 参考等位基因校验/左对齐） ---
    norm_cmd = [
        bcftools_path, "norm",
        "--multiallelics", "-any",
        "--fasta-ref", ref_fa,
        "--check-ref", "s",
        vcf_gz_path,
        "-Oz",
        "-o", norm_vcf_gz_path,
    ]
    print("[信息] 正在用 bcftools norm 进行规范化（--multiallelics -any, --check-ref s）…")
    proc_norm = subprocess.run(norm_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if proc_norm.returncode != 0:
        raise RuntimeError(
            "bcftools norm 失败\n命令: " + " ".join(norm_cmd) +
            "\nSTDOUT:\n" + (proc_norm.stdout or "") +
            "\nSTDERR:\n" + (proc_norm.stderr or "")
        )
    if not os.path.exists(norm_vcf_gz_path) or os.path.getsize(norm_vcf_gz_path) == 0:
        raise RuntimeError(f"未找到规范化后的 .vcf.gz 文件或文件大小为 0：{norm_vcf_gz_path}")

    # --- 规范化后建立 tbi 索引 ---
    index_cmd_2 = [bcftools_path, "index", "-t", "--threads", str(plink_threads), norm_vcf_gz_path]
    print("[信息] 正在用 bcftools 为规范化后的 .vcf.gz 建立 tbi 索引…")
    proc_index2 = subprocess.run(index_cmd_2, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if proc_index2.returncode != 0:
        raise RuntimeError(
            "bcftools index（规范化 VCF）失败\n命令: " + " ".join(index_cmd_2) +
            "\nSTDOUT:\n" + (proc_index2.stdout or "") +
            "\nSTDERR:\n" + (proc_index2.stderr or "")
        )
    if not os.path.exists(norm_tbi_path) or os.path.getsize(norm_tbi_path) == 0:
        raise RuntimeError(f"未找到规范化后的 tbi 索引或文件大小为 0：{norm_tbi_path}")

    # --- 4) 删除所有中间文件（导出阶段及重命名阶段） ---
    for f in tmp_files:
        try:
            if os.path.exists(f):
                os.remove(f)
        except Exception as e:
            print(f"[警告] 删除中间文件失败：{f}，错误：{e}")

    abs_norm_vcf_gz = os.path.abspath(norm_vcf_gz_path)
    print(f"[完成] 生成规范化后的压缩 VCF：{abs_norm_vcf_gz}")
    print(f"[完成] 生成规范化后的索引：{os.path.abspath(norm_tbi_path)}")

    # 若需要，删除 ./tmp 文件夹（将一并删除输出文件）
    if delete_tmp:
        print("[警告] delete_tmp=True：将删除 ./tmp 文件夹以及其中的输出文件。")
        try:
            shutil.rmtree(out_dir)
            print("[信息] 已删除 ./tmp 文件夹。")
        except Exception as e:
            print(f"[警告] 删除 ./tmp 失败：{e}")

    return abs_norm_vcf_gz
