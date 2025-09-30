# 模块说明：提供将 PLINK 二进制（.bed/.bim/.fam）导出为 bgzip 压缩的 VCF（.vcf.gz）并在 ./tmp 目录下建立 .tbi 索引的工具函数。

import os
import shutil
import subprocess

def plink_to_vcf_raw(
    bed_prefix: str,
    plink2_path: str = "/home/b/b37974/plink2",
    bcftools_path: str = "/home/b/b37974/bcftools/bcftools",
    threads: int = 6,
    delete_tmp: bool = False,
) -> str:
    """
    将 PLINK 格式（.bed/.bim/.fam）转换为 **bgzip 压缩的 VCF（.vcf.gz）**，并使用 bcftools 为该 .vcf.gz 建立 **tbi** 索引。

    变更点（与旧版相比）
    --------------------
    1) 直接由 plink2 导出 `.vcf.gz`（`--export vcf bgz`），不再生成未压缩 `.vcf`。
    2) 输出文件默认写入当前工作目录下的 `./tmp` 文件夹（若不存在则自动创建）。
    3) 新增参数 `delete_tmp` 控制是否在流程完成后删除 `./tmp` 文件夹（默认 False）。
       **注意**：若设为 True，将在函数结束时删除 `./tmp`，届时本次生成的 `.vcf.gz` 与 `.tbi` 也会被删除。

    参数
    ----
    bed_prefix : str
        PLINK 输入前缀（不含后缀），例如 "/path/to/data"。
    plink2_path : str, 默认 "/home/b/b37974/plink2"
        plink2 可执行文件路径；若在 $PATH 中可直接写 "plink2"。
    bcftools_path : str, 默认 "/home/b/b37974/bcftools/bcftools"
        bcftools 可执行文件路径。
    threads : int, 默认 6
        并行线程数。
    delete_tmp : bool, 默认 False
        是否在完成后删除 `./tmp` 文件夹（**会删除本次输出文件**）。

    返回
    ----
    str
        生成的 `.vcf.gz` 文件的绝对路径（若 `delete_tmp=True` 则返回的路径随后将失效）。
    """
    import os
    import shutil
    import subprocess

    # --- 校验可执行文件 ---
    if not (os.path.isfile(plink2_path) or shutil.which(plink2_path)):
        raise FileNotFoundError(f"找不到 plink2：{plink2_path}")
    if not (os.path.isfile(bcftools_path) or shutil.which(bcftools_path)):
        raise FileNotFoundError(f"找不到 bcftools：{bcftools_path}")

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

    vcf_gz_path = f"{out_prefix}.vcf.gz"
    tbi_path = f"{vcf_gz_path}.tbi"

    # --- 组装 plink2 命令（直接导出 bgzip 压缩 VCF） ---
    export_cmd = [
        plink2_path,
        "--bfile", bed_prefix,
        "--export", "vcf", "bgz",
        "--threads", str(plink_threads),
        "--out", out_prefix,
    ]

    print("[信息] 正在用 plink2 导出 bgzip 压缩 VCF（.vcf.gz）…")
    proc = subprocess.run(export_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if proc.returncode != 0:
        raise RuntimeError(
            "plink2 导出失败\n命令: " + " ".join(export_cmd) +
            "\nSTDOUT:\n" + (proc.stdout or "") +
            "\nSTDERR:\n" + (proc.stderr or "")
        )

    if not os.path.exists(vcf_gz_path) or os.path.getsize(vcf_gz_path) == 0:
        raise RuntimeError(f"未找到导出的 .vcf.gz 文件或文件大小为 0：{vcf_gz_path}")

    # --- 为 .vcf.gz 建 tbi 索引 ---
    index_cmd = [bcftools_path, "index", "-t", "--threads", str(plink_threads), vcf_gz_path]
    print("[信息] 正在用 bcftools 为 .vcf.gz 建立 tbi 索引…")
    proc_index = subprocess.run(index_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if proc_index.returncode != 0:
        raise RuntimeError(
            "bcftools index 失败\n命令: " + " ".join(index_cmd) +
            "\nSTDOUT:\n" + (proc_index.stdout or "") +
            "\nSTDERR:\n" + (proc_index.stderr or "")
        )

    if not os.path.exists(tbi_path) or os.path.getsize(tbi_path) == 0:
        raise RuntimeError(f"未找到生成的 tbi 索引或文件大小为 0：{tbi_path}")

    abs_vcf_gz = os.path.abspath(vcf_gz_path)
    print(f"[完成] 生成压缩 VCF：{abs_vcf_gz}")
    print(f"[完成] 生成索引：{os.path.abspath(tbi_path)}")

    # 若需要，删除 ./tmp 文件夹（将一并删除输出文件）
    if delete_tmp:
        print("[警告] delete_tmp=True：将删除 ./tmp 文件夹以及其中的输出文件。")
        try:
            shutil.rmtree(out_dir)
            print("[信息] 已删除 ./tmp 文件夹。")
        except Exception as e:
            print(f"[警告] 删除 ./tmp 失败：{e}")

    return abs_vcf_gz

# new function to norm input vcf.gz
