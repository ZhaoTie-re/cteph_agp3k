"""
LD 比较分析工具集（tommo_ld_tools.py）

本工具集用于针对关注变体（focus variants），从日本人群参考数据库 ToMMo 中提取其 LD（连锁不平衡）信息，
并与用户提供的基因型数据中的 LD 计算结果进行对比分析。主要功能包括：

1. 提取 ToMMo LD 数据：
   - 通过 tabix 在 ToMMo 提供的 .tsv.gz 文件中查询指定的 focus loci；
   - 构建 ToMMo 中每个 focus 变体与其关联变体的 r² 结果表；
   - 生成匹配日志信息和 pickle 格式的 LD 字典。

2. 基于 Plink2 计算用户样本中的 LD：
   - 使用样本的 bed/bim/fam 文件，按 case/control 分组；
   - 对每个 focus 变体，与 ToMMo 中对应的 variation2_id 进行 LD 计算；
   - 识别并标注未出现在 BIM 文件中的变体、PLINK 输出缺失变体等情况；
   - 结果以 pickle 格式保存，包括 case/control 各自的 LD 表和日志信息。

3. ToMMo 与用户 LD 结果合并：
   - 对每个 focus 变体，根据 variation2_id 将 ToMMo 的 r² 值合并到用户 LD 结果中；
   - 插入新列 TOMMO_R2（在 UNPHASED_R2 列前）；
   - 标注未匹配记录或说明备注。

4. 可视化输出 PDF：
   - 对每个 focus 变体绘制散点图，x轴为 ToMMo R²，y轴为用户计算的 R²；
   - 添加斜率为1的参考线；
   - 标注 ToMMo 可用配对数 与 用户中可计算配对数；
   - 若无 ToMMo 数据或用户缺失所有关联变体，图中显示注释说明。

适用场景：
- 验证公开数据库（如 ToMMo）中的 LD 模式是否适用于用户自身的群体；
- 对比平台、群体、深度等对 LD 模式的影响；
- 制作适用于报告和论文的 LD 对比图。

使用建议：
- 建议在高质量的 WGS 数据或深度充足的 Array 数据基础上运行；
- 推荐使用高线程数并行运行 LD 计算（默认线程数为 4）；
- 可根据需要调整 LD 窗口大小（kb）和最小 r² 阈值。

作者：ZHAO TIE
"""
import os
import gzip
import pickle
import warnings
import subprocess
import pandas as pd
from collections import defaultdict
import concurrent.futures
from concurrent.futures import ProcessPoolExecutor
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

def run_subprocess(cmd, **kwargs):
    subprocess.run(cmd, **kwargs)

def extract_ld_from_tommo_for_focus_loci(
    focus_loci_path: str,
    tommo_ld_dir: str,
    output_chr_prefix: str = 'tommo-54kjpn-20230828-GRCh38-autosome-'
) -> tuple[dict[str, pd.DataFrame], pd.DataFrame]:
    """
    从 TOMMO LD 文件中提取关注变体的 LD 信息。

    参数:
    - focus_loci_path (str): 包含需要关注的变体的文件路径，每行格式为 chr:pos:ref:alt。
    - tommo_ld_dir (str): 包含 TOMMO LD 按染色体分割的 tsv.gz 和 .tbi 文件的目录。
    - output_chr_prefix (str): TOMMO LD 文件前缀（默认适用于 GRCh38 数据）。

    返回:
    - result_dict (dict): 以变体 ID 为键（chr:pos:ref:alt），对应匹配到的 LD 数据（DataFrame）为值。
    - log_df (DataFrame): 包含每个变体的匹配信息日志（是否命中、匹配数、是否异常等）。
    """

    # Step 1. 读取 focus loci 文件，并按染色体分组，准备从 TOMMO LD 数据中提取信息
    chr_locus_dict = defaultdict(list)
    with open(focus_loci_path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            chrom, pos, ref, alt = line.split(":")
            chr_locus_dict[chrom].append((pos, ref, alt, line))  # 保存原始 ID

    result = {}          # 存储最终 TOMMO LD 匹配结果
    log_records = []     # 存储日志信息，便于后续展示与追踪

    # Step 2. 遍历每条染色体下的变体，从 TOMMO LD 文件中查询匹配记录
    for chrom, loci in chr_locus_dict.items():
        ld_file = os.path.join(tommo_ld_dir, f"{output_chr_prefix}{chrom}-plink-r2.tsv.gz")
        if not os.path.exists(ld_file):
            print(f"[警告] 染色体 {chrom} 的 TOMMO LD 文件未找到: {ld_file}")
            continue

        for pos, ref, alt, full_id in loci:
            region = f"{chrom}:{pos}-{pos}"  # tabix 查询范围
            matched_rows = []

            try:
                cmd = ["tabix", ld_file, region]
                output = subprocess.check_output(cmd, text=True).splitlines()

                for row in output:
                    fields = row.strip().split("\t")
                    # 检查是否为该变体本身（可能是变体1或变体2）
                    if (fields[1] == chrom and fields[2] == pos and
                        fields[3] == ref and fields[4] == alt):
                        # skip if variation1_id == variation2_id
                        if fields[0] == fields[6]:
                            continue
                        matched_rows.append(fields)
                    elif (fields[7] == chrom and fields[8] == pos and
                          fields[9] == ref and fields[10] == alt):
                        # skip if variation1_id == variation2_id
                        if fields[0] == fields[6]:
                            continue
                        matched_rows.append(fields)

                if matched_rows:
                    df = pd.DataFrame(matched_rows, columns=[
                        "variation1_id", "variation1_chromosome", "variation1_position",
                        "variation1_reference", "variation1_alternative", "variation1_maf",
                        "variation2_id", "variation2_chromosome", "variation2_position",
                        "variation2_reference", "variation2_alternative", "variation2_maf", "r2"
                    ])
                    # 更新 variation1_id 和 variation2_id 列（使用新变量）
                    v1_chrom = df["variation1_chromosome"]
                    v1_pos = df["variation1_position"]
                    v1_ref = df["variation1_reference"]
                    v1_alt = df["variation1_alternative"]
                    v2_chrom = df["variation2_chromosome"]
                    v2_pos = df["variation2_position"]
                    v2_ref = df["variation2_reference"]
                    v2_alt = df["variation2_alternative"]

                    df["variation1_id"] = [f"{c}:{p}:{r}:{a}" for c, p, r, a in zip(v1_chrom, v1_pos, v1_ref, v1_alt)]
                    df["variation2_id"] = [f"{c}:{p}:{r}:{a}" for c, p, r, a in zip(v2_chrom, v2_pos, v2_ref, v2_alt)]

                    result[full_id] = df
                else:
                    # 构造空 DataFrame（无匹配记录时）
                    result[full_id] = pd.DataFrame(columns=[
                        "variation1_id", "variation1_chromosome", "variation1_position",
                        "variation1_reference", "variation1_alternative", "variation1_maf",
                        "variation2_id", "variation2_chromosome", "variation2_position",
                        "variation2_reference", "variation2_alternative", "variation2_maf", "r2"
                    ])

                # 记录日志
                log_records.append({
                    "ID": full_id,
                    "CHR": chrom,
                    "POS": pos,
                    "Match_Found": len(matched_rows) > 0,
                    "Num_Matches": len(matched_rows)
                })

            except subprocess.CalledProcessError as e:
                # tabix 报错处理
                print(f"[错误] tabix 查询失败：{region}\n{e}")
                result[full_id] = pd.DataFrame()
                log_records.append({
                    "ID": full_id,
                    "CHR": chrom,
                    "POS": pos,
                    "Match_Found": False,
                    "Num_Matches": 0,
                    "Error": str(e)
                })

    # 将 result 字典保存为 pickle 文件
    output_pkl_path = os.path.abspath("tommo_ld_dict.pkl")
    with open(output_pkl_path, "wb") as f:
        pickle.dump(result, f)

    return output_pkl_path, pd.DataFrame(log_records) # type: ignore



def calculate_ld_for_focus_loci_by_group(
    focus_loci_path: str,
    bed_prefix: str,
    case_prefix: str,
    plink2_path: str = "/home/b/b37974/plink2_alpha6/plink2",
    ld_window_r2: float = 0.2,
    ld_window_kb: int = 1000,
    threads: int = 4,
    strip_chr_prefix: bool = False
) -> tuple[dict[str, dict[str, pd.DataFrame]], pd.DataFrame]:
    """
    使用 plink2 分别计算 case/control 中关注变体的 LD。

    参数:
    - focus_loci_path (str): 文件路径，每行一个关注变体（格式：chr:pos:ref:alt）
    - bed_prefix (str): plink 数据前缀，包含 .bed/.bim/.fam
    - case_prefix (str): 用于区分 case/control 的个体前缀，例如 PHOM 表示 case，其他为 control
    - plink2_path (str): plink2 执行路径，默认 /home/b/b37974/plink2
    - ld_window_r2 (float): 最小 r2 值（默认 0.2）
    - ld_window_kb (int): LD 窗口大小（kb，默认 1000）
    - threads (int): plink2 使用的线程数（默认 4）
    - strip_chr_prefix (bool): 是否移除 locus 中的 chr 前缀（默认 False）

    返回:
    - tuple: (result_dict, log_df)
      每个变体对应一个子字典 {'case': case_ld_df, 'control': ctrl_ld_df}，以及日志 DataFrame
    """
    import tempfile

    log_records = []

    # Step 1. 读取关注变体列表
    with open(focus_loci_path) as f:
        loci = [line.strip() for line in f if line.strip() and not line.startswith("#")]

    result = {}

    with tempfile.TemporaryDirectory() as tmpdir:
        case_plink = os.path.join(tmpdir, "case_data")
        ctrl_plink = os.path.join(tmpdir, "ctrl_data")

        fam_path = f"{bed_prefix}.fam"
        case_keep_path = os.path.join(tmpdir, "case.keep")
        ctrl_remove_path = os.path.join(tmpdir, "ctrl.remove")

        # 读取 .fam 文件，提取 sample IDs
        case_samples = []
        ctrl_samples = []
        with open(fam_path) as fam_file:
            for line in fam_file:
                fields = line.strip().split()
                if not fields:
                    continue
                if fields[0].startswith(case_prefix):
                    case_samples.append((fields[0], fields[1]))
                else:
                    ctrl_samples.append((fields[0], fields[1]))

        # 写入 case.keep 和 ctrl.remove 文件（plink keep/remove: FID IID）
        with open(case_keep_path, "w") as f:
            for fid, iid in case_samples:
                f.write(f"{fid} {iid}\n")
        with open(ctrl_remove_path, "w") as f:
            for fid, iid in case_samples:
                f.write(f"{fid} {iid}\n")

        # 使用线程池并行提取 case 和 control 的独立 plink 文件
        cmds = [
            [
                plink2_path, "--bfile", bed_prefix,
                "--keep", case_keep_path,
                "--make-bed", "--out", case_plink,
                "--threads", str(threads),
            ],
            [
                plink2_path, "--bfile", bed_prefix,
                "--remove", ctrl_remove_path,
                "--make-bed", "--out", ctrl_plink,
                "--threads", str(threads),
            ]
        ]

        with concurrent.futures.ThreadPoolExecutor(max_workers=2) as executor:
            futures = []
            for cmd in cmds:
                group = 'case' if '--keep' in cmd else 'ctrl'
                log_path = os.path.abspath(f"{group}_plink2.command.out")
                futures.append(
                    executor.submit(
                        run_subprocess, cmd, check=True,
                        stdout=open(log_path, "a"),
                        stderr=subprocess.STDOUT
                    )
                )
            for future in concurrent.futures.as_completed(futures):
                future.result()

        for locus in loci:
            locus_for_ld_snp = locus
            if strip_chr_prefix and locus.startswith("chr"):
                locus_for_ld_snp = locus[3:]

            chrom, pos, ref, alt = locus.split(":")
            snp_id = locus_for_ld_snp

            out_case = os.path.join(tmpdir, "case_ld")
            out_ctrl = os.path.join(tmpdir, "ctrl_ld")

            cmd_case = [
                plink2_path,
                "--bfile", case_plink,
                "--r2-unphased",
                "--ld-snp", snp_id,
                "--ld-window-kb", str(ld_window_kb),
                "--ld-window-r2", str(ld_window_r2),
                "--out", out_case,
                "--threads", str(threads),
            ]
            cmd_ctrl = [
                plink2_path,
                "--bfile", ctrl_plink,
                "--r2-unphased",
                "--ld-snp", snp_id,
                "--ld-window-kb", str(ld_window_kb),
                "--ld-window-r2", str(ld_window_r2),
                "--out", out_ctrl,
                "--threads", str(threads),
            ]

            try:
                log_path_case = os.path.abspath("case_plink2.command.out")
                with open(log_path_case, "a") as log_f:
                    subprocess.run(cmd_case, stdout=log_f, stderr=log_f, check=True)
                df_case = pd.read_csv(f"{out_case}.vcor", sep='\s+')
                df_case = df_case[df_case['ID_A'] != df_case['ID_B']]
                log_records.append({
                    "ID": locus,
                    "CHR": chrom,
                    "POS": pos,
                    "Group": "case",
                    "Match_Found": len(df_case) > 0,
                    "Num_Matches": len(df_case),
                    "CMD_Status": "Success"
                })
            except subprocess.CalledProcessError as e:
                df_case = pd.DataFrame()
                # 检查 SNP 是否存在于 .bim 文件中
                bim_path = f"{case_plink}.bim"
                with open(bim_path) as bim_f:
                    bim_ids = set(line.split()[1] for line in bim_f)
                if snp_id not in bim_ids:
                    status = "Failed: SNP not found in BIM file"
                else:
                    status = f"Failed: {e}"
                log_records.append({
                    "ID": locus,
                    "CHR": chrom,
                    "POS": pos,
                    "Group": "case",
                    "Match_Found": False,
                    "Num_Matches": 0,
                    "CMD_Status": status
                })

            try:
                log_path_ctrl = os.path.abspath("ctrl_plink2.command.out")
                with open(log_path_ctrl, "a") as log_f:
                    subprocess.run(cmd_ctrl, stdout=log_f, stderr=log_f, check=True)
                df_ctrl = pd.read_csv(f"{out_ctrl}.vcor", sep='\s+')
                df_ctrl = df_ctrl[df_ctrl['ID_A'] != df_ctrl['ID_B']]
                log_records.append({
                    "ID": locus,
                    "CHR": chrom,
                    "POS": pos,
                    "Group": "control",
                    "Match_Found": len(df_ctrl) > 0,
                    "Num_Matches": len(df_ctrl),
                    "CMD_Status": "Success"
                })
            except subprocess.CalledProcessError as e:
                df_ctrl = pd.DataFrame()
                bim_path = f"{ctrl_plink}.bim"
                with open(bim_path) as bim_f:
                    bim_ids = set(line.split()[1] for line in bim_f)
                if snp_id not in bim_ids:
                    status = "Failed: SNP not found in BIM file"
                else:
                    status = f"Failed: {e}"
                log_records.append({
                    "ID": locus,
                    "CHR": chrom,
                    "POS": pos,
                    "Group": "control",
                    "Match_Found": False,
                    "Num_Matches": 0,
                    "CMD_Status": status
                })

            result[locus] = {
                "case": df_case,
                "control": df_ctrl
            }

    return result, pd.DataFrame(log_records)


def compute_ld_between_focus_and_tommo_linked_variants(
    tommo_dict: dict,
    bed_prefix: str,
    case_prefix: str,
    plink2_path: str = "/home/b/b37974/plink2_alpha6/plink2",
    threads: int = 4,
    ld_window_kb: int = 2000,
    ld_window_r2: float = 0.0
) -> dict:
    """
    计算每个关注变体（focus variant）与其在 ToMMo 数据中已知的 LD 关联变体（linked variants）
    在我们自己的 case/control WGS 数据中是否也存在连锁不平衡（LD）关系。

    【函数功能】
    - 本函数用于验证 ToMMo 提供的 LD 信息是否在我们的实际样本数据中也成立。
    - 每个 focus 变体与其 ToMMo 关联变体进行 plink2 LD 计算；
    - 对于 LD 缺失的变体执行反向查询补全；
    - 最终返回每个 focus 变体在 case/control 群体中的 LD 结果表。

    【参数说明】
    - tommo_dict: 字典类型，每个键是一个关注变体（focus_id），值是该变体对应的 ToMMo LD 信息（DataFrame），需包含 ID_B、CHROM_B、POS_B 等字段；
    - bed_prefix: 输入的 plink 文件前缀（.bed/.bim/.fam）；
    - case_prefix: 用于区分 case/control 样本 ID 的前缀（如 "PHOM"）；
    - plink2_path: plink2 可执行文件的完整路径；
    - threads: 同时运行的最大子进程数（即多进程并发度）；
    - ld_window_kb: plink LD 计算的窗口大小，单位为 kb；
    - ld_window_r2: plink 输出 LD 的最小 r² 阈值。

    【返回值】
    返回一个嵌套字典 result_dict，结构如下：
        result_dict[focus_id]["case"] = DataFrame (plink2 LD 结果 + 补全信息)
        result_dict[focus_id]["control"] = DataFrame
        result_dict[focus_id]["note"] = 说明未能计算 LD 的情况

    每个结果表包含如下列：
        - CHROM_A, POS_A, ID_A：目标变体位置信息
        - CHROM_B, POS_B, ID_B：关联变体位置信息
        - UNPHASED_R2：plink 计算出的 r²
        - NOTE：若 LD 结果缺失，提供缺失原因说明，例如：
            - "ID_B missing from Genotype Data"
            - "ID_B missing from PLINK output"
            - 若反向补全成功，NOTE 留空
    """
    import hashlib
    result_dict = {}

    # ===================== 样本分组：case/control 拆分 =====================
    # 创建临时目录用于中间文件
    tmpdir = os.path.abspath("tmp")
    os.makedirs(tmpdir, exist_ok=True)

    # 生成 case/control plink 文件名
    case_plink = os.path.join(tmpdir, "case_data")
    ctrl_plink = os.path.join(tmpdir, "ctrl_data")
    fam_path = f"{bed_prefix}.fam"
    case_keep_path = os.path.join(tmpdir, "case.keep")
    ctrl_remove_path = os.path.join(tmpdir, "ctrl.remove")

    # ---------- 读取 .fam 文件，按 case_prefix 拆分样本 ----------
    case_samples = []
    ctrl_samples = []
    with open(fam_path) as fam_file:
        for line in fam_file:
            fields = line.strip().split()
            if not fields:
                continue
            if fields[0].startswith(case_prefix):
                case_samples.append((fields[0], fields[1]))
            else:
                ctrl_samples.append((fields[0], fields[1]))

    # ---------- 写入 keep/remove 文件 ----------
    with open(case_keep_path, "w") as f:
        for fid, iid in case_samples:
            f.write(f"{fid} {iid}\n")
    with open(ctrl_remove_path, "w") as f:
        for fid, iid in case_samples:
            f.write(f"{fid} {iid}\n")

    # ---------- plink2 拆分出 case/control plink 文件 ----------
    cmds = [
        [
            plink2_path, "--bfile", bed_prefix,
            "--keep", case_keep_path,
            "--make-bed", "--out", case_plink,
            "--threads", str(threads),
        ],
        [
            plink2_path, "--bfile", bed_prefix,
            "--remove", ctrl_remove_path,
            "--make-bed", "--out", ctrl_plink,
            "--threads", str(threads),
        ]
    ]
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as executor:
        futures = []
        for cmd in cmds:
            executor.submit(run_subprocess, cmd, check=True)
        executor.shutdown(wait=True)

    # ===================== 有效变体读取（BIM 文件） =====================
    bim_path = f"{bed_prefix}.bim"
    with open(bim_path) as f_bim:
        valid_vids = set(line.strip().split()[1] for line in f_bim)

    # ===================== 辅助函数定义 =====================

    # --- append_missing_rows: 用于补全缺失的 variation2_id 记录，并标注 NOTE ---
    # NOTE 列区分两种情况：
    #   - "ID_B missing from Genotype Data"：ID_B 不在 BIM 文件（基因型数据）中
    #   - "ID_B missing from PLINK output"：ID_B 在 BIM 文件中，但未出现在 plink2 LD 结果中
    def append_missing_rows(df, id_a, missing_ids, missing_ids_from_bim):
        """
        对于缺失的 variation2_id，补全行，并在 NOTE 列进行区分标注。
        参数:
          df: 当前 DataFrame
          id_a: focus 变体 ID
          missing_ids: 需要补全的 ID_B 列表
          missing_ids_from_bim: 不在 BIM 文件中的 ID_B 列表
        返回: 增补后的 DataFrame
        """
        chrom_a, pos_a = id_a.split(":")[0].replace("chr", ""), id_a.split(":")[1]
        rows = []
        for vid in missing_ids:
            try:
                chrom_b, pos_b = vid.split(":")[0].replace("chr", ""), vid.split(":")[1]
            except Exception:
                chrom_b, pos_b = None, None
            if vid in missing_ids_from_bim:
                note = "ID_B missing from Genotype Data"  # 不在 BIM 文件
            else:
                note = "ID_B missing from PLINK output"   # 在 BIM 文件但未被 PLINK 输出
            rows.append({
                '#CHROM_A': chrom_a,
                'POS_A': pos_a,
                'ID_A': id_a,
                'CHROM_B': chrom_b,
                'POS_B': pos_b,
                'ID_B': vid,
                'UNPHASED_R2': None,
                'NOTE': note
            })
        if df.empty:
            return pd.DataFrame(rows)
        else:
            df = df.copy()
            if 'NOTE' not in df.columns:
                df['NOTE'] = None
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", FutureWarning)
                df = pd.concat([df, pd.DataFrame(rows)], ignore_index=True)
            return df


    # ===================== 并行运行每个 focus_id 的 LD 计算 =====================
    with ProcessPoolExecutor(max_workers=threads) as executor:
        future_to_focus = {
            executor.submit(
                process_focus_id, focus_id, df, case_plink, ctrl_plink, valid_vids, tmpdir, plink2_path, ld_window_kb, ld_window_r2
            ): focus_id
            for focus_id, df in tommo_dict.items()
        }
        for future in concurrent.futures.as_completed(future_to_focus):
            focus_id, result = future.result()
            result_dict[focus_id] = result

    # 保存为 pickle 文件，不返回内存对象
    output_path = os.path.abspath("focus_vs_tommo_ld_result.pkl")
    with open(output_path, "wb") as f:
        pickle.dump(result_dict, f)

    return output_path # type: ignore


# ===================== 提升为模块级函数: process_focus_id =====================
def process_focus_id(
    focus_id,
    df,
    case_plink,
    ctrl_plink,
    valid_vids,
    tmpdir,
    plink2_path,
    ld_window_kb,
    ld_window_r2
):
    """
    针对单个 focus_id：
      - 提取 TOMMO 提供的 variation2_id
      - 过滤出 BIM 文件中存在的
      - 分别在 case/control 下运行 plink2 计算 LD
      - 补全缺失记录，标准化输出 DataFrame
    返回: (focus_id, {"case": df, "control": df})
    """
    import pandas as pd
    import subprocess
    import warnings
    import os
    # --- append_missing_rows: 用于补全缺失的 variation2_id 记录，并标注 NOTE ---
    def append_missing_rows(df, id_a, missing_ids, missing_ids_from_bim):
        chrom_a, pos_a = id_a.split(":")[0].replace("chr", ""), id_a.split(":")[1]
        rows = []
        for vid in missing_ids:
            try:
                chrom_b, pos_b = vid.split(":")[0].replace("chr", ""), vid.split(":")[1]
            except Exception:
                chrom_b, pos_b = None, None
            if vid in missing_ids_from_bim:
                note = "ID_B missing from Genotype Data"
            else:
                note = "ID_B missing from PLINK output"
            rows.append({
                '#CHROM_A': chrom_a,
                'POS_A': pos_a,
                'ID_A': id_a,
                'CHROM_B': chrom_b,
                'POS_B': pos_b,
                'ID_B': vid,
                'UNPHASED_R2': None,
                'NOTE': note
            })
        if df.empty:
            return pd.DataFrame(rows)
        else:
            df = df.copy()
            if 'NOTE' not in df.columns:
                df['NOTE'] = None
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", FutureWarning)
                df = pd.concat([df, pd.DataFrame(rows)], ignore_index=True)
            return df

    # ---------- TOMMO 为空时直接返回空表 ----------
    if df.empty:
        empty_df = pd.DataFrame(columns=['#CHROM_A', 'POS_A', 'ID_A', 'CHROM_B', 'POS_B', 'ID_B', 'UNPHASED_R2', 'NOTE'])
        return focus_id, {
            "case": empty_df.copy(),
            "control": empty_df.copy(),
            "note": "No ToMMo LD Record"
        }
    # variation2_id 列不存在时直接返回空表
    if "variation2_id" not in df.columns:
        empty_df = pd.DataFrame(columns=['#CHROM_A', 'POS_A', 'ID_A', 'CHROM_B', 'POS_B', 'ID_B', 'UNPHASED_R2', 'NOTE'])
        return focus_id, {"case": empty_df.copy(), "control": empty_df.copy()}
    # ---------- TOMMO 变体列表处理 ----------
    variation2_ids = df["variation2_id"].dropna().unique().tolist()
    # 去除自身
    variation2_ids = [vid for vid in variation2_ids if vid != focus_id]
    # 保留原始 TOMMO variation2_id 列表（用于最终完整性输出）
    original_vids = variation2_ids.copy()
    # 仅保留 BIM 文件中存在的
    filtered_vids = [vid for vid in variation2_ids if vid in valid_vids]
    # 记录哪些 variation2_id 不在 BIM 文件（基因型数据）中
    missing_vids_from_bim = [vid for vid in variation2_ids if vid not in valid_vids]
    # ---------- 若全部 variation2_id 不在 BIM，直接返回缺失标注空表 ----------
    if not filtered_vids:
        empty_df = pd.DataFrame(columns=[
            '#CHROM_A', 'POS_A', 'ID_A', 'CHROM_B', 'POS_B', 'ID_B', 'UNPHASED_R2', 'NOTE'
        ])
        return focus_id, {
            "case": append_missing_rows(empty_df.copy(), focus_id, variation2_ids, missing_vids_from_bim),
            "control": append_missing_rows(empty_df.copy(), focus_id, variation2_ids, missing_vids_from_bim),
            "note": "No PLINK run: all variation2_id missing in BIM"
        }
    # ---------- 写入 variation2_id 列表文件（可选，当前仅保留为备查） ----------
    with open(os.path.join(tmpdir, "focus_ld_snps.txt"), "w") as f:
        for vid in filtered_vids:
            f.write(f"{vid}\n")
    # ---------- 运行 plink2 计算 LD ----------
    ld_snps_list = list(set(filtered_vids))
    ld_snps_str = ",".join([focus_id] + ld_snps_list)
    focus_tag = focus_id.replace(":", "_")
    out_case = os.path.join(tmpdir, f"{focus_tag}_case_ld")
    out_ctrl = os.path.join(tmpdir, f"{focus_tag}_ctrl_ld")
    # 使用 threads=1 防止子进程竞争
    cmd_case = [
        plink2_path,
        "--bfile", case_plink,
        "--r2-unphased",
        "--ld-snps", ld_snps_str,
        "--ld-window-kb", str(ld_window_kb),
        "--ld-window-r2", str(ld_window_r2),
        "--out", out_case,
        "--threads", "1",
    ]
    cmd_ctrl = [
        plink2_path,
        "--bfile", ctrl_plink,
        "--r2-unphased",
        "--ld-snps", ld_snps_str,
        "--ld-window-kb", str(ld_window_kb),
        "--ld-window-r2", str(ld_window_r2),
        "--out", out_ctrl,
        "--threads", "1",
    ]
    # ---------- 运行 plink2，并读取结果（case） ----------
    try:
        subprocess.run(cmd_case, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        df_case_all = pd.read_csv(f"{out_case}.vcor", sep='\s+')
        original_df_case = df_case_all.copy()
        df_case_all = df_case_all[df_case_all['ID_A'] != df_case_all['ID_B']]
        # 仅保留 ID_B 在 filtered_vids 且 ID_A == focus_id 的行
        df_case_all = df_case_all[(df_case_all['ID_B'].isin(filtered_vids)) & (df_case_all['ID_A'] == focus_id)]
    except Exception:
        df_case_all = pd.DataFrame()
        original_df_case = pd.DataFrame()
    # ---------- 运行 plink2，并读取结果（control） ----------
    try:
        subprocess.run(cmd_ctrl, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        df_ctrl_all = pd.read_csv(f"{out_ctrl}.vcor", sep='\s+')
        original_df_ctrl = df_ctrl_all.copy()
        df_ctrl_all = df_ctrl_all[df_ctrl_all['ID_A'] != df_ctrl_all['ID_B']]
        df_ctrl_all = df_ctrl_all[(df_ctrl_all['ID_B'].isin(filtered_vids)) & (df_ctrl_all['ID_A'] == focus_id)]
    except Exception:
        df_ctrl_all = pd.DataFrame()
        original_df_ctrl = pd.DataFrame()

    # --- to_final_df: 规范化最终 DataFrame，补全缺失，并标准化列顺序 ---
    def to_final_df(df, id_a, original_vids, original_df):
        if df.empty:
            df2 = pd.DataFrame(columns=['#CHROM_A', 'POS_A', 'ID_A', 'CHROM_B', 'POS_B', 'ID_B', 'UNPHASED_R2', 'NOTE'])
        else:
            df2 = df.copy()
            df2['#CHROM_A'] = df2['ID_A'].apply(lambda x: x.split(":")[0].replace("chr", ""))
            df2['POS_A'] = df2['ID_A'].apply(lambda x: x.split(":")[1])
            df2['CHROM_B'] = df2['ID_B'].apply(lambda x: x.split(":")[0].replace("chr", ""))
            df2['POS_B'] = df2['ID_B'].apply(lambda x: x.split(":")[1])
            df2['UNPHASED_R2'] = df2['R2'] if 'R2' in df2.columns else (df2['UNPHASED_R2'] if 'UNPHASED_R2' in df2.columns else None)
            df2 = df2[['#CHROM_A', 'POS_A', 'ID_A', 'CHROM_B', 'POS_B', 'ID_B', 'UNPHASED_R2']]
            df2['NOTE'] = None
        existing_ids = set(df2['ID_B'])
        missing_ids = [vid for vid in original_vids if vid not in existing_ids]
        df2 = append_missing_rows(df2, id_a, missing_ids, missing_ids_from_bim=missing_vids_from_bim)
        df2 = df2[['#CHROM_A', 'POS_A', 'ID_A', 'CHROM_B', 'POS_B', 'ID_B', 'UNPHASED_R2', 'NOTE']]
        if not original_df.empty:
            needed_keys = set(zip(
                df2.loc[df2["NOTE"] == "ID_B missing from PLINK output", "ID_B"],
                df2.loc[df2["NOTE"] == "ID_B missing from PLINK output", "ID_A"]
            ))
            reverse_lookup_dict = {
                (row['ID_A'], row['ID_B']): row.get('R2') if 'R2' in row else row.get('UNPHASED_R2')
                for _, row in original_df.iterrows()
                if (row['ID_A'], row['ID_B']) in needed_keys
            }
            def try_reverse_lookup(row):
                if row['NOTE'] == "ID_B missing from PLINK output":
                    key = (row['ID_B'], row['ID_A'])
                    if key in reverse_lookup_dict:
                        return pd.Series({
                            "UNPHASED_R2": reverse_lookup_dict[key],
                            "NOTE": None
                        })
                return pd.Series({
                    "UNPHASED_R2": row["UNPHASED_R2"],
                    "NOTE": row["NOTE"]
                })
            df2[["UNPHASED_R2", "NOTE"]] = df2.apply(try_reverse_lookup, axis=1)
        return df2

    case_df = to_final_df(df_case_all, focus_id, original_vids, original_df_case)
    ctrl_df = to_final_df(df_ctrl_all, focus_id, original_vids, original_df_ctrl)
    return focus_id, {
        "case": case_df,
        "control": ctrl_df
    }



def merge_tommo_and_focus_ld_dict(
    tommo_ld_dict_path: str,
    focus_ld_dict_path: str,
    select: str = "control"
) -> dict[str, dict[str, pd.DataFrame | str]]:
    """
    合并 ToMMo 的 LD 结果与我们自己计算的 case/control LD 结果。

    参数:
    - tommo_ld_dict_path (str): 指向 tommo_ld_dict.pkl 文件的路径，格式为 {focus_id: DataFrame}
    - focus_ld_dict_path (str): 指向 focus_ld_dict.pkl 文件的路径，格式为 {focus_id: {'control': df, 'case': df, 'note': str (可选)}}
    - select (str): 选择合并 'control' 或 'case' 的 LD 表（默认 'control'）

    返回:
    - merged_dict (dict): 字典格式为 {focus_id: {'control': 合并后的 DataFrame, 'note': str (可选)}}
    """

    with open(tommo_ld_dict_path, "rb") as f:
        tommo_ld_dict = pickle.load(f)
    with open(focus_ld_dict_path, "rb") as f:
        focus_ld_dict = pickle.load(f)

    merged_dict = {}

    for focus_id in tommo_ld_dict:
        tommo_df = tommo_ld_dict.get(focus_id, pd.DataFrame())
        focus_entry = focus_ld_dict.get(focus_id, {})
        focus_df = focus_entry.get(select, pd.DataFrame())
        note = focus_entry.get("note") if "note" in focus_entry else None

        # 如果两个都是空表，返回空 DataFrame
        if tommo_df.empty and focus_df.empty:
            merged = pd.DataFrame(columns=focus_df.columns.tolist())
        else:
            merged = focus_df.copy()
            # 若 UNPHASED_R2 存在，则插入 TOMMO_R2 在其之前；否则插入到最后
            insert_loc = merged.columns.get_loc("UNPHASED_R2") if "UNPHASED_R2" in merged.columns else len(merged.columns)
            merged.insert(insert_loc, "TOMMO_R2", None)

            # 进行值填充：focus_df 的 ID_B == tommo_df 的 variation2_id
            if not tommo_df.empty and "variation2_id" in tommo_df.columns and "r2" in tommo_df.columns:
                tommo_r2_map = tommo_df.set_index("variation2_id")["r2"].to_dict()
                merged["TOMMO_R2"] = merged["ID_B"].map(tommo_r2_map)
            else:
                merged["TOMMO_R2"] = None

        merged_dict[focus_id] = {select: merged}
        if note is not None:
            merged_dict[focus_id]["note"] = note

    # 保存为 pickle 文件并返回路径
    output_path = os.path.abspath("merged_tommo_focus_ld.pkl")
    with open(output_path, "wb") as f:
        pickle.dump(merged_dict, f)
    return output_path # type: ignore



def plot_ld_comparison_from_merged_dict(
    merge_pkl_path: str,
    select: str = "control"
) -> str:
    """
    从合并后的 TOMMO vs 我们的 LD 结果中，绘制每个 key 的对比图，并生成 PDF。

    参数:
    - merge_pkl_path (str): 指向合并字典的 .pkl 文件路径
    - select (str): 选择绘制 'control' 或 'case'（默认 'control'）

    返回:
    - pdf_path (str): 生成的 PDF 文件路径
    """
    with open(merge_pkl_path, "rb") as f:
        merged_dict = pickle.load(f)

    pdf_path = os.path.abspath("tommo_vs_focus_ld_scatter.pdf")
    
    plt.style.use('default')  # 使用默认样式
    
    with PdfPages(pdf_path) as pdf:
        for focus_id, entry in merged_dict.items():
            df = entry.get(select, pd.DataFrame())
            note = entry.get("note")

            # 转换 TOMMO_R2 和 UNPHASED_R2 为数值型
            if "TOMMO_R2" in df.columns:
                df["TOMMO_R2"] = pd.to_numeric(df["TOMMO_R2"], errors="coerce")
            if "UNPHASED_R2" in df.columns:
                df["UNPHASED_R2"] = pd.to_numeric(df["UNPHASED_R2"], errors="coerce")

            fig, ax = plt.subplots(figsize=(6, 6))
            ax.set_title(f"{focus_id}\nToMMo LD $R^2$ vs User Genotype LD $R^2$", fontsize=12, loc="center", pad=20)

            if note:
                # 新的 note 说明逻辑
                if note == "No ToMMo LD Record":
                    display_note = "No LD information found in ToMMo for this variant."
                elif note == "No PLINK run: all variation2_id missing in BIM":
                    display_note = "ToMMo found LD variants, but none were found in the user's genotype data."
                else:
                    display_note = note
                ax.text(0.5, 0.5, display_note, fontsize=12, ha="center", va="center", wrap=True)
                ax.set_axis_off()
            else:
                df_valid = df.dropna(subset=["TOMMO_R2", "UNPHASED_R2"])
                x = df_valid["TOMMO_R2"]
                y = df_valid["UNPHASED_R2"]

                ax.scatter(x, y, alpha=0.7)
                ax.plot([0, 1], [0, 1], ls="--", color="red")
                ax.set_xlim(0, 1)
                ax.set_ylim(0, 1)
                ax.set_xlabel("ToMMo LD $R^2$")
                ax.set_ylabel(f"User Genotype LD $R^2$ ({select.capitalize()})")
                ax.set_aspect('equal', adjustable='box')

                ax.text(0.05, 0.95, f"ToMMo LD pairs available: {df['TOMMO_R2'].notna().sum()}", fontsize=10, transform=ax.transAxes, verticalalignment='top')
                ax.text(0.05, 0.90, f"LD pairs computed from user genotypes: {df['UNPHASED_R2'].notna().sum()}", fontsize=10, transform=ax.transAxes, verticalalignment='top')

            pdf.savefig(fig)
            plt.close(fig)

    return pdf_path