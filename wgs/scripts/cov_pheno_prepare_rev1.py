#!/usr/bin/env python3
"""
从 FAM、JHR 信息（Excel）与 BBJ 投影 sscore 生成 GWAS 所需的 phenotype 与 covariate 文件。

功能概述
--------
- 仅保留 (#FID, IID) 存在于 FAM 的 sscore 观测。
- 与 info.xlsx 通过 IID==ID 左连接合并，重命名：'Age at time of DNA collection'->'AGE'，'Sex'->'SEX'。
- 生成 pheno_df（#FID, IID, PHENO1），PHENO1 由 IID 前缀判定：病例=2，对照=1。
- 生成 cov_df（#FID, IID, AGE, AGE_Z, SEX, PC1_AVG~PC10_AVG）。SEX：M→1，F→2，其它→0；AGE_Z 为 AGE 的 Z 分数（仅对非缺失计算，缺失保持 NA）。
- 生成两种 cov_df：① 全量样本但不含 AGE/AGE_Z；② 仅包含 AGE 与 AGE_Z 非缺失样本，且保留 AGE/AGE_Z。

作者: ZHAO TIE
"""
from typing import Tuple
import argparse
import logging

def prepare_pheno_cov(
    fam_path: str,
    info_path: str,
    bbj_sscore_path: str,
    case_prefix: str,
    output_dir: str | None = None,
    output_prefix: str | None = None,
) -> Tuple[str, str]:
    """
    功能: 基于 FAM、JHR 信息表 (Excel) 与 BBJ 投影的 sscore 文件，构建 GWAS 所需的
    phenotype (pheno_df) 与 covariate (cov_df) 文件，并写出为**tab 分隔的 CSV** 文件，
    返回写出文件的路径 (pheno_path, cov_path)。

    输入参数:
    - fam_path: FAM 文件路径。无表头，前两列分别为 FID 与 IID，示例: `0000063134\t0000063134\t0\t0\t1\t1`。
    - info_path: Excel 文件路径，仅提取列: 'ID', 'Age at time of DNA collection', 'Sex'。
    - bbj_sscore_path: BBJ 投影 sscore 文件路径，包含列名: '#FID', 'IID', 'PHENO1',
      'ALLELE_CT', 'NAMED_ALLELE_DOSAGE_SUM', 'PC1_AVG'~'PC10_AVG'。
    - case_prefix: 病例样本在 IID 中的前缀标识，例如 'PHOM'。若 IID 以此前缀开头，则 PHENO1 置为 2，否则置为 1。
    - output_dir: 输出目录。默认为当前工作目录（os.getcwd()）。
    - output_prefix: 输出文件名前缀。默认取自 `bbj_sscore_path` 的文件名前缀(去掉扩展名)。

    处理流程:
    1) 仅保留 `bbj_sscore` 中 (#FID, IID) 组合出现在 FAM 前两列组合中的观测。
    2) 以筛选后的 `bbj_sscore` 为基准，按 `IID` = `info.ID` 左连接合并信息表；
       并将 'Age at time of DNA collection' 重命名为 'AGE'，'Sex' 重命名为 'SEX'。
    3) 生成 `pheno_df`: 提取并按顺序排列 `#FID`, `IID`, `PHENO1` 列；
       其中 `PHENO1` 根据 IID 是否以 `case_prefix` 开头重新编码: case=2, control=1。
    4) 生成 `cov_df`: 提取并按顺序排列
       `#FID`, `IID`, `AGE`, `AGE_Z`, `SEX`, `PC1_AVG`...`PC10_AVG`，并完成:
         - `SEX` 映射: M -> 1, F -> 2, 其他 -> 0 (大小写不敏感，仅取首字母)。
         - 在 `AGE` 列之后、`SEX` 之前插入一列 `AGE_Z`，对所有非 NA 的 AGE 做 Z-score 变换；
           AGE 缺失的观测其 AGE_Z 依然为 NA。同时将缺失 AGE 的样本 `#FID`、`IID` 写出为
           `missing_age_samples.tsv` (tab 分隔)。
    5) 写出以下文件（均为 tab 分隔 CSV）：
       - pheno_df：<prefix>.pheno_df.csv
       - cov_df（仅非缺失 AGE/AGE_Z）：<prefix>.cov_df.csv
       - cov_df（全量样本、无 AGE/AGE_Z 列）：<prefix>.cov_df.no_age.csv

    返回:
    - (pheno_path, cov_path): 两个文件的绝对路径。

    注意:
    - 为了避免字符串前导零丢失，FID/IID/ID 一律按字符串读取与处理。
    - 本函数对 Excel 的 'Sex' 字段做宽松映射: 取首字母并大写后映射到 {M,F}；不匹配则记为 0。
    """
    import os
    import pandas as pd
    import numpy as np

    # ---- 解析输出目录/文件前缀 ----
    if output_dir is None:
        output_dir = os.getcwd()
    os.makedirs(output_dir, exist_ok=True)

    if output_prefix is None:
        stem = os.path.basename(bbj_sscore_path)
        # 去除多重扩展名，例如 .sscore 或 .txt.gz 等
        for ext in [".gz", ".bgz", ".bz2", ".xz", ".zip", ".sscore", ".txt", ".tsv"]:
            if stem.endswith(ext):
                stem = stem[: -len(ext)]
        output_prefix = stem

    # ---- 读取 FAM (无表头) ----
    fam = pd.read_csv(
        fam_path,
        sep=r"\s+",
        header=None,
        usecols=[0, 1],
        names=["FID", "IID"],
        dtype={"FID": str, "IID": str},
        engine="python",
    )

    # ---- 读取 sscore ----
    sscore = pd.read_csv(
        bbj_sscore_path,
        sep="\t",
        dtype={"#FID": str, "IID": str},
        engine="python",
    )

    # ---- 仅保留 (#FID, IID) 在 fam 中存在的观测 ----
    fam_pairs = fam.assign(_key=fam["FID"] + "\t" + fam["IID"]) ["_key"]
    fam_set = set(fam_pairs.tolist())
    sscore["_key"] = sscore["#FID"].astype(str) + "\t" + sscore["IID"].astype(str)
    sscore_f = sscore[sscore["_key"].isin(fam_set)].drop(columns=["_key"]).reset_index(drop=True)

    # ---- 读取并精简 info (Excel) ----
    info = pd.read_excel(
        info_path,
        dtype={"ID": str},
        engine=None,
    )
    keep_cols = ["ID", "Age at time of DNA collection", "Sex"]
    missing_cols = [c for c in keep_cols if c not in info.columns]
    if missing_cols:
        raise ValueError(f"info 缺失所需列: {missing_cols}")
    info_s = info[keep_cols].rename(
        columns={"Age at time of DNA collection": "AGE", "Sex": "SEX"}
    )

    # ---- 合并: sscore 基准，IID 对应 info.ID ----
    merged = sscore_f.merge(info_s, how="left", left_on="IID", right_on="ID")
    # 合并后不再需要 ID 列
    if "ID" in merged.columns:
        merged = merged.drop(columns=["ID"])  # IID 已在左侧

    # ---- 生成 pheno_df (#FID, IID, PHENO1) ----
    # 根据 IID 前缀重编码 PHENO1: case=2, control=1
    iid_str = merged["IID"].astype(str)
    is_case = iid_str.str.startswith(str(case_prefix))
    pheno = merged[["#FID", "IID"]].copy()
    pheno["PHENO1"] = np.where(is_case, 2, 1).astype(int)

    # ---- 生成 cov_df (#FID, IID, AGE, AGE_Z, SEX, PC1_AVG..PC10_AVG) ----
    # AGE: 数值化
    merged["AGE"] = pd.to_numeric(merged["AGE"], errors="coerce")

    # 计算 AGE_Z (仅对非 NA 参与均值/标准差计算；默认 ddof=1)
    age_valid = merged["AGE"].dropna()
    if len(age_valid) > 1:
        age_mean = age_valid.mean()
        age_std = age_valid.std()  # ddof=1
        merged["AGE_Z"] = (merged["AGE"] - age_mean) / age_std
    else:
        merged["AGE_Z"] = np.nan

    # SEX: 宽松映射
    def _map_sex(x):
        if pd.isna(x):
            return 0
        s = str(x).strip().upper()
        if len(s) > 0:
            c = s[0]
            if c == "M":
                return 1
            if c == "F":
                return 2
        return 0

    merged["SEX"] = merged["SEX"].apply(_map_sex).astype(int)

    # PCs: 确保存在并排序
    pc_cols = [f"PC{i}_AVG" for i in range(1, 11)]
    for c in pc_cols:
        if c not in merged.columns:
            raise ValueError(f"sscore 缺失列: {c}")

    cov_cols = ["#FID", "IID", "AGE", "AGE_Z", "SEX"] + pc_cols
    cov = merged[cov_cols].copy()

    # 构建两种 cov_df：
    # 1) 全量样本但不含 AGE/AGE_Z 列
    cov_no_age = cov.drop(columns=["AGE", "AGE_Z"], errors="ignore")
    # 2) 仅保留 AGE 与 AGE_Z 非缺失的样本（保留 AGE/AGE_Z 列）
    cov_with_age = cov.dropna(subset=["AGE", "AGE_Z"]).copy()

    # ---- 写出缺失 AGE 的样本列表 (#FID, IID) ----
    missing_age = cov[cov["AGE"].isna()][["#FID", "IID"]].copy()
    missing_age_path = os.path.join(output_dir, f"{output_prefix}.missing_age_samples.csv")
    if not missing_age.empty:
        missing_age.to_csv(missing_age_path, sep="\t", index=False)
    else:
        # 如果无缺失，仍写出空表头，便于审计
        missing_age.to_csv(missing_age_path, sep="\t", index=False)

    # ---- 写出 pheno/cov 文件 ----
    pheno_path = os.path.join(output_dir, f"{output_prefix}.pheno_df.csv")
    cov_path = os.path.join(output_dir, f"{output_prefix}.cov_df.csv")  # 带 AGE/AGE_Z，且二者均非缺失
    cov_no_age_path = os.path.join(output_dir, f"{output_prefix}.cov_df.no_age.csv")  # 全量样本但无 AGE/AGE_Z

    pheno.to_csv(pheno_path, sep="\t", index=False)
    cov_with_age.to_csv(cov_path, sep="\t", index=False)
    cov_no_age.to_csv(cov_no_age_path, sep="\t", index=False)

    logging.info(f"cov_df(no_age) 写出: {os.path.abspath(cov_no_age_path)}")

    # 额外文件 cov_no_age 写出，但不在返回值中
    return os.path.abspath(pheno_path), os.path.abspath(cov_path)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="基于 FAM/INFO/BBJ sscore 生成 GWAS 的 phenotype 与 covariate（tab-delimited CSV）"
    )
    parser.add_argument("--fam_path", required=True, help="FAM 文件路径（无表头，前两列为 FID/IID）")
    parser.add_argument("--info_path", required=True, help="信息表 Excel 路径（需包含 ID、Age at time of DNA collection、Sex）")
    parser.add_argument("--bbj_sscore_path", required=True, help="BBJ 投影 sscore 文件路径")
    parser.add_argument("--case_prefix", required=True, help="病例样本 IID 前缀（如 PHOM）")
    parser.add_argument("--output_dir", default=None, help="输出目录（默认: 当前工作目录）")
    parser.add_argument("--output_prefix", default=None, help="输出文件名前缀（默认取自 sscore 文件名）")
    parser.add_argument("--log_level", default="INFO", choices=["DEBUG","INFO","WARNING","ERROR","CRITICAL"], help="日志级别")

    args = parser.parse_args()

    logging.basicConfig(
        level=getattr(logging, args.log_level.upper(), logging.INFO),
        format="[%(levelname)s] %(message)s"
    )

    try:
        pheno_fp, cov_fp = prepare_pheno_cov(
            fam_path=args.fam_path,
            info_path=args.info_path,
            bbj_sscore_path=args.bbj_sscore_path,
            case_prefix=args.case_prefix,
            output_dir=args.output_dir,
            output_prefix=args.output_prefix,
        )
        logging.info(f"pheno_df 写出: {pheno_fp}")
        logging.info(f"cov_df 写出: {cov_fp}")
    except Exception as e:
        logging.exception("执行失败")
        raise SystemExit(1)
