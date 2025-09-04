# CTEPH AGP3K 脚本说明文档

本文档介绍 `/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/script` 目录下脚本的功能与用途，适用于慢性血栓栓塞性肺动脉高压（CTEPH）基因组数据的质量控制、人群结构分析及关联分析预处理。

## 目录结构与脚本功能

### 1. 数据预处理与质量控制

- **`sample.qc.hetero.miss.rev1.py`**  
    样本杂合性与缺失度质控，支持自定义阈值，筛选高质量样本。

- **`variant.qc.miss.py`**  
    变异位点缺失率与MAF质控，生成报告与可视化。

- **`variant.parameter.py`**  
    变异参数（MQ、VQSLOD）分布分析，辅助阈值设定。

- **`hwe.check.rev1.py`**  
    Hardy-Weinberg平衡检验，支持病例/对照分组，输出可视化报告。

### 2. 测序深度与基因型质量参数优化

- **`dp.parameter.py`**  
    测序深度（DP）分析，区分不同测序深度样本，辅助参数选择。

- **`dp.parameter.vis.py`**  
    DP参数可视化，生成优化图表。

- **`gq.parameter.py`**  
    基因型质量（GQ）分析，评估不同阈值对数据保留率影响。

- **`gq.parameter.vis.rev3.py`**  
    GQ参数可视化，输出优化报告。

### 3. 人群结构分析

- **`bbj_prepare.py`**  
    BBJ参考数据预处理，导出VCF并进行MAF过滤。

- **`bbj_pca.rev2.py`**  
    BBJ主成分分析，样本ID更新、LD修剪，生成特征向量。

- **`bbj_projection.py`**  
    研究样本投影至BBJ主成分空间，评估人群结构。

- **`bbj_projection_keep.py`**  
    投影结果筛选与可视化，定义样本保留标准。

- **`pca.run.py`**  
    研究数据主成分分析，区分病例/对照，输出可视化图表。

### 4. 亲缘关系检查

- **`kinship.pi_hat.rev1.py`**  
    计算PI_HAT评估亲缘关系，移除高度相关样本，优先保留病例。

### 5. 数据整合与最终准备

- **`clean_gt.py`**  
    整合质控结果，合并筛选列表，生成最终基因型数据。

- **`update.fam.pheno.py`**  
    更新PLINK FAM文件与表型信息，生成协变量文件。

## 推荐使用流程

1. **初始质控**  
     样本与变异质控  
     ```bash
     python sample.qc.hetero.miss.rev1.py --bed_prefix input --threads 8
     python variant.qc.miss.py --bed_prefix input --threads 8 --out output
     ```

2. **参数优化**  
     DP与GQ参数分析  
     ```bash
     python dp.parameter.py --chrom 1 --vcf_path input.vcf.gz --metadata metadata.csv
     python gq.parameter.py --chrom 1 --vcf_path input.vcf.gz --metadata metadata.csv
     ```

3. **HWE检验**  
     ```bash
     python hwe.check.rev1.py --bed_prefix input --mode case_control --out hwe_filtered
     ```

4. **人群结构分析**  
     BBJ数据准备、PCA与投影  
     ```bash
     python bbj_prepare.py --bbj_bed_prefix bbj_input --threads 8
     python bbj_pca.rev2.py --bbj_bed_prefix bbj_prepared --output_prefix bbj_pca
     python bbj_projection.py --my_bed_prefix my_data --bbj_bed_prefix bbj_pca
     python bbj_projection_keep.py --projection_sscore_path projection.sscore
     ```

5. **亲缘关系检查**  
     ```bash
     python kinship.pi_hat.rev1.py --bed_prefix qc_data --sample_keep projection_keep.txt --prune_in prune.in
     ```

6. **最终数据准备**  
     数据清理与表型更新  
     ```bash
     python clean_gt.py --bed_prefix qc_data --projection_keep keep.txt --kinship_remove remove.txt
     python update.fam.pheno.py --info_file clinical.xlsx --bed_prefix final_data --case_prefix PHOM --out final
     ```

## 输出文件类型

- 质控报告与可视化图表
- 样本与变异筛选列表
- PCA特征值、特征向量及图表
- 完整质控后的PLINK基因型数据

## 依赖环境

- **Python包**：pandas, numpy, matplotlib, seaborn, pysam
- **外部工具**：PLINK2, bcftools
- **系统要求**：Linux，建议内存≥16GB

## 注意事项

- 所有脚本支持命令行参数，使用 `--help` 查看详细用法
- 建议按流程顺序执行，后续步骤依赖前面输出
- 大数据集建议多线程处理（`--threads` 参数）
- 定期备份中间结果，避免重复计算

---

*最后更新：2025年2月*  
*作者：ZHAO TIE*  
*说明：本脚本集基于旧版本，不再维护。部分函数已集成至新工具集，建议保留。*
