# CTEPH-AGP3K WGS Data Processing Pipeline

## 概述 (Overview)

本目录包含了 CTEPH-AGP3K 项目的全基因组测序 (WGS) 数据处理管道。该管道使用 Nextflow 进行流程管理，包含从原始 VCF 文件处理到质量控制、PCA 分析、BBJ 投影等完整的基因组数据预处理步骤。

## 主要文件说明 (Main Files)

### 核心流程文件
- **`select.nf`**: 主要的 Nextflow 管道脚本，定义了完整的 WGS 数据处理流程
- **`nextflow.config`**: Nextflow 配置文件，包含各个进程的计算资源分配设置
- **`scripts/`**: Python 脚本目录，包含各种质量控制和数据处理脚本

### 配置文件
- **`nextflow.config`**: 定义了每个进程的 SLURM 集群资源配置

## 流程步骤详解 (Pipeline Steps)

### 阶段 1: 基础数据过滤 (Basic Filtering)

#### 1. selectPASS (`01.pass/`)
- **功能**: 从原始 VQSR VCF 文件中提取 PASS 质量的变异位点
- **输入**: VQSR 处理后的 VCF 文件 (按染色体分割)
- **输出**: 仅包含 PASS 标记变异的 VCF 文件
- **工具**: bcftools view

#### 2. filter_ac1 (`02.mac1/`)
- **功能**: 过滤等位基因计数 ≥1 的变异 (去除单态位点)
- **输入**: PASS 过滤后的 VCF 文件
- **输出**: MAC≥1 的 VCF 文件
- **工具**: bcftools view

#### 3. set_ids (`03.set_ids/`)
- **功能**: 为变异位点设置标准化的 ID 格式
- **输入**: MAC≥1 的 VCF 文件
- **输出**: 具有标准化 ID 的 VCF 文件
- **工具**: bcftools annotate

### 阶段 2: 注释和标准化 (Annotation & Normalization)

#### 4. add_gt_AF (`04.add_gt_AF/`)
- **功能**: 添加基因型和等位基因频率信息
- **输入**: 标准化 ID 的 VCF 文件
- **输出**: 包含 AF 信息的 VCF 文件
- **工具**: bcftools +fill-tags

#### 5. gt_norm (`05.gt_norm/`)
- **功能**: 基因型标准化处理
- **输入**: 带 AF 信息的 VCF 文件
- **输出**: 标准化的 VCF 文件
- **工具**: bcftools norm

#### 6. variant_filter (`06.variant_filter/`)
- **功能**: 应用变异位点过滤条件
- **输入**: 标准化的 VCF 文件
- **输出**: 过滤后的 VCF 文件

### 阶段 3: 质量控制准备 (QC Preparation)

#### 7. gt_qc (`07.gt_qc/`)
- **功能**: 基因型质量控制预处理
- **输入**: 过滤后的 VCF 文件
- **输出**: QC 预处理的 VCF 文件

#### 8. gt_norm_2 (`08.gt_norm_2/`)
- **功能**: 二次基因型标准化
- **输入**: QC 预处理的 VCF 文件
- **输出**: 二次标准化的 VCF 文件

#### 9. het_qt_qc (`09.het_qt_qc/`)
- **功能**: 杂合子质量控制
- **输入**: 二次标准化的 VCF 文件
- **输出**: 杂合子 QC 的 VCF 文件

### 阶段 4: 格式转换 (Format Conversion)

#### 10. vcf2bed (`10.vcf2bed/`)
- **功能**: 将 VCF 文件转换为 PLINK 格式 (BED/BIM/FAM)
- **输入**: 杂合子 QC 的 VCF 文件
- **输出**: PLINK 二进制文件集
- **工具**: plink2

### 阶段 5: 样本质量控制 (Sample QC)

#### 11. RunSampleQC (`11.run_sample_qc/`)
- **功能**: 执行全面的样本质量控制分析
- **脚本**: `scripts/sample_qc_main.py`
- **主要指标**:
  - 样本缺失率 (Sample Missing Rate)
  - 杂合度检测 (Heterozygosity)
  - 深度分布分析 (Depth Distribution)
  - 亲缘关系检测 (PI_HAT)
- **输出**: 样本 QC 统计结果和可视化图表

#### 12. rmMAF0orVMISS1 (`12.rm_maf0_vmiss1/`)
- **功能**: 移除 MAF=0 或 VMISS=1 的变异位点
- **输入**: PLINK 格式文件
- **输出**: 过滤后的 PLINK 文件
- **工具**: plink2

### 阶段 6: 变异质量控制 (Variant QC)

#### 13. RunVariantQC (`13.run_variant_qc/`)
- **功能**: 执行变异位点质量控制
- **脚本**: `scripts/variant_qc_main.py`
- **主要指标**:
  - 变异缺失率 (Variant Missing Rate)
  - Hardy-Weinberg 平衡检验 (HWE)
  - MAF 分层分析
- **输出**: 变异 QC 统计结果和通过的变异列表

### 阶段 7: 主成分分析 (PCA Analysis)

#### 14. RunPCA (`14.run_pca/`)
- **功能**: 执行主成分分析
- **脚本**: `scripts/pca_qc_main.py`
- **输入**: 质量控制后的 PLINK 文件
- **输出**: PCA 结果文件和可视化图表

### 阶段 8: BBJ 数据库投影分析 (BBJ Projection)

#### 15. BBJSampleKeep (`15.bbj_sample_keep/`)
- **功能**: 根据 BBJ 投影结果筛选保留的样本
- **脚本**: `scripts/bbj_sample_keep_main.py`
- **输入**: BBJ 投影结果
- **输出**: 筛选后的样本列表

#### 16. rmMAF0orVMISS1_repeat (`16.rm_maf0_vmiss1_repeat/`)
- **功能**: 重复执行 MAF/VMISS 过滤 (基于更新的样本集)
- **输入**: 更新样本集的 PLINK 文件
- **输出**: 再次过滤的 PLINK 文件

### 阶段 9: ToMMo 面板比较分析 (ToMMo Panel Analysis)

#### 17. ToMMoPanelCompare (`17.tommo_panel_compare/`)
- **功能**: 与 ToMMo 60KJPN 面板进行比较分析
- **脚本**: `scripts/panel_compare_main.py`
- **输入**: 质量控制后的数据
- **输出**: 面板比较结果

#### 18. ToMMoPanelThr (`18.tommo_panel_thr/`)
- **功能**: 应用 ToMMo 面板比较的阈值筛选
- **脚本**: `scripts/panel_thr_main.py`
- **输入**: 面板比较结果
- **输出**: 阈值筛选结果

#### 19. ToMMoPanelFilter (`19.tommo_panel_filter/`)
- **功能**: 基于 ToMMo 面板比较结果进行最终过滤
- **脚本**: `scripts/panel_filter_main.py`
- **输入**: 阈值筛选结果
- **输出**: 最终过滤的数据集

### 阶段 10: 最终数据准备 (Final Data Preparation)

#### 20. CovPhenoPrep (`20.cov_pheno_prepare/`)
- **功能**: 准备协变量和表型数据用于下游 GWAS 分析
- **脚本**: `scripts/cov_pheno_prepare_rev1.py`
- **输入**: 最终筛选的样本和 BBJ 投影结果
- **输出**: 
  - `cteph_agp3k.bbj.projection.cov_df.csv` - 协变量数据框
  - `cteph_agp3k.bbj.projection.pheno_df.csv` - 表型数据框
  - `cteph_agp3k.bbj.projection.cov_df.no_age.csv` - 无年龄信息的协变量
  - `cteph_agp3k.bbj.projection.missing_age_samples.csv` - 缺失年龄信息的样本

#### 21. MissBiasFilter (`21.miss_bias_filter/`)
- **功能**: 缺失偏倚过滤
- **脚本**: `scripts/miss_bias_main.py`
- **输入**: 最终数据集
- **输出**: 缺失偏倚过滤后的数据

## 核心 Python 脚本说明 (Core Python Scripts)

### 样本质量控制模块
- **`sample_qc_main.py`**: 样本 QC 主脚本
- **`sample_qc_pipeline.py`**: 样本 QC 流程管理
- **`sample_qc_calculator.py`**: 样本 QC 指标计算
- **`sample_qc_flags.py`**: 样本 QC 标记生成

### 变异质量控制模块
- **`variant_qc_main.py`**: 变异 QC 主脚本
- **`variant_qc_calculator.py`**: 变异 QC 指标计算
- **`variant_qc_flags.py`**: 变异 QC 标记生成
- **`variant_qc_rm_maf0_vmiss1.py`**: MAF/VMISS 过滤脚本

### PCA 分析模块
- **`pca_qc_main.py`**: PCA 主脚本
- **`pca_qc_tools.py`**: PCA 分析工具

### BBJ 投影分析模块
- **`bbj_prepare_main.py`**: BBJ 数据准备
- **`bbj_pca_main.py`**: BBJ PCA 分析
- **`bbj_projection_main.py`**: BBJ 投影分析
- **`bbj_projection_tools.py`**: BBJ 投影工具
- **`bbj_sample_keep_main.py`**: BBJ 样本筛选

### ToMMo 面板分析模块
- **`panel_compare_main.py`**: 面板比较主脚本
- **`panel_compare_tools.py`**: 面板比较工具
- **`panel_thr_main.py`**: 面板阈值处理
- **`panel_filter_main.py`**: 面板过滤主脚本

### 其他工具模块
- **`cov_pheno_prepare_rev1.py`**: 协变量和表型数据准备
- **`miss_bias_main.py`**: 缺失偏倚分析主脚本
- **`miss_bias_tools.py`**: 缺失偏倚分析工具
- **`random_plink_subset.py`**: 随机 PLINK 子集选择

## 运行说明 (Execution Instructions)

### 环境要求
- Nextflow
- SLURM 集群环境
- bcftools
- plink2
- Python 3.x with pandas, matplotlib, numpy
- Singularity (for containerized tools)

### 运行命令
```bash
# 在 wgs 目录下运行完整管道
nextflow run select.nf -c nextflow.config

# 从特定步骤开始运行 (示例)
nextflow run select.nf -c nextflow.config --entry <process_name>
```

### 重要参数设置
管道的关键参数在 `select.nf` 文件顶部定义：
- `params.filePath`: 输入 VQSR VCF 文件路径
- `params.samplelist`: 样本列表文件
- `params.tommodir`: ToMMo 数据库路径
- `params.outdir`: 输出目录

## 输出文件结构 (Output Structure)

每个处理步骤都有对应的输出目录 (`01.pass/` 到 `21.miss_bias_filter/`)，包含：
- 处理后的数据文件 (VCF/PLINK 格式)
- 质量控制统计结果
- 可视化图表 (PNG/PDF)
- 日志文件

## 质量控制标准 (QC Standards)

### 样本 QC 阈值
- 样本缺失率: < 10%
- PI_HAT (亲缘关系): < 0.2
- 杂合度: 5SD 范围内
- 测序深度: 稳健 Z 值 > -3.0

### 变异 QC 阈值
- 变异缺失率: < 5%
- HWE p-value: 
  - 常见变异 (MAF>0.05): 控制组 > 1e-6, 病例组 > 1e-10
  - 低频变异 (0.01≤MAF≤0.05): 控制组 > 1e-6
  - 稀有变异 (MAF<0.01): 无 HWE 过滤

## 注意事项 (Important Notes)

1. **计算资源**: 该管道需要大量计算资源，建议在高性能集群上运行
2. **存储空间**: 中间文件会占用大量存储空间，确保有足够的磁盘空间
3. **依赖关系**: 确保所有必需的软件工具已正确安装和配置
4. **日志监控**: 定期检查 `.nextflow.log` 文件以监控运行状态

## 联系信息 (Contact)

- 作者: ZHAO TIE
- 项目: CTEPH-AGP3K
- 更新日期: 2025年9月

---

此管道为 CTEPH-AGP3K 项目的核心数据处理流程，经过充分测试和优化，适用于大规模 WGS 数据的质量控制和预处理。