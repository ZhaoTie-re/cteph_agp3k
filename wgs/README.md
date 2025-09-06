# CTEPH-AGP3K WGS 数据处理管道

## 概述 (Overview)

本目录包含了 CTEPH-AGP3K 项目的全基因组测序 (WGS) 数据处理管道。该管道使用 Nextflow 进行流程管理，包含从原始 VCF 文件处理到质量控制、PCA 分析、BBJ 投影等完整的基因组数据预处理步骤。

## 目录结构 (Directory Structure)

```
wgs/
├── select.nf                    # 主要的 Nextflow 管道脚本
├── nextflow.config              # Nextflow 配置文件
├── README.md                    # 本文档
└── scripts/                     # Python 脚本目录
    ├── sample_qc_main.py        # 样本质量控制主脚本
    ├── sample_qc_pipeline.py    # 样本QC流程管理
    ├── sample_qc_calculator.py  # 样本QC指标计算
    ├── sample_qc_flags.py       # 样本QC标记生成
    ├── variant_qc_main.py       # 变异质量控制主脚本
    ├── variant_qc_calculator.py # 变异QC指标计算
    ├── variant_qc_flags.py      # 变异QC标记生成
    ├── variant_qc_rm_maf0_vmiss1.py # MAF/VMISS过滤脚本
    ├── pca_qc_main.py           # PCA分析主脚本
    ├── pca_qc_tools.py          # PCA分析工具
    ├── bbj_prepare_main.py      # BBJ数据准备脚本
    ├── bbj_pca_main.py          # BBJ PCA分析脚本
    ├── bbj_projection_main.py   # BBJ投影分析主脚本
    ├── bbj_projection_tools.py  # BBJ投影分析工具
    ├── bbj_sample_keep_main.py  # BBJ样本筛选脚本
    ├── panel_compare_main.py    # ToMMo面板比较主脚本
    ├── panel_compare_tools.py   # ToMMo面板比较工具
    ├── panel_thr_main.py        # ToMMo面板阈值处理脚本
    ├── panel_filter_main.py     # ToMMo面板过滤主脚本
    ├── cov_pheno_prepare_rev1.py # 协变量和表型数据准备脚本
    ├── miss_bias_main.py        # 缺失偏倚分析主脚本
    ├── miss_bias_tools.py       # 缺失偏倚分析工具
    ├── random_plink_subset.py   # 随机PLINK子集选择工具
    ├── hwe.json                 # HWE检验配置文件
    ├── panel_select_config.json # 面板选择配置文件
    └── test.ipynb               # 测试notebook
```

## 主要文件说明 (Main Files)

### 核心流程文件
- **`select.nf`**: 主要的 Nextflow 管道脚本，定义了完整的 WGS 数据处理流程
- **`nextflow.config`**: Nextflow 配置文件，包含各个进程的计算资源分配设置
- **`scripts/`**: Python 脚本目录，包含各种质量控制和数据处理脚本

### 配置文件
- **`nextflow.config`**: 定义了每个进程的 SLURM 集群资源配置
- **`scripts/hwe.json`**: Hardy-Weinberg 平衡检验的参数配置
- **`scripts/panel_select_config.json`**: ToMMo 面板选择的配置参数

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

**样本质控方法详解：**

1. **样本缺失率 (SMISS)**  
   统计每个样本的基因型缺失比例，绘制分布图。

2. **杂合度检测 (Heterozygosity F)**  
   计算每个样本的杂合性 F 值，结合均值±5SD进行异常检测。

3. **测序深度分布 (Mean DP)**  
   统计每个样本的平均测序深度，使用稳健Z分数识别低深度样本。

4. **亲缘关系检测 (PI_HAT)**  
   - 计算所有样本两两之间的 PI_HAT 值，筛选高亲缘关系对（PI_HAT > 0.2）。
   - 绘制 PI_HAT 分布直方图，区分 case-case、control-control、case-control。
   - 构建亲缘网络图，采用最小加权顶点覆盖算法（MWVC）自动识别需删除的样本节点，优先保留病例（case），并结合缺失率权重。
   - 删除样本流程：  
     1. 对所有高 PI_HAT 的样本对，构建网络；
     2. 节点权重 = case样本权重高，且缺失率低者优先保留；
     3. MWVC算法自动选择需删除的样本，生成删除名单；
     4. 在最终 QC 标记文件中，PASS_PI_HAT 为 False 的样本即为被删除对象。

5. **综合 QC 标记输出**  
   - 生成 `sample_qc_flags.csv`，包含每个样本的各项 QC 标记（PASS_SMISS, PASS_MEAN_DP, PASS_HET_F, PASS_PI_HAT）。
   - 可视化 QC 关系图表。

### 阶段 6: 变异质量控制 (Variant QC)

#### 13. RunVariantQC (`13.run_variant_qc/`)
- **功能**: 执行变异位点质量控制
- **脚本**: `scripts/variant_qc_main.py`
- **主要指标**:
  - 变异缺失率 (Variant Missing Rate)
  - Hardy-Weinberg 平衡检验 (HWE)
  - MAF 分层分析
- **输出**: 变异 QC 统计结果和通过的变异列表

**变异质控方法详解：**

1. **变异缺失率 (VMISS)**  
   统计每个位点的缺失率，绘制分布图。

2. **MAF分层分析**  
   按 MAF 分层，分别统计各层的 QC 指标。

3. **Hardy-Weinberg 平衡检验 (HWE)**  
   - 对病例组和对照组分别进行 HWE 检验，采用不同阈值（常见变异、低频变异、稀有变异）。
   - 输出 HWE 检验结果和可视化图表。

4. **最终变异筛选**  
   - 结合缺失率、HWE、MAF等多项指标，生成通过 QC 的变异列表。
   - 输出 QC 标记文件和筛选建议。

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

### 阶段 9: ToMMo 面板比较分析 (ToMMo Panel Analysis)

#### 17. ToMMoPanelCompare (`17.tommo_panel_compare/`)
- **功能**: 与 ToMMo 60KJPN 面板进行比较分析
- **脚本**: `scripts/panel_compare_main.py`
- **输入**: 质量控制后的数据
- **输出**: 面板比较结果

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

### 1. 样本质量控制模块 (Sample QC Module)

#### `sample_qc_main.py`
- **功能**: 样本质量控制的主控脚本
- **主要任务**:
  - 协调整个样本QC流程
  - 调用各种QC计算器和标记生成器
  - 生成最终的样本QC报告
- **输入**: PLINK格式文件 (BED/BIM/FAM)
- **输出**: 样本QC统计结果、可视化图表、通过QC的样本列表

#### `sample_qc_pipeline.py`
- **功能**: 样本QC流程管理器
- **主要任务**:
  - 定义QC流程的执行顺序
  - 管理不同QC步骤之间的数据传递
  - 处理异常情况和错误恢复
- **特点**: 模块化设计，便于维护和扩展

#### `sample_qc_calculator.py`
- **功能**: 样本QC指标计算引擎
- **计算指标**:
  - 样本缺失率 (Sample Missing Rate)
  - 杂合度 (Heterozygosity Rate)
  - 测序深度分布 (Depth Distribution)
  - 亲缘关系系数 (PI_HAT)
  - F系数 (Inbreeding Coefficient)
- **算法**: 使用稳健统计方法识别异常样本

#### `sample_qc_flags.py`
- **功能**: 样本QC标记生成器
- **标记类型**:
  - `HIGH_MISSING`: 高缺失率样本
  - `OUTLIER_HET`: 杂合度异常样本
  - `RELATED`: 高亲缘关系样本
  - `LOW_DEPTH`: 低测序深度样本
- **输出**: 样本标记文件和过滤建议

### 2. 变异质量控制模块 (Variant QC Module)

#### `variant_qc_main.py`
- **功能**: 变异质量控制的主控脚本
- **主要任务**:
  - 执行变异位点的质量评估
  - 应用多层次过滤标准
  - 生成变异QC报告和统计图表
- **输入**: 经过初步过滤的VCF文件
- **输出**: 变异QC统计、通过QC的变异列表

#### `variant_qc_calculator.py`
- **功能**: 变异QC指标计算引擎
- **计算指标**:
  - 变异缺失率 (Variant Missing Rate)
  - Hardy-Weinberg平衡检验 (HWE test)
  - 等位基因频率 (Allele Frequency)
  - 变异质量评分 (Variant Quality Scores)
- **特殊处理**: 按MAF分层进行不同的QC标准

#### `variant_qc_flags.py`
- **功能**: 变异QC标记生成器
- **标记类型**:
  - `HIGH_MISSING`: 高缺失率变异
  - `HWE_FAIL`: HWE检验失败变异
  - `LOW_QUAL`: 低质量变异
  - `MONO`: 单态变异
- **过滤策略**: 基于变异类型和频率的差异化过滤

#### `variant_qc_rm_maf0_vmiss1.py`
- **功能**: 专门用于移除MAF=0或VMISS=1的变异
- **算法**: 高效识别和移除无信息变异位点
- **优化**: 针对大规模数据集的内存优化算法

### 3. 主成分分析模块 (PCA Analysis Module)

#### `pca_qc_main.py`
- **功能**: 主成分分析的主控脚本
- **分析步骤**:
  - LD修剪 (Linkage Disequilibrium Pruning)
  - PCA计算 (Principal Component Analysis)
  - 人群分层检测 (Population Stratification Detection)
  - 异常样本识别 (Outlier Detection)
- **输出**: PC坐标、载荷矩阵、可视化图表

### 4. BBJ投影分析模块 (BBJ Projection Module)

#### `bbj_prepare_main.py`
- **功能**: BBJ数据库投影分析的数据准备
- **主要任务**:
  - 提取与BBJ数据库共同的变异位点
  - 标准化等位基因编码
  - 处理缺失数据
- **输出**: 标准化的投影分析输入文件

#### `bbj_pca_main.py`
- **功能**: BBJ数据库的PCA分析
- **分析内容**:
  - 在BBJ参考数据上执行PCA
  - 生成BBJ人群的主成分空间
  - 为投影分析准备参考坐标系
- **参考数据**: 使用BBJ数据库的高质量样本

#### `bbj_projection_main.py`
- **功能**: BBJ投影分析主脚本
- **投影算法**:
  - 将研究样本投影到BBJ主成分空间
  - 计算投影坐标和置信区间
  - 识别人群异常样本
- **质控标准**: 基于BBJ人群分布的异常值检测

#### `bbj_projection_tools.py`
- **功能**: BBJ投影分析工具集
- **工具函数**:
  - 投影坐标计算算法
  - 人群归属判断
  - 投影质量评估
  - 结果可视化工具

#### `bbj_sample_keep_main.py`
- **功能**: 基于BBJ投影结果的样本筛选
- **筛选标准**:
  - 投影距离阈值
  - 人群归属一致性
  - 投影质量评分
- **输出**: 通过BBJ投影QC的样本列表

### 5. ToMMo面板分析模块 (ToMMo Panel Analysis Module)

#### `panel_compare_main.py`
- **功能**: 与ToMMo 60KJPN面板的比较分析主脚本
- **比较内容**:
  - 等位基因频率比较
  - 变异位点覆盖度分析
  - 人群特异性变异识别
- **参考数据**: ToMMo 60KJPN全基因组参考面板

#### `panel_compare_tools.py`
- **功能**: ToMMo面板比较分析工具集
- **工具函数**:
  - 频率相关性计算
  - 异常频率变异检测
  - 比较结果可视化
  - 统计显著性检验

#### `panel_thr_main.py`
- **功能**: ToMMo面板比较结果的阈值处理
- **阈值设定**:
  - 频率差异阈值
  - 相关系数阈值
  - 覆盖度阈值
- **输出**: 阈值筛选结果和建议过滤列表

#### `panel_filter_main.py`
- **功能**: 基于ToMMo面板比较的最终过滤
- **过滤策略**:
  - 多维度综合评估
  - 保守过滤策略
  - 敏感性分析
- **输出**: 最终通过面板QC的变异列表

### 6. 协变量和表型准备模块 (Covariate & Phenotype Module)

#### `cov_pheno_prepare_rev1.py`
- **功能**: 为下游GWAS分析准备协变量和表型数据
- **数据整合**:
  - 样本信息匹配
  - BBJ投影结果整合
  - 临床信息处理
  - 缺失数据处理
- **输出文件**:
  - `cteph_agp3k.bbj.projection.cov_df.csv`: 完整协变量矩阵
  - `cteph_agp3k.bbj.projection.pheno_df.csv`: 表型数据矩阵
  - `cteph_agp3k.bbj.projection.cov_df.no_age.csv`: 无年龄协变量矩阵
  - `cteph_agp3k.bbj.projection.missing_age_samples.csv`: 年龄缺失样本列表

### 7. 缺失偏倚分析模块 (Missing Bias Analysis Module)

#### `miss_bias_main.py`
- **功能**: 缺失偏倚分析的主控脚本
- **分析内容**:
  - 系统性缺失模式检测
  - 病例对照缺失差异分析
  - 缺失偏倚影响评估
- **统计方法**: 使用卡方检验和Fisher精确检验

#### `miss_bias_tools.py`
- **功能**: 缺失偏倚分析工具集
- **工具函数**:
  - 缺失模式识别算法
  - 偏倚统计量计算
  - 缺失数据可视化
  - 偏倚校正建议

### 8. 辅助工具模块 (Utility Tools)

#### `random_plink_subset.py`
- **功能**: 随机选择PLINK数据子集
- **应用场景**:
  - 测试数据生成
  - 计算资源优化
  - 方法验证
- **算法**: 保持LD结构的分层随机抽样

### 9. 配置文件说明 (Configuration Files)

#### `hwe.json`
- **功能**: Hardy-Weinberg平衡检验的参数配置
- **参数设置**:
  - 不同MAF层级的HWE阈值
  - 病例组和对照组的差异化标准
  - 检验方法选择

#### `panel_select_config.json`
- **功能**: ToMMo面板选择和比较的配置参数
- **配置内容**:
  - 面板数据路径
  - 比较算法参数
  - 过滤阈值设定

## 脚本依赖关系 (Script Dependencies)

```
sample_qc_main.py
├── sample_qc_pipeline.py
├── sample_qc_calculator.py
└── sample_qc_flags.py

variant_qc_main.py
├── variant_qc_calculator.py
├── variant_qc_flags.py
└── variant_qc_rm_maf0_vmiss1.py

pca_qc_main.py
└── pca_qc_tools.py

bbj_projection_main.py
├── bbj_prepare_main.py
├── bbj_pca_main.py
├── bbj_projection_tools.py
└── bbj_sample_keep_main.py

panel_compare_main.py
├── panel_compare_tools.py
├── panel_thr_main.py
└── panel_filter_main.py

miss_bias_main.py
└── miss_bias_tools.py
```

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

---