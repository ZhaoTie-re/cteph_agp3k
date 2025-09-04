# CTEPH-AGP3K 基因型一致性调优分析流程 (Tuning Concordance Rev1)

## 概述 (Overview)

本模块是 CTEPH-AGP3K 项目中**最核心的质量控制和参数优化**流程之一，专门用于评估和优化 WGS (全基因组测序) 与 Array (芯片) 数据之间的基因型一致性。通过系统性地测试不同的质量控制参数组合，为项目提供最优的基因型过滤标准，确保下游分析的准确性和可靠性。

**🔥 重要性说明**: 这是一个持续更新和优化的核心流程，直接影响整个项目的数据质量和分析结果的可信度。

## 项目背景与意义 (Background & Significance)

### 科学问题
在大规模基因组学研究中，经常需要整合多个平台的数据：
- **WGS 数据**: 提供全基因组覆盖，但成本较高，样本量相对较小
- **Array 数据**: 覆盖已知重要位点，成本较低，样本量较大
- **数据整合挑战**: 不同平台之间存在系统性差异，需要严格的质量控制

### 解决方案
通过 **WGS vs Array 基因型一致性分析**：
1. **验证数据质量**: 确保两个平台在重叠位点上的基因型调用准确性
2. **优化过滤参数**: 找到最佳的质量控制阈值组合
3. **量化系统误差**: 识别和纠正平台间的系统性偏差
4. **建立质控标准**: 为项目制定统一的数据质量标准

## 目录结构 (Directory Structure)

```
tuning.concordance.rev1/
├── README.md                                            # 本说明文档
├── tuning.concordance.rev1.nf                           # Nextflow 主流程脚本
├── nextflow.config                                      # 计算资源配置文件
├── scripts/                                             # 核心分析脚本目录
│   ├── extract_vcf_format.py                            # VCF FORMAT字段提取
│   ├── evaluate_genotype_concordance_and_vmiss.rev3.py  # 一致性评估核心算法
│   ├── genotype_concordance_vmiss.summary.rev1.py       # 结果汇总分析
│   ├── merge_concordance_results.rev1.py                # 染色体结果合并
│   └── analysis_vis/                                    # 可视化分析脚本
│       ├── vis.tuning.gt.auto.chr.ipynb                 # 参数调优可视化
│       ├── vis.cross.platform.auto.chr.ipynb            # 跨平台比较可视化
│       ├── vis.auto.chr.ipynb                           # 自动化分析可视化
│       ├── check.ipynb                                  # 结果验证脚本
│       ├── merged_summary_15x.csv                       # 15x深度汇总结果
│       ├── merged_summary_30x.csv                       # 30x深度汇总结果
│       └── merged_summary_all.csv                       # 全体样本汇总结果
├── results/                                             # 分析结果目录
│   ├── 01.prepare_format_matrix/                        # 数据预处理结果
│   ├── 02.concordance_vmiss_calculator/                 # 一致性计算结果
│   ├── 03.concordance_vmiss_summary/                    # 染色体级别汇总
│   ├── 04.merge_concordance_vmiss/                      # 跨染色体合并结果
│   ├── 05.merge_concordance_vmiss_summary/              # 最终汇总结果
│   └── tmp/                                             # 临时文件
└── work/                                                # Nextflow 工作目录
```

## 核心分析流程 (Core Analysis Workflow)

### 第一阶段: 数据预处理 (Data Preprocessing)

#### Process 1: FormatMatrixPrepare
**功能**: 从 WGS 和 Array VCF 文件中提取共享样本和变异的 FORMAT 字段矩阵

**输入数据**:
- **WGS VCF**: `wgs/06.variant_filter/{chr}.pass.mac1.vfilter.vcf.gz`
- **Array VCF**: `array/07.z3.re.select.ids.vcf/cteph_agp3k.ajsa.qc.rechr.norm.reselect.ids.vcf.gz`

**核心算法** (`extract_vcf_format.py`):
```python
# 1. 并行提取两个VCF文件的变异ID和样本ID
# 2. 计算共享变异和共享样本集合
# 3. 使用bcftools批量提取FORMAT字段矩阵
# 4. 输出压缩的TSV格式矩阵文件
```

**输出文件**:
```
{chr}.shared_samples.txt           # 共享样本列表
{chr}.shared_variant_ids.txt       # 共享变异位点列表
{chr}.wgs.GT.tsv.gz               # WGS基因型矩阵
{chr}.wgs.DP.tsv.gz               # WGS测序深度矩阵
{chr}.wgs.GQ.tsv.gz               # WGS基因型质量矩阵
{chr}.wgs.AF.tsv.gz               # WGS等位基因频率矩阵
{chr}.array.GT.tsv.gz             # Array基因型矩阵(真值)
```

**技术特点**:
- **内存优化**: 分块处理大型矩阵，避免内存溢出
- **并行加速**: 多线程并行提取，显著提升处理速度
- **格式标准化**: 统一的矩阵格式便于后续分析

### 第二阶段: 参数网格搜索 (Parameter Grid Search)

#### 参数空间定义
我们系统性地测试以下参数组合：

```nextflow
dp_channel = Channel.from(1..30)           # 测序深度阈值: 1-30
gq_channel = Channel.of(10, 20, 30)        # 基因型质量阈值: 10, 20, 30
laf_channel = Channel.of(0.0, 0.1, 0.15, 0.2, 0.25, 0.3)  # 低AF阈值
haf_channel = Channel.of(1.0, 0.9, 0.85, 0.8, 0.75, 0.7)  # 高AF阈值

# 总参数组合数: 30 × 3 × 6 × 6 = 3,240 种组合
```

**参数含义**:
- **DP (Depth)**: 测序深度阈值，低于此值的基因型被标记为缺失
- **GQ (Genotype Quality)**: 基因型质量阈值，低于此值的基因型被标记为缺失
- **LAF (Low Allele Fraction)**: 低等位基因频率阈值，用于过滤杂合子中的偏斜等位基因比例
- **HAF (High Allele Fraction)**: 高等位基因频率阈值，与LAF配合使用

#### Process 2: ConcordanceVmissCalculator
**功能**: 对每个参数组合和每个染色体计算基因型一致性和缺失率

**核心算法** (`evaluate_genotype_concordance_and_vmiss.rev3.py`):

```python
# 主要分析步骤:
1. 加载WGS和Array的GT、DP、GQ、AF矩阵
2. 根据当前参数组合应用质量过滤:
   - DP < threshold → 标记为缺失
   - GQ < threshold → 标记为缺失  
   - LAF < AF < HAF (杂合子) → 标记为缺失
3. 计算一致性指标:
   - 混淆矩阵 (Confusion Matrix)
   - 一致性率 (Concordance Rate) 
   - 非参考一致性 (NRC, Non-Reference Concordance)
4. 计算缺失率指标:
   - VMISS: 每个变异位点的样本缺失率
   - SMISS: 每个样本的变异缺失率
5. 按样本测序深度分组分析 (15x vs 30x)
```

**输出结果**:
```python
# 保存为pickle文件，包含:
results_dict = {
    parameter_combination: (vmiss_df, smiss_df, confusion_df)
}
```

#### Process 3: ConcordanceVmissSummary  
**功能**: 对每个参数组合的结果进行统计汇总

**汇总指标**:
```python
# 计算各种汇总统计:
- 总体一致性率
- 分基因型一致性率 (0/0, 0/1, 1/1)
- 平均VMISS和SMISS
- 各测序深度组的表现对比
- 敏感性和特异性指标
```

### 第三阶段: 结果整合与分析 (Results Integration & Analysis)

#### Process 4: MergeConcordanceVmiss
**功能**: 将同一参数组合下所有染色体的结果进行合并

```python
# 合并策略:
1. 按参数组合(DP, GQ, LAF, HAF)分组
2. 将chr1-chr22的结果文件合并
3. 计算全基因组水平的统计指标
4. 生成参数组合的综合评估报告
```

#### Process 5: MergeConcordanceVmissSummary
**功能**: 生成最终的参数优化报告

**最终输出**:
- **全基因组一致性矩阵**: 所有参数组合的表现评估
- **最优参数推荐**: 基于多个指标的最优参数组合
- **平台差异分析**: WGS vs Array的系统性差异评估

## 关键Python脚本详解 (Key Scripts Analysis)

### 1. extract_vcf_format.py
**作用**: VCF文件FORMAT字段的高效提取工具

**核心功能**:
```python
def extract_shared_data():
    """
    1. 并行提取两个VCF的变异ID和样本ID
    2. 计算交集得到共享变异和样本
    3. 使用bcftools query高效提取FORMAT字段
    4. 输出标准化的矩阵格式
    """
```

**技术亮点**:
- **内存优化**: 流式处理，不一次性加载全部数据
- **并行加速**: ThreadPoolExecutor并行处理
- **错误处理**: 完善的异常处理和数据验证

### 2. evaluate_genotype_concordance_and_vmiss.rev3.py
**作用**: 一致性评估的核心算法引擎

**核心算法**:
```python
def evaluate_concordance_batch(batch_data, filters):
    """
    批量处理基因型一致性评估
    
    参数:
    - batch_data: 基因型、DP、GQ、AF数据批次
    - filters: 质量过滤参数(DP, GQ, LAF, HAF)
    
    返回:
    - confusion_matrix: 混淆矩阵
    - vmiss_stats: 变异缺失统计
    - smiss_stats: 样本缺失统计
    """
    
    # 1. 应用质量过滤
    wgs_filtered = apply_quality_filters(wgs_data, filters)
    
    # 2. 计算一致性
    confusion = calculate_confusion_matrix(wgs_filtered, array_data)
    
    # 3. 计算缺失率
    vmiss = calculate_variant_missing_rate(wgs_filtered)
    smiss = calculate_sample_missing_rate(wgs_filtered)
    
    return confusion, vmiss, smiss
```

**算法优势**:
- **可扩展性**: 支持新增质量控制指标
- **高效性**: 向量化计算，处理速度快
- **稳健性**: 处理各种边界情况和异常数据

### 3. genotype_concordance_vmiss.summary.rev1.py
**作用**: 结果汇总和统计分析

**核心指标计算**:
```python
def calculate_summary_metrics(results_dict):
    """
    计算综合评估指标
    
    指标包括:
    1. 总体一致性率 (Overall Concordance Rate)
    2. 非参考一致性率 (Non-Reference Concordance)  
    3. 敏感性 (Sensitivity)
    4. 特异性 (Specificity)
    5. F1分数 (F1 Score)
    6. 平均VMISS/SMISS
    """
```

## 可视化分析模块 (Visualization Analysis)

### 核心可视化脚本

#### 1. vis.tuning.gt.auto.chr.ipynb
**功能**: 参数调优的交互式可视化分析

**主要图表**:
```python
# 1. 参数热图 (Parameter Heatmaps)
# - DP vs GQ的一致性热图
# - LAF vs HAF的一致性热图
# - 多维参数的综合评估

# 2. 趋势分析图 (Trend Analysis)
# - 一致性随参数变化的趋势
# - 数据保留率 vs 一致性的权衡
# - 不同测序深度的表现对比

# 3. 最优参数识别
# - Pareto最优解可视化
# - 多目标优化结果展示
```

#### 2. vis.cross.platform.auto.chr.ipynb  
**功能**: 跨平台数据质量比较分析

**比较维度**:
```python
# 1. 平台系统差异分析
# - WGS vs Array的基因型分布差异
# - 不同MAF区间的平台表现
# - 染色体特异性差异模式

# 2. 质量控制效果评估
# - 过滤前后的数据质量对比
# - 不同参数设置的效果展示
# - 质量改善的量化分析
```

#### 3. vis.auto.chr.ipynb
**功能**: 自动化分析报告生成

**自动化内容**:
```python
# 1. 质量控制报告
# - 数据质量评估总结
# - 异常样本/变异识别
# - 推荐的质控参数

# 2. 一致性分析报告  
# - 平台间一致性评估
# - 问题位点和样本标记
# - 数据整合建议
```

## 重要输出文件解读 (Key Output Interpretation)

### 1. 汇总统计文件
```csv
# merged_summary_all.csv 示例结构
DP,GQ,LAF,HAF,CHR,PLATFORM,TOTAL_VARIANTS,CONCORDANT_VARIANTS,CONCORDANCE_RATE,NRC,VMISS_MEAN,SMISS_MEAN
8,20,0.25,0.75,merged,ALL,1500000,1425000,0.95,0.92,0.08,0.05
...
```

**字段含义**:
- **CONCORDANCE_RATE**: 总体一致性率
- **NRC**: 非参考基因型一致性率  
- **VMISS_MEAN**: 平均变异缺失率
- **SMISS_MEAN**: 平均样本缺失率

### 2. 混淆矩阵文件
```python
# confusion_matrix 结构
{
    'overall': {
        '0/0_0/0': count,  # 双平台都是纯合子参考型
        '0/0_0/1': count,  # WGS参考型，Array杂合子
        '0/0_1/1': count,  # WGS参考型，Array纯合子变异型
        # ... 其他组合
    },
    'by_platform': {
        '15x': {...},  # 15x测序深度样本的混淆矩阵
        '30x': {...}   # 30x测序深度样本的混淆矩阵
    }
}
```

## 参数优化策略 (Parameter Optimization Strategy)

### 多目标优化框架

我们的参数优化考虑以下目标：

```python
# 目标函数设计
objectives = {
    'maximize': [
        'concordance_rate',      # 最大化一致性率
        'data_retention_rate',   # 最大化数据保留率
        'sensitivity',           # 最大化敏感性
        'specificity'           # 最大化特异性
    ],
    'minimize': [
        'vmiss_rate',           # 最小化变异缺失率
        'smiss_rate',           # 最小化样本缺失率
        'false_positive_rate'   # 最小化假阳性率
    ]
}
```

### Pareto最优解选择

```python
def find_pareto_optimal_parameters():
    """
    基于多目标优化寻找Pareto最优参数组合
    
    考虑因素:
    1. 一致性率 vs 数据保留率的权衡
    2. 不同测序深度样本的平衡表现
    3. 计算复杂度和实际应用的可行性
    """
```

### 推荐参数设置

基于当前分析结果，我们的初步推荐：

```python
# 保守设置 (高质量，低保留率)
conservative_params = {
    'DP': 15,
    'GQ': 30, 
    'LAF': 0.25,
    'HAF': 0.75
}

# 平衡设置 (中等质量，中等保留率)
balanced_params = {
    'DP': 10,
    'GQ': 20,
    'LAF': 0.2, 
    'HAF': 0.8
}

# 宽松设置 (保留更多数据)
liberal_params = {
    'DP': 6,
    'GQ': 10,
    'LAF': 0.15,
    'HAF': 0.85
}
```

## 运行说明 (Execution Instructions)

### 环境要求

```bash
# 必需软件环境
- Nextflow (≥ 20.0)
- SLURM 集群环境
- Conda 环境: cteph_geno_pro

# Python依赖包
- pandas ≥ 1.3.0
- numpy ≥ 1.20.0
- pickle (标准库)
- multiprocessing (标准库)
- concurrent.futures (标准库)

# 外部工具
- bcftools ≥ 1.10
```

### 运行步骤

#### 1. 数据准备验证
```bash
# 验证输入文件存在性
ls ${params.wgsDir}/*.vcf.gz
ls ${params.arrayDir}/*.vcf.gz
ls ${params.infoDir}/wgs_array_dp.csv
```

#### 2. 启动完整流程
```bash
cd /LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/tuning.concordance.rev1

# 测试运行 (仅chr1和chr22)
nextflow run tuning.concordance.rev1.nf -c nextflow.config --chr_subset "chr1,chr22"

# 完整运行 (全部染色体)
nextflow run tuning.concordance.rev1.nf -c nextflow.config
```

#### 3. 监控运行状态
```bash
# 检查Nextflow日志
tail -f .nextflow.log

# 查看进程状态
nextflow log

# 检查结果输出
ls -la results/*/
```

#### 4. 结果分析
```bash
# 启动Jupyter进行可视化分析
jupyter notebook scripts/analysis_vis/
```

### 资源配置优化

```nextflow
// nextflow.config中的资源分配
process {
    withName: 'FormatMatrixPrepare' {
        clusterOptions = '--rsc p=1:t=8:c=4:m=18284M'  // 内存密集型
    }
    withName: 'ConcordanceVmissCalculator' {
        clusterOptions = '--rsc p=1:t=4:c=2:m=9142M'   // 计算密集型
    }
    // ... 其他进程配置
}
```

## 质量控制与验证 (Quality Control & Validation)

### 数据完整性检查

```python
# 自动化质量检查脚本
def validate_analysis_results():
    """
    验证分析结果的完整性和一致性
    
    检查项目:
    1. 所有参数组合是否都有结果
    2. 染色体结果是否完整
    3. 统计指标是否在合理范围内
    4. 数据格式是否标准化
    """
    
    # 检查结果文件完整性
    check_file_completeness()
    
    # 验证统计指标合理性
    validate_statistical_metrics()
    
    # 检查数据格式一致性
    verify_data_formats()
```

### 结果可信度评估

```python
# 可信度评估指标
reliability_metrics = {
    'sample_size_adequacy': 'check_minimum_sample_size()',
    'variant_count_sufficiency': 'check_variant_coverage()', 
    'technical_replicates_consistency': 'check_replicate_concordance()',
    'batch_effect_assessment': 'evaluate_batch_effects()'
}
```

## 故障排除指南 (Troubleshooting Guide)

### 常见问题与解决方案

#### 1. 内存不足错误
```bash
# 问题: Java heap space / Out of memory
# 解决: 增加内存分配
export NXF_OPTS='-Xms2g -Xmx8g'

# 或修改nextflow.config
process.memory = '16 GB'
```

#### 2. 文件权限问题
```bash
# 问题: Permission denied
# 解决: 检查文件权限
chmod 755 scripts/*.py
chmod 644 input_files/*.vcf.gz
```

#### 3. bcftools相关错误
```bash
# 问题: bcftools命令失败
# 解决: 验证VCF文件完整性
bcftools index -f input.vcf.gz
bcftools view -H input.vcf.gz | wc -l
```

#### 4. Python模块导入错误
```bash
# 问题: ModuleNotFoundError
# 解决: 激活正确的conda环境
conda activate cteph_geno_pro
conda list | grep pandas
```

### 性能优化建议

```python
# 1. 并行度调优
optimal_threads = min(available_cores, data_size_dependent_threads)

# 2. 内存使用优化
batch_size = calculate_optimal_batch_size(available_memory, data_size)

# 3. I/O优化
use_compression = True  # 启用压缩减少磁盘I/O
cache_intermediate_results = True  # 缓存中间结果
```

## 扩展与定制 (Extensions & Customization)

### 新增质量控制指标

```python
# 扩展质量控制参数示例
def add_new_qc_metric():
    """
    添加新的质量控制指标
    
    可扩展指标:
    1. PL (Phred-scaled likelihoods)
    2. AD (Allelic depth)
    3. VAF (Variant allele fraction)
    4. 自定义复合指标
    """
    
    # 在evaluate_genotype_concordance_and_vmiss.py中添加
    new_filters = {
        'PL_threshold': 50,
        'AD_ratio_threshold': 0.2,
        'VAF_deviation_threshold': 0.1
    }
```

### 适配其他项目

```python
# 项目适配指南
def adapt_to_new_project():
    """
    适配新项目的步骤:
    
    1. 修改输入文件路径
    2. 调整参数搜索空间
    3. 自定义评估指标
    4. 更新可视化模板
    """
    
    # 配置文件模板
    project_config = {
        'input_paths': {...},
        'parameter_space': {...},
        'evaluation_metrics': {...},
        'output_formats': {...}
    }
```

## 持续更新计划 (Continuous Update Plan)

### 当前版本 (Rev1) 的主要特性
- ✅ 基础一致性评估算法
- ✅ 多参数网格搜索
- ✅ 分测序深度分析
- ✅ 基础可视化功能

### 计划中的更新 (Rev2)
- 🔄 **算法优化**: 更高效的一致性计算算法
- 🔄 **参数空间扩展**: 增加更多质量控制维度
- 🔄 **机器学习集成**: 基于ML的参数优化
- 🔄 **实时监控**: 分析进度的实时可视化

## 技术文档与引用 (Documentation & Citation)

### 方法学说明
```
本分析流程基于以下核心原理:
1. 基于共享变异位点的直接基因型比较
2. 多维质量控制参数的网格搜索优化
3. 分层分析策略(按测序深度分组)
4. 综合评估指标的多目标优化
```

### 引用格式
```bibtex
@misc{cteph_concordance_tuning,
  title={CTEPH-AGP3K基因型一致性调优分析流程},
  author={ZHAO TIE and CTEPH-AGP3K Consortium},
  year={2025},
  url={https://github.com/ZhaoTie-re/cteph_agp3k},
  note={Version Rev1, 持续更新中}
}
```

## 联系与支持 (Contact & Support)

### 开发团队
- **开发者**: ZHAO TIE

### 更新日志
- **2025-06**: Rev1版本发布，基础功能完成
- **2025-09**: 持续优化中，添加新的分析模块

---

**⚠️ 重要提醒**: 
1. 这是一个**活跃开发中的核心模块**，请定期检查更新
2. 在生产环境使用前，请务必在测试数据上验证结果
3. 任何问题或建议请及时反馈给开发团队
4. 建议在重要分析前备份当前版本和参数设置

**🎯 项目目标**: 为CTEPH-AGP3K项目提供最可靠、最优化的基因型质量控制标准，确保研究结果的科学性和可重现性。