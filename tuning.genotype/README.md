# CTEPH-AGP3K Genotype Quality Parameter Tuning

## 概述 (Overview)

本目录包含了 CTEPH-AGP3K 项目中用于基因型质量参数调优的脚本和分析流程。主要用于评估和优化基因型质量控制参数，特别是针对 Depth of Coverage (DP) 和 Genotype Quality (GQ) 阈值的设定，以确保下游分析使用高质量的基因型数据。

## 目录结构 (Directory Structure)

```
tuning.genotype/
├── README.md                    # 本说明文档
├── dp.gq.check.nf               # Nextflow 主流程脚本
├── nextflow.config              # Nextflow 配置文件
├── dp.parameter.ipynb           # DP 参数分析 Jupyter notebook
├── gq.parameter.ipynb           # GQ 参数分析 Jupyter notebook  
├── sample.dp.check.ipynb        # 样本信息准备 notebook
├── chr21.22.PAR.dp.check.pdf    # 合并的 DP 分析报告
├── chr21.22.PAR.gq.check.pdf    # 合并的 GQ 分析报告
├── tmp/                         # 临时输出文件目录
│   ├── chr21.dp.check.csv       # DP 统计结果文件
│   ├── chr21.dp.check.pdf       # DP 可视化报告
│   ├── chr21.gq.check.csv       # GQ 统计结果文件
│   ├── chr21.gq.check.pdf       # GQ 可视化报告
│   ├── chr22.*                  # chr22 对应文件
│   └── PAR.*                    # PAR 区域对应文件
└── work/                        # Nextflow 工作目录
```

## 主要功能 (Main Functions)

### 1. 深度覆盖度 (DP) 参数调优
- 评估不同 DP 阈值下的变异保留情况
- 分析 30x 和 15x 测序深度样本的表现差异
- 生成 DP 阈值优化建议

### 2. 基因型质量 (GQ) 参数调优  
- 评估不同 GQ 阈值下的基因型质量分布
- 分析不同测序深度样本的 GQ 表现
- 提供 GQ 过滤标准优化方案

### 3. 可视化分析
- 生成阈值-变异保留数量曲线图
- 比较不同样本组的参数表现
- 提供直观的参数选择依据

## 核心文件详解 (Core Files)

### 流程控制文件

#### dp.gq.check.nf
Nextflow 主流程脚本，包含四个核心进程：

```nextflow
process dp_check {
    // 计算不同 DP 阈值下的变异保留数量
    // 分别针对 30x 和 15x 样本组
}

process gq_check {
    // 计算不同 GQ 阈值下的变异保留数量
    // 分别针对 30x 和 15x 样本组
}

process dp_check_vis {
    // 生成 DP 参数的可视化报告
}

process gq_check_vis {
    // 生成 GQ 参数的可视化报告
}
```

**关键参数设置**:
- **输入目录**: `params.vcfDir = '.../wgs/02.mac1'`
- **分析范围**: chr21, chr22, PAR (用于快速测试)
- **脚本目录**: `params.scriptDir = '.../script'`
- **输出目录**: `params.outdir = '.../tuning.genotype'`

#### nextflow.config
配置文件，定义计算资源分配：
```nextflow-config
process {
    withName: 'dp_check' {
        clusterOptions = '--rsc p=1:t=2:c=1:m=4571M'
    }
    withName: 'gq_check' {
        clusterOptions = '--rsc p=1:t=2:c=1:m=4571M'
    }
}
```

### 核心分析脚本

#### dp.parameter.py
深度覆盖度参数分析脚本：

**功能特性**:
- **阈值范围**: DP 1-30 的全面评估
- **并行处理**: 使用 ProcessPoolExecutor 提高效率
- **样本分组**: 30x 和 15x 测序深度样本分别分析
- **输出格式**: CSV 格式的统计结果

**核心算法**:
```python
def process_vcf(vcf_path, dp_thresholds, sample_codes):
    # 遍历 VCF 记录
    # 检查指定样本的 DP 值
    # 统计满足阈值的变异数量
```

#### gq.parameter.py
基因型质量参数分析脚本：

**功能特性**:
- **阈值范围**: GQ 值的全面评估
- **并行处理**: 高效的多进程分析
- **样本分组**: 30x 和 15x 样本独立分析
- **质量评估**: 基因型调用置信度分析

#### 可视化脚本

#### dp.parameter.vis.py
DP 参数可视化脚本：

**可视化元素**:
- **线性图**: DP 阈值 vs 保留变异数量
- **参考线**: 总变异数量水平线
- **分组比较**: 30x vs 15x 样本表现
- **网格线**: 便于读数的辅助线

#### gq.parameter.vis.rev3.py
GQ 参数可视化脚本：

**图表特性**:
- **趋势分析**: GQ 阈值变化对数据保留的影响
- **质量评估**: 不同 GQ 阈值的适用性
- **比较分析**: 多样本组的表现对比

## 分析流程 (Analysis Workflow)

### 阶段 1: 数据准备 (Data Preparation)

#### 1.1 样本信息整理
```python
# sample.dp.check.ipynb 中的核心步骤
- 读取 WGS 样本列表
- 筛选 CTEPH 样本 (PHOM 前缀)
- 匹配 JHRP4 面板信息
- 生成样本元数据文件
```

#### 1.2 输入数据验证
- **VCF 文件**: `chr*.pass.mac1.vcf.gz`
- **样本分组**: 30x 和 15x 测序深度
- **分析范围**: chr21, chr22, PAR (快速测试)

### 阶段 2: 参数分析 (Parameter Analysis)

#### 2.1 DP 参数评估
```bash
# 对每个染色体执行 DP 分析
python dp.parameter.py \
    --chrom chr21 \
    --vcf_path chr21.pass.mac1.vcf.gz \
    --metadata cteph_jhrp4_info.csv \
    --output_csv chr21.dp.check.csv
```

**分析维度**:
- **阈值范围**: DP 1-30
- **样本分组**: 30x vs 15x
- **统计指标**: 每个阈值下的变异保留数量

#### 2.2 GQ 参数评估
```bash
# 对每个染色体执行 GQ 分析
python gq.parameter.py \
    --chrom chr21 \
    --vcf_path chr21.pass.mac1.vcf.gz \
    --metadata cteph_jhrp4_info.csv \
    --output_csv chr21.gq.check.csv
```

**分析维度**:
- **质量评估**: 基因型调用置信度
- **阈值优化**: 最佳 GQ 切点确定
- **样本比较**: 不同测序深度的表现

### 阶段 3: 可视化报告 (Visualization)

#### 3.1 DP 可视化
```bash
python dp.parameter.vis.py \
    --chrom chr21 \
    --count_file chr21.count.txt \
    --dp_summary chr21.dp.check.csv \
    --output_pdf chr21.dp.check.pdf
```

#### 3.2 GQ 可视化
```bash
python gq.parameter.vis.rev3.py \
    --chrom chr21 \
    --count_file chr21.count.txt \
    --gq_summary chr21.gq.check.csv \
    --output_pdf chr21.gq.check.pdf
```

## 运行说明 (Execution Instructions)

### 环境要求
- **Nextflow**: 工作流管理
- **SLURM**: 集群任务调度
- **Conda 环境**: `cteph_geno_pro`
- **必需工具**:
  - pysam (VCF 文件处理)
  - pandas (数据分析)
  - matplotlib (可视化)
  - numpy (数值计算)

### 运行步骤

#### 1. 准备样本信息
```bash
# 在 Jupyter 环境中运行
jupyter notebook sample.dp.check.ipynb
```

#### 2. 执行参数分析
```bash
cd /LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/tuning.genotype
nextflow run dp.gq.check.nf -c nextflow.config
```

#### 3. 查看结果
```bash
# 检查临时文件
ls tmp/

# 查看合并报告
open chr21.22.PAR.dp.check.pdf
open chr21.22.PAR.gq.check.pdf
```

### 自定义参数

#### 修改染色体范围
```nextflow
// 在 dp.gq.check.nf 中修改
Channel
    .from((1..22).collect { "chr${it}" } + ["PAR"])  // 全基因组分析
    .set { chr_ch_dp }
```

#### 调整阈值范围
```python
# 在 dp.parameter.py 中修改
dp_thresholds = list(range(1, 51))  // 扩展到 DP 1-50
```

## 输出文件解读 (Output Interpretation)

### CSV 统计文件
```csv
# dp.check.csv 示例
,1,2,3,4,5,...,30
30x,15234,14832,14256,13845,13234,...,8956
15x,12456,12034,11567,11234,10845,...,7234
```

**数据含义**:
- **行**: 样本组 (30x, 15x)
- **列**: DP 阈值 (1-30)
- **值**: 满足阈值的变异数量

### PDF 可视化报告

#### DP 报告内容
1. **趋势曲线**: DP 阈值 vs 保留变异数
2. **样本比较**: 30x vs 15x 表现差异
3. **参考线**: 总变异数量基准线
4. **统计信息**: 关键数值标注

#### GQ 报告内容
1. **质量分布**: GQ 阈值影响分析
2. **保留率曲线**: 数据损失评估
3. **分组对比**: 不同测序深度比较
4. **推荐阈值**: 最优参数建议

## 参数优化建议 (Parameter Optimization)

### DP 阈值选择原则

#### 1. 平衡原则
- **数据保留**: 保持足够的变异数量
- **质量保证**: 确保测序深度可靠性
- **样本公平**: 兼顾不同测序深度样本

#### 2. 推荐设置
- **30x 样本**: DP ≥ 10
- **15x 样本**: DP ≥ 8
- **保守设置**: DP ≥ 6 (全样本通用)

### GQ 阈值选择原则

#### 1. 质量标准
- **高质量**: GQ ≥ 30 (推荐)
- **中等质量**: GQ ≥ 20 (可接受)
- **宽松过滤**: GQ ≥ 10 (特殊情况)

#### 2. 应用场景
- **GWAS 分析**: 建议 GQ ≥ 30
- **变异发现**: 可用 GQ ≥ 20
- **探索性分析**: 可用 GQ ≥ 10

## 质量控制检查 (Quality Control)

### 数据一致性检查
```bash
# 检查样本数量一致性
wc -l cteph_jhrp4_info.csv

# 验证 VCF 文件完整性
bcftools index -s input.vcf.gz
```

### 结果验证
```python
# 验证统计结果的合理性
- 检查阈值递增时数量递减
- 验证 30x > 15x 的一般规律
- 确认极值的合理性
```

## 故障排除 (Troubleshooting)

### 常见问题

#### 1. 内存不足
```bash
# 增加内存配置
clusterOptions = '--rsc p=1:t=2:c=1:m=9142M'
```

#### 2. 文件路径错误
```bash
# 检查文件存在性
ls /LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/wgs/02.mac1/
```

#### 3. 样本信息不匹配
```python
# 验证样本代码一致性
vcf_samples = set(vcf_file.header.samples)
meta_samples = set(cteph_jhrp4['sample_code'])
print(f"Common samples: {len(vcf_samples & meta_samples)}")
```

### 日志分析
```bash
# 检查 Nextflow 日志
tail -f .nextflow.log

# 查看具体进程日志
cat work/*/*/.command.log
```

## 扩展应用 (Extended Applications)

### 1. 全基因组分析
- 修改染色体范围至 chr1-22+X+Y
- 调整计算资源配置
- 增加并行度设置

### 2. 其他质量参数
- **AD** (Allelic Depth): 等位基因深度
- **AB** (Allele Balance): 等位基因平衡
- **VAF** (Variant Allele Frequency): 变异等位基因频率

### 3. 多项目适配
- 修改样本元数据格式
- 调整阈值范围
- 自定义可视化样式

## 技术细节 (Technical Details)

### 并行处理优化
```python
# 使用多进程提高效率
with ProcessPoolExecutor(max_workers=2) as executor:
    futures = [executor.submit(process_vcf, vcf_path, thresholds, samples) 
               for samples in sample_sets]
```

### 内存管理
```python
# 逐条记录处理避免内存溢出
for record in vcf_in.fetch():
    # 处理单个变异记录
    # 避免全部载入内存
```

### 数据结构优化
```python
# 使用字典和集合提高查询效率
sample_codes = set(sample_list)  # O(1) 查询
counts = {threshold: 0 for threshold in thresholds}  # 高效计数
```

## 版本信息 (Version Info)

- **创建日期**: 2025年9月
- **适用版本**: CTEPH-AGP3K v1.0
- **依赖版本**:
  - Nextflow: 20.0+
  - Python: 3.7+
  - pysam: 0.16+
  - pandas: 1.3+
  - matplotlib: 3.5+

## 联系信息 (Contact)

- **开发者**: ZHAO TIE
- **项目**: CTEPH-AGP3K
- **更新**: 此项目为就参数调优项目，不再维护

---

此基因型质量参数调优模块为 CTEPH-AGP3K 项目提供了科学的参数选择依据，确保后续分析使用最优的质量控制标准。通过系统性的阈值评估和可视化分析，为项目的数据质量把控提供了重要支撑。