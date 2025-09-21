# CTEPH AGP3K 蛋白质关联分析

## 项目概述

本项目是慢性血栓栓塞性肺动脉高压（CTEPH）研究中的蛋白质关联分析模块，主要用于分析基因变异与蛋白质表达水平之间的关联关系。该分析是CTEPH AGP3K项目审查分析的重要组成部分。

## 项目结构

```
01.pro_assoc/
├── README.md                           # 本文档
├── pro_assoc_main.ipynb               # 主要分析脚本（Jupyter Notebook）
├── gt_mx.tsv                          # 基因型矩阵文件
├── gt_mx.tsv.summary.json             # 基因型矩阵摘要信息
├── gt_counts_summary.tsv              # 基因型计数统计
├── variant_meta.tsv                   # 变体元数据表
├── protein_genotype_boxplots.ADD.pdf  # 加性模型箱线图
├── protein_genotype_boxplots.DOM.pdf  # 显性模型箱线图
├── protein_genotype_boxplots.REC.pdf  # 隐性模型箱线图
└── *.log                              # 各种日志文件
```

## 分析流程

### 1. 数据准备
- **蛋白质表达数据**: 来源于Soma扫描技术的蛋白质组学数据
- **基因型数据**: 来源于全基因组测序（WGS）数据，经过质控和过滤
- **GWAS结果**: 包含显著关联变异的汇总统计结果

### 2. 主要分析步骤

#### 2.1 变异筛选
- 从GWAS结果中提取包含'ADD'（加性模型）的显著变异
- 过滤得到目标变异ID列表

#### 2.2 样本ID处理
- 从蛋白质表达数据中提取样本ID
- 清理样本ID格式（去除`_day*`后缀）
- 确保样本ID与基因型数据匹配

#### 2.3 基因型数据提取
- 使用PLINK2从二进制格式文件中提取目标变异的基因型信息
- 生成基因型矩阵（`gt_mx.tsv`）
- 统计基因型分布情况

#### 2.4 元数据整合
- 整合变异信息、基因注释和蛋白质信息
- 生成综合的变体元数据表（`variant_meta.tsv`）

#### 2.5 关联分析可视化
- 针对每个变异，生成蛋白质表达水平的箱线图
- 分别使用三种遗传模型：加性（ADD）、显性（DOM）、隐性（REC）
- 对蛋白质表达数据进行log2转换以改善分布

## 输出文件说明

### 数据文件
- **gt_mx.tsv**: 基因型矩阵，行为样本，列为变异，值为基因型编码（0/1/2）
- **gt_counts_summary.tsv**: 每个变异的基因型计数统计
  - `n0/n1/n2`: 基因型为0/1/2的样本数
  - `nmiss`: 缺失基因型的样本数
  - `ncalled/nsamples`: 成功调用/总样本数
- **variant_meta.tsv**: 变异元数据表
  - `ID`: 变异ID
  - `rsID`: dbSNP ID
  - `Gene`: 相关基因
  - `SeqID`: 蛋白质序列ID
  - `UniProt`: UniProt蛋白质ID

### 可视化文件
- **protein_genotype_boxplots.*.pdf**: 蛋白质表达箱线图
  - 每个变异对应一页图表
  - 显示不同基因型组别的蛋白质表达分布
  - 包含统计显著性检验结果

## 技术依赖

### 软件工具
- **PLINK2**: 基因型数据处理
- **bcftools**: VCF文件操作
- **Python 3.x**: 主要分析环境

### Python包
- `pandas`: 数据处理和分析
- `numpy`: 数值计算
- `matplotlib`: 基础绘图
- `seaborn`: 统计图表
- `scipy`: 统计分析

## 使用说明

### 1. 环境准备
确保已安装所需的软件工具和Python包。

### 2. 运行分析
打开`pro_assoc_main.ipynb`，按顺序执行各个cell：

```python
# 1. 导入模块和设置路径
import pro_assoc_tools

# 2. 设置数据路径
pro_ex_path = "蛋白质表达数据路径"
plink_prefix = "PLINK文件前缀"
gwas_path = "GWAS结果路径"

# 3. 执行分析流程
# - 变异筛选
# - 基因型提取
# - 数据整合
# - 可视化生成
```

### 3. 结果解读
- 查看生成的PDF文件了解蛋白质表达模式
- 检查统计摘要了解数据质量
- 分析元数据表确认变异-基因-蛋白质映射关系

## 重要参数

### 基因型提取参数
- `plink_threads`: PLINK运行线程数（默认8）
- `variant_ids`: 目标变异ID列表
- `sample_ids`: 目标样本ID列表

### 可视化参数
- `log2_transform`: 是否对蛋白质表达进行log2转换（推荐True）
- `model`: 遗传模型（ADD/DOM/REC）
- `preview_n`: 预览显示的记录数（默认10）

## 注意事项

1. **数据匹配**: 确保蛋白质数据、基因型数据和GWAS结果中的样本ID一致
2. **内存使用**: 基因型矩阵可能较大，注意内存使用情况
3. **文件路径**: 所有路径必须使用绝对路径
4. **质控检查**: 运行前检查输入数据的质量和完整性

## 联系信息

*作者: ZHAO TIE*

---
*最后更新: 2025年9月21日*