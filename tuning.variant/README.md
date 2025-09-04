# CTEPH-AGP3K Variant Parameter Tuning

## 概述 (Overview)

本目录包含了 CTEPH-AGP3K 项目中用于变异参数调优的脚本和分析流程。主要用于评估和验证变异质量控制参数的设置，特别是针对 Mapping Quality (MQ) 和 Variant Quality Score Recalibration (VQSLOD) 阈值的优化。

## 目录结构 (Directory Structure)

```
tuning.variant/
├── README.md                    # 本说明文档
├── mq.vq.check.nf              # Nextflow 主流程脚本
├── variant.parameter.ipynb      # Jupyter notebook 分析文件
├── merged_output.pdf            # 合并后的结果报告
├── tmp/                         # 临时输出文件目录
│   ├── chr1.mq.vq.check.pdf    # 各染色体的分析结果
│   ├── chr2.mq.vq.check.pdf
│   ├── ...
│   └── PAR.mq.vq.check.pdf
└── work/                        # Nextflow 工作目录
```

## 主要文件说明 (Main Files)

### 核心流程文件
- **`mq.vq.check.nf`**: Nextflow 主流程脚本，用于批量分析各染色体的变异质量参数
- **`variant.parameter.ipynb`**: Jupyter notebook，用于合并 PDF 报告和进一步的数据分析
- **`merged_output.pdf`**: 合并所有染色体分析结果的综合报告

### 依赖脚本
- **`../script/variant.parameter.py`**: Python 分析脚本，负责具体的变异参数分析和可视化

## 分析流程 (Analysis Workflow)

### 1. 数据输入 (Input Data)
- **输入目录**: `/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/wgs/05.gt_norm`
- **输入文件格式**: `{chr}.pass.mac1.setid.gt_af.norm.vcf.gz`
- **染色体范围**: chr1-chr22 + PAR (伪常染色体区域)

### 2. 参数分析 (Parameter Analysis)

#### 2.1 Mapping Quality (MQ) 分析
- **默认阈值**: 58.75
- **分析内容**: 
  - MQ 值分布直方图
  - 超过阈值的变异数量和百分比
  - MQ 值范围统计

#### 2.2 Variant Quality Score Recalibration (VQSLOD) 分析
- **默认阈值**: 10
- **分析内容**:
  - VQSLOD 值分布直方图
  - 超过阈值的变异数量和百分比
  - VQSLOD 值范围统计

### 3. 可视化输出 (Visualization Output)
每个染色体生成独立的 PDF 报告，包含：
- MQ 分布图 (上图)
- VQSLOD 分布图 (下图)
- 统计信息标注
- 阈值线标记

## 核心脚本详解 (Core Scripts)

### mq.vq.check.nf
Nextflow 流程脚本，主要功能：
```nextflow
process mq_vq_check {
    // 对每个染色体进行 MQ 和 VQSLOD 分析
    // 调用 variant.parameter.py 脚本
    // 生成染色体特异性的 PDF 报告
}
```

**关键参数**:
- `--MQ 58.75`: Mapping Quality 阈值
- `--VQSLOD 10`: Variant Quality Score Recalibration 阈值
- `--chrom`: 染色体标识
- `--vcf_path`: 输入 VCF 文件路径
- `--output_pdf`: 输出 PDF 文件名

### variant.parameter.py
Python 分析脚本，主要功能：

#### 数据提取
```bash
bcftools query -f "%CHROM\t%POS\t%REF\t%ALT\t%INFO/MQ\t%INFO/VQSLOD\n" input.vcf.gz
```

#### 统计分析
- 读取 MQ 和 VQSLOD 值
- 计算分布统计
- 应用阈值过滤
- 生成百分比统计

#### 可视化生成
- 双子图布局 (2行1列)
- 直方图分布展示
- 阈值线标记
- 统计信息注释

### variant.parameter.ipynb
Jupyter notebook，主要功能：

#### PDF 合并
```python
# 自动识别所有染色体的 PDF 文件
# 按染色体顺序排序 (chr1-chr22, PAR)
# 合并为单一 PDF 报告
```

## 运行说明 (Execution Instructions)

### 环境要求
- Nextflow
- SLURM 集群环境
- Conda 环境: `cteph_geno_pro`
- 必需工具:
  - bcftools
  - Python 3.x
  - pandas, numpy, matplotlib
  - PyPDF2

### 运行命令

#### 1. 执行 Nextflow 流程
```bash
cd /LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/tuning.variant
nextflow run mq.vq.check.nf
```

#### 2. 合并 PDF 报告
```bash
# 在 Jupyter notebook 或 Python 环境中运行
python -c "exec(open('variant.parameter.ipynb').read())"
```

### 参数配置
可以通过修改 `mq.vq.check.nf` 中的参数来调整阈值：
```bash
python ${scriptDir}/variant.parameter.py \
    --chrom ${chr} \
    --vcf_path ${vcf} \
    --MQ 58.75 \        # 可调整 MQ 阈值
    --VQSLOD 10 \       # 可调整 VQSLOD 阈值
    --output_pdf ${pdf}
```

## 输出文件说明 (Output Files)

### 临时文件 (tmp/)
- **`chr*.mq.vq.check.pdf`**: 各染色体的独立分析报告
- **`chr*.MQ.VQSLOD.txt`**: 提取的原始数据文件 (临时)

### 最终输出
- **`merged_output.pdf`**: 合并所有染色体的综合报告

## 分析结果解读 (Result Interpretation)

### MQ (Mapping Quality) 评估
- **高质量阈值**: MQ > 58.75
- **评估指标**:
  - 通过阈值的变异比例
  - MQ 值分布特征
  - 极值情况检查

### VQSLOD (Variant Quality Score Recalibration) 评估
- **高质量阈值**: VQSLOD > 10
- **评估指标**:
  - 通过阈值的变异比例
  - VQSLOD 值分布特征
  - 质量得分的一致性

### 质量控制建议
1. **MQ 阈值调整**: 根据分布特征适当调整阈值
2. **VQSLOD 阈值优化**: 平衡变异保留率和质量要求
3. **染色体差异**: 关注不同染色体间的质量差异
4. **批次效应**: 检查是否存在系统性偏差

## 应用场景 (Use Cases)

### 1. 阈值优化
- 评估当前 MQ 和 VQSLOD 阈值的合理性
- 比较不同阈值设置对变异保留率的影响
- 为下游分析确定最优参数

### 2. 质量评估
- 检查 VQSR 处理后的数据质量
- 识别潜在的技术问题或批次效应
- 验证数据预处理的有效性

### 3. 报告生成
- 为项目团队提供可视化的质量评估报告
- 支持数据质量的技术文档编制
- 辅助同行评议和发表准备

## 技术细节 (Technical Details)

### 计算资源配置
```nextflow
process mq_vq_check {
    executor 'slurm'
    queue 'gr10478b'
    time '36h'
    conda '/home/b/b37974/anaconda3/envs/cteph_geno_pro'
}
```

### 并行处理
- 染色体级别并行处理
- 每个染色体独立分析
- 减少内存占用和计算时间

### 数据处理流程
1. **数据提取**: 使用 bcftools 提取 INFO 字段
2. **统计计算**: pandas 进行数据处理和统计
3. **可视化**: matplotlib 生成高质量图表
4. **报告合并**: PyPDF2 整合多个报告

## 故障排除 (Troubleshooting)

### 常见问题
1. **内存不足**: 调整进程资源配置
2. **文件权限**: 检查输入文件的访问权限
3. **依赖缺失**: 确保所有 Python 包已安装
4. **路径错误**: 验证输入文件路径的正确性

### 日志检查
- `.nextflow.log`: Nextflow 执行日志
- `work/`: 详细的进程执行日志

## 版本信息 (Version Info)

- **创建日期**: 2025年9月
- **适用版本**: CTEPH-AGP3K v1.0
- **兼容性**: Nextflow 20.0+, Python 3.7+

## 联系信息 (Contact)

- **开发者**: ZHAO TIE
- **项目**: CTEPH-AGP3K
- **更新**: 此项目为就参数调优项目，不再维护

---

此变异参数调优模块为 CTEPH-AGP3K 项目提供了重要的质量控制评估功能，确保下游分析使用高质量的变异数据。