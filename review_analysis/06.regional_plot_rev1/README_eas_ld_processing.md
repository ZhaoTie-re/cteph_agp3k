# EAS LD数据处理脚本使用说明

## 文件说明

## 计算资源限制导致`regional_plot_prepare.ipynb`针对1000G的EAS计算受限的时候，可以选用这个脚本作为备选

- `run_eas_ld_processing.sh`: SLURM批处理脚本，用于处理EAS 1000G LD数据

## 使用方法

### 1. 检查必要文件

确保以下文件存在：
```bash
cd /LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/review_analysis/06.regional_plot_rev1

# 检查LD JSON文件
ls -lh cteph_agp3k.lowfreq_common.ld_matrices_by_lead.summary.json

# 检查EAS bed文件
ls -lh eas_all.{bed,bim,fam}

# 检查交集JSON文件
ls -lh eas_ld_tmp/variant_intersection_summary.json
```

### 2. 提交任务

```bash
cd /LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/review_analysis/06.regional_plot_rev1

# 提交SLURM任务
sbatch run_eas_ld_processing.sh
```

### 3. 监控任务

```bash
# 查看任务状态
squeue -u $USER

# 查看实时输出
tail -f eas_ld_processing_<JOB_ID>.out

# 查看错误日志
tail -f eas_ld_processing_<JOB_ID>.err
```

### 4. 检查结果

任务成功完成后，会生成：
- `eas_ld/` - 包含每个lead variant的EAS LD数据
- `eas_ld/<prefix>.eas_ld_sum_stat_summary.json` - 摘要JSON文件
- `eas_result_summary.txt` - 摘要文件路径（方便后续使用）

## 资源配置

当前SLURM配置：
- **队列**: gr10478b
- **时间限制**: 168小时（7天）
- **CPU**: 8核
- **内存**: 36568M (~36GB)
- **并行任务数**: 1（顺序处理，避免内存问题）

## 常见问题

### 问题1: 进程被杀死（killed）

**原因**: 内存不足（OOM）

**解决方案**:
1. 脚本已经添加了 `--chr` 参数来限制plink2只读取特定染色体
2. 设置了 `--memory 16000` 限制plink2内存使用
3. 使用顺序处理（n_jobs=1）而非并行处理

如果仍然失败，可以尝试：
- 增加SLURM的内存申请：修改 `#SBATCH --rsc` 中的 `m=` 参数
- 进一步减少plink2的内存限制：修改代码中的 `--memory` 参数

### 问题2: 函数返回None

**原因**: 输入文件不存在或plink执行失败

**解决方案**:
1. 检查日志中的错误信息
2. 确认所有输入文件路径正确
3. 确认plink2和plink1.9可执行文件存在且有执行权限

### 问题3: 交集为空

**原因**: EAS 1000G数据与研究数据没有共同变体

**解决方案**:
1. 检查变体ID格式是否一致
2. 检查染色体编号格式（chr16 vs 16）
3. 检查参考等位基因是否匹配

## 修改记录

- 2024-12-02: 
  - 添加 `--chr` 参数解决内存不足问题
  - 改为顺序处理避免并行导致的内存峰值
  - 优化错误处理和日志输出
  - 添加详细的文件检查
