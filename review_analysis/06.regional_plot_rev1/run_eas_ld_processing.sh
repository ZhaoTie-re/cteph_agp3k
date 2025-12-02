#!/bin/bash
#SBATCH --job-name=eas_ld_process
#SBATCH --output=eas_ld_processing_%j.out
#SBATCH --error=eas_ld_processing_%j.err
#SBATCH -p gr10478b
#SBATCH -t 168:0:0
#SBATCH --rsc p=1:t=16:c=8:m=36568M

# EAS LD数据处理任务
# 用于处理 regional_plot_prepare.ipynb 中第6个cell的任务

set -e  # 遇到错误立即退出
set -u  # 使用未定义变量时报错

# 打印任务信息
echo "=========================================="
echo "EAS LD数据处理任务"
echo "=========================================="
echo "作业ID: $SLURM_JOB_ID"
echo "节点: $SLURM_NODELIST"
echo "开始时间: $(date)"
echo "工作目录: $(pwd)"
echo "=========================================="
echo ""

# 设置工作目录
WORK_DIR="/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/review_analysis/06.regional_plot_rev1"
cd "$WORK_DIR"

# 创建日志目录
mkdir -p logs

# 设置Python路径
MODULE_DIR="/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/review_analysis/script"
export PYTHONPATH="${MODULE_DIR}:${PYTHONPATH:-}"

# 输入文件路径
LD_JSON="${WORK_DIR}/cteph_agp3k.lowfreq_common.ld_matrices_by_lead.summary.json"
EAS_BED_PREFIX="${WORK_DIR}/eas_all"
INTERSECTION_JSON="${WORK_DIR}/eas_ld_tmp/variant_intersection_summary.json"

# 检查必要文件是否存在
echo "检查输入文件..."
if [ ! -f "$LD_JSON" ]; then
    echo "错误: LD JSON文件不存在: $LD_JSON"
    exit 1
fi
echo "✓ LD JSON文件: $LD_JSON"

if [ ! -f "${EAS_BED_PREFIX}.bed" ]; then
    echo "错误: EAS BED文件不存在: ${EAS_BED_PREFIX}.bed"
    exit 1
fi
echo "✓ EAS BED文件: ${EAS_BED_PREFIX}.bed"

if [ ! -f "$INTERSECTION_JSON" ]; then
    echo "错误: 交集JSON文件不存在: $INTERSECTION_JSON"
    exit 1
fi
echo "✓ 交集JSON文件: $INTERSECTION_JSON"
echo ""

# 运行Python脚本
echo "=========================================="
echo "开始处理EAS LD数据..."
echo "=========================================="
echo ""

python3 << 'EOF'
import sys
import os

# 添加模块路径
MODULE_DIR = "/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/review_analysis/script"
if MODULE_DIR not in sys.path:
    sys.path.insert(0, MODULE_DIR)

# 导入模块
import reginal_plot_tools_rev1
import importlib
importlib.reload(reginal_plot_tools_rev1)

# 设置参数（使用环境变量中的路径）
ld_json = os.getenv('LD_JSON', "cteph_agp3k.lowfreq_common.ld_matrices_by_lead.summary.json")
eas_bed_prefix = "/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/review_analysis/06.regional_plot_rev1/eas_all"
intersection_json = "/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/review_analysis/06.regional_plot_rev1/eas_ld_tmp/variant_intersection_summary.json"
output_dir = None  # 使用默认的eas_ld目录
n_jobs = 1  # 使用1个任务，避免内存问题（顺序处理）

print(f"参数设置:")
print(f"  - LD JSON: {ld_json}")
print(f"  - EAS BED前缀: {eas_bed_prefix}")
print(f"  - 交集JSON: {intersection_json}")
print(f"  - 并行任务数: {n_jobs}")
print("")

# 再次检查文件是否存在（Python层面）
print("检查输入文件（Python）...")
files_ok = True
if not os.path.exists(ld_json):
    print(f"✗ LD JSON文件不存在: {ld_json}")
    files_ok = False
else:
    print(f"✓ LD JSON文件存在: {ld_json}")

if not os.path.exists(f"{eas_bed_prefix}.bed"):
    print(f"✗ EAS BED文件不存在: {eas_bed_prefix}.bed")
    files_ok = False
else:
    print(f"✓ EAS BED文件存在: {eas_bed_prefix}.bed")

if not os.path.exists(intersection_json):
    print(f"✗ 交集JSON文件不存在: {intersection_json}")
    files_ok = False
else:
    print(f"✓ 交集JSON文件存在: {intersection_json}")

if not files_ok:
    print("")
    print("✗ 存在缺失文件，终止执行")
    sys.exit(1)

print("")

# 运行处理
try:
    eas_result_summary = reginal_plot_tools_rev1.process_ld_matrices_with_eas(
        ld_json_file=ld_json,
        eas_bed_prefix=eas_bed_prefix,
        intersection_json=intersection_json,
        output_dir=output_dir,
        n_jobs=n_jobs
    )
    
    print("")
    print("=" * 60)
    
    # 检查返回值
    if eas_result_summary is None:
        print("✗ EAS LD数据处理失败！")
        print("=" * 60)
        print("函数返回了 None，说明处理过程中遇到了错误。")
        print("请检查上面的日志信息获取详细错误原因。")
        print("")
        print("常见原因:")
        print("  1. 输入文件不存在或路径错误")
        print("  2. plink2/plink1.9 执行失败（内存不足、文件格式错误等）")
        print("  3. 变体交集为空")
        sys.exit(1)
    else:
        print("✓ EAS LD数据处理完成！")
        print("=" * 60)
        print(f"输出摘要文件: {eas_result_summary}")
        
        # 将结果路径保存到文件，方便后续使用
        with open("eas_result_summary.txt", "w") as f:
            f.write(eas_result_summary)
        print(f"结果路径已保存到: eas_result_summary.txt")
    
except Exception as e:
    print("")
    print("=" * 60)
    print("✗ 处理失败！")
    print("=" * 60)
    print(f"错误信息: {str(e)}")
    import traceback
    traceback.print_exc()
    sys.exit(1)

EOF

EXIT_CODE=$?

echo ""
echo "=========================================="
echo "任务完成"
echo "=========================================="
echo "结束时间: $(date)"
echo "退出码: $EXIT_CODE"
echo ""

if [ $EXIT_CODE -eq 0 ]; then
    echo "✓ 任务成功完成！"
    echo ""
    echo "输出文件位置:"
    echo "  - EAS LD数据: eas_ld/"
    echo "  - 摘要文件: $(cat eas_result_summary.txt 2>/dev/null || echo '未找到')"
else
    echo "✗ 任务失败，退出码: $EXIT_CODE"
    echo "请查看日志文件获取详细信息"
fi

echo ""
echo "=========================================="
