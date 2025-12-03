#!/bin/zsh
#SBATCH --job-name=vmiss_analysis_all
#SBATCH --output=vmiss_analysis_all_%j.out
#SBATCH --error=vmiss_analysis_all_%j.err
#SBATCH -p gr10478b
#SBATCH -t 168:0:0
#SBATCH --rsc p=1:t=64:c=32:m=146272M

# 并行运行所有4种模式的VMISS阈值分析

# 设置工作目录
cd /LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/tuning.vmiss

# 激活conda环境
source $(conda info --base)/etc/profile.d/conda.sh
conda activate cteph_geno_pro

# 定义输入文件路径
export VARIANT_METRICS_FILE="./vmiss_results/cteph_agp3k.sqc.rm_maf0_vmiss1_variant_metrics.tsv"

# 检查输入文件是否存在
if [[ ! -f ${VARIANT_METRICS_FILE} ]]; then
    echo "错误：输入文件不存在: ${VARIANT_METRICS_FILE}"
    exit 1
fi

echo "============================================================"
echo "并行运行所有4种模式的VMISS阈值分析"
echo "输入文件: ${VARIANT_METRICS_FILE}"
echo "开始时间: $(date)"
echo "============================================================"

# 定义运行单个模式的函数
run_analysis_mode() {
    local MODE=$1
    local OUTPUT_DIR="./vmiss_analysis_mode${MODE}"
    
    echo ""
    echo "----------------------------------------"
    echo "启动模式${MODE}分析..."
    echo "输出目录: ${OUTPUT_DIR}"
    echo "----------------------------------------"
    
    mkdir -p ${OUTPUT_DIR}
    
    python3 << EOF
import sys
import os
import importlib

sys.path.insert(0, '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/tuning.vmiss/scripts')

if 'tuning_vmiss_tools' in sys.modules:
    importlib.reload(sys.modules['tuning_vmiss_tools'])
else:
    import tuning_vmiss_tools

from tuning_vmiss_tools import analyze_vmiss_thresholds

variant_metrics_file = "${VARIANT_METRICS_FILE}"
output_dir = "${OUTPUT_DIR}"
mode = ${MODE}

print(f"\n{'='*60}")
print(f"模式{mode}分析")
print(f"输入: {variant_metrics_file}")
print(f"输出: {output_dir}")
print(f"{'='*60}\n")

try:
    # 根据模式设置参数
    if mode == 1:
        # 模式1：所有变体，单一VMISS
        results = analyze_vmiss_thresholds(
            variant_metrics_file=variant_metrics_file,
            output_dir=output_dir,
            analysis_mode=1,
            hist_step=0.01,
            cdf_step=0.01,
            bin_position='right',
            knee_curve='concave',
            knee_direction='increasing',
            knee_S=1.0,
            knee_weight_x=10.0,
            knee_weight_y=1.0,
            dpi=600,
            chunksize=50000,
            n_threads=8
        )
    elif mode == 2:
        # 模式2：所有变体，分15X和30X
        results = analyze_vmiss_thresholds(
            variant_metrics_file=variant_metrics_file,
            output_dir=output_dir,
            analysis_mode=2,
            hist_step=0.01,
            cdf_step_15x=0.01,
            cdf_step_30x=0.01,
            bin_position='right',
            knee_curve='concave',
            knee_direction='increasing',
            knee_S=1.0,
            knee_weight_x=10.0,
            knee_weight_y=1.0,
            dpi=600,
            chunksize=50000,
            n_threads=8
        )
    elif mode == 3:
        # 模式3：三组MAF，单一VMISS
        results = analyze_vmiss_thresholds(
            variant_metrics_file=variant_metrics_file,
            output_dir=output_dir,
            analysis_mode=3,
            maf_thresholds=(0.01, 0.05),
            hist_step=0.01,
            cdf_step=0.01,
            bin_position='right',
            knee_curve='concave',
            knee_direction='increasing',
            knee_S=1.0,
            knee_weight_x=10.0,
            knee_weight_y=1.0,
            dpi=600,
            chunksize=50000,
            n_threads=8
        )
    elif mode == 4:
        # 模式4：三组MAF，分15X和30X
        results = analyze_vmiss_thresholds(
            variant_metrics_file=variant_metrics_file,
            output_dir=output_dir,
            analysis_mode=4,
            maf_thresholds=(0.01, 0.05),
            hist_step=0.01,
            cdf_step_15x=0.01,
            cdf_step_30x=0.01,
            bin_position='right',
            knee_curve='concave',
            knee_direction='increasing',
            knee_S=1.0,
            knee_weight_x=10.0,
            knee_weight_y=1.0,
            dpi=600,
            chunksize=50000,
            n_threads=8
        )
    
    print(f"\n{'='*60}")
    print(f"✓ 模式{mode}分析完成！")
    print(f"{'='*60}")
    print(f"JSON: {results['json_path']}")
    print(f"图片: {results['plot_path']}")
    print(f"\n拐点检测结果:")
    for group_key, knee_info in results['knee_points'].items():
        if knee_info.get('knee_found'):
            print(f"  [{group_key}]")
            print(f"    阈值: {knee_info['knee_vmiss']:.4f}")
            print(f"    累积: {knee_info['knee_cumulative_percentage']:.2f}%")
            print(f"    变体数: {knee_info['knee_cumulative_count']:,}/{knee_info['total_variants']:,}")
        else:
            print(f"  [{group_key}]: {knee_info.get('message', '未找到拐点')}")
    print(f"{'='*60}\n")
    
    # 写入完成标记
    with open(os.path.join(output_dir, 'DONE'), 'w') as f:
        f.write(f"模式{mode}完成\n")
    
except Exception as e:
    print(f"\n{'='*60}")
    print(f"✗ 模式{mode}分析失败！")
    print(f"{'='*60}")
    print(f"错误: {str(e)}")
    import traceback
    traceback.print_exc()
    sys.exit(1)

EOF
    
    local EXIT_CODE=$?
    if [[ ${EXIT_CODE} -eq 0 ]]; then
        echo "✓ 模式${MODE}完成 ($(date))"
    else
        echo "✗ 模式${MODE}失败 ($(date))"
    fi
    
    return ${EXIT_CODE}
}

# 并行运行所有4种模式
echo ""
echo "开始并行运行4种模式..."
echo ""

# 导出函数以便后台进程可以使用
export -f run_analysis_mode
export VARIANT_METRICS_FILE

# 启动4个后台任务
run_analysis_mode 1 > vmiss_mode1_$$.log 2>&1 &
PID1=$!
echo "模式1启动 (PID: $PID1)"

run_analysis_mode 2 > vmiss_mode2_$$.log 2>&1 &
PID2=$!
echo "模式2启动 (PID: $PID2)"

run_analysis_mode 3 > vmiss_mode3_$$.log 2>&1 &
PID3=$!
echo "模式3启动 (PID: $PID3)"

run_analysis_mode 4 > vmiss_mode4_$$.log 2>&1 &
PID4=$!
echo "模式4启动 (PID: $PID4)"

echo ""
echo "所有模式已启动，等待完成..."
echo "监控日志文件："
echo "  - vmiss_mode1_$$.log"
echo "  - vmiss_mode2_$$.log"
echo "  - vmiss_mode3_$$.log"
echo "  - vmiss_mode4_$$.log"
echo ""

# 等待所有后台任务完成
wait $PID1
EXIT1=$?
wait $PID2
EXIT2=$?
wait $PID3
EXIT3=$?
wait $PID4
EXIT4=$?

# 汇总结果
echo ""
echo "============================================================"
echo "所有模式运行完成"
echo "结束时间: $(date)"
echo "============================================================"
echo ""
echo "执行结果："
[[ ${EXIT1} -eq 0 ]] && echo "  ✓ 模式1: 成功" || echo "  ✗ 模式1: 失败 (退出码: ${EXIT1})"
[[ ${EXIT2} -eq 0 ]] && echo "  ✓ 模式2: 成功" || echo "  ✗ 模式2: 失败 (退出码: ${EXIT2})"
[[ ${EXIT3} -eq 0 ]] && echo "  ✓ 模式3: 成功" || echo "  ✗ 模式3: 失败 (退出码: ${EXIT3})"
[[ ${EXIT4} -eq 0 ]] && echo "  ✓ 模式4: 成功" || echo "  ✗ 模式4: 失败 (退出码: ${EXIT4})"
echo ""
echo "输出目录："
echo "  - ./vmiss_analysis_mode1/"
echo "  - ./vmiss_analysis_mode2/"
echo "  - ./vmiss_analysis_mode3/"
echo "  - ./vmiss_analysis_mode4/"
echo ""

# 检查是否全部成功
if [[ ${EXIT1} -eq 0 && ${EXIT2} -eq 0 && ${EXIT3} -eq 0 && ${EXIT4} -eq 0 ]]; then
    echo "🎉 所有模式分析成功完成！"
    exit 0
else
    echo "⚠️  部分模式分析失败，请检查日志文件"
    exit 1
fi
