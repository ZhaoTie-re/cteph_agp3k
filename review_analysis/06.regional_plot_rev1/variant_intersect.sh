#!/bin/bash
#SBATCH --job-name=eas_variant_intersect
#SBATCH --output=eas_variant_intersect_%j.log
#SBATCH --error=eas_variant_intersect_%j.err
#SBATCH -p gr10478b
#SBATCH -t 168:0:0
#SBATCH --rsc p=1:t=16:c=8:m=36568M

# Script to find variant intersections between EAS bed/bim files and lead variant summary stats
# 为每个lead variant计算与EAS数据的交集，生成变体ID列表

# 默认参数
BED_PREFIX="/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/review_analysis/06.regional_plot_rev1/eas_all"
JSON_FILE="cteph_agp3k.lowfreq_common.ld_matrices_by_lead.summary.json"
OUTPUT_DIR="eas_ld_tmp"

# 解析命令行参数
while [[ $# -gt 0 ]]; do
    case $1 in
        --bed-prefix)
            BED_PREFIX="$2"
            shift 2
            ;;
        --json-file)
            JSON_FILE="$2"
            shift 2
            ;;
        --output-dir)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        *)
            echo "Unknown option: $1"
            echo "Usage: $0 [--bed-prefix PATH] [--json-file PATH] [--output-dir PATH]"
            exit 1
            ;;
    esac
done

# 检查必要文件
BIM_FILE="${BED_PREFIX}.bim"

if [ ! -f "$BIM_FILE" ]; then
    echo "Error: BIM file not found: $BIM_FILE"
    exit 1
fi

if [ ! -f "$JSON_FILE" ]; then
    echo "Error: JSON file not found: $JSON_FILE"
    exit 1
fi

# 创建输出目录
mkdir -p "$OUTPUT_DIR"

echo "=========================================="
echo "EAS Variant Intersection Pipeline"
echo "=========================================="
echo "BED prefix: $BED_PREFIX"
echo "BIM file: $BIM_FILE"
echo "JSON file: $JSON_FILE"
echo "Output directory: $OUTPUT_DIR"
echo "=========================================="

# 步骤1: 提取EAS变体ID（使用awk节省内存，只读取第2列）
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Step 1: Extracting EAS variant IDs from BIM file..."
EAS_VARIANTS_FILE="${OUTPUT_DIR}/eas_all_variants.txt"

awk '{print $2}' "$BIM_FILE" > "$EAS_VARIANTS_FILE"
EAS_VARIANT_COUNT=$(wc -l < "$EAS_VARIANTS_FILE")
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Found $EAS_VARIANT_COUNT EAS variants"

# 步骤2: 排序EAS变体ID（用于快速查找）
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Step 2: Sorting EAS variant IDs..."
sort "$EAS_VARIANTS_FILE" -o "${EAS_VARIANTS_FILE}.sorted"
mv "${EAS_VARIANTS_FILE}.sorted" "$EAS_VARIANTS_FILE"

# 步骤3: 解析JSON并处理每个lead variant
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Step 3: Processing lead variants..."

# 提取lead variants列表
LEAD_VARIANTS=$(python3 -c "
import json
import sys

with open('$JSON_FILE', 'r') as f:
    data = json.load(f)

if 'per_lead_outputs' in data:
    for lead_var in data['per_lead_outputs'].keys():
        print(lead_var)
")

# 初始化结果JSON
RESULT_JSON="${OUTPUT_DIR}/variant_intersection_summary.json"
echo "{" > "$RESULT_JSON"
echo "  \"created_at\": \"$(date '+%Y-%m-%d %H:%M:%S')\"," >> "$RESULT_JSON"
echo "  \"eas_bed_prefix\": \"$BED_PREFIX\"," >> "$RESULT_JSON"
echo "  \"source_json\": \"$JSON_FILE\"," >> "$RESULT_JSON"
echo "  \"n_eas_variants_total\": $EAS_VARIANT_COUNT," >> "$RESULT_JSON"
echo "  \"per_lead_intersections\": {" >> "$RESULT_JSON"

FIRST_LEAD=1
TOTAL_COMMON=0
LEAD_COUNT=0

# 处理每个lead variant
for LEAD_VAR in $LEAD_VARIANTS; do
    LEAD_COUNT=$((LEAD_COUNT + 1))
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Processing lead variant $LEAD_COUNT: $LEAD_VAR"
    
    # 获取该lead variant的sum_stat_tsv文件路径
    SUM_STAT_FILE=$(python3 -c "
import json
with open('$JSON_FILE', 'r') as f:
    data = json.load(f)
print(data['per_lead_outputs']['$LEAD_VAR']['sum_stat_tsv'])
")
    
    if [ ! -f "$SUM_STAT_FILE" ]; then
        echo "  Warning: Sum stat file not found: $SUM_STAT_FILE"
        continue
    fi
    
    # 提取sum_stat中的SNPID（跳过表头）
    SAFE_VAR_NAME=$(echo "$LEAD_VAR" | sed 's/:/_/g' | sed 's/>/_/g' | sed 's/</_/g')
    SUMSTAT_VARIANTS="${OUTPUT_DIR}/${SAFE_VAR_NAME}.sumstat_variants.txt"
    
    tail -n +2 "$SUM_STAT_FILE" | cut -f1 > "$SUMSTAT_VARIANTS"
    SUMSTAT_COUNT=$(wc -l < "$SUMSTAT_VARIANTS")
    
    # 排序sum_stat变体（用于comm命令）
    sort "$SUMSTAT_VARIANTS" -o "${SUMSTAT_VARIANTS}.sorted"
    
    # 使用comm命令找交集（两个文件都已排序）
    INTERSECT_FILE="${OUTPUT_DIR}/${SAFE_VAR_NAME}.intersect_variants.txt"
    comm -12 "$EAS_VARIANTS_FILE" "${SUMSTAT_VARIANTS}.sorted" > "$INTERSECT_FILE"
    COMMON_COUNT=$(wc -l < "$INTERSECT_FILE")
    
    # 检查lead variant是否在交集中
    LEAD_IN_INTERSECT="false"
    if grep -qxF "$LEAD_VAR" "$INTERSECT_FILE"; then
        LEAD_IN_INTERSECT="true"
    fi
    
    TOTAL_COMMON=$((TOTAL_COMMON + COMMON_COUNT))
    
    echo "  - Sum stat variants: $SUMSTAT_COUNT"
    echo "  - Common variants: $COMMON_COUNT"
    echo "  - Lead variant in intersection: $LEAD_IN_INTERSECT"
    
    # 清理临时文件
    rm -f "$SUMSTAT_VARIANTS" "${SUMSTAT_VARIANTS}.sorted"
    
    # 添加到JSON（处理逗号）
    if [ $FIRST_LEAD -eq 0 ]; then
        echo "    ," >> "$RESULT_JSON"
    fi
    FIRST_LEAD=0
    
    echo "    \"$LEAD_VAR\": {" >> "$RESULT_JSON"
    echo "      \"lead_variant\": \"$LEAD_VAR\"," >> "$RESULT_JSON"
    echo "      \"sum_stat_file\": \"$SUM_STAT_FILE\"," >> "$RESULT_JSON"
    echo "      \"intersect_file\": \"$INTERSECT_FILE\"," >> "$RESULT_JSON"
    echo "      \"n_sumstat_variants\": $SUMSTAT_COUNT," >> "$RESULT_JSON"
    echo "      \"n_common_variants\": $COMMON_COUNT," >> "$RESULT_JSON"
    echo "      \"lead_variant_in_intersection\": $LEAD_IN_INTERSECT" >> "$RESULT_JSON"
    echo -n "    }" >> "$RESULT_JSON"
done

# 完成JSON
echo "" >> "$RESULT_JSON"
echo "  }," >> "$RESULT_JSON"
echo "  \"summary\": {" >> "$RESULT_JSON"
echo "    \"n_lead_variants\": $LEAD_COUNT," >> "$RESULT_JSON"
echo "    \"n_common_variants_total\": $TOTAL_COMMON" >> "$RESULT_JSON"
echo "  }" >> "$RESULT_JSON"
echo "}" >> "$RESULT_JSON"

echo "=========================================="
echo "Processing Complete!"
echo "=========================================="
echo "Total lead variants processed: $LEAD_COUNT"
echo "Total common variants: $TOTAL_COMMON"
echo "Output directory: $OUTPUT_DIR"
echo "Summary JSON: $RESULT_JSON"
echo "=========================================="

# 清理EAS变体列表（可选，如果想保留则注释掉）
# rm -f "$EAS_VARIANTS_FILE"

exit 0
