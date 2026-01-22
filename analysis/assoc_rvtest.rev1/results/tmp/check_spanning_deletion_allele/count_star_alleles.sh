#!/usr/bin/env zsh

# 脚本功能：计数 .bim 文件中第二列（Variant ID）以 ":*" 结尾的行数
# 用途：统计 Spanning Deletion Alleles

# 检查参数
if [[ -z "$1" ]]; then
    echo "Usage: $0 <path_to_bim_file>"
    exit 1
fi

BIM_FILE="$1"

# 检查文件是否存在
if [[ ! -f "$BIM_FILE" ]]; then
    echo "Error: File '$BIM_FILE' not found."
    exit 1
fi

echo "Processing bim file: $BIM_FILE"
echo "Counting IDs in column 2 ending with ':*'..."

# 使用 awk 进行处理
# $2 代表第二列
# ~ /:\*$/ 匹配以 :* 结尾的字符串
count=$(awk '$2 ~ /:\*$/ {sum++} END {print sum+0}' "$BIM_FILE")

echo "------------------------------------------------"
echo "Count of spanning deletion alleles (*): $count"
echo "------------------------------------------------"
