#!/bin/bash
# 02_check_duplicate_iid.sh

# 列出所有fam文件
echo "Checking fam files:"
ls -1 *.fam
echo ""

# 收集所有fam文件的IID（第二列）
all_iids=$(cat *.fam | awk '{print $2}')

# 统计总数和唯一数
total_count=$(echo "$all_iids" | wc -l)
unique_count=$(echo "$all_iids" | sort -u | wc -l)

echo "Total IIDs across all fam files: $total_count"
echo "Unique IIDs: $unique_count"

# 检查是否有重复
if [ $total_count -ne $unique_count ]; then
    echo ""
    echo "ERROR: Found duplicate IIDs across different fam files!"
    echo "Duplicate IIDs:"
    echo "$all_iids" | sort | uniq -d
    echo ""
    echo "Showing which files contain duplicates:"
    for iid in $(echo "$all_iids" | sort | uniq -d); do
        echo "  IID: $iid found in:"
        grep -l "\s$iid\s" *.fam | sed 's/^/    /'
    done
    exit 1
else
    echo ""
    echo "SUCCESS: No duplicate IIDs found across fam files."
fi
