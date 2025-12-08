params.ArrayRAWPath = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/array/00.raw_data_ph'
params.SampleSelectList = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/array_rev1/info/cteph_agp3k.v4.ls'
params.NagasakiPipelinePath = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/nagasaki_pipeline'
params.ChrRenameFile = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/Info/chr_rename.txt'
params.OutputDir = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/array_rev1/results'
params.Plink2Path = '/home/b/b37974/plink2'
params.PlinkPath = '/home/b/b37974/plink'
params.SampleInfo = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/info/cteph_agp3k_jhrpv4.rev1.xlsx'
params.SampleIDColumn = 'ID'
params.SampleSexColumn = 'Sex'
params.OutcomeColumn = 'OUTCOME2'
params.CaseValue = 'CTEPH'
// MAF filtering parameters
params.MAF_Threshold = 0.01          // MAF threshold for filtering
params.MAF_FilterMode = 'ALL'        // Options: 'ALL', 'CASE', 'CTRL'
// HWE filtering parameters
params.HWE_FilterMode = 'CTRL|CASE'  // Options: 'CTRL', 'CASE', 'CTRL|CASE'
params.HWE_CTRL_Threshold = 1e-6     // HWE threshold for controls
params.HWE_CASE_Threshold = 1e-10    // HWE threshold for cases (only used if FilterMode includes 'CASE')

// 识别成配对的 .bed, .bim, .fam 文件
Channel
    .fromFilePairs("${params.ArrayRAWPath}/*.{bed,bim,fam}", size: 3) { file ->
        // 提取基础文件名（去掉扩展名）
        file.name.replaceAll(/\.(bed|bim|fam)$/, '')
    }
    .map { prefix, files ->
        // 只保留文件，不保留prefix: [bed, bim, fam]
        def bed = files.find { it.name.endsWith('.bed') }
        def bim = files.find { it.name.endsWith('.bim') }
        def fam = files.find { it.name.endsWith('.fam') }
        tuple(bed, bim, fam)
    }
    .set { plink_files_ch }


// 样本选择
process SampleSelect {
    executor 'slurm'
    queue 'gr10478b'
    time '12h'
    tag "${bed.simpleName}"
    publishDir "${params.OutputDir}/01.sample_select", mode: 'symlink'

    input:
    tuple path(bed), path(bim), path(fam) from plink_files_ch
    path(sample_list) from params.SampleSelectList

    output:
    tuple path("*.selected.bed"), path("*.selected.bim"), path("*.selected.fam") into selected_files
    path "*.selected.log"
    path "sample_list.formatted.txt"

    script:
    def prefix = bed.baseName
    """
    # 格式化样本列表：添加列名并复制为两列（FID和IID）
    echo -e "#FID\\tIID" > sample_list.formatted.txt
    awk '{print \$1"\\t"\$1}' ${sample_list} >> sample_list.formatted.txt
    
    # 选择样本
    ${params.Plink2Path} \
        --bfile ${prefix} \
        --keep sample_list.formatted.txt \
        --make-bed \
        --out ${prefix}.selected \
        --threads 4
    """
}

// 样本QC: 去除call-rate低于99%的样本
process SampleQC {
    executor 'slurm'
    queue 'gr10478b'
    time '36h'
    tag "${bed.simpleName}"
    publishDir "${params.OutputDir}/02.sample_qc", mode: 'symlink'
    
    input:
    tuple path(bed), path(bim), path(fam) from selected_files

    output:
    tuple path("*.sqc.bed"), path("*.sqc.bim"), path("*.sqc.fam") into sqc_out, sqc_out_2
    path "*.smiss"
    path "*.sqc.log"

    script:
    def prefix = bed.baseName
    """
    # 计算样本缺失率
    ${params.Plink2Path} \
        --bfile ${prefix} \
        --missing sample-only \
        --out ${prefix} \
        --threads 4

    # 过滤call-rate < 99%的样本 (即缺失率 > 0.01)
    ${params.Plink2Path} \
        --bfile ${prefix} \
        --mind 0.01 \
        --make-bed \
        --out ${prefix}.sqc \
        --threads 4
    """
}

// 将所有fam文件收集到一起，准备检查重复
sqc_out
    .map { bed, bim, fam -> fam }
    .collect()
    .set { all_fam_files }

// 检查多个fam文件之间的IID重复
process CheckDuplicateIID {
    executor 'local'
    tag "Checking IID duplicates across all fam files"
    
    input:
    path(fam_files) from all_fam_files

    output:
    val(true) into dup_check_done

    script:
    """
    #!/bin/bash
    
    # 列出所有fam文件
    echo "Checking fam files:"
    ls -1 *.fam
    echo ""
    
    # 收集所有fam文件的IID（第二列）
    all_iids=\$(cat *.fam | awk '{print \$2}')
    
    # 统计总数和唯一数
    total_count=\$(echo "\$all_iids" | wc -l)
    unique_count=\$(echo "\$all_iids" | sort -u | wc -l)
    
    echo "Total IIDs across all fam files: \$total_count"
    echo "Unique IIDs: \$unique_count"
    
    # 检查是否有重复
    if [ \$total_count -ne \$unique_count ]; then
        echo ""
        echo "ERROR: Found duplicate IIDs across different fam files!"
        echo "Duplicate IIDs:"
        echo "\$all_iids" | sort | uniq -d
        echo ""
        echo "Showing which files contain duplicates:"
        for iid in \$(echo "\$all_iids" | sort | uniq -d); do
            echo "  IID: \$iid found in:"
            grep -l "\\s\$iid\\s" *.fam | sed 's/^/    /'
        done
        exit 1
    else
        echo ""
        echo "SUCCESS: No duplicate IIDs found across fam files."
    fi
    """
}

// 如果检查通过，继续处理
sqc_out_2
    .combine(dup_check_done)
    .map { bed, bim, fam, check -> tuple(bed, bim, fam) }
    .set { checked_files }

// 转换为VCF格式并标准化
process ConvertToVCF {
    executor 'slurm'
    queue 'gr10478b'
    time '24h'
    tag "${bed.simpleName}"
    publishDir "${params.OutputDir}/03.vcf_convert", mode: 'symlink'
    
    input:
    tuple path(bed), path(bim), path(fam) from checked_files
    path(chr_rename) from params.ChrRenameFile

    output:
    tuple path("*.norm.vcf.gz"), path("*.norm.vcf.gz.tbi") into vcf_out
    path "*.norm.vcf.gz.log"

    script:
    def prefix = bed.baseName
    """
    # 优化的流式处理：过滤染色体 -> PLINK转VCF -> 重命名染色体 -> 标准化 -> 填充ID
    # 使用管道减少I/O操作
    
    # 步骤0: 过滤只保留标准染色体（1-22, X, Y, MT）并清理非标准等位基因
    echo "=== VCF Conversion QC Report ===" > ${prefix}.norm.vcf.gz.log
    echo "Processing: ${prefix}" >> ${prefix}.norm.vcf.gz.log
    echo "Date: \$(date)" >> ${prefix}.norm.vcf.gz.log
    echo "" >> ${prefix}.norm.vcf.gz.log
    
    echo "Step 0: Input variants (before filtering)" >> ${prefix}.norm.vcf.gz.log
    wc -l ${prefix}.bim | awk '{printf "  Total variants: %'"'"'d\\n", \$1}' >> ${prefix}.norm.vcf.gz.log
    cut -f1 ${prefix}.bim | sort | uniq -c | awk '{printf "  Chr %s: %'"'"'d variants\\n", \$2, \$1}' >> ${prefix}.norm.vcf.gz.log
    echo "" >> ${prefix}.norm.vcf.gz.log
    
    # 统计问题变异
    echo "  Quality issues detected:" >> ${prefix}.norm.vcf.gz.log
    # 使用更简单的匹配方法
    missing_allele_count=\$(awk '\$5 == "." || \$6 == "."' ${prefix}.bim | wc -l)
    non_standard_count=\$(awk '\$5 !~ /^[ATCG]\$/ || \$6 !~ /^[ATCG]\$/' ${prefix}.bim | wc -l)
    non_acgt_count=\$((non_standard_count))
    true_non_acgt_count=\$((non_acgt_count - missing_allele_count))
    
    printf "    - Non-ACGT alleles (indels, etc.): %'d\\n" \${true_non_acgt_count} >> ${prefix}.norm.vcf.gz.log
    printf "    - Missing allele information (. in REF/ALT): %'d\\n" \${missing_allele_count} >> ${prefix}.norm.vcf.gz.log
    printf "    - Total problematic variants: %'d\\n" \${non_acgt_count} >> ${prefix}.norm.vcf.gz.log
    echo "" >> ${prefix}.norm.vcf.gz.log
    
    ${params.Plink2Path} \
        --bfile ${prefix} \
        --chr 1-22,X,Y,MT \
        --snps-only just-acgt \
        --make-bed \
        --out temp_filtered_${prefix} \
        --threads 8
    
    echo "Step 1: After chromosome filtering (1-22,X,Y,MT) and ACGT-only SNPs" >> ${prefix}.norm.vcf.gz.log
    wc -l temp_filtered_${prefix}.bim | awk '{printf "  Remaining variants: %'"'"'d\\n", \$1}' >> ${prefix}.norm.vcf.gz.log
    
    # 计算被过滤掉的变体数量
    original_count=\$(wc -l < ${prefix}.bim)
    filtered_count=\$(wc -l < temp_filtered_${prefix}.bim)
    removed_count=\$((original_count - filtered_count))
    printf "  Removed variants: %'d\\n" \${removed_count} >> ${prefix}.norm.vcf.gz.log
    echo "  Removal rate: \$(awk "BEGIN {printf \\"%.2f%%\\", (\${removed_count}/\${original_count})*100}")" >> ${prefix}.norm.vcf.gz.log
    echo "" >> ${prefix}.norm.vcf.gz.log
    
    # 检查过滤后是否还有缺失等位基因
    remaining_missing=\$(awk '\$5 == "." || \$6 == "."' temp_filtered_${prefix}.bim | wc -l)
    if [ \${remaining_missing} -gt 0 ]; then
        echo "  WARNING: \${remaining_missing} variants with missing alleles remain after PLINK2 filtering" >> ${prefix}.norm.vcf.gz.log
        echo "  These will be removed manually..." >> ${prefix}.norm.vcf.gz.log
        
        # 手动移除缺失等位基因的变异
        awk '\$5 != "." && \$6 != "."' temp_filtered_${prefix}.bim > temp_filtered_${prefix}.bim.clean
        
        # 提取这些变异的ID用于过滤
        awk '{print \$2}' temp_filtered_${prefix}.bim.clean > variants_to_keep.txt
        
        ${params.Plink2Path} \
            --bfile temp_filtered_${prefix} \
            --extract variants_to_keep.txt \
            --make-bed \
            --out temp_filtered_clean_${prefix} \
            --threads 8
        
        # 替换文件
        mv temp_filtered_clean_${prefix}.bed temp_filtered_${prefix}.bed
        mv temp_filtered_clean_${prefix}.bim temp_filtered_${prefix}.bim
        mv temp_filtered_clean_${prefix}.fam temp_filtered_${prefix}.fam
        
        # 更新统计
        filtered_count=\$(wc -l < temp_filtered_${prefix}.bim)
        removed_count=\$((original_count - filtered_count))
        
        printf "  After removing missing alleles: %'d variants\\n" \${filtered_count} >> ${prefix}.norm.vcf.gz.log
        printf "  Total removed: %'d\\n" \${removed_count} >> ${prefix}.norm.vcf.gz.log
        
        rm -f variants_to_keep.txt temp_filtered_${prefix}.bim.clean
    fi
    
    # 步骤2: PLINK转VCF（写入磁盘，因为plink不支持stdout）
    echo "" >> ${prefix}.norm.vcf.gz.log
    echo "Step 2: PLINK binary to VCF conversion" >> ${prefix}.norm.vcf.gz.log
    
    ${params.PlinkPath} \
        --bfile temp_filtered_${prefix} \
        --recode vcf-iid bgz \
        --out temp_${prefix}
    
    vcf_count=\$(bcftools view -H temp_${prefix}.vcf.gz | wc -l)
    printf "  VCF variants: %'d\\n" \${vcf_count} >> ${prefix}.norm.vcf.gz.log
    echo "  Format: VCF 4.2, bgzipped" >> ${prefix}.norm.vcf.gz.log
    echo "  Reason: Convert PLINK binary format to VCF for downstream processing" >> ${prefix}.norm.vcf.gz.log
    
    # 步骤3: 重命名染色体
    echo "" >> ${prefix}.norm.vcf.gz.log
    echo "Step 3: Chromosome renaming (numeric to chr-prefix)" >> ${prefix}.norm.vcf.gz.log
    
    bcftools annotate \
        --rename-chrs ${chr_rename} \
        --threads 8 \
        -Oz \
        -o temp_renamed_${prefix}.vcf.gz \
        temp_${prefix}.vcf.gz || { echo "Error: bcftools annotate (rename) failed" >> ${prefix}.norm.vcf.gz.log; exit 1; }
    
    bcftools index -t temp_renamed_${prefix}.vcf.gz
    
    renamed_count=\$(bcftools view -H temp_renamed_${prefix}.vcf.gz | wc -l)
    printf "  Variants after renaming: %'d\\n" \${renamed_count} >> ${prefix}.norm.vcf.gz.log
    echo "  Chromosomes: 1-22 → chr1-chr22, 23 → chrX, 24 → chrY, 26 → chrM" >> ${prefix}.norm.vcf.gz.log
    echo "  Reason: Match reference genome chromosome naming convention (GRCh38)" >> ${prefix}.norm.vcf.gz.log
    
    # 步骤4: 标准化并严格检查REF一致性
    echo "" >> ${prefix}.norm.vcf.gz.log
    echo "Step 4: VCF normalization with reference genome" >> ${prefix}.norm.vcf.gz.log
    echo "  Using --check-ref s (skip variants with REF mismatch)" >> ${prefix}.norm.vcf.gz.log
    
    # 首先记录重命名后的变体数
    before_norm_count=\${renamed_count}
    
    # 执行标准化，使用-s参数跳过REF不匹配的变体
    bcftools norm \
        --multiallelics -any \
        --fasta-ref ${params.NagasakiPipelinePath}/data/hs38DH.fa \
        --check-ref s \
        --threads 8 \
        -Oz \
        -o temp_norm_${prefix}.vcf.gz \
        temp_renamed_${prefix}.vcf.gz \
        2> ${prefix}.norm_warnings.txt || { echo "Error: bcftools norm failed" >> ${prefix}.norm.vcf.gz.log; exit 1; }
    
    bcftools index -t temp_norm_${prefix}.vcf.gz
    
    # 统计标准化后的变体数
    norm_count=\$(bcftools view -H temp_norm_${prefix}.vcf.gz | wc -l)
    total_norm_change=\$((norm_count - before_norm_count))
    
    printf "  Variants before normalization: %'d\\n" \${before_norm_count} >> ${prefix}.norm.vcf.gz.log
    printf "  Variants after normalization: %'d\\n" \${norm_count} >> ${prefix}.norm.vcf.gz.log
    
    if [ \${total_norm_change} -eq 0 ]; then
        echo "  Net change: 0 (no variants added or removed)" >> ${prefix}.norm.vcf.gz.log
    elif [ \${total_norm_change} -gt 0 ]; then
        printf "  Net change: +%'d variants (from multiallelic splitting)\\n" \${total_norm_change} >> ${prefix}.norm.vcf.gz.log
    else
        removed=\$((-total_norm_change))
        printf "  Net change: -%'d variants removed\\n" \${removed} >> ${prefix}.norm.vcf.gz.log
    fi
    echo "" >> ${prefix}.norm.vcf.gz.log
    
    # 解析bcftools norm的输出统计
    echo "  Normalization summary (from bcftools norm):" >> ${prefix}.norm.vcf.gz.log
    cat ${prefix}.norm_warnings.txt >> ${prefix}.norm.vcf.gz.log
    echo "" >> ${prefix}.norm.vcf.gz.log
    
    # 提取关键统计信息
    if [ -f ${prefix}.norm_warnings.txt ]; then
        # 提取Lines统计行
        lines_stats=\$(grep "^Lines" ${prefix}.norm_warnings.txt 2>/dev/null || echo "")
        if [ ! -z "\$lines_stats" ]; then
            # 解析各个数值: total/split/joined/realigned/removed/skipped
            # 提取冒号后的数字部分
            numbers_part=\$(echo "\$lines_stats" | sed 's/.*: *//')
            total_lines=\$(echo "\$numbers_part" | awk -F'/' '{print \$1}')
            split_lines=\$(echo "\$numbers_part" | awk -F'/' '{print \$2}')
            joined_lines=\$(echo "\$numbers_part" | awk -F'/' '{print \$3}')
            realigned_lines=\$(echo "\$numbers_part" | awk -F'/' '{print \$4}')
            removed_lines=\$(echo "\$numbers_part" | awk -F'/' '{print \$5}')
            skipped_lines=\$(echo "\$numbers_part" | awk -F'/' '{print \$6}')
            
            echo "  Detailed breakdown:" >> ${prefix}.norm.vcf.gz.log
            printf "    - Total variants processed: %'d\\n" \${total_lines} >> ${prefix}.norm.vcf.gz.log
            printf "    - Multiallelic sites split: %'d\\n" \${split_lines} >> ${prefix}.norm.vcf.gz.log
            printf "    - Variants joined: %'d\\n" \${joined_lines} >> ${prefix}.norm.vcf.gz.log
            printf "    - Indels realigned: %'d\\n" \${realigned_lines} >> ${prefix}.norm.vcf.gz.log
            printf "    - Variants removed: %'d\\n" \${removed_lines} >> ${prefix}.norm.vcf.gz.log
            printf "    - Variants skipped: %'d\\n" \${skipped_lines} >> ${prefix}.norm.vcf.gz.log
            echo "" >> ${prefix}.norm.vcf.gz.log
        fi
        
        # 提取REF/ALT统计行
        ref_alt_stats=\$(grep "^REF/ALT" ${prefix}.norm_warnings.txt 2>/dev/null || echo "")
        if [ ! -z "\$ref_alt_stats" ]; then
            # 提取冒号后的数字部分
            ref_numbers_part=\$(echo "\$ref_alt_stats" | sed 's/.*: *//')
            ref_alt_total=\$(echo "\$ref_numbers_part" | awk -F'/' '{print \$1}')
            ref_alt_modified=\$(echo "\$ref_numbers_part" | awk -F'/' '{print \$2}')
            ref_alt_added=\$(echo "\$ref_numbers_part" | awk -F'/' '{print \$3}')
            
            echo "  REF/ALT allele adjustments:" >> ${prefix}.norm.vcf.gz.log
            printf "    - Total variants: %'d\\n" \${ref_alt_total} >> ${prefix}.norm.vcf.gz.log
            printf "    - REF/ALT swapped (to match reference): %'d\\n" \${ref_alt_modified} >> ${prefix}.norm.vcf.gz.log
            printf "    - ALT alleles added: %'d\\n" \${ref_alt_added} >> ${prefix}.norm.vcf.gz.log
            
            if [ "\${ref_alt_modified}" != "0" ] && [ "\${ref_alt_total}" != "0" ] && [ ! -z "\${ref_alt_total}" ]; then
                swap_rate=\$(awk "BEGIN {printf \\"%.2f%%\\", (\${ref_alt_modified}/\${ref_alt_total})*100}")
                echo "    - Swap rate: \${swap_rate}" >> ${prefix}.norm.vcf.gz.log
                echo "    - Reason: VCF REF allele differs from reference genome" >> ${prefix}.norm.vcf.gz.log
                echo "    - Action: REF and ALT were swapped, genotypes flipped (0↔1)" >> ${prefix}.norm.vcf.gz.log
                echo "    - Note: Biological meaning preserved (e.g., AA stays AA)" >> ${prefix}.norm.vcf.gz.log
            fi
            echo "" >> ${prefix}.norm.vcf.gz.log
        fi
    fi
    
    echo "  Operations performed:" >> ${prefix}.norm.vcf.gz.log
    echo "    1. Split multiallelic sites into biallelic records" >> ${prefix}.norm.vcf.gz.log
    echo "    2. Left-align and normalize indels" >> ${prefix}.norm.vcf.gz.log
    echo "    3. Swap REF/ALT when VCF REF != reference genome (--check-ref s)" >> ${prefix}.norm.vcf.gz.log
    echo "    4. Flip genotypes accordingly to preserve biological meaning" >> ${prefix}.norm.vcf.gz.log
    
    # 步骤5: 设置ID
    echo "" >> ${prefix}.norm.vcf.gz.log
    echo "Step 5: Set variant IDs" >> ${prefix}.norm.vcf.gz.log
    
    bcftools annotate \
        --set-id '%CHROM:%POS:%REF:%ALT' \
        --threads 8 \
        -Oz \
        -o ${prefix}.norm.vcf.gz \
        temp_norm_${prefix}.vcf.gz || { echo "Error: bcftools annotate (set-id) failed" >> ${prefix}.norm.vcf.gz.log; exit 1; }
    
    # 创建索引
    bcftools index --threads 8 -t ${prefix}.norm.vcf.gz
    
    final_count=\$(bcftools view -H ${prefix}.norm.vcf.gz | wc -l)
    printf "  Final variants: %'d\\n" \${final_count} >> ${prefix}.norm.vcf.gz.log
    echo "  ID format: CHROM:POS:REF:ALT (e.g., chr1:12345:A:G)" >> ${prefix}.norm.vcf.gz.log
    echo "  Reason: Unique, reproducible variant identifiers" >> ${prefix}.norm.vcf.gz.log
    
    # 最终统计汇总
    echo "" >> ${prefix}.norm.vcf.gz.log
    echo "======================================" >> ${prefix}.norm.vcf.gz.log
    echo "SUMMARY" >> ${prefix}.norm.vcf.gz.log
    echo "======================================" >> ${prefix}.norm.vcf.gz.log
    printf "Input variants (Step 0):     %'d\\n" \${original_count} >> ${prefix}.norm.vcf.gz.log
    printf "After filtering (Step 1):    %'d (removed: %'d)\\n" \${filtered_count} \${removed_count} >> ${prefix}.norm.vcf.gz.log
    printf "After VCF conversion (Step 2): %'d\\n" \${vcf_count} >> ${prefix}.norm.vcf.gz.log
    printf "After chr rename (Step 3):   %'d\\n" \${renamed_count} >> ${prefix}.norm.vcf.gz.log
    if [ \${total_norm_change} -ge 0 ]; then
        printf "After normalization (Step 4): %'d (change: +%'d)\\n" \${norm_count} \${total_norm_change} >> ${prefix}.norm.vcf.gz.log
    else
        abs_change=\$((-total_norm_change))
        printf "After normalization (Step 4): %'d (change: -%'d)\\n" \${norm_count} \${abs_change} >> ${prefix}.norm.vcf.gz.log
    fi
    printf "Final output (Step 5):       %'d\\n" \${final_count} >> ${prefix}.norm.vcf.gz.log
    echo "" >> ${prefix}.norm.vcf.gz.log
    total_removed=\$((original_count - final_count))
    retention_rate=\$(awk "BEGIN {printf \\"%.2f%%\\", (\${final_count}/\${original_count})*100}")
    printf "Total variants removed: %'d\\n" \${total_removed} >> ${prefix}.norm.vcf.gz.log
    echo "Retention rate: \${retention_rate}" >> ${prefix}.norm.vcf.gz.log
    echo "" >> ${prefix}.norm.vcf.gz.log
    echo "Main processing steps:" >> ${prefix}.norm.vcf.gz.log
    echo "  1. Chromosome filtering: Removed non-standard chromosomes" >> ${prefix}.norm.vcf.gz.log
    echo "  2. Allele filtering: Removed non-ACGT alleles" >> ${prefix}.norm.vcf.gz.log
    echo "  3. Normalization: REF/ALT swapped where needed, multiallelic splitting" >> ${prefix}.norm.vcf.gz.log
    echo "" >> ${prefix}.norm.vcf.gz.log
    echo "Output files:" >> ${prefix}.norm.vcf.gz.log
    echo "  - ${prefix}.norm.vcf.gz (normalized VCF)" >> ${prefix}.norm.vcf.gz.log
    echo "  - ${prefix}.norm.vcf.gz.tbi (tabix index)" >> ${prefix}.norm.vcf.gz.log
    echo "  - ${prefix}.norm.vcf.gz.log (this log file)" >> ${prefix}.norm.vcf.gz.log
    echo "  - ${prefix}.norm_warnings.txt (normalization warnings)" >> ${prefix}.norm.vcf.gz.log
    
    # 清理临时文件
    rm -f temp_*${prefix}.*
    """
}

// 将VCF转换回PLINK格式
process VCFtoPLINK {
    executor 'slurm'
    queue 'gr10478b'
    time '12h'
    tag "${vcf.simpleName}"
    publishDir "${params.OutputDir}/04.plink_final", mode: 'symlink'
    
    input:
    tuple path(vcf), path(tbi) from vcf_out
    path(sample_info) from params.SampleInfo

    output:
    tuple path("*.bed"), path("*.bim"), path("*.fam") into plink_final
    path "*.log"
    path "*.update_sex.txt"

    script:
    def prefix = vcf.baseName.replaceAll(/\.vcf$/, '')
    def id_col = params.SampleIDColumn
    def sex_col = params.SampleSexColumn
    """
    # 执行Python脚本读取Excel并创建性别更新文件
    python3 << 'EOF'
import pandas as pd
import sys

# 读取Excel文件中的样本信息
print("Reading sample information from Excel file...")
try:
    df = pd.read_excel("${sample_info}")
    print(f"Total records in Excel: {len(df)}")
    
    # 检查必需的列是否存在
    if "${id_col}" not in df.columns or "${sex_col}" not in df.columns:
        print(f"Error: Required columns '${id_col}' and/or '${sex_col}' not found in Excel file")
        print(f"Available columns: {', '.join(df.columns)}")
        sys.exit(1)
    
    # 提取ID和Sex列
    sex_data = df[["${id_col}", "${sex_col}"]].copy()
    
    # 转换性别编码: F->2 (female), M->1 (male)
    sex_mapping = {'F': '2', 'M': '1', 'f': '2', 'm': '1'}
    sex_data['Sex_Code'] = sex_data["${sex_col}"].map(sex_mapping)
    
    # 检查是否有未识别的性别代码
    unknown_sex = sex_data[sex_data['Sex_Code'].isna()]
    if len(unknown_sex) > 0:
        print(f"Warning: {len(unknown_sex)} samples with unknown sex codes:")
        print(unknown_sex["${sex_col}"].value_counts())
        print("These will be set as unknown (0)")
        sex_data['Sex_Code'] = sex_data['Sex_Code'].fillna('0')
    
    # 创建PLINK格式的性别更新文件: FID IID Sex
    # FID和IID都使用样本ID
    with open("${prefix}.update_sex.txt", "w") as f:
        for _, row in sex_data.iterrows():
            sample_id = str(row["${id_col}"])
            sex_code = row['Sex_Code']
            f.write(f"{sample_id}\\t{sample_id}\\t{sex_code}\\n")
    
    print(f"Sex information written for {len(sex_data)} samples")
    print(f"  Females (F): {(sex_data['Sex_Code'] == '2').sum()}")
    print(f"  Males (M): {(sex_data['Sex_Code'] == '1').sum()}")
    print(f"  Unknown: {(sex_data['Sex_Code'] == '0').sum()}")
    
except Exception as e:
    print(f"Error reading Excel file: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)
EOF
    
    # 转换VCF到PLINK格式
    # 使用--double-id确保FID和IID相同
    # --split-par b38: 处理X染色体伪常染色体区域（PAR）
    # --update-sex: 使用从Excel读取的性别信息
    ${params.Plink2Path} \
        --vcf ${vcf} \
        --double-id \
        --split-par b38 \
        --update-sex ${prefix}.update_sex.txt \
        --make-bed \
        --out ${prefix} \
        --threads 8
    
    # 验证转换结果
    echo "VCF to PLINK conversion completed" > ${prefix}.conversion.log
    echo "Date: \$(date)" >> ${prefix}.conversion.log
    echo "" >> ${prefix}.conversion.log
    
    # 统计变异数
    variant_count=\$(wc -l < ${prefix}.bim)
    sample_count=\$(wc -l < ${prefix}.fam)
    
    printf "Variants: %'d\\n" \${variant_count} >> ${prefix}.conversion.log
    printf "Samples: %'d\\n" \${sample_count} >> ${prefix}.conversion.log
    echo "" >> ${prefix}.conversion.log
    
    # 统计性别信息
    echo "Sex distribution:" >> ${prefix}.conversion.log
    male_count=\$(awk '\$5==1' ${prefix}.fam | wc -l)
    female_count=\$(awk '\$5==2' ${prefix}.fam | wc -l)
    unknown_count=\$(awk '\$5==0' ${prefix}.fam | wc -l)
    printf "  Males: %'d\\n" \${male_count} >> ${prefix}.conversion.log
    printf "  Females: %'d\\n" \${female_count} >> ${prefix}.conversion.log
    printf "  Unknown: %'d\\n" \${unknown_count} >> ${prefix}.conversion.log
    echo "" >> ${prefix}.conversion.log
    
    # 显示前几行样本ID以确认FID=IID和性别
    echo "Sample ID format (first 5 samples):" >> ${prefix}.conversion.log
    head -5 ${prefix}.fam | awk '{sex=\$5; if(sex==1) sex_str="Male"; else if(sex==2) sex_str="Female"; else sex_str="Unknown"; print "  FID: " \$1 "  IID: " \$2 "  Sex: " sex_str " (" \$5 ")"}' >> ${prefix}.conversion.log
    echo "" >> ${prefix}.conversion.log
    
    echo "Output files:" >> ${prefix}.conversion.log
    echo "  - ${prefix}.bed (binary genotype file)" >> ${prefix}.conversion.log
    echo "  - ${prefix}.bim (variant information)" >> ${prefix}.conversion.log
    echo "  - ${prefix}.fam (sample information, FID=IID, with sex)" >> ${prefix}.conversion.log
    echo "  - ${prefix}.update_sex.txt (sex information from Excel)" >> ${prefix}.conversion.log
    echo "" >> ${prefix}.conversion.log
    echo "Sex encoding: 1=Male, 2=Female, 0=Unknown" >> ${prefix}.conversion.log
    """
}

// 收集所有PLINK文件用于合并
plink_final
    .collect()
    .set { all_plink_files }

// 提取共同变体并合并所有PLINK文件
process ExtractCommonVariantsAndMerge {
    executor 'slurm'
    queue 'gr10478b'
    time '24h'
    tag "Merging all PLINK files"
    publishDir "${params.OutputDir}/05.merged", mode: 'symlink'
    
    input:
    path(plink_files) from all_plink_files

    output:
    tuple path("cteph_agp3k.ajsa.sqc.norm.bed"), path("cteph_agp3k.ajsa.sqc.norm.bim"), path("cteph_agp3k.ajsa.sqc.norm.fam") into merged_plink
    path "merge.log"

    script:
    """
    # 提取所有.bim文件的路径
    bim_files=\$(ls -1 *.bim)
    
    echo "=== PLINK File Merging Report ===" > merge.log
    echo "Date: \$(date)" >> merge.log
    echo "" >> merge.log
    
    # 统计每个文件的变体数
    echo "Input files:" >> merge.log
    file_count=0
    for bim in \${bim_files}; do
        file_count=\$((file_count + 1))
        prefix=\${bim%.bim}
        variant_count=\$(wc -l < \${bim})
        sample_count=\$(wc -l < \${prefix}.fam)
        printf "  %d. %s: %'d variants, %'d samples\\n" \${file_count} \${prefix} \${variant_count} \${sample_count} >> merge.log
    done
    echo "" >> merge.log
    
    # 提取每个.bim文件的第2列（variant ID）
    echo "Extracting variant IDs from all files..." >> merge.log
    for bim in \${bim_files}; do
        awk '{print \$2}' \${bim} | sort > \${bim}.variants.txt
    done
    
    # 找到所有文件共有的变体（交集）
    echo "Finding common variants across all files..." >> merge.log
    
    # 获取第一个文件作为起点
    first_bim=\$(echo "\${bim_files}" | head -1)
    cp \${first_bim}.variants.txt common_variants.txt
    
    # 依次与其他文件求交集
    for bim in \${bim_files}; do
        if [ "\${bim}" != "\${first_bim}" ]; then
            comm -12 common_variants.txt \${bim}.variants.txt > temp_common.txt
            mv temp_common.txt common_variants.txt
        fi
    done
    
    common_count=\$(wc -l < common_variants.txt)
    printf "Common variants found: %'d\\n" \${common_count} >> merge.log
    echo "" >> merge.log
    
    # 为每个PLINK文件提取共同变体
    echo "Extracting common variants from each file..." >> merge.log
    first_filtered=""
    file_num=0
    for bim in \${bim_files}; do
        prefix=\${bim%.bim}
        
        ${params.Plink2Path} \\
            --bfile \${prefix} \\
            --extract common_variants.txt \\
            --make-bed \\
            --out \${prefix}.common \\
            --threads 4
        
        # 验证提取后的变体数
        extracted_count=\$(wc -l < \${prefix}.common.bim)
        printf "  %s: %'d variants extracted\\n" \${prefix} \${extracted_count} >> merge.log
        
        # 记录第一个文件，其他文件添加到merge列表
        file_num=\$((file_num + 1))
        if [ \${file_num} -eq 1 ]; then
            first_filtered=\${prefix}.common
        else
            echo "\${prefix}.common" >> merge_list.txt
        fi
    done
    echo "" >> merge.log
    
    # 检查文件数量并执行相应操作
    if [ \${file_num} -gt 1 ]; then
        # 多个文件，执行合并
        echo "Merging \${file_num} PLINK files using PLINK 1.9..." >> merge.log
        echo "Base file: \${first_filtered}" >> merge.log
        echo "Files to merge with base:" >> merge.log
        cat merge_list.txt >> merge.log
        echo "" >> merge.log
        
        ${params.PlinkPath} \\
            --bfile \${first_filtered} \\
            --merge-list merge_list.txt \\
            --keep-allele-order \\
            --allow-extra-chr \\
            --make-bed \\
            --out cteph_agp3k.ajsa.sqc.norm.plink1
        
        echo "Merge completed successfully with PLINK 1.9" >> merge.log
        echo "" >> merge.log
        
        # 使用PLINK2转换为PLINK2格式
        echo "Converting merged file to PLINK2 format..." >> merge.log
        ${params.Plink2Path} \\
            --bfile cteph_agp3k.ajsa.sqc.norm.plink1 \\
            --make-bed \\
            --out cteph_agp3k.ajsa.sqc.norm \\
            --threads 8
        
        echo "Conversion to PLINK2 format completed" >> merge.log
        
        # 清理PLINK1中间文件
        rm -f cteph_agp3k.ajsa.sqc.norm.plink1.bed
        rm -f cteph_agp3k.ajsa.sqc.norm.plink1.bim
        rm -f cteph_agp3k.ajsa.sqc.norm.plink1.fam
        rm -f cteph_agp3k.ajsa.sqc.norm.plink1.log
        
    elif [ \${file_num} -eq 1 ]; then
        # 只有一个文件，使用PLINK2转换
        echo "Only one file found, converting to PLINK2 format..." >> merge.log
        ${params.Plink2Path} \\
            --bfile \${first_filtered} \\
            --make-bed \\
            --out cteph_agp3k.ajsa.sqc.norm \\
            --threads 8
        
        echo "Conversion to PLINK2 format completed" >> merge.log
    else
        echo "ERROR: No files found to merge!" >> merge.log
        exit 1
    fi
    
    # 最终统计
    echo "" >> merge.log
    echo "======================================" >> merge.log
    echo "MERGE SUMMARY" >> merge.log
    echo "======================================" >> merge.log
    
    final_variant_count=\$(wc -l < cteph_agp3k.ajsa.sqc.norm.bim)
    final_sample_count=\$(wc -l < cteph_agp3k.ajsa.sqc.norm.fam)
    
    printf "Input files: %'d\\n" \${file_count} >> merge.log
    printf "Common variants: %'d\\n" \${common_count} >> merge.log
    printf "Final variants in merged file: %'d\\n" \${final_variant_count} >> merge.log
    printf "Total samples in merged file: %'d\\n" \${final_sample_count} >> merge.log
    echo "" >> merge.log
    
    # 统计合并后的性别分布
    echo "Sex distribution in merged file:" >> merge.log
    male_count=\$(awk '\$5==1' cteph_agp3k.ajsa.sqc.norm.fam | wc -l)
    female_count=\$(awk '\$5==2' cteph_agp3k.ajsa.sqc.norm.fam | wc -l)
    unknown_count=\$(awk '\$5==0' cteph_agp3k.ajsa.sqc.norm.fam | wc -l)
    printf "  Males: %'d\\n" \${male_count} >> merge.log
    printf "  Females: %'d\\n" \${female_count} >> merge.log
    printf "  Unknown: %'d\\n" \${unknown_count} >> merge.log
    echo "" >> merge.log
    
    # 统计合并后的表型分布（注意：此时表型尚未更新）
    echo "Phenotype distribution in merged file (before update):" >> merge.log
    case_count=\$(awk '\$6==2' cteph_agp3k.ajsa.sqc.norm.fam | wc -l)
    ctrl_count=\$(awk '\$6==1' cteph_agp3k.ajsa.sqc.norm.fam | wc -l)
    missing_pheno=\$(awk '\$6==-9 || \$6==0' cteph_agp3k.ajsa.sqc.norm.fam | wc -l)
    printf "  Cases (pheno=2): %'d\\n" \${case_count} >> merge.log
    printf "  Controls (pheno=1): %'d\\n" \${ctrl_count} >> merge.log
    printf "  Missing/Unknown: %'d\\n" \${missing_pheno} >> merge.log
    echo "  Note: Phenotype will be updated in the next process step" >> merge.log
    echo "" >> merge.log
    
    echo "Output files:" >> merge.log
    echo "  - cteph_agp3k.ajsa.sqc.norm.bed (merged binary genotype file)" >> merge.log
    echo "  - cteph_agp3k.ajsa.sqc.norm.bim (merged variant information)" >> merge.log
    echo "  - cteph_agp3k.ajsa.sqc.norm.fam (merged sample information)" >> merge.log
    echo "  - merge.log (this log file)" >> merge.log
    
    # 清理临时文件
    rm -f *.common.* *.variants.txt common_variants.txt merge_list.txt
    """
}

// 检查并去除Multiallelic variants
process RemoveMultiallelicVariants {
    executor 'slurm'
    queue 'gr10478b'
    time '12h'
    tag "Removing multiallelic variants"
    publishDir "${params.OutputDir}/06.biallelic_only", mode: 'symlink'
    
    input:
    tuple path(bed), path(bim), path(fam) from merged_plink

    output:
    tuple path("*.biallelic.bed"), path("*.biallelic.bim"), path("*.biallelic.fam") into biallelic_plink
    path "multiallelic_check.log"
    path "multiallelic_variants.txt" optional true

    script:
    def prefix = bed.baseName
    """
    echo "=== Multiallelic Variants Check Report ===" > multiallelic_check.log
    echo "Date: \$(date)" >> multiallelic_check.log
    echo "" >> multiallelic_check.log
    
    # 统计输入文件信息
    input_variant_count=\$(wc -l < ${prefix}.bim)
    input_sample_count=\$(wc -l < ${prefix}.fam)
    
    printf "Input file: %s\\n" "${prefix}" >> multiallelic_check.log
    printf "Input variants: %'d\\n" \${input_variant_count} >> multiallelic_check.log
    printf "Input samples: %'d\\n" \${input_sample_count} >> multiallelic_check.log
    echo "" >> multiallelic_check.log
    
    # 识别multiallelic variants（相同chr:pos的变体）
    echo "Identifying multiallelic variants (same chr:pos)..." >> multiallelic_check.log
    
    # 提取chr和pos（bim文件的第1列和第4列），统计重复
    awk '{print \$1":"\$4}' ${prefix}.bim | sort | uniq -c | awk '\$1 > 1 {print \$2}' > multiallelic_sites.txt
    
    multiallelic_site_count=\$(wc -l < multiallelic_sites.txt)
    printf "Multiallelic sites found: %'d\\n" \${multiallelic_site_count} >> multiallelic_check.log
    echo "" >> multiallelic_check.log
    
    if [ \${multiallelic_site_count} -gt 0 ]; then
        # 提取所有multiallelic位点的详细信息
        echo "Extracting details of multiallelic variants..." >> multiallelic_check.log
        
        # 创建multiallelic variant ID列表
        touch multiallelic_variants.txt
        while IFS= read -r site; do
            chr=\$(echo "\$site" | cut -d: -f1)
            pos=\$(echo "\$site" | cut -d: -f2)
            # 找到所有该位点的变体ID
            awk -v chr="\$chr" -v pos="\$pos" '\$1==chr && \$4==pos {print \$2}' ${prefix}.bim >> multiallelic_variants.txt
        done < multiallelic_sites.txt
        
        multiallelic_variant_count=\$(wc -l < multiallelic_variants.txt)
        printf "Total variants at multiallelic sites: %'d\\n" \${multiallelic_variant_count} >> multiallelic_check.log
        echo "" >> multiallelic_check.log
        
        # 显示前10个multiallelic位点的示例
        echo "Example multiallelic sites (first 10):" >> multiallelic_check.log
        head -10 multiallelic_sites.txt | while IFS= read -r site; do
            chr=\$(echo "\$site" | cut -d: -f1)
            pos=\$(echo "\$site" | cut -d: -f2)
            echo "  Site: \$site" >> multiallelic_check.log
            awk -v chr="\$chr" -v pos="\$pos" '\$1==chr && \$4==pos {printf "    - ID: %s, ALT: %s, REF: %s\\n", \$2, \$5, \$6}' ${prefix}.bim >> multiallelic_check.log
        done
        echo "" >> multiallelic_check.log
        
        # 使用PLINK去除multiallelic variants
        echo "Removing multiallelic variants with PLINK2..." >> multiallelic_check.log
        
        ${params.Plink2Path} \
            --bfile ${prefix} \
            --exclude multiallelic_variants.txt \
            --make-bed \
            --out ${prefix}.biallelic \
            --threads 8
        
        # 统计过滤后的结果
        output_variant_count=\$(wc -l < ${prefix}.biallelic.bim)
        removed_count=\$((input_variant_count - output_variant_count))
        
        printf "Variants after removal: %'d\\n" \${output_variant_count} >> multiallelic_check.log
        printf "Variants removed: %'d\\n" \${removed_count} >> multiallelic_check.log
        removal_rate=\$(awk "BEGIN {printf \\"%.2f%%\\", (\${removed_count}/\${input_variant_count})*100}")
        echo "Removal rate: \${removal_rate}" >> multiallelic_check.log
        
    else
        echo "No multiallelic variants detected!" >> multiallelic_check.log
        echo "All variants are biallelic. Copying files without modification..." >> multiallelic_check.log
        
        # 如果没有multiallelic variants，直接复制文件
        cp ${prefix}.bed ${prefix}.biallelic.bed
        cp ${prefix}.bim ${prefix}.biallelic.bim
        cp ${prefix}.fam ${prefix}.biallelic.fam
        
        printf "Output variants: %'d (unchanged)\\n" \${input_variant_count} >> multiallelic_check.log
    fi
    
    echo "" >> multiallelic_check.log
    echo "======================================" >> multiallelic_check.log
    echo "SUMMARY" >> multiallelic_check.log
    echo "======================================" >> multiallelic_check.log
    
    final_variant_count=\$(wc -l < ${prefix}.biallelic.bim)
    final_sample_count=\$(wc -l < ${prefix}.biallelic.fam)
    
    printf "Input variants: %'d\\n" \${input_variant_count} >> multiallelic_check.log
    printf "Output variants: %'d\\n" \${final_variant_count} >> multiallelic_check.log
    printf "Samples (unchanged): %'d\\n" \${final_sample_count} >> multiallelic_check.log
    echo "" >> multiallelic_check.log
    
    echo "Output files:" >> multiallelic_check.log
    echo "  - ${prefix}.biallelic.bed (binary genotype file, biallelic only)" >> multiallelic_check.log
    echo "  - ${prefix}.biallelic.bim (variant information, biallelic only)" >> multiallelic_check.log
    echo "  - ${prefix}.biallelic.fam (sample information, unchanged)" >> multiallelic_check.log
    echo "  - multiallelic_check.log (this log file)" >> multiallelic_check.log
    if [ \${multiallelic_site_count} -gt 0 ]; then
        echo "  - multiallelic_variants.txt (list of removed variant IDs)" >> multiallelic_check.log
    fi
    
    # 清理临时文件
    rm -f multiallelic_sites.txt
    """
}

// Variant QC 步骤1: 计算变异缺失率并过滤
process VariantQC_MissingRate {
    executor 'slurm'
    queue 'gr10478b'
    time '12h'
    tag "Variant missing rate QC"
    publishDir "${params.OutputDir}/07.variant_qc_missing", mode: 'symlink'
    
    input:
    tuple path(bed), path(bim), path(fam) from biallelic_plink

    output:
    tuple path("*.vmiss_qc.bed"), path("*.vmiss_qc.bim"), path("*.vmiss_qc.fam") into vmiss_qc_out
    path "*.vmiss"
    path "*.vmiss_qc.log"
    path "variant_missing_qc.log"

    script:
    def prefix = bed.baseName
    """
    echo "=== Variant Missing Rate QC Report ===" > variant_missing_qc.log
    echo "Date: \$(date)" >> variant_missing_qc.log
    echo "" >> variant_missing_qc.log
    
    # 统计输入文件信息
    input_variant_count=\$(wc -l < ${prefix}.bim)
    input_sample_count=\$(wc -l < ${prefix}.fam)
    
    printf "Input file: %s\\n" "${prefix}" >> variant_missing_qc.log
    printf "Input variants: %'d\\n" \${input_variant_count} >> variant_missing_qc.log
    printf "Input samples: %'d\\n" \${input_sample_count} >> variant_missing_qc.log
    echo "" >> variant_missing_qc.log
    
    # 统计样本性别分布（用于说明性染色体过滤策略）
    echo "Sample sex distribution:" >> variant_missing_qc.log
    male_count=\$(awk '\$5==1' ${prefix}.fam | wc -l)
    female_count=\$(awk '\$5==2' ${prefix}.fam | wc -l)
    unknown_count=\$(awk '\$5==0' ${prefix}.fam | wc -l)
    printf "  Males: %'d\\n" \${male_count} >> variant_missing_qc.log
    printf "  Females: %'d\\n" \${female_count} >> variant_missing_qc.log
    printf "  Unknown: %'d\\n" \${unknown_count} >> variant_missing_qc.log
    echo "" >> variant_missing_qc.log
    
    # 统计每条染色体的变异数
    echo "Variants per chromosome (before filtering):" >> variant_missing_qc.log
    cut -f1 ${prefix}.bim | sort | uniq -c | awk '{printf "  Chr %s: %'"'"'d variants\\n", \$2, \$1}' >> variant_missing_qc.log
    echo "" >> variant_missing_qc.log
    
    # 计算变异缺失率
    echo "Step 1: Calculating variant missing rates..." >> variant_missing_qc.log
    ${params.Plink2Path} \
        --bfile ${prefix} \
        --missing variant-only \
        --out ${prefix} \
        --threads 8
    
    # 分析vmiss文件
    total_variants=\$(tail -n +2 ${prefix}.vmiss | wc -l)
    echo "  Total variants analyzed: \${total_variants}" >> variant_missing_qc.log
    echo "" >> variant_missing_qc.log
    
    # 统计不同缺失率范围的变异数
    echo "Missing rate distribution:" >> variant_missing_qc.log
    awk 'NR>1 {
        missing_rate = \$5;
        if (missing_rate == 0) perfect++;
        else if (missing_rate <= 0.001) very_low++;
        else if (missing_rate <= 0.005) low++;
        else if (missing_rate <= 0.01) moderate++;
        else if (missing_rate <= 0.05) high++;
        else very_high++;
    }
    END {
        printf "  Missing rate = 0:          %'"'"'d (%.2f%%)\\n", perfect, (perfect/NR)*100;
        printf "  0 < missing rate ≤ 0.001:  %'"'"'d (%.2f%%)\\n", very_low, (very_low/NR)*100;
        printf "  0.001 < missing rate ≤ 0.005: %'"'"'d (%.2f%%)\\n", low, (low/NR)*100;
        printf "  0.005 < missing rate ≤ 0.01: %'"'"'d (%.2f%%)\\n", moderate, (moderate/NR)*100;
        printf "  0.01 < missing rate ≤ 0.05: %'"'"'d (%.2f%%)\\n", high, (high/NR)*100;
        printf "  Missing rate > 0.05:       %'"'"'d (%.2f%%)\\n", very_high, (very_high/NR)*100;
    }' ${prefix}.vmiss >> variant_missing_qc.log
    echo "" >> variant_missing_qc.log
    
    # 统计将被过滤的变异数（missing rate > 0.01，注意是严格大于）
    # PLINK2 --geno 使用 > 而非 >=，即保留 missing_rate <= threshold 的变异
    to_remove=\$(awk 'NR>1 && \$5 > 0.01' ${prefix}.vmiss | wc -l)
    printf "Variants with missing rate > 0.01: %'d (%.2f%%)\\n" \${to_remove} \$(awk "BEGIN {printf \\"%.2f\\", (\${to_remove}/\${total_variants})*100}") >> variant_missing_qc.log
    echo "" >> variant_missing_qc.log
    
    # PLINK2性染色体过滤说明
    echo "Note on sex chromosome filtering:" >> variant_missing_qc.log
    echo "  PLINK2 --geno filter handles sex chromosomes intelligently:" >> variant_missing_qc.log
    echo "  - Autosomes (chr1-22): Missing rate calculated across all samples" >> variant_missing_qc.log
    echo "  - X chromosome: Missing rate calculated only for female samples" >> variant_missing_qc.log
    echo "  - Y chromosome: Missing rate calculated only for male samples" >> variant_missing_qc.log
    echo "  - MT (mitochondrial): Missing rate calculated across all samples" >> variant_missing_qc.log
    echo "  This prevents incorrect filtering due to sex-specific biology." >> variant_missing_qc.log
    echo "" >> variant_missing_qc.log
    
    # 过滤缺失率 > 0.01 的变异
    echo "Step 2: Filtering variants with missing rate > 0.01..." >> variant_missing_qc.log
    ${params.Plink2Path} \
        --bfile ${prefix} \
        --geno 0.01 \
        --make-bed \
        --out ${prefix}.vmiss_qc \
        --threads 8
    
    # 统计过滤后的结果
    output_variant_count=\$(wc -l < ${prefix}.vmiss_qc.bim)
    removed_count=\$((input_variant_count - output_variant_count))
    
    printf "Variants after filtering: %'d\\n" \${output_variant_count} >> variant_missing_qc.log
    printf "Variants removed: %'d\\n" \${removed_count} >> variant_missing_qc.log
    removal_rate=\$(awk "BEGIN {printf \\"%.2f%%\\", (\${removed_count}/\${input_variant_count})*100}")
    echo "Removal rate: \${removal_rate}" >> variant_missing_qc.log
    retention_rate=\$(awk "BEGIN {printf \\"%.2f%%\\", (\${output_variant_count}/\${input_variant_count})*100}")
    echo "Retention rate: \${retention_rate}" >> variant_missing_qc.log
    echo "" >> variant_missing_qc.log
    
    # 统计每条染色体过滤后的变异数
    echo "Variants per chromosome (after filtering):" >> variant_missing_qc.log
    cut -f1 ${prefix}.vmiss_qc.bim | sort | uniq -c | awk '{printf "  Chr %s: %'"'"'d variants\\n", \$2, \$1}' >> variant_missing_qc.log
    echo "" >> variant_missing_qc.log
    
    # 详细统计每条染色体的过滤情况
    echo "Filtering details by chromosome:" >> variant_missing_qc.log
    for chr in \$(cut -f1 ${prefix}.bim | sort -u); do
        before=\$(awk -v c="\$chr" '\$1==c' ${prefix}.bim | wc -l)
        after=\$(awk -v c="\$chr" '\$1==c' ${prefix}.vmiss_qc.bim | wc -l)
        removed=\$((before - after))
        if [ \${before} -gt 0 ]; then
            pct=\$(awk "BEGIN {printf \\"%.2f\\", (\${removed}/\${before})*100}")
            printf "  Chr %s: %'d → %'d (removed: %'d, %.2f%%)\\n" "\$chr" \${before} \${after} \${removed} \${pct} >> variant_missing_qc.log
        fi
    done
    echo "" >> variant_missing_qc.log
    
    # 最终统计汇总
    echo "======================================" >> variant_missing_qc.log
    echo "SUMMARY" >> variant_missing_qc.log
    echo "======================================" >> variant_missing_qc.log
    printf "Input variants: %'d\\n" \${input_variant_count} >> variant_missing_qc.log
    printf "Output variants: %'d\\n" \${output_variant_count} >> variant_missing_qc.log
    printf "Variants removed: %'d\\n" \${removed_count} >> variant_missing_qc.log
    printf "Samples (unchanged): %'d\\n" \${input_sample_count} >> variant_missing_qc.log
    echo "" >> variant_missing_qc.log
    echo "Quality control threshold:" >> variant_missing_qc.log
    echo "  Maximum allowed missing rate: 1% (0.01, strict inequality)" >> variant_missing_qc.log
    echo "  Filter behavior: PLINK2 --geno uses > (not >=)" >> variant_missing_qc.log
    echo "    - Variants with missing_rate > 0.01 are EXCLUDED" >> variant_missing_qc.log
    echo "    - Variants with missing_rate <= 0.01 are RETAINED" >> variant_missing_qc.log
    echo "  Rationale: Variants with >1% missing genotypes may indicate" >> variant_missing_qc.log
    echo "             technical issues or poor probe performance" >> variant_missing_qc.log
    echo "" >> variant_missing_qc.log

    
    echo "Output files:" >> variant_missing_qc.log
    echo "  - ${prefix}.vmiss_qc.bed (binary genotype file, QC-passed)" >> variant_missing_qc.log
    echo "  - ${prefix}.vmiss_qc.bim (variant information, QC-passed)" >> variant_missing_qc.log
    echo "  - ${prefix}.vmiss_qc.fam (sample information, unchanged)" >> variant_missing_qc.log
    echo "  - ${prefix}.vmiss (variant missing rate statistics)" >> variant_missing_qc.log
    echo "  - ${prefix}.vmiss_qc.log (PLINK2 log file)" >> variant_missing_qc.log
    echo "  - variant_missing_qc.log (this report)" >> variant_missing_qc.log
    """
}

// 更新表型信息并拆分染色体区域
process UpdatePhenotypeAndSplitChromosomes {
    executor 'slurm'
    queue 'gr10478b'
    time '24h'
    tag "Update phenotype and split chromosomes"
    publishDir "${params.OutputDir}/08.pheno_update_chrsplit", mode: 'symlink'
    
    input:
    tuple path(bed), path(bim), path(fam) from vmiss_qc_out
    path(sample_info) from params.SampleInfo
    
    output:
    tuple val("autosomes"), path("*.autosomes.bed"), path("*.autosomes.bim"), path("*.autosomes.fam") optional true into chr_autosomes, chr_autosomes_2
    tuple val("chrX"), path("*.chrX.bed"), path("*.chrX.bim"), path("*.chrX.fam") optional true into chr_x
    tuple val("chrY"), path("*.chrY.bed"), path("*.chrY.bim"), path("*.chrY.fam") optional true into chr_y
    tuple val("PAR1"), path("*.PAR1.bed"), path("*.PAR1.bim"), path("*.PAR1.fam") optional true into chr_par1
    tuple val("PAR2"), path("*.PAR2.bed"), path("*.PAR2.bim"), path("*.PAR2.fam") optional true into chr_par2
    tuple val("chrM"), path("*.chrM.bed"), path("*.chrM.bim"), path("*.chrM.fam") optional true into chr_m
    path "pheno_update.log"
    path "*.update_pheno.txt"
    path "chr_split_summary.log"
    
    script:
    def prefix = bed.baseName
    def id_col = params.SampleIDColumn
    def outcome_col = params.OutcomeColumn
    def case_val = params.CaseValue
    """
    echo "=== Phenotype Update and Chromosome Split Report ===" > pheno_update.log
    echo "Date: \$(date)" >> pheno_update.log
    echo "" >> pheno_update.log
    
    # ===== 步骤1: 读取Excel并创建表型更新文件 =====
    echo "Step 1: Reading phenotype information from Excel..." >> pheno_update.log
    
    python3 << 'EOF'
import pandas as pd
import sys

print("Reading Excel file for phenotype information...")
try:
    df = pd.read_excel("${sample_info}")
    print(f"Total records in Excel: {len(df)}")
    
    # 检查必需的列是否存在
    if "${id_col}" not in df.columns:
        print(f"Error: Required column '${id_col}' not found in Excel file")
        print(f"Available columns: {', '.join(df.columns)}")
        sys.exit(1)
    
    if "${outcome_col}" not in df.columns:
        print(f"Error: Required column '${outcome_col}' not found in Excel file")
        print(f"Available columns: {', '.join(df.columns)}")
        sys.exit(1)
    
    # 提取ID和分组列
    pheno_data = df[["${id_col}", "${outcome_col}"]].copy()
    
    # 转换表型编码: case->2, control->1, missing/other->-9
    pheno_data['Pheno_Code'] = pheno_data["${outcome_col}"].apply(
        lambda x: '2' if str(x).strip().upper() == "${case_val}".upper() else 
                  ('1' if pd.notna(x) and str(x).strip() != '' else '-9')
    )
    
    # 统计表型分布
    case_count = (pheno_data['Pheno_Code'] == '2').sum()
    ctrl_count = (pheno_data['Pheno_Code'] == '1').sum()
    missing_count = (pheno_data['Pheno_Code'] == '-9').sum()
    
    print(f"Phenotype distribution:")
    print(f"  Cases (${case_val}): {case_count}")
    print(f"  Controls (other non-missing): {ctrl_count}")
    print(f"  Missing/Unknown: {missing_count}")
    
    # 创建PLINK格式的表型更新文件: FID IID Phenotype
    with open("${prefix}.update_pheno.txt", "w") as f:
        for _, row in pheno_data.iterrows():
            sample_id = str(row["${id_col}"])
            pheno_code = row['Pheno_Code']
            f.write(f"{sample_id}\\t{sample_id}\\t{pheno_code}\\n")
    
    print(f"Phenotype file created: ${prefix}.update_pheno.txt")
    
except Exception as e:
    print(f"Error reading Excel file: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)
EOF
    
    # 记录Python脚本输出到log
    echo "  Python output logged above" >> pheno_update.log
    echo "" >> pheno_update.log
    
    # ===== 步骤2: 更新PLINK fam文件的表型信息 =====
    echo "Step 2: Updating phenotype information in PLINK files..." >> pheno_update.log
    
    ${params.Plink2Path} \\
        --bfile ${prefix} \\
        --pheno ${prefix}.update_pheno.txt \\
        --make-bed \\
        --out ${prefix}.pheno \\
        --threads 8
    
    # 验证表型更新
    echo "  Phenotype update completed" >> pheno_update.log
    case_in_fam=\$(awk '\$6==2' ${prefix}.pheno.fam | wc -l)
    ctrl_in_fam=\$(awk '\$6==1' ${prefix}.pheno.fam | wc -l)
    missing_in_fam=\$(awk '\$6==-9 || \$6==0' ${prefix}.pheno.fam | wc -l)
    
    printf "  Cases in updated FAM: %'d\\n" \${case_in_fam} >> pheno_update.log
    printf "  Controls in updated FAM: %'d\\n" \${ctrl_in_fam} >> pheno_update.log
    printf "  Missing in updated FAM: %'d\\n" \${missing_in_fam} >> pheno_update.log
    echo "" >> pheno_update.log
    
    # ===== 步骤3: 分析染色体分布 =====
    echo "Step 3: Analyzing chromosome distribution..." >> pheno_update.log
    
    # 统计每条染色体的变异数，显示时添加"Chr"前缀
    echo "Chromosome distribution in input file:" >> pheno_update.log
    cut -f1 ${prefix}.pheno.bim | sort | uniq -c | \\
        awk '{printf "  Chr %s: %'"'"'d variants\\n", \$2, \$1}' >> pheno_update.log
    echo "" >> pheno_update.log
    
    # ===== 步骤4: 拆分染色体区域 =====
    echo "Step 4: Splitting genotype data by chromosomal regions..." >> pheno_update.log
    echo "" > chr_split_summary.log
    
    # 数据格式: 染色体标记为 1-22, X, Y, MT, PAR1, PAR2 (无chr前缀)
    echo "  Chromosome format: numeric without 'chr' prefix (1,2,...,22,X,Y,MT,PAR1,PAR2)" >> pheno_update.log
    echo "" >> pheno_update.log
    
    # 检查是否存在常染色体 (1-22)
    autosome_count=\$(awk '\$1 ~ /^([1-9]|1[0-9]|2[0-2])\$/' ${prefix}.pheno.bim | wc -l)
    if [ \${autosome_count} -gt 0 ]; then
        echo "  Extracting autosomes (1-22)..." >> pheno_update.log
        ${params.Plink2Path} \\
            --bfile ${prefix}.pheno \\
            --chr 1-22 \\
            --make-bed \\
            --out ${prefix}.autosomes \\
            --threads 8 2>&1 | tee -a chr_split_summary.log
        
        if [ -f ${prefix}.autosomes.bim ]; then
            variants=\$(wc -l < ${prefix}.autosomes.bim)
            printf "    Autosomes: %'d variants extracted\\n" \${variants} >> pheno_update.log
        fi
    else
        echo "    Autosomes: No variants found (skipped)" >> pheno_update.log
    fi
    
    # 检查是否存在X染色体 (标记为 "X"，不含PAR区域)
    # 注意: 如果数据中PAR1和PAR2是独立标记的，X染色体应该已经排除了PAR区域
    chrx_count=\$(awk '\$1 == "X"' ${prefix}.pheno.bim | wc -l)
    if [ \${chrx_count} -gt 0 ]; then
        echo "  Extracting chrX (non-PAR, marked as 'X')..." >> pheno_update.log
        ${params.Plink2Path} \\
            --bfile ${prefix}.pheno \\
            --chr X \\
            --make-bed \\
            --out ${prefix}.chrX \\
            --threads 8 2>&1 | tee -a chr_split_summary.log
        
        if [ -f ${prefix}.chrX.bim ]; then
            variants=\$(wc -l < ${prefix}.chrX.bim)
            printf "    chrX (non-PAR): %'d variants extracted\\n" \${variants} >> pheno_update.log
        fi
    else
        echo "    chrX: No variants found (skipped)" >> pheno_update.log
    fi
    
    # 检查是否存在Y染色体 (标记为 "Y"，不含PAR区域)
    # 注意: 如果数据中PAR1和PAR2是独立标记的，Y染色体应该已经排除了PAR区域
    chry_count=\$(awk '\$1 == "Y"' ${prefix}.pheno.bim | wc -l)
    if [ \${chry_count} -gt 0 ]; then
        echo "  Extracting chrY (non-PAR, marked as 'Y')..." >> pheno_update.log
        ${params.Plink2Path} \\
            --bfile ${prefix}.pheno \\
            --chr Y \\
            --make-bed \\
            --out ${prefix}.chrY \\
            --threads  8 2>&1 | tee -a chr_split_summary.log
        
        if [ -f ${prefix}.chrY.bim ]; then
            variants=\$(wc -l < ${prefix}.chrY.bim)
            printf "    chrY (non-PAR): %'d variants extracted\\n" \${variants} >> pheno_update.log
        fi
    else
        echo "    chrY: No variants found (skipped)" >> pheno_update.log
    fi
    
    # 检查是否存在PAR1区域 (直接标记为 "PAR1")
    par1_count=\$(awk '\$1 == "PAR1"' ${prefix}.pheno.bim | wc -l)
    if [ \${par1_count} -gt 0 ]; then
        echo "  Extracting PAR1 region (marked as 'PAR1')..." >> pheno_update.log
        ${params.Plink2Path} \\
            --bfile ${prefix}.pheno \\
            --chr PAR1 \\
            --make-bed \\
            --out ${prefix}.PAR1 \\
            --threads 8 2>&1 | tee -a chr_split_summary.log
        
        if [ -f ${prefix}.PAR1.bim ]; then
            variants=\$(wc -l < ${prefix}.PAR1.bim)
            printf "    PAR1: %'d variants extracted\\n" \${variants} >> pheno_update.log
        fi
    else
        echo "    PAR1: No variants found (skipped)" >> pheno_update.log
    fi
    
    # 检查是否存在PAR2区域 (直接标记为 "PAR2")
    par2_count=\$(awk '\$1 == "PAR2"' ${prefix}.pheno.bim | wc -l)
    if [ \${par2_count} -gt 0 ]; then
        echo "  Extracting PAR2 region (marked as 'PAR2')..." >> pheno_update.log
        ${params.Plink2Path} \\
            --bfile ${prefix}.pheno \\
            --chr PAR2 \\
            --make-bed \\
            --out ${prefix}.PAR2 \\
            --threads 8 2>&1 | tee -a chr_split_summary.log
        
        if [ -f ${prefix}.PAR2.bim ]; then
            variants=\$(wc -l < ${prefix}.PAR2.bim)
            printf "    PAR2: %'d variants extracted\\n" \${variants} >> pheno_update.log
        fi
    else
        echo "    PAR2: No variants found (skipped)" >> pheno_update.log
    fi
    
    # 检查是否存在线粒体染色体 (标记为 "MT")
    chrm_count=\$(awk '\$1 == "MT"' ${prefix}.pheno.bim | wc -l)
    if [ \${chrm_count} -gt 0 ]; then
        echo "  Extracting chrM (mitochondrial)..." >> pheno_update.log
        ${params.Plink2Path} \\
            --bfile ${prefix}.pheno \\
            --chr chrM \\
            --make-bed \\
            --out ${prefix}.chrM \\
            --threads 8 2>&1 | tee -a chr_split_summary.log
        
        if [ -f ${prefix}.chrM.bim ]; then
            variants=\$(wc -l < ${prefix}.chrM.bim)
            printf "    chrM: %'d variants extracted\\n" \${variants} >> pheno_update.log
        fi
    else
        echo "    chrM: No variants found (skipped)" >> pheno_update.log
    fi
    
    echo "" >> pheno_update.log
    
    # ===== 最终统计汇总 =====
    echo "======================================" >> pheno_update.log
    echo "SUMMARY" >> pheno_update.log
    echo "======================================" >> pheno_update.log
    
    input_variants=\$(wc -l < ${prefix}.pheno.bim)
    input_samples=\$(wc -l < ${prefix}.pheno.fam)
    
    printf "Input data:\\n" >> pheno_update.log
    printf "  Total variants: %'d\\n" \${input_variants} >> pheno_update.log
    printf "  Total samples: %'d\\n" \${input_samples} >> pheno_update.log
    printf "  Cases: %'d\\n" \${case_in_fam} >> pheno_update.log
    printf "  Controls: %'d\\n" \${ctrl_in_fam} >> pheno_update.log
    printf "  Missing: %'d\\n" \${missing_in_fam} >> pheno_update.log
    echo "" >> pheno_update.log
    
    printf "Chromosomal regions created:\\n" >> pheno_update.log
    for region in autosomes chrX chrY PAR1 PAR2 chrM; do
        if [ -f ${prefix}.\${region}.bim ]; then
            variants=\$(wc -l < ${prefix}.\${region}.bim)
            printf "  %s: %'d variants\\n" \${region} \${variants} >> pheno_update.log
        else
            printf "  %s: Not created (no variants)\\n" \${region} >> pheno_update.log
        fi
    done
    echo "" >> pheno_update.log
    
    echo "Output files:" >> pheno_update.log
    echo "  - ${prefix}.update_pheno.txt (phenotype update file)" >> pheno_update.log
    echo "  - ${prefix}.pheno.{bed,bim,fam} (updated phenotype files)" >> pheno_update.log
    echo "  - ${prefix}.<region>.{bed,bim,fam} (split chromosome files)" >> pheno_update.log
    echo "  - pheno_update.log (this report)" >> pheno_update.log
    echo "  - chr_split_summary.log (PLINK2 split logs)" >> pheno_update.log
    echo "" >> pheno_update.log
    
    echo "Note: Only regions with variants present in the input data are created." >> pheno_update.log
    echo "Missing regions indicate no variants were found for that chromosomal location." >> pheno_update.log
    
    # 清理临时文件
    rm -f ${prefix}.pheno.bed ${prefix}.pheno.bim ${prefix}.pheno.fam
    """
}

// Variant QC 步骤2: 计算AAF, MAF和HWE (仅针对常染色体)
process VariantQC_CalculateMetrics {
    executor 'slurm'
    queue 'gr10478b'
    time '12h'
    tag "Variant AAF/MAF/HWE calculation for autosomes"
    publishDir "${params.OutputDir}/09.variant_qc_stats", mode: 'symlink'
    
    input:
    tuple val(region), path(bed), path(bim), path(fam) from chr_autosomes
    
    output:
    tuple val(region), path("*.variant_qc_sum.tsv") into variant_qc_stats
    path "variant_qc_calculation.log"
    // path "*.all.afreq"
    // path "*.case.afreq"
    // path "*.ctrl.afreq"
    // path "*.all.hardy"
    // path "*.case.hardy"
    // path "*.ctrl.hardy"
    
    script:
    def prefix = bed.baseName
    """
    echo "=== Variant QC Statistics Calculation Report ===" > variant_qc_calculation.log
    echo "Date: \$(date)" >> variant_qc_calculation.log
    echo "Region: ${region}" >> variant_qc_calculation.log
    echo "" >> variant_qc_calculation.log
    
    # 统计输入文件信息
    input_variant_count=\$(wc -l < ${prefix}.bim)
    input_sample_count=\$(wc -l < ${prefix}.fam)
    case_count=\$(awk '\$6==2' ${prefix}.fam | wc -l)
    ctrl_count=\$(awk '\$6==1' ${prefix}.fam | wc -l)
    missing_pheno_count=\$(awk '\$6==-9 || \$6==0' ${prefix}.fam | wc -l)
    
    printf "Input file: %s\\n" "${prefix}" >> variant_qc_calculation.log
    printf "Total variants: %'d\\n" \${input_variant_count} >> variant_qc_calculation.log
    printf "Total samples: %'d\\n" \${input_sample_count} >> variant_qc_calculation.log
    printf "  Cases (pheno=2): %'d\\n" \${case_count} >> variant_qc_calculation.log
    printf "  Controls (pheno=1): %'d\\n" \${ctrl_count} >> variant_qc_calculation.log
    printf "  Missing phenotype: %'d\\n" \${missing_pheno_count} >> variant_qc_calculation.log
    echo "" >> variant_qc_calculation.log
    
    # ===== 步骤1: 计算AAF (ALL) =====
    echo "Step 1: Calculating ALT allele frequency (AAF) for ALL samples..." >> variant_qc_calculation.log
    ${params.Plink2Path} \\
        --bfile ${prefix} \\
        --freq \\
        --out ${prefix}.all \\
        --threads 8
    
    echo "  AAF (ALL) calculation completed" >> variant_qc_calculation.log
    echo "" >> variant_qc_calculation.log
    
    # ===== 步骤2: 计算AAF (CASE) =====
    echo "Step 2: Calculating ALT allele frequency (AAF) for CASE samples..." >> variant_qc_calculation.log
    ${params.Plink2Path} \\
        --bfile ${prefix} \\
        --keep <(awk '\$6==2 {print \$1, \$2}' ${prefix}.fam) \\
        --freq \\
        --out ${prefix}.case \\
        --threads 8
    
    echo "  AAF (CASE) calculation completed" >> variant_qc_calculation.log
    echo "" >> variant_qc_calculation.log
    
    # ===== 步骤3: 计算AAF (CTRL) =====
    echo "Step 3: Calculating ALT allele frequency (AAF) for CTRL samples..." >> variant_qc_calculation.log
    ${params.Plink2Path} \\
        --bfile ${prefix} \\
        --keep <(awk '\$6==1 {print \$1, \$2}' ${prefix}.fam) \\
        --freq \\
        --out ${prefix}.ctrl \\
        --threads 8
    
    echo "  AAF (CTRL) calculation completed" >> variant_qc_calculation.log
    echo "" >> variant_qc_calculation.log
    
    # ===== 步骤4: 计算HWE (ALL) =====
    echo "Step 4: Calculating Hardy-Weinberg equilibrium for ALL samples..." >> variant_qc_calculation.log
    ${params.Plink2Path} \\
        --bfile ${prefix} \\
        --hardy \\
        --out ${prefix}.all \\
        --threads 8
    
    echo "  HWE (ALL) calculation completed" >> variant_qc_calculation.log
    echo "" >> variant_qc_calculation.log
    
    # ===== 步骤5: 计算HWE (CASE) =====
    echo "Step 5: Calculating Hardy-Weinberg equilibrium for CASE samples..." >> variant_qc_calculation.log
    ${params.Plink2Path} \\
        --bfile ${prefix} \\
        --keep <(awk '\$6==2 {print \$1, \$2}' ${prefix}.fam) \\
        --hardy \\
        --out ${prefix}.case \\
        --threads 8
    
    echo "  HWE (CASE) calculation completed" >> variant_qc_calculation.log
    echo "" >> variant_qc_calculation.log
    
    # ===== 步骤6: 计算HWE (CTRL) =====
    echo "Step 6: Calculating Hardy-Weinberg equilibrium for CTRL samples..." >> variant_qc_calculation.log
    ${params.Plink2Path} \\
        --bfile ${prefix} \\
        --keep <(awk '\$6==1 {print \$1, \$2}' ${prefix}.fam) \\
        --hardy \\
        --out ${prefix}.ctrl \\
        --threads 8
    
    echo "  HWE (CTRL) calculation completed" >> variant_qc_calculation.log
    echo "" >> variant_qc_calculation.log
    
    # ===== 步骤7: 合并所有统计结果 =====
    echo "Step 7: Merging all statistics into summary file..." >> variant_qc_calculation.log
    
    python3 << 'EOF'
import pandas as pd
import sys

print("Reading PLINK2 output files...")

try:
    # 读取AAF文件 (PLINK2 --freq 输出格式)
    # 列: #CHROM ID REF ALT ALT_FREQS OBS_CT
    afreq_all = pd.read_csv("${prefix}.all.afreq", sep="\\t")
    afreq_case = pd.read_csv("${prefix}.case.afreq", sep="\\t")
    afreq_ctrl = pd.read_csv("${prefix}.ctrl.afreq", sep="\\t")
    
    print(f"  AAF ALL: {len(afreq_all)} variants")
    print(f"  AAF CASE: {len(afreq_case)} variants")
    print(f"  AAF CTRL: {len(afreq_ctrl)} variants")
    
    # 读取HWE文件 (PLINK2 --hardy 输出格式)
    # 列: #CHROM ID A1 AX HETX HOMX1 HOMX2 HET_A1 HOMXAX P
    hardy_all = pd.read_csv("${prefix}.all.hardy", sep="\\t")
    hardy_case = pd.read_csv("${prefix}.case.hardy", sep="\\t")
    hardy_ctrl = pd.read_csv("${prefix}.ctrl.hardy", sep="\\t")
    
    print(f"  HWE ALL: {len(hardy_all)} variants")
    print(f"  HWE CASE: {len(hardy_case)} variants")
    print(f"  HWE CTRL: {len(hardy_ctrl)} variants")
    
    # 提取需要的列并重命名
    # AAF = ALT allele frequency (即ALT_FREQS列)
    afreq_all_sub = afreq_all[["ID", "ALT_FREQS"]].rename(columns={"ALT_FREQS": "AAF_ALL"})
    afreq_case_sub = afreq_case[["ID", "ALT_FREQS"]].rename(columns={"ALT_FREQS": "AAF_CASE"})
    afreq_ctrl_sub = afreq_ctrl[["ID", "ALT_FREQS"]].rename(columns={"ALT_FREQS": "AAF_CTRL"})
    
    # HWE = P value (即P列)
    hardy_all_sub = hardy_all[["ID", "P"]].rename(columns={"P": "HWE_ALL"})
    hardy_case_sub = hardy_case[["ID", "P"]].rename(columns={"P": "HWE_CASE"})
    hardy_ctrl_sub = hardy_ctrl[["ID", "P"]].rename(columns={"P": "HWE_CTRL"})
    
    # 合并所有AAF数据
    merged = afreq_all_sub.merge(afreq_case_sub, on="ID", how="outer")
    merged = merged.merge(afreq_ctrl_sub, on="ID", how="outer")
    
    # 计算MAF (Minor Allele Frequency)
    # MAF = min(AAF, 1-AAF)
    merged["MAF_ALL"] = merged["AAF_ALL"].apply(lambda x: min(x, 1-x) if pd.notna(x) else None)
    merged["MAF_CASE"] = merged["AAF_CASE"].apply(lambda x: min(x, 1-x) if pd.notna(x) else None)
    merged["MAF_CTRL"] = merged["AAF_CTRL"].apply(lambda x: min(x, 1-x) if pd.notna(x) else None)
    
    # 合并HWE数据
    merged = merged.merge(hardy_all_sub, on="ID", how="outer")
    merged = merged.merge(hardy_case_sub, on="ID", how="outer")
    merged = merged.merge(hardy_ctrl_sub, on="ID", how="outer")
    
    # 重新排列列顺序
    merged = merged[["ID", "AAF_ALL", "AAF_CASE", "AAF_CTRL", 
                     "MAF_ALL", "MAF_CASE", "MAF_CTRL",
                     "HWE_ALL", "HWE_CASE", "HWE_CTRL"]]
    
    # 重命名ID列为VARIANT_ID
    merged = merged.rename(columns={"ID": "VARIANT_ID"})
    
    # 解析VARIANT_ID格式 (CHROM:POS:REF:ALT) 用于排序
    print("\\nParsing VARIANT_ID for sorting...")
    merged[["CHROM", "POS", "REF", "ALT"]] = merged["VARIANT_ID"].str.split(":", expand=True)
    
    # 转换POS为整数
    merged["POS"] = merged["POS"].astype(int)
    
    # 提取染色体编号用于排序 (去除chr前缀)
    # chr1 -> 1, chr2 -> 2, ..., chr22 -> 22, chrX -> 23, chrY -> 24, chrM -> 25
    def get_chr_sort_key(chrom):
        chrom_clean = chrom.replace("chr", "").replace("Chr", "").replace("CHR", "")
        if chrom_clean.isdigit():
            return int(chrom_clean)
        elif chrom_clean.upper() == "X":
            return 23
        elif chrom_clean.upper() == "Y":
            return 24
        elif chrom_clean.upper() in ["M", "MT"]:
            return 25
        elif chrom_clean.upper() == "PAR1":
            return 26
        elif chrom_clean.upper() == "PAR2":
            return 27
        else:
            return 99  # 其他未知染色体排在最后
    
    merged["CHR_SORT_KEY"] = merged["CHROM"].apply(get_chr_sort_key)
    
    # 按照染色体编号(主升序)和位置(次升序)排序
    merged = merged.sort_values(by=["CHR_SORT_KEY", "POS"])
    
    # 删除临时排序列
    merged = merged.drop(columns=["CHROM", "POS", "REF", "ALT", "CHR_SORT_KEY"])
    
    print(f"  Sorted by chromosome (chr1-chr22) and position")
    
    # 保存为TSV文件
    output_file = "${prefix}.variant_qc_sum.tsv"
    merged.to_csv(output_file, sep="\\t", index=False, na_rep="NA")
    
    print(f"  Summary file created: {output_file}")
    print(f"  Total variants in summary: {len(merged)}")
    
    # 统计一些基本信息
    print("\\nStatistics summary:")
    print(f"  AAF_ALL range: [{merged['AAF_ALL'].min():.4f}, {merged['AAF_ALL'].max():.4f}]")
    print(f"  MAF_ALL range: [{merged['MAF_ALL'].min():.4f}, {merged['MAF_ALL'].max():.4f}]")
    print(f"  HWE_ALL range: [{merged['HWE_ALL'].min():.4e}, {merged['HWE_ALL'].max():.4e}]")
    
    # 统计MAF和AAF不同的变异数 (ALT allele是major allele的情况)
    # MAF != AAF 意味着 AAF > 0.5，即ALT allele是major allele
    maf_aaf_diff_all = (merged['MAF_ALL'] != merged['AAF_ALL']).sum()
    maf_aaf_diff_case = (merged['MAF_CASE'] != merged['AAF_CASE']).sum()
    maf_aaf_diff_ctrl = (merged['MAF_CTRL'] != merged['AAF_CTRL']).sum()
    
    print(f"\\nVariants where MAF != AAF (ALT allele is major allele, AAF > 0.5):")
    print(f"  ALL: {maf_aaf_diff_all} variants ({maf_aaf_diff_all/len(merged)*100:.2f}%)")
    print(f"  CASE: {maf_aaf_diff_case} variants ({maf_aaf_diff_case/len(merged)*100:.2f}%)")
    print(f"  CTRL: {maf_aaf_diff_ctrl} variants ({maf_aaf_diff_ctrl/len(merged)*100:.2f}%)")
    print("  Note: MAF = min(AAF, 1-AAF), so MAF != AAF when ALT is the major allele")
    
    # 统计HWE显著偏离的变异数 (P < 1e-6)
    hwe_fail_all = (merged['HWE_ALL'] < 1e-6).sum()
    hwe_fail_case = (merged['HWE_CASE'] < 1e-6).sum()
    hwe_fail_ctrl = (merged['HWE_CTRL'] < 1e-6).sum()
    
    print(f"\\nHWE violations (P < 1e-6):")
    print(f"  ALL: {hwe_fail_all} variants ({hwe_fail_all/len(merged)*100:.2f}%)")
    print(f"  CASE: {hwe_fail_case} variants ({hwe_fail_case/len(merged)*100:.2f}%)")
    print(f"  CTRL: {hwe_fail_ctrl} variants ({hwe_fail_ctrl/len(merged)*100:.2f}%)")
    
except Exception as e:
    print(f"Error processing files: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)
EOF
    
    echo "  Statistics merged successfully" >> variant_qc_calculation.log
    echo "" >> variant_qc_calculation.log
    
    # ===== 最终统计汇总 =====
    echo "======================================" >> variant_qc_calculation.log
    echo "SUMMARY" >> variant_qc_calculation.log
    echo "======================================" >> variant_qc_calculation.log
    
    printf "Region: %s\\n" "${region}" >> variant_qc_calculation.log
    printf "Total variants: %'d\\n" \${input_variant_count} >> variant_qc_calculation.log
    printf "Samples used:\\n" >> variant_qc_calculation.log
    printf "  ALL: %'d\\n" \${input_sample_count} >> variant_qc_calculation.log
    printf "  CASE: %'d\\n" \${case_count} >> variant_qc_calculation.log
    printf "  CTRL: %'d\\n" \${ctrl_count} >> variant_qc_calculation.log
    echo "" >> variant_qc_calculation.log
    
    echo "Statistics calculated:" >> variant_qc_calculation.log
    echo "  1. AAF (ALT allele frequency): Frequency of the ALT allele (bim column 5)" >> variant_qc_calculation.log
    echo "  2. MAF (Minor allele frequency): min(AAF, 1-AAF)" >> variant_qc_calculation.log
    echo "  3. HWE (Hardy-Weinberg equilibrium): P-value for deviation from HWE" >> variant_qc_calculation.log
    echo "" >> variant_qc_calculation.log
    echo "Important notes:" >> variant_qc_calculation.log
    echo "  - AAF = Frequency of ALT allele (can be >0.5 if ALT is the major allele)" >> variant_qc_calculation.log
    echo "  - MAF = Frequency of the less common allele (always ≤0.5)" >> variant_qc_calculation.log
    echo "  - When AAF > 0.5: ALT is the major allele, MAF = 1 - AAF" >> variant_qc_calculation.log
    echo "  - When AAF ≤ 0.5: ALT is the minor allele, MAF = AAF" >> variant_qc_calculation.log
    echo "  - MAF != AAF indicates ALT allele is more common than REF allele" >> variant_qc_calculation.log
    echo "" >> variant_qc_calculation.log
    echo "Stratification:" >> variant_qc_calculation.log
    echo "  - ALL: All samples included" >> variant_qc_calculation.log
    echo "  - CASE: Only case samples (phenotype = 2)" >> variant_qc_calculation.log
    echo "  - CTRL: Only control samples (phenotype = 1)" >> variant_qc_calculation.log
    echo "" >> variant_qc_calculation.log
    
    echo "Output files:" >> variant_qc_calculation.log
    echo "  - ${prefix}.variant_qc_sum.tsv (main summary file)" >> variant_qc_calculation.log
    echo "      Columns: VARIANT_ID, AAF_ALL, AAF_CASE, AAF_CTRL," >> variant_qc_calculation.log
    echo "               MAF_ALL, MAF_CASE, MAF_CTRL," >> variant_qc_calculation.log
    echo "               HWE_ALL, HWE_CASE, HWE_CTRL" >> variant_qc_calculation.log
    echo "  - ${prefix}.all.afreq (PLINK2 frequency output, ALL)" >> variant_qc_calculation.log
    echo "  - ${prefix}.case.afreq (PLINK2 frequency output, CASE)" >> variant_qc_calculation.log
    echo "  - ${prefix}.ctrl.afreq (PLINK2 frequency output, CTRL)" >> variant_qc_calculation.log
    echo "  - ${prefix}.all.hardy (PLINK2 HWE output, ALL)" >> variant_qc_calculation.log
    echo "  - ${prefix}.case.hardy (PLINK2 HWE output, CASE)" >> variant_qc_calculation.log
    echo "  - ${prefix}.ctrl.hardy (PLINK2 HWE output, CTRL)" >> variant_qc_calculation.log
    echo "  - variant_qc_calculation.log (this report)" >> variant_qc_calculation.log
    """
}

// Variant QC 步骤3: 根据MAF和HWE过滤变异
process VariantQC_FilterByMAFandHWE {
    executor 'slurm'
    queue 'gr10478b'
    time '12h'
    tag "Filtering variants by MAF and HWE"
    publishDir "${params.OutputDir}/10.variant_qc_filtered", mode: 'symlink'
    
    input:
    tuple val(region), path(tsv) from variant_qc_stats
    tuple val(region2), path(bed), path(bim), path(fam) from chr_autosomes_2
    
    output:
    tuple val(region), path("cteph_agp3k.array.sqc.vqc.bed"), path("cteph_agp3k.array.sqc.vqc.bim"), path("cteph_agp3k.array.sqc.vqc.fam") into qc_filtered_plink
    path "*.variants_to_remove.txt"
    path "variant_filter.log"
    
    script:
    def prefix = bed.baseName
    def output_prefix = "cteph_agp3k.array.sqc.vqc"
    def maf_mode = params.MAF_FilterMode ?: "ALL"
    def maf_threshold = params.MAF_Threshold ?: 0.01
    def hwe_mode = params.HWE_FilterMode ?: "CTRL|CASE"
    def hwe_ctrl_threshold = params.HWE_CTRL_Threshold ?: 1e-6
    def hwe_case_threshold = params.HWE_CASE_Threshold ?: 1e-10
    """
    echo "=== Variant Filtering by MAF and HWE Report ===" > variant_filter.log
    echo "Date: \$(date)" >> variant_filter.log
    echo "Region: ${region}" >> variant_filter.log
    echo "" >> variant_filter.log
    
    # 统计输入文件信息
    input_variant_count=\$(wc -l < ${prefix}.bim)
    input_sample_count=\$(wc -l < ${prefix}.fam)
    
    printf "Input file: %s\\n" "${prefix}" >> variant_filter.log
    printf "Total variants: %'d\\n" \${input_variant_count} >> variant_filter.log
    printf "Total samples: %'d\\n" \${input_sample_count} >> variant_filter.log
    echo "" >> variant_filter.log
    
    # 过滤条件
    echo "Filtering criteria:" >> variant_filter.log
    echo "  MAF filter mode: ${maf_mode}" >> variant_filter.log
    echo "  MAF threshold: MAF_${maf_mode} < ${maf_threshold}" >> variant_filter.log
    echo "  HWE filter mode: ${hwe_mode}" >> variant_filter.log
    if [[ "${hwe_mode}" == *"CTRL"* ]]; then
        echo "  HWE CTRL threshold: HWE_CTRL < ${hwe_ctrl_threshold}" >> variant_filter.log
    fi
    if [[ "${hwe_mode}" == *"CASE"* ]]; then
        echo "  HWE CASE threshold: HWE_CASE < ${hwe_case_threshold}" >> variant_filter.log
    fi
    echo "" >> variant_filter.log
    
    # 使用Python提取需要移除的变异ID
    echo "Step 1: Identifying variants to remove..." >> variant_filter.log
    
    python3 << 'EOF' 2>&1 | tee -a variant_filter.log
import pandas as pd
import sys

print("Reading variant QC summary file...")

try:
    # 读取TSV文件
    df = pd.read_csv("${tsv}", sep="\\t")
    print(f"  Total variants in TSV: {len(df)}")
    
    # 检查所需的列是否存在
    maf_col = "MAF_${maf_mode}"
    hwe_mode = "${hwe_mode}"
    
    if maf_col not in df.columns:
        print(f"Error: Column '{maf_col}' not found in TSV file")
        print(f"Available columns: {', '.join(df.columns)}")
        sys.exit(1)
    
    # 检查HWE相关列
    use_hwe_ctrl = "CTRL" in hwe_mode
    use_hwe_case = "CASE" in hwe_mode
    
    if use_hwe_ctrl and "HWE_CTRL" not in df.columns:
        print(f"Error: Column 'HWE_CTRL' not found in TSV file")
        print(f"Available columns: {', '.join(df.columns)}")
        sys.exit(1)
    
    if use_hwe_case and "HWE_CASE" not in df.columns:
        print(f"Error: Column 'HWE_CASE' not found in TSV file")
        print(f"Available columns: {', '.join(df.columns)}")
        sys.exit(1)
    
    # 应用过滤条件
    maf_threshold = float("${maf_threshold}")
    hwe_ctrl_threshold = float("${hwe_ctrl_threshold}")
    hwe_case_threshold = float("${hwe_case_threshold}")
    
    # MAF过滤
    maf_fail = df[maf_col] < maf_threshold
    
    # HWE过滤（根据模式）
    hwe_fail = pd.Series([False] * len(df), index=df.index)
    hwe_ctrl_fail = pd.Series([False] * len(df), index=df.index)
    hwe_case_fail = pd.Series([False] * len(df), index=df.index)
    
    if use_hwe_ctrl:
        hwe_ctrl_fail = df["HWE_CTRL"] < hwe_ctrl_threshold
        hwe_fail = hwe_fail | hwe_ctrl_fail
    
    if use_hwe_case:
        hwe_case_fail = df["HWE_CASE"] < hwe_case_threshold
        hwe_fail = hwe_fail | hwe_case_fail
    
    # 合并所有过滤条件
    to_remove = df[maf_fail | hwe_fail]
    
    print(f"\\nFiltering results:")
    print(f"  Variants failing MAF filter ({maf_col} < {maf_threshold}): {maf_fail.sum()} ({maf_fail.sum()/len(df)*100:.2f}%)")
    
    if use_hwe_ctrl:
        print(f"  Variants failing HWE_CTRL filter (< {hwe_ctrl_threshold}): {hwe_ctrl_fail.sum()} ({hwe_ctrl_fail.sum()/len(df)*100:.2f}%)")
    
    if use_hwe_case:
        print(f"  Variants failing HWE_CASE filter (< {hwe_case_threshold}): {hwe_case_fail.sum()} ({hwe_case_fail.sum()/len(df)*100:.2f}%)")
    
    print(f"  Variants failing any filter (to remove): {len(to_remove)} ({len(to_remove)/len(df)*100:.2f}%)")
    print(f"  Variants passing all filters (to keep): {len(df) - len(to_remove)} ({(len(df)-len(to_remove))/len(df)*100:.2f}%)")
    
    # 保存需要移除的变异ID（不包含列名）
    to_remove["VARIANT_ID"].to_csv("${prefix}.variants_to_remove.txt", index=False, header=False)
    
    print(f"\\nVariant IDs to remove saved to: ${prefix}.variants_to_remove.txt")
    
    # 统计失败组合
    print(f"\\nFailure breakdown:")
    
    if use_hwe_ctrl and use_hwe_case:
        # 三个过滤器都启用
        maf_only = maf_fail & ~hwe_ctrl_fail & ~hwe_case_fail
        hwe_ctrl_only = ~maf_fail & hwe_ctrl_fail & ~hwe_case_fail
        hwe_case_only = ~maf_fail & ~hwe_ctrl_fail & hwe_case_fail
        maf_hwe_ctrl = maf_fail & hwe_ctrl_fail & ~hwe_case_fail
        maf_hwe_case = maf_fail & ~hwe_ctrl_fail & hwe_case_fail
        hwe_both = ~maf_fail & hwe_ctrl_fail & hwe_case_fail
        all_three = maf_fail & hwe_ctrl_fail & hwe_case_fail
        
        print(f"  MAF only: {maf_only.sum()} variants ({maf_only.sum()/len(df)*100:.2f}%)")
        print(f"  HWE_CTRL only: {hwe_ctrl_only.sum()} variants ({hwe_ctrl_only.sum()/len(df)*100:.2f}%)")
        print(f"  HWE_CASE only: {hwe_case_only.sum()} variants ({hwe_case_only.sum()/len(df)*100:.2f}%)")
        print(f"  MAF + HWE_CTRL: {maf_hwe_ctrl.sum()} variants ({maf_hwe_ctrl.sum()/len(df)*100:.2f}%)")
        print(f"  MAF + HWE_CASE: {maf_hwe_case.sum()} variants ({maf_hwe_case.sum()/len(df)*100:.2f}%)")
        print(f"  HWE_CTRL + HWE_CASE: {hwe_both.sum()} variants ({hwe_both.sum()/len(df)*100:.2f}%)")
        print(f"  All three filters: {all_three.sum()} variants ({all_three.sum()/len(df)*100:.2f}%)")
    
    elif use_hwe_ctrl:
        # 只有MAF和HWE_CTRL
        maf_only = maf_fail & ~hwe_ctrl_fail
        hwe_ctrl_only = ~maf_fail & hwe_ctrl_fail
        both = maf_fail & hwe_ctrl_fail
        
        print(f"  MAF only: {maf_only.sum()} variants ({maf_only.sum()/len(df)*100:.2f}%)")
        print(f"  HWE_CTRL only: {hwe_ctrl_only.sum()} variants ({hwe_ctrl_only.sum()/len(df)*100:.2f}%)")
        print(f"  Both MAF + HWE_CTRL: {both.sum()} variants ({both.sum()/len(df)*100:.2f}%)")
    
    elif use_hwe_case:
        # 只有MAF和HWE_CASE
        maf_only = maf_fail & ~hwe_case_fail
        hwe_case_only = ~maf_fail & hwe_case_fail
        both = maf_fail & hwe_case_fail
        
        print(f"  MAF only: {maf_only.sum()} variants ({maf_only.sum()/len(df)*100:.2f}%)")
        print(f"  HWE_CASE only: {hwe_case_only.sum()} variants ({hwe_case_only.sum()/len(df)*100:.2f}%)")
        print(f"  Both MAF + HWE_CASE: {both.sum()} variants ({both.sum()/len(df)*100:.2f}%)")
    
    else:
        # 只有MAF过滤
        print(f"  MAF only: {maf_fail.sum()} variants ({maf_fail.sum()/len(df)*100:.2f}%)")
    
except Exception as e:
    print(f"Error processing TSV file: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)
EOF
    
    echo "" >> variant_filter.log
    
    # 统计要移除的变异数
    variants_to_remove=\$(wc -l < ${prefix}.variants_to_remove.txt)
    printf "Variants identified for removal: %'d\\n" \${variants_to_remove} >> variant_filter.log
    echo "" >> variant_filter.log
    
    # 使用PLINK2移除这些变异
    if [ \${variants_to_remove} -gt 0 ]; then
        echo "Step 2: Removing filtered variants with PLINK2..." >> variant_filter.log
        
        ${params.Plink2Path} \\
            --bfile ${prefix} \\
            --exclude ${prefix}.variants_to_remove.txt \\
            --make-bed \\
            --out ${output_prefix} \\
            --threads 8
        
        # 统计过滤后的结果
        output_variant_count=\$(wc -l < ${output_prefix}.bim)
        removed_count=\$((input_variant_count - output_variant_count))
        
        echo "" >> variant_filter.log
        printf "Variants after filtering: %'d\\n" \${output_variant_count} >> variant_filter.log
        printf "Variants removed: %'d\\n" \${removed_count} >> variant_filter.log
        removal_rate=\$(awk "BEGIN {printf \\"%.2f%%\\", (\${removed_count}/\${input_variant_count})*100}")
        echo "Removal rate: \${removal_rate}" >> variant_filter.log
        retention_rate=\$(awk "BEGIN {printf \\"%.2f%%\\", (\${output_variant_count}/\${input_variant_count})*100}")
        echo "Retention rate: \${retention_rate}" >> variant_filter.log
        
    else
        echo "Step 2: No variants to remove (all passed filters)" >> variant_filter.log
        echo "  Copying input files to output..." >> variant_filter.log
        
        cp ${prefix}.bed ${output_prefix}.bed
        cp ${prefix}.bim ${output_prefix}.bim
        cp ${prefix}.fam ${output_prefix}.fam
        
        output_variant_count=\${input_variant_count}
        echo "  All variants retained: \${output_variant_count}" >> variant_filter.log
    fi
    
    echo "" >> variant_filter.log
    
    # ===== 最终统计汇总 =====
    echo "======================================" >> variant_filter.log
    echo "SUMMARY" >> variant_filter.log
    echo "======================================" >> variant_filter.log
    
    printf "Region: %s\\n" "${region}" >> variant_filter.log
    printf "Input variants: %'d\\n" \${input_variant_count} >> variant_filter.log
    printf "Output variants: %'d\\n" \${output_variant_count} >> variant_filter.log
    
    if [ \${variants_to_remove} -gt 0 ]; then
        removed_count=\$((input_variant_count - output_variant_count))
        printf "Variants removed: %'d\\n" \${removed_count} >> variant_filter.log
        removal_rate=\$(awk "BEGIN {printf \\"%.2f%%\\", (\${removed_count}/\${input_variant_count})*100}")
        retention_rate=\$(awk "BEGIN {printf \\"%.2f%%\\", (\${output_variant_count}/\${input_variant_count})*100}")
        echo "  Removal rate: \${removal_rate}" >> variant_filter.log
        echo "  Retention rate: \${retention_rate}" >> variant_filter.log
    else
        echo "Variants removed: 0 (all passed filters)" >> variant_filter.log
    fi
    
    printf "Samples (unchanged): %'d\\n" \${input_sample_count} >> variant_filter.log
    echo "" >> variant_filter.log
    
    echo "Quality control filters applied:" >> variant_filter.log
    echo "  1. MAF filter: MAF_${maf_mode} < ${maf_threshold}" >> variant_filter.log
    echo "     Rationale: Variants with very low minor allele frequency may have" >> variant_filter.log
    echo "                insufficient statistical power for association testing" >> variant_filter.log
    echo "" >> variant_filter.log
    
    filter_num=2
    if [[ "${hwe_mode}" == *"CTRL"* ]]; then
        echo "  \${filter_num}. HWE_CTRL filter: HWE_CTRL < ${hwe_ctrl_threshold}" >> variant_filter.log
        echo "     Rationale: Significant deviation from Hardy-Weinberg equilibrium in" >> variant_filter.log
        echo "                controls may indicate genotyping errors or population stratification" >> variant_filter.log
        echo "     Note: Using controls avoids removing disease-associated variants" >> variant_filter.log
        echo "" >> variant_filter.log
        filter_num=\$((filter_num + 1))
    fi
    
    if [[ "${hwe_mode}" == *"CASE"* ]]; then
        echo "  \${filter_num}. HWE_CASE filter: HWE_CASE < ${hwe_case_threshold}" >> variant_filter.log
        echo "     Rationale: Extreme deviation from HWE in cases may indicate" >> variant_filter.log
        echo "                genotyping errors or technical artifacts" >> variant_filter.log
        echo "" >> variant_filter.log
    fi
    
    echo "Filter logic: Variants are REMOVED if they fail ANY criterion" >> variant_filter.log
    echo "              (i.e., OR logic among all active filters)" >> variant_filter.log
    echo "" >> variant_filter.log
    echo "Active filter modes:" >> variant_filter.log
    echo "  - MAF mode: ${maf_mode} (filtering based on MAF_${maf_mode})" >> variant_filter.log
    echo "  - HWE mode: ${hwe_mode} (filtering based on ${hwe_mode} samples)" >> variant_filter.log
    echo "" >> variant_filter.log
    
    echo "Output files:" >> variant_filter.log
    echo "  - ${output_prefix}.bed (filtered binary genotype file)" >> variant_filter.log
    echo "  - ${output_prefix}.bim (filtered variant information)" >> variant_filter.log
    echo "  - ${output_prefix}.fam (sample information, unchanged)" >> variant_filter.log
    echo "  - ${prefix}.variants_to_remove.txt (list of removed variant IDs)" >> variant_filter.log
    echo "  - variant_filter.log (this report)" >> variant_filter.log
    """
}




