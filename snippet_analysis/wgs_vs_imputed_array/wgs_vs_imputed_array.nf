nextflow.enable.dsl=1

params.ImputedArrayPath = '/LARGE0/gr10478/project/Pulmonary_Hypertension/DATA_SOURCE/impute20251029/JSA'
params.SciptsDir = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/snippet_analysis/wgs_vs_imputed_array/scripts'
params.OutDir = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/snippet_analysis/wgs_vs_imputed_array/results'
params.SampleSelectList = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/array_rev1/info/cteph_agp3k.v4.ls'
params.NagasakiPipelinePath = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/nagasaki_pipeline'
params.Bcftools = '/home/b/b37974/bcftools/bcftools'
params.Plink = '/home/b/b37974/plink'
params.Plink2 = '/home/b/b37974/plink2'
params.SampleInfo = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/info/cteph_agp3k_jhrpv4.rev1.xlsx'
params.SampleIDColumn = 'ID'
params.SampleSexColumn = 'Sex'
params.OutcomeColumn = 'OUTCOME2'
params.CaseValue = 'CTEPH'

// Create channel for chr1-chr22 imputed array BCF files
imputed_array_ch = Channel
    .from(1..22)
    .map { chr -> 
        def chrName = "chr${chr}"
        def bcfFile = file("${params.ImputedArrayPath}/impute.${chrName}.imputed.bcf")
        return tuple(chrName, bcfFile)
    }

process indexBcf {
    tag "${chr}"
    executor 'slurm'
    queue 'gr10478b'
    time '12h'
    
    publishDir "${params.OutDir}/00.raw", mode: 'copy'
    
    input:
    set val(chr), path(bcf) from imputed_array_ch
    
    output:
    set val(chr), path(bcf), path("${bcf}.csi") into indexed_array_ch
    
    script:
    """
    ${params.Bcftools} index -f -c --threads 4 ${bcf}
    """
}

process normalizeBcf {
    tag "${chr}"
    executor 'slurm'
    queue 'gr10478b'
    time '12h'
    
    publishDir "${params.OutDir}/01.normalized", mode: 'copy'
    
    input:
    set val(chr), path(bcf), path(csi) from indexed_array_ch
    
    output:
    set val(chr), path("${chr}.normalized.vcf.gz"), path("${chr}.normalized.vcf.gz.tbi") into normalized_ch, normalized_ch_2
    
    script:
    """
    ${params.Bcftools} view \
        -S ${params.SampleSelectList} \
        --force-samples \
        --threads 8 \
        ${bcf} \
    | ${params.Bcftools} norm \
        --multiallelics -any \
        --fasta-ref ${params.NagasakiPipelinePath}/data/hs38DH.fa \
        --check-ref s \
        --threads 8 \
        -Ou \
    | ${params.Bcftools} annotate \
        --set-id '%CHROM:%POS:%REF:%ALT' \
        --threads 8 \
        -Oz \
        -o ${chr}.normalized.vcf.gz
    
    ${params.Bcftools} index -t --threads 8 ${chr}.normalized.vcf.gz
    """
}

process convertToPlink {
    tag "${chr}"
    executor 'slurm'
    queue 'gr10478b'
    time '12h'
    
    publishDir "${params.OutDir}/02.plink", mode: 'copy'
    
    input:
    set val(chr), path(vcf), path(tbi) from normalized_ch
    
    output:
    set val(chr), path("${chr}.bed"), path("${chr}.bim"), path("${chr}.fam") into plink_ch
    
    script:
    """
    # Convert VCF to PLINK format
    ${params.Plink2} \
        --vcf ${vcf} \
        --make-bed \
        --double-id \
        --out ${chr} \
        --threads 8
    
    # Update sex and pheno information using Python
    python3 << 'EOF'
import pandas as pd
import sys

# Read sample info from Excel
sample_info = pd.read_excel("${params.SampleInfo}")

# Read FAM file
fam = pd.read_csv("${chr}.fam", sep='\\s+', header=None, 
                  names=['FID', 'IID', 'PID', 'MID', 'SEX', 'PHENO'])

# Create mapping dictionaries
sex_map = dict(zip(sample_info['${params.SampleIDColumn}'], sample_info['${params.SampleSexColumn}']))
outcome_map = dict(zip(sample_info['${params.SampleIDColumn}'], sample_info['${params.OutcomeColumn}']))

# Update SEX (1=male, 2=female, 0=unknown)
def convert_sex(iid):
    sex = sex_map.get(iid, 0)
    if pd.isna(sex):
        return 0
    sex_str = str(sex).upper()
    if sex_str in ['M', 'MALE', '1']:
        return 1
    elif sex_str in ['F', 'FEMALE', '2']:
        return 2
    else:
        return 0

# Update PHENO (2=case, 1=control, 0=missing)
def convert_pheno(iid):
    outcome = outcome_map.get(iid, 0)
    if pd.isna(outcome):
        return 0
    if str(outcome) == '${params.CaseValue}':
        return 2
    else:
        return 1

fam['SEX'] = fam['IID'].apply(convert_sex)
fam['PHENO'] = fam['IID'].apply(convert_pheno)

# Write updated FAM file
fam.to_csv("${chr}.fam", sep='\\t', header=False, index=False)
EOF
    """
}

plink_ch
    .toList()
    .map { entries ->
        def chr_order = (1..22).collect { "chr${it}" }
        def filtered = entries.findAll { it[0] in chr_order }
        def sorted = filtered.sort { a, b ->
            Integer.parseInt(a[0].replace("chr", "")) <=> Integer.parseInt(b[0].replace("chr", ""))
        }
        def paths = sorted.collect { it[1].toString().replaceAll(/\.bed$/, '') }
        return paths.join('\n')
    }
    .set { pmerge_list_content }

process writePmergeList {
    executor 'slurm'
    queue 'gr10478b'
    time '1h'
    publishDir "${params.OutDir}/02.plink", mode: 'copy'

    input:
    val(list_text) from pmerge_list_content

    output:
    file("pmerge_lst.txt") into pmerge_lst_ch

    script:
    """
    echo "${list_text}" > pmerge_lst.txt
    """
}

process pmerge {
    executor 'slurm'
    queue 'gr10478b'
    time '36h'
    tag "pmerge"

    publishDir "${params.OutDir}/02.plink", mode: 'symlink'

    input:
    file(pmerge_list) from pmerge_lst_ch

    output:
    tuple file("*.bed"), file("*.bim"), file("*.fam") into pmerge_out

    script:
    out_prefix = "cteph_agp3k.imputed_array"
    """
    export PATH=/home/b/b37974/:$PATH
    plink2 \
        --pmerge-list ${pmerge_list} bfile \
        --threads 8 \
        --make-bed \
        --out ${out_prefix}
    """
}

// Variant QC 步骤1: 计算变异缺失率并过滤
process VariantQC_MissingRate {
    executor 'slurm'
    queue 'gr10478b'
    time '12h'
    tag "Variant missing rate QC"
    publishDir "${params.OutDir}/03.variant_qc_missing", mode: 'symlink'
    
    input:
    tuple path(bed), path(bim), path(fam) from pmerge_out

    output:
    tuple val("autosomes"), path("*.vmiss_qc.bed"), path("*.vmiss_qc.bim"), path("*.vmiss_qc.fam") into vmiss_qc_out
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
    ${params.Plink2} \
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
    ${params.Plink2} \
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


// Variant QC 步骤2: 计算AAF, MAF和HWE (仅针对常染色体)
process VariantQC_CalculateMetrics {
    executor 'slurm'
    queue 'gr10478b'
    time '24h'
    tag "Variant AAF/MAF/HWE calculation for autosomes"
    publishDir "${params.OutDir}/04.variant_stats", mode: 'symlink'
    
    input:
    tuple val(region), path(bed), path(bim), path(fam) from vmiss_qc_out
    
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
    ${params.Plink2} \\
        --bfile ${prefix} \\
        --freq \\
        --out ${prefix}.all \\
        --threads 12
    
    echo "  AAF (ALL) calculation completed" >> variant_qc_calculation.log
    echo "" >> variant_qc_calculation.log
    
    # ===== 步骤2: 计算AAF (CASE) =====
    echo "Step 2: Calculating ALT allele frequency (AAF) for CASE samples..." >> variant_qc_calculation.log
    ${params.Plink2} \\
        --bfile ${prefix} \\
        --keep <(awk '\$6==2 {print \$1, \$2}' ${prefix}.fam) \\
        --freq \\
        --out ${prefix}.case \\
        --threads 12
    
    echo "  AAF (CASE) calculation completed" >> variant_qc_calculation.log
    echo "" >> variant_qc_calculation.log
    
    # ===== 步骤3: 计算AAF (CTRL) =====
    echo "Step 3: Calculating ALT allele frequency (AAF) for CTRL samples..." >> variant_qc_calculation.log
    ${params.Plink2} \\
        --bfile ${prefix} \\
        --keep <(awk '\$6==1 {print \$1, \$2}' ${prefix}.fam) \\
        --freq \\
        --out ${prefix}.ctrl \\
        --threads 12
    
    echo "  AAF (CTRL) calculation completed" >> variant_qc_calculation.log
    echo "" >> variant_qc_calculation.log
    
    # ===== 步骤4: 计算HWE (ALL) =====
    echo "Step 4: Calculating Hardy-Weinberg equilibrium for ALL samples..." >> variant_qc_calculation.log
    ${params.Plink2} \\
        --bfile ${prefix} \\
        --hardy \\
        --out ${prefix}.all \\
        --threads 12
    
    echo "  HWE (ALL) calculation completed" >> variant_qc_calculation.log
    echo "" >> variant_qc_calculation.log
    
    # ===== 步骤5: 计算HWE (CASE) =====
    echo "Step 5: Calculating Hardy-Weinberg equilibrium for CASE samples..." >> variant_qc_calculation.log
    ${params.Plink2} \\
        --bfile ${prefix} \\
        --keep <(awk '\$6==2 {print \$1, \$2}' ${prefix}.fam) \\
        --hardy \\
        --out ${prefix}.case \\
        --threads 12
    
    echo "  HWE (CASE) calculation completed" >> variant_qc_calculation.log
    echo "" >> variant_qc_calculation.log
    
    # ===== 步骤6: 计算HWE (CTRL) =====
    echo "Step 6: Calculating Hardy-Weinberg equilibrium for CTRL samples..." >> variant_qc_calculation.log
    ${params.Plink2} \\
        --bfile ${prefix} \\
        --keep <(awk '\$6==1 {print \$1, \$2}' ${prefix}.fam) \\
        --hardy \\
        --out ${prefix}.ctrl \\
        --threads 12
    
    echo "  HWE (CTRL) calculation completed" >> variant_qc_calculation.log
    echo "" >> variant_qc_calculation.log
    
    # ===== 步骤7: 合并所有统计结果 (内存优化版本) =====
    echo "Step 7: Merging all statistics into summary file (memory-optimized)..." >> variant_qc_calculation.log
    
    python3 << 'EOF'
import pandas as pd
import sys
import gc
import os

print("Reading PLINK2 output files with memory optimization...")

try:
    # 设置内存优化参数
    chunk_size = 50000  # 每次读取的行数
    
    # 检查文件大小
    def check_file_size(filename):
        if os.path.exists(filename):
            size_mb = os.path.getsize(filename) / (1024 * 1024)
            print(f"  File {filename}: {size_mb:.1f} MB")
            return size_mb
        return 0
    
    print("\\nFile sizes:")
    total_size = 0
    for suffix in ['all.afreq', 'case.afreq', 'ctrl.afreq', 'all.hardy', 'case.hardy', 'ctrl.hardy']:
        filename = "${prefix}." + suffix
        size = check_file_size(filename)
        total_size += size
    print(f"  Total input files: {total_size:.1f} MB")
    
    # 使用更小的chunk_size处理大文件
    if total_size > 1000:  # 如果总文件大小超过1GB
        chunk_size = 20000
        print(f"  Large files detected, using smaller chunk size: {chunk_size}")
    
    # 分批读取和合并文件
    print("\\nReading AAF files...")
    
    # 读取AAF文件
    afreq_all = pd.read_csv("${prefix}.all.afreq", sep="\\t", dtype={'ID': str, 'ALT_FREQS': 'float32'})
    print(f"  AAF ALL: {len(afreq_all)} variants")
    
    afreq_case = pd.read_csv("${prefix}.case.afreq", sep="\\t", dtype={'ID': str, 'ALT_FREQS': 'float32'})
    print(f"  AAF CASE: {len(afreq_case)} variants")
    
    afreq_ctrl = pd.read_csv("${prefix}.ctrl.afreq", sep="\\t", dtype={'ID': str, 'ALT_FREQS': 'float32'})
    print(f"  AAF CTRL: {len(afreq_ctrl)} variants")
    
    # 提取需要的列
    afreq_all_sub = afreq_all[["ID", "ALT_FREQS"]].rename(columns={"ALT_FREQS": "AAF_ALL"})
    afreq_case_sub = afreq_case[["ID", "ALT_FREQS"]].rename(columns={"ALT_FREQS": "AAF_CASE"})
    afreq_ctrl_sub = afreq_ctrl[["ID", "ALT_FREQS"]].rename(columns={"ALT_FREQS": "AAF_CTRL"})
    
    # 清理内存
    del afreq_all, afreq_case, afreq_ctrl
    gc.collect()
    
    print("\\nReading HWE files...")
    
    # 读取HWE文件
    hardy_all = pd.read_csv("${prefix}.all.hardy", sep="\\t", dtype={'ID': str, 'P': 'float32'})
    print(f"  HWE ALL: {len(hardy_all)} variants")
    
    hardy_case = pd.read_csv("${prefix}.case.hardy", sep="\\t", dtype={'ID': str, 'P': 'float32'})
    print(f"  HWE CASE: {len(hardy_case)} variants")
    
    hardy_ctrl = pd.read_csv("${prefix}.ctrl.hardy", sep="\\t", dtype={'ID': str, 'P': 'float32'})
    print(f"  HWE CTRL: {len(hardy_ctrl)} variants")
    
    # 提取需要的列
    hardy_all_sub = hardy_all[["ID", "P"]].rename(columns={"P": "HWE_ALL"})
    hardy_case_sub = hardy_case[["ID", "P"]].rename(columns={"P": "HWE_CASE"})
    hardy_ctrl_sub = hardy_ctrl[["ID", "P"]].rename(columns={"P": "HWE_CTRL"})
    
    # 清理内存
    del hardy_all, hardy_case, hardy_ctrl
    gc.collect()
    
    print("\\nMerging data...")
    
    # 逐步合并数据
    merged = afreq_all_sub.merge(afreq_case_sub, on="ID", how="outer")
    del afreq_all_sub, afreq_case_sub
    gc.collect()
    
    merged = merged.merge(afreq_ctrl_sub, on="ID", how="outer")
    del afreq_ctrl_sub
    gc.collect()
    
    # 计算MAF
    print("\\nCalculating MAF...")
    merged["MAF_ALL"] = merged["AAF_ALL"].apply(lambda x: min(x, 1-x) if pd.notna(x) else None)
    merged["MAF_CASE"] = merged["AAF_CASE"].apply(lambda x: min(x, 1-x) if pd.notna(x) else None)
    merged["MAF_CTRL"] = merged["AAF_CTRL"].apply(lambda x: min(x, 1-x) if pd.notna(x) else None)
    
    # 合并HWE数据
    print("\\nMerging HWE data...")
    merged = merged.merge(hardy_all_sub, on="ID", how="outer")
    del hardy_all_sub
    gc.collect()
    
    merged = merged.merge(hardy_case_sub, on="ID", how="outer")
    del hardy_case_sub
    gc.collect()
    
    merged = merged.merge(hardy_ctrl_sub, on="ID", how="outer")
    del hardy_ctrl_sub
    gc.collect()
    
    # 重新排列列顺序
    merged = merged[["ID", "AAF_ALL", "AAF_CASE", "AAF_CTRL", 
                     "MAF_ALL", "MAF_CASE", "MAF_CTRL",
                     "HWE_ALL", "HWE_CASE", "HWE_CTRL"]]
    
    # 重命名ID列为VARIANT_ID
    merged = merged.rename(columns={"ID": "VARIANT_ID"})
    
    print("\\nParsing VARIANT_ID for sorting...")
    
    # 更内存友好的排序方法
    def get_sort_key(variant_id):
        try:
            parts = variant_id.split(":")
            chrom = parts[0].replace("chr", "").replace("Chr", "").replace("CHR", "")
            pos = int(parts[1])
            
            if chrom.isdigit():
                chr_num = int(chrom)
            elif chrom.upper() == "X":
                chr_num = 23
            elif chrom.upper() == "Y":
                chr_num = 24
            elif chrom.upper() in ["M", "MT"]:
                chr_num = 25
            else:
                chr_num = 99
                
            return (chr_num, pos)
        except:
            return (99, 999999999)
    
    merged["sort_key"] = merged["VARIANT_ID"].apply(get_sort_key)
    merged = merged.sort_values("sort_key")
    merged = merged.drop("sort_key", axis=1)
    
    print(f"  Sorted by chromosome and position")
    
    # 保存为TSV文件
    output_file = "${prefix}.variant_qc_sum.tsv"
    merged.to_csv(output_file, sep="\\t", index=False, na_rep="NA")
    
    print(f"  Summary file created: {output_file}")
    print(f"  Total variants in summary: {len(merged)}")
    
    # 统计一些基本信息
    if len(merged) > 0:
        print("\\nStatistics summary:")
        print(f"  AAF_ALL range: [{merged['AAF_ALL'].min():.4f}, {merged['AAF_ALL'].max():.4f}]")
        print(f"  MAF_ALL range: [{merged['MAF_ALL'].min():.4f}, {merged['MAF_ALL'].max():.4f}]")
        print(f"  HWE_ALL range: [{merged['HWE_ALL'].min():.4e}, {merged['HWE_ALL'].max():.4e}]")
        
        # 统计MAF和AAF不同的变异数
        maf_aaf_diff_all = (merged['MAF_ALL'] != merged['AAF_ALL']).sum()
        maf_aaf_diff_case = (merged['MAF_CASE'] != merged['AAF_CASE']).sum()
        maf_aaf_diff_ctrl = (merged['MAF_CTRL'] != merged['AAF_CTRL']).sum()
        
        print(f"\\nVariants where MAF != AAF (ALT allele is major allele, AAF > 0.5):")
        print(f"  ALL: {maf_aaf_diff_all} variants ({maf_aaf_diff_all/len(merged)*100:.2f}%)")
        print(f"  CASE: {maf_aaf_diff_case} variants ({maf_aaf_diff_case/len(merged)*100:.2f}%)")
        print(f"  CTRL: {maf_aaf_diff_ctrl} variants ({maf_aaf_diff_ctrl/len(merged)*100:.2f}%)")
        
        # 统计HWE显著偏离的变异数
        hwe_fail_all = (merged['HWE_ALL'] < 1e-6).sum()
        hwe_fail_case = (merged['HWE_CASE'] < 1e-6).sum()
        hwe_fail_ctrl = (merged['HWE_CTRL'] < 1e-6).sum()
        
        print(f"\\nHWE violations (P < 1e-6):")
        print(f"  ALL: {hwe_fail_all} variants ({hwe_fail_all/len(merged)*100:.2f}%)")
        print(f"  CASE: {hwe_fail_case} variants ({hwe_fail_case/len(merged)*100:.2f}%)")
        print(f"  CTRL: {hwe_fail_ctrl} variants ({hwe_fail_ctrl/len(merged)*100:.2f}%)")
    
    print("\\nMemory cleanup completed")
    
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

// Variant QC 步骤3: 按染色体拆分variant stats用于并行注释
process VariantQC_SplitStatsByChr {
    tag "Split stats for ${region}"
    
    input:
    tuple val(region), path(variant_stats) from variant_qc_stats
    
    output:
    tuple val(region), path("*_stats.tsv") into variant_qc_stats_by_chr mode flatten
    
    script:
    """
    python3 << 'EOF'
import pandas as pd
import gc

print("Splitting variant stats by chromosome with memory optimization...")

# 使用分块读取大文件
chunk_size = 100000  # 每次读取10万行
chr_files = {}
header = None

print("Reading input file in chunks...")
chunk_count = 0

for chunk in pd.read_csv("${variant_stats}", sep="\\t", 
                          dtype={'VARIANT_ID': str}, 
                          chunksize=chunk_size):
    chunk_count += 1
    
    # 保存表头
    if header is None:
        header = chunk.columns.tolist()
    
    # 提取染色体信息
    chunk['CHR'] = chunk['VARIANT_ID'].str.split(':', expand=True)[0]
    
    # 按染色体分组
    for chr_name, group in chunk.groupby('CHR', sort=False):
        # 初始化文件
        if chr_name not in chr_files:
            output_file = f"{chr_name}_stats.tsv"
            chr_files[chr_name] = output_file
            # 写入表头
            group.drop('CHR', axis=1).to_csv(output_file, sep='\\t', index=False, mode='w')
        else:
            # 追加数据（不写表头）
            output_file = chr_files[chr_name]
            group.drop('CHR', axis=1).to_csv(output_file, sep='\\t', index=False, mode='a', header=False)
    
    # 清理内存
    del chunk
    gc.collect()
    
    if chunk_count % 10 == 0:
        print(f"  Processed {chunk_count * chunk_size:,} variants...")

print(f"\\nSplit completed:")
for chr_name, output_file in chr_files.items():
    line_count = sum(1 for _ in open(output_file)) - 1  # 减去表头
    print(f"  {chr_name}: {line_count:,} variants -> {output_file}")

print(f"\\nTotal chromosomes: {len(chr_files)}")

EOF
    """
}

// 将输出转换为正确的channel格式，从文件名提取染色体信息
// 提取 [chr, tsv] 格式用于join
variant_qc_stats_by_chr
    .flatMap { region, files ->
        def fileList = files instanceof List ? files : [files]
        fileList.collect { file ->
            def chr = file.baseName.replaceAll(/_stats$/, '')
            [chr, file]
        }
    }
    .set { variant_qc_stats_by_chr_only }

// normalized_ch_2 格式: [chr, vcf.gz, vcf.gz.tbi]
// variant_qc_stats_by_chr_only 格式: [chr, tsv]
// join后格式: [chr, vcf.gz, vcf.gz.tbi, tsv]
normalized_ch_2
    .join(variant_qc_stats_by_chr_only)
    .set { chr_vcf_stats_joined }

// chr_vcf_stats_joined.view()

// Variant QC 步骤3a: 按染色体并行添加Marker类型和R2信息
process VariantQC_AddMarkerInfo_ByChr {
    executor 'slurm'
    queue 'gr10478b'
    time '24h'
    tag "Add Marker info for ${chr}"
    
    input:
    tuple val(chr), path(vcf_file), path(vcf_idx), path(chr_stats) from chr_vcf_stats_joined
    
    output:
    tuple val(chr), file("${chr}.annotated.tsv") into variant_qc_annotated_by_chr_only
    
    when:
    true
    
    script:
    """
    python3 << 'EOF'
import pandas as pd
import subprocess
import sys
import gc

print("Processing chromosome: ${chr}")

# 使用分块读取避免内存溢出
chunk_size = 50000  # 每次处理5万个变体
all_results = []

print(f"  Reading input file in chunks (chunk_size={chunk_size})...")

chunk_num = 0
for chunk in pd.read_csv("${chr_stats}", sep="\\t", dtype={'VARIANT_ID': str}, chunksize=chunk_size):
    chunk_num += 1
    print(f"\\n  Processing chunk {chunk_num} ({len(chunk)} variants)...")
    
    # 提取位置信息
    chunk['POS'] = chunk['VARIANT_ID'].str.split(':', expand=True)[1].astype(int)
    chunk['REF'] = chunk['VARIANT_ID'].str.split(':', expand=True)[2]
    chunk['ALT'] = chunk['VARIANT_ID'].str.split(':', expand=True)[3]
    
    # 初始化结果列
    chunk['Marker_Type'] = 'NA'
    chunk['IMPUTED_R2'] = float('nan')
    
    # 批量查询：获取该chunk覆盖的位置范围
    min_pos = chunk['POS'].min()
    max_pos = chunk['POS'].max()
    
    print(f"    Batch querying VCF for region ${chr}:{min_pos}-{max_pos}...")
    
    # 使用 chr:pos-pos 格式批量查询整个范围
    tabix_cmd = ['tabix', '${vcf_file}', f"${chr}:{min_pos}-{max_pos}"]
    
    try:
        result = subprocess.run(tabix_cmd, capture_output=True, text=True, check=False)
        
        if result.returncode != 0:
            print(f"    WARNING: tabix returned error code {result.returncode}")
            print(f"    stderr: {result.stderr}")
        
        # 解析VCF输出，建立位置->VCF行的索引
        print(f"    Parsing VCF output and building position index...")
        vcf_index = {}  # key: (pos, ref, alt), value: (marker_type, r2_value)
        
        if result.stdout.strip():
            for vcf_line in result.stdout.strip().split('\\n'):
                fields = vcf_line.split('\\t')
                if len(fields) < 8:
                    continue
                
                vcf_pos = int(fields[1])
                vcf_ref = fields[3]
                vcf_alt = fields[4]
                info_field = fields[7]
                
                # 检查 TYPED 和 IMPUTED flags
                has_typed = 'TYPED' in info_field
                has_imputed = 'IMPUTED' in info_field
                
                if has_typed and has_imputed:
                    marker_type = 'TYPED-IMPUTED'
                elif has_typed:
                    marker_type = 'TYPED'
                elif has_imputed:
                    marker_type = 'IMPUTED'
                else:
                    marker_type = 'UNKNOWN'
                
                # 提取 R2 值（精确匹配 R2=，不匹配 ER2=）
                r2_value = float('nan')
                for field in info_field.split(';'):
                    if field.startswith('R2='):
                        try:
                            r2_value = float(field.split('=', 1)[1])
                        except (ValueError, IndexError):
                            pass
                        break
                
                # 存入索引（key为位置+REF+ALT的组合）
                key = (vcf_pos, vcf_ref, vcf_alt)
                vcf_index[key] = (marker_type, r2_value)
        
        print(f"    VCF index built: {len(vcf_index)} unique variants")
        
        # 使用索引快速匹配变体
        print(f"    Annotating variants...")
        annotated_count = 0
        not_found_count = 0
        
        for idx, row in chunk.iterrows():
            pos = row['POS']
            ref = row['REF']
            alt = row['ALT']
            
            key = (pos, ref, alt)
            if key in vcf_index:
                marker_type, r2_value = vcf_index[key]
                chunk.at[idx, 'Marker_Type'] = marker_type
                chunk.at[idx, 'IMPUTED_R2'] = r2_value
                annotated_count += 1
            else:
                not_found_count += 1
        
        print(f"    Chunk {chunk_num} results:")
        print(f"      Annotated: {annotated_count}/{len(chunk)} ({annotated_count/len(chunk)*100:.1f}%)")
        print(f"      Not found: {not_found_count}/{len(chunk)} ({not_found_count/len(chunk)*100:.1f}%)")
        
    except Exception as e:
        print(f"    ERROR processing chunk {chunk_num}: {e}")
        import traceback
        traceback.print_exc()
    
    # 删除临时列
    chunk = chunk.drop(['POS', 'REF', 'ALT'], axis=1)
    
    # 保存该chunk的结果
    all_results.append(chunk)
    
    # 清理内存
    del chunk
    gc.collect()

# 合并所有chunks
print(f"\\n  Merging all chunks...")
final_result = pd.concat(all_results, ignore_index=True)
print(f"  Total variants processed: {len(final_result)}")

# 统计总体结果
total_annotated = (final_result['Marker_Type'] != 'NA').sum()
total_not_found = (final_result['Marker_Type'] == 'NA').sum()

print(f"\\n=== Final Statistics ===")
print(f"  Total annotated: {total_annotated}/{len(final_result)} ({total_annotated/len(final_result)*100:.1f}%)")
print(f"  Total not found: {total_not_found}/{len(final_result)} ({total_not_found/len(final_result)*100:.1f}%)")

# 保存注释后的文件
output_file = "${chr}.annotated.tsv"
final_result.to_csv(output_file, sep='\\t', index=False, na_rep='NA')
print(f"\\n  Output: {output_file}")

EOF
    """
}

// 将所有染色体的注释文件收集到一起，按染色体顺序排序
variant_qc_annotated_by_chr_only
    .toSortedList { a, b -> 
        // 提取染色体编号进行排序
        def chrA = a[0].replaceAll(/chr/, '').toInteger()
        def chrB = b[0].replaceAll(/chr/, '').toInteger()
        chrA <=> chrB
    }
    .map { list -> list.collect { chr, file -> file } }  // 只取文件部分
    .set { all_annotated_files }

// all_annotated_files.view()

// Variant QC 步骤3b: 合并所有染色体的注释结果
process VariantQC_MergeAnnotations {
    executor 'slurm'
    queue 'gr10478b'
    time '24h'
    tag "Merge annotations for autosomes"
    publishDir "${params.OutDir}/05.variant_stats_final", mode: 'symlink'
    
    input:
    file(annotated_files) from all_annotated_files
    
    output:
    file("*.variant_qc_final.tsv") into variant_qc_final
    file "marker_info_merge.log"
    
    script:
    """
    echo "=== Marker Type and R2 Annotation Merge Report ===" > marker_info_merge.log
    echo "Date: \$(date)" >> marker_info_merge.log
    echo "Region: autosomes (chr1-chr22)" >> marker_info_merge.log
    echo "" >> marker_info_merge.log
    
    # 按照chr1-chr22的顺序合并TSV文件
    echo "Merging chromosome files in order (chr1-chr22)..."
    
    output_file="cteph_agp3k.imputed_array.vmiss_qc.variant_qc_final.tsv"
    
    # 使用Python合并文件（文件已按chr1-chr22顺序排列）
    python3 << 'EOF'
import sys

print("Merging annotated chromosome files...")

# 获取所有输入文件，Nextflow已经按chr顺序传入
input_files = "${annotated_files}".split()
print(f"Total input files: {len(input_files)}")

if len(input_files) == 0:
    print("ERROR: No input files found!")
    sys.exit(1)

output_file = "cteph_agp3k.imputed_array.vmiss_qc.variant_qc_final.tsv"
variant_count = 0

with open("marker_info_merge.log", "a") as log:
    log.write("\\nFile processing order:\\n")
    
    # 第一个文件：保留表头
    first_file = input_files[0]
    print(f"  Processing file 1/{len(input_files)} (with header): {first_file}")
    log.write(f"  File 1: {first_file}\\n")
    
    with open(first_file, 'r') as f:
        lines = f.readlines()
        with open(output_file, 'w') as out:
            out.writelines(lines)
        variant_count += len(lines) - 1  # 减去表头
        log.write(f"    Variants: {len(lines) - 1}\\n")
    
    # 其余文件：跳过表头
    for i, input_file in enumerate(input_files[1:], start=2):
        print(f"  Processing file {i}/{len(input_files)} (without header): {input_file}")
        log.write(f"  File {i}: {input_file}\\n")
        
        with open(input_file, 'r') as f:
            lines = f.readlines()[1:]  # 跳过表头
            with open(output_file, 'a') as out:
                out.writelines(lines)
            variant_count += len(lines)
            log.write(f"    Variants: {len(lines)}\\n")
    
    log.write(f"\\nTotal merged variants: {variant_count}\\n")

print(f"\\nMerge completed: {variant_count} variants")
print(f"Output file: {output_file}")

EOF
    
    echo "" >> marker_info_merge.log
    echo "======================================" >> marker_info_merge.log
    echo "SUMMARY" >> marker_info_merge.log
    echo "======================================" >> marker_info_merge.log
    echo "Merge completed at: \$(date)" >> marker_info_merge.log
    echo "Merge method: Sequential concatenation by chromosome order (chr1-chr22)" >> marker_info_merge.log
    echo "  - First file (chr1): included with header" >> marker_info_merge.log
    echo "  - Remaining files (chr2-chr22): appended without headers" >> marker_info_merge.log
    """
}
