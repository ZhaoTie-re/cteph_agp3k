params.ArrayRAWPath = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v5/array/raw_data_ph'
params.SampleSelectList = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v5/info/cteph_agp3k.v5.ls'
params.NagasakiPipelinePath = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/nagasaki_pipeline'
params.ToMMoVCF = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/ToMMo_60KJPN/tommo-60kjpn-20240904-GRCh38-snvindel-af-autosome.norm.vcf.gz'
params.ToMMoIndex = params.ToMMoVCF + '.tbi'
params.ChrRenameFile = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/Info/chr_rename.txt'
params.OutputDir = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v5/array/results'
params.ScriptsDir = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v5/array/scripts'
params.Plink2Path = '/home/b/b37974/plink2'
params.PlinkPath = '/home/b/b37974/plink'
params.TabixPath = '/home/b/b37974/htslib-1.9/tabix'
params.SampleInfo = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v5/info/ph_agp3k_combined_final_20251112.xlsx'
params.SampleIDColumn = 'Genome Scan ID'
params.SampleSexColumn = 'Sex'
params.OutcomeColumn = 'Outcome'
params.CaseValue = 'PH'
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
    ${params.ScriptsDir}/01_sample_select.sh \\
        ${params.Plink2Path} \\
        ${sample_list} \\
        ${prefix}
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
    path "*.sqc.log" into sample_qc_logs

    script:
    def prefix = bed.baseName
    """
    ${params.ScriptsDir}/02_sample_qc.sh \\
        ${params.Plink2Path} \\
        ${prefix}
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
    ${params.ScriptsDir}/02_check_duplicate_iid.sh
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
    ${params.ScriptsDir}/03_convert_to_vcf.sh \\
        ${params.Plink2Path} \\
        ${params.PlinkPath} \\
        ${params.NagasakiPipelinePath} \\
        ${chr_rename} \\
        ${prefix}
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
    ${params.ScriptsDir}/04_vcf_to_plink.sh \\
        ${params.Plink2Path} \\
        ${sample_info} \\
        "${id_col}" \\
        "${sex_col}" \\
        ${vcf} \\
        ${prefix}
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
    path "merge.log" into merge_log

    script:
    """
    ${params.ScriptsDir}/05_merge_plink.sh \\
        ${params.Plink2Path} \\
        ${params.PlinkPath}
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
    path "multiallelic_check.log" into multiallelic_log
    path "multiallelic_variants.txt" optional true

    script:
    def prefix = bed.baseName
    """
    ${params.ScriptsDir}/06_remove_multiallelic.sh \\
        ${params.Plink2Path} \\
        ${prefix}
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
    path "variant_missing_qc.log" into variant_miss_log

    script:
    def prefix = bed.baseName
    """
    ${params.ScriptsDir}/07_variant_qc_missing.sh \\
        ${params.Plink2Path} \\
        ${prefix}
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
    ${params.ScriptsDir}/08_update_pheno_split.sh \\
        ${params.Plink2Path} \\
        ${sample_info} \\
        "${id_col}" \\
        "${outcome_col}" \\
        "${case_val}" \\
        ${prefix}
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
    ${params.ScriptsDir}/09_variant_qc_metrics.sh \\
        ${params.Plink2Path} \\
        "${region}" \\
        ${prefix}
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
    tuple val(region), path("cteph_agp3k.array.sqc.vqc.bed"), path("cteph_agp3k.array.sqc.vqc.bim"), path("cteph_agp3k.array.sqc.vqc.fam") into qc_filtered_plink, qc_filtered_plink_2
    path "*.variants_to_remove.txt"
    path "variant_filter.log" into variant_filter_logs
    
    script:
    def prefix = bed.baseName
    def output_prefix = "cteph_agp3k.array.sqc.vqc"
    def maf_mode = params.MAF_FilterMode ?: "ALL"
    def maf_threshold = params.MAF_Threshold ?: 0.01
    def hwe_mode = params.HWE_FilterMode ?: "CTRL|CASE"
    def hwe_ctrl_threshold = params.HWE_CTRL_Threshold ?: 1e-6
    def hwe_case_threshold = params.HWE_CASE_Threshold ?: 1e-10
    """
    ${params.ScriptsDir}/10_variant_qc_filter.sh \\
        ${params.Plink2Path} \\
        "${region}" \\
        ${prefix} \\
        ${tsv} \\
        "${maf_mode}" \\
        "${maf_threshold}" \\
        "${hwe_mode}" \\
        "${hwe_ctrl_threshold}" \\
        "${hwe_case_threshold}"
    """
}

// 对比TOMMO数据的AF (Autosomes only)
process CompareWithToMMo {
    executor 'slurm'
    queue 'gr10478b'
    time '24h'
    tag "Comparing with ToMMo"
    publishDir "${params.OutputDir}/11.tommo_compare", mode: 'symlink'
    
    input:
    tuple val(region), path(bed), path(bim), path(fam) from qc_filtered_plink
    path(tommo_vcf) from params.ToMMoVCF
    path(tommo_index) from params.ToMMoIndex
    
    output:
    tuple val(region), path("*.variant_qc_with_tommo.tsv.gz"), path("*.variant_qc_with_tommo.tsv.gz.tbi") into tommo_compare_out
    path("*.report.txt")
    path("*.panel_compare.pdf")
    
    script:
    def prefix = bed.baseName
    """
    source activate cteph_geno_pro
    python3 ${params.ScriptsDir}/panel_compare_main.py \\
        --bed_prefix ${prefix} \\
        --output_prefix ${prefix} \\
        --threads 32 \\
        --tommo_vcf_path ${tommo_vcf} \\
        --chunk_size 50000 \\
        --max_workers 4 \\
        --tabix_path ${params.TabixPath} \\
        --grouping_metric MAF_ALL \\
        --grouping_threshold 0.05
    """
}

// 从ToMMo对比结果中提取最终变异列表(TOMMO_FILTER==PASS & 95% Central Range of AAF_CTRL - TOMMO_AAF)
process FilterToMMoCentralRange {
    executor 'slurm'
    queue 'gr10478b'
    time '12h'
    tag "Filtering ToMMo Central Range"
    publishDir "${params.OutputDir}/12.tommo_central_filter", mode: 'symlink'

    input:
    tuple val(region), path(tsv), path(tbi) from tommo_compare_out

    output:
    tuple val(region), path("*.tommo_pass_95pct.ids.txt") into tommo_central_ids
    path "*.tommo_pass_95pct.plot.pdf"
    path "*.tommo_pass_95pct.report.txt" into tommo_filter_reports

    script:
    def prefix = tsv.simpleName.replaceAll(/\.variant_qc_with_tommo$/, '')
    """
    source activate cteph_geno_pro
    python3 ${params.ScriptsDir}/12_filter_tommo_central_range.py \\
        --input_tsv ${tsv} \\
        --output_prefix ${prefix}
    """
}

// Join the channels for the final extraction step
tommo_central_ids
    .join(qc_filtered_plink_2)
    .set { extract_input }

// 根据ToMMo过滤结果提取最终PLINK文件
process ExtractToMMoCentralVariants {
    executor 'slurm'
    queue 'gr10478b'
    time '12h'
    tag "Extracting ToMMo Central Variants"
    publishDir "${params.OutputDir}/13.final_variants", mode: 'symlink'

    input:
    tuple val(region), path(ids), path(bed), path(bim), path(fam) from extract_input

    output:
    tuple val(region), path("*.tommo.bed"), path("*.tommo.bim"), path("*.tommo.fam") into final_plink
    path "*.log"

    script:
    def prefix = bed.baseName
    """
    ${params.Plink2Path} \\
        --bfile ${prefix} \\
        --extract ${ids} \\
        --make-bed \\
        --out ${prefix}.tommo \\
        --threads 8
    """
}

//使用plink2将最终的plink文件转换为vcf格式（vcf.gz）并索引（bcftools index, tbi）
process FinalPLINKtoVCF {
    executor 'slurm'
    queue 'gr10478b'
    time '12h'
    tag "Final PLINK to VCF conversion"
    publishDir "${params.OutputDir}/14.final_vcf", mode: 'symlink'

    input:
    tuple val(region), path(bed), path(bim), path(fam) from final_plink
    path(chr_rename) from params.ChrRenameFile

    output:
    tuple val(region), path("*.tommo.vcf.gz"), path("*.tommo.vcf.gz.tbi") into final_vcf_out
    path "*.log" into final_vcf_logs

    script:
    def prefix = bed.baseName
    """
    ${params.Plink2Path} \\
        --bfile ${prefix} \\
        --recode vcf bgz id-paste=iid \\
        --out ${prefix}.temp \\
        --threads 8

    bcftools annotate \\
        --rename-chrs ${chr_rename} \\
        -O z \\
        -o ${prefix}.tommo.vcf.gz \\
        --threads 8 \\
        ${prefix}.temp.vcf.gz

    bcftools index \\
        --tbi \\
        --threads 8 \\
        ${prefix}.tommo.vcf.gz
    """
}


