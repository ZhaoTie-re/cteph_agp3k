params.gtPath = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/wgs'
params.snpeffDir = '/LARGE1/gr10478/platform/JHRPv4/workspace/pipeline/output/snpEff.v4.index'
params.scriptDir = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/analysis/assoc_rvtest.rev1/scripts'
params.outDir = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/analysis/assoc_rvtest.rev1/results'
params.performNorm = false
params.ac_threshold = 2
params.num_var_threshold = 3
params.plink2 = '/home/b/b37974/plink2_alpha6/plink2'
params.remove_samples = "${params.scriptDir}/sample_rm.ls"

process FilterMAC {
    executor 'slurm'
    queue 'gr10478b'
    time '12h'

    publishDir "${params.outDir}/00.pre_step", mode: 'symlink'

    input:
    val gtPath from params.gtPath

    output:
    tuple file("*.bed"), file("*.bim"), file("*.fam") into mac_filtered_ch
    file('*.log')

    script:
    input_bfile = "${gtPath}/18.tommo_panel_filter/cteph_agp3k.rare"
    out_bfile = "cteph_agp3k.rare.mac${params.ac_threshold}"
    """
    ${params.plink2} --bfile ${input_bfile} \
        --mac ${params.ac_threshold} \
        --make-bed \
        --out ${out_bfile} \
        --threads 32
    """
}

process RemoveSamples {
    executor 'slurm'
    queue 'gr10478b'
    time '12h'

    publishDir "${params.outDir}/00.pre_step", mode: 'symlink'

    input:
    tuple file(bed), file(bim), file(fam) from mac_filtered_ch

    output:
    tuple file("*.bed"), file("*.bim"), file("*.fam") into ready_plink_ch
    file('*.log')

    script:
    base = bed.baseName
    out_bfile = "${base}.rm_samples"
    """
    # Generate FID IID (same as IID) for plink remove
    awk '{print \$1, \$1}' ${params.remove_samples} > sample_rm.fid_iid.txt

    ${params.plink2} --bfile ${base} \
        --remove sample_rm.fid_iid.txt \
        --make-bed \
        --out ${out_bfile} \
        --threads 32
    """
}

process RVtestPrepare {
    executor 'slurm'
    queue 'gr10478b'
    time '36h'

    publishDir "${params.outDir}/01.rvtest_prepare", mode: 'symlink'

    input:
    tuple file(bed), file(bim), file(fam) from ready_plink_ch
    val gtPath from params.gtPath

    output:
    file('*.log')
    tuple file("*.vcf.gz"), file("*.vcf.gz.tbi") into rvtest_vcf_ch
    tuple file("*.pheno_df.csv"), file("*.cov_df.no_age.csv") into (rvtest_pheno_covar_ch, sens_pheno_covar_ch)
    file("refFlat.hg38.nochr.txt.gz") into (rvtest_refflat_ch, sens_refflat_ch)

    script:
    bed_prefix = bed.baseName
    pheno = "${gtPath}/19.cov_pheno_prepare/cteph_agp3k.bbj.projection.pheno_df.csv"
    covar = "${gtPath}/19.cov_pheno_prepare/cteph_agp3k.bbj.projection.cov_df.no_age.csv"
    refflat = "/LARGE0/gr10478/b37974/Pulmonary_Hypertension/anno_raw/refFlat.hg38.txt.gz"
    norm_flag = params.performNorm ? "--norm" : ""
    """
    source activate cteph_geno_pro
    python ${params.scriptDir}/rvtest_prepare_main.py \
        --bed-prefix ${bed_prefix} \
        --pheno-path ${pheno} \
        --covar-path ${covar} \
        --refflat-path ${refflat} \
        --threads 32 \
        --verbose \
        --log-file rvtest_prepare.log \
        ${norm_flag}
    """
}

process SnpEffAnnotate {
    executor 'slurm'
    queue 'gr10478b'
    time '48h'

    publishDir "${params.outDir}/02.snpeff_annotate", mode: 'symlink'

    input:
    val snpeffDir from params.snpeffDir
    tuple file(vcf), file(vcf_tbi) from rvtest_vcf_ch

    output:
    file('*.log')
    file('*.tsv') into snpeff_annotated_stats
    tuple file("*.snpeff.vcf.gz"), file("*.snpeff.vcf.gz.tbi") into snpeff_annotated_vcf_ch

    script:
    """
    source activate cteph_geno_pro
    python ${params.scriptDir}/snpeff_anno_main.py \
        --vcf-path ${vcf} \
        --snpeff-dir ${snpeffDir} \
        --parallel \
        --max-workers 16 \
        --threads 16 \
        --keep-cache \
    """
}

process PlotSnpEffStats {
    // executor 'slurm'
    // queue 'gr10478b'
    time '1h'

    publishDir "${params.outDir}/02.snpeff_annotate", mode: 'symlink'

    input:
    file stats_file from snpeff_annotated_stats

    output:
    file("*.pdf")

    script:
    """
    source activate cteph_geno_pro
    python ${params.scriptDir}/plot_snpeff_stats.py \
        --input ${stats_file} \
        --output snpeff_impact_effect_dist
    """
}

impact_filters = Channel.from([
    ["impact_low_moderate_high", "LOW MODERATE HIGH"],
    ["impact_moderate_high", "MODERATE HIGH"],
    ["impact_high", "HIGH"]
])

process InfoFilter {
    // executor 'slurm'
    // queue 'gr10478b'
    time '48h'
    tag "${filter_tag}"

    publishDir "${params.outDir}/03.info_filter", mode: 'symlink'

    input:
    tuple val(filter_tag), val(filter_values), file(vcf), file(vcf_tbi) from impact_filters.combine(snpeff_annotated_vcf_ch)

    output:
    // file('*.log')
    file('*.json')
    tuple val(filter_tag), file("${out_prefix}.vcf.gz"), file("${out_prefix}.vcf.gz.tbi") into (info_filtered_vcf_ch, sensitivity_vcf_ch)

    script:
    out_prefix="cteph_agp3k.rare.${filter_tag}"
    """
    source activate cteph_geno_pro
    python ${params.scriptDir}/info_filter_main.py \
        --input ${vcf} \
        --info-key impact \
        --values ${filter_values} \
        --out-prefix ${out_prefix} \
        --threads 16 \
        --check-chr-prefix \
        --no-keep-chr-prefix \
        --verbose
    """
}

test_methods = Channel.from([
    ["skato", "--kernel skato"],
    ["burden", "--burden cmc"]
])

rvtest_run_inputs = test_methods
    .combine(info_filtered_vcf_ch)
    .combine(rvtest_pheno_covar_ch)
    .combine(rvtest_refflat_ch)

process RVtestRun {
    // executor 'slurm'
    // queue 'gr10478b'
    time '48h'
    tag "${filter_tag}_${method_tag}"

    publishDir "${params.outDir}/04.rvtest_run/${filter_tag}/${method_tag}", mode: 'symlink'

    input:
    tuple val(method_tag), val(method_opt), 
          val(filter_tag), file(vcf), file(vcf_tbi), 
          file(pheno), file(covar), 
          file(refflat) from rvtest_run_inputs

    output:
    file("*.log")
    tuple val(filter_tag), val(method_tag), file("*.assoc") into rvtest_result_ch

    script:
    out_prefix="cteph_agp3k.rare.${filter_tag}.${method_tag}"
    """
    export PATH=/home/b/b37974/rvtests/executable/:$PATH
    rvtest \
        --inVcf ${vcf} \
        --pheno ${pheno} --pheno-name pheno1 \
        --covar ${covar} --covar-name sex,pc1_avg,pc2_avg,pc3_avg,pc4_avg,pc5_avg,pc6_avg,pc7_avg,pc8_avg,pc9_avg,pc10_avg \
        --geneFile ${refflat} \
        --out ${out_prefix} \
        --noweb \
        --numThread 4 \
        ${method_opt}
    """
}

process RVtestPostProcess {
    // executor 'slurm'
    // queue 'gr10478b'
    time '1h'
    tag "${filter_tag}_${method_tag}"

    publishDir "${params.outDir}/05.post_process/${filter_tag}/${method_tag}", mode: 'symlink'

    input:
    tuple val(filter_tag), val(method_tag), file(assoc_file) from rvtest_result_ch

    output:
    tuple val(filter_tag), val(method_tag), file("*.fdr.assoc") into rvtest_final_ch

    script:
    out_file = "${assoc_file.baseName}.filtered.fdr.assoc"
    """
    source activate cteph_geno_pro
    python ${params.scriptDir}/rvtest_post_process.py \
        --input ${assoc_file} \
        --output ${out_file} \
        --num-var-threshold ${params.num_var_threshold}
    """
}

process RVtestVisualization {
    // executor 'slurm'
    // queue 'gr10478b'
    time '1h'
    tag "${filter_tag}_${method_tag}"

    publishDir "${params.outDir}/06.visualization/${filter_tag}/${method_tag}", mode: 'symlink'

    input:
    tuple val(filter_tag), val(method_tag), file(assoc_file) from rvtest_final_ch

    output:
    tuple val(filter_tag), val(method_tag), file("*.png"), file("*.pdf") into rvtest_plots_ch

    script:
    out_prefix = "${assoc_file.baseName}"
    title = "RVTest: ${filter_tag} - ${method_tag}"
    """
    source activate cteph_geno_pro
    python ${params.scriptDir}/plot_manhattan_qq.py \
        --input ${assoc_file} \
        --output-prefix ${out_prefix} \
        --title "${title}"
    """
}

// Sensitivity Check Data Preparation
summary_file_path = "${params.gtPath}/18.tommo_panel_filter/cteph_agp3k.variant_qc_summary.variant_qc_with_tommo.tsv.summary_filter.tsv"
sensitivity_summary_ch = Channel.value(file(summary_file_path))

process SensitivityDataPrepare {
    // executor 'slurm'
    // queue 'gr10478b'
    tag "${filter_tag}"
    publishDir "${params.outDir}/07.sensitivity_check/00.data_prepare", mode: 'symlink'
    time '12h' 

    input:
    tuple val(filter_tag), file(vcf), file(tbi) from sensitivity_vcf_ch
    file summary_file from sensitivity_summary_ch

    output:
    tuple val(filter_tag), file("*.stat1.vcf.gz"), file("*.stat1.vcf.gz.tbi") into sensitivity_stat1_ch
    tuple val(filter_tag), file("*.stat1_stat2.vcf.gz"), file("*.stat1_stat2.vcf.gz.tbi") into sensitivity_stat12_ch
    file "*.json"

    script:
    out_prefix = vcf.name.replace(".vcf.gz", "")
    """
    bash ${params.scriptDir}/run_sensitivity_prepare.sh \
        ${vcf} \
        ${summary_file} \
        ${out_prefix} \
        ${task.cpus} \
        ${params.scriptDir}
    """
}

// -----------------------------------------------------------
// Sensitivity Analysis Steps
// -----------------------------------------------------------

sensitivity_test_methods = Channel.from([
    ["skato", "--kernel skato"],
    ["burden", "--burden cmc"]
])

// Merge Stat1 and Stat12 channels for processing
// Update filter_tag to include sensitivity subset name
sensitivity_combined_vcf_ch = sensitivity_stat1_ch
    .map { tag, vcf, tbi -> [ "${tag}.stat1", vcf, tbi ] }
    .mix( sensitivity_stat12_ch.map { tag, vcf, tbi -> [ "${tag}.stat1_stat2", vcf, tbi ] } )

sensitivity_run_inputs = sensitivity_test_methods
    .combine(sensitivity_combined_vcf_ch)
    .combine(sens_pheno_covar_ch)
    .combine(sens_refflat_ch)

process SensitivityRVtestRun {
    // executor 'slurm'
    // queue 'gr10478b'
    time '48h'
    tag "${filter_tag}_${method_tag}"

    publishDir "${params.outDir}/07.sensitivity_check/01.rvtest_run/${filter_tag}/${method_tag}", mode: 'symlink'

    input:
    tuple val(method_tag), val(method_opt), 
          val(filter_tag), file(vcf), file(vcf_tbi), 
          file(pheno), file(covar), 
          file(refflat) from sensitivity_run_inputs

    output:
    file("*.log")
    tuple val(filter_tag), val(method_tag), file("*.assoc") into sensitivity_result_ch

    script:
    out_prefix="cteph_agp3k.rare.${filter_tag}.${method_tag}"
    """
    export PATH=/home/b/b37974/rvtests/executable/:$PATH
    rvtest \
        --inVcf ${vcf} \
        --pheno ${pheno} --pheno-name pheno1 \
        --covar ${covar} --covar-name sex,pc1_avg,pc2_avg,pc3_avg,pc4_avg,pc5_avg,pc6_avg,pc7_avg,pc8_avg,pc9_avg,pc10_avg \
        --geneFile ${refflat} \
        --out ${out_prefix} \
        --noweb \
        --numThread 4 \
        ${method_opt}
    """
}

process SensitivityRVtestPostProcess {
    // executor 'slurm'
    // queue 'gr10478b'
    time '1h'
    tag "${filter_tag}_${method_tag}"

    publishDir "${params.outDir}/07.sensitivity_check/02.post_process/${filter_tag}/${method_tag}", mode: 'symlink'

    input:
    tuple val(filter_tag), val(method_tag), file(assoc_file) from sensitivity_result_ch

    output:
    tuple val(filter_tag), val(method_tag), file("*.fdr.assoc") into sensitivity_final_ch

    script:
    out_file = "${assoc_file.baseName}.filtered.fdr.assoc"
    """
    source activate cteph_geno_pro
    python ${params.scriptDir}/rvtest_post_process.py \
        --input ${assoc_file} \
        --output ${out_file} \
        --num-var-threshold ${params.num_var_threshold}
    """
}

process SensitivityRVtestVisualization {
    // executor 'slurm'
    // queue 'gr10478b'
    time '1h'
    tag "${filter_tag}_${method_tag}"

    publishDir "${params.outDir}/07.sensitivity_check/03.visualization/${filter_tag}/${method_tag}", mode: 'symlink'

    input:
    tuple val(filter_tag), val(method_tag), file(assoc_file) from sensitivity_final_ch

    output:
    tuple val(filter_tag), val(method_tag), file("*.png"), file("*.pdf") into sensitivity_plots_ch

    script:
    out_prefix = "${assoc_file.baseName}"
    title = "Sensitivity: ${filter_tag} - ${method_tag}"
    """
    source activate cteph_geno_pro
    python ${params.scriptDir}/plot_manhattan_qq.py \
        --input ${assoc_file} \
        --output-prefix ${out_prefix} \
        --title "${title}"
    """
}




