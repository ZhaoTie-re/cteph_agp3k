params.gtPath = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/wgs'
params.snpeffDir = '/LARGE1/gr10478/platform/JHRPv4/workspace/pipeline/output/snpEff.v4.index'
params.scriptDir = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/analysis/assoc_rvtest/scripts'
params.outDir = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/analysis/assoc_rvtest/results'

params.numvarThr = 3 // 每个基因至少包含的变异数阈值

process RVtestPrepare {
    executor 'slurm'
    queue 'gr10478b'
    time '36h'

    publishDir "${params.outDir}/01.rvtest_prepare", mode: 'symlink'

    input:
    val gtPath from params.gtPath

    output:
    file('*.log')
    tuple file("*.vcf.gz"), file("*.vcf.gz.tbi") into rvtest_vcf_ch
    tuple file("*.pheno_df.csv"), file("*.cov_df.no_age.csv") into rvtest_pheno_covar_ch
    file("refFlat.hg38.nochr.txt.gz") into rvtest_refflat_ch

    script:
    bed_prefix = "${gtPath}/19.tommo_panel_filter/cteph_agp3k.rare"
    pheno = "${gtPath}/20.cov_pheno_prepare/cteph_agp3k.bbj.projection.pheno_df.csv"
    covar = "${gtPath}/20.cov_pheno_prepare/cteph_agp3k.bbj.projection.cov_df.no_age.csv"
    refflat = "/LARGE0/gr10478/b37974/Pulmonary_Hypertension/anno_raw/refFlat.hg38.txt.gz"
    """
    source activate cteph_geno_pro
    python ${params.scriptDir}/rvtest_prepare_main.py \
        --bed-prefix ${bed_prefix} \
        --pheno-path ${pheno} \
        --covar-path ${covar} \
        --refflat-path ${refflat} \
        --threads 32 \
        --verbose
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

process InfoFilter {
    executor 'slurm'
    queue 'gr10478b'
    time '48h'

    publishDir "${params.outDir}/03.info_filter", mode: 'symlink'

    input:
    tuple file(vcf), file(vcf_tbi) from snpeff_annotated_vcf_ch

    output:
    file('*.log')
    tuple file("${out_prefix}.vcf.gz"), file("${out_prefix}.vcf.gz.tbi") into info_filtered_vcf_ch, info_filtered_vcf_ch2

    script:
    out_prefix='cteph_agp3k.rare.impact_moderate_high'
    """
    source activate cteph_geno_pro
    python ${params.scriptDir}/info_filter_main.py \
        --input ${vcf} \
        --info-key impact \
        --values MODERATE HIGH \
        --out-prefix ${out_prefix} \
        --threads 16 \
        --check-chr-prefix \
        --no-keep-chr-prefix \
        --verbose
    """
}

process AddMacInfo {
    executor 'slurm'
    queue 'gr10478b'
    time '6h'

    publishDir "${params.outDir}/03.info_filter", mode: 'symlink'

    input:
    tuple file(vcf), file(vcf_tbi) from info_filtered_vcf_ch2

    output:
    tuple file("*.vcf.gz"), file("*.vcf.gz.tbi") into info_filtered_vcf_mac_ch

    script:
    """
    source activate cteph_geno_pro
    bcftools +fill-tags ${vcf} -Oz -o cteph_agp3k.rare.impact_moderate_high.mac.vcf.gz --threads 16 -- -t AF,AC,AN,MAF
    tabix -p vcf cteph_agp3k.rare.impact_moderate_high.mac.vcf.gz
    """
}

macmin_opts = Channel.from([
    ["no_macmin", ""],               // 不指定--siteMACMin
    ["macmin_5",  "--siteMACMin 5"], // 指定为5
    ["macmin_10", "--siteMACMin 10"] // 指定为10
])

process RVtestRun {
    executor 'slurm'
    queue 'gr10478b'
    time '48h'
    tag "${macmin_tag}"

    publishDir "${params.outDir}/04.rvtest_run/${macmin_tag}", mode: 'symlink'

    input:
    tuple val(macmin_tag), val(macmin_opt) from macmin_opts
    tuple file(vcf), file(vcf_tbi) from info_filtered_vcf_ch
    tuple file(pheno), file(covar) from rvtest_pheno_covar_ch
    file(refflat) from rvtest_refflat_ch

    output:
    file("*.log")
    tuple val(macmin_tag), file("cteph_agp3k.rare.${macmin_tag}.SkatO.assoc") into skato_result_ch, skato_result_ch2, skato_result_ch3

    script:
    """
    export PATH=/home/b/b37974/rvtests/executable/:$PATH
    rvtest \
        --inVcf ${vcf} \
        --pheno ${pheno} --pheno-name pheno1 \
        --covar ${covar} --covar-name sex,pc1_avg,pc2_avg,pc3_avg,pc4_avg,pc5_avg,pc6_avg,pc7_avg,pc8_avg,pc9_avg,pc10_avg \
        --geneFile ${refflat} \
        --kernel skato \
        --out cteph_agp3k.rare.${macmin_tag} \
        --noweb \
        --numThread 16 \
        ${macmin_opt}
    """
}

process JsonManifest {
    executor 'slurm'
    queue 'gr10478b'
    time '24h'
    tag "json manifest"

    publishDir "${params.outDir}/05.json_manifest", mode: 'symlink'

    input:
    tuple val(tag_no_macmin), file(assoc_no_macmin) from skato_result_ch.filter { it[0] == 'no_macmin' }
    tuple val(tag_macmin_5), file(assoc_macmin_5) from skato_result_ch2.filter { it[0] == 'macmin_5' }
    tuple val(tag_macmin_10), file(assoc_macmin_10) from skato_result_ch3.filter { it[0] == 'macmin_10' }
    tuple file(vcf), file(vcf_tbi) from info_filtered_vcf_mac_ch

    output:
    file("*.json") into json_manifest_ch

    script:
    """
    source activate cteph_geno_pro
    python ${params.scriptDir}/create_file_manifest.py \
        --vcf ${vcf} \
        --assoc-no-macmin ${assoc_no_macmin} \
        --assoc-macmin-5 ${assoc_macmin_5} \
        --assoc-macmin-10 ${assoc_macmin_10} \
        --results-dir ${params.outDir} \
        --scripts-dir ${params.scriptDir} \
        --output cteph_agp3k_file_manifest.json
    """
}

process SummaryStats {
    executor 'slurm'
    queue 'gr10478b'
    time '24h'
    tag "summary stats"

    publishDir "${params.outDir}/06.summary_stats", mode: 'symlink'

    input:
    file(json_manifest) from json_manifest_ch
    val numvar_thr from params.numvarThr

    output:
    file('*.tsv')
    file('*.csv')
    tuple file("*.allele_counts.tsv.gz"), file("*.allele_counts.tsv.gz.tbi") into allele_counts_tsv_ch
    file('*.updated.json') into updated_json_manifest_ch

    script:
    """
    source activate cteph_geno_pro
    python ${params.scriptDir}/summary_main.py \
        --json-path ${json_manifest} \
        --num-var-thr ${numvar_thr} \
        --threads 8 \
        --verbose
    """
}
