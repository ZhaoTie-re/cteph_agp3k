params.gtPath = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/wgs'
params.snpeffDir = '/LARGE1/gr10478/platform/JHRPv4/workspace/pipeline/output/snpEff.v4.index'
params.scriptDir = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/analysis/assoc_rvtest_rev1/scripts'
params.outDir = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/analysis/assoc_rvtest_rev1/results'

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

process AddInfo {
    executor 'slurm'
    queue 'gr10478b'
    time '6h'

    publishDir "${params.outDir}/03.add_info", mode: 'symlink'

    input:
    tuple file(vcf), file(vcf_tbi) from snpeff_annotated_vcf_ch

    output:
    tuple file("*.vcf.gz"), file("*.vcf.gz.tbi") into info_filtered_vcf_mac_ch

    script:
    """
    source activate cteph_geno_pro
    bcftools +fill-tags ${vcf} -Oz -o cteph_agp3k.rare.snpeff.info.vcf.gz --threads 16 -- -t AF,AC,AN,MAF
    tabix -p vcf cteph_agp3k.rare.snpeff.info.vcf.gz
    """
}

