params.gtPath = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/wgs'
params.tommoPath = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/ToMMo_60KJPN'
params.jhrpv4Path = '/LARGE1/gr10478/platform/JHRPv4/workspace/pipeline/output/VQSR.v4'
params.scriptDir = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/analysis/assoc_plink2/scripts'
params.outDir = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/analysis/assoc_plink2/results'

// Define models to test
model_ch = Channel.of('additive', 'dominant', 'recessive')
gtPath_ch = Channel.value(params.gtPath)

// Combine the model and gtPath channels to create a channel of tuples
assoc_input_ch = gtPath_ch.combine(model_ch)

process AssocPlink2Model {
    // executor 'slurm'
    // queue 'gr10478b'
    // time '6d'
    tag "assoc_plink2: ${model}"

    publishDir "${params.outDir}/01.assoc_result/${model}", mode: 'symlink'

    input:
    tuple val(gtPath), val(model) from assoc_input_ch

    output:
    file("*.log")
    tuple val(model), file("*.glm.logistic") \
          into assoc_plink2_out1, assoc_plink2_out2, assoc_plink2_out3

    script:
    BED_FILE_PR = "${gtPath}/19.sep_maf/cteph_agp3k.clean.rm_mono.not_rare"
    PHENO_FILE = "${gtPath}/18.clean_gt/pheno.txt"
    PHENO_NAME = "PHENO1"
    COVAR_FILE = "${gtPath}/18.clean_gt/covariance.txt"
    COVAR_NAME = "SEX, PC1_AVG-PC10_AVG"
    OUT_PR = "cteph_agp3k.sex.10pc.${model}"

    MODEL_OPT = model == 'additive' ? '--glm' : "--glm ${model}"

    """
    export PATH=/home/b/b37974/:$PATH
    plink2 \\
        --bfile ${BED_FILE_PR} \\
        --pheno ${PHENO_FILE} \\
        --pheno-name ${PHENO_NAME} \\
        --covar ${COVAR_FILE} \\
        --covar-name ${COVAR_NAME} \\
        ${MODEL_OPT} omit-ref no-firth hide-covar \\
        --out ${OUT_PR} \\
        --ci 0.95 \\
        --threads 8
    """
}

process AssocPlink2Summary {
    // executor 'slurm'
    // queue 'gr10478b'
    // time '6d'
    tag "assoc_plink2: summary"

    publishDir "${params.outDir}/02.summary_vis", mode: 'symlink'

    input:
    tuple val(tagADD), file(glmADD) from assoc_plink2_out1.filter { it[0] == 'additive' }
    tuple val(tagDOM), file(glmDOM) from assoc_plink2_out2.filter { it[0] == 'dominant' }
    tuple val(tagREC), file(glmREC) from assoc_plink2_out3.filter { it[0] == 'recessive' }
    val(gtPath) from params.gtPath
    val(tommoPath) from params.tommoPath
    val(jhrpv4Path) from params.jhrpv4Path

    output:
    file("*.csv")

    script:
    BED_FILE_PR = "${gtPath}/19.sep_maf/cteph_agp3k.clean.rm_mono.not_rare"
    TOMMO_VCF = "${tommoPath}/tommo-60kjpn-20240904-GRCh38-snvindel-af-autosome.norm.vcf.gz"
    TOMMO_VCF_TBI = "${tommoPath}/tommo-60kjpn-20240904-GRCh38-snvindel-af-autosome.norm.vcf.gz.tbi"
    """
    source activate cteph_geno_pro
    python ${params.scriptDir}/gwas_summary.plink2.py \\
        --bed_prefix ${BED_FILE_PR} \\
        --tommo_vcf_file ${TOMMO_VCF} \\
        --jhrp4_vcf_path ${jhrpv4Path} \\
        --jhrp4_vcf_prefix all.VQSR3 \\
        --p_threshold 5e-8 \\
        --add_path ${glmADD} \\
        --dom_path ${glmDOM} \\
        --rec_path ${glmREC} \\
    """
}
