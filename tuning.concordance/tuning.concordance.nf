params.scriptDir = "/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/script"
params.infoDir = "/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/info"
params.wgsDir = "/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/wgs"
params.arrayDir = "/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/array"

params.tmpDir = "/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/tuning.concordance/tmp"
params.outDir = "/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/tuning.concordance"

min_dp_channel = Channel.from(1..30)
min_gq_channel = Channel.of(10, 20, 30)

process TrueCallPrepare {

    executor 'slurm'
    queue 'gr10478b'
    time '6d'

    tag "Prepare true & call VCFs"

    publishDir "${params.outDir}/01.prepare_call_true", mode: 'symlink'

    input:
    val(arrayDir) from params.arrayDir
    val(wgsDir) from params.wgsDir

    output:
    tuple file("array_common_reheadered.vcf.gz"), file("array_common_reheadered.vcf.gz.tbi") into true_vcf
    tuple file("wgs_common_reheadered.vcf.gz"), file("wgs_common_reheadered.vcf.gz.tbi") into call_vcf, call_vcf_2
    file('common_ids.txt') into common_ids
    file('common_samples.txt')

    script:
    array_vcf = "${arrayDir}/07.z3.re.select.ids.vcf/cteph_agp3k.ajsa.qc.rechr.norm.reselect.ids.vcf.gz"
    wgs_vcf = "${wgsDir}/08.variant_filter/chr21.pass.mac1.het.smiss.v_filtered.vcf.gz"
    """
    chmod +x ${params.scriptDir}/concordance.vcf.prepare.sh
    ${params.scriptDir}/concordance.vcf.prepare.sh ${array_vcf} ${wgs_vcf}
    """
}

min_dp_channel
    .combine(min_gq_channel)
    .set { combined_channel }

process TrueCallConcordance {

    executor 'slurm'
    queue 'gr10478b'
    time '6d'

    tag "${dp} ${gq} Concordance"

    publishDir "${params.tmpDir}", mode: 'symlink'

    input:
    tuple file(true_vcf), file(true_vcf_tbi) from true_vcf
    tuple file(call_vcf), file(call_vcf_tbi) from call_vcf
    tuple val(dp), val(gq) from combined_channel

    output:
    tuple val(dp), val(gq), file("*.genotype_concordance_contingency_metrics") into gt_contingency
    tuple val(dp), val(gq), file("*.genotype_concordance_detail_metrics") into gt_detail
    tuple val(dp), val(gq), file("*.genotype_concordance_summary_metrics") into gt_summary

    script:
    """
    chmod +x ${params.scriptDir}/concordance.dp.gq.sh
    ${params.scriptDir}/concordance.dp.gq.sh ${true_vcf} ${call_vcf} ${dp} ${gq}
    """
}

process ConcordanceSummary {

    executor 'slurm'
    queue 'gr10478b'
    time '6d'

    tag "${dp} ${gq} Concordance summary"

    publishDir "${params.outDir}/02.concordance_summary", mode: 'symlink'

    input:
    tuple val(dp), val(gq), file(detail_mx) from gt_detail
    val(infoDir) from params.infoDir
    file(common_ids) from common_ids

    output:
    tuple val(dp), val(gq), file("DP${dp}_GQ${gq}_dict.pkl") into concordance_summary

    script:
    group_info = "${infoDir}/wgs_array_dp.csv"
    """
    source activate cteph_geno_pro
    python ${params.scriptDir}/stats_rev4.py --detail_mx ${detail_mx} --group_info ${group_info} --DP ${dp} --GQ ${gq} --common_ids ${common_ids}
    """
}

process DPGQMatrix {

    executor 'slurm'
    queue 'gr10478b'
    time '6d'

    tag "DP & GQ matrix"

    publishDir "${params.outDir}/03.dp_gq_matrix", mode: 'symlink'

    input:
    val(wgsDir) from params.wgsDir

    output:
    tuple file("wgs_gt_matrix.txt"), file("wgs_gq_matrix.txt"), file("wgs_dp_matrix.txt") into wgs_gq_dp_matrix

    script:
    wgs_vcf = "${wgsDir}/08.variant_filter/chr21.pass.mac1.het.smiss.v_filtered.vcf.gz"
    """
    bcftools view -v snps --threads 4 ${wgs_vcf} | bcftools query -f "[%GT\t]\n" > wgs_gt_matrix.txt
    bcftools view -v snps --threads 4 ${wgs_vcf} | bcftools query -f "[%GQ\t]\n" > wgs_gq_matrix.txt
    bcftools view -v snps --threads 4 ${wgs_vcf} | bcftools query -f "[%DP\t]\n" > wgs_dp_matrix.txt
    """
}

min_dp_channel_2 = Channel.from(1..30)
min_gq_channel_2 = Channel.of(10, 20, 30)

min_dp_channel_2
    .combine(min_gq_channel_2)
    .set { combined_channel_2 }

process DPGQVmiss {

    executor 'slurm'
    queue 'gr10478b'
    time '6d'

    tag "${dp} ${gq} Vmiss"

    publishDir "${params.outDir}/04.dp_gq_vmiss", mode: 'symlink'

    input:
    tuple val(dp), val(gq) from combined_channel_2
    tuple file(gt_mx), file(gq_mx), file(dp_mx) from wgs_gq_dp_matrix

    output:
    tuple val(dp), val(gq), file("DP${dp}_GQ${gq}_vmissing_rate.pkl") into vmiss_dp_gq

    script:
    """
    source activate cteph_geno_pro
    python ${params.scriptDir}/vmiss.dp.gq.py --DP ${dp} --GQ ${gq} --GT_MX ${gt_mx} --GQ_MX ${gq_mx} --DP_MX ${dp_mx}
    """
}

process DPGQMatrix_CallVCF {

    executor 'slurm'
    queue 'gr10478b'
    time '6d'

    tag "DP & GQ matrix"

    publishDir "${params.outDir}/05.call_vcf_dp_gq_matrix", mode: 'symlink'

    input:
    tuple file(call_vcf), file(call_vcf_tbi) from call_vcf_2

    output:
    tuple file("call_wgs_gt_matrix.txt"), file("call_wgs_gq_matrix.txt"), file("call_wgs_dp_matrix.txt") into call_wgs_gq_dp_matrix

    script:
    """
    bcftools query -f "[%GT\t]\n" ${call_vcf} > call_wgs_gt_matrix.txt
    bcftools query -f "[%GQ\t]\n" ${call_vcf} > call_wgs_gq_matrix.txt
    bcftools query -f "[%DP\t]\n" ${call_vcf} > call_wgs_dp_matrix.txt
    """
}

min_dp_channel_3 = Channel.from(1..30)
min_gq_channel_3 = Channel.of(10, 20, 30)

min_dp_channel_3
    .combine(min_gq_channel_3)
    .set { combined_channel_3 }

process DPGQVmiss_CallVCF {

    executor 'slurm'
    queue 'gr10478b'
    time '6d'

    tag "${dp} ${gq} Vmiss"

    publishDir "${params.outDir}/06.call_vcf_dp_gq_vmiss", mode: 'symlink'

    input:
    tuple val(dp), val(gq) from combined_channel_3
    tuple file(gt_mx), file(gq_mx), file(dp_mx) from call_wgs_gq_dp_matrix

    output:
    tuple val(dp), val(gq), file("DP${dp}_GQ${gq}_vmissing_rate.pkl") into call_vmiss_dp_gq

    script:
    """
    source activate cteph_geno_pro
    python ${params.scriptDir}/vmiss.dp.gq.py --DP ${dp} --GQ ${gq} --GT_MX ${gt_mx} --GQ_MX ${gq_mx} --DP_MX ${dp_mx}
    """
}