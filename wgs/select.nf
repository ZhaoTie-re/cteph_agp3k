params.filePath = '/LARGE1/gr10478/platform/JHRPv4/workspace/pipeline/output/VQSR.v4'
params.infoPath = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/info'
params.samplelist = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/info/cteph_agp3k_wgs_ids.txt'
params.tommodir = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/ToMMo_60KJPN'
params.sifdir = '/home/b/b37974/simg'
params.scriptDir = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/script'
params.outdir = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/wgs'

params.cteph_agp3k_main = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k'
params.bbj_main = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/BBJ_genome_b38'

Channel
    .from((1..22).collect { "chr${it}" } + ["PAR"])
    .set { chr_ch }

process selectPASS {
    executor 'slurm'
    queue 'gr10478b'
    time '36h'
    tag "${chr}"

    publishDir "${params.outdir}/01.pass", mode: 'symlink'

    input:
    val(filePath) from params.filePath
    val(chr) from chr_ch
    path(samplelist) from params.samplelist

    output:
    tuple chr, file(pass_vcf), file(pass_vcf_tbi) into selectPASS_out

    script:
    vcf = "${filePath}/all.VQSR3.${chr}.vcf.gz"
    pass_vcf = "${chr}.pass.vcf.gz"
    pass_vcf_tbi = "${chr}.pass.vcf.gz.tbi"

    """
    bcftools view ${vcf} --threads 2 -S ${samplelist} --force-samples -Ou | bcftools view --threads 2 -f"PASS" -Oz -o ${pass_vcf}
    bcftools index --threads 2 -t ${pass_vcf}
    """
}

process filter_ac1 {
    executor 'slurm'
    queue 'gr10478b'
    time '36h'
    tag "${chr}"

    publishDir "${params.outdir}/02.mac1", mode: 'symlink'

    input:
    tuple chr, file(vcf), file(vcf_tbi) from selectPASS_out

    output:
    tuple chr, file(mac1_vcf), file(mac1_vcf_tbi) into filter_ac1_out

    script:
    mac1_vcf = "${chr}.pass.mac1.vcf.gz"
    mac1_vcf_tbi = "${chr}.pass.mac1.vcf.gz.tbi"

    """
    bcftools view ${vcf} --threads 2 --min-ac 1:nref -Oz -o ${mac1_vcf}
    bcftools index --threads 2 -t ${mac1_vcf}
    """
}

process set_ids {
    executor 'slurm'
    queue 'gr10478b'
    time '36h'
    tag "${chr}"

    publishDir "${params.outdir}/03.set_ids", mode: 'symlink'

    input:
    tuple chr, file(vcf), file(vcf_tbi) from filter_ac1_out

    output:
    tuple chr, file(setid_vcf), file(setid_vcf_tbi) into set_ids_out

    script:
    setid_vcf = "${chr}.pass.mac1.setid.vcf.gz"
    setid_vcf_tbi = "${chr}.pass.mac1.setid.vcf.gz.tbi"

    """
    bcftools annotate --set-id '%CHROM:%POS:%REF:%ALT' ${vcf} -Oz -o ${setid_vcf}
    bcftools index --threads 2 -t ${setid_vcf}
    """
}

process add_gt_AF {
    executor 'slurm'
    queue 'gr10478b'
    time '36h'
    tag "${chr}"

    publishDir "${params.outdir}/04.add_gt_AF", mode: 'symlink'

    input:
    tuple chr, file(vcf), file(vcf_tbi) from set_ids_out
    val(sifdir) from params.sifdir

    output:
    tuple chr, file(gt_af_vcf), file(gt_af_vcf_tbi) into add_gt_af_out

    script:
    gatk_sif = "${sifdir}/gatk_latest.sif"
    gt_af_vcf = "${chr}.pass.mac1.setid.gt_af.vcf.gz"
    gt_af_vcf_tbi = "${chr}.pass.mac1.setid.gt_af.vcf.gz.tbi"

    """
    singularity exec --bind /LARGE0:/LARGE0 ${gatk_sif} gatk VariantAnnotator \
        -R /LARGE0/gr10478/b37974/Pulmonary_Hypertension/nagasaki_pipeline/data/hs38DH.fa \
        -V ${vcf} \
        -O ${gt_af_vcf} \
        -A AlleleFraction \
        --create-output-variant-index true \
        --java-options "-XX:ParallelGCThreads=4"
    """
}

process gt_norm {
    executor 'slurm'
    queue 'gr10478b'
    time '36h'
    tag "${chr}"

    publishDir "${params.outdir}/05.gt_norm", mode: 'symlink'

    input:
    tuple chr, file(tmp_vcf), file(tmp_vcf_tbi) from add_gt_af_out

    output:
    tuple chr, file(gt_norm_vcf), file(gt_norm_vcf_tbi) into gt_norm_out

    script:
    gt_norm_vcf = "${chr}.pass.mac1.setid.gt_af.norm.vcf.gz"
    gt_norm_vcf_tbi = "${chr}.pass.mac1.setid.gt_af.norm.vcf.gz.tbi"
    """
    bcftools view ${tmp_vcf} | sed 's/nan/NaN/g' | bgzip > ${gt_norm_vcf}
    bcftools +setGT ${gt_norm_vcf} -- -t a -n u
    bcftools index --threads 4 -t ${gt_norm_vcf}
    """
}

gt_norm_out
    .filter { chr, vcf, vcf_tbi -> chr != "PAR" }
    .set {autosome_vcf_ch}

process variant_filter {
    executor 'slurm'
    queue 'gr10478b'
    time '36h'
    tag "${chr}"

    publishDir "${params.outdir}/06.variant_filter", mode: 'symlink'

    input:
    tuple chr, file(vcf), file(vcf_tbi) from autosome_vcf_ch

    output:
    tuple chr, file(filtered_vcf), file(filtered_vcf_tbi) into filtered_vcf_out

    script:
    filtered_vcf = "${chr}.pass.mac1.vfilter.vcf.gz"
    filtered_vcf_tbi = "${chr}.pass.mac1.vfilter.vcf.gz.tbi"

    """
    bcftools view ${vcf} --threads 8 -i 'VQSLOD > 10 & MQ > 58.75' -Oz -o ${filtered_vcf}
    bcftools index --threads 8 -t ${filtered_vcf}
    """
}

process gt_qc {
    executor 'slurm'
    queue 'gr10478b'
    time '36h'
    tag "${chr}"

    publishDir "${params.outdir}/07.gt_qc", mode: 'symlink'

    input:
    tuple chr, file(vcf), file(vcf_tbi) from filtered_vcf_out
    val(sifdir) from params.sifdir

    output:
    tuple chr, file(gt_qc_vcf), file(gt_qc_vcf_tbi) into gt_qc_out

    script:
    gatk_sif = "${sifdir}/gatk_latest.sif"
    gt_qc_vcf = "${chr}.pass.mac1.vfilter.gt_qc.vcf.gz"
    gt_qc_vcf_tbi = "${chr}.pass.mac1.vfilter.gt_qc.vcf.gz.tbi"

    """
    singularity exec --bind /LARGE0:/LARGE0 ${gatk_sif} gatk --java-options "-Xmx1G" VariantFiltration \
        -R /LARGE0/gr10478/b37974/Pulmonary_Hypertension/nagasaki_pipeline/data/hs38DH.fa \
        -V ${vcf} \
        -O ${gt_qc_vcf} \
        --genotype-filter-name "LowGQ" \
        --genotype-filter-expression "GQ<20" \
        --genotype-filter-name "LowDP" \
        --genotype-filter-expression "DP<9" \
        --set-filtered-genotype-to-no-call true \
        --create-output-variant-index true
    """
}

process gt_norm_2 {
    executor 'slurm'
    queue 'gr10478b'
    time '36h'
    tag "${chr}"

    publishDir "${params.outdir}/08.gt_norm_2", mode: 'symlink'

    input:
    tuple chr, file(vcf), file(vcf_tbi) from gt_qc_out

    output:
    tuple chr, file(gt_norm_vcf_2), file(gt_norm_vcf_tbi_2) into gt_norm_out_2

    script:
    gt_norm_vcf_2 = "${chr}.pass.mac1.vfilter.gt_qc.norm.vcf.gz"
    gt_norm_vcf_tbi_2 = "${chr}.pass.mac1.vfilter.gt_qc.norm.vcf.gz.tbi"
    """
    bcftools view ${vcf} | sed 's/nan/NaN/g' | bgzip > ${gt_norm_vcf_2}
    bcftools index --threads 4 -t ${gt_norm_vcf_2}
    """
}

process het_gt_qc {
    executor 'slurm'
    queue 'gr10478b'
    time '36h'
    tag "${chr}"

    publishDir "${params.outdir}/09.het_qt_qc", mode: 'symlink'

    input:
    tuple chr, file(vcf), file(vcf_tbi) from gt_norm_out_2
    val(sifdir) from params.sifdir

    output:
    tuple chr, file(het_qc_vcf), file(het_qc_vcf_tbi) into het_qc_out

    script:
    gatk_sif = "${sifdir}/gatk_latest.sif"
    het_qc_vcf = "${chr}.pass.mac1.vfilter.gt_qc.het_qc.vcf.gz"
    het_qc_vcf_tbi = "${chr}.pass.mac1.vfilter.gt_qc.het_qc.vcf.gz.tbi"
    """
    singularity exec --bind /LARGE0:/LARGE0 ${gatk_sif} gatk --java-options "-Xmx1G" VariantFiltration \
        -R /LARGE0/gr10478/b37974/Pulmonary_Hypertension/nagasaki_pipeline/data/hs38DH.fa \
        -V ${vcf} \
        -O ${het_qc_vcf} \
        --genotype-filter-name "ABB_outlier" \
        --genotype-filter-expression "isHet == 1 && (AF < 0.25 || AF > 0.75)" \
        --genotype-filter-name "ABB_NaN" \
        --genotype-filter-expression "AF == 'NaN'" \
        --set-filtered-genotype-to-no-call true \
        --create-output-variant-index true
    """
}

// process autosome_vcf2bed {
//     executor 'slurm'
//     queue 'gr10478b'
//     time '36h'
//     tag "${chr}"

//     publishDir "${params.outdir}/06.autosome_vcf2bed", mode: 'symlink'

//     input:
//     tuple chr, file(vcf), file(vcf_tbi) from autosome_vcf_ch

//     output:
//     tuple chr, file("*.bed"), file("*.bim"), file("*.fam") into bed_ch

//     script:
//     bed_prefix = "${chr}.pass.mac1"
//     """
//     export PATH=/home/b/b37974/:$PATH
//     plink \
//         --vcf ${vcf} \
//         --make-bed \
//         --keep-allele-order \
//         --double-id \
//         --out ${bed_prefix}
//     """
// }

// bed_ch
//     .toList()
//     .map { entries ->
//         def chr_order = (1..22).collect { "chr${it}" }
//         def filtered = entries.findAll { it[0] in chr_order }
//         def sorted = filtered.sort { a, b ->
//             Integer.parseInt(a[0].replace("chr", "")) <=> Integer.parseInt(b[0].replace("chr", ""))
//         }
//         def paths = sorted.collect { it[1].toString().replaceAll(/\.bed$/, '') }
//         return paths.join('\n')
//     }
//     .set { pmerge_list_content }

// process autosome_writePmergeList {
//     executor 'slurm'
//     queue 'gr10478b'
//     time '1h'
//     publishDir "${params.outdir}/06.autosome_vcf2bed", mode: 'copy'

//     input:
//     val(list_text) from pmerge_list_content

//     output:
//     file("pmerge_lst.txt") into pmerge_lst_ch

//     script:
//     """
//     echo "${list_text}" > pmerge_lst.txt
//     """
// }

// process autosome_pmerge {
//     executor 'slurm'
//     queue 'gr10478b'
//     time '36h'
//     tag "pmerge"

//     publishDir "${params.outdir}/06.autosome_vcf2bed", mode: 'symlink'

//     input:
//     file(pmerge_list) from pmerge_lst_ch

//     output:
//     tuple file("*.bed"), file("*.bim"), file("*.fam") into pmerge_out

//     script:
//     out_prefix = "cteph_agp3k.all_chr.no_qc"
//     """
//     export PATH=/home/b/b37974/:$PATH
//     plink2 \
//         --pmerge-list ${pmerge_list} bfile \
//         --threads 16 \
//         --make-bed \
//         --out ${out_prefix}
//     """
// }

// process sample_qc_miss_het {
//     executor 'slurm'
//     queue 'gr10478b'
//     time '36h'
//     tag "sample_qc_miss_het"

//     publishDir "${params.outdir}/06.autosome_vcf2bed", mode: 'symlink'

//     input:
//     tuple file(bed), file(bim), file(fam) from pmerge_out

//     output:
//     file('*.pdf')
//     file('*.txt') into pass_miss_het_out

//     script:
//     bed_prefix = bed.baseName
//     """
//     source activate cteph_geno_pro
//     python ${params.scriptDir}/sample.qc.hetero.miss.py \
//         --bed_prefix ${bed_prefix} \
//         --threads 8 \
//         --miss_threshold 0.05 \
//         --het_threshold 5sd \
//     """
// }

// gt_norm_out_2
//     .filter { chr, vcf, vcf_tbi -> chr != "PAR" }
//     .set {autosome_vcf_ch_2}

// process select_pass_qc_samples {
//     executor 'slurm'
//     queue 'gr10478b'
//     time '36h'
//     tag "${chr}"

//     publishDir "${params.outdir}/07.select_pass_qc_samples", mode: 'symlink'

//     input:
//     tuple chr, file(vcf), file(vcf_tbi) from autosome_vcf_ch_2
//     file(pass_qc_sample) from pass_miss_het_out

//     output:
//     tuple chr, file(pass_vcf), file(pass_vcf_tbi) into vcf_sample_qc_out

//     script:
//     pass_vcf = "${chr}.pass.mac1.het.smiss.vcf.gz"
//     pass_vcf_tbi = "${chr}.pass.mac1.het.smiss.vcf.gz.tbi"

//     """
//     bcftools view ${vcf} --threads 4 -S ${pass_qc_sample} --force-samples -Oz -o ${pass_vcf}
//     bcftools index --threads 4 -t ${pass_vcf}
//     """
// }

// process variant_filter {
//     executor 'slurm'
//     queue 'gr10478b'
//     time '36h'
//     tag "${chr}"

//     publishDir "${params.outdir}/08.variant_filter", mode: 'symlink'

//     input:
//     tuple chr, file(vcf), file(vcf_tbi) from vcf_sample_qc_out

//     output:
//     tuple chr, file(filtered_vcf), file(filtered_vcf_tbi) into filtered_vcf_out

//     script:
//     filtered_vcf = "${chr}.pass.mac1.het.smiss.v_filtered.vcf.gz"
//     filtered_vcf_tbi = "${chr}.pass.mac1.het.smiss.v_filtered.vcf.gz.tbi"

//     """
//     bcftools view ${vcf} --threads 8 -i 'VQSLOD > 10 & MQ > 58.75' -Oz -o ${filtered_vcf}
//     bcftools index --threads 8 -t ${filtered_vcf}
//     """
// }

// process gt_qc {
//     executor 'slurm'
//     queue 'gr10478b'
//     time '36h'
//     tag "${chr}"

//     publishDir "${params.outdir}/09.gt_qc", mode: 'symlink'

//     input:
//     tuple chr, file(vcf), file(vcf_tbi) from filtered_vcf_out
//     val(sifdir) from params.sifdir

//     output:
//     tuple chr, file(gt_qc_vcf), file(gt_qc_vcf_tbi) into gt_qc_out

//     script:
//     gatk_sif = "${sifdir}/gatk_latest.sif"
//     gt_qc_vcf = "${chr}.pass.mac1.het.smiss.v_filtered.gt_qc.vcf.gz"
//     gt_qc_vcf_tbi = "${chr}.pass.mac1.het.smiss.v_filtered.gt_qc.vcf.gz.tbi"

//     """
//     singularity exec --bind /LARGE0:/LARGE0 ${gatk_sif} gatk --java-options "-Xmx1G" VariantFiltration \
//         -R /LARGE0/gr10478/b37974/Pulmonary_Hypertension/nagasaki_pipeline/data/hs38DH.fa \
//         -V ${vcf} \
//         -O ${gt_qc_vcf} \
//         --genotype-filter-name "LowGQ" \
//         --genotype-filter-expression "GQ<20" \
//         --genotype-filter-name "LowDP" \
//         --genotype-filter-expression "DP<8" \
//         --set-filtered-genotype-to-no-call true \
//         --create-output-variant-index true
//     """
// }

// process filter_ac1_2 {
//     executor 'slurm'
//     queue 'gr10478b'
//     time '36h'
//     tag "${chr}"

//     publishDir "${params.outdir}/09.gt_qc", mode: 'symlink'

//     input:
//     tuple chr, file(vcf), file(vcf_tbi) from gt_qc_out

//     output:
//     tuple chr, file(mac1_vcf), file(mac1_vcf_tbi) into filter_ac1_2_out

//     script:
//     mac1_vcf = "${chr}.pass.mac1.het.smiss.v_filtered.gt_qc.mac1.vcf.gz"
//     mac1_vcf_tbi = "${chr}.pass.mac1.het.smiss.v_filtered.gt_qc.mac1.vcf.gz.tbi"

//     """
//     bcftools view ${vcf} --threads 4 --min-ac 1:nref -Oz -o ${mac1_vcf}
//     bcftools index --threads 4 -t ${mac1_vcf}
//     """
// }

// process vcf2bed {
//     executor 'slurm'
//     queue 'gr10478b'
//     time '36h'
//     tag "${chr}"

//     publishDir "${params.outdir}/10.vcf2bed", mode: 'symlink'

//     input:
//     tuple chr, file(vcf), file(vcf_tbi) from filter_ac1_2_out

//     output:
//     tuple chr, file("*.bed"), file("*.bim"), file("*.fam") into bed_ch_qc_out

//     script:
//     bed_prefix = "${chr}.pass.mac1.het.smiss.v_filtered.gt_qc"
//     """
//     export PATH=/home/b/b37974/:$PATH
//     plink \
//         --vcf ${vcf} \
//         --make-bed \
//         --keep-allele-order \
//         --double-id \
//         --out ${bed_prefix}
//     """
// }

// bed_ch_qc_out
//     .toList()
//     .map { entries ->
//         def chr_order = (1..22).collect { "chr${it}" }
//         def filtered = entries.findAll { it[0] in chr_order }
//         def sorted = filtered.sort { a, b ->
//             Integer.parseInt(a[0].replace("chr", "")) <=> Integer.parseInt(b[0].replace("chr", ""))
//         }
//         def paths = sorted.collect { it[1].toString().replaceAll(/\.bed$/, '') }
//         return paths.join('\n')
//     }
//     .set { pmerge_list_content_qc_out }

// process writePmergeList {
//     executor 'slurm'
//     queue 'gr10478b'
//     time '1h'
//     publishDir "${params.outdir}/10.vcf2bed", mode: 'copy'

//     input:
//     val(list_text) from pmerge_list_content_qc_out

//     output:
//     file("pmerge_lst.txt") into pmerge_lst_qc_out

//     script:
//     """
//     echo "${list_text}" > pmerge_lst.txt
//     """
// }

// process pmerge {
//     executor 'slurm'
//     queue 'gr10478b'
//     time '36h'
//     tag "pmerge"

//     publishDir "${params.outdir}/11.pmerge", mode: 'symlink'

//     input:
//     file(pmerge_list) from pmerge_lst_qc_out

//     output:
//     tuple file("*.bed"), file("*.bim"), file("*.fam") into after_qc_ch

//     script:
//     out_prefix = "cteph_agp3k.s_qc.gt_qc"
//     """
//     export PATH=/home/b/b37974/:$PATH
//     plink2 \
//         --pmerge-list ${pmerge_list} bfile \
//         --threads 4 \
//         --make-bed \
//         --out ${out_prefix}
//     """
// }

// process update_fam_pheno {

//     executor 'slurm'
//     queue 'gr10478b'
//     time '36h'
//     tag "update_fam_pheno"

//     publishDir "${params.outdir}/12.update_fam_pheno", mode: 'symlink'

//     input:
//     tuple file(bed), file(bim), file(fam) from after_qc_ch
//     val(infoPath) from params.infoPath

//     output:
//     file('*.log')
//     tuple file("*.bed"), file("*.bim"), file("*.fam") into update_fam_pheno_out

//     script:
//     bed_prefix = bed.baseName
//     out_prefix = "cteph_agp3k.s_qc.gt_qc.pheno"
//     info_df = "${infoPath}/cteph_agp3k_jhrpv4.xlsx"
//     """
//     source activate cteph_geno_pro
//     python ${params.scriptDir}/update.fam.pheno.py \
//         --info_file ${info_df} \
//         --bed_prefix ${bed_prefix} \
//         --case_prefix PHOM \
//         --out ${out_prefix}
//     """
// }

// process variant_qc_miss {
//     executor 'slurm'
//     queue 'gr10478b'
//     time '36h'
//     tag "variant_qc_miss"

//     publishDir "${params.outdir}/13.variant_qc_miss", mode: 'symlink'

//     input:
//     tuple file(bed), file(bim), file(fam) from update_fam_pheno_out

//     output:
//     file('*.log')
//     file('*.pdf')
//     file('*.png')
//     file('*.txt')
//     file('*.vmiss')
//     tuple file("*.bed"), file("*.bim"), file("*.fam") into variant_qc_miss_out

//     script:
//     bed_prefix = bed.baseName
//     out_prefix = "cteph_agp3k.s_qc.gt_qc.v_qc"
//     """
//     source activate cteph_geno_pro
//     python ${params.scriptDir}/variant.qc.miss.py \
//         --bed_prefix ${bed_prefix} \
//         --threads 4 \
//         --maf_threshold 0.01 \
//         --vmiss_threshold 0.04 \
//         --out ${out_prefix}
//     """
// }


// process hwe_check {
//     executor 'slurm'
//     queue 'gr10478b'
//     time '36h'
//     tag "hwe_qc"

//     publishDir "${params.outdir}/14.hwe_qc", mode: 'symlink'

//     input:
//     tuple file(bed), file(bim), file(fam) from variant_qc_miss_out

//     output:
//     file('*.log') 
//     file('*.pdf')
//     file('*.png')
//     file('*.afreq')
//     file('*.hardy')
//     tuple file("pass_hwe_ids.txt"), file("*.maf") into pass_hwe_ids_out
//     tuple file("*.bed"), file("*.bim"), file("*.fam") into hwe_qc_out, qc_done_ch, qc_done_ch_2, qc_done_ch_3

//     script:
//     bed_prefix = bed.baseName
//     out_prefix = "cteph_agp3k.s_qc.gt_qc.v_qc.hwe"
//     """
//     source activate cteph_geno_pro
//     python ${params.scriptDir}/hwe.check.rev1.py \
//         --bed_prefix ${bed_prefix} \
//         --threads 4 \
//         --plt_maf 0.01 \
//         --plt_hwe_control 1e-6 \
//         --plt_hwe_case 1e-10 \
//         --mode case_control \
//         --control_hwe_threshold 1e-6 \
//         --case_hwe_threshold 1e-10 \
//         --out ${out_prefix}
//     """
// }

// process pca {
//     executor 'slurm'
//     queue 'gr10478b'
//     time '36h'
//     tag "pca"

//     publishDir "${params.outdir}/15.pca", mode: 'symlink'

//     input:
//     tuple file(bed), file(bim), file(fam) from hwe_qc_out
//     tuple file(pass_hwe_ids), file(maf) from pass_hwe_ids_out
//     val(infoPath) from params.infoPath

//     output:
//     file('*.pdf')
//     file('*.log')
//     file('*.pca.eigenvec')
//     file('*.pca.eigenval')
//     tuple file("*.no_high_ld.bed"), file("*.no_high_ld.bim"), file("*.no_high_ld.fam") into no_high_ld_prune_bed_ch
//     tuple file("*.prune.in"), file("*.prune.out") into prune_in_out

//     script:
//     bed_prefix = bed.baseName
//     high_ld = "${infoPath}/high-LD-regions-hg38-GRCh38.txt"
//     """
//     source activate cteph_geno_pro
//     python ${params.scriptDir}/pca.run.py \
//         --hwe_pass_ids ${pass_hwe_ids} \
//         --check_maf ${maf} \
//         --maf_selection 0.01 \
//         --bed_prefix ${bed_prefix} \
//         --high_ld ${high_ld} \
//         --threads 8 \
//         --case_prefix PHOM \
//         --case_name CTEPH \
//         --control_name AGP3K
//     """
// }

// process bbj_prepare {
//     executor 'slurm'
//     queue 'gr10478b'
//     time '36h'
//     tag "bbj_prepare"

//     publishDir "${params.cteph_agp3k_main}/bbj_projection/01.bbj_prepare", mode: 'symlink'

//     input:
//     val(bbj_main) from params.bbj_main

//     output:
//     file('*.log')
//     tuple file("*.setid.bed"), file("*.setid.bim"), file("*.setid.fam") into bbj_prepare_out

//     script:
//     bed_prefix = "${bbj_main}/02.sample_qc/NewOE13_Auto.id.b38.sqc"
//     """
//     source activate cteph_geno_pro
//     python ${params.scriptDir}/bbj_prepare.py \
//         --bbj_bed_prefix ${bed_prefix} \
//         --threads 24 \
//     """
// }

// process bbj_pca {
//     executor 'slurm'
//     queue 'gr10478b'
//     time '36h'
//     tag "bbj_pca"

//     publishDir "${params.cteph_agp3k_main}/bbj_projection/02.bbj_pca", mode: 'symlink'

//     input:
//     tuple file(bed), file(bim), file(fam) from bbj_prepare_out
//     val(infoPath) from params.infoPath

//     output:
//     file('*.log')
//     file('*.pdf')
//     tuple file("*.no_high_ld.bed"), file("*.no_high_ld.bim"), file("*.no_high_ld.fam") into bbj_no_highld_bed
//     tuple file("*.no_high_ld.prune.in"), file("*.no_high_ld.prune.out") into bbj_prune_out
//     tuple file("*.no_high_ld.prune.pca.eigenvec"), file("*.no_high_ld.prune.pca.eigenval") into bbj_pca_out
//     tuple file("*.no_high_ld.prune.pca.acount"), file("*.no_high_ld.prune.pca.eigenvec.allele") into bbj_pca_acount_out 
    
//     script:
//     bed_prefix = bed.baseName
//     high_ld = "${infoPath}/high-LD-regions-hg38-GRCh38.txt"
//     """
//     source activate cteph_geno_pro
//     python ${params.scriptDir}/bbj_pca.rev2.py \
//         --bbj_bed_prefix ${bed_prefix} \
//         --output_prefix bbj.b38.auto.sqc.vqc.norm \
//         --high_ld ${high_ld} \
//         --threads 16
//     """
// }

// process bbj_projection {
//     executor 'slurm'
//     queue 'gr10478b'
//     time '36h'
//     tag "bbj_projection"

//     publishDir "${params.cteph_agp3k_main}/bbj_projection/03.bbj_projection", mode: 'symlink'

//     input:
//     tuple file(my_bed), file(my_bim), file(my_fam) from qc_done_ch
//     tuple file(bbj_bed), file(bbj_bim), file(bbj_fam) from bbj_no_highld_bed
//     tuple file(bbj_prune_in), file(bbj_prune_out) from bbj_prune_out
//     tuple file(bbj_pca_acount), file(bbj_pca_eigenvec_allele) from bbj_pca_acount_out

//     output:
//     file('*.log')
//     file('*.txt')
//     tuple file("*.sscore"), file("*.sscore.vars") into bbj_projection_out, bbj_projection_out_2

//     script:
//     my_bed_prefix = my_bed.baseName
//     bbj_bed_prefix = bbj_bed.baseName
//     my_prefix_out = "cteph_agp3k"
//     bbj_prefix_out = "bbj"
//     """
//     source activate cteph_geno_pro
//     python ${params.scriptDir}/bbj_projection.py \
//         --my_bed_prefix ${my_bed_prefix} \
//         --bbj_bed_prefix ${bbj_bed_prefix} \
//         --bbj_prune_in ${bbj_prune_in} \
//         --my_prefix_out ${my_prefix_out} \
//         --bbj_prefix_out ${bbj_prefix_out} \
//         --bbj_pca_acount ${bbj_pca_acount} \
//         --bbj_pca_eigenvec_allele ${bbj_pca_eigenvec_allele} \
//         --threads 16
//     """
// }

// process bbj_projection_sample_keep {
//     executor 'slurm'
//     queue 'gr10478b'
//     time '36h'
//     tag "bbj_projection_sample_keep"

//     publishDir "${params.outdir}/16.bbj_projection_sample_keep", mode: 'symlink'

//     input:
//     tuple file(sscore), file(sscore_vars) from bbj_projection_out

//     output:
//     file('*.png')
//     file('*.subset.keep.txt') into bbj_projection_keep_out, bbj_projection_keep_out_2

//     script:
//     """
//     source activate cteph_geno_pro
//     python ${params.scriptDir}/bbj_projection_keep.py \
//         --projection_sscore_path ${sscore} \
//         --rect_xlim -0.025 0.016 \
//         --rect_ylim -0.025 0.025
//     """
// }

// process kinship_sample_remove {
//     executor 'slurm'
//     queue 'gr10478b'
//     time '36h'
//     tag "kinship_sample_remove"

//     publishDir "${params.outdir}/17.kinship_sample_remove", mode: 'symlink'

//     input:
//     tuple file(my_bed), file(my_bim), file(my_fam) from qc_done_ch_2
//     file(sample_keep) from bbj_projection_keep_out
//     tuple file(my_prune_in), file(my_prune_out) from prune_in_out

//     output:
//     file('*.log')
//     file('*.pdf')
//     file('*.csv')
//     file('*.to_remove.txt') into kinship_sample_remove_out

//     script:
//     bed_prefix = my_bed.baseName
//     """
//     source activate cteph_geno_pro
//     python ${params.scriptDir}/kinship.pi_hat.rev1.py \
//         --bed_prefix ${bed_prefix} \
//         --sample_keep ${sample_keep} \
//         --prune_in ${my_prune_in} \
//         --subset_prefix projection_keep \
//         --out_prefix kinship_check \
//         --case_prefix PHOM \
//         --pi_hat_threshold 0.20 \
//         --threads 8 
//     """
// }

// process clean_gt {
//     executor 'slurm'
//     queue 'gr10478b'
//     time '36h'
//     tag "clean_gt"

//     publishDir "${params.outdir}/18.clean_gt", mode: 'symlink'

//     input:
//     tuple file(my_bed), file(my_bim), file(my_fam) from qc_done_ch_3
//     tuple file(sscore), file(sscore_vars) from bbj_projection_out_2
//     file(sample_keep) from bbj_projection_keep_out_2
//     file(sample_remove) from kinship_sample_remove_out

//     output:
//     file('*.log')
//     tuple file('pheno.txt'), file('covariance.txt') into pheno_cov_out
//     tuple file("*.bed"), file("*.bim"), file("*.fam") into clean_gt_out

//     script:
//     bed_prefix = my_bed.baseName
//     """
//     source activate cteph_geno_pro
//     python ${params.scriptDir}/clean_gt.py \
//         --bed_prefix ${bed_prefix} \
//         --projection_keep ${sample_keep} \
//         --kinship_remove ${sample_remove} \
//         --sscore ${sscore} \
//         --out_prefix cteph_agp3k.clean \
//         --threads 8
//     """
// }




