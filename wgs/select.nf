params.filePath = '/LARGE1/gr10478/platform/JHRPv4/workspace/pipeline/output/VQSR.v4'
params.infoPath = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/info'
params.samplelist = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/info/cteph_agp3k_wgs_ids.txt'
params.tommodir = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/ToMMo_60KJPN'
params.sifdir = '/home/b/b37974/simg'
params.scriptDir = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/wgs/scripts'
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
        --genotype-filter-expression "DP<8" \
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
    gt_norm_vcf_2 = "${chr}.pass.mac1.vfilter.gt_qc.norm2.vcf.gz"
    gt_norm_vcf_tbi_2 = "${chr}.pass.mac1.vfilter.gt_qc.norm2.vcf.gz.tbi"
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
        --genotype-filter-expression "isHet == 1 && (AF < 0.2 || AF > 0.8)" \
        --genotype-filter-name "ABB_NaN" \
        --genotype-filter-expression "AF == 'NaN'" \
        --set-filtered-genotype-to-no-call true \
        --create-output-variant-index true
    """
}

process vcf2bed {
    executor 'slurm'
    queue 'gr10478b'
    time '36h'
    tag "${chr}"

    publishDir "${params.outdir}/10.vcf2bed", mode: 'symlink'

    input:
    tuple chr, file(vcf), file(vcf_tbi) from het_qc_out

    output:
    tuple chr, file("*.bed"), file("*.bim"), file("*.fam") into bed_ch

    script:
    bed_prefix = "${chr}.pass.gt_qc"
    """
    export PATH=/home/b/b37974/:$PATH
    plink \
        --vcf ${vcf} \
        --make-bed \
        --keep-allele-order \
        --double-id \
        --out ${bed_prefix}
    """
}

bed_ch
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
    publishDir "${params.outdir}/10.vcf2bed", mode: 'copy'

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

    publishDir "${params.outdir}/10.vcf2bed", mode: 'symlink'

    input:
    file(pmerge_list) from pmerge_lst_ch

    output:
    tuple file("*.bed"), file("*.bim"), file("*.fam") into pmerge_out

    script:
    out_prefix = "cteph_agp3k.autosome.gt_qc"
    """
    export PATH=/home/b/b37974/:$PATH
    plink2 \
        --pmerge-list ${pmerge_list} bfile \
        --threads 16 \
        --make-bed \
        --out ${out_prefix}
    """
}

process RunSampleQC {
    executor 'slurm'
    queue 'gr10478b'
    time '36h'
    tag "RunSampleQC"

    publishDir "${params.outdir}/11.run_sample_qc", mode: 'symlink'

    input:
    tuple file(bed), file(bim), file(fam) from pmerge_out
    val(infoPath) from params.infoPath

    output:
    file('.command.log')
    file('*.pdf')
    tuple file('*.sample_qc_flags.csv'), file('*.sample_qc_summary.csv'), file('*.pi_hat.csv') into sample_qc_out
    tuple file("*.bed"), file("*.bim"), file("*.fam") into update_fam_out

    script:
    bed_prefix = bed.baseName
    info_df = "${infoPath}/cteph_agp3k_jhrpv4.xlsx"
    high_ld = "${infoPath}/high-LD-regions-hg38-GRCh38.txt"
    """
    source activate cteph_geno_pro
    python ${params.scriptDir}/sample_qc_main.py \
        --info_file ${info_df} \
        --bed_prefix ${bed_prefix} \
        --high_ld ${high_ld} \
        --case_prefix PHOM \
        --out_prefix cteph_agp3k \
        --script_path ${params.scriptDir} \
        --pi_threshold 0.2 \
        --dp_robust_z_threshold -3 \
        --het_threshold 5sd \
        --smiss_threshold 0.1 \
        --threads 16
    """
}

process FinishSampleQC {
    executor 'local'
    queue 'gr10478b'
    time '6h'
    tag "FinishSampleQC"

    publishDir "${params.outdir}/11.run_sample_qc", mode: 'symlink'

    input:
    tuple file(bed), file(bim), file(fam) from update_fam_out
    tuple file(sample_qc_flags), file(sample_qc_summary), file(pi_hat) from sample_qc_out

    output:
    file('*.log')
    tuple file("*.bed"), file("*.bim"), file("*.fam") into finish_sample_qc_out

    script:
    bed_prefix = bed.baseName
    """
    source activate cteph_geno_pro
    python ${params.scriptDir}/sample_qc_pipeline.py \
        --sample_qc_flags ${sample_qc_flags} \
        --bed_prefix ${bed_prefix} \
        --out_prefix cteph_agp3k.sqc \
        --mode meandp \
        --include_pass_pi_hat
    """
}   

process rmMAF0orVMISS1 {
    executor 'slurm'
    queue 'gr10478b'
    time '36h'
    tag "rmMAF0orVMISS1"

    publishDir "${params.outdir}/12.rm_maf0_vmiss1", mode: 'symlink'

    input:
    tuple file(bed), file(bim), file(fam) from finish_sample_qc_out

    output:
    file('*.log')
    tuple file("*.bed"), file("*.bim"), file("*.fam") into rm_maf0_vmiss1_out
    tuple file('maf0_or_vmiss1_variants.variant_ids.tsv'), file('maf0_or_vmiss1_variants.with_flags.tsv') into maf0_or_vmiss1_ids_tsv

    script:
    bed_prefix = bed.baseName
    output_prefix = "cteph_agp3k.sqc"
    """
    source activate cteph_geno_pro
    python ${params.scriptDir}/variant_qc_rm_maf0_vmiss1.py \
        --script_path ${params.scriptDir} \
        --bed_prefix ${bed_prefix} \
        --output_prefix ${output_prefix} \
        --threads 32
    """
}

process RunVariantQC {
    executor 'slurm'
    queue 'gr10478b'
    time '36h'
    tag "RunVariantQC"

    publishDir "${params.outdir}/13.run_variant_qc", mode: 'symlink'

    input:
    tuple file(bed), file(bim), file(fam) from rm_maf0_vmiss1_out

    output:
    file('*.pdf')
    file('*.png')
    file('*.log')
    tuple file('vmiss_pass_variants.tsv'), file('hwe_pass_variants.tsv'), file('pass_variants.tsv') into pass_variants_out
    tuple file("*.bed"), file("*.bim"), file("*.fam") into variant_qc_out, variant_qc_out_2

    script:
    bed_prefix = bed.baseName
    """
    source activate cteph_geno_pro
    python ${params.scriptDir}/variant_qc_main.py \
        --script_path ${params.scriptDir} \
        --bed_prefix ${bed_prefix} \
        --output_prefix cteph_agp3k.sqc.vqc \
        --threads 32 \
        --vmiss_threshold 0.04 \
        --hwe_json ${params.scriptDir}/hwe.json
    """
    }

process RunPCA {
    executor 'slurm'
    queue 'gr10478b'
    time '36h'
    tag "RunPCA"

    publishDir "${params.outdir}/14.run_pca", mode: 'symlink'

    input:
    tuple file(bed), file(bim), file(fam) from variant_qc_out
    val(infoPath) from params.infoPath

    output:
    file('*.pdf')
    file('*.log')
    tuple file('*.prune.in'), file('*.prune.out') into pca_prune_out
    tuple file('*.eigenvec'), file('*.eigenval'), file('*.eigenvec.allele') into pca_out
    tuple file('*.no_high_ld.bed'), file('*.no_high_ld.bim'), file('*.no_high_ld.fam') into no_high_ld_out

    script:
    bed_prefix = bed.baseName
    high_ld = "${infoPath}/high-LD-regions-hg38-GRCh38.txt"
    """
    source activate cteph_geno_pro
    python ${params.scriptDir}/pca_qc_main.py \
        --bed_prefix ${bed_prefix} \
        --high_ld ${high_ld} \
        --output_prefix cteph_agp3k \
        --threads 32 \
        --maf_threshold 0.05 \
        --case_prefix PHOM \
        --case_name CTEPH \
        --control_name AGP3K \
        --plink2_path /home/b/b37974/plink2
    """
}

process PrepareBBJ {
    executor 'slurm'
    queue 'gr10478b'
    time '36h'
    tag "PrepareBBJ"

    publishDir "${params.cteph_agp3k_main}/bbj_projection/01.bbj_prepare", mode: 'symlink'

    input:
    val(bbj_main) from params.bbj_main

    output:
    file('*.log')
    tuple file("*.setid.bed"), file("*.setid.bim"), file("*.setid.fam") into bbj_prepare_out

    script:
    bed_prefix = "${bbj_main}/02.sample_qc/NewOE13_Auto.id.b38.sqc"
    """
    source activate cteph_geno_pro
    python ${params.scriptDir}/bbj_prepare_main.py \
        --bbj_bed_prefix ${bed_prefix} \
        --maf_threshold 0.05 \
        --threads 64 \
    """
}

process RunBBJPCA {
    executor 'slurm'
    queue 'gr10478b'
    time '36h'
    tag "RunBBJPCA"

    publishDir "${params.cteph_agp3k_main}/bbj_projection/02.bbj_pca", mode: 'symlink'

    input:
    tuple file(bed), file(bim), file(fam) from bbj_prepare_out
    val(infoPath) from params.infoPath

    output:
    file('*.log')
    file('*.pdf')
    tuple file("*.prune.in"), file("*.prune.out") into bbj_pca_prune_out
    tuple file("*.no_high_ld.bed"), file("*.no_high_ld.bim"), file("*.no_high_ld.fam") into bbj_no_high_ld_out
    tuple file("*.eigenvec"), file("*.eigenval"), file("*.eigenvec.allele"), file("*.acount") into bbj_pca_out
    
    script:
    bed_prefix = bed.baseName
    high_ld = "${infoPath}/high-LD-regions-hg38-GRCh38.txt"
    """
    source activate cteph_geno_pro
    python ${params.scriptDir}/bbj_pca_main.py \
        --bbj_bed_prefix ${bed_prefix} \
        --output_prefix bbj.maf \
        --high_ld ${high_ld} \
        --threads 32
    """
}

process RunBBJProjection {
    executor 'slurm'
    queue 'gr10478b'
    time '36h'
    tag "RunBBJProjection"

    publishDir "${params.cteph_agp3k_main}/bbj_projection/03.bbj_projection", mode: 'symlink'

    input:
    tuple file(my_bed), file(my_bim), file(my_fam) from no_high_ld_out
    tuple file(bbj_bed), file(bbj_bim), file(bbj_fam) from bbj_no_high_ld_out
    tuple file(bbj_prune_in), file(bbj_prune_out) from bbj_pca_prune_out
    tuple file(bbj_eigenvec), file(bbj_eigenval), file(bbj_eigenvec_allele), file(bbj_acount) from bbj_pca_out

    output:
    file('*.log')
    file('*.pdf')
    tuple file('*.sscore'), file('*.sscore.vars') into bbj_projection_out, bbj_projection_out_2

    script:
    my_bed_prefix = my_bed.baseName
    bbj_bed_prefix = bbj_bed.baseName
    """
    source activate cteph_geno_pro
    python ${params.scriptDir}/bbj_projection_main.py \
        --my_bed_prefix ${my_bed_prefix} \
        --bbj_bed_prefix ${bbj_bed_prefix} \
        --bbj_prune_in ${bbj_prune_in} \
        --my_prefix_out cteph_agp3k \
        --bbj_prefix_out bbj \
        --bbj_pca_acount ${bbj_acount} \
        --bbj_pca_eigenvec_allele ${bbj_eigenvec_allele} \
        --threads 32 \
        --case_prefix PHOM \
        --case_name CTEPH \
        --control_name AGP3K
    """
}

process BBJSampleKeep {
    executor 'slurm'
    queue 'gr10478b'
    time '36h'
    tag "BBJSampleKeep"

    publishDir "${params.outdir}/15.bbj_sample_keep", mode: 'symlink'

    input:
    tuple file(bed), file(bim), file(fam) from variant_qc_out_2
    tuple file(sscore), file(sscore_vars) from bbj_projection_out

    output:
    file('*.log')
    file('*.png')
    file('*.txt')
    tuple file("*.bed"), file("*.bim"), file("*.fam") into bbj_sample_keep_out

    script:
    bed_prefix = bed.baseName
    """
    source activate cteph_geno_pro
    python ${params.scriptDir}/bbj_sample_keep_main.py \
        --sscore_path ${sscore} \
        --case_prefix PHOM \
        --case_name CTEPH \
        --control_name AGP3K \
        --prefix_out cteph_agp3k \
        --rect_xlim -0.028 0.016 \
        --rect_ylim -0.024 0.033 \
        --bed_prefix ${bed_prefix} \
        --threads 32
    """
}

process rmMAF0orVMISS1_repeat {
    executor 'slurm'
    queue 'gr10478b'
    time '36h'
    tag "rmMAF0orVMISS1_repeat"

    publishDir "${params.outdir}/16.rm_maf0_vmiss1_repeat", mode: 'symlink'

    input:
    tuple file(bed), file(bim), file(fam) from bbj_sample_keep_out

    output:
    file('*.log')
    tuple file("*.bed"), file("*.bim"), file("*.fam") into rm_maf0_vmiss1_repeat_out, rm_maf0_vmiss1_repeat_out_2
    tuple file('maf0_or_vmiss1_variants.variant_ids.tsv'), file('maf0_or_vmiss1_variants.with_flags.tsv')

    script:
    bed_prefix = bed.baseName
    output_prefix = "cteph_agp3k.sqc.vqc.bbj_sample_keep"
    """
    source activate cteph_geno_pro
    python ${params.scriptDir}/variant_qc_rm_maf0_vmiss1.py \
        --script_path ${params.scriptDir} \
        --bed_prefix ${bed_prefix} \
        --output_prefix ${output_prefix} \
        --threads 32
    """
}

process ToMMoPanelCompare {
    executor 'slurm'
    queue 'gr10478b'
    time '36h'
    tag "ToMMoPanelCompare"

    publishDir "${params.outdir}/17.tommo_panel_compare", mode: 'symlink'

    input:
    tuple file(bed), file(bim), file(fam) from rm_maf0_vmiss1_repeat_out
    val(tommodir) from params.tommodir

    output:
    file('*.pdf')
    file('*.variant_qc_summary.variant_qc_with_tommo.tsv') into tommo_panel_compare_out

    script:
    bed_prefix = bed.baseName
    tommo_vcf_path = "${tommodir}/tommo-60kjpn-20240904-GRCh38-snvindel-af-autosome.norm.vcf.gz"
    """
    source activate cteph_geno_pro
    python ${params.scriptDir}/panel_compare_main.py \
        --bed_prefix ${bed_prefix} \
        --tommo_vcf_path ${tommo_vcf_path} \
        --threads 6 \
        --chunk_size 500000 \
        --max_workers 10 \
        --output_prefix cteph_agp3k \
        --regions_chunk_lines 500000
    """
}

process ToMMoPanelThr {
    executor 'slurm'
    queue 'gr10478b'
    time '36h'
    tag "ToMMoPanelThr"

    publishDir "${params.outdir}/18.tommo_panel_thr", mode: 'symlink'

    input:
    file(variant_qc_with_tommo) from tommo_panel_compare_out

    output:
    file('*.pdf')
    file('*.tsv.gz')
    tuple file('knee_variants.rare.tsv'), file('knee_variants.lowfreq.tsv'), file('knee_variants.common.tsv') into tommo_panel_thr_variant_out
    file('manifest.json') into tommo_panel_thr_out

    script:
    """
    source activate cteph_geno_pro
    python ${params.scriptDir}/panel_thr_main.py \
        --variant_qc_with_tommo ${variant_qc_with_tommo} \
        --chunk_size 500000 \
        --knee_weight_y_map '{"rare": 1.0, "lowfreq": 4.0, "common": 4.0}'
    """
}

process ToMMoPanelFilter {
    executor 'slurm'
    queue 'gr10478b'
    time '36h'
    tag "ToMMoPanelFilter"

    publishDir "${params.outdir}/19.tommo_panel_filter", mode: 'symlink'

    input:
    file(manifest) from tommo_panel_thr_out
    tuple file(bed), file(bim), file(fam) from rm_maf0_vmiss1_repeat_out_2

    output:
    file('*.log')
    file('*.json')
    file('*.tsv')
    file('*.bed')
    file('*.bim')
    file('*.fam') into final_sample_out
    tuple file('*.lowfreq_common.bed'), file('*.lowfreq_common.bim'), file('*.lowfreq_common.fam') into lowfreq_common_out

    script:
    bed_prefix = bed.baseName
    """
    source activate cteph_geno_pro
    python ${params.scriptDir}/panel_filter_main.py \
        --manifest_path ${manifest} \
        --config_json   ${params.scriptDir}/panel_select_config.json \
        --bed_prefix    ${bed_prefix} \
        --out_prefix    cteph_agp3k \
        --threads 6 --chunk_size 500000 --max_workers 10 \
        --keep_tmp \
        --merge_low_common \
        --plink2_path /home/b/b37974/plink2 \
        --save_out_map
    """
}

sample_ch = final_sample_out.map { it[0] }

process CovPhenoPrepare {
    executor 'slurm'
    queue 'gr10478b'
    time '36h'
    tag "CovPhenoPrepare"

    publishDir "${params.outdir}/20.cov_pheno_prepare", mode: 'symlink'

    input:
    file(sample) from sample_ch
    val(infoPath) from params.infoPath
    tuple file(sscore), file(sscore_vars) from bbj_projection_out_2

    output:
    file('*.missing_age_samples.csv')
    tuple file('*.pheno_df.csv'), file('*.cov_df.csv'), file('*.cov_df.no_age.csv') into cov_pheno_out

    script:
    info_df = "${infoPath}/cteph_agp3k_jhrpv4.xlsx"
    """
    source activate cteph_geno_pro
    python ${params.scriptDir}/cov_pheno_prepare_rev1.py \
        --info_path ${info_df} \
        --fam_path ${sample} \
        --bbj_sscore_path ${sscore} \
        --case_prefix PHOM 
    """
}

process MissBiasFilter {
    executor 'slurm'
    queue 'gr10478b'
    time '36h'
    tag "MissBiasFilter"

    publishDir "${params.outdir}/21.miss_bias_filter", mode: 'symlink'

    input:
    tuple file(bed), file(bim), file(fam) from lowfreq_common_out

    output:
    file('*.png')
    file('*.log')
    file('*.missing')
    tuple file('*.bed'), file('*.bim'), file('*.fam') into miss_bias_out

    script:
    bed_prefix = bed.baseName
    """
    source activate cteph_geno_pro
    python ${params.scriptDir}/miss_bias_main.py \
        --bed_prefix ${bed_prefix} \
        --out_prefix_run "cteph_agp3k.missing_bias" \
        --use_midp \
        --color_by_variant \
        --fdr_threshold 0.05 \
        --threads 32 \
        --out_prefix_remove "cteph_agp3k.lowfreq_common.rm_q_lt_0.05"
    """
}