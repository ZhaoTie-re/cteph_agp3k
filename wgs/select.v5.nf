params.JHRPv5 = '/LARGE1/gr10478/platform/JHRPv5/workspace/pipeline/output/annovar.v5/annovar_concat'
params.OutDir = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v5/wgs/results'
params.ScriptDir = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v5/wgs/scripts'
params.SifDir = '/home/b/b37974/simg'
params.SampleList = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v5/info/cteph_agp3k.v5.ls'
params.InfoPath = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v5/info/ph_agp3k_combined_final_20251112.xlsx'
params.HihgLD = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v5/info/high-LD-regions-hg38-GRCh38_modified.txt'
params.NagasakiPipelinePath = '/LARGE0/gr10478/b37974/Pulmonary_Hypertension/nagasaki_pipeline'
params.IdCol = 'ID JHRPv5'
params.PlatformCol = 'WGS'
params.TargetDPCol = 'Target_DP'
params.MeanDPCol = 'DP'    
params.SexCol = 'Sex'
params.GroupCol = 'Outcome'

// Create channel for all chromosomes including X, Y, and PAR
Channel
    .from((1..22).collect { "chr${it}" } + ["chrX", "chrY", "PAR"])
    .set { chr_ch }

process selectPassNormSetID {
    executor 'slurm'
    queue 'gr10478b'
    time '36h'
    tag "${chr}"

    publishDir "${params.OutDir}/01.pass_norm_setid", mode: 'symlink'

    input:
    val(chr) from chr_ch
    path(samplelist) from params.SampleList

    output:
    tuple val(chr), file(output_vcf), file(output_vcf_tbi) into selected_normalized_vcf_out

    script:
    input_vcf = "${params.JHRPv5}/annovar.v5.${chr}.hg38_multianno.vcf.gz"
    output_vcf = "${chr}.selected.pass.norm_split.setid.vcf.gz"
    output_vcf_tbi = "${chr}.selected.pass.norm_split.setid.vcf.gz.tbi"

    """
    # Step 1: Select samples and filter PASS variants
    # Step 2: Filter variants with MAC >= 1 (remove monomorphic sites, both all-ref and all-alt)
    # Step 3: Normalize and split multiallelic variants:
    #         --multiallelics -any: split multiallelic sites into biallelic records
    #         --fasta-ref: use reference genome for left-alignment and normalization
    #         --check-ref s: set/fix variants with REF allele mismatches (e=exit, w=warn, x=exclude, s=set/fix)
    # Step 4: Set variant IDs to CHROM:POS:REF:ALT format
    bcftools view ${input_vcf} --threads 2 -S ${samplelist} --force-samples -Ou | \
        bcftools view --threads 2 -f"PASS" -Ou | \
        bcftools view --threads 2 --min-ac 1:minor -Ou | \
        bcftools norm \
            --multiallelics -any \
            --fasta-ref ${params.NagasakiPipelinePath}/data/hs38DH.fa \
            --check-ref s \
            --threads 2 \
            -Ou | \
        bcftools annotate --set-id '%CHROM:%POS:%REF:%ALT' -Oz -o ${output_vcf}
    
    bcftools index --threads 2 -t ${output_vcf}
    """
}

process filterVariantQuality {
    executor 'slurm'
    queue 'gr10478b'
    time '36h'
    tag "${chr}"

    publishDir "${params.OutDir}/02.vqc_filter", mode: 'symlink'

    input:
    tuple val(chr), file(vcf), file(vcf_tbi) from selected_normalized_vcf_out

    output:
    tuple val(chr), file(filtered_vcf), file(filtered_vcf_tbi) into vqc_filtered_out

    script:
    filtered_vcf = "${chr}.vqc.vcf.gz"
    filtered_vcf_tbi = "${chr}.vqc.vcf.gz.tbi"

    """
    # Filter variants based on quality metrics:
    # VQSLOD > 10: Variant Quality Score Log-Odds (confidence in variant call)
    # MQ > 58.75: Mapping Quality (alignment quality of reads supporting the variant)
    bcftools view ${vcf} --threads 8 -i 'VQSLOD > 10 & MQ > 58.75' -Oz -o ${filtered_vcf}
    bcftools index --threads 8 -t ${filtered_vcf}
    """
}

process addAFAndNormalizeGT {
    executor 'slurm'
    queue 'gr10478b'
    time '36h'
    tag "${chr}"

    publishDir "${params.OutDir}/03.af_gtnorm", mode: 'symlink'

    input:
    tuple val(chr), file(vcf), file(vcf_tbi) from vqc_filtered_out
    val(sifdir) from params.SifDir

    output:
    tuple val(chr), file(output_vcf), file(output_vcf_tbi) into af_gtnorm_out

    script:
    gatk_sif = "${sifdir}/gatk_latest.sif"
    tmp_vcf = "${chr}.vqc.af.tmp.vcf.gz"
    output_vcf = "${chr}.vqc.af.gtnorm.vcf.gz"
    output_vcf_tbi = "${chr}.vqc.af.gtnorm.vcf.gz.tbi"

    """
    # Step 1: Add AlleleFraction annotation using GATK VariantAnnotator
    #         AlleleFraction (AF): Allele balance for heterozygous calls (ratio of ALT allele depth to total depth)
    singularity exec --bind /LARGE0:/LARGE0 ${gatk_sif} gatk VariantAnnotator \
        -R ${params.NagasakiPipelinePath}/data/hs38DH.fa \
        -V ${vcf} \
        -O ${tmp_vcf} \
        -A AlleleFraction \
        --create-output-variant-index true \
        --java-options "-XX:ParallelGCThreads=4"
    
    # Step 2: Normalize genotypes - replace 'nan' with 'NaN' for VCF spec compliance
    # Step 3: Unphase and sort all genotypes with bcftools +setGT
    #         -t a: target all genotypes
    #         -n u: unphase genotype (change | to /) and sort by allele (e.g., 1|0 becomes 0/1)
    # Step 4: Index the final normalized VCF
    bcftools view ${tmp_vcf} --threads 4 | sed 's/nan/NaN/g' | bgzip > ${output_vcf}
    bcftools +setGT ${output_vcf} -- -t a -n u
    bcftools index --threads 4 -t ${output_vcf}
    """
}

process splitAndTagBySex {
    executor 'slurm'
    queue 'gr10478b'
    time '12h'
    tag "${chr}"
    
    publishDir "${params.OutDir}/04.split_by_sex", mode: 'symlink'
    
    input:
    tuple val(chr), file(vcf), file(tbi) from af_gtnorm_out
    path(info_file) from params.InfoPath
    
    output:
    tuple val(chr), file("*.vcf.gz"), file("*.vcf.gz.tbi") into split_sex_raw
    
    script:
    """
    # Setup python script
    python3 ${params.ScriptDir}/get_samples_by_sex.py \
        --excel ${info_file} \
        --id-col "${params.IdCol}" \
        --sex-col "${params.SexCol}" \
        --male-out male_samples.txt \
        --female-out female_samples.txt

    if [[ "${chr}" == "chrX" ]]; then
        # Male
        bcftools view --threads 4 --force-samples -S male_samples.txt ${vcf} -Oz -o ${chr}.1.MaleChrX.vcf.gz
        bcftools index --threads 4 -t ${chr}.1.MaleChrX.vcf.gz
        # Female
        bcftools view --threads 4 --force-samples -S female_samples.txt ${vcf} -Oz -o ${chr}.2.FemaleChrX.vcf.gz
        bcftools index --threads 4 -t ${chr}.2.FemaleChrX.vcf.gz
        
    elif [[ "${chr}" == "chrY" ]]; then
        bcftools view --threads 4 --force-samples -S male_samples.txt ${vcf} -Oz -o ${chr}.MaleChrY.vcf.gz
        bcftools index --threads 4 -t ${chr}.MaleChrY.vcf.gz
    elif [[ "${chr}" == "PAR" ]]; then
        ln -s ${vcf} ${chr}.PseudoautosomalChr.vcf.gz
        ln -s ${tbi} ${chr}.PseudoautosomalChr.vcf.gz.tbi
    else
        # Autosomes
        ln -s ${vcf} ${chr}.AutosomalChr.vcf.gz
        ln -s ${tbi} ${chr}.AutosomalChr.vcf.gz.tbi
    fi
    """
}

split_sex_out = split_sex_raw
    .transpose()
    .map { chr, vcf, tbi ->
        def tag = ""
        def vcfName = vcf.name
        if (vcfName.contains("MaleChrX")) tag = "MaleChrX"
        else if (vcfName.contains("FemaleChrX")) tag = "FemaleChrX"
        else if (vcfName.contains("MaleChrY")) tag = "MaleChrY"
        else if (vcfName.contains("PseudoautosomalChr")) tag = "PseudoautosomalChr"
        else tag = "AutosomalChr"
        
        [tag, chr, vcf, tbi]
    }

process genotypeQC {
    executor 'slurm'
    queue 'gr10478b'
    time '36h'
    tag "${tag}:${chr}"

    publishDir "${params.OutDir}/05.genotype_qc", mode: 'symlink'

    input:
    tuple val(tag), val(chr), file(vcf), file(tbi) from split_sex_out
    val(sifdir) from params.SifDir

    output:
    tuple val(tag), val(chr), file(final_vcf), file(final_tbi) into genotype_qc_out

    script:
    gatk_sif = "${sifdir}/gatk_latest.sif"
    tmp_vcf = "${chr}.${tag}.gt_qc.tmp.vcf.gz"
    final_vcf = "${chr}.${tag}.gt_qc.norm.vcf.gz"
    final_tbi = "${chr}.${tag}.gt_qc.norm.vcf.gz.tbi"
    
    if (tag == 'AutosomalChr' || tag == 'PseudoautosomalChr' || tag == 'FemaleChrX') {
        """
        # Step 1: Normalization (sed nan -> NaN) and re-indexing
        bcftools view ${vcf} --threads 4 | sed 's/nan/NaN/g' | bgzip > ${tmp_vcf}
        bcftools index --threads 4 -t ${tmp_vcf}

        # Step 2: GATK VariantFiltration for Autosomes/PAR/FemaleX
        # DP < 8, GQ < 20, Heterozygous AB outlier
        singularity exec --bind /LARGE0:/LARGE0 ${gatk_sif} gatk --java-options "-Xmx4G" VariantFiltration \
            -R ${params.NagasakiPipelinePath}/data/hs38DH.fa \
            -V ${tmp_vcf} \
            -O ${final_vcf} \
            --genotype-filter-name "LowGQ" \
            --genotype-filter-expression "GQ < 20" \
            --genotype-filter-name "LowDP" \
            --genotype-filter-expression "DP < 8" \
            --genotype-filter-name "ABB_outlier" \
            --genotype-filter-expression "isHet == 1 && (AF < 0.2 || AF > 0.8)" \
            --genotype-filter-name "ABB_NaN" \
            --genotype-filter-expression "AF == 'NaN'" \
            --set-filtered-genotype-to-no-call true \
            --create-output-variant-index true

        # Step 3: Delete tmp_vcf
        rm ${tmp_vcf}*
        """
    } else {
        """
        # Step 1: Normalization (sed nan -> NaN) and re-indexing
        bcftools view ${vcf} --threads 4 | sed 's/nan/NaN/g' | bgzip > ${tmp_vcf}
        bcftools index --threads 4 -t ${tmp_vcf}

        # Step 2: GATK VariantFiltration for MaleChrX/MaleChrY
        # DP < 4
        singularity exec --bind /LARGE0:/LARGE0 ${gatk_sif} gatk --java-options "-Xmx4G" VariantFiltration \
            -R ${params.NagasakiPipelinePath}/data/hs38DH.fa \
            -V ${tmp_vcf} \
            -O ${final_vcf} \
            --genotype-filter-name "LowDP" \
            --genotype-filter-expression "DP < 4" \
            --set-filtered-genotype-to-no-call true \
            --create-output-variant-index true

        # Step 3: Delete tmp_vcf
        rm ${tmp_vcf}*
        """
    }
}

process convertToBed {
    executor 'slurm'
    queue 'gr10478b'
    time '12h'
    tag "${tag}:${chr}"

    publishDir "${params.OutDir}/06.vcf2bed", mode: 'symlink'

    input:
    tuple val(tag), val(chr), file(vcf), file(tbi) from genotype_qc_out

    output:
    tuple val(tag), val(chr), file("*.bed"), file("*.bim"), file("*.fam") into bed_out

    script:
    bed_prefix = "${chr}.${tag}.gt_qc.norm"
    SPLIT_PAR_OPT = tag == 'PseudoautosomalChr' ? '--split-par hg38' : ''
    """
    export PATH=/home/b/b37974/:$PATH
    
    # Generate sample info file for PLINK2 using the python script
    python3 ${params.ScriptDir}/update_sample_info.py \
        --excel ${params.InfoPath} \
        --id-col "${params.IdCol}" \
        --sex-col "${params.SexCol}" \
        --pheno-col "${params.GroupCol}" \
        --out sample_info.txt

    # Run PLINK2 with sex and phenotype update
    # --update-sex: expects FID IID SEX (cols 1, 2, 3)
    # --pheno: loads phenotype from file (cols 1, 2, 4) which is then written to .fam by --make-bed
    # Note: plink2 uses --pheno instead of --update-pheno
    # --split-par hg38: Only added for PseudoautosomalChr to handle PAR regions correctly
    plink2 \
        --vcf ${vcf} \
        --double-id \
        --update-sex sample_info.txt col-num=3 \
        --pheno sample_info.txt \
        --pheno-col-nums 4 \
        ${SPLIT_PAR_OPT} \
        --make-bed \
        --out ${bed_prefix} \
        --threads 4
    """
}

process mergeGenotypes {
    executor 'slurm'
    queue 'gr10478b'
    time '24h'
    tag "${tag}"

    publishDir "${params.OutDir}/07.merged_bed", mode: 'symlink'

    input:
    tuple val(tag), val(chrs), file(beds), file(bims), file(fams) from bed_out.groupTuple()

    output:
    tuple val(tag), file("${output_bed}"), file("${output_bim}"), file("${output_fam}") into merged_bed_out

    script:
    output_bed = "${tag}.merged.bed"
    output_bim = "${tag}.merged.bim"
    output_fam = "${tag}.merged.fam"
    
    """
    export PATH=/home/b/b37974/:$PATH
    
    if [ "${tag}" == "AutosomalChr" ]; then
        rm -f merge_list.txt
        for f in *.bed; do
            echo "\${f%.bed}" >> merge_list.txt
        done
        
        count=\$(wc -l < merge_list.txt)
        
        if [ "\$count" -eq "1" ]; then
            PREFIX=\$(cat merge_list.txt)
            plink2 --bfile \$PREFIX --make-bed --out ${tag}.merged --threads 4
        else
            plink2 \
                --pmerge-list merge_list.txt bfile \
                --make-bed \
                --out ${tag}.merged \
                --threads 4
        fi
    else
        IN_BED=\$(ls *.bed | head -n 1)
        IN_BIM=\$(ls *.bim | head -n 1)
        IN_FAM=\$(ls *.fam | head -n 1)
        
        cp \$IN_BED ${output_bed}
        cp \$IN_BIM ${output_bim}
        cp \$IN_FAM ${output_fam}
    fi
    """
}
