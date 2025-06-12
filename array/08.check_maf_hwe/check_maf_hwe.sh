# usr/bin/zsh
# check MAF and HWE for plink files

bed_prefix='/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/array/07.z2.norm.vcf2bed/cteph_agp3k.ajsa.qc.rechr.norm'

# check MAF & HWE for autosomes
/home/b/b37974/plink2 \
    --bfile ${bed_prefix}\
    --freq \
    --hardy \
    --out cteph_agp3k.ajsa.qc.norm.auto_chr

# get cteph & agp3k samples from ${bed_prefix}.fam
# cteph samples iid begin with 'PHOM' and agp3k samples iid begin without 'PHOM'
awk '{if ($2 ~ /^PHOM/) {print $2}}' ${bed_prefix}.fam > cteph_samples.txt
awk '{if ($2 !~ /^PHOM/) {print $2}}' ${bed_prefix}.fam > agp3k_samples.txt

# check MAF & HWE for cteph samples
/home/b/b37974/plink2 \
    --bfile ${bed_prefix}\
    --keep cteph_samples.txt \
    --freq \
    --hardy \
    --out cteph.ajsa.qc.norm.auto_chr

# check MAF & HWE for agp3k samples
/home/b/b37974/plink2 \
    --bfile ${bed_prefix}\
    --keep agp3k_samples.txt \
    --freq \
    --hardy \
    --out agp3k.ajsa.qc.norm.auto_chr
