# usr/bin/zsh

export PATH=/home/b/b37974/:$PATH

vcf='/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/wgs/02.mac1/chr22.pass.mac1.vcf.gz'

plink \
    --vcf ${vcf} \
    --make-bed \
    --keep-allele-order \
    --double-id \
    --out test