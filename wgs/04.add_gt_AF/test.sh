# usr/bin/zsh

export PATH=/home/b/b37974/:$PATH

vcf='/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/wgs/04.add_gt_AF/chr22.pass.mac1.setid.gt_af.vcf.gz'

plink \
    --vcf ${vcf} \
    --make-bed \
    --keep-allele-order \
    --double-id \
    --out test