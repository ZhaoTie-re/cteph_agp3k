# usr/bin/zsh
# filter plink files by MAF and HWE

bed_prefix='/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/array/07.seperate_chr/cteph_agp3k.ajsa.mind.geno.rm_mav.auto_chr'

# filter by MAF > 0.01 and HWE p-value > 1e-6
/home/b/b37974/plink2 \
    --bfile ${bed_prefix} \
    --maf 0.01 \
    --hwe 1e-6 \
    --make-bed \
    --out cteph_agp3k.ajsa.mind.geno.rm_mav.auto_chr.maf_hwe_qc