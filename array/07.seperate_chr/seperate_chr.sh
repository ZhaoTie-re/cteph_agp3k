#usr/bin/zsh 
# seperate plink files by autosome, sex chromosomes, and mitochondrial

bed_prefix='/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/array/06.rm_mav/cteph_agp3k.ajsa.mind.geno.rm_mav'

# Separate autosomes
/home/b/b37974/plink2 --bfile ${bed_prefix} --chr 1-22 --make-bed --out cteph_agp3k.ajsa.mind.geno.rm_mav.auto_chr

# Separate sex chromosomes
/home/b/b37974/plink2 --bfile ${bed_prefix} --chr X,Y,XY --make-bed --out  cteph_agp3k.ajsa.mind.geno.rm_mav.sex_chr

# Separate mitochondrial chromosome
/home/b/b37974/plink2 --bfile ${bed_prefix} --chr MT --make-bed --out cteph_agp3k.ajsa.mind.geno.rm_mav.mito
