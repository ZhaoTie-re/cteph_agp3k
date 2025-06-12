#usr/bin/zsh 
# convert bed to vcf

# This bed2vcf.sh is used for all variants (before filtering)
# This script converts a bed file to a vcf file, renames the chromosomes, and normalizes the vcf file.
# It is used to prepare the data for further analysis.

# Initialize environment
export PATH=/home/b/b37974/:$PATH

bed_prefix='/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/array/07.seperate_chr/cteph_agp3k.ajsa.mind.geno.rm_mav.auto_chr'

plink \
    --bfile ${bed_prefix} \
    --recode vcf-iid bgz \
    --out cteph_agp3k.ajsa.mind.geno.rm_mav.auto_chr

# index the vcf file
tabix -p vcf cteph_agp3k.ajsa.mind.geno.rm_mav.auto_chr.vcf.gz

# rename the chr names in the vcf file
bcftools annotate --rename-chrs /LARGE0/gr10478/b37974/Pulmonary_Hypertension/BBJ_genome_b38/03.maf_norm/rename_chr.txt \
    cteph_agp3k.ajsa.mind.geno.rm_mav.auto_chr.vcf.gz \
    -Oz \
    -o cteph_agp3k.ajsa.mind.geno.rm_mav.auto_chr.rechr.vcf.gz

# index the renamed vcf file
bcftools index -t cteph_agp3k.ajsa.mind.geno.rm_mav.auto_chr.rechr.vcf.gz

# remove the old vcf file
rm cteph_agp3k.ajsa.mind.geno.rm_mav.auto_chr.vcf.gz
rm cteph_agp3k.ajsa.mind.geno.rm_mav.auto_chr.vcf.gz.tbi

# normalize the vcf file
bcftools norm \
    --multiallelics -any \
    --fasta-ref /LARGE0/gr10478/b37974/Pulmonary_Hypertension/nagasaki_pipeline/data/hs38DH.fa \
    --check-ref s cteph_agp3k.ajsa.mind.geno.rm_mav.auto_chr.rechr.vcf.gz \
    -Oz \
    -o cteph_agp3k.ajsa.qc.rechr.norm.vcf.gz
# index the normalized vcf file
bcftools index -t cteph_agp3k.ajsa.qc.rechr.norm.vcf.gz

# remove the old vcf file
rm cteph_agp3k.ajsa.mind.geno.rm_mav.auto_chr.rechr.vcf.gz
rm cteph_agp3k.ajsa.mind.geno.rm_mav.auto_chr.rechr.vcf.gz.tbi