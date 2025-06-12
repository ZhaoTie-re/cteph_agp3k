#!/bin/bash
#SBATCH -p gr10478b
#SBATCH -t 48:0:0
#SBATCH --rsc p=1:t=4:c=4:m=16384M
#SBATCH --job-name=re.select.ids

#This script is used to select the after QC filtering snp ids in array vcf file

# Initialize environment
export PATH=/home/b/b37974/:$PATH

array_vcf='/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/array/07.z1.norm.all/cteph_agp3k.ajsa.qc.rechr.norm.vcf.gz'

# using bcftools to extract the select_ids from array vcf vcf
bcftools view -i'ID=@selected_variants.txt' $array_vcf -Oz -o cteph_agp3k.ajsa.qc.rechr.norm.reselect.ids.vcf.gz --threads 4
# index the vcf files
tabix -p vcf cteph_agp3k.ajsa.qc.rechr.norm.reselect.ids.vcf.gz
