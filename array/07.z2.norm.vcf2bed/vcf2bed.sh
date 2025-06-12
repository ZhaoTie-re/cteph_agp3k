# usr/bin/zsh
# convert bed to vcf

# This bed2vcf.sh is used for all variants (after filtering)
# This script converts vcf.gz file to plink format (bed/bim/fam) file, keeps the original ref and alt alleles.

# Initialize environment
export PATH=/home/b/b37974/:$PATH

# Define the input and output file paths
input_vcf='/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/array/07.z1.norm.all/cteph_agp3k.ajsa.qc.rechr.norm.vcf.gz'

output_prefix='cteph_agp3k.ajsa.qc.rechr.norm'

# Convert the VCF file to PLINK format (bed/bim/fam)
plink2 \
    --vcf "${input_vcf}" \
    --make-bed \
    --out "${output_prefix}" \
    --keep-allele-order \
    --threads 16