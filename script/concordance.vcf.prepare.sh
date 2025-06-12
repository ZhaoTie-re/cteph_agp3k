# usr/bin/zsh 
# prepare data for concordance

array_vcf=$1
wgs_vcf=$2
PICARD_SIG="/home/b/b37974/simg/picard.sif"

# get the common variants (ID) between the two vcf files (list) using bcftools
bcftools query -f '%ID\n' $array_vcf | sort > array_ids.txt
bcftools query -f '%ID\n' $wgs_vcf | sort > wgs_ids.txt
comm -12 array_ids.txt wgs_ids.txt > common_ids.txt

# get common samples
bcftools query -l $array_vcf | sort > array_samples.txt
bcftools query -l $wgs_vcf | sort > wgs_samples.txt
comm -12 array_samples.txt wgs_samples.txt > common_samples.txt

# extract common variants AND common samples
bcftools view -i'ID=@common_ids.txt' -S common_samples.txt $array_vcf -Oz -o array_common.vcf.gz --threads 4
bcftools view -i'ID=@common_ids.txt' -S common_samples.txt $wgs_vcf -Oz -o wgs_common.vcf.gz --threads 4

# index the vcf files
tabix -p vcf array_common.vcf.gz
tabix -p vcf wgs_common.vcf.gz

# create a reference genome dictionary for the vcf files

singularity exec -B /LARGE0/gr10478/b37974:/LARGE0/gr10478/b37974 ${PICARD_SIG} java -XX:ParallelGCThreads=2 -jar /home/b/b37974/anaconda3/envs/picard_env/share/picard-3.3.0-0/picard.jar CreateSequenceDictionary \
    --REFERENCE /LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_array_gwas/data_src/nagasaki_pipeline/hs38DH.fa \
    --OUTPUT reference.dict

# update the VCF sequence dictionary with the reference genome dictionary
singularity exec -B /LARGE0/gr10478/b37974:/LARGE0/gr10478/b37974 ${PICARD_SIG} java -XX:ParallelGCThreads=2 -jar /home/b/b37974/anaconda3/envs/picard_env/share/picard-3.3.0-0/picard.jar UpdateVcfSequenceDictionary \
    --INPUT array_common.vcf.gz \
    --OUTPUT array_common_reheadered.vcf.gz \
    --SEQUENCE_DICTIONARY reference.dict

singularity exec -B /LARGE0/gr10478/b37974:/LARGE0/gr10478/b37974 ${PICARD_SIG} java -XX:ParallelGCThreads=2 -jar /home/b/b37974/anaconda3/envs/picard_env/share/picard-3.3.0-0/picard.jar UpdateVcfSequenceDictionary \
    --INPUT wgs_common.vcf.gz \
    --OUTPUT wgs_common_reheadered.vcf.gz \
    --SEQUENCE_DICTIONARY reference.dict

# index the reheadered vcf files
tabix -p vcf array_common_reheadered.vcf.gz
tabix -p vcf wgs_common_reheadered.vcf.gz