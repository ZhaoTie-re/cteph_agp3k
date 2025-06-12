#usr/bin/zsh
# rm MAV (multi-allelic variants) from plink files

# extract MAV from .log file
log_file='/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/array/05.merge_asa_jsa/cteph_agp3k.ajsa.mind.geno.log'
bed_prefix='/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/array/05.merge_asa_jsa/cteph_agp3k.ajsa.mind.geno'

grep "Warning: Variants" "$log_file" | awk -F"'" '{print $2"\n"$4}' > mav_variants.txt

# remove MAV from plink files
/home/b/b37974/plink2 \
    --bfile "$bed_prefix" \
    --exclude mav_variants.txt \
    --make-bed \
    --out cteph_agp3k.ajsa.mind.geno.rm_mav