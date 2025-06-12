# usr/bin/zsh


# extract common snps from illuminaASA24v1.mind.geno and illuminaJSA24v1.mind.geno

asa_prefix='/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/array/04.variant_qc/illuminaASA24v1.mind.geno'
jsa_prefix='/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/array/04.variant_qc/illuminaJSA24v1.mind.geno'

asa_bim="${asa_prefix}.bim"
jsa_bim="${jsa_prefix}.bim"
common_snps="common_snps.txt"

awk 'NR==FNR {snp[$2]; next} $2 in snp {print $2}' "$asa_bim" "$jsa_bim" > "$common_snps"

# using common snps to merge illuminaASA24v1.mind.geno and illuminaJSA24v1.mind.geno based on plink2

/home/b/b37974/plink \
    --bfile "$asa_prefix" \
    --bmerge "$jsa_prefix" \
    --extract "$common_snps" \
    --make-bed \
    --out cteph_agp3k.ajsa.mind.geno