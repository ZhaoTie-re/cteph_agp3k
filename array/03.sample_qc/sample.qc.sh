# usr/bin/zsh /home/hdd1/Projects/02.sample_qc/sample.qc.sh
# sample QC for illuminaASA24v1.cteph.agp3k & illuminaASA24v1.cteph.agp3k separately

asa_prefix='/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/array/02.raw_data_cteph/illuminaASA24v1.cteph.agp3k'
jsa_prefix='/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/array/02.raw_data_cteph/illuminaJSA24v1.cteph.agp3k'

# sample QC for asa
/home/b/b37974/plink2 \
    --bfile ${asa_prefix} \
    --mind 0.01 \
    --make-bed \
    --out illuminaASA24v1.mind

/home/b/b37974/plink2 \
    --bfile ${jsa_prefix} \
    --mind 0.01 \
    --make-bed \
    --out illuminaJSA24v1.mind


