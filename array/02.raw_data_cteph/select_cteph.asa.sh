# usr/bin/zsh /home/hdd1/Projects/01.raw_data_cteph/select_cteph.asa.sh
# select CTEPH patients and part of AGP3K from asa

sample_list='/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/array/02.raw_data_cteph/asa.cteph.agp3k_two_columns.ls'
bed_prefix='/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k/array/01.asa.snp.rm/illuminaASA24v1'

# select target samples
/home/b/b37974/plink2 \
    --bfile ${bed_prefix} \
    --keep ${sample_list} \
    --make-bed \
    --out illuminaASA24v1.cteph.agp3k

