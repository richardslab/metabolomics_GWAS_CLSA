#!/bin/bash
#PBS -o logs/
#PBS -e logs/
#PBS -l walltime=4:00:00
#PBS -l nodes=1:ppn=4
#PBS -l mem=24G
#PBS -l vmem=24G
#PBS -t 1-1091%20

ldpaneldir=/home/richards/yiheng.chen/scratch/project1_2_metabolomics_GWAS_CLSA/codes/data/genotype_data_8299
GWAS_F=/home/richards/yiheng.chen/scratch/project1_2_metabolomics_GWAS_CLSA/codes/GWAS_related_results/metabo_GWAS
output_dir=/home/richards/yiheng.chen/scratch/project1_2_metabolomics_GWAS_CLSA/codes/GWAS_related_results/COJO_results/metabolites

cd ${GWAS_F}
FILE=$(ls *.gz | sed -n ${PBS_ARRAYID}p)

## make directory for each metabolite
mkdir ${output_dir}/${FILE}_clumpped;

for m in {1..22}; do
 /scratch/richards/yiheng.chen/GCTA/gcta_1.93.2beta/gcta64 --bfile ${ldpaneldir}/clsa_imp_${m}_v3_8299 --chr ${m} --maf 0.01 --cojo-p 5e-8 --cojo-wind 5000 --cojo-collinear 0.9 --cojo-slct --thread-num 4 --cojo-file <(zcat ${GWAS_F}/${FILE} | awk 'BEGIN{OFS="\t"; print "SNP","A1","A2","freq","b", "se","p", "N"} NR>1 {print $2,$4,$5,$7,$8,$9,$10,$6}') --out ${output_dir}/${FILE}_clumpped/${FILE}_chr${m};
done;
