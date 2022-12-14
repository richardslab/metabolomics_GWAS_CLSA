#!/bin/bash
#PBS -N clsa.metabolomics.gwas.gcta.v3
#PBS -o logs/
#PBS -e logs/
#PBS -l walltime=5:00:00
#PBS -l nodes=1:ppn=4
#PBS -l mem=24G
#PBS -l vmem=24G
#PBS -t 1-309%10

cd ${analysis_dir}/GWAS_related_results/metabolites_ratio_phenotype

file_dir=${analysis_dir}/data/europ_grm
output_dir=${analysis_dir}/GWAS_related_results/metabo_ratio_GWAS

FILE=$(ls *.txt | sed -n ${PBS_ARRAYID}p)

gcta64 --mbfile ${file_dir}/geno_clsa_v3_bfiles.txt --grm-sparse ${file_dir}/sp_grm_european_clsa_v3_ldpruned --fastGWA-mlm --pheno ${FILE} --qcovar ${file_dir}/qcovarCol_clsa_europ_unrelated_noHeader.txt --covar ${file_dir}/CovarCol_clsa_europ_unrelated_noHeader.txt --threads 4 --out ${output_dir}/${FILE}_gwas_v3;

gzip ${output_dir}/${FILE}_gwas_v3.fastGWA
