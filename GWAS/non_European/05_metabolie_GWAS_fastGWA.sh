#!/bin/bash
#PBS -N clsa.metabolomics.gwas.gcta.v3
#PBS -o logs
#PBS -e logs
#PBS -l walltime=1:00:00
#PBS -l nodes=1:ppn=6
#PBS -l mem=36G
#PBS -l vmem=36G
#PBS -t 1-971%10

pop=south_asian
work_dir=work_dir=/home/richards/yiheng.chen/scratch/project1_2_metabolomics_GWAS_CLSA/codes/non_EUR_GWAS/non_European_ancestry_GWAS_results
program_dir=/scratch/richards/yiheng.chen/project1_2_metabolomics_GWAS_CLSA/codes/program

cd $work_dir/$pop/phenotype/metabolites

FILE=$(ls *.txt | sed -n ${PBS_ARRAYID}p)

${program_dir}/gcta64 --mbfile $work_dir/$pop/genotype/clsa_imp_dat/geno_clsa_v3_bfiles.txt \
--grm-sparse $work_dir/$pop/genotype/grm/sp_grm_${pop}_clsa_ldpruned --fastGWA-mlm --pheno ${FILE} --qcovar $work_dir/$pop/phenotype/qcovarCol_clsa_${pop}_noHeader.txt \
--covar $work_dir/$pop/phenotype/CovarCol_clsa_${pop}_noHeader.txt --threads 6 --out $work_dir/$pop/GWAS/${FILE}_gwas_v3;
