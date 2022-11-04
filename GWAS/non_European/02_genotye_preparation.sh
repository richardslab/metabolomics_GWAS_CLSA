#!/bin/bash
#PBS -N clsa.bgen.to.bfiles
#PBS -o logs/
#PBS -e logs/
#PBS -l walltime=6:00:00
#PBS -l nodes=1:ppn=6
#PBS -l mem=36G
#PBS -l vmem=36G
## filter cutoff: maf=0.05; imputation info score: >0.3

pop=south_asian
work_dir=work_dir=/home/richards/yiheng.chen/scratch/project1_2_metabolomics_GWAS_CLSA/codes/non_EUR_GWAS/non_European_ancestry_GWAS_results
cd /scratch/richards/yiheng.chen/project1_2_metabolomics_GWAS_CLSA/codes

for i in {1..22}; {
        plink2 --bgen ./data/clsa_imputed_genotype_data/clsa_imp_${i}_v3.bgen ref-first \
        --sample ./data/clsa_imputed_genotype_data/clsa_imp_v3.sample --maf 0.05 --mach-r2-filter 0.3 2.0 \
        --keep $work_dir/${pop}/phenotype/sqc_file_${pop}_ID.txt \
        --make-bed \
        --out $work_dir/${pop}/genotype/clsa_imp_dat/clsa_imp_${pop}_${i}_v3
}
