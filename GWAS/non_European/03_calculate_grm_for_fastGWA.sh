#!/bin/bash
#PBS -N clsa.bgen.to.bfiles
#PBS -o logs/
#PBS -e logs/
#PBS -l walltime=1:00:00
#PBS -l nodes=1:ppn=4
#PBS -l mem=24G
#PBS -l vmem=24G

pop=south_asian
work_dir=/home/richards/yiheng.chen/scratch/project1_2_metabolomics_GWAS_CLSA/codes/non_EUR_GWAS/non_European_ancestry_GWAS_results

########## prune the list of SNPs used for calculating grm
cd /scratch/richards/yiheng.chen/project1_2_metabolomics_GWAS_CLSA/codes

plink \
--bfile $work_dir/$pop/genotype/clsa_gen_v3_${pop}_nosex \
--exclude ./data/high-LD-regions-hg19.txt \
--geno 0.05 \
--hwe 0.000001 \
--indep-pairwise 50 5 0.6 \
--keep-allele-order \
--maf 0.001 \
--out $work_dir/$pop/genotype/grm/clsa_gen_v3_nosexchr_${pop}.pruned \
--write-snplist

################# calculate grm
gcta64 \
--bfile $work_dir/$pop/genotype/clsa_gen_v3_${pop}_nosex \
--autosome \
--make-grm \
--exclude $work_dir/$pop/genotype/grm/clsa_gen_v3_nosexchr_${pop}.pruned.prune.out \
--thread-num 4 \
--out $work_dir/$pop/genotype/grm/grm_${pop}_clsa_ldpruned

######### makse sparse grm 
gcta64 --grm $work_dir/$pop/genotype/grm/grm_${pop}_clsa_ldpruned --make-bK-sparse 0.05 --out $work_dir/$pop/genotype/grm/sp_grm_${pop}_clsa_ldpruned
