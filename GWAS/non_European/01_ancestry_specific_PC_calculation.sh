## this is an example of codes for South Asian sub-population. We performed the same analyses on African and East Asian sub-populations.
#!/bin/bash
#PBS -N clsa_pca
#PBS -o ./logs/
#PBS -e ./logs/
#PBS -l walltime=1:00:00
#PBS -l nodes=1:ppn=16
#PBS -l mem=24G
#PBS -l vmem=24G

pop=south_asian
work_dir=/home/richards/yiheng.chen/scratch/project1_2_metabolomics_GWAS_CLSA/codes/non_EUR_GWAS/non_European_ancestry_GWAS_results

###### extract ancestry-specific genotype data
cd /scratch/richards/yiheng.chen/project1_2_metabolomics_GWAS_CLSA/codes

plink --bfile ./data/clsa_imputed_genotype_data/clsa_gen_v3 --chr 1-22 \
--keep $work_dir/$pop/phenotype/sqc_file_${pop}_ID.txt \
--make-bed \
--out $work_dir/$pop/genotype/clsa_gen_v3_${pop}_nosex

## check the related individuals
king -b $work_dir/$pop/genotype/clsa_gen_v3_${pop}_nosex.bed --unrelated --degree 2
## no one needs to be removed within each ancestry

###### snp filter for PC calculation
plink --bfile $work_dir/$pop/genotype/clsa_gen_v3_${pop}_nosex --chr 1-22 \
--exclude ./data/high-LD-regions-hg19.txt \
--snps-only just-acgt \
--geno 0.05 \
--hwe 0.000001 \
--maf 0.05 \
--keep-allele-order \
--make-bed \
--out $work_dir/$pop/genotype/clsa_gen_v3_${pop}_nosex_filtered

#identify palindromic snps
awk '($5~/[AT]/ && $6 ~/[AT]/) || ($5~/[GC]/ && $6~/[GC]/) { print $2 }' $work_dir/$pop/genotype/clsa_gen_v3_${pop}_nosex_filtered.bim > $work_dir/$pop/genotype/clsa_gen_v3_${pop}_nosex_filtered.palin;

####### get  independent snps from non-sex chromosome 
plink --bfile $work_dir/$pop/genotype/clsa_gen_v3_${pop}_nosex_filtered \
--indep-pairwise 1000 kb 5 0.1 \
--out $work_dir/$pop/genotype/clsa_gen_v3_${pop}_nosex_filtered_ld_checked

## create genotype data for specific ancestry
plink --bfile $work_dir/$pop/genotype/clsa_gen_v3_${pop}_nosex_filtered \
--exclude $work_dir/$pop/genotype/clsa_gen_v3_${pop}_nosex_filtered.palin \
--extract $work_dir/$pop/genotype/clsa_gen_v3_${pop}_nosex_filtered_ld_checked.prune.in \
--make-bed \
--out $work_dir/$pop/genotype/clsa_gen_v3_pruned_forPC_${pop}

### run PC calcultion
plink --bfile $work_dir/$pop/genotype/clsa_gen_v3_pruned_forPC_${pop} \
--pca \
--out $work_dir/$pop/genotype/clsa_gen_v3_pruned_forPC_${pop}_pca
