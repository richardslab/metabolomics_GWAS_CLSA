## the following codes filter all imputed SNPs with two critiera: maf=0.001; imputation info score: >0.3
cd /scratch/richards/yiheng.chen/project1_2_metabolomics_GWAS_CLSA/codes

for i in {1..22}; {
        plink2 --bgen ./data/clsa_imputed_genotype_data/clsa_imp_${i}_v3.bgen ref-first \
        --sample ./data/clsa_imputed_genotype_data/clsa_imp_v3.sample --maf 0.001 --mach-r2-filter 0.3 2.0 --make-bed \
        --out ./data/clsa_imputed_genotype_processed/clsa_imp_${i}_v3
}

##### generate a list of unrelated individuals use king
## merge plink files for all autosomes
plink --bfile ./data/clsa_imputed_genotype_data/clsa_gen_v3 \
--chr 1-22 --make-bed --out ./data/clsa_genotype_processed/clsa_gen_v3_nosexchr
## use KING to generated a list of unrelated individuals
king -b ./data/clsa_genotype_processed/clsa_gen_v3_nosexchr.bed --unrelated --degree 2

## get  genotying snps from non-sex chromosome 

plink --bfile ./data/clsa_imputed_genotype_data/clsa_gen_v3 --chr 1-22 --make-bed \
--out ./data/clsa_genotype_processed_grm/clsa_gen_v3_nosexchr

########## prune the list of SNPs used for calculating grm
plink \
--bfile ./data/clsa_genotype_processed_grm/clsa_gen_v3_nosexchr \
--exclude ./data/high-LD-regions-hg19.txt \
--geno 0.05 \
--hwe 0.000001 \
--indep-pairwise 50 5 0.6 \
--keep ./data/clsa_v3_europ_id_keep \
--keep-allele-order \
--maf 0.001 \
--out ./data/LD_pruned_snp_for_grm/clsa_gen_v3_nosexchr.europ.pruned \
--write-snplist

################# calculate grm
gcta64 \
--bfile ./data/clsa_genotype_processed_grm/clsa_gen_v3_nosexchr \
--autosome \
--make-grm \
--exclude ./data/LD_pruned_snp_for_grm/clsa_gen_v3_nosexchr.europ.pruned.prune.out \
--thread-num 4 \
--out ./data/europ_grm/grm_european_clsa_ldpruned

######### makse sparse grm 
gcta64 --grm ./data/europ_grm/grm_european_clsa_ldpruned --make-bK-sparse 0.05 --out ./data/europ_grm/sp_grm_european_clsa_v3_ldpruned
