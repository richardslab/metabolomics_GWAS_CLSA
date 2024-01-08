##### the following script extract eQTL associations for mQTLs
library(data.table)
library(dplyr)
library(stringr)
library(tidyr)
library(biomaRt)

clsa_metabolite_cojo_SNPs_filtered<-read.csv("./results/clsa_metabolites_cojo_formatted_rsid_formatted_filtered.csv",row.names = 1)
clsa_metabolite_ratios_cojo_SNPs_formatted_rsid_filtered<-read.csv("./results/clsa_metabolite_ratios_cojo_SNPs_formatted_rsid_bonf_filtered_enzyme_transporter.csv",row.names = 1)
clsa_metabolite_and_ratio_cojo_SNPs_filtered<-rbind(clsa_metabolite_cojo_SNPs_filtered[,c(1:11,18)], clsa_metabolite_ratios_cojo_SNPs_formatted_rsid_filtered[,c(1:11,18)])

#### get gene annotations
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
ensembl = useMart(biomart="ENSEMBL_MART_ENSEMBL",
                  dataset="hsapiens_gene_ensembl", 
                  host="uswest.ensembl.org",
                  ensemblRedirect = FALSE)

############################## eQTL  ############################## 
################ whole blood eqtl ################ 
Whole_Blood.v8.signif_variant_gene_pairs<-fread("./data/GTEX_eqtl_database/eqtls_EUR/Whole_Blood.v8.EUR.signif_pairs.txt.gz")
clsa_metabolite_and_ratio_cojo_SNPs_filtered<-clsa_metabolite_and_ratio_cojo_SNPs_filtered %>% rowwise() %>%
  mutate(variant_id_1 = paste0(str_split(SNP_ID,":",4)[[1]][1],"_",str_split(SNP_ID,":",4)[[1]][2],"_",str_split(SNP_ID,":",4)[[1]][3],"_",str_split(SNP_ID,":",4)[[1]][4],"_b38"),
         variant_id_2 = paste0(str_split(SNP_ID,":",4)[[1]][1],"_",str_split(SNP_ID,":",4)[[1]][2],"_",str_split(SNP_ID,":",4)[[1]][4],"_",str_split(SNP_ID,":",4)[[1]][3],"_b38"))

clsa_metabolite_and_ratio_cojo_SNPs_filtered_whole_blood_eqtl<-left_join(clsa_metabolite_and_ratio_cojo_SNPs_filtered,Whole_Blood.v8.signif_variant_gene_pairs, by=c("variant_id_1"="variant_id"))
clsa_metabolite_and_ratio_cojo_SNPs_filtered_whole_blood_eqtl_no_match<-clsa_metabolite_and_ratio_cojo_SNPs_filtered_whole_blood_eqtl %>% filter(is.na(phenotype_id)==T)
clsa_metabolite_and_ratio_cojo_SNPs_filtered_whole_blood_eqtl_no_match1<-left_join(clsa_metabolite_and_ratio_cojo_SNPs_filtered_whole_blood_eqtl_no_match,Whole_Blood.v8.signif_variant_gene_pairs, by=c("variant_id_2"="variant_id"))## no matching found with variant_id_2
clsa_metabolite_and_ratio_cojo_SNPs_filtered_whole_blood_eqtl_filtered<-clsa_metabolite_and_ratio_cojo_SNPs_filtered_whole_blood_eqtl %>% filter(is.na(phenotype_id)!=T)

### get the gene symbol for the eqtl
clsa_metabolite_and_ratio_cojo_SNPs_filtered_whole_blood_eqtl_filtered$gene_id_2 <- sub("[.][0-9]*","",clsa_metabolite_and_ratio_cojo_SNPs_filtered_whole_blood_eqtl_filtered$phenotype_id)
genes <- clsa_metabolite_and_ratio_cojo_SNPs_filtered_whole_blood_eqtl_filtered$gene_id_2
gene_IDs <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),
                  values = genes, mart= mart)
clsa_metabolite_and_ratio_cojo_SNPs_filtered_whole_blood_eqtl_filtered$gene_symbol<-apply(clsa_metabolite_and_ratio_cojo_SNPs_filtered_whole_blood_eqtl_filtered,1,function(x) gene_IDs$hgnc_symbol[gene_IDs$ensembl_gene_id == x[26]][1])
write.csv(clsa_metabolite_and_ratio_cojo_SNPs_filtered_whole_blood_eqtl_filtered,"./analysis_dir/expression_relevance_check/eqtl/clsa_metabolite_and_ratio_cojo_SNPs_whole_blood_eqtl_filtered.csv")

################ tissue-specific eqtl ################
eqtl_list<-c("Adipose_Subcutaneous.v8.EUR.signif_pairs.txt.gz", "Brain_Spinal_cord_cervical_c-1.v8.EUR.signif_pairs.txt.gz","Nerve_Tibial.v8.EUR.signif_pairs.txt.gz",
             "Adipose_Visceral_Omentum.v8.EUR.signif_pairs.txt.gz", "Brain_Substantia_nigra.v8.EUR.signif_pairs.txt.gz", "Ovary.v8.EUR.signif_pairs.txt.gz",
             "Adrenal_Gland.v8.EUR.signif_pairs.txt.gz", "Breast_Mammary_Tissue.v8.EUR.signif_pairs.txt.gz", "Pancreas.v8.EUR.signif_pairs.txt.gz",
             "Artery_Aorta.v8.EUR.signif_pairs.txt.gz", "Cells_Cultured_fibroblasts.v8.EUR.signif_pairs.txt.gz", "Pituitary.v8.EUR.signif_pairs.txt.gz",
             "Artery_Coronary.v8.EUR.signif_pairs.txt.gz", "Cells_EBV-transformed_lymphocytes.v8.EUR.signif_pairs.txt.gz", "Prostate.v8.EUR.signif_pairs.txt.gz",
             "Artery_Tibial.v8.EUR.signif_pairs.txt.gz", "Colon_Sigmoid.v8.EUR.signif_pairs.txt.gz", "Skin_Not_Sun_Exposed_Suprapubic.v8.EUR.signif_pairs.txt.gz",
             "Brain_Amygdala.v8.EUR.signif_pairs.txt.gz", "Colon_Transverse.v8.EUR.signif_pairs.txt.gz", "Skin_Sun_Exposed_Lower_leg.v8.EUR.signif_pairs.txt.gz",
             "Brain_Anterior_cingulate_cortex_BA24.v8.EUR.signif_pairs.txt.gz","Esophagus_Gastroesophageal_Junction.v8.EUR.signif_pairs.txt.gz", "Small_Intestine_Terminal_Ileum.v8.EUR.signif_pairs.txt.gz",
             "Brain_Caudate_basal_ganglia.v8.EUR.signif_pairs.txt.gz", "Esophagus_Mucosa.v8.EUR.signif_pairs.txt.gz","Spleen.v8.EUR.signif_pairs.txt.gz",
             "Brain_Cerebellar_Hemisphere.v8.EUR.signif_pairs.txt.gz","Esophagus_Muscularis.v8.EUR.signif_pairs.txt.gz", "Stomach.v8.EUR.signif_pairs.txt.gz",
             "Brain_Cerebellum.v8.EUR.signif_pairs.txt.gz","Heart_Atrial_Appendage.v8.EUR.signif_pairs.txt.gz","Testis.v8.EUR.signif_pairs.txt.gz",
             "Brain_Cortex.v8.EUR.signif_pairs.txt.gz", "Heart_Left_Ventricle.v8.EUR.signif_pairs.txt.gz", "Thyroid.v8.EUR.signif_pairs.txt.gz",
             "Brain_Frontal_Cortex_BA9.v8.EUR.signif_pairs.txt.gz", "Kidney_Cortex.v8.EUR.signif_pairs.txt.gz", "Uterus.v8.EUR.signif_pairs.txt.gz",
             "Brain_Hippocampus.v8.EUR.signif_pairs.txt.gz", "Liver.v8.EUR.signif_pairs.txt.gz", "Vagina.v8.EUR.signif_pairs.txt.gz",
             "Brain_Hypothalamus.v8.EUR.signif_pairs.txt.gz", "Lung.v8.EUR.signif_pairs.txt.gz",
             "Brain_Nucleus_accumbens_basal_ganglia.v8.EUR.signif_pairs.txt.gz",  "Minor_Salivary_Gland.v8.EUR.signif_pairs.txt.gz",
             "Brain_Putamen_basal_ganglia.v8.EUR.signif_pairs.txt.gz", "Muscle_Skeletal.v8.EUR.signif_pairs.txt.gz")

for (i in 1:length(eqtl_list)) {
  tissue<-as.character(str_split(eqtl_list[i],".v8")[[1]][1])
  clsa_metabolite_and_ratio_cojo_SNPs_filtered_1<-clsa_metabolite_and_ratio_cojo_SNPs_filtered %>% rowwise() %>%mutate(variant_id_1 = paste0(str_split(SNP_ID,":",4)[[1]][1],"_",str_split(SNP_ID,":",4)[[1]][2],"_",str_split(SNP_ID,":",4)[[1]][3],"_",str_split(SNP_ID,":",4)[[1]][4],"_b38"),variant_id_2 = paste0(str_split(SNP_ID,":",4)[[1]][1],"_",str_split(SNP_ID,":",4)[[1]][2],"_",str_split(SNP_ID,":",4)[[1]][4],"_",str_split(SNP_ID,":",4)[[1]][3],"_b38"))
  tissue_specific_signif_variant_gene_pairs<-fread(paste0("./data/GTEX_eqtl_database/eqtls_EUR/",eqtl_list[i]))
  clsa_metabolite_and_ratio_cojo_SNPs_filtered_etql<-left_join(clsa_metabolite_and_ratio_cojo_SNPs_filtered_1,tissue_specific_signif_variant_gene_pairs, by=c("variant_id_1"="variant_id"))
  clsa_metabolite_and_ratio_cojo_SNPs_filtered_etql_no_match<-clsa_metabolite_and_ratio_cojo_SNPs_filtered_etql %>% filter(is.na(phenotype_id)==T)
  clsa_metabolite_and_ratio_cojo_SNPs_filtered_etql_no_match_checked<-left_join(clsa_metabolite_and_ratio_cojo_SNPs_filtered_1,tissue_specific_signif_variant_gene_pairs, by=c("variant_id_2"="variant_id"))## no matching found with variant_id_2
  ## checked, variant_id_1 is enough for matching to variant_id in eqtl database
  clsa_metabolite_and_ratio_cojo_SNPs_filtered_etql_no_match_check<-clsa_metabolite_and_ratio_cojo_SNPs_filtered_etql_no_match_checked%>%filter(is.na(phenotype_id)!=T)
  if (nrow(clsa_metabolite_and_ratio_cojo_SNPs_filtered_etql_no_match_check)==0) {
    clsa_metabolite_and_ratio_cojo_SNPs_filtered_etql_filtered<-clsa_metabolite_and_ratio_cojo_SNPs_filtered_etql %>% filter(is.na(phenotype_id)!=T)
    clsa_metabolite_and_ratio_cojo_SNPs_filtered_etql_filtered$gene_id_2 <-sub("[.][0-9]*","",clsa_metabolite_and_ratio_cojo_SNPs_filtered_etql_filtered$phenotype_id)
    genes <- clsa_metabolite_and_ratio_cojo_SNPs_filtered_etql_filtered$gene_id_2
    gene_IDs <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),values = genes, mart= mart)
    clsa_metabolite_and_ratio_cojo_SNPs_filtered_etql_filtered$gene_symbol<-apply(clsa_metabolite_and_ratio_cojo_SNPs_filtered_etql_filtered,1,function(x) gene_IDs$hgnc_symbol[gene_IDs$ensembl_gene_id == x[26]][1])
    write.csv(clsa_metabolite_and_ratio_cojo_SNPs_filtered_etql_filtered,paste0("./analysis_dir/expression_relevance_check/eqtl/clsa_metabolite_and_ratio_cojo_SNPs_",tissue,"_eqtl_filtered.csv"))
  }
  else {clsa_metabolite_and_ratio_cojo_SNPs_filtered_etql<-rbind(clsa_metabolite_and_ratio_cojo_SNPs_filtered_etql,clsa_metabolite_and_ratio_cojo_SNPs_filtered_etql_no_match_check)
  clsa_metabolite_and_ratio_cojo_SNPs_filtered_etql_filtered<-clsa_metabolite_and_ratio_cojo_SNPs_filtered_etql %>% filter(is.na(phenotype_id)!=T)
  clsa_metabolite_and_ratio_cojo_SNPs_filtered_etql_filtered$gene_id_2 <- sub("[.][0-9]*","",clsa_metabolite_and_ratio_cojo_SNPs_filtered_etql_filtered$phenotype_id)
  genes <- clsa_metabolite_and_ratio_cojo_SNPs_filtered_etql_filtered$gene_id_2
  gene_IDs <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),values = genes, mart= mart)
  clsa_metabolite_and_ratio_cojo_SNPs_filtered_etql_filtered$gene_symbol<-apply(clsa_metabolite_and_ratio_cojo_SNPs_filtered_etql_filtered,1,function(x) gene_IDs$hgnc_symbol[gene_IDs$ensembl_gene_id == x[26]][1])
  write.csv(clsa_metabolite_and_ratio_cojo_SNPs_filtered_etql_filtered,paste0("./analysis_dir/expression_relevance_check/eqtl/clsa_metabolite_and_ratio_cojo_SNPs_",tissue,"_eqtl_filtered.csv"))
  }
}

###### combine all tissue-specific eqtl to one file
file_list<-c("clsa_metabolite_and_ratio_cojo_SNPs_Adipose_Subcutaneous_eqtl_filtered.csv", "clsa_metabolite_and_ratio_cojo_SNPs_Esophagus_Mucosa_eqtl_filtered.csv",
             "clsa_metabolite_and_ratio_cojo_SNPs_Adipose_Visceral_Omentum_eqtl_filtered.csv", "clsa_metabolite_and_ratio_cojo_SNPs_Esophagus_Muscularis_eqtl_filtered.csv",
             "clsa_metabolite_and_ratio_cojo_SNPs_Adrenal_Gland_eqtl_filtered.csv", "clsa_metabolite_and_ratio_cojo_SNPs_Heart_Atrial_Appendage_eqtl_filtered.csv",
             "clsa_metabolite_and_ratio_cojo_SNPs_Artery_Aorta_eqtl_filtered.csv", "clsa_metabolite_and_ratio_cojo_SNPs_Heart_Left_Ventricle_eqtl_filtered.csv",
             "clsa_metabolite_and_ratio_cojo_SNPs_Artery_Coronary_eqtl_filtered.csv", "clsa_metabolite_and_ratio_cojo_SNPs_Kidney_Cortex_eqtl_filtered.csv",
             "clsa_metabolite_and_ratio_cojo_SNPs_Artery_Tibial_eqtl_filtered.csv", "clsa_metabolite_and_ratio_cojo_SNPs_Liver_eqtl_filtered.csv",
             "clsa_metabolite_and_ratio_cojo_SNPs_Brain_Amygdala_eqtl_filtered.csv", "clsa_metabolite_and_ratio_cojo_SNPs_Lung_eqtl_filtered.csv",
             "clsa_metabolite_and_ratio_cojo_SNPs_Brain_Anterior_cingulate_cortex_BA24_eqtl_filtered.csv",   "clsa_metabolite_and_ratio_cojo_SNPs_Minor_Salivary_Gland_eqtl_filtered.csv",
             "clsa_metabolite_and_ratio_cojo_SNPs_Brain_Caudate_basal_ganglia_eqtl_filtered.csv", "clsa_metabolite_and_ratio_cojo_SNPs_Muscle_Skeletal_eqtl_filtered.csv",
             "clsa_metabolite_and_ratio_cojo_SNPs_Brain_Cerebellar_Hemisphere_eqtl_filtered.csv", "clsa_metabolite_and_ratio_cojo_SNPs_Nerve_Tibial_eqtl_filtered.csv",
             "clsa_metabolite_and_ratio_cojo_SNPs_Brain_Cerebellum_eqtl_filtered.csv", "clsa_metabolite_and_ratio_cojo_SNPs_Ovary_eqtl_filtered.csv",
             "clsa_metabolite_and_ratio_cojo_SNPs_Brain_Cortex_eqtl_filtered.csv", "clsa_metabolite_and_ratio_cojo_SNPs_Pancreas_eqtl_filtered.csv",
             "clsa_metabolite_and_ratio_cojo_SNPs_Brain_Frontal_Cortex_BA9_eqtl_filtered.csv", "clsa_metabolite_and_ratio_cojo_SNPs_Pituitary_eqtl_filtered.csv",
             "clsa_metabolite_and_ratio_cojo_SNPs_Brain_Hippocampus_eqtl_filtered.csv", "clsa_metabolite_and_ratio_cojo_SNPs_Prostate_eqtl_filtered.csv",
             'clsa_metabolite_and_ratio_cojo_SNPs_Brain_Hypothalamus_eqtl_filtered.csv', "clsa_metabolite_and_ratio_cojo_SNPs_Skin_Not_Sun_Exposed_Suprapubic_eqtl_filtered.csv",
             "clsa_metabolite_and_ratio_cojo_SNPs_Brain_Nucleus_accumbens_basal_ganglia_eqtl_filtered.csv",  "clsa_metabolite_and_ratio_cojo_SNPs_Skin_Sun_Exposed_Lower_leg_eqtl_filtered.csv",
             "clsa_metabolite_and_ratio_cojo_SNPs_Brain_Putamen_basal_ganglia_eqtl_filtered.csv", "clsa_metabolite_and_ratio_cojo_SNPs_Small_Intestine_Terminal_Ileum_eqtl_filtered.csv",
             "clsa_metabolite_and_ratio_cojo_SNPs_Brain_Spinal_cord_cervical_c-1_eqtl_filtered.csv", "clsa_metabolite_and_ratio_cojo_SNPs_Spleen_eqtl_filtered.csv",
             "clsa_metabolite_and_ratio_cojo_SNPs_Brain_Substantia_nigra_eqtl_filtered.csv", "clsa_metabolite_and_ratio_cojo_SNPs_Stomach_eqtl_filtered.csv",
             "clsa_metabolite_and_ratio_cojo_SNPs_Breast_Mammary_Tissue_eqtl_filtered.csv", "clsa_metabolite_and_ratio_cojo_SNPs_Testis_eqtl_filtered.csv",
             "clsa_metabolite_and_ratio_cojo_SNPs_Cells_Cultured_fibroblasts_eqtl_filtered.csv", "clsa_metabolite_and_ratio_cojo_SNPs_Thyroid_eqtl_filtered.csv",
             "clsa_metabolite_and_ratio_cojo_SNPs_Cells_EBV-transformed_lymphocytes_eqtl_filtered.csv", "clsa_metabolite_and_ratio_cojo_SNPs_Uterus_eqtl_filtered.csv",
             "clsa_metabolite_and_ratio_cojo_SNPs_Colon_Sigmoid_eqtl_filtered.csv", "clsa_metabolite_and_ratio_cojo_SNPs_Vagina_eqtl_filtered.csv",
             "clsa_metabolite_and_ratio_cojo_SNPs_Colon_Transverse_eqtl_filtered.csv", "clsa_metabolite_and_ratio_cojo_SNPs_whole_blood_eqtl_filtered.csv",
             "clsa_metabolite_and_ratio_cojo_SNPs_Esophagus_Gastroesophageal_Junction_eqtl_filtered.csv")
all_eqtl_dat<-data.frame(
  SNP=character(),
  SNP_ID=character(),
  gene_id=character(),
  tss_distance=integer(),
  ma_samples=integer(),
  ma_count=integer(), 
  maf=numeric(),
  pval_nominal=numeric(),
  slope=numeric(),
  slope_se=numeric(),
  pval_nominal_threshold=numeric(),
  min_pval_nominal=numeric(),
  pval_beta=numeric(),
  gene_id_2=character(),
  gene_symbol=character())

for (file in file_list){
  eqtl_dat<-read.csv(paste0("./analysis_dir/expression_relevance_check/eqtl/",file), row.names = 1)
  tissue<-str_split(str_split(file,"SNPs_")[[1]][2], "_eqtl")[[1]][1]
  eqtl_dat_sub <-eqtl_dat %>% select(SNP,SNP_ID,phenotype_id:gene_symbol)
  eqtl_dat_sub$tissue<-rep(tissue, nrow(eqtl_dat_sub))
  all_eqtl_dat<-rbind(all_eqtl_dat,eqtl_dat_sub)
}
all_eqtl_dat_distinct<-all_eqtl_dat %>% distinct()

write.csv(all_eqtl_dat_distinct,"./analysis_dir/expression_relevance_check/eqtl/all_eqtl_dat_all_49_tissues.csv")
