library(data.table)
library(dplyr)
library(stringr)
library(tidyr)

clsa_metabolite_cojo_SNPs_filtered<-read.csv("./results/clsa_metabolites_cojo_formatted_rsid_formatted_filtered.csv",row.names = 1)
clsa_metabolite_ratios_cojo_SNPs_formatted_rsid_filtered<-read.csv("./results/clsa_metabolite_ratios_cojo_SNPs_formatted_rsid_bonf_filtered_enzyme_transporter.csv",row.names = 1)
clsa_metabolite_and_ratio_cojo_SNPs_filtered<-rbind(clsa_metabolite_cojo_SNPs_filtered[,c(1:11,18)], clsa_metabolite_ratios_cojo_SNPs_formatted_rsid_filtered[,c(1:11,18)])

#### get gene annotations
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
ensembl = useMart(biomart="ENSEMBL_MART_ENSEMBL",
                  dataset="hsapiens_gene_ensembl", 
                  host="uswest.ensembl.org",
                  ensemblRedirect = FALSE)

############################## sqtl  ############################## 
################ whole blood sqtl ################ 
#combine the SNPs from metabolites and metabolite ratios

Whole_Blood.v8.EUR.signif_variant_gene_pairs<-fread("./data/GTEX_eqtl_database/sqtls_EUR/Whole_Blood.v8.EUR.signif_pairs.txt.gz")

clsa_metabolite_and_ratio_cojo_SNPs_filtered<-clsa_metabolite_and_ratio_cojo_SNPs_filtered %>% rowwise() %>%
  mutate(variant_id_1 = paste0(str_split(SNP_ID,":",4)[[1]][1],"_",str_split(SNP_ID,":",4)[[1]][2],"_",str_split(SNP_ID,":",4)[[1]][3],"_",str_split(SNP_ID,":",4)[[1]][4],"_b38"),
         variant_id_2 = paste0(str_split(SNP_ID,":",4)[[1]][1],"_",str_split(SNP_ID,":",4)[[1]][2],"_",str_split(SNP_ID,":",4)[[1]][4],"_",str_split(SNP_ID,":",4)[[1]][3],"_b38"))

clsa_metabolite_and_ratio_cojo_SNPs_filtered_whole_blood_sqtl<-left_join(clsa_metabolite_and_ratio_cojo_SNPs_filtered,Whole_Blood.v8.EUR.signif_variant_gene_pairs, by=c("variant_id_1"="variant_id"))
clsa_metabolite_and_ratio_cojo_SNPs_filtered_whole_blood_sqtl_no_match<-clsa_metabolite_and_ratio_cojo_SNPs_filtered_whole_blood_sqtl %>% filter(is.na(phenotype_id)==T)
clsa_metabolite_and_ratio_cojo_SNPs_filtered_whole_blood_sqtl_no_match<-left_join(clsa_metabolite_and_ratio_cojo_SNPs_filtered_whole_blood_sqtl_no_match,Whole_Blood.v8.EUR.signif_variant_gene_pairs, by=c("variant_id_2"="variant_id"))## no matching found with variant_id_2
clsa_metabolite_and_ratio_cojo_SNPs_filtered_whole_blood_sqtl_filtered<-clsa_metabolite_and_ratio_cojo_SNPs_filtered_whole_blood_sqtl %>% filter(is.na(phenotype_id)!=T)

### get the gene symbol for the sqtl
clsa_metabolite_and_ratio_cojo_SNPs_filtered_whole_blood_sqtl_filtered<-clsa_metabolite_and_ratio_cojo_SNPs_filtered_whole_blood_sqtl_filtered%>%rowwise()%>%mutate(gene_id=str_split(phenotype_id,":")[[1]][5])
clsa_metabolite_and_ratio_cojo_SNPs_filtered_whole_blood_sqtl_filtered$gene_id_2 <- sub("[.][0-9]*","",clsa_metabolite_and_ratio_cojo_SNPs_filtered_whole_blood_sqtl_filtered$gene_id)
genes <- clsa_metabolite_and_ratio_cojo_SNPs_filtered_whole_blood_sqtl_filtered$gene_id_2
gene_IDs <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),
                  values = genes, mart= mart)
clsa_metabolite_and_ratio_cojo_SNPs_filtered_whole_blood_sqtl_filtered$gene_symbol<-apply(clsa_metabolite_and_ratio_cojo_SNPs_filtered_whole_blood_sqtl_filtered,1,function(x) gene_IDs$hgnc_symbol[gene_IDs$ensembl_gene_id == x[27]][1])
write.csv(clsa_metabolite_and_ratio_cojo_SNPs_filtered_whole_blood_sqtl_filtered,"./analysis_dir/expression_relevance_check/sqtl/clsa_metabolite_and_ratio_cojo_SNPs_whole_blood_sqtl_filtered.csv")

################ tissue-specific sqtl ################
sqtl_list<-c("Adipose_Subcutaneous.v8.EUR.signif_pairs.txt.gz", "Colon_Transverse.v8.EUR.signif_pairs.txt.gz", "Lung.v8.EUR.signif_pairs.txt.gz",
             "Adipose_Visceral_Omentum.v8.EUR.signif_pairs.txt.gz", "Esophagus_Gastroesophageal_Junction.v8.EUR.signif_pairs.txt.gz", "Minor_Salivary_Gland.v8.EUR.signif_pairs.txt.gz",
             "Adrenal_Gland.v8.EUR.signif_pairs.txt.gz", "Esophagus_Mucosa.v8.EUR.signif_pairs.txt.gz", "Muscle_Skeletal.v8.EUR.signif_pairs.txt.gz",
             "Artery_Coronary.v8.EUR.signif_pairs.txt.gz", "Esophagus_Muscularis.v8.EUR.signif_pairs.txt.gz","Nerve_Tibial.v8.EUR.signif_pairs.txt.gz",
             "Artery_Tibial.v8.EUR.signif_pairs.txt.gz", "Heart_Atrial_Appendage.v8.EUR.signif_pairs.txt.gz","Ovary.v8.EUR.signif_pairs.txt.gz",
             "Cells_Cultured_fibroblasts.v8.EUR.signif_pairs.txt.gz", "Heart_Left_Ventricle.v8.EUR.signif_pairs.txt.gz", "Pancreas.v8.EUR.signif_pairs.txt.gz",
             "Cells_EBV-transformed_lymphocytes.v8.EUR.signif_pairs.txt.gz", "Kidney_Cortex.v8.EUR.signif_pairs.txt.gz",
             "Colon_Sigmoid.v8.EUR.signif_pairs.txt.gz","Liver.v8.EUR.signif_pairs.txt.gz")

for (i in 1:length(sqtl_list)) {
  tissue<-as.character(str_split(sqtl_list[i],".v8")[[1]][1])
  clsa_metabolite_and_ratio_cojo_SNPs_filtered_1<-clsa_metabolite_and_ratio_cojo_SNPs_filtered %>% rowwise() %>%mutate(variant_id_1 = paste0(str_split(SNP_ID,":",4)[[1]][1],"_",str_split(SNP_ID,":",4)[[1]][2],"_",str_split(SNP_ID,":",4)[[1]][3],"_",str_split(SNP_ID,":",4)[[1]][4],"_b38"),variant_id_2 = paste0(str_split(SNP_ID,":",4)[[1]][1],"_",str_split(SNP_ID,":",4)[[1]][2],"_",str_split(SNP_ID,":",4)[[1]][4],"_",str_split(SNP_ID,":",4)[[1]][3],"_b38"))
  tissue_specific_signif_variant_gene_pairs<-fread(paste0("./data/GTEX_eqtl_database/sqtls_EUR/",sqtl_list[i]))
  clsa_metabolite_and_ratio_cojo_SNPs_filtered_sqtl<-left_join(clsa_metabolite_and_ratio_cojo_SNPs_filtered_1,tissue_specific_signif_variant_gene_pairs, by=c("variant_id_1"="variant_id"))
  clsa_metabolite_and_ratio_cojo_SNPs_filtered_sqtl_no_match<-clsa_metabolite_and_ratio_cojo_SNPs_filtered_sqtl %>% filter(is.na(phenotype_id)==T)
  clsa_metabolite_and_ratio_cojo_SNPs_filtered_sqtl_no_match_checked<-left_join(clsa_metabolite_and_ratio_cojo_SNPs_filtered_1,tissue_specific_signif_variant_gene_pairs, by=c("variant_id_2"="variant_id"))## no matching found with variant_id_2
  ## checked, variant_id_1 is enough for matching to variant_id in sqtl database
  clsa_metabolite_and_ratio_cojo_SNPs_filtered_sqtl_no_match_check<-clsa_metabolite_and_ratio_cojo_SNPs_filtered_sqtl_no_match_checked%>%filter(is.na(phenotype_id)!=T)
  if (nrow(clsa_metabolite_and_ratio_cojo_SNPs_filtered_sqtl_no_match_check)==0) {
    clsa_metabolite_and_ratio_cojo_SNPs_filtered_sqtl_filtered<-clsa_metabolite_and_ratio_cojo_SNPs_filtered_sqtl %>% filter(is.na(phenotype_id)!=T)
    clsa_metabolite_and_ratio_cojo_SNPs_filtered_sqtl_filtered<-clsa_metabolite_and_ratio_cojo_SNPs_filtered_sqtl_filtered%>%rowwise()%>%mutate(gene_id=str_split(phenotype_id,":")[[1]][5])
    clsa_metabolite_and_ratio_cojo_SNPs_filtered_sqtl_filtered$gene_id_2 <-sub("[.][0-9]*","",clsa_metabolite_and_ratio_cojo_SNPs_filtered_sqtl_filtered$gene_id)
    genes <- clsa_metabolite_and_ratio_cojo_SNPs_filtered_sqtl_filtered$gene_id_2
    gene_IDs <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),values = genes, mart= mart)
    clsa_metabolite_and_ratio_cojo_SNPs_filtered_sqtl_filtered$gene_symbol<-apply(clsa_metabolite_and_ratio_cojo_SNPs_filtered_sqtl_filtered,1,function(x) gene_IDs$hgnc_symbol[gene_IDs$ensembl_gene_id == x[27]][1])
    write.csv(clsa_metabolite_and_ratio_cojo_SNPs_filtered_sqtl_filtered,paste0("./analysis_dir/expression_relevance_check/sqtl/clsa_metabolite_and_ratio_cojo_SNPs_",tissue,"_sqtl_filtered.csv"))
  }
  else {clsa_metabolite_and_ratio_cojo_SNPs_filtered_sqtl<-rbind(clsa_metabolite_and_ratio_cojo_SNPs_filtered_sqtl,clsa_metabolite_and_ratio_cojo_SNPs_filtered_sqtl_no_match_check)
  clsa_metabolite_and_ratio_cojo_SNPs_filtered_sqtl_filtered<-clsa_metabolite_and_ratio_cojo_SNPs_filtered_sqtl %>% filter(is.na(phenotype_id)!=T)
  clsa_metabolite_and_ratio_cojo_SNPs_filtered_sqtl_filtered<-clsa_metabolite_and_ratio_cojo_SNPs_filtered_sqtl_filtered%>%rowwise()%>%mutate(gene_id=str_split(phenotype_id,":")[[1]][5])
  clsa_metabolite_and_ratio_cojo_SNPs_filtered_sqtl_filtered$gene_id_2 <- sub("[.][0-9]*","",clsa_metabolite_and_ratio_cojo_SNPs_filtered_sqtl_filtered$gene_id)
  genes <- clsa_metabolite_and_ratio_cojo_SNPs_filtered_sqtl_filtered$gene_id_2
  gene_IDs <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),values = genes, mart= mart)
  clsa_metabolite_and_ratio_cojo_SNPs_filtered_sqtl_filtered$gene_symbol<-apply(clsa_metabolite_and_ratio_cojo_SNPs_filtered_sqtl_filtered,1,function(x) gene_IDs$hgnc_symbol[gene_IDs$ensembl_gene_id == x[27]][1])
  write.csv(clsa_metabolite_and_ratio_cojo_SNPs_filtered_sqtl_filtered,paste0("./analysis_dir/expression_relevance_check/sqtl/clsa_metabolite_and_ratio_cojo_SNPs_",tissue,"_sqtl_filtered.csv"))
  }
}

###### combine all tissue-specific sqtl to one file
file_list<-c("clsa_metabolite_and_ratio_cojo_SNPs_Adipose_Subcutaneous_sqtl_filtered.csv", "clsa_metabolite_and_ratio_cojo_SNPs_Heart_Atrial_Appendage_sqtl_filtered.csv",
             "clsa_metabolite_and_ratio_cojo_SNPs_Adipose_Visceral_Omentum_sqtl_filtered.csv", "clsa_metabolite_and_ratio_cojo_SNPs_Heart_Left_Ventricle_sqtl_filtered.csv",
             "clsa_metabolite_and_ratio_cojo_SNPs_Adrenal_Gland_sqtl_filtered.csv","clsa_metabolite_and_ratio_cojo_SNPs_Kidney_Cortex_sqtl_filtered.csv",
             "clsa_metabolite_and_ratio_cojo_SNPs_Artery_Coronary_sqtl_filtered.csv", "clsa_metabolite_and_ratio_cojo_SNPs_Liver_sqtl_filtered.csv",
             "clsa_metabolite_and_ratio_cojo_SNPs_Artery_Tibial_sqtl_filtered.csv", "clsa_metabolite_and_ratio_cojo_SNPs_Lung_sqtl_filtered.csv",
             "clsa_metabolite_and_ratio_cojo_SNPs_Cells_Cultured_fibroblasts_sqtl_filtered.csv", "clsa_metabolite_and_ratio_cojo_SNPs_Minor_Salivary_Gland_sqtl_filtered.csv",
             "clsa_metabolite_and_ratio_cojo_SNPs_Cells_EBV-transformed_lymphocytes_sqtl_filtered.csv", "clsa_metabolite_and_ratio_cojo_SNPs_Muscle_Skeletal_sqtl_filtered.csv",
             "clsa_metabolite_and_ratio_cojo_SNPs_Colon_Sigmoid_sqtl_filtered.csv", "clsa_metabolite_and_ratio_cojo_SNPs_Nerve_Tibial_sqtl_filtered.csv",
             "clsa_metabolite_and_ratio_cojo_SNPs_Colon_Transverse_sqtl_filtered.csv", "clsa_metabolite_and_ratio_cojo_SNPs_Ovary_sqtl_filtered.csv",
             "clsa_metabolite_and_ratio_cojo_SNPs_Esophagus_Gastroesophageal_Junction_sqtl_filtered.csv",  "clsa_metabolite_and_ratio_cojo_SNPs_Pancreas_sqtl_filtered.csv",
             "clsa_metabolite_and_ratio_cojo_SNPs_Esophagus_Mucosa_sqtl_filtered.csv","clsa_metabolite_and_ratio_cojo_SNPs_whole_blood_sqtl_filtered.csv",
             "clsa_metabolite_and_ratio_cojo_SNPs_Esophagus_Muscularis_sqtl_filtered.csv")
all_sqtl_dat<-data.frame(
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
  sqtl_dat<-read.csv(paste0("./analysis_dir/expression_relevance_check/sqtl/",file), row.names = 1)
  tissue<-str_split(str_split(file,"SNPs_")[[1]][2], "_sqtl")[[1]][1]
  sqtl_dat_sub <-sqtl_dat%>% select(SNP,SNP_ID,phenotype_id:gene_symbol)
  sqtl_dat_sub$tissue<-rep(tissue, nrow(sqtl_dat_sub))
  all_sqtl_dat<-rbind(all_sqtl_dat,sqtl_dat_sub)
}
all_sqtl_dat_distinct<-all_sqtl_dat %>% distinct()

write.csv(all_sqtl_dat_distinct,"./analysis_dir/expression_relevance_check/sqtl/all_sqtl_dat_EUR_23_tissues.csv")
