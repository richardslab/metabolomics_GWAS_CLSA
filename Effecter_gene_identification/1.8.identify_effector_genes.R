library(data.table)
library(dplyr)
library(stringr)

clsa_metabolite_cojo_SNPs_formatted_rsid_filtered<-read.csv("./results/clsa_metabolites_cojo_formatted_rsid_formatted_filtered.csv",row.name=1)
clsa_metabolite_cojo_SNPs_formatted_rsid_filtered1<-clsa_metabolite_cojo_SNPs_formatted_rsid_filtered %>% rowwise() %>% mutate(metabo_snpid_pair=paste0(metabolite,SNP_ID))

clsa_metabolite_ratios_cojo_SNPs_formatted_rsid_filtered<-read.csv("./results/clsa_metabolite_ratios_cojo_SNPs_formatted_rsid_bonf_filtered_enzyme_transporter.csv",row.name=1)
clsa_metabolite_ratios_cojo_SNPs_formatted_rsid_filtered1<-clsa_metabolite_ratios_cojo_SNPs_formatted_rsid_filtered %>% rowwise() %>% mutate(metabo_snpid_pair=paste0(metabolite,SNP_ID))

##metabo id full name match
metabo_id_full_name_match<-rbind(clsa_metabolite_cojo_SNPs_formatted_rsid_filtered1[,c(11,21)], clsa_metabolite_ratios_cojo_SNPs_formatted_rsid_filtered1[,c(11,25)])
metabo_id_full_name_match_unique<-metabo_id_full_name_match %>% distinct()

### load the file containing all gene within 1Mb --> a list of snp gene pairs with metabolism-related genes
all_IVs_nearby_genes<-read.csv("./results/all_IVs_nearby_genes.csv",row.names = 1)
clsa_metabolite_and_ratio_cojo_SNPs_filtered_SNPID<-rbind(clsa_metabolite_cojo_SNPs_formatted_rsid_filtered1[,c(1:11,18)], clsa_metabolite_ratios_cojo_SNPs_formatted_rsid_filtered1[,c(1:11,18)])
rsid_snpid_match<-clsa_metabolite_and_ratio_cojo_SNPs_filtered_SNPID %>% select(SNP,SNP_ID) %>% distinct()

############# all the snps with biological relevance ############# 
metabo_metabolism_related_genes<-read.csv("./analysis_dir/biological_relevance_check/metabolite_cojo_snp_metabolism_related_gene_checked.csv", row.names = 1)
metabo_ratio_metabolism_related_genes<-read.csv("./analysis_dir/biological_relevance_check/metabolite_ratio_cojo_snp_metabolism_related_gene_checked.csv", row.names = 1)
## add rsID back
metabo_metabolism_related_genes<-left_join(metabo_metabolism_related_genes,rsid_snpid_match, by=c("SNP_ID"="SNP_ID"))
metabo_ratio_metabolism_related_genes<-left_join(metabo_ratio_metabolism_related_genes,rsid_snpid_match, by=c("SNP_ID"="SNP_ID"))

metabo_metabolism_related_genes<-metabo_metabolism_related_genes %>% filter(gene_symbol !="" & is.na(gene_symbol) == F) %>% rowwise() %>% mutate(test_metabo = paste0(metabolite, SNP_ID), test_gene = paste0(SNP_ID, gene_symbol), test_metabo_gene=paste0(metabolite,SNP_ID, gene_symbol))
metabo_ratio_metabolism_related_genes<-metabo_ratio_metabolism_related_genes %>% filter(gene_symbol !="" & is.na(gene_symbol) == F) %>% rowwise() %>% mutate(test_metabo = paste0(metabolite, SNP_ID), test_gene = paste0(SNP_ID, gene_symbol), test_metabo_gene=paste0(metabolite,SNP_ID, gene_symbol))

biological_rela_metabo_list<-metabo_metabolism_related_genes %>% filter(HMDB_gene_pass == "Y" | KEGG_gene_pass == "Y" | PubChem_gene_pass == "Y") %>% select(test_metabo, test_gene, test_metabo_gene, SNP, gene_symbol, SNP_ID, metabolite)
biological_rela_metabo_ratio_list<-metabo_ratio_metabolism_related_genes %>% filter(HMDB_gene_pass == "Y" | KEGG_gene_pass == "Y" | PubChem_gene_pass == "Y") %>% select(test_metabo, test_gene, test_metabo_gene, SNP, gene_symbol, SNP_ID, metabolite)
tier1_list<-rbind(biological_rela_metabo_list, biological_rela_metabo_ratio_list)

clsa_metabolite_and_ratio_cojo_SNPs_filtered_SNPID1<-clsa_metabolite_and_ratio_cojo_SNPs_filtered_SNPID%>% rowwise() %>% 
  mutate(test_metabo= paste0(metabolite,SNP_ID))

tier1_list_with_sig_associations<-tier1_list %>% rowwise() %>% mutate(sig_metabo_snp_asso=ifelse(test_metabo %in% clsa_metabolite_and_ratio_cojo_SNPs_filtered_SNPID1$test_metabo, "Y","N"))%>%filter(sig_metabo_snp_asso == "Y")

tier1_list_with_sig_associations1<-tier1_list_with_sig_associations%>%rowwise()%>%mutate(SNP_fix=ifelse(is.na(str_extract(SNP, "rs"))==F, SNP, SNP_ID))
tier1_list_with_sig_associations1$full_name<-apply(tier1_list_with_sig_associations1,1,function(x) metabo_id_full_name_match_unique[metabo_id_full_name_match_unique$metabolite == x[7],]$full_name[1])
biological_relevant_genes<-tier1_list_with_sig_associations1 %>%rowwise() %>% mutate(snp_metabo=paste0(SNP_fix,full_name))%>% select(snp_metabo, gene_symbol) 
biological_relevant_genes1<-biological_relevant_genes%>% group_by(snp_metabo) %>% summarise(biological_genes=paste(unique(gene_symbol),collapse="; "))

write.csv(biological_relevant_genes1, "./analysis_dir/biological_relevance_check/biological_relevance_genes_all.csv")

############## eqtl / sqtl genes #############
### load the file containing the eqtl-related gene
all_eqtl_dat_all_tissues<-read.csv("./analysis_dir/expression_relevance_check/eqtl/all_eqtl_dat_all_49_tissues.csv",row.names = 1)
all_eqtl_dat_all_tissues<-all_eqtl_dat_all_tissues %>% filter(gene_symbol !="" & is.na(gene_symbol) == F) %>% rowwise() %>% mutate(test_snp_gene = paste0(SNP_ID, gene_symbol))
### load the file containing the sqtl-related gene
all_sqtl_dat_all_tissues<-read.csv("./analysis_dir/expression_relevance_check/sqtl/all_sqtl_dat_EUR_23_tissues.csv",row.names = 1)
all_sqtl_dat_all_tissues<-all_sqtl_dat_all_tissues %>% filter(gene_symbol !="" & is.na(gene_symbol) == F) %>% rowwise() %>% mutate(test_snp_gene = paste0(SNP_ID, gene_symbol))
eqtl_sqtl_list<-c(all_eqtl_dat_all_tissues$test_snp_gene, all_sqtl_dat_all_tissues$test_snp_gene)

## mark the nearby genes (of the SNPs that are eqtl or sqtl)
clsa_metabolite_and_ratio_cojo_SNPs_filtered_SNPID_gene<-left_join(clsa_metabolite_and_ratio_cojo_SNPs_filtered_SNPID, all_IVs_nearby_genes[,c(4,10)], by=c("SNP_ID"="SNP"))
clsa_metabolite_and_ratio_cojo_SNPs_filtered_SNPID_gene1<-clsa_metabolite_and_ratio_cojo_SNPs_filtered_SNPID_gene %>% rowwise() %>% mutate(test_metabo = paste0(metabolite, SNP_ID), test_gene = paste0(SNP_ID, gene_symbol), test_metabo_gene=paste0(metabolite,SNP_ID, gene_symbol))

clsa_metabolite_and_ratio_cojo_SNPs_filtered_SNPID_gene2<-clsa_metabolite_and_ratio_cojo_SNPs_filtered_SNPID_gene1 %>% rowwise() %>% mutate(eqtl_check = ifelse(gene_symbol %in% all_eqtl_dat_all_tissues[all_eqtl_dat_all_tissues$SNP_ID == SNP_ID,]$gene_symbol, "Y", "N"))
clsa_metabolite_and_ratio_cojo_SNPs_filtered_SNPID_gene3<-clsa_metabolite_and_ratio_cojo_SNPs_filtered_SNPID_gene2 %>% rowwise() %>% mutate(sqtl_check = ifelse(gene_symbol %in% all_sqtl_dat_all_tissues[all_sqtl_dat_all_tissues$SNP_ID == SNP_ID,]$gene_symbol, "Y", "N"))
clsa_metabolite_and_ratio_cojo_SNPs_filtered_SNPID_gene4<- clsa_metabolite_and_ratio_cojo_SNPs_filtered_SNPID_gene3 %>%rowwise()%>%mutate(SNP_fix=ifelse(is.na(str_extract(SNP, "rs"))==F, SNP, SNP_ID))
clsa_metabolite_and_ratio_cojo_SNPs_filtered_SNPID_gene5<-clsa_metabolite_and_ratio_cojo_SNPs_filtered_SNPID_gene4 %>%filter(eqtl_check == "Y" | sqtl_check == "Y")

expression_relavant_genes_metabo<-clsa_metabolite_and_ratio_cojo_SNPs_filtered_SNPID_gene5 %>%  filter(!(str_detect(metabolite, "_")))%>% select(SNP_fix,gene_symbol) %>% group_by(SNP_fix) %>% summarise(expression_genes=paste(unique(gene_symbol),collapse="; "))
expression_relavant_genes_metabo_ratios<-clsa_metabolite_and_ratio_cojo_SNPs_filtered_SNPID_gene5 %>%  filter(str_detect(metabolite, "_"))%>%select(SNP_fix,gene_symbol) %>% group_by(SNP_fix) %>% summarise(expression_genes=paste(unique(gene_symbol),collapse="; "))

write.csv(expression_relavant_genes_metabo, "./analysis_dir/expression_relevance_check/expression_relavant_genes_metabo.csv")
write.csv(expression_relavant_genes_metabo_ratios, "./analysis_dir/expression_relevance_check/expression_relavant_genes_metabo_ratios.csv")

###### include eqtl and sqtl coloc results
eqtl_coloc_results_metabo_names_pph4_sig<-read.csv("./restuls/coloc_results/Summary_eqtl_coloc_results.csv", row.names = 1)
sqtl_coloc_results_metabo_names_pph4_sig<-read.csv("./restuls/coloc_results/Summary_sqtl_coloc_results.csv", row.names = 1)
clsa_metabolite_and_ratio_cojo_SNPs_filtered_SNPID_gene6<-clsa_metabolite_and_ratio_cojo_SNPs_filtered_SNPID_gene5%>%rowwise()%>%mutate(test_ID = paste0(metabolite,SNP),
                                                                                                                                        test_ID_fixed = paste0(metabolite,str_replace_all(SNP_ID,":","_")))
##check nearby gene eqtl coloc
clsa_metabolite_and_ratio_cojo_SNPs_filtered_SNPID_gene7<-clsa_metabolite_and_ratio_cojo_SNPs_filtered_SNPID_gene6 %>% rowwise() %>% 
  mutate(eqtl_coloc_check = ifelse((eqtl_check =="Y" & (gene_symbol %in% eqtl_coloc_results_metabo_names_pph4_sig_previous_run[eqtl_coloc_results_metabo_names_pph4_sig_previous_run$test_ID == test_ID,]$gene_symbol) | (gene_symbol %in% eqtl_coloc_results_metabo_names_pph4_sig_additional_run[eqtl_coloc_results_metabo_names_pph4_sig_additional_run$testID_fixed_test == test_ID_fixed,]$gene_symbol)), "Y", "N"))

##check nearby gene sqtl coloc
clsa_metabolite_and_ratio_cojo_SNPs_filtered_SNPID_gene8<-clsa_metabolite_and_ratio_cojo_SNPs_filtered_SNPID_gene7 %>% rowwise() %>% 
  mutate(sqtl_coloc_check = ifelse((sqtl_check =="Y" & (gene_symbol %in% sqtl_coloc_results_metabo_names_pph4_sig_previous_run[sqtl_coloc_results_metabo_names_pph4_sig_previous_run$test_ID == test_ID,]$gene_symbol) | (gene_symbol %in% sqtl_coloc_results_metabo_names_pph4_sig_additional_run[sqtl_coloc_results_metabo_names_pph4_sig_additional_run$testID_fixed_test == test_ID_fixed,]$gene_symbol)), "Y", "N"))

clsa_metabolite_and_ratio_cojo_SNPs_filtered_SNPID_gene9<-clsa_metabolite_and_ratio_cojo_SNPs_filtered_SNPID_gene8 %>%filter(eqtl_coloc_check == "Y" | sqtl_coloc_check == "Y")
expression_relavant_colocalized_genes_metabo<-clsa_metabolite_and_ratio_cojo_SNPs_filtered_SNPID_gene9 %>%  filter(!(str_detect(metabolite, "_")))%>% select(SNP_fix,gene_symbol) %>% group_by(SNP_fix) %>% summarise(expression_genes=paste(unique(gene_symbol),collapse="; "))
expression_relavant_colocalized_genes_metabo_ratios<-clsa_metabolite_and_ratio_cojo_SNPs_filtered_SNPID_gene9 %>%  filter(str_detect(metabolite, "_"))%>%select(SNP_fix,gene_symbol) %>% group_by(SNP_fix) %>% summarise(expression_genes=paste(unique(gene_symbol),collapse="; "))

write.csv(expression_relavant_colocalized_genes_metabo, "./analysis_dir/expression_relevance_check/expression_relavant_colocalized_genes_metabo.csv")
write.csv(expression_relavant_colocalized_genes_metabo_ratios, "./analysis_dir/expression_relevance_check/expression_relavant_colocalized_genes_metabo_ratios.csv")

############# filter for tier 1 snp-metabolites associated and the SNPs are associated with eqtl or sqtl of the gene that are metabolize that metabolites
# Potential effector genes are corresponding to each SNP-metabolite pairs --> the gene are metabolism-related and nearby snp and associated (coloc) with the eqtl / sqtl of the snp
#clsa_metabolite_and_ratio_cojo_SNPs_filtered_SNPID_gene9 include snp with colocalized cis-eqtl or cis-sqtl
clsa_metabolite_and_ratio_cojo_SNPs_filtered_SNPID_gene1_credible<-clsa_metabolite_and_ratio_cojo_SNPs_filtered_SNPID_gene9 %>% filter(test_metabo_gene %in% tier1_list_with_sig_associations1$test_metabo_gene)
clsa_metabolite_and_ratio_cojo_SNPs_filtered_SNPID_gene1_credible$full_name<-apply(clsa_metabolite_and_ratio_cojo_SNPs_filtered_SNPID_gene1_credible,1,function(x) metabo_id_full_name_match_unique[metabo_id_full_name_match_unique$metabolite == x[11],]$full_name[1])
clsa_metabolite_and_ratio_cojo_SNPs_filtered_SNPID_gene1_credible<-clsa_metabolite_and_ratio_cojo_SNPs_filtered_SNPID_gene1_credible%>%rowwise()%>%mutate(SNP_fix=ifelse(is.na(str_extract(SNP, "rs"))==F, SNP, SNP_ID))
potential_effector_genes<-clsa_metabolite_and_ratio_cojo_SNPs_filtered_SNPID_gene1_credible %>%rowwise() %>% mutate(snp_metabo=paste0(SNP_fix,full_name))%>% select(SNP_fix, full_name,snp_metabo, gene_symbol) 
potential_effector_genes1<-potential_effector_genes%>% group_by(SNP_fix,full_name) %>% summarise(effector_genes=paste(unique(gene_symbol),collapse="; "))
write.csv(potential_effector_genes1, "./analysis_dir/potential_effector_genes_with_eqtl_sqtl_coloc.csv")
