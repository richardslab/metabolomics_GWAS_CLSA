library(data.table)
library(dplyr)
library(stringr)
library(tidyr)

clsa_metabolite_cojo_SNPs_filtered<-read.csv("./results/clsa_metabolites_cojo_formatted_rsid_formatted_filtered.csv",row.names = 1)
clsa_metabolite_ratio_cojo_SNPs_filtered<-read.csv("./results/clsa_metabolite_ratios_cojo_SNPs_formatted_rsid_bonf_filtered_enzyme_transporter.csv",row.names = 1)

all_IVs_nearby_genes<-read.csv("./results/all_IVs_nearby_genes.csv",row.names = 1)

matabolites_annotation<-read.csv("./data/clsa_phenotypes/clsa_metabo_annotation_dat")

##################### the gene involved in the metabolism of the metabolites #######################
###Step 1. add nearby genes to each SNPs
all_IVs_nearby_genes_formated_with_metabolite<-left_join(clsa_metabolite_cojo_SNPs_filtered,all_IVs_nearby_genes, by=c("SNP_ID"="SNP"))
###Step 2. get the HMDB and PubChem ID for all metabolites
all_IVs_nearby_genes_formated_with_metabolite_2<-left_join(all_IVs_nearby_genes_formated_with_metabolite,matabolites_annotation_July2[,c(2,17,20)], by=c("metabolite"="metabo_ID"))
all_IVs_nearby_genes_formated_with_metabolite_2$PUBCHEM_curated<-as.character(all_IVs_nearby_genes_formated_with_metabolite_2$PUBCHEM_curated)
all_IVs_nearby_genes_formated_with_metabolite_2_filtered<-all_IVs_nearby_genes_formated_with_metabolite_2 %>% filter((HMDB_curated !="" & HMDB_curated != 0)|(PUBCHEM_curated !="" & PUBCHEM_curated != 0))
###Step 3. prepare the data frame for bioloigcally relevance check
all_IVs_nearby_genes_formated_with_metabolite_3<-all_IVs_nearby_genes_formated_with_metabolite_2_filtered %>% select(SNP_ID, metabolite, CHR, POS, ncbi_gene_id, ensembl_gene_id, gene_symbol, HMDB_curated, PUBCHEM_curated)
all_IVs_nearby_genes_formated_with_metabolite_3$ncbi_gene_id<-as.character(all_IVs_nearby_genes_formated_with_metabolite_3$ncbi_gene_id)
## for the metabolite with more than 1 HMDB_ID --> split they into separate rows
all_IVs_nearby_genes_formated_with_metabolite_4<-all_IVs_nearby_genes_formated_with_metabolite_3 %>% mutate(HMDB_curated = strsplit(as.character(HMDB_curated), ",")) %>% unnest(HMDB_curated)
###Step 4. check if any of the nearby genes belong to metabolites' HMDB-related proteins/genes
HMDB_genes_clsa<-read.csv("./analysis_dir/biological_relevance_check/clsa_metabolites_HMDB_genes_with_metabo_id.csv", row.names = 1)
HMDB_genes_clsa$Entrez_ID<-as.character(HMDB_genes_clsa$Entrez_ID)
all_IVs_nearby_genes_formated_with_metabolite_4$HMDB_gene_pass<-apply(all_IVs_nearby_genes_formated_with_metabolite_4, 1, function(x) ifelse(x[5] %in% HMDB_genes_clsa[HMDB_genes_clsa$metabolite_ID == x[2],]$Entrez_ID, "Y", "N"))
all_IVs_nearby_genes_formated_with_metabolite_4$hmdb.protein.name<-apply(all_IVs_nearby_genes_formated_with_metabolite_4, 1, function(x)
  ifelse(x[10]== "Y", as.character(HMDB_genes_clsa[HMDB_genes_clsa$metabolite_ID == x[2] & HMDB_genes_clsa$Entrez_ID ==x[5],]$hmdb.protein.name), "NA"))
###Step 5. check if the nearby genes belong to metabolites' KEGG-related proteins/genes
KEGG_genes_clsa<-read.csv("./analysis_dir/biological_relevance_check/CLSA_metabolites_KEGG_pathway_genes_with_metabo_id.csv", row.names = 1)
KEGG_genes_clsa$ncbi.gene.id<-as.character(KEGG_genes_clsa$ncbi.gene.id)
all_IVs_nearby_genes_formated_with_metabolite_4$KEGG_gene_pass<-apply(all_IVs_nearby_genes_formated_with_metabolite_4, 1, function(x) ifelse(x[5] %in% KEGG_genes_clsa[KEGG_genes_clsa$metabolite_ID == x[2],]$ncbi.gene.id, "Y", "N"))
all_IVs_nearby_genes_formated_with_metabolite_4$kegg_ncbi_gene<-apply(all_IVs_nearby_genes_formated_with_metabolite_4, 1, function(x)
  ifelse(x[12]== "Y", as.character(KEGG_genes_clsa[KEGG_genes_clsa$metabolite_ID == x[2] & KEGG_genes_clsa$ncbi.gene.id ==x[5],]$ncbi.gene.desc), "NA"))
###Step 8. check if the nearby genes belong to metabolites' PubChem literature co-occurance proteins/genes
##### extract the Chemical-Gene Co-Occurrences in Literature from PubChem
metabolite_gene_co_occurancce_PubChem_clsa<-read.csv("./analysis_dir/biological_relevance_check/pubchem_metabo_gene_cooccurance_database_clsa.csv",row.names = 1)
metabolite_gene_co_occurancce_PubChem_clsa$entrez_ID<-as.character(metabolite_gene_co_occurancce_PubChem_clsa$entrez_ID)
class(metabolite_gene_co_occurancce_PubChem_clsa$PubChem_ID)<-as.character(metabolite_gene_co_occurancce_PubChem_clsa$PubChem_ID)
## for the metabolite with more than 1 Pubchem ID --> split they into separate rows
all_IVs_nearby_genes_formated_with_metabolite_5<-all_IVs_nearby_genes_formated_with_metabolite_4 %>% mutate(PUBCHEM_curated = strsplit(as.character(PUBCHEM_curated), ";")) %>% unnest(PUBCHEM_curated)
all_IVs_nearby_genes_formated_with_metabolite_5$PubChem_gene_pass<-apply(all_IVs_nearby_genes_formated_with_metabolite_5, 1, function(x) ifelse(x[5] %in% metabolite_gene_co_occurancce_PubChem_clsa[metabolite_gene_co_occurancce_PubChem_clsa$PubChem_ID == x[9],]$entrez_ID, "Y", "N"))
all_IVs_nearby_genes_formated_with_metabolite_5$PubChem_related_gene<-apply(all_IVs_nearby_genes_formated_with_metabolite_5, 1, function(x)
  ifelse(x[14]== "Y", as.character(metabolite_gene_co_occurancce_PubChem_clsa[metabolite_gene_co_occurancce_PubChem_clsa$PubChem_ID == x[9] & metabolite_gene_co_occurancce_PubChem_clsa$entrez_ID ==x[5],]$Gene_Symbol), "NA"))
all_IVs_nearby_genes_formated_with_metabolite_5$PubChem_lit_evidence<-apply(all_IVs_nearby_genes_formated_with_metabolite_5, 1, function(x)
  ifelse(x[14]== "Y", as.character(metabolite_gene_co_occurancce_PubChem_clsa[metabolite_gene_co_occurancce_PubChem_clsa$PubChem_ID == x[9] & metabolite_gene_co_occurancce_PubChem_clsa$entrez_ID ==x[5],]$evidence), "NA"))
all_IVs_nearby_genes_formated_with_metabolite_5$full_name<-apply(all_IVs_nearby_genes_formated_with_metabolite_5, 1, function(x) matabolites_annotation_July2[matabolites_annotation_July2$metabo_ID ==x[2],]$CHEMICAL_NAME[1])

write.csv(all_IVs_nearby_genes_formated_with_metabolite_5, "./analysis_dir/biological_relevance_check/metabolite_cojo_snp_metabolism_related_gene_checked.csv")

################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
##################### the gene involved in the metabolism of the metabolite pairs ####################### 
###Step 1. add nearby genes to each SNPs
all_IVs_nearby_genes_formated_with_metabolite_ratios<-left_join(clsa_metabolite_ratio_cojo_SNPs_filtered[,1:18],all_IVs_nearby_genes, by=c("SNP_ID"="SNP"))
###Step 2. split the metabolite pairs to two metabolite and then get the PubChem ID for all metabolites
all_IVs_nearby_genes_formated_with_metabolite_ratios_1<-all_IVs_nearby_genes_formated_with_metabolite_ratios %>% mutate(metabo= str_split(as.character(metabolite), "_")) %>% unnest(metabo)
all_IVs_nearby_genes_formated_with_metabolite_ratios_2<-left_join(all_IVs_nearby_genes_formated_with_metabolite_ratios_1,matabolites_annotation_July2[,c(2,17,20)], by=c("metabo"="metabo_ID"))
all_IVs_nearby_genes_formated_with_metabolite_ratios_2$PUBCHEM_curated<-as.character(all_IVs_nearby_genes_formated_with_metabolite_ratios_2$PUBCHEM_curated)
all_IVs_nearby_genes_formated_with_metabolite_ratios_2_filtered<-all_IVs_nearby_genes_formated_with_metabolite_ratios_2 %>% filter((HMDB_curated !="" & HMDB_curated != 0)|(PUBCHEM_curated !="" & PUBCHEM_curated != 0))
###Step 3. prepare effector genes
all_IVs_nearby_genes_formated_with_metabolite_ratios_3<-all_IVs_nearby_genes_formated_with_metabolite_ratios_2_filtered %>% select(SNP_ID, metabolite, CHR, POS, ncbi_gene_id, ensembl_gene_id, gene_symbol, HMDB_curated, PUBCHEM_curated, metabo)
all_IVs_nearby_genes_formated_with_metabolite_ratios_3$ncbi_gene_id<-as.character(all_IVs_nearby_genes_formated_with_metabolite_ratios_3$ncbi_gene_id)
## for the metabolite with more than 1 HMDB_ID --> split they into separate rows
all_IVs_nearby_genes_formated_with_metabolite_ratios_4<-all_IVs_nearby_genes_formated_with_metabolite_ratios_3 %>% mutate(HMDB_curated = strsplit(as.character(HMDB_curated), ",")) %>% unnest(HMDB_curated)
###Step 4. check if any of the nearby genes belong to metabolites' HMDB-related proteins/genes
HMDB_genes_clsa<-read.csv("./analysis_dir/biological_relevance_check/clsa_metabolites_HMDB_genes_with_metabo_id.csv", row.names = 1)
HMDB_genes_clsa$Entrez_ID<-as.character(HMDB_genes_clsa$Entrez_ID)
all_IVs_nearby_genes_formated_with_metabolite_ratios_4$HMDB_gene_pass<-apply(all_IVs_nearby_genes_formated_with_metabolite_ratios_4, 1, function(x) ifelse(x[5] %in% HMDB_genes_clsa[HMDB_genes_clsa$metabolite_ID == x[10],]$Entrez_ID, "Y", "N"))
all_IVs_nearby_genes_formated_with_metabolite_ratios_4$hmdb.protein.name<-apply(all_IVs_nearby_genes_formated_with_metabolite_ratios_4, 1, function(x) 
  ifelse(x[11]== "Y", as.character(HMDB_genes_clsa[HMDB_genes_clsa$metabolite_ID == x[10] & HMDB_genes_clsa$Entrez_ID ==x[5],]$hmdb.protein.name), "NA"))
###Step 5. check if the nearby genes belong to metabolites' KEGG-related proteins/genes
KEGG_genes_clsa<-read.csv("./analysis_dir/biological_relevance_check/CLSA_metabolites_KEGG_pathway_genes_with_metabo_id.csv", row.names = 1)
KEGG_genes_clsa$ncbi.gene.id<-as.character(KEGG_genes_clsa$ncbi.gene.id)
all_IVs_nearby_genes_formated_with_metabolite_ratios_4$KEGG_gene_pass<-apply(all_IVs_nearby_genes_formated_with_metabolite_ratios_4, 1, function(x) ifelse(x[5] %in% KEGG_genes_clsa[KEGG_genes_clsa$metabolite_ID == x[10],]$ncbi.gene.id, "Y", "N"))
all_IVs_nearby_genes_formated_with_metabolite_ratios_4$kegg_ncbi_gene<-apply(all_IVs_nearby_genes_formated_with_metabolite_ratios_4, 1, function(x) 
  ifelse(x[13]== "Y", as.character(KEGG_genes_clsa[KEGG_genes_clsa$metabolite_ID == x[10] & KEGG_genes_clsa$ncbi.gene.id ==x[5],]$ncbi.gene.desc), "NA"))
###Step 6. check if the nearby genes belong to metabolites' PubChem literature co-occurance proteins/genes
##### extract the Chemical-Gene Co-Occurrences in Literature from PubChem
metabolite_gene_co_occurancce_PubChem_clsa<-read.csv("./analysis_dir/biological_relevance_check/pubchem_metabo_gene_cooccurance_database_clsa.csv",row.names = 1)
metabolite_gene_co_occurancce_PubChem_clsa$entrez_ID<-as.character(metabolite_gene_co_occurancce_PubChem_clsa$entrez_ID)
class(metabolite_gene_co_occurancce_PubChem_clsa$PubChem_ID)<-as.character(metabolite_gene_co_occurancce_PubChem_clsa$PubChem_ID)
## for the metabolite with more than 1 Pubchem ID --> split they into separate rows
all_IVs_nearby_genes_formated_with_metabolite_ratios_5<-all_IVs_nearby_genes_formated_with_metabolite_ratios_4 %>% mutate(PUBCHEM_curated = strsplit(as.character(PUBCHEM_curated), ";")) %>% unnest(PUBCHEM_curated)
all_IVs_nearby_genes_formated_with_metabolite_ratios_5$PubChem_gene_pass<-apply(all_IVs_nearby_genes_formated_with_metabolite_ratios_5, 1, function(x) ifelse(x[5] %in% metabolite_gene_co_occurancce_PubChem_clsa[metabolite_gene_co_occurancce_PubChem_clsa$PubChem_ID == x[9],]$entrez_ID, "Y", "N"))
all_IVs_nearby_genes_formated_with_metabolite_ratios_5$PubChem_lit_evidence<-apply(all_IVs_nearby_genes_formated_with_metabolite_ratios_5, 1, function(x) 
  ifelse(x[15]== "Y", as.character(metabolite_gene_co_occurancce_PubChem_clsa[metabolite_gene_co_occurancce_PubChem_clsa$PubChem_ID == x[9] & metabolite_gene_co_occurancce_PubChem_clsa$entrez_ID ==x[5],]$evidence), "NA"))
all_IVs_nearby_genes_formated_with_metabolite_ratios_5$PubChem_related_gene<-apply(all_IVs_nearby_genes_formated_with_metabolite_ratios_5, 1, function(x) 
  ifelse(x[15]== "Y", as.character(metabolite_gene_co_occurancce_PubChem_clsa[metabolite_gene_co_occurancce_PubChem_clsa$PubChem_ID == x[9] & metabolite_gene_co_occurancce_PubChem_clsa$entrez_ID ==x[5],]$Gene_Symbol), "NA"))

all_IVs_nearby_genes_formated_with_metabolite_ratios_5$full_name<-apply(all_IVs_nearby_genes_formated_with_metabolite_ratios_5, 1, function(x) matabolites_annotation_July2[matabolites_annotation_July2$metabo_ID ==x[10],]$CHEMICAL_NAME[1])

write.csv(all_IVs_nearby_genes_formated_with_metabolite_ratios_5, "./analysis_dir/biological_relevance_check/metabolite_ratio_cojo_snp_metabolism_related_gene_checked.csv")
