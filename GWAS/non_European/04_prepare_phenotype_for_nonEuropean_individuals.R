library(data.table)
library(dplyr)
library(tidyverse)

setwd("./analysis_dir")
CLSA_NORMDATAALL<-read.csv("./clsa_phenotypes/clsa_batch_norm_metabo_dat",na.strings=c("","NA"," ","Metabolite_not_called_in_this_set", "NaN"))

## metabolomics ADM_GWAS_ID
Metabolon_client_identifier<-read.csv("./clsa_phenotypes/linking_file")
## baseline data
Baseline_data_v5<-read.csv("./clsa_phenotypes/clsa_baseline_data", header=TRUE,stringsAsFactors=FALSE)
##metabolites annotation
matabolites_annotation<-read.csv("./clsa_phenotypes/clsa_metabo_annotation_dat")
#metabolomics metadata
ID_metadat<-read.csv("./clsa_phenotypes/clsa_metabo_meta_dat")
##subset baseline data
Baseline_data_v5_sub<-Baseline_data_v5 %>% select("ADM_GWAS3_COM","entity_id","SEX_ASK_COM", "AGE_NMBR_COM","BLD_ALC24_HR_COM","BLD_ALC24_MIN_COM", "BLD_FD24_HR_COM", "BLD_FD24_MIN_COM")
colnames(Baseline_data_v5_sub)<-c("ADM_GWAS_COM","entity_id","SEX_ASK_COM", "AGE_NMBR_COM","BLD_ALC24_HR_COM","BLD_ALC24_MIN_COM", "BLD_FD24_HR_COM", "BLD_FD24_MIN_COM")

##remove the row without ADM_GWAS_COM
Baseline_data_v5_sub_withADMid<-Baseline_data_v5_sub[is.na(Baseline_data_v5_sub$ADM_GWAS_COM)==FALSE,]
##merge data
dat1<-merge(CLSA_NORMDATAALL, Metabolon_client_identifier, by="PARENT_SAMPLE_NAME")
dat2<-left_join(dat1, Baseline_data_v5_sub_withADMid, by="ADM_GWAS_COM")
dat3<-left_join(dat2, ID_metadat, by="PARENT_SAMPLE_NAME")

##### remove the metabolites have no measurements in >50% of samples
## get the list of metabolites with > 50% missing measurements
dat3_filtered_metaboMeasurement<-dat3 %>% select(X35:X999926111)
metabolite_list<-c()
for(i in names(dat3_filtered_metaboMeasurement)){
  if(sum(is.na(dat3_filtered_metaboMeasurement[[i]]))/length(dat3_filtered_metaboMeasurement[[i]])<0.5){
    metabolite_list<-append(metabolite_list,i)
  }
} ## 1091 metabolites pass this filter
dat4<-dat3 %>% select(PARENT_SAMPLE_NAME,metabolite_list,BLD_FD24_HR_COM,ADM_GWAS_COM)

#### check if there is any individual with over 50% missing
dat4_filtered<-dat4 %>% mutate(missingMeasurement_perce = rowSums(is.na(across(X35:X999926111)))/1091) %>% filter(missingMeasurement_perce<0.5) ## all the individual have >50% measurements

#### get the individual with hour to last meal/drink information
dat5 <-dat4_filtered %>% filter(BLD_FD24_HR_COM>=0)

################################# genomic information
fam_dat<-fread("./clsa_imputed_genotype_data/clsa_gen_v3.fam")
sample_dat<-fread("./clsa_imputed_genotype_data/clsa_imp_v3.sample")
sqc_file<-fread("./clsa_imputed_genotype_data/clsa_sqc_v3.txt")

################ merge metabolomics data and genomic information
### list of individuals without metabolomics measurement or did not pass previous filter
ID_list<-data.frame(cbind(Metabolon_client_identifier$ADM_GWAS_COM, Metabolon_client_identifier$ADM_GWAS_COM))
## get ID for individuals pass previous filter
ID_list_compl<-ID_list[complete.cases(ID_list),]
colnames(ID_list_compl)<-c("FID","IID")
ID_list_no_metabo<-sample_dat %>% select(ID_1, ID_2) %>% filter(!(ID_1 %in% ID_list_compl$FID) | !(ID_1 %in% dat5$ADM_GWAS_COM )) 
colnames(ID_list_no_metabo)<-c("FID","IID") ### 17445 individuals can be removed from GWAS

## prepare a new sample file
sample_dat_sub<-sample_dat %>% filter((ID_1 %in% ID_list_compl$FID) & (ID_1 %in% dat5$ADM_GWAS_COM ))

#### prepare covariates data for gwas (non-european indiv)
pop_matching_dat<-data.frame(index=c(1,2,3),pop=c("south_asian","east_asian", "black"))

work_dir="../GWAS_related_results/metabolites_phenotype_nonEurop"

###### generate PCs for specific Ancestry
for (i in 1:3){
   pop=pop_matching_dat$pop[i]
   ##read in the ancestry specific PC results
   pc_path=paste0(work_dir,pop,"/genotype/clsa_gen_v3_pruned_forPC_",pop,"_pca.eigenvec")
   pca_res<-read_table(pc_path, col_names=F)
   pca_res <- pca_res[,-1]
   names(pca_res)[1] <- "ind"
   names(pca_res)[2:ncol(pca_res)] <- paste0("PC", 1:(ncol(pca_res)-1))
   ##prepare covariables for fastGWA
   covarCol_data<-sqc_file_all_ancestries%>%filter(pca.cluster.id == i)
   covarCol_data1<-left_join(covarCol_data, Baseline_data_v5_sub[,c(1,4)])#add age info
   covarCol_data2<-left_join(covarCol_data1, pca_res, by=c("ADM_GWAS_COM"="ind"))
   dat5_sub<-dat5%>%select(ADM_GWAS_COM, BLD_FD24_HR_COM)## select BLD_FD24_HR_COM, BLD_FD24_HR_COM columns
   covarCol_data3<-left_join(covarCol_data2,dat5_sub) 
   covarCol_data3_sub<-covarCol_data3%>%select(ADM_GWAS_COM, batch, chromosomal.sex, pca.cluster.id, PC1:PC5, Sex, AGE_NMBR_COM, BLD_FD24_HR_COM)
   covarCol_data3_sub1<-cbind(covarCol_data3_sub$ADM_GWAS_COM, covarCol_data3_sub)
   colnames(covarCol_data3_sub1)<-c("FID","IID","batch","chromosomal.sex","pca.cluster.id","PC1", "PC2","PC3","PC4","PC5","Sex","Age", "BLD_FD24_HR_COM")
   qcovarCol_clsa_pop<-covarCol_data3_sub1 %>% select(FID, IID, PC1, PC2, PC3, PC4, PC5, Age, BLD_FD24_HR_COM)
   CovarCol_clsa_pop<-covarCol_data3_sub1 %>% select(FID, IID, batch, Sex)
   write.table(qcovarCol_clsa_pop,paste0(work_dir, pop,"/phenotype/qcovarCol_clsa_",pop,"_noHeader.txt"),sep="\t",quote = FALSE, row.names = FALSE, col.names = FALSE)
   write.table(CovarCol_clsa_pop,paste0(work_dir, pop,"/phenotype/CovarCol_clsa_",pop,"_noHeader.txt"),sep="\t",quote = FALSE, row.names = FALSE, col.names = FALSE)
   ##check the metabolites present in over 50% indiv of each ancestry
   dat6 <-dat5%>% filter(ADM_GWAS_COM %in% covarCol_data$ADM_GWAS_COM)%>%select(metabolite_list)
   metabolite_list_totest<-c()
   for(i in names(dat6)){
     if(sum(is.na(dat6[[i]]))/length(dat6[[i]])<0.5){
       metabolite_list_totest<-append(metabolite_list_totest,i)
     }
   }
   ###prepapre phenotype files
   post_outlier_removal_mean_sd<-data.frame(
       metabolite=character(),
       n_metabo_nonNA=integer(),
       n_metabo_post_3sd=integer(),
       mean=numeric(),
       sd=numeric(),
       standardized_mean=numeric(),
       standardized_sd=numeric())
   for (metabo in metabolite_list_totest) {
          dat_sub<-dat5 %>% filter(ADM_GWAS_COM %in% covarCol_data$ADM_GWAS_COM) %>% select(ADM_GWAS_COM,metabo)
          dat_sub$log_metabo<-apply(dat_sub, 1, function(x) log(as.numeric(x[2])))##natural log transform the data
          ##remove the measurements that 3 SD away
          mean_metabo<-mean(dat_sub$log_metabo,na.rm = TRUE)
          sd_metabo<-sd(dat_sub$log_metabo,na.rm = TRUE)
          dat_sub_filter<-dat_sub %>% filter(abs((log_metabo-mean_metabo)/sd_metabo)<=3)
          #check n indiv for each meetabolites
          dat_sub_filter1<-dat_sub_filter%>%filter(is.na(log_metabo)==F)
          n_metabo_before_filter<-dat_sub%>%filter(is.na(log_metabo)==F)%>%nrow()
          n_metabo<-nrow(dat_sub_filter1)
          ##standardize data
          mean_metabo_post<-mean(dat_sub_filter$log_metabo, na.rm = TRUE)
          sd_metabo_post<-sd(dat_sub_filter$log_metabo, na.rm = TRUE)
          dat_sub_filter$std_metabo<-apply(dat_sub_filter, 1, function(x) (as.numeric(x[3])-mean_metabo_post)/sd_metabo_post)
          phenoCol<-left_join(covarCol_data3_sub1[,1:2],dat_sub_filter[,c(1,4)],by=c("FID"="ADM_GWAS_COM"))
          colnames(phenoCol)<-c("FID","IID","Standardized_metabolites_level")
          standardized_mean<-mean(phenoCol$Standardized_metabolites_level, na.rm = T)
          standardized_sd<-sd(phenoCol$Standardized_metabolites_level, na.rm = T)
          tmp_dat<-data.frame(metabo, n_metabo_before_filter, n_metabo, mean_metabo_post, sd_metabo_post, standardized_mean, standardized_sd)
          colnames(tmp_dat)<-c("metabolite", "n_metabo_nonNA","n_metabo_post_3sd","mean","sd", "standardized_mean", "standardized_sd")
          post_outlier_removal_mean_sd<-rbind(post_outlier_removal_mean_sd,tmp_dat)
          write.table(phenoCol, paste0(work_dir,pop,"/phenotype/metabolites/phenoCol_",metabo,"_clsa_",pop,"_3sd.txt"),sep="\t",quote = FALSE,row.names = FALSE, col.names = FALSE)
  }
  write.csv(post_outlier_removal_mean_sd, paste0(work_dir,pop,"/phenotype/",pop,"_post_outlier_metabolites_removal_mean_sd.csv"))
}
