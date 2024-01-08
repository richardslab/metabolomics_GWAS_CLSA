library(data.table)
library(LDlinkR)
library(tidyr)
library(stringr)
library(Matching)
library(dplyr)
library(stringr)
library("coloc")
library(arrow)

args<-commandArgs(TRUE)

## load the SNP-metabo-gene pairs with significant eQTL associations
clsa_metabolite_cojo_SNPs_filtered_SNPID_gene5_eqtl_ALL1<-read.csv("./results/clsa_metabolite_cojo_SNPs_filtered_SNPID_gene5_eqtl_ALL1.csv",row.names = 1)

### Metabolites GWAS: get the SNP info within +/-500kb of the leading snps using terminal(awk).<all_seg_metabo.txt> file contains all subsetted summary stats for all metabolites
metabo_gwas_dat<-fread("./restuls/coloc_results/all_seg_metabo.txt")
metabo_gwas_dat1<-metabo_gwas_dat%>%select(V4,V1:V3,V5:V14)
colnames(metabo_gwas_dat1)<-c("testID_fixed","testID","metabolite","snp","CHR","SNP_ID","POS_hg38","A1","A2","N","FREQ","BETA","SE","P")

metabo_gwas_dat_formatted<-copy(metabo_gwas_dat1)
## Flip the allele to make sure the allele frequency is for minor allele
metabo_gwas_dat_formatted$FREQ <- ifelse(metabo_gwas_dat1$FREQ < 0.5, metabo_gwas_dat1$FREQ, 1-metabo_gwas_dat1$FREQ)
metabo_gwas_dat_formatted$A1 <- ifelse(metabo_gwas_dat1$FREQ < 0.5, as.character(metabo_gwas_dat1$A1), as.character(metabo_gwas_dat1$A2))
metabo_gwas_dat_formatted$A2 <- ifelse(metabo_gwas_dat1$FREQ < 0.5, as.character(metabo_gwas_dat1$A2), as.character(metabo_gwas_dat1$A1))
metabo_gwas_dat_formatted$BETA <- ifelse(metabo_gwas_dat1$FREQ < 0.5, metabo_gwas_dat1$BETA, (-1)*metabo_gwas_dat1$BETA)

## the 3 information needed to load the eqtl summary stats files
coloc_result<-data.frame(
  testID_fixed=character(),
  mqtl=character(),
  metabo=character(),
  phenotype_id=character(),
  gene_symbol=character(),
  tissue=character(),
  nsnps=integer(),
  PP.H0.abf=numeric(),
  PP.H1.abf=numeric(),
  PP.H2.abf=numeric(),
  PP.H3.abf=numeric(),
  PP.H4.abf=numeric())

i=args[1]
set<-clsa_metabolite_cojo_SNPs_filtered_SNPID_gene5_eqtl_ALL1[i,]
testID_fixed_test<-as.character(set$testID_fixed)
mqtl<-set$SNP.x
metabo<-set$metabolite
phenotype_id_sig<-set$phenotype_id
gene_id<-set$gene_id
gene_symbol<-set$gene_symbol
tissue<-set$tissue
pos_hg38<-set$POS_hg38
pos_hg38_lower<-ifelse(pos_hg38-500000>0,pos_hg38-500000,0)
pos_hg38_upper<-pos_hg38+500000
file_name<-paste0(set$tissue,".v8.EUR.allpairs.chr",set$CHR,".parquet")
df<- read_parquet(paste0("./data/GTEx_Analysis_v8_QTLs/GTEx_Analysis_v8_EUR_eQTL_all_associations/",file_name), as_tibble = TRUE)
df_sub<-df%>%filter(phenotype_id == phenotype_id_sig)
df_sub_in_range<-df_sub %>% rowwise()%>%mutate(POS_hg38=as.integer(strsplit(variant_id,"_")[[1]][2]),
                                               SNP_ID_1=paste0(strsplit(variant_id,"_")[[1]][1],":",strsplit(variant_id,"_")[[1]][2],":",strsplit(variant_id,"_")[[1]][3],":",strsplit(variant_id,"_")[[1]][4]),
                                               SNP_ID_2=paste0(strsplit(variant_id,"_")[[1]][1],":",strsplit(variant_id,"_")[[1]][2],":",strsplit(variant_id,"_")[[1]][4],":",strsplit(variant_id,"_")[[1]][3]))%>%
  filter(POS_hg38>=pos_hg38_lower & POS_hg38<pos_hg38_upper)
##metabo summary stats
metabo_gwas_dat_formatted_sub<-metabo_gwas_dat_formatted%>%filter(testID_fixed == testID_fixed_test)
## merge by positions
summary_stats_merged<-inner_join(metabo_gwas_dat_formatted_sub,df_sub_in_range, by=c("POS_hg38"="POS_hg38"))
summary_stats_merged_filtered<-summary_stats_merged%>%filter((SNP_ID == SNP_ID_1 | SNP_ID == SNP_ID_2))%>%filter(is.na(slope)==F & is.na(slope_se)==F & is.na(maf)==F & is.na(pval_nominal)==F)%>%
  filter(FREQ>0.05 & maf>0.05)
summary_stats_merged_filtered$SNP_ID<-as.character(summary_stats_merged_filtered$SNP_ID)
## run coloc
d1=summary_stats_merged_filtered[,1:14]
d2=summary_stats_merged_filtered[,c(6,15:26)]
n_metabo<-max(d1$N)
n_eqtl<-max(d2$ma_count)
coloc.res <- coloc.abf(dataset1=list(pvalues=d1$P, 
                                     beta=d1$BETA, 
                                     varbeta=(d1$SE)^2,
                                     snp=d1$SNP_ID, 
                                     MAF=d1$FREQ, 
                                     N=n_metabo,
                                     type="quant"), 
                       dataset2=list(beta=d2$slope, 
                                     varbeta=(d2$slope_se)^2, 
                                     snp=d2$SNP_ID, 
                                     MAF=d2$maf,
                                     pvalues=d2$pval_nominal,
                                     N=n_eqtl,
                                     type="quant"))
temp_result<-data.frame(testID_fixed_test, mqtl,metabo,phenotype_id_sig,gene_id,
                        gene_symbol,tissue,t(coloc.res$summary))

write.csv(summary_stats_merged_filtered,paste0("./restuls/coloc_results/eqtl/summary_stats/",i,"_summary_stats_",testID_fixed_test,"_",gene_symbol,"_",tissue,"_",".csv"))
write.csv(temp_result,paste0("./restuls/coloc_results/eqtl/",i,"_coloc_result_set",testID_fixed_test,"_",gene_symbol,"_",tissue,"_",".csv"))
