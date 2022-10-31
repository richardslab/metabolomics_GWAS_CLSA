library(data.table)
library(ggplot2)
library(dplyr)
library(ggpubr)

setwd("/home/richards/yiheng.chen/scratch/project1_2_metabolomics_GWAS_CLSA/codes/plotting")
metabolite_variance_explained_results_summary<-fread("metabolite_variance_explained_results_summary.txt")
colnames(metabolite_variance_explained_results_summary)<-c("metabolite","test","variance_explained_by_genetics","se")

clsa_metabolite_cojo_SNPs_filtered1 <- read.csv("clsa_metabolite_cojo_SNPs_with_loci_May2022.csv", row.names = 1)

### annotate the metabolites (variance explained dataset)
matabolites_annotation_July2<-read.csv("CLSA_COMBINED_ANNOTATIONTABLEALL_v1_modi_July2021.csv")
metabo_superpathway_list<-matabolites_annotation_July2 %>% rowwise() %>% mutate(SUPER_PATHWAY_mod=ifelse(SUPER_PATHWAY == "", "Unknown", SUPER_PATHWAY)) %>% select(metabo_ID, SUPER_PATHWAY_mod) %>% distinct()
metabolite_variance_explained_results_summary$full_name <- apply(metabolite_variance_explained_results_summary, 1, function(x) matabolites_annotation_July2$CHEMICAL_NAME[matabolites_annotation_July2$metabo_ID == x[1]][1])
metabolite_variance_explained_results_summary$SUPER_PATHWAY_mod <- apply(metabolite_variance_explained_results_summary, 1, function(x) metabo_superpathway_list$SUPER_PATHWAY_mod[metabo_superpathway_list$metabo_ID == x[1]][1])
## keep the metabolites with independent associations
metabolite_variance_explained_results_summary_filtered1<-metabolite_variance_explained_results_summary %>% filter(metabolite %in% clsa_metabolite_cojo_SNPs_filtered1$metabolite)

#### get the number of associated SNPs for each metabolites
#####################
CLSA_LeadSNPs_dat_sub<-clsa_metabolite_cojo_SNPs_filtered1 %>%select(SNP, A1, A2,FREQ, BETA, SE, P,N,CHR,POS,metabolite,SNP_ID,locus)
colnames(CLSA_LeadSNPs_dat_sub) <- c("SNP", "A1", "A2","FREQ", "BETA", "SE", "P","N","CHR","POS","metabolite","SNP_ID","locus")
CLSA_LeadSNPs_dat_sub$P<-as.numeric(CLSA_LeadSNPs_dat_sub$P)

### plot n metabolites per snp
n_metaoblite_per_snp<-CLSA_LeadSNPs_dat_sub %>% select(locus,metabolite)%>% group_by(locus)%>%summarise(n_metabolite=length(unique(metabolite)))
n_loci_per_metabo<-CLSA_LeadSNPs_dat_sub %>% select(locus,metabolite)%>% group_by(metabolite)%>%summarise(n_loci=length(unique(locus)))
hist(n_loci_per_metabo$n_loci)
median(n_loci_per_metabo$n_loci)#2

n_loci_per_metabo1<-left_join(n_loci_per_metabo,metabo_superpathway_list, by=c("metabolite"="metabo_ID"))
#n_snp_per_metabo1_filtered<-n_snp_per_metabo1 %>%filter(SUPER_PATHWAY_mod !="Unknown" & SUPER_PATHWAY_mod != "Partially Characterized Molecules")

metabolite_variance_explained_n_loci<-left_join(metabolite_variance_explained_results_summary_filtered1, n_loci_per_metabo1, by="metabolite")

###### final revision (remove background grids)
p_herta_n_loci_all<-ggplot(metabolite_variance_explained_n_loci, aes(x=n_loci, y=variance_explained_by_genetics)) +
  geom_point(position = "jitter", size=5, color="grey") +
  geom_smooth(method = "lm", aes(x=n_loci, y=variance_explained_by_genetics),level=0.95) +
  scale_x_continuous(breaks=c(0, 2, 4, 6, 8,10))+
  theme_bw()+
  labs(x="Number of loci associated per metabolite", y="Heritability explained by assayed genetic variants")+
  theme(axis.text.x = element_text(size=25, angle = 20, hjust=1,colour = "black"), axis.text.y = element_text(size=25,colour = "black"),
        axis.title=element_text(size=30),  legend.text=element_text(size=20),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black"))+
  stat_cor(method = "spearman", label.x = 6, label.y = 0.8, size=10)


tiff("Fig2C_heritability_n_loci_metabolite_all_metabolites_Oct162022.tiff", width = 850, height =750, units = "px")
p_herta_n_loci_all
dev.off()

##get exact p-value for correlation
cor.test(metabolite_variance_explained_n_loci$n_loci, metabolite_variance_explained_n_loci$variance_explained_by_genetics, method="spearman") -> tmp
tmp$estimate
#0.3584618 
tmp$p.value
#2.383618e-22


