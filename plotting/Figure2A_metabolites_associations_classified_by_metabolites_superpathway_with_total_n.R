library(dplyr)
library(data.table)
library(ggplot2)
setwd("./analysis_dir/plotting")
matabolites_annotation_July2<-read.csv("CLSA_COMBINED_ANNOTATIONTABLEALL_v1_modi_July20201.csv")
metabo_superpathway_list<-matabolites_annotation_July2 %>% rowwise() %>% mutate(SUPER_PATHWAY_mod=ifelse(SUPER_PATHWAY == "", "Unknown", SUPER_PATHWAY)) %>% select(metabo_ID,  CHEMICAL_NAME, SUPER_PATHWAY_mod) %>% distinct()

####### plot the proportion of metabolites tested by SUPER-PATHWAY
metabolite_list_totest1<-read.csv("metabolites_to_be_gwas.csv")
colnames(metabolite_list_totest1)<-c("id","metabolite")
metabolite_list_totest2<-left_join(metabolite_list_totest1, metabo_superpathway_list, by=c("metabolite"="metabo_ID"))
n_superpathway_tested_gwas<-metabolite_list_totest2 %>% select(SUPER_PATHWAY_mod,metabolite) %>% group_by(SUPER_PATHWAY_mod) %>% summarise(n_superpathway = length(unique(metabolite))) 

####### plot the proportion of metabolites that have novel associations (independent)
### read in the summary of association checked 
summary_of_checked_associations<-read.csv("clsa_metabolite_cojo_SNPs_formatted_rsid_filtered_novelty_checked_bonf.csv", row.names = 1)

####### plot the proportion of metabolites have independent significant IVs 
clsa_metabolite_cojo_SNPs_filtered1<-left_join(summary_of_checked_associations, metabo_superpathway_list, by=c("metabolite"="metabo_ID"))
n_superpathway_with_IVs<-clsa_metabolite_cojo_SNPs_filtered1 %>%select(SUPER_PATHWAY_mod,metabolite) %>% group_by(SUPER_PATHWAY_mod) %>% summarise(n_superpathway = length(unique(metabolite))) 
sum(n_superpathway_with_IVs$n_superpathway)#690 metabolites have associations

## annotate novel metabolites with novel IVs (for associations pass bonferroni correction)
n_superpathway_with_novel_IVs<-clsa_metabolite_cojo_SNPs_filtered1 %>% filter(pass_bonferroni ==1 &novelty_check == "novel")%>%select(SUPER_PATHWAY_mod,metabolite) %>% group_by(SUPER_PATHWAY_mod) %>% summarise(n_superpathway = length(unique(metabolite)))
sum(n_superpathway_with_novel_IVs$n_superpathway)#313 metabolites have novel associations

#### join all number
#n_superpathway_tested_gwas: all metabolites (in each pathway-tested in GWAS) -->n_superpathway.x
#n_superpathway_with_IVs: all metabolites with GWS IVs -->n_superpathway.y
#n_superpathway_with_novel_IVs: all metabolites with novel GWS IVs -->n_superpathway
all_superpathway_n<-left_join(n_superpathway_tested_gwas, n_superpathway_with_IVs, by=c("SUPER_PATHWAY_mod"))
all_superpathway_n<-left_join(all_superpathway_n, n_superpathway_with_novel_IVs, by=c("SUPER_PATHWAY_mod"))
all_superpathway_n<-all_superpathway_n %>% rowwise() %>% mutate(percent_metabo_with_IV = n_superpathway.y/n_superpathway.x*100,
                                                                super_path_n_metabo= paste0(SUPER_PATHWAY_mod, " (", n_superpathway.x ,")"))

colnames(all_superpathway_n)<-c("SUPER_PATHWAY", "N_testd_metabolites_in_GWAS", "N_metabolites_with_GWS_SNPs", "N_metabolites_with_novel_SNPs", "percent_of_tested_metabolites_with_GWS_SNPs", "super_pathway(n)")### for a percentage stacking bar

all_superpathway_n1<-all_superpathway_n %>% rowwise() %>% mutate(proportion_metabo_with_novel_IVs = N_metabolites_with_novel_SNPs ,
                                                                 proportion_metabo_with_no_IVs = N_testd_metabolites_in_GWAS - N_metabolites_with_GWS_SNPs,
                                                                 proportion_metabo_witout_novel_IV = N_metabolites_with_GWS_SNPs - proportion_metabo_with_novel_IVs)

all_superpathway_n1_melt<-melt(all_superpathway_n1[,c(6,7:9)])
all_superpathway_n1_melt1<-left_join(all_superpathway_n1_melt, all_superpathway_n1[,c(6,2)], by=c("super_pathway(n)"))
colnames(all_superpathway_n1_melt1)<-c("super_pathway(n)","Class","value","N_testd_metabolites_in_GWAS")

all_superpathway_n1_melt1$Class <- factor(all_superpathway_n1_melt1$Class, levels = c("proportion_metabo_with_no_IVs", "proportion_metabo_witout_novel_IV", "proportion_metabo_with_novel_IVs"))

p_bar<-ggplot(all_superpathway_n1_melt1, aes(fill=Class, x=reorder(`super_pathway(n)`, -N_testd_metabolites_in_GWAS), y=value))+ 
  geom_bar(position="stack", stat="identity")+
  theme_bw() +
  theme(axis.text.x = element_text(size=25, angle = 20, hjust=1,colour = "black"), axis.text.y = element_text(size=25,colour = "black"),
        axis.title=element_text(size=30), legend.text=element_text(size=25), legend.position = "top", legend.title=element_text(size=25),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black"))+
  scale_fill_manual(labels = c("Metabolites without associations", "Metabolites with associations (not novel)","Metabolites with novel associations"), values=c("#b2df8a","#1f78b4","#a6cee3"))+
  labs(x="Super pathways (N metabolites tested in GWAS)", y="Counts", size=10)+
  guides(fill = guide_legend(title = "Classes", title.position = "top", ncol=1))

tiff("Fig2A.tiff", width = 1050, height = 1000, units = "px")
p_bar
dev.off()
