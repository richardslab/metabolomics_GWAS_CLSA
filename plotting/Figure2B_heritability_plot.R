library(data.table)
library(ggplot2)
library(dplyr)

setwd("/home/richards/yiheng.chen/scratch/project1_2_metabolomics_GWAS_CLSA/codes/plotting")
metabolite_variance_explained_results_summary<-fread("metabolite_variance_explained_results_summary.txt")
colnames(metabolite_variance_explained_results_summary)<-c("metabolite","test","variance_explained_by_genetics","se")
clsa_metabolite_cojo_SNPs_filtered<-read.csv("clsa_metabolites_cojo_formatted_rsid_formatted_filtered_May2022.csv",row.name=1)

### annotate the metabolites
matabolites_annotation_July2<-read.csv("CLSA_COMBINED_ANNOTATIONTABLEALL_v1_modi_July2021.csv")
metabo_superpathway_list<-matabolites_annotation_July2 %>% rowwise() %>% mutate(SUPER_PATHWAY_mod=ifelse(SUPER_PATHWAY == "", "Unknown", SUPER_PATHWAY)) %>% select(metabo_ID, SUPER_PATHWAY_mod) %>% distinct()
metabolite_variance_explained_results_summary$full_name <- apply(metabolite_variance_explained_results_summary, 1, function(x) matabolites_annotation_July2$CHEMICAL_NAME[matabolites_annotation_July2$metabo_ID == x[1]][1])
metabolite_variance_explained_results_summary$SUPER_PATHWAY_mod <- apply(metabolite_variance_explained_results_summary, 1, function(x) metabo_superpathway_list$SUPER_PATHWAY_mod[metabo_superpathway_list$metabo_ID == x[1]][1])

## remove the Unknown and partially characterized molecules
metabolite_variance_explained_results_summary_filtered1<-metabolite_variance_explained_results_summary %>% group_by(SUPER_PATHWAY_mod) %>% mutate(max_heri=max(variance_explained_by_genetics), bar_width= 0.02*n())
group.colors <- c(`Amino Acid` = "#fc8d62", Lipid = "#66c2a5",`Cofactors and Vitamins`="#e78ac3",Xenobiotics="#1f78b4", Nucleotide="#ffd92f", 
                  Carbohydrate="#8da0cb",Peptide="#e5c494", Energy="#a6d854", `Unknown`="#bc80bd", `Partially Characterized Molecules`="#fdb462")
metabolite_variance_explained_results_summary_filtered1$SUPER_PATHWAY_mod<-factor(metabolite_variance_explained_results_summary_filtered1$SUPER_PATHWAY_mod, 
                                                                                  levels=c("Amino Acid", "Lipid","Cofactors and Vitamins","Xenobiotics", "Nucleotide", "Carbohydrate","Peptide", "Energy", "Unknown","Partially Characterized Molecules"))

## do not plot unknown or partially characterized metabolites
metabolite_variance_explained_results_summary_filtered2<-metabolite_variance_explained_results_summary_filtered1%>%filter(!(SUPER_PATHWAY_mod %in% c("Partially Characterized Molecules", "Unknown")))

p_violin<-ggplot(metabolite_variance_explained_results_summary_filtered2, aes(x=reorder(SUPER_PATHWAY_mod, -max_heri), variance_explained_by_genetics, fill=SUPER_PATHWAY_mod)) +
  geom_violin(scale = "count",width = 3) +
  scale_fill_manual(values=group.colors)+
  theme_bw() +
  theme(axis.text.x = element_text(size=26, angle = 20, hjust=1.2,colour = "black"), axis.text.y = element_text(size=25,colour = "black"),
        axis.title=element_text(size=30),  legend.text=element_text(size=20), legend.title = element_text(size=25),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black"))+
  labs(x="Metabolites ranked by estimated heritability", y="Heritability explained by assayed genetic variants", size=100)+
  guides(fill = guide_legend(title = "Super Pathways", ncol=1))+
  stat_summary(fun = "median",
               geom = "crossbar", 
               width = 0.5,
               colour = "red") +
  geom_hline(yintercept=0.199774, linetype=2, col = 'blue', size=1.2)

tiff("Fig2B_metabo_hertiability_by_superpathway_violin_v4_median_Oct162022.tiff", width = 1400, height = 880, units = "px")
p_violin
dev.off()
