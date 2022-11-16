library(tidyverse)
library(data.table)
library(dplyr)
library(ggplot2)
library(tidyr)
library(RColorBrewer)

setwd("./analysis_dir/plotting")
matabolites_annotation_July2<-read.csv("CLSA_COMBINED_ANNOTATIONTABLEALL_v1_modi_July2021")
metabo_superpathway_list<-matabolites_annotation_July2 %>% rowwise() %>% mutate(SUPER_PATHWAY_mod=ifelse(SUPER_PATHWAY == "", "Unknown", SUPER_PATHWAY)) %>% select(metabo_ID, SUPER_PATHWAY_mod) %>% distinct()

clsa_metabolite_cojo_SNPs_filtered1 <- read.csv("clsa_metabolites_cojo_formatted_rsid_formatted_filtered_May2022.csv", row.names = 1)

##### mahattan plots for the SNPs around the leading snps 
##get the SNP info within +/-500kb of the leading snps using terminal(awk) -- both metabolite and metabolite ratios
metabo_gwas_dat<-fread("all_seg_metabo.txt")
colnames(metabo_gwas_dat)<-c("testID", "metabolite","SNP","testID_fixed", "CHR", "SNP_ID","POS","A1","A2","N","FREQ","BETA","SE","P")
metabo_gwas_dat_metabolite <- metabo_gwas_dat %>% filter(metabolite %in% clsa_metabolite_cojo_SNPs_filtered1$metabolite)
metabo_gwas_dat_metabolite1<-left_join(metabo_gwas_dat_metabolite, metabo_superpathway_list, by=c("metabolite"="metabo_ID"))
metabo_gwas_dat_metabolite1_filtered <-metabo_gwas_dat_metabolite1 %>% filter(FREQ > 0.01) %>% mutate(P.fix=ifelse(P==0, 5e-324, P)) %>% select(SNP_ID, CHR, POS, P.fix, SUPER_PATHWAY_mod)

#### read in a GWAS for background (grey scale peaks)
background_gwas<-fread("X999925957_EUR_GWAS.fastGWA.gz")
background_gwas_filtered<-background_gwas %>%filter(AF1 > 0.05 & !(SNP %in% metabo_gwas_dat_metabolite1_filtered$SNP_ID)) %>% select(SNP, CHR, POS, P) %>% rowwise() %>% mutate(SUPER_PATHWAY_mod = NA)
colnames(background_gwas_filtered)<-c("SNP_ID", "CHR", "POS", "P.fix", "SUPER_PATHWAY_mod")

### combine the background gwas and leading snp gwas
gwas_together<-rbind(metabo_gwas_dat_metabolite1_filtered, background_gwas_filtered)

## prepare a common X axis
don_all_metabolite <- gwas_together %>% 
  # Compute chromosome size
  group_by(CHR) %>% 
  summarise(chr_len=max(POS)) %>% 
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(as.numeric(chr_len))-chr_len) %>%
  select(-chr_len) %>%
  # Add this info to the initial dataset
  left_join(gwas_together, ., by=c("CHR"="CHR")) %>%
  # Add a cumulative position of each SNP
  arrange(CHR, POS) %>%
  mutate(BPcum=POS+tot)

axisdf <- don_all_metabolite %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

color_match=c("Lipid"="#66c2a5", "Amino Acid"="#fc8d62", "Carbohydrate"="#8da0cb",
              "Cofactors and Vitamins"="#e78ac3", "Energy"="#a6d854", "Nucleotide"="#ffd92f",
              "Peptide"="#e5c494","Xenobiotics"="#1f78b4", "Unknown"="#bc80bd", "Partially Characterized Molecules"="#fdb462")

p1=ggplot(data=don_all_metabolite, aes(x=BPcum, y=-log10(P.fix))) +
  # Show all points
  geom_point(data=subset(don_all_metabolite, P.fix>5e-8), aes(color=as.factor(CHR)), alpha=0.8, size=2) +
  scale_color_manual(values = rep(c("#636363", "#bdbdbd"), 11)) +
  new_scale_color() +
  geom_point(data=subset(don_all_metabolite, P.fix<5e-8), aes(color=as.factor(SUPER_PATHWAY_mod)), alpha=0.8, size=2) +
  scale_color_manual(values = color_match) +
  # custom X axis:
  scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center) +
  scale_y_continuous(breaks = seq(0, 400, 50),  label=seq(0, 400, 50)) +   # remove space between plot area and x axis
  # Custom the theme:
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    plot.margin = margin(2, 1, 1, 1, "cm"),
    axis.text.x = element_text(size=18),
    axis.text.y = element_text(size=18),
    axis.title=element_text(size=30),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank())+
  labs(x="Chromosome", size=65) +
  ylab(expression(-log[10](P)))

tiff("Fig1_clsa_metabo_cojo_gws_plot_all.tiff", units="px", width=2800, height=700)
p1
dev.off()
