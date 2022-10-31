library(ggplot2)
library(dplyr)
library(stringr)
library(data.table)

setwd("/home/richards/yiheng.chen/scratch/project1_2_metabolomics_GWAS_CLSA/codes/plotting")
clsa_metabolite_cojo_SNPs_formatted_rsid_filtered<-read.csv("clsa_metabolite_cojo_SNPs_with_loci_May2022.csv",row.name=1)
clsa_metabolite_ratios_cojo_SNPs_formatted_rsid_filtered<-read.csv("clsa_metabolite_ratios_cojo_SNPs_formatted_rsid_bonf_filtered_enzyme_transporter_with_locus_May2022.csv",row.name=1)

mr_result<-read.csv("mr_results_for_plot_with_beta_or.csv")

mr_result_pleiotrpy_check<-read.csv("mr_results_pleiotropy_reverse_causation_check.csv")
mr_result_pleiotrpy_check_sub<-mr_result_pleiotrpy_check%>%select(Outcomes, `Metabolite.or.metabolite.ratios`,pass_metabo_pleiotropy_check, pass_reverse_causation_check)%>%distinct()
mr_result1<-left_join(mr_result, mr_result_pleiotrpy_check_sub, by=c("full_name"="Metabolite.or.metabolite.ratios", "traits"="Outcomes"))

mr_result_filtered<-mr_result1%>%filter(pass_metabo_pleiotropy_check =="Y" & pass_reverse_causation_check =="Y")
mr_result_metabo<-mr_result_filtered%>%filter(full_name %in% clsa_metabolite_cojo_SNPs_formatted_rsid_filtered$full_name)
mr_result_metabo_ratio<-mr_result_filtered%>%filter(full_name %in% clsa_metabolite_ratios_cojo_SNPs_formatted_rsid_filtered$full_name)

unique_trait_list<-unique(mr_result_filtered$traits)

## for the immune-related traits
mr_result_filtered_immune <- mr_result_filtered[mr_result_filtered$Category=="Immune-related",]
mr_result_filtered_immune$name_plot <- paste0(mr_result_filtered_immune$full_name,"__",mr_result_filtered_immune$traits)

p1<-ggplot(mr_result_filtered_immune, aes(x = log(Beta_Odds.ratio), y = reorder(name_plot,-Beta_Odds.ratio))) +
  geom_errorbarh(aes(xmin = log(CIL), xmax = log(CIU)), color = "#1b9e77", height = 0, size=1.2) +
  geom_point(color = "#1b9e77",size=2.5) +
  facet_grid(traits~.,scales = "free",space = "free_y") +
  theme_bw() +
  theme(axis.text = element_text(size = 15,color="black"),
        axis.title = element_text(size = 16),
        strip.text = element_text(size = 17,face="bold"),
        axis.text.y=element_text(size=16,color="black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black")) +
  geom_vline(xintercept = 0,col = "grey",lty = 2) +
  ylab("") +
  scale_y_discrete(labels = function(x) unlist(lapply(str_split(x,"__"), function(y) stringr::str_wrap(y[1], width = 30))))+
  xlab(stringr::str_wrap('OR (95% CI) per unit* increase of exposure', width = 40)) +
  scale_x_continuous(breaks = c(log(0.2),log(0.5),log(1),log(2),log(5),log(10)),
                     labels = c(0.2,0.5,1,2,5,10))

tiff("mr_result_filtered_immune_Oct162022.tiff", width = 1300, height = 2700, units = "px", res=150)
p1
dev.off()

## for the metabolism-related traits (OR trait)
mr_result_filtered_metabo <- mr_result_filtered[mr_result_filtered$Category=="Metabolism-related" & mr_result_filtered$traits !="BMI",]
mr_result_filtered_metabo$name_plot <- paste0(mr_result_filtered_metabo$full_name,"__",mr_result_filtered_metabo$traits)
tiff("mr_result_filtered_metabolism_OR_Oct162022.tiff", width = 1300, height = 1000, units = "px", res=150)
ggplot(mr_result_filtered_metabo, aes(x = log(Beta_Odds.ratio), y = reorder(name_plot,-Beta_Odds.ratio))) +
  geom_errorbarh(aes(xmin = log(CIL), xmax = log(CIU)), color = "#d95f02", height = 0, size=1.2) +
  geom_point(color = "#d95f02", size=2.5) +
  facet_grid(traits~.,scales = "free",space = "free_y") +
  theme_bw() +
  theme(axis.text = element_text(size = 15,color="black"),
        axis.title = element_text(size = 16),
        strip.text = element_text(size = 17,face="bold"),
        axis.text.y=element_text(size=16,color="black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black")) +
  geom_vline(xintercept = 0,col = "grey",lty = 2) +
  ylab("") +
  scale_y_discrete(labels = function(x) unlist(lapply(str_split(x,"__"), function(y) stringr::str_wrap(y[1], width = 30))))+
  xlab(stringr::str_wrap('OR (95% CI) per unit* increase of exposure', width = 40)) +
  scale_x_continuous(breaks = c(log(0.5),log(0.7),log(1),log(1.5),log(2)),
                     labels = c(0.5,0.7,1,1.5,2))
dev.off()

## for the metabolism-related traits (beta trait) -->BMI
mr_result_filtered_bmi <- mr_result_filtered[mr_result_filtered$traits =="BMI",]
mr_result_filtered_bmi$name_plot <- paste0(mr_result_filtered_bmi$full_name,"__",mr_result_filtered_bmi$traits)
tiff("mr_result_filtered_metabolism_beta_Oct162022.tiff", width = 1300, height = 300, units = "px", res=150)
ggplot(mr_result_filtered_bmi, aes(x = Beta_Odds.ratio, y = reorder(name_plot,-Beta_Odds.ratio))) +
  geom_errorbarh(aes(xmin =CIL, xmax = CIU), color = "#d95f02", height = 0,size=1.2) +
  geom_point(color = "#d95f02",size=2.5) +
  facet_grid(traits~.,scales = "free",space = "free_y") +
  theme_bw() +
  theme(axis.text = element_text(size = 15,color="black"),
        axis.title = element_text(size = 16),
        strip.text = element_text(size = 17,face="bold"),
        axis.text.y=element_text(size=16,color="black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black")) +
  geom_vline(xintercept = 0,col = "grey",lty = 2) +
  ylab("") +
  scale_y_discrete(labels = function(x) unlist(lapply(str_split(x,"__"), function(y) stringr::str_wrap(y[1], width = 30))))+
  xlab(stringr::str_wrap('SD effect size (95% CI) per unit* increase of exposure', width = 40))+
  scale_x_continuous(breaks = c(-0.1,-0.05,0,0.05,0.1,0.2),
                     labels = c(-0.1,-0.05,0,0.05,0.1,0.2))
dev.off()


## for the aging-related traits (OR trait)
mr_result_filtered_aging <- mr_result_filtered[mr_result_filtered$Category=="Aging-related" & mr_result_filtered$traits !="eBMD",]
mr_result_filtered_aging$name_plot <- paste0(mr_result_filtered_aging$full_name,"__",mr_result_filtered_aging$traits)
tiff("mr_result_filtered_aging_OR_Oct162022.tiff", width = 1300, height = 400, units = "px", res=150)
ggplot(mr_result_filtered_aging, aes(x = log(Beta_Odds.ratio), y = reorder(name_plot,-Beta_Odds.ratio))) +
  geom_errorbarh(aes(xmin = log(CIL), xmax = log(CIU)), color = "#7570b3", height = 0,size=1.2) +
  geom_point(color = "#7570b3",size=2.5) +
  facet_grid(traits~.,scales = "free",space = "free_y") +
  theme_bw() +
  theme(axis.text = element_text(size = 15,color="black"),
        axis.title = element_text(size = 16),
        strip.text = element_text(size = 17,face="bold"),
        axis.text.y=element_text(size=16,color="black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black")) +
  geom_vline(xintercept = 0,col = "grey",lty = 2) +
  ylab("") +
  scale_y_discrete(labels = function(x) unlist(lapply(str_split(x,"__"), function(y) stringr::str_wrap(y[1], width = 30))))+
  xlab(stringr::str_wrap('OR (95% CI) per unit* increase of exposure', width = 40)) +
  scale_x_continuous(breaks = c(log(0.6),log(0.7),log(0.8),log(0.9),log(1),log(1.1)),
                     labels = c(0.6,0.7, 0.8,0.9,1,1.1))
dev.off()


## for the aging-related traits (beta trait) -->eBMD
mr_result_filtered_eBMD <- mr_result_filtered[mr_result_filtered$traits =="eBMD",]
#change the (N(1) + N(8))-acetylspermidine --> (N(1)+N(8))-acetylspermidine
mr_result_filtered_eBMD<-mr_result_filtered_eBMD%>%rowwise()%>%mutate(full_name_modi=ifelse(full_name == "N-acetylputrescine / (N(1) + N(8))-acetylspermidine", "N-acetylputrescine / (N(1)+N(8))-acetylspermidine", full_name))

mr_result_filtered_eBMD$name_plot <- paste0(mr_result_filtered_eBMD$full_name_modi,"__",mr_result_filtered_eBMD$traits)
tiff("mr_result_filtered_aging_beta_Oct162022.tiff", width = 1300, height = 1250, units = "px", res=150)
ggplot(mr_result_filtered_eBMD, aes(x = Beta_Odds.ratio, y = reorder(name_plot,-Beta_Odds.ratio))) +
  geom_errorbarh(aes(xmin =CIL, xmax = CIU), color = "#7570b3", height = 0,size=1.2) +
  geom_point(color = "#7570b3",size=2.5) +
  facet_grid(traits~.,scales = "free",space = "free_y") +
  theme_bw() +
  theme(axis.text = element_text(size = 15,color="black"),
        axis.title = element_text(size = 16),
        strip.text = element_text(size = 17,face="bold"),
        axis.text.y=element_text(size=16,color="black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black")) +
  geom_vline(xintercept = 0,col = "grey",lty = 2) +
  ylab("") +
  scale_y_discrete(labels = function(x) unlist(lapply(str_split(x,"__"), function(y) stringr::str_wrap(y[1], width = 30))))+
  xlab(stringr::str_wrap('SD effect size (95% CI) per unit* increase of exposure', width = 40))+
  scale_x_continuous(breaks = c(-0.05,-0.02,0, 0.02,0.05,0.1),
                     labels = c(-0.05,-0.02,0, 0.02, 0.05,0.1))
dev.off()




