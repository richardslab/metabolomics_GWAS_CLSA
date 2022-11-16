library(data.table)
library(ggplot2)
library(dplyr)
library(ggpubr)

setwd("./analysis_dir/plotting")
clsa_metabolite_cojo_SNPs_filtered1 <- read.csv("clsa_metabolite_cojo_SNPs_with_loci", row.names = 1)

#### plot n metabolite per locus (pleitropy)
n_metaoblite_per_locus<-clsa_metabolite_cojo_SNPs_filtered1 %>% select(locus,metabolite)%>% group_by(locus)%>%summarise(n_metabolite=length(unique(metabolite)))
n_metaoblite_per_locus<-data.frame(n_metaoblite_per_locus)
n_metaoblite_per_locus1 <- n_metaoblite_per_locus[order(-n_metaoblite_per_locus$n_metabolite),]
median(n_metaoblite_per_locus1$n_metabolite)#2
n_metaoblite_per_locus1$lous_index<-seq(1,length(n_metaoblite_per_locus1$locus), 1)

p_pleio_line<-ggplot(data=n_metaoblite_per_locus2, aes(x=lous_index/248*100, y=log2(n_metabolite+1))) + geom_line(color="blue", size=2) +
  scale_y_continuous(breaks = unique(sort(log2(n_metaoblite_per_locus2$n_metabolite+1)))[seq(1,length(unique(n_metaoblite_per_locus2$n_metabolite)),3)], 
                     labels =unique(sort(n_metaoblite_per_locus1$n_metabolite))[seq(1,length(unique(n_metaoblite_per_locus1$n_metabolite)),3)])+
  xlab(seq(1,100,10))+
  theme_bw() +
  theme(axis.text.x = element_text(size=25, angle = 20, hjust=1, colour="black"), axis.text.y = element_text(size=20, colour="black"),
        axis.title=element_text(size=30),  legend.text=element_text(size=20), legend.title = element_text(size=25),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black"))+
  labs(x="Percentage of all loci", y="Number of associated metabolites", size=100)

tiff("Fig2D.tiff", width = 800, height = 750, units = "px")
p_pleio_line
dev.off()
