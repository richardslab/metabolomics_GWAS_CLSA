library(circlize)
library(tidyr)
library(dplyr)
library(data.table)

setwd("./analysis_dir/plotting")
clsa_metabolite_ratios_cojo_SNPs_formatted_rsid_bonf_filtered_enzyme_transporter<-read.csv("clsa_metabolite_ratios_cojo_SNPs_formatted_rsid_bonf_filtered_enzyme_transporter_with_locus",row.names = 1)

#add super pathway info for metabo1 and metabo2
CLSA_metabo_annotation<-read.csv("CLSA_COMBINED_ANNOTATIONTABLEALL_v1_modi_July2021.csv")
clsa_metabolite_ratios_cojo_SNPs_formatted_rsid_bonf_filtered_enzyme_transporter$SUPER_PATHWAY_mod1<-apply(clsa_metabolite_ratios_cojo_SNPs_formatted_rsid_bonf_filtered_enzyme_transporter,1,function(x) CLSA_metabo_annotation$SUPER_PATHWAY[CLSA_metabo_annotation$metabo_ID == x[21]])
clsa_metabolite_ratios_cojo_SNPs_formatted_rsid_bonf_filtered_enzyme_transporter$SUPER_PATHWAY_mod2<-apply(clsa_metabolite_ratios_cojo_SNPs_formatted_rsid_bonf_filtered_enzyme_transporter,1,function(x) CLSA_metabo_annotation$SUPER_PATHWAY[CLSA_metabo_annotation$metabo_ID == x[22]])

clsa_metabolite_ratio_cojo_SNPs_filtered1_sub<-clsa_metabolite_ratios_cojo_SNPs_formatted_rsid_bonf_filtered_enzyme_transporter%>% select(metabo1, metabo2, SUPER_PATHWAY_mod1,SUPER_PATHWAY_mod2, metabolite, P, metabo1_name, metabo2_name)
colnames(clsa_metabolite_ratio_cojo_SNPs_filtered1_sub)<-c("metabo1_name", "metabo2_name", "SUPER_PATHWAY_mod1","SUPER_PATHWAY_mod2", "metabolite", "P", "full_name_mod1", "full_name_mod2")
## keep the smallest P for each metabolite pairs
clsa_metabolite_ratio_cojo_SNPs_filtered1_sub1<-clsa_metabolite_ratio_cojo_SNPs_filtered1_sub %>% rowwise() %>% mutate(fix.p=ifelse(P <= 0, 5e-324, P)) %>% 
  mutate(neg_log10_P = -log10(fix.p))%>%group_by(metabolite)%>% slice(which.max(neg_log10_P))

plotting_set2<-clsa_metabolite_ratio_cojo_SNPs_filtered1_sub1%>%ungroup()%>%select(full_name_mod1,full_name_mod2,neg_log10_P)
## convert long table to wide table

plotting_set_wide <- spread(plotting_set2, full_name_mod1, neg_log10_P, fill=0)
plotting_set_wide1<-data.frame(plotting_set_wide, check.names = F)
rownames(plotting_set_wide1)<-plotting_set_wide1$full_name_mod2
plotting_set_wide2<-plotting_set_wide1[,2:59]

## create the group match table
metabo1_set<-clsa_metabolite_ratios_cojo_SNPs_formatted_rsid_bonf_filtered_enzyme_transporter %>%select(metabo1_name, SUPER_PATHWAY_mod1) %>% distinct()
colnames(metabo1_set)<-c("full_name_mod","super_pathway")
metabo2_set<-clsa_metabolite_ratios_cojo_SNPs_formatted_rsid_bonf_filtered_enzyme_transporter %>%select(metabo2_name, SUPER_PATHWAY_mod2) %>% distinct()
colnames(metabo2_set)<-c("full_name_mod","super_pathway")


metabo_set<-rbind(metabo1_set,metabo2_set) %>%distinct()


group = structure(as.vector(metabo_set$super_pathway), names = as.vector(metabo_set$full_name_mod))

# set grid.col based on the pvalue
grid.col = structure(c(rep("#fc8d62", length(metabo_set[metabo_set$super_pathway=="Amino Acid",]$full_name_mod)),
                       rep("#8da0cb", length(metabo_set[metabo_set$super_pathway=="Carbohydrate",]$full_name_mod)),
                       rep("#e78ac3", length(metabo_set[metabo_set$super_pathway=="Cofactors and Vitamins",]$full_name_mod)),
                       rep("#a6d854", length(metabo_set[metabo_set$super_pathway=="Energy",]$full_name_mod)),
                       rep("#66c2a5", length(metabo_set[metabo_set$super_pathway=="Lipid",]$full_name_mod)),
                       rep("#ffd92f", length(metabo_set[metabo_set$super_pathway=="Nucleotide",]$full_name_mod)),
                       rep("#e5c494", length(metabo_set[metabo_set$super_pathway=="Peptide",]$full_name_mod)),
                       rep("#1f78b4", length(metabo_set[metabo_set$super_pathway=="Xenobiotics",]$full_name_mod))),
                     names = c(as.vector(metabo_set[metabo_set$super_pathway=="Amino Acid",]$full_name_mod),
                               as.vector(metabo_set[metabo_set$super_pathway=="Carbohydrate",]$full_name_mod),
                               as.vector(metabo_set[metabo_set$super_pathway=="Cofactors and Vitamins",]$full_name_mod),
                               as.vector(metabo_set[metabo_set$super_pathway=="Energy",]$full_name_mod),
                               as.vector(metabo_set[metabo_set$super_pathway=="Lipid",]$full_name_mod),
                               as.vector(metabo_set[metabo_set$super_pathway=="Nucleotide",]$full_name_mod),
                               as.vector(metabo_set[metabo_set$super_pathway=="Peptide",]$full_name_mod),
                               as.vector(metabo_set[metabo_set$super_pathway=="Xenobiotics",]$full_name_mod))
)

pathway_color_match<-data.frame(c("#fc8d62","#8da0cb","#e78ac3","#a6d854","#66c2a5","#ffd92f","#e5c494","#1f78b4"),
                                c("Amino Acid","Carbohydrate","Cofactors and Vitamins","Energy","Lipid","Nucleotide","Peptide","Xenobiotics"),
                                c(1,2,3,4,5,6,7,8))
colnames(pathway_color_match)<-c("color","super_pathway","index")

#################################
plotting_set2$width<-rep(2,nrow(plotting_set2))
plotting_set3<-plotting_set2 %>% rowwise() %>% mutate(grey_scale = ifelse(neg_log10_P>=100, 20,
                                                                          ifelse(neg_log10_P>=50&neg_log10_P<100, 15,10)))
#include a column for line
plotting_set4<-left_join(plotting_set3,metabo_set,by=c("full_name_mod1"="full_name_mod"))
plotting_set5<-left_join(plotting_set4,pathway_color_match,by="super_pathway")

col_fun = colorRamp2(c(10,15,20), c("#bdbdbd","#636363","#252525"), transparency = 0.1)
link_fun=colorRamp2(c(1,2,3,4,5,6,7,8),
                    c("#fc8d62","#8da0cb","#e78ac3","#a6d854","#66c2a5","#ffd92f","#e5c494","#1f78b4"))

## remove 4 metabolite names with shorter form
plotting_set6<-plotting_set5%>%rowwise()%>%mutate(full_name_mod1_new=ifelse(full_name_mod1 == "N-acetylglucosamine/N-acetylgalactosamine", "GlcNAc/alpha-GalNAc", 
                                                                        ifelse(full_name_mod1 =="linoleoyl-arachidonoyl-glycerol (18:2/20:4) [1]*","diacylglycerol 1",
                                                                               ifelse(full_name_mod1 =="linoleoyl-arachidonoyl-glycerol (18:2/20:4) [2]*","diacylglycerol 2",
                                                                                      ifelse(full_name_mod1 =="5-acetylamino-6-formylamino-3-methyluracil","AFMU", 
                                                                                             ifelse(full_name_mod1 == "oleoyl-linoleoyl-glycerol (18:1/18:2) [2]", "diacylglycerol 3",
                                                                                                    ifelse(full_name_mod1 =="flavin adenine dinucleotide (FAD)", "flavin adenine dinucleotide", 
                                                                                                           ifelse(full_name_mod1 =="(N(1) + N(8))-acetylspermidine", "(N(1)+N(8))-acetylspermidine",
                                                                                                                  ifelse(full_name_mod1 =="S-adenosylhomocysteine (SAH)", "S-adenosylhomocysteine",
                                                                                                                         ifelse(full_name_mod1 =="5-methylthioadenosine (MTA)", "5-methylthioadenosine",as.character(full_name_mod1)))))))))),
                                                  full_name_mod2_new=ifelse(full_name_mod2 == "N-acetylglucosamine/N-acetylgalactosamine", "GlcNAc/alpha-GalNAc", 
                                                                            ifelse(full_name_mod2 =="linoleoyl-arachidonoyl-glycerol (18:2/20:4) [1]*","diacylglycerol 1",
                                                                                   ifelse(full_name_mod2 =="linoleoyl-arachidonoyl-glycerol (18:2/20:4) [2]*","diacylglycerol 2",
                                                                                          ifelse(full_name_mod2 =="5-acetylamino-6-formylamino-3-methyluracil","AFMU",
                                                                                                 ifelse(full_name_mod2 == "oleoyl-linoleoyl-glycerol (18:1/18:2) [2]", "diacylglycerol 3",
                                                                                                        ifelse(full_name_mod2 =="flavin adenine dinucleotide (FAD)", "flavin adenine dinucleotide", 
                                                                                                               ifelse(full_name_mod2 =="(N(1) + N(8))-acetylspermidine", "(N(1)+N(8))-acetylspermidine",
                                                                                                                      ifelse(full_name_mod2 =="S-adenosylhomocysteine (SAH)", "S-adenosylhomocysteine",
                                                                                                                             ifelse(full_name_mod2 =="5-methylthioadenosine (MTA)", "5-methylthioadenosine",as.character(full_name_mod2)))))))))))


metabo_set_new<-metabo_set%>%rowwise()%>%mutate(full_name_mod = ifelse(full_name_mod == "N-acetylglucosamine/N-acetylgalactosamine", "GlcNAc/alpha-GalNAc", 
                                                                       ifelse(full_name_mod =="linoleoyl-arachidonoyl-glycerol (18:2/20:4) [1]*","diacylglycerol 1",
                                                                              ifelse(full_name_mod =="linoleoyl-arachidonoyl-glycerol (18:2/20:4) [2]*","diacylglycerol 2",
                                                                                     ifelse(full_name_mod =="5-acetylamino-6-formylamino-3-methyluracil","AFMU", 
                                                                                            ifelse(full_name_mod =="oleoyl-linoleoyl-glycerol (18:1/18:2) [2]", "diacylglycerol 3", 
                                                                                                   ifelse(full_name_mod =="flavin adenine dinucleotide (FAD)", "flavin adenine dinucleotide", 
                                                                                                          ifelse(full_name_mod =="(N(1) + N(8))-acetylspermidine", "(N(1)+N(8))-acetylspermidine",
                                                                                                                 ifelse(full_name_mod =="S-adenosylhomocysteine (SAH)", "S-adenosylhomocysteine",
                                                                                                                        ifelse(full_name_mod =="5-methylthioadenosine (MTA)", "5-methylthioadenosine",as.character(full_name_mod)))))))))))

group1 = structure(as.vector(metabo_set_new$super_pathway), names = as.vector(metabo_set_new$full_name_mod))

plotting_set6$full_name_mod1_new<-factor(plotting_set6$full_name_mod1_new)
plotting_set6$full_name_mod2_new<-factor(plotting_set6$full_name_mod2_new)

grid.col1 = structure(c(rep("#fc8d62", length(metabo_set_new[metabo_set_new$super_pathway=="Amino Acid",]$full_name_mod)),
                       rep("#8da0cb", length(metabo_set_new[metabo_set_new$super_pathway=="Carbohydrate",]$full_name_mod)),
                       rep("#e78ac3", length(metabo_set_new[metabo_set_new$super_pathway=="Cofactors and Vitamins",]$full_name_mod)),
                       rep("#a6d854", length(metabo_set_new[metabo_set_new$super_pathway=="Energy",]$full_name_mod)),
                       rep("#66c2a5", length(metabo_set_new[metabo_set_new$super_pathway=="Lipid",]$full_name_mod)),
                       rep("#ffd92f", length(metabo_set_new[metabo_set_new$super_pathway=="Nucleotide",]$full_name_mod)),
                       rep("#e5c494", length(metabo_set_new[metabo_set_new$super_pathway=="Peptide",]$full_name_mod)),
                       rep("#1f78b4", length(metabo_set_new[metabo_set_new$super_pathway=="Xenobiotics",]$full_name_mod))),
                     names = c(as.vector(metabo_set_new[metabo_set_new$super_pathway=="Amino Acid",]$full_name_mod),
                               as.vector(metabo_set_new[metabo_set_new$super_pathway=="Carbohydrate",]$full_name_mod),
                               as.vector(metabo_set_new[metabo_set_new$super_pathway=="Cofactors and Vitamins",]$full_name_mod),
                               as.vector(metabo_set_new[metabo_set_new$super_pathway=="Energy",]$full_name_mod),
                               as.vector(metabo_set_new[metabo_set_new$super_pathway=="Lipid",]$full_name_mod),
                               as.vector(metabo_set_new[metabo_set_new$super_pathway=="Nucleotide",]$full_name_mod),
                               as.vector(metabo_set_new[metabo_set_new$super_pathway=="Peptide",]$full_name_mod),
                               as.vector(metabo_set_new[metabo_set_new$super_pathway=="Xenobiotics",]$full_name_mod))
)

circos.clear()
tiff("Fig3_chord_diagram.tiff", units="px", width=2500, height=2700, res=300)
par(cex = 0.7, mar = c(0, 0, 0, 0))
chordDiagram(plotting_set6[,c(9,10,4)], group = group1, grid.col=grid.col1, annotationTrack = "grid",col=col_fun(plotting_set6[,5]),
             preAllocateTracks = list(track.height = 0.3), link.border=link_fun(plotting_set6$index), link.lwd = 1)#link.border = "white"

circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
}, bg.border = NA) # here set bg.border to NA is important
dev.off()


