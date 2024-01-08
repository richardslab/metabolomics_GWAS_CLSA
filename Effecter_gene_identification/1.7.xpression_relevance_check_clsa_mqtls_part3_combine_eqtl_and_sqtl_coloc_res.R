###### the following script merge the colocalization results and filter for gene-metabolite pairs with PP.H4>0.8
library(dplyr)
library(stringr)
library(data.table)
library(tidyr)

matabolites_annotation<-read.csv("./data/clsa_phenotypes/clsa_metabo_annotation_dat")

## combine all coloc results for metabolites/ratios and eqtl
file_list_metabo <- list.files("./restuls/coloc_results/eqtl/")
dataset<-NULL
for (file in file_list_metabo){
  file_path<-paste0("./restuls/coloc_results/eqtl/",file)
  if (!exists("dataset")){
    dataset <- read.csv(file_path, header=TRUE, row.names = 1)
  }
  
  # if the merged dataset does exist, append to it
  if (exists("dataset")){
    temp_dataset <-read.csv(file_path, header=TRUE, row.names = 1)
    dataset<-rbind(dataset, temp_dataset)
    rm(temp_dataset)
  }
}

dataset$full_name<- apply(dataset, 1, function(x) matabolites_annotation$CHEMICAL_NAME[matabolites_annotation$metabo_ID == x[3]])
eqtl_coloc_results_pph4_sig<-dataset%>%filter(PP.H4.abf>0.8)

write.csv(eqtl_coloc_results_pph4_sig,"./restuls/coloc_results/Summary_eqtl_coloc_results.csv")

#########################################################################
###### combine all coloc results for metabolites/ratios and sqtl ########
file_list_metabo <- list.files("./restuls/coloc_results/sqtl")
dataset1<-NULL
for (file in file_list_metabo){
  file_path<-paste0("./restuls/coloc_results/sqtl/",file)
  if (!exists("dataset1")){
    dataset1 <- read.csv(file_path, header=TRUE, row.names = 1)
  }
  
  # if the merged dataset1 does exist, append to it
  if (exists("dataset1")){
    temp_dataset1 <-read.csv(file_path, header=TRUE, row.names = 1)
    dataset1<-rbind(dataset1, temp_dataset1)
    rm(temp_dataset1)
  }
}

dataset1$full_name<- apply(dataset1, 1, function(x) matabolites_annotation$CHEMICAL_NAME[matabolites_annotation$metabo_ID == x[3]])
sqtl_coloc_results_pph4_sig<-dataset1%>%filter(PP.H4.abf>0.8)

write.csv(sqtl_coloc_results_pph4_sig,"./restuls/coloc_results/Summary_sqtl_coloc_results.csv")
