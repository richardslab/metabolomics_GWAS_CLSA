######## The following script extracts all biological information for metabolites
library(dplyr)
library(xml2)

####### metabolites annotation (containing Metabolon ID, HMDB ID, PubChem ID, and metabolite names)
matabolites_annotation<-read.csv("./data/clsa_phenotypes/clsa_metabo_annotation_dat")

####### load the biological information
### HMDB ###
metabolites_in_test<-read.csv("./results/metabo_with_mqtl.csv",row.names = 1)
metabolites_in_test$HMDB_ID<-apply(metabolites_in_test,1,function(x) matabolites_annotation$HMDB_curated[matabolites_annotation$metabo_ID == x[1]][1])
## remove the metabolites without HMDB ID
metabolites_in_test_filtered<-metabolites_in_test %>% filter(HMDB_ID!="")
## split the metabolites with more than 1 HMDB ID
metabolites_in_test_filtered_2<-metabolites_in_test_filtered %>% separate(HMDB_ID, c("HMDB_1","HMDB_2","HMDB_3","HMDB_4"))
metabolites_in_test_filtered_2_melt<-melt(metabolites_in_test_filtered_2,id="metabolite_list")
metabolites_matching_HMDBID<-metabolites_in_test_filtered_2_melt %>%filter(is.na(value) !=TRUE) %>%select(metabolite_list,value)
colnames(metabolites_matching_HMDBID)<-c("metabolite_ID","HMDB_ID")

#############################################################################
###### Load and process biological information related to metabolites #######
### load gene and unipro information
gene_info <- fread("./data/Homo_sapiens.gene_info.gz", header=TRUE, sep="\t")
uniprot_info <- fread("./data/HUMAN_9606_idmapping.dat__Gene_Name.Gene_Synonym.GeneID.dat", header=FALSE, col.name=c("uniprotId", "type", "value"))
### HMDB ###
serum.xml <- xmlParse(./data/serum_metabolites.xml")
rootnode <- xmlRoot(serum.xml)
xmlIndex <- sapply(1:xmlSize(rootnode), function(i) {xmlValue(rootnode[[i]][["accession"]]) })

xmlIndex_table<-data.frame(
  index=integer(),
  hmdb.ID=character(),
  hmdb.secondary.ID=character())

for(i in 1:xmlSize(rootnode)){
  for (j in 1:xmlSize(rootnode[[i]][["secondary_accessions"]])){
    TMP_index <- data.frame(i,
                      xmlValue(rootnode[[i]][["accession"]]),
                      xmlValue(rootnode[[i]][["secondary_accessions"]][j][["accession"]]))
    names(TMP_index) <- names(xmlIndex_table)
    xmlIndex_table <- rbind(xmlIndex_table, TMP_index)
  }
}
### KEGG ###
kg2keggEntrez <- fread("./data/knownToKeggEntrez.txt.gz", header=FALSE, col.names = c("knownGeneId", "keggHsaId", "keggHsaId.EntrezGeneId"))
keggPathway <- fread("./data/keggPathway.txt.gz", header=FALSE, col.names = c("enst.id", "entrezGeneId", "keggHsaId"))
entrezGene2Hsa <- fread("./data/entrezId-hsa.txt", header=FALSE, col.names = c("entrezGeneId", "keggHsaId"))

##########################################################
###### query biological information for metabolites ######
query.list_clsa<-as.character(unique(metabolites_matching_HMDBID$HMDB_ID))
########################################################################
#### Query HMDB metabolites (associated proteins and genes) ############
DF1 <- data.frame(
  query_HMDB_ID=character(),
  hmdb.accession=character(),
  hmdb.secondary.accession=character(),
  hmdb.status=character(),
  hmdb.name=character(),
  hmdb.protein.name=character(),
  hmdb.uniprot.name=character(),
  hmdb.prtoein.genename=character(),
  Entrez_ID=integer()
)

for (query in query.list_clsa) {
  if (length(unique(xmlIndex_table[which(xmlIndex_table$hmdb.ID == query),1]))!=0) {
    i<-unique(xmlIndex_table[which(xmlIndex_table$hmdb.ID == query),1])
  }
  else {
    i<-unique(xmlIndex_table[which(xmlIndex_table$hmdb.secondary.ID == query),1])
  }
  if (length(i) > 0) {
    # If current metabolite in query list...
    message(xmlValue(rootnode[[i]][["accession"]]))
    # "Iterate over all pathways for the found query metabolite...
    if (xmlSize(rootnode[[i]][["protein_associations"]]) >= 1) {
      for (j in 1:xmlSize(rootnode[[i]][["protein_associations"]])) {
        # Extract proteim name, unipro id, and gene name
        protein_name <-
          xmlValue(rootnode[[i]][["protein_associations"]][[j]][["name"]])
        unipro_id <-
          xmlValue(rootnode[[i]][["protein_associations"]][[j]][["uniprot_id"]])
        gene_name <-
          xmlValue(rootnode[[i]][["protein_associations"]][[j]][["gene_name"]])
        TMP <- data.frame(query,
                          xmlValue(rootnode[[i]][["accession"]]),
                          xmlValue(rootnode[[i]][["secondary_accessions"]]),
                          xmlValue(rootnode[[i]][["status"]]),
                          xmlValue(rootnode[[i]][["name"]]),
                          protein_name,
                          unipro_id,
                          gene_name,
                          ifelse(length(uniprot_info[uniprotId == unipro_id & type == "GeneID", value])>0, uniprot_info[uniprotId == unipro_id & type == "GeneID", value], NA))
        names(TMP) <- names(DF1)
        DF1 <- rbind(DF1, TMP)
      }
    }
  }
  else{ next }
}

DF1_clsa<-left_join(DF1,metabolites_matching_HMDBID, by=c("query_HMDB_ID"="HMDB_ID"))
write.csv(DF1_clsa, "./analysis_dir/biological_relevance_check/clsa_metabolites_HMDB_genes_with_metabo_id.csv")

##############################################################################################################
########## Query KEGG pathway ID, pathway name, and genes that are associated with metabolites ############
DF <- data.frame(
  query_HMDB_ID=character(),
  hmdb.accession=character(),
  hmdb.secondary.accession=character(),
  hmdb.status=character(),
  hmdb.name=character(),
  kegg.name=character(),
  kegg.hsa.id=character(),
  ncbi.gene.id=character(),
  ncbi.gene.symbol=character(),
  ncbi.gene.desc=character(),
  kegg.id=character(),
  kegg.name=character()
)

for (query in query.list_clsa) {
  if (length(unique(xmlIndex_table[which(xmlIndex_table$hmdb.ID == query),1]))!=0) {
    i<-unique(xmlIndex_table[which(xmlIndex_table$hmdb.ID == query),1])
  }
  else {
    i<-unique(xmlIndex_table[which(xmlIndex_table$hmdb.secondary.ID == query),1])
  }
  # If current metabolite in query list...
  if (length(i)>0){
    message(xmlValue(rootnode[[i]][["accession"]]))
    # "Iterate over all pathways for the found query metabolite...
    for (j in 1:xmlSize(rootnode[[i]][["biological_properties"]][["pathways"]])) {
      # Extract KEGG name and KEGG id
      keggid <-
        xmlValue(rootnode[[i]][["biological_properties"]][["pathways"]][[j]][["kegg_map_id"]])
      keggname <-
        xmlValue(rootnode[[i]][["biological_properties"]][["pathways"]][[j]][["name"]])
      # Ensure KEGG id is present (missing when entry is for SMPDB/Pathwhiz only)
      if (!is.na(keggid) & keggid != "") {
        # Convert the cross-species map ("map#####") to the human map ("hsa#####").
        hsaId <- paste0("hsa", substr(keggid, 4, 8))
        # Find all instances of the human map in the UCSC annotation contaning Entrez Gene Ids
        # For each found entry cut out the just the Entrez Gene Id.
        for (geneid in unique(entrezGene2Hsa[keggHsaId == hsaId, entrezGeneId])) {
          TMP1 <- data.frame(query,
                             xmlValue(rootnode[[i]][["accession"]]),
                             xmlValue(rootnode[[i]][["secondary_accessions"]]),
                             xmlValue(rootnode[[i]][["status"]]),
                             xmlValue(rootnode[[i]][["name"]]),
                             keggname,
                             hsaId,
                             geneid,
                             ifelse(length(gene_info[GeneID == as.numeric(geneid), Symbol])>0, gene_info[GeneID == as.numeric(geneid), Symbol], NA),
                             ifelse(length(gene_info[GeneID == as.numeric(geneid), description])>0, gene_info[GeneID == as.numeric(geneid), description], NA),
                             keggid,
                             keggname)
          names(TMP1) <- names(DF)
          DF <- rbind(DF, TMP1)
        }
      }
    }
  }
  else{ next }
}
DF_clsa<-left_join(DF,metabolites_matching_HMDBID, by=c("query_HMDB_ID"="HMDB_ID"))
write.csv(DF_clsa, "./analysis_dir/biological_relevance_check/CLSA_metabolites_KEGG_pathway_genes_with_metabo_id.csv")

############################################################################################
##### extract the Chemical-Gene Co-Occurrences in Literature from PubChem ##################
matabolites_annotation_filtered<-matabolites_annotation %>% filter(PUBCHEM != "" & PUBCHEM != 0)
metabo_id_pubchem_formatted<-matabolites_annotation_filtered %>%select(metabo_ID,PUBCHEM_curated) %>% separate(PUBCHEM_curated,c("PubChem1","PubChem2","PubChem3","PubChem4"))
metabo_id_pubchem_list<-melt(metabo_id_pubchem_formatted, na.rm=T, id.vars="metabo_ID") %>% select(metabo_ID, value) %>%filter(value!="")
Pubchem_list<-unique(metabo_id_pubchem_list$value)

gene_info<-data.frame(PubChem=character(),
                      category_name=character(),
                      gene=character(),
                      evidence=character())
for (metabolite in Pubchem_list){
  PubChem_ID<-as.character(metabolite)
  url_chem<-paste0("https://pubchem.ncbi.nlm.nih.gov/link_db/link_db_server.cgi?format=XML&type=ChemicalGeneSymbolNeighbor&operation=GetAllLinks&id_1=",PubChem_ID)
  chem_lit<- xml2::read_xml(url_chem)
  n_records<-chem_lit%>%xml_children()%>%length()
  message<-chem_lit%>%xml_children()%>%`[`(1)%>%`[`(1)%>%xml_text()
  if (n_records >1 & message !="NotFound"){
    for (n in 1:n_records){
      cat_name<-chem_lit%>%xml_children()%>%`[`(n)%>%xml_children%>%`[`(2)%>%xml_name()
      gene_name<-chem_lit%>%xml_children()%>%`[`(n)%>%xml_children%>%`[`(2)%>%xml_text()
      evidence<-chem_lit%>%xml_children()%>%`[`(n)%>%xml_children%>%`[`(3)%>%xml_text()
      gene_info_intermediate<-data.frame(PubChem_ID,cat_name,gene_name,evidence)
      gene_info<-rbind(gene_info,gene_info_intermediate)
    }}
  else{next}
}

write.csv(gene_info, "metabolite_gene_cooccurancce_PubChem_clsa.csv")
metabolite_gene_co_occurancce_PubChem_clsa<-read.csv("./analysis_dir/biological_relevance_check/metabolite_gene_cooccurancce_PubChem_clsa.csv", row.names = 1)

##get gene names for the enzymes
library(KEGGREST)
metabolite_gene_co_occurancce_PubChem_clsa$gene_name<-as.character(metabolite_gene_co_occurancce_PubChem_clsa$gene_name)
enzyme_list<-metabolite_gene_co_occurancce_PubChem_clsa %>% select(gene_name) %>% filter(str_detect(gene_name, "EC:*"))
enzyme_list_for_check<-unique(unlist(enzyme_list))
enzyme_gene_matching<-data.frame(
  enzyme=character(),
  gene=character()
)
for (enzyme in enzyme_list_for_check) {
  if (length(tryCatch(keggGet(tolower(enzyme)), error=function(e) NULL))==0) {next}
  else{
    enzyme_info<-keggGet(tolower(enzyme))
    enzyme_dat<-data.frame(enzyme_info[[1]]$GENES)
    if(nrow(enzyme_dat)>0){
      colnames(enzyme_dat)<-"genes"
      HSA_genes<-enzyme_dat %>% filter(str_detect(genes, "HSA:"))##extract human genes for the enzyme only
      if (nrow(HSA_genes)>0){
        HSA_gene_list<-strsplit(as.character(HSA_genes$genes), " ")
        for(i in HSA_gene_list){
          gene<-gsub("(?<=\\()[^()]*(?=\\))(*SKIP)(*F)|.", "", i, perl=T)
          print(enzyme)
          print(gene)
          gene_interm<-data.frame(enzyme, gene)
          enzyme_gene_matching<-rbind(enzyme_gene_matching, gene_interm)
        }
      }
      else{
        gene_interm<-data.frame(enzyme, NA)
        colnames(gene_interm)<-c("enzyme", "gene")
        enzyme_gene_matching<-rbind(enzyme_gene_matching, gene_interm)
      }
    }
    else{next}
    }
}
enzyme_gene_matching_filtered1<-enzyme_gene_matching %>% filter(!is.na(gene) & gene != "")
write.csv(enzyme_gene_matching_filtered1, "./analysis_dir/biological_relevance_check/enzyme_gene_matching_filtered_clsa.csv")

enzyme_gene_matching_filtered<-read.csv("./analysis_dir/biological_relevance_check/enzyme_gene_matching_filtered_clsa.csv", header=T, row.names = 1)
enzyme_gene_matching_filtered$enzyme<-as.character(enzyme_gene_matching_filtered$enzyme)
enzyme_gene_matching_filtered$gene<-as.character(enzyme_gene_matching_filtered$gene)
metabolite_gene_co_occurancce_PubChem_clsa_1<-left_join(metabolite_gene_co_occurancce_PubChem_clsa, enzyme_gene_matching_filtered,by=c("gene_name"="enzyme"))

metabolite_gene_co_occurancce_PubChem_clsa_1$Gene_Symbol<-apply(metabolite_gene_co_occurancce_PubChem_clsa_1, 1, 
                                                                function(x) ifelse(is.na(x[5]), toupper(x[3]), x[5]))
gene_info <- fread("Homo_sapiens.gene_info.gz", header=TRUE, sep="\t")
metabolite_gene_co_occurancce_PubChem_clsa_1$entrez_ID<-apply(metabolite_gene_co_occurancce_PubChem_clsa_1,1, function(x) gene_info[gene_info$Symbol == toupper(x[6]),]$GeneID[1])
write.csv(metabolite_gene_co_occurancce_PubChem_clsa_1,"./analysis_dir/biological_relevance_check/pubchem_metabo_gene_cooccurance_database_clsa.csv")
