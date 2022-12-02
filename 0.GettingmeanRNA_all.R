
setwd("/mnt/nocode/juliana/PID/ThemeII")

#Import RNAseq MAtrix

dfcountsEC <-
  readRDS(
    "/mnt/nocode/juliana/PID/Chromatin_States_Unpublished/results/2022/RNAseq_matrix/dfcountsEC.rds"
  )

#Format the matrix to add celltype metadata column

dfcountsEC$GENE_ID<-rownames(dfcountsEC)
#----------------------------------------------------------------------------

dfcountsEC_list<-dfcountsEC$Gene_ID
#Two genes are not in the list of the RNAseq
# ##Still dont know wich ones
# #Two genes are not in the list of the RNAseq
# unique(dfcountsEC_list[dfcountsEC_list%in%gsub("\\..*", "",dfcountsEC_RNA_EC$GENE_ID)==FALSE])
# # ENSG00000289620 ENSG00000289007


#----------------------------------------------------------------------------
library(tidyverse)

dfcountsECRNAEC_long <-
  pivot_longer(
    data = dfcountsEC,
    cols = 1:length(dfcountsEC) - 1,
    names_to = "SAMPLE_NAME",
    values_to = "EC"
  )

#Add metadata column with celltype info

#Import metadata
all_bluprint_metadata <-
  read.delim(
    "~/../../mnt/nocode/juliana/PID/Chromatin_States_Blueprint/20160816.data",
    sep = "\t",
    header = TRUE
  )

metadata<-
  all_bluprint_metadata %>% 
  select("BIOMATERIAL_TYPE",
  "SAMPLE_NAME","CELL_TYPE",
  "TISSUE_TYPE","DONOR_SEX",
  "DONOR_ID") %>% 
  unique

dfcountsECRNAEC_long<-left_join(dfcountsECRNAEC_long,metadata)

write.table(dfcountsECRNAEC_long, "OUTPUT/RNAEC_long.txt", sep="\t")


Counts_sex <- dfcountsECRNAEC_long %>%
  select(
    "BIOMATERIAL_TYPE",
    "SAMPLE_NAME",
    "CELL_TYPE",
    "TISSUE_TYPE",
    "DONOR_SEX",
    "DONOR_ID"
  ) %>%
  unique %>%
  count(CELL_TYPE, DONOR_SEX)

Counts_sex <-pivot_wider(Counts_sex,names_from = "DONOR_SEX",values_from = "n")
write.table(Counts_sex, "OUTPUT/Counts_sex_RNAseq.txt", sep="\t")

#Get the order form here
Genes<-readRDS("extData/gene_Anotations38_lcincluded.rds") %>% data.frame()
Genes$Gene_names <- ifelse(
  Genes$gene_name == "",
  as.character(Genes$ID),
  as.character(Genes$gene_name)
)


colnames(Genes)[4]<-"GENE_ID"
Genes$GENE_ID<-gsub("\\..*", "", Genes$GENE_ID)
dfcountsECRNAEC_long$GENE_ID<-gsub("\\..*", "", dfcountsECRNAEC_long$GENE_ID)

Genes_annot<-Genes %>% data.frame()%>%
  select(GENE_ID, Gene_names, gene_type)%>%
  unique()

dfcountsECRNAEC_long<-left_join(dfcountsECRNAEC_long,Genes_annot)

write.table(dfcountsECRNAEC_long, "OUTPUT/RNAEC_long.txt", sep="\t")

ECAvg <-
  dfcountsECRNAEC_long %>% group_by(Gene_names,CELL_TYPE) %>% summarise(
    VALS = paste0(EC, collapse = ", "),
    MEAN = mean(EC),
    SD = sd(EC),
    GENE_ID=paste(unique(GENE_ID), collapse=", "),
    gene_type= unique(gene_type))

saveRDS(ECAvg, "OUTPUT/ECAvg.rds")
write.table(ECAvg, "OUTPUT/ECAvg.txt", sep="\t")

