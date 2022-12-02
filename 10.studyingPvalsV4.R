
#------------------------------------------------------------------------------------
#                           IMPORTING PVALS
#------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------
#               Building Matrix Disease*Pattern filling with NA
#------------------------------------------------------------------------------------
EnrichmentEnhparalel_2_100 <- readRDS("/mnt/nocode/juliana/PID/ThemeII/OUTPUT/EnrichmentEnhparalel_2_100.rds")
colnames(EnrichmentEnhparalel_2_100)
# [1] "DISEASE"            "PATTERN"            "AEinPATTERN"        "TotalAEInd"        
# [5] "AEinPATTERNToSNPs"  "TotalAEIndToSNPs"   "AEinPATTERNoutSNPs" "TotalAEIndoutSNPs" 
# [9] "p.value"       

GWAS_results<-EnrichmentEnhparalel_2_100
summary(GWAS_results$p.value)

library(tidyverse)
GWAS.sign<-GWAS_results %>% filter(p.value<0.05)

summary(GWAS.sign$AEinPATTERNToSNPs)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.000   3.000   4.000   6.691   8.000 123.000 

GWAS.sign<-GWAS.sign %>% filter(AEinPATTERNToSNPs>=mean(AEinPATTERNToSNPs))
summary(GWAS.sign$AEinPATTERNToSNPs)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 7.00    8.00   11.00   13.37   15.00  123.00

GWAS.sign<-GWAS.sign %>% separate(DISEASE, into=c("DISEASE", "PUBMED.ID"), sep="_")

GWAS.collapsed<-GWAS.sign %>% group_by(DISEASE,PATTERN) %>% 
  summarise(sig.PUBMED.IDs=paste0(unique(PUBMED.ID), collapse = ", "),
            N.PUBMED.IDs= length(unique(PUBMED.ID)))
#------------------------------------------------------------------------------------
Patterns <- read.table("OUTPUT/Table_PatternsAE.txt", sep="\t")
Patterns.names <- colnames(Patterns)
Patterns <- Patterns[,c(29,1:27)]
Patterns$matrix <- apply(Patterns[,2:28], 1, paste, collapse = ",")
Patterns <- Patterns[, c("Patron", "matrix")]
colnames(Patterns)[1] <- "PATTERN"

GWAS.collapsed <- left_join(GWAS.collapsed,Patterns)
GWAS.collapsed <- GWAS.collapsed %>% 
  separate(matrix, into=Patterns.names[1:27], sep=",")

GWAS.collapsed <-data.frame(GWAS.collapsed[,1:4], 
           apply(GWAS.collapsed[,5:length(GWAS.collapsed)], 2, as.numeric))

GWAS.collapsed.disease <- 
  GWAS.collapsed %>% group_by(DISEASE) %>%
  dplyr::select(-c("sig.PUBMED.IDs","PATTERN","N.PUBMED.IDs"))%>%
  summarise(across(everything(), ~ sum(., is.na(.), 0)))
#------------------------------------------------------------------------------------

GWAS.disease.binary<-ifelse(GWAS.collapsed.disease[,2:28]==0,0,1)
rownames(GWAS.disease.binary)<-GWAS.collapsed.disease$DISEASE

library(pheatmap)
pheatmap(
  t(GWAS.disease.binary),
  color = c("white", "#4682B4"),
  border_color = "black",
  scale = "none",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  treeheight_col = 50,
  treeheight_row = 25,
  show_colnames = TRUE,
  fontsize_col = 4,
  fontsize_row = 3,
)

#650 2000
#"#ad2487"