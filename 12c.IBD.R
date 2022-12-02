
#------------------------------------------------------------------------------------
#                           IMPORTING PVALS
#------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------
#               Building Matrix Disease*Pattern filling with NA
#------------------------------------------------------------------------------------

GWAS_filename<-read.delim("OUTPUT/GWAS_results_0001.txt")
colnames(GWAS_filename)
#[1] "disease"           "EnhancerGroup"     "padj"              "AEinPATTERN"       "AEinPATTERNToSNPs"   

library(tidyverse)

GWAS.sign<-GWAS_filename

GWAS.sign<-GWAS.sign %>% filter(AEinPATTERNToSNPs>=mean(AEinPATTERNToSNPs))
summary(GWAS.sign$AEinPATTERNToSNPs)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 14.00   16.00   20.00   23.23   25.25  102.00 

GWAS.sign<-GWAS.sign %>% separate(disease, into=c("DISEASE", "PUBMED.ID"), sep="_")

GWAS.collapsed<-GWAS.sign %>% group_by(DISEASE,EnhancerGroup) %>% 
  summarise(sig.PUBMED.IDs=paste0(unique(PUBMED.ID), collapse = ", "),
            Nsig.PUBMED.IDs= length(unique(PUBMED.ID)))



# IBD.traits <- c(
# "Ulcerative colitis",
# "Crohn's disease",
# "Pyoderma gangrenosum in inflammatory bowel disease"
# )

IBD.traits <- c(
  "Ulcerative colitis")


IBD.enrich <- GWAS.collapsed[GWAS.collapsed$DISEASE %in% IBD.traits,]
#------------------------------------------------------------------------------------
Patterns <- read.table("OUTPUT/Table_PatternsAE.txt", sep="\t")
Patterns.names <- colnames(Patterns)
Patterns <- Patterns[,c(29,1:27)]
Patterns$matrix <- apply(Patterns[,2:28], 1, paste, collapse = ",")
Patterns <- Patterns[, c("Patron", "matrix")]
colnames(Patterns)[1] <- "EnhancerGroup"

IBD.enrich <- left_join(IBD.enrich,Patterns)

annot<-IBD.enrich[, c("DISEASE", "EnhancerGroup", "matrix")]
annot$values<-1

annot <-pivot_wider(annot, names_from = DISEASE, values_from = values)


annotlefth<-annot[, 1:2]
annotlefth <-annotlefth %>% 
  separate(matrix, into=Patterns.names[1:27], sep=",")

annotlefth <-data.frame(annotlefth[,1], 
           apply(annotlefth[,2:length(annotlefth)], 2, as.numeric))

rownames(annotlefth)<-annotlefth$EnhancerGroup

numericAnnotLeft<-annotlefth[,2:length(annotlefth)]

numericAnnotLeft<-numericAnnotLeft[,colSums(numericAnnotLeft)>0]

numericAnnotLeft <- numericAnnotLeft[c(
  "effmemCD8T",
  "germcenterB|naiveB",
  "naiveB",
  "mDC",
  "imDC|mDC",
  "mneut|segneutBM",
  "eryth",
  "endprog|endumbprolif",
  "endprog|endumbprolif|mesenBM"
), c(
  "effmemCD8T",
  "germcenterB",
  "naiveB",
  "mDC",
  "imDC",
  "segneutBM",
  "mneut",
  "eryth",
  "mesenBM",
  "endumbprolif",
  "endprog"
)]


library(pheatmap)
pheatmap(
  numericAnnotLeft,
  color = c("#642194", "yellow"),
  border_color = "black",
  scale = "none",
  cellwidth = 10 ,
  cellheight = 10,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  treeheight_col = 50,
  treeheight_row = 25,
  show_colnames = TRUE,
  fontsize_col = 10,
  fontsize_row = 10,
)


annotrigth<-annot[, c(1, 3:5)]
colnames(annotrigth)<- c("Pattern", "crohn disease", 
                         "pyoderma gangrenosum", "ulcerative colitis")
annotrigth <- data.frame(annotrigth[, 1],
                         apply(annotrigth[, 2:length(annotrigth)], 2, as.numeric))
rownames(annotrigth)<-annotrigth$Pattern
annotrigth[is.na(annotrigth)]<-0



annotrigth <- annotrigth[c(
  "effmemCD8T",
  "germcenterB|naiveB",
  "naiveB",
  "mDC",
  "imDC|mDC",
  "mneut|segneutBM",
  "eryth",
  "endprog|endumbprolif",
  "endprog|endumbprolif|mesenBM"
), c(2, 4, 3)]

pheatmap(
  annotrigth,
  color = c("#642194", "yellow"),
  cellwidth = 10 ,
  cellheight = 10,
  border_color = "black",
  scale = "none",
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  treeheight_col = 50,
  treeheight_row = 25,
  show_colnames = TRUE,
  fontsize_col = 10,
  fontsize_row = 10,
)
