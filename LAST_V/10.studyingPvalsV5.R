
#------------------------------------------------------------------------------------
#                           IMPORTING PVALS
#------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------
#               Building Matrix Disease*Pattern filling with NA
#------------------------------------------------------------------------------------

GWAS_filename<-read.delim("OUTPUT/GWAS_results_0001.txt")
colnames(GWAS_filename)
#[1] "disease"           "EnhancerGroup"     "padj"              "AEinPATTERN"       "AEinPATTERNToSNPs"   


summary(GWAS_filename$padj)
#0.000e+00 0.000e+00 3.560e-09 6.891e-06 1.913e-06 9.874e-05 

library(tidyverse)

summary(GWAS_filename$AEinPATTERNToSNPs)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 3.00    8.00   11.00   13.66   16.00  102.00

GWAS.sign<-GWAS_filename

GWAS.sign<-GWAS.sign %>% filter(AEinPATTERNToSNPs>=mean(AEinPATTERNToSNPs))
summary(GWAS.sign$AEinPATTERNToSNPs)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 14.00   16.00   20.00   23.23   25.25  102.00 

GWAS.sign<-GWAS.sign %>% separate(disease, into=c("DISEASE", "PUBMED.ID"), sep="_")

GWAS.collapsed<-GWAS.sign %>% group_by(DISEASE,EnhancerGroup) %>% 
  summarise(sig.PUBMED.IDs=paste0(unique(PUBMED.ID), collapse = ", "),
            N.PUBMED.IDs= length(unique(PUBMED.ID)))
#------------------------------------------------------------------------------------
Patterns <- read.table("OUTPUT/Table_PatternsAE.txt", sep="\t")
Patterns.names <- colnames(Patterns)
Patterns <- Patterns[,c(29,1:27)]
Patterns$matrix <- apply(Patterns[,2:28], 1, paste, collapse = ",")
Patterns <- Patterns[, c("Patron", "matrix")]
colnames(Patterns)[1] <- "EnhancerGroup"

GWAS.collapsed <- left_join(GWAS.collapsed,Patterns)
GWAS.collapsed <- GWAS.collapsed %>% 
  separate(matrix, into=Patterns.names[1:27], sep=",")

GWAS.collapsed <-data.frame(GWAS.collapsed[,1:4], 
           apply(GWAS.collapsed[,5:length(GWAS.collapsed)], 2, as.numeric))


colnames(GWAS.collapsed)
GWAS.collapsed.disease <- 
  GWAS.collapsed %>% group_by(DISEASE) %>%
  dplyr::select(-c("sig.PUBMED.IDs","EnhancerGroup","N.PUBMED.IDs"))%>%
  summarise(across(everything(), ~ sum(., is.na(.), 0)))
#------------------------------------------------------------------------------------
#GWAS.disease.nveces<-data.frame(GWAS.collapsed.disease[,2:28])
GWAS.disease.binary<-ifelse(GWAS.collapsed.disease[,2:28]==0,0,1)
rownames(GWAS.disease.binary)<-GWAS.collapsed.disease$DISEASE

GWAS.disease.binary <-
  GWAS.disease.binary[, colnames(GWAS.disease.binary) %in%  c("endprog", "mesenBM", "endumbprolif") ==
                        FALSE]
GWAS.disease.binary<-GWAS.disease.binary[rowSums(GWAS.disease.binary)>0,]


write.table(GWAS.disease.binary, "GWAS.disease.nveces.txt", sep="\t")
write.table(GWAS.disease.binary, "GWAS.disease.binary.txt", sep="\t")

library(pheatmap)
pheatmap(
  t(GWAS.disease.binary),
  color = c("#cae2eb", "#7da4b3", "#2e6073"),
  border_color = "black",
  scale = "none",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  treeheight_col = 50,
  treeheight_row = 25,
  show_colnames = TRUE,
  fontsize_col = 4,
  fontsize_row = 4,
)

#650 2000
#"#ad2487"

#------------------------------------------------------------------------------------
# Dividing the heathmap by lineage
#------------------------------------------------------------------------------------

GWAS.disease.binary<-read.table("OUTPUT/GWAS.disease.nveces.txt", sep="\t")
       
#"megkaryo""eryth" "meos" "NK" 
Mono.lineage<-c("M2","mono","imDC", "M1", "M0", "mDC", "osteo")
T.lineage<-c("CD8T", "effmemCD8T", "CD4T")
B.lineage<-c("CD38negnaiveB", "germcenterB", "cswitmB", "naiveB", "plasma")
neut.linage<-c("neutmetmyelo", "neutmyelo","mneut","segneutBM", "bfneut")
lymph.linage<-c(T.lineage,B.lineage,"NK" )
  
Mono.lineage.DF<-GWAS.disease.binary[,Mono.lineage]
Mono.lineage.DF<-Mono.lineage.DF[rowSums(Mono.lineage.DF)>0, 
                                 colSums(Mono.lineage.DF)>0 ]

neut.linage.DF<-GWAS.disease.binary[,neut.linage]
neut.linage.DF<-neut.linage.DF[rowSums(neut.linage.DF)>0, 
                               colSums(neut.linage.DF)>0 ]

B.lineage.DF<-GWAS.disease.binary[,B.lineage]
B.lineage.DF<-B.lineage.DF[rowSums(B.lineage.DF)>0, 
                               colSums(B.lineage.DF)>0 ]

lymph.lineage.DF<-GWAS.disease.binary[,lymph.linage]
lymph.lineage.DF<-lymph.lineage.DF[rowSums(lymph.lineage.DF)>0, 
                           colSums(lymph.lineage.DF)>0 ]
# apply(lymph.lineage.DF, 2, as.character)
# sapply(lymph.lineage.DF, as.character)

#col = c("white", "#cae2eb", "#7da4b3", "#2e6073")
# col = c("#f5dcea", "#c284a7", "#cc2b89", "#6a269e")
# col = c("#DCE319FF", "#3CBB75FF", "#287D8EFF", "#440154FF")
col = c("#fff382", "#3CBB75FF", "#287D8EFF", "#440154FF")
#"#cc2b89"
# "#cc2b59",
names(col)=c("0", "1","2","3")
  


pheatmap(
  Mono.lineage.DF,#600, 750
  cellheight = 7,
  cellwidth = 7,
  legend_breaks = c(0,1,2,3),
  color = col,
  border_color = "white",
  # border_color = "black",
  scale = "none",
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  # treeheight_col = 50,
  # treeheight_row = 25,
  show_colnames = TRUE,
  fontsize_col = 7,
  fontsize_row = 7,
)

#annoate number of patherns significant by disease
#------------------------------------------------------------------------------------
# Only monocytes 
#------------------------------------------------------------------------------------

Mono.lineage.DF
colnames(GWAS.disease.binary)
rownames(GWAS.disease.binary)
GWAS.mono.disease<-GWAS.disease.binary[rownames(GWAS.disease.binary) %in% rownames(Mono.lineage.DF),]
GWAS.mono.disease<-GWAS.mono.disease[,colSums(GWAS.mono.disease)>0]
GWAS.mono.disease<-GWAS.mono.disease[,colnames(GWAS.mono.disease) %in% 
                                         colnames(Mono.lineage.DF)==FALSE]

write(rownames(Mono.lineage.DF), "OUTPUT/Mono.lineage.diseases.txt", sep= "\n")
write(colnames(Mono.lineage.DF), "OUTPUT/Mono.lineage.cells.txt", sep= "\n")

##Manually organized
levels.Mono.diseases<-read.delim2("OUTPUT/Mono.lineage.diseases.txt", sep= "\n", header = FALSE)

library(dplyr)

GWAS.mono.disease <-
  GWAS.mono.disease[match(as.character(levels.Mono.diseases$V1),
                          rownames(GWAS.mono.disease)), ]

levels.Mono.cells<-read.delim2("OUTPUT/Mono.lineage.cells.txt", sep= "\n", header = FALSE)

Mono.lineage.DF <-Mono.lineage.DF[,as.character(levels.Mono.cells$V1)]
colnames(Mono.lineage.DF)

Mono.lineage.DF <-Mono.lineage.DF[match(as.character(levels.Mono.diseases$V1),
                                         rownames(Mono.lineage.DF)), ]
colnames(Mono.lineage.DF)

write(colnames(GWAS.mono.disease), "OUTPUT/diseseMonoothercells.txt", sep= "\n")
levels.other.cells<-read.delim2("OUTPUT/diseseMonoothercells.txt", sep= "\n", header = FALSE)

GWAS.mono.disease <-GWAS.mono.disease[,as.character(levels.other.cells$V1)]
colnames(GWAS.mono.disease)

pheatmap(
  cellheight = 7,
  cellwidth = 7,
  GWAS.mono.disease,
  legend_breaks = c(0, 1, 2, 3),
  color = col,
  border_color = "white",
  # border_color = "black",
  scale = "none",
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  # treeheight_col = 50,
  # treeheight_row = 25,
  show_colnames = TRUE,
  fontsize_col = 7,
  fontsize_row = 7,
)


#------------------------------------------------------------------------------------
# Dividing the heathmap by diseease group
#------------------------------------------------------------------------------------

oncolo.traits<-paste0(c("chemotherapy", "cancer", "carcinoma","leukemia", 
                        "lymphoma", "melanoma", "myeloma", "Melanoma"), collapse = "|")
oncolo.traits<-rownames(GWAS.disease.binary)[grepl(oncolo.traits, rownames(GWAS.disease.binary))]

genetic.environm.traits<-paste0(c("ratio", "level", "concentration", "counts", "mean", "Mean", "length", 
                       "trait", "measurement", "biomarker", "tissue", "age", "antibody positivity","density",
                       "triglycerides", "cholesterol", "pressure","index", "Atrial fibrillation", 
                       "Allergic sensitization","responsiveness", "Cholesterol", 
                       "time to first tooth eruption",
                       "interval", "Response to statin therapy", "Smoking behavior",
                       "Smoking cessation in chronic obstructive pulmonary disease",
                       "Temperament", "Triglycerides", "Verbal declarative memory",
                       "White blood cell types", "Word reading" ), collapse = "|")
genetic.environm.traits<-rownames(GWAS.disease.binary)[grepl(genetic.environm.traits, rownames(GWAS.disease.binary))]

clinical.traits<-rownames(GWAS.disease.binary)
clinical.traits <-disease.traits[disease.traits %in% genetic.environm.traits == FALSE &
                   disease.traits %in% oncolo.traits == FALSE]


oncolo.traits.DF<-GWAS.disease.binary[rownames(GWAS.disease.binary)%in%oncolo.traits,]
oncolo.traits.DF<-oncolo.traits.DF[,colSums(oncolo.traits.DF)!=0]

pheatmap(
  t(oncolo.traits.DF),
  color = c("#cae2eb", "#7da4b3", "#2e6073"),
  border_color = "black",
  scale = "none",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  treeheight_col = 50,
  treeheight_row = 25,
  show_colnames = TRUE,
  fontsize_col = 7,
  fontsize_row = 7,
)


clinical.traits.DF<-GWAS.disease.binary[rownames(GWAS.disease.binary)%in%clinical.traits,]
clinical.traits.DF<-clinical.traits.DF[,colSums(clinical.traits.DF)!=0]

pheatmap(
  t(clinical.traits.DF),
  color = c("#cae2eb", "#7da4b3", "#2e6073"),
  border_color = "black",
  scale = "none",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  treeheight_col = 50,
  treeheight_row = 25,
  show_colnames = TRUE,
  fontsize_col = 7,
  fontsize_row = 7,
)



genetic.environm.traits.DF <-
  GWAS.disease.binary[rownames(GWAS.disease.binary) %in% genetic.environm.traits, ]
genetic.environm.traits.DF <-
  genetic.environm.traits.DF[, colSums(genetic.environm.traits.DF) != 0]

pheatmap(
  t(genetic.environm.traits.DF),
  color = c("#cae2eb", "#7da4b3", "#2e6073"),
  border_color = "black",
  scale = "none",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  treeheight_col = 50,
  treeheight_row = 25,
  show_colnames = TRUE,
  fontsize_col = 7,
  fontsize_row = 7,
)
