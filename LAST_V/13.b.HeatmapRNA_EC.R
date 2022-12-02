#-----------------------------------------------------------------------
#Heatmap  with  genes expression
#-----------------------------------------------------------------------
#Import  genes average expression
#---------------------------------------------------------------------
setwd("/mnt/nocode/juliana/")
#---------------------------------------------------------------------
dfcountsECRNAEC_long<-read.table("OUTPUT/RNAEC_long.txt", sep="\t")

EnrichedGenes<-read.table("OUTPUT/IBD.GENES.DF.txt", sep="\t")
EnrichedGenes<-read.table("OUTPUT/Lipid.GENESandlc.DF.txt", sep="\t")
EnrichedGenesCounts<-dfcountsECRNAEC_long %>% filter( dfcountsECRNAEC_long$Gene_names %in% EnrichedGenes$gene_name)

ECAvg <-
  EnrichedGenesCounts %>% group_by(Gene_names,CELL_TYPE) %>% summarise(
    VALS = paste0(EC, collapse = ", "),
    MEAN = mean(EC),
    SD = sd(EC),
    GENE_ID=paste(unique(GENE_ID), collapse=", "),
    gene_type= unique(gene_type))

EC.Matrix<-pivot_wider(ECAvg, id_cols = Gene_names, names_from = CELL_TYPE,
              values_from = MEAN) %>% data.frame()

rownames(EC.Matrix)<-EC.Matrix$Gene_names

annotationrow<-ECAvg %>% select(Gene_names, gene_type)%>%unique()%>%data.frame()
row.names(annotationrow)<-annotationrow$Gene_names
annotationrow<-annotationrow%>% select(gene_type)

EC.Matrix <- 
  EC.Matrix %>% 
  mutate(Gene_names = factor(Gene_names,levels = rownames(annotationrow))) %>% 
  arrange(Gene_names)

EC.Matrix<-EC.Matrix %>% select(-c("Gene_names", "NA."))

#-----------------------------------------------------------------------
#Plot heatmap
#---------------------------------------------------------------------
celltype_abb<-read.table("OUTPUT/celltype_Abb.txt", sep="\t")
colannot<-data.frame(variable=colnames(EC.Matrix))
colannot<-left_join(colannot,celltype_abb)
colnames(EC.Matrix)<-colannot$acronym

celltypes<-read.delim("INPUTS/celltype_data2.txt")
EC.Matrix<-EC.Matrix%>% select(celltypes$PATTERN_NAME)


library(pheatmap)

pheatmap(EC.Matrix , scale= "row",
         annotation_row = annotationrow,cluster_rows =TRUE,
         cluster_cols = FALSE,
         fontsize = 7, fontsize_row = 7,fontsize_col = 7, border_color = "white",
         cellwidth = 7,cellheight = 7)





