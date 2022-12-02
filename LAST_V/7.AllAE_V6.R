#-----------------------------------------------------------------
#               IMPORTING M1 and M2 ENHANCERS 
#   (conserved in the 3 replicates of each cell type)
#------------------------------------------------------------------
#Which enhacers are GWAS enriched
setwd("/mnt/nocode/juliana/PID/ThemeII")

# #Chromstates conserved in 75% of the replicates
conserved<-readRDS("OUTPUT/ConservedDF.rds")

# IDS of AE in at least 2 samples
AllAE<-readRDS("OUTPUT/BActiveE.rds")

##Only AE
All_dataAE<-conserved[rownames(conserved)%in%AllAE,]
#Maybe some enhancers are not conserved in the celltype?
#27 celltypes

##Rename if AE or no
data_modif<-as.matrix(All_dataAE)
data_modif[data_modif %in% c(1,2,3,7,5)] <- 0#"No_E", 5 "PoisE"
data_modif[data_modif %in% c(4,6)] <- 1#"AE
colnames(data_modif)<-gsub("\\.", "", colnames(data_modif))
#--------------------------------------------------------------------------------------------------------------
library(plyr)
Patterns<-plyr::count(df = data.frame(data_modif), vars = colnames(data_modif))
Patterns<-Patterns[order(Patterns$freq, decreasing = TRUE),]
length(Patterns$freq)#107562

colnames(Patterns)<-c("endprog", "M2","bfneut","mono","megkaryo", "CD38negnaiveB","CD4T","CD8T","cswitmB","NK",
  "effmemCD8T", "endumbprolif", "eryth","germcenterB","imDC", "M1", "M0","mDC","meos", 
  "mneut","mesenBM", "naiveB","neutmetmyelo","neutmyelo","osteo","plasma","segneutBM","freq")

patron<-function(x){
  x<-which(x==1)
  return(x)
}

Patrons<-list()
for (i in 1:length(Patterns$freq)){
  patrons<-apply(Patterns[i,], 1, patron)
  patrons<-paste(colnames(Patterns[i,])[patrons], collapse = "|")
  Patrons[[i]]<-patrons
  print(patrons)
}

Patterns$Patron<-do.call(rbind,Patrons)
write.table(Patterns, "OUTPUT/Table_PatternsAE.txt", sep="\t")
#--------------------------------------------------------------------------
Mylist<-data.frame(data_modif)
Milista<-as.matrix(Mylist)
colnames(Milista)<-c("endprog", "M2","bfneut","mono","megkaryo", "CD38negnaiveB","CD4T","CD8T","cswitmB","NK",
                      "effmemCD8T", "endumbprolif", "eryth","germcenterB","imDC", "M1", "M0","mDC","meos", 
                      "mneut","mesenBM", "naiveB","neutmetmyelo","neutmyelo","osteo","plasma","segneutBM")

patron2 <- function(y) {
  patrons <- patron(y)
  patrons <- colnames(Milista)[patrons]
  patrons <- paste(patrons, collapse = "|")
  return(patrons)
}

Patrons<-apply(Milista, 1, patron2)
Mylist2<-cbind(Milista,Patrons)
Mylist <-data.frame(Mylist2)
head(Mylist)

CtypeChrStates <- Mylist
CtypeChrStates$regionID <- rownames(CtypeChrStates)
CtypeChrStates$coord <- rownames(CtypeChrStates)

library(GenomicRanges); library(tidyr)

CtypeChrStates <-
  tidyr::separate(CtypeChrStates, coord, c("chr", "region"))
CtypeChrStates$end <- as.numeric(CtypeChrStates$region) * 200
CtypeChrStates$start <- CtypeChrStates$end - 200


gr_CtypeChrStates<- makeGRangesFromDataFrame(
  CtypeChrStates,
  keep.extra.columns = TRUE,
  ignore.strand = TRUE,
  seqinfo = NULL,
  seqnames.field = c(
    "seqnames",
    "seqname",
    "chromosome",
    "chrom",
    "chr",
    "chromosome_name",
    "seqid"
  ),
  start.field = "start",
  end.field = "end",
  strand.field = NULL,
  starts.in.df.are.0based = TRUE
)

saveRDS(gr_CtypeChrStates, "OUTPUT/gr_AE_list.rds")

#---------------------------------------------------------------------------
Patterns<-read.table("OUTPUT/Table_PatternsAE.txt", sep="\t")

library(superheat)
library(RColorBrewer)
# par(margin(3,5,5,10))

#Figure with heatmap on the left showing the celltypes in the
#first top 100 patterns and barplot in the right showing the 
#abundance of the pattern in the enhancer dataset

data<-Patterns[2:100,1:27]
rownames(data)<-Patterns$Patron[2:100]

# rownames(data)<-NULL
superheat(
  data,
  col.dendrogram = TRUE,
  smooth.heat = FALSE,
  # set heatmap color map
  # heat.pal = c("white", "#542788"),
  # heat.na.col = "white",
  
  # grid line colors
  grid.vline.col = "lightgrey",
  grid.hline.col = "lightgrey",
  
  # right plot: HDI
  yr = Patterns$freq[2:100] / 10000,
  yr.plot.type = "bar",
  yr.axis.name = "Number of AE bins\n(*10^4)",
  yr.plot.size = 0.6,
  yr.bar.col = "#542788",
  yr.obs.col = rep("#542788", 99),
  yr.axis.name.angle = 90,
  
  # # left labels
  left.label = "none",
  # left.label.size = 1,
  # left.label.text.size = 3,
  # left.label.col = adjustcolor("grey", alpha.f = 0.3),
  
  # bottom labels
  bottom.label.size = 0.5,
  bottom.label.text.size = 3,
  bottom.label.col = "white",
  bottom.label.text.angle = 90,
  bottom.label.text.alignment = "right"
)

#Biplot Clustering rows and columns of the patterns

data<-Patterns[2:100,1:27]
rownames(data) <- Patterns$Patron[2:100]

superheat(
  data,
  col.dendrogram = TRUE,
  row.dendrogram = TRUE,
  smooth.heat = FALSE,
  
  # grid line colors
  grid.vline.col = "lightgrey",
  grid.hline.col = "lightgrey",
  
  # # left labels
  # left.label="none",
  left.label.size = 0.5,
  left.label.text.size = 3,
  # left.label.col = adjustcolor("grey", alpha.f = 0.3),
  
  # bottom labels
  bottom.label.size = 0.5,
  bottom.label.text.size = 3,
  bottom.label.col = "white",
  bottom.label.text.angle = 90,
  bottom.label.text.alignment = "right"

)


library(pheatmap)

#Pattern frequency plot
anotcols <- data.frame(Freq = Patterns[2:100, "freq"] / 10000)

bars<-anotcols
bars$Pattern<-Patterns[2:100,"Patron"]
bars$Pattern<-factor(bars$Pattern, 
                     levels= rev(bars[order(bars$Freq),"Pattern"]))

bars$color<-"#e8b83f"
bars[bars$Freq>=2, "color"]<-"#de6c2a"
bars[bars$Freq>=4, "color"]<-"#728c1b"
bars[bars$Freq>=6, "color"]<-"#178562"


library(ggplot2)
ggplot(bars, aes(x=Pattern, y=Freq, fill=color)) +
  geom_bar(stat="identity", size=0)+theme_classic()+
  scale_fill_manual(values=rev(c("#e8b83f","#de6c2a", "#728c1b", "#178562")))+
  theme(axis.text.x = element_text(angle = 45,  hjust=0.98), #element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 12))

anotcols[anotcols<2]<-"level1"
anotcols[anotcols>6]<-"level4"
anotcols[anotcols<4]<-"level2"
anotcols[anotcols<6]<-"level3"
rownames(anotcols) <- Patterns$Patron[2:100]
anotcols$Freq<-factor(anotcols$Freq, levels= c("level1", "level2","level3", "level4"))

ann_colors = c("lightgrey", "#8c1b2c", "#2accde", "#178562")
names(ann_colors)<-c("level1", "level2","level3", "level4")
ann_colors<-list(Freq=ann_colors)

#biplot
#1000 * 370 fig dimensions
pheatmap(
  t(data),
  color = c("white", "#ad2487"),
  border_color = "white",
  scale = "none",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  treeheight_col = 30,
  treeheight_row = 25,
  show_colnames = FALSE,
  fontsize_col = 8,
  annotation_col = anotcols,
  annotation_colors = ann_colors
)

pheatmap
# clustering_distance_rows = "euclidean", 
# clustering_distance_cols = "euclidean", 
# clustering_method = "complete"

str(hclust(dist(data, method = "euclidean"), method="complete"))
cellclust<-hclust(dist(t(data), method = "euclidean"), method="complete")
cellclust<-cellclust$labels[c(cellclust$order)]
write.table(cellclust, "OUTPUT/cellclust.txt", sep="\t")

patternclust<-hclust(dist(data, method = "euclidean"), method="complete")
patternclust<-patternclust$labels[c(patternclust$order)]
write.table(patternclust, "OUTPUT/patternclust.txt", sep="\t")

#--------------------------------------------------------------------------
#                     PLOT ENHANCERS BY CELLTYPE
#--------------------------------------------------------------------------

gr_CtypeChrStates<-readRDS("OUTPUT/gr_AE_list.rds")
library(tidyverse)
CtypeChrStates.DF<-gr_CtypeChrStates %>% data.frame 


All.enhancers<-lapply(CtypeChrStates.DF[,6:32], 
       function(x) sum(as.numeric(as.character(x))))

All.enhancers<-do.call(rbind, All.enhancers)


specific.enhancers <-
  lapply(CtypeChrStates.DF[, 6:32], function(x)
    as.numeric(as.character(x)))

specific.enhancers <-
  do.call(cbind, specific.enhancers)

CtypeChrStates.DF$specific <- rowSums(specific.enhancers) == 1

sumas<-CtypeChrStates.DF %>% filter(specific == TRUE)  
sumas<-lapply(sumas[,6:32], 
       function(x) sum(as.numeric(as.character(x))))
spec.enhancers<-do.call(rbind, sumas)
spec.enhancers<-data.frame(celltype=rownames(spec.enhancers),all=All.enhancers, specific=spec.enhancers)
spec.enhancers.long<-spec.enhancers
spec.enhancers.long$all<-spec.enhancers.long$all-spec.enhancers.long$specific
spec.enhancers.long<-pivot_longer(spec.enhancers, cols = 2:3, names_to = "type", values_to = "NEnhancers")

level<-spec.enhancers.long[order(spec.enhancers.long$NEnhancers, decreasing = TRUE),] %>%
  select(celltype)%>% unique

ggplot(spec.enhancers.long, 
       aes(fill=factor(type, levels = c( "all", "specific")),
                y=NEnhancers, x=factor(celltype, levels = level$celltype ))) + 
  geom_bar(position="stack", stat="identity") + theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c("#194d82", "#8a214e"), name = c("all", "specific"))

