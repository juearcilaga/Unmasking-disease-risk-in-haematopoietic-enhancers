#-----------------------------------------------------------------
#               IMPORTING M1 and M2 ENHANCERS 
#   (conserved in the 3 replicates of each cell type)
#------------------------------------------------------------------
#Which enhacers are GWAS enriched
setwd("/mnt/nocode/juliana/PID/ThemeII")

#Chromstates conserved in 75% of the replicates
conserved<-readRDS("OUTPUT/ConservedDF.rds")

# AE in at least 2 samples
AllAE<-readRDS("OUTPUT/BActiveE.rds")

##Only AE
All_dataAE<-conserved[rownames(conserved)%in%AllAE,]
#Maybe some enhancers are not conserved in the celltype?

##Rename if AE or no
data_modif<-as.matrix(All_dataAE)
data_modif[data_modif %in% c(1,2,3,7,5)] <- 0#"No_E", 5 "PoisE"
data_modif[data_modif %in% c(4,6)] <- 1#"AE

colnames(data_modif)<-gsub("\\.", "", colnames(data_modif))

library(plyr)
Patterns<-plyr::count(df = data.frame(data_modif), vars = colnames(data_modif))
Patterns<-Patterns[order(Patterns$freq, decreasing = TRUE),]
length(Patterns$freq)#12901
colnames(Patterns)<-c("meos", "mDC", "imDC", "osteo", "mono",
                      "M0", "M2", "M1", "mneut", "CD38negnaiveB",
                      "cswitmB", "naiveB","CD4T","effmemCD8T", "mesenBM",
                      "endprog","freq")

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

Patterns<-read.table("OUTPUT/Table_PatternsAE.txt", sep="\t")
library(superheat)
library(RColorBrewer)

par(margin(3,5,5,10))
data<-Patterns[2:50,1:16]
rownames(data)<-Patterns$Patron[2:50]
# rownames(data)<-NULL
  superheat(data,
            col.dendrogram=TRUE,
            smooth.heat=FALSE,
            # set heatmap color map
            # heat.pal = c("white", "#542788"),
            # heat.na.col = "white",
        
            # grid line colors
            grid.vline.col = "lightgrey",
            grid.hline.col = "lightgrey",
            
            # right plot: HDI
            yr = Patterns$freq[2:50]/10000,
            yr.plot.type = "bar",
            yr.axis.name = "Number of AE bins\n(*10^4)",
            yr.plot.size = 0.6,
            yr.bar.col = "#542788",
            yr.obs.col = rep("#542788", 49),
            yr.axis.name.angle = 90,
            
            # # left labels
            left.label="none",
            # left.label.size = 1,
            # left.label.text.size = 3,
            # left.label.col = adjustcolor("grey", alpha.f = 0.3),
            
            # bottom labels
            bottom.label.size = 0.5,
            bottom.label.text.size = 3,
            bottom.label.col = "white",
            bottom.label.text.angle = 90,
            bottom.label.text.alignment = "right")

  data<-Patterns[2:100,1:16]
  rownames(data)<-Patterns$Patron[2:100]
  
  superheat(data,
            col.dendrogram=TRUE,
            row.dendrogram = TRUE,
            smooth.heat=FALSE,
            
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
            bottom.label.text.alignment = "right")
library(pheatmap)

anotcols<-data.frame(Freq=Patterns[2:100, "freq"]/10000)
rownames(anotcols)<-Patterns$Patron[2:100]
ann_colors = list(
  Freq=c("#e8b83f","#de6c2a", "#728c1b", "#178562"))

pheatmap(t(data),color = c("white", "#ad2487"), border_color= "white",
         scale="none", cluster_rows=TRUE, cluster_cols = TRUE, treeheight_col=30,
         treeheight_row = 25, show_colnames=FALSE, fontsize_col = 8, 
         annotation_col=anotcols, annotation_colors=ann_colors)

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

bars<-anotcols
bars$Pattern<-rownames(bars)
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
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 12))
#--------------------------------------------------------------------------

#####################################
Mylist<-data.frame(data_modif)

Milista<-as.matrix(Mylist)
colnames(Milista)<-c("meos", "mDC", "imDC", "osteo", "mono",
                    "M0", "M2", "M1", "mneut", "CD38negnaiveB",
                    "cswitmB", "naiveB","CD4T","effmemCD8T", "mesenBM",
                    "endprog")

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
