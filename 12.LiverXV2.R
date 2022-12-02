setwd("/mnt/nocode/juliana/PID/ThemeII")

library(dplyr)
library(tidyverse)
library(tidyr)
library(parallel)
library(MASS)

chrlabels <-
  read.table("INPUTS/labelsChr.txt", header = TRUE, sep = "\t")

#Import chromatin states merged
All_data <- readRDS("INPUTS/states_220422_BACKGROUND_merged.rds")
dim(All_data)#15441337      108
head(All_data)

#NR1H3
#Chromosome 11: 47,248,300-47,269,033 forward strand.
#GRCh38:CM000673.2

start <- 47248300 - 5000
end <- 47269033 + 5000

Chrom <- rownames(All_data)[grepl("chr11", rownames(All_data))]
Region <- Chrom[ceiling(start / 200):ceiling(end / 200)]

RegionStates <- All_data[rownames(All_data) %in% Region, ]

# #Import metadata
metadata <- read.table("OUTPUT/METADATA_paper.txt", sep = "\t")


#-----------------------------------------------------------------------------------------
###Consensus chromatin states 
#-----------------------------------------------------------------------------------------
# 75% of the samples should  be consistent
Sum <- function(x) {
  sum(x) >= round(length(x) * 0.75)
}

FindConserved <- function(chromstate) {
  # chromstate<-chrlabels$SymbState[1]
  celltype.chromstate <- celltype == chromstate
  ifelse(as.numeric(dim(data.frame(celltype))[2]) > 1,
         Suma <- apply(celltype.chromstate, 1, Sum),
         Suma <- celltype.chromstate)
  
  celltype2 <- data.frame(celltype, Suma)
  celltype2$conserved <- 0
  celltype2[celltype2$Suma == TRUE, "conserved"] <-  chromstate
  print(table(celltype2$conserved))
  conserved <-
    data.frame(Conserved = as.numeric(celltype2$conserved),
               row.names = rownames(celltype.chromstate))
  colnames(conserved) <- chrlabels$NameState[as.numeric(chromstate)]
  return(conserved)
}

Conserved = list()

for (cell in unique(metadata$CELL_TYPE)) {
  print(cell)
  celltype <-
    RegionStates[, colnames(RegionStates) %in%
                   metadata[metadata$CELL_TYPE == cell, "SAMPLE_NAME"]]
  
  Conservedlist<-mclapply(chrlabels$SymbState, FindConserved, mc.cores = 35 )
  Conserved[[as.character(cell)]] <- do.call(cbind, Conservedlist)
}

ConservedDF<-mclapply(Conserved, rowSums, mc.cores = 35 )

ConservedDF2<-data.frame(do.call(cbind, ConservedDF))

saveRDS(ConservedDF2, "OUTPUT/ConservedList_LiverX.rds")

#---------------------------------------------------------------------------------
# PLOT
#---------------------------------------------------------------------------------
# ConservedDF2<-readRDS("OUTPUT/ConservedList_LiverX.rds")

celltypes_Abb<-read.delim("OUTPUT/celltype_Abb.txt")
celltypes<-read.delim("INPUTS/celltype_data2.txt")

library(RColorBrewer)
library(scales)


Conserved.Long <- ConservedDF2

Conserved.Long$coord <- rownames(Conserved.Long)

Conserved.Long <-
  pivot_longer(
    Conserved.Long,
    cols = 1:length(Conserved.Long) - 1,
    names_to = "CELL_TYPE",
    values_to = "COLOR"
  )

Conserved.Long <- transform(Conserved.Long, ID = as.numeric(factor(CELL_TYPE)))

Conserved.Long <- separate(Conserved.Long, col = "coord", into = c("chr", "region"))

mydata = data.frame(
  variable = Conserved.Long$CELL_TYPE,
  pos = Conserved.Long$ID,
  start = (as.numeric(Conserved.Long$region) * 200) - 199,
  end = as.numeric(Conserved.Long$region) * 200,
  SymbState = Conserved.Long$COLOR
)

mydata <-left_join(mydata, celltypes_Abb)
mydata$acronym<-factor(mydata$acronym, levels=rev(unique(celltypes$PATTERN_NAME)))

mydata$SymbState <- as.character(mydata$SymbState)
chrlabels$SymbState <- as.character(chrlabels$SymbState)
chrlabels$NameState <- as.character(chrlabels$NameState)
chrlabels[8, ] <- c("0", "Not conserved")
chrlabels[9, ] <- c("8", "Gene")

mydata <- left_join(mydata, chrlabels)


chrcolor <- brewer.pal(n = 7, name = 'Dark2')
chrcolor[6] <- c("grey")
chrcolor[8] <- "#3f4963"
chrcolor[9] <- "black"

level_order <- c("gene", "Not conserved","TransElo","HetChrom",
                 "PolyComb","AProm","AeloE","AconE","PoisE")
names(chrcolor) <- rev(level_order)

#Gene coordninates (genes in enhancers snps regions)
lipid.GENES.DF<-read.table("OUTPUT/Lipid.GENES.DF.txt", sep="\t")

InRegionGenes<-lipid.GENES.DF[lipid.GENES.DF$seqnames=="chr11",]

InRegionGenes<-  
  data.frame(
    variable = InRegionGenes$gene_name,
    pos = max(mydata$pos) +1,
    start = as.numeric(InRegionGenes$gene_start),
    end = as.numeric(InRegionGenes$gene_end),
    SymbState ="8",
    NameState = "gene",
    acronym = InRegionGenes$gene_name
    )

mydata.genes<- rbind(mydata,InRegionGenes)

ggplot(
  mydata.genes,
  aes(
    x = start,
    xend = end,
    y =  acronym,
    yend =  acronym,
    color = NameState
  ),
  size = 3
) +
  geom_segment(size = 3.5) +
  scale_color_manual(values = chrcolor) +
  theme_bw() +
  xlab("Position in chr 11") +
  ggtitle("Conserved chromatin states in NR1H3") +
  ylab("cell type") + theme(text = element_text(size = 12),
                            axis.text = element_text(color = "black"))

# saveRDS(mydata, "Region_chrstates_bycelltype.rds")
# write.table(mydata, "Region_chrstates_bycelltype.txt", sep = "\t")
# 639 * 417

#-----------------------------------------------------------
# by sample
#-----------------------------------------------------------
ConservedDF.longer2<-RegionStates
ConservedDF.longer2$coord<-rownames(ConservedDF.longer2)
ConservedDF.longer2<-pivot_longer(ConservedDF.longer2, cols = 1:length(ConservedDF.longer2)-1, names_to = "SAMPLE_NAME", values_to = "COLOR")
ConservedDF.longer2<-left_join(ConservedDF.longer2, metadata)
ConservedDF.longer2<-with(ConservedDF.longer2, ConservedDF.longer2[order(CELL_TYPE),])

ConservedDF.longer2<-transform(ConservedDF.longer2, ID = as.numeric(factor(SAMPLE_NAME)))
ConservedDF.longer2<-separate(ConservedDF.longer2,col = "coord", into = c("chr", "region"))


mydata3 = data.frame(
  variable = ConservedDF.longer2$SAMPLE_NAME,
  # pos = ConservedDF.longer2$ID,
  start = (as.numeric(ConservedDF.longer2$region)*200)-199,
  end = as.numeric(ConservedDF.longer2$region)*200,
  SymbState=ConservedDF.longer2$COLOR,
  celltype=ConservedDF.longer2$CELL_TYPE
)

mydata3$SymbState<-as.character(mydata3$SymbState)
mydata3<-left_join(mydata3, chrlabels)

M012<-unique(metadata[grepl("macrophage", metadata$CELL_TYPE), "SAMPLE_NAME"])

mydata3.genes<-mydata3[mydata3$variable %in% M012, ]


InRegionGenes.sample<-data.frame(InRegionGenes[,c(1,3:5)], 
                                 celltype = "NA", 
                                 NameState = InRegionGenes[,6])

mydata3.genes<- rbind(mydata3.genes,InRegionGenes.sample)

# rm(InRegionGenes.sample)

positions<-mydata3.genes%>%dplyr::select(variable, celltype) %>% unique()
positions$pos<-0
positions[positions$celltype!="NA", "pos"]<-2:(length(positions[positions$celltype!="NA", "pos"])+1)
positions[positions$celltype=="NA", "pos"]<- 1

mydata3.genes<-left_join(mydata3.genes, positions)
# mydata3.genes<-with(mydata3.genes, mydata3.genes[order(celltype),])
mydata3.genes$variable<-factor(mydata3.genes$variable, 
                               levels=unique(mydata3.genes$variable))

mydata3.genes%>%dplyr::select(variable, celltype,pos) %>% unique()

coloraxis<-c(rep("#6b2520",7), rep("#3a5431",7), 
             rep("#1f5387",7), 
             rep("black",length(positions[positions$celltype=="NA", "pos"])))
names(coloraxis)<-positions$variable

ggplot(
  mydata3.genes,
  aes(
    x = start,
    xend = end,
    y = variable,
    yend =  variable,
    color = NameState), size = 3) + 
  geom_segment(size = 3.5) +
  scale_color_manual(values = chrcolor )+
  theme_bw()+ 
  xlab("Chromosomename") +
  # ggtitle("Conserved chromatin states by sample")+ 
  ylab("Cell type")+
  theme(text = element_text(size = 12),
        axis.text.y = element_text(colour = coloraxis))

#Size 639*317

