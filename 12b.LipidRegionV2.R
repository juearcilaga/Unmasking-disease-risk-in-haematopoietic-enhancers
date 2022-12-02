#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
setwd("/mnt/nocode/juliana/PID/ThemeII")
library(dplyr)
library(tidyverse)
library(tidyr)
library(parallel)
library(MASS)
library(RColorBrewer)
library(scales)
#-----------------------------------------------------------------------------------------
#Import Genes and SNPS related to lipid metabolism
#-----------------------------------------------------------------------------------------
#Gene coordninates (genes in enhancers snps regions)
lipid.GENES.DF<-read.table("OUTPUT/Lipid.GENES.DF.txt", sep="\t")

#Enhancer coordinates (enhancers snps regions)
SNPsToEnhM0123.summary.bychr<-read.table("OUTPUT/SNPsToEnhM0123.summary.bychr.txt", sep="\t")

#Snps coordinates
SNPsToEnhM0123.summary<-read.table("OUTPUT/SNPsToEnhM0123summary.txt", sep="\t")

##Chromatin state labels for the plot
chrlabels <-read.table("INPUTS/labelsChr.txt", header = TRUE, sep = "\t")

#Import chromatin states merged
All_data <- readRDS("INPUTS/states_220422_BACKGROUND_merged.rds")

# #Import metadata
metadata <- read.table("OUTPUT/METADATA_paper.txt", sep = "\t")
metadata2<-read.table("../SecondYearAssesment/RESULTADOS/celltype_data.txt", sep = "\t", header = TRUE)

FixLabels <- function(vector.of.names) {
  fixed.names <- gsub(" ", ".", vector.of.names)
  fixed.names <- gsub("-", ".", fixed.names)
  fixed.names <- gsub(",", ".", fixed.names)
  fixed.names <- gsub("\\(", ".", fixed.names)
  fixed.names <- gsub("\\)", ".", fixed.names)
  fixed.names <- gsub("\\.+", ".", fixed.names)
  return(fixed.names)
}

metadata2$CELL_TYPE<-FixLabels(metadata2$CELL_TYPE)
metadata$CELL_TYPE<-factor(FixLabels(metadata$CELL_TYPE), levels=metadata2$CELL_TYPE)
metadata<-with(metadata, metadata[order(CELL_TYPE),])
#-----------------------------------------------------------------------------------------
# Lipid region
#Chromosome 8:#19954801 - 20012600
#GRCh38:
#-----------------------------------------------------------------------------------------
# start <- 19954801
# end <- 20012600
# chromosome<-"chr8"

#Only regions with genes??c("11","19","22","8","X")

#Set parameters for 1:length(SNPsToEnhM0123.summary.bychr$start20k)
##Falta el 5 chr  6
SetofEnh<-5
chromosome<-as.character(SNPsToEnhM0123.summary.bychr$chr[SetofEnh])
start<-SNPsToEnhM0123.summary.bychr$start20k[SetofEnh]
end <- SNPsToEnhM0123.summary.bychr$end20k[SetofEnh]

Chrom <- rownames(All_data)[grepl(chromosome, rownames(All_data))]
Region <- Chrom[ceiling(start / 200):ceiling(end / 200)]
RegionStates <- All_data[rownames(All_data) %in% Region, ]

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

filename<-sprintf("OUTPUT/ConservedList_Lipids_%s.rds", chromosome)
saveRDS(ConservedDF2, filename)
#-----------------------------------------------------------------------------------------
# Formating data for the plot
#-----------------------------------------------------------------------------------------
# chromosome<- XXXX
# filename<-sprintf("OUTPUT/ConservedList_Lipids_%s.rds", chromosome)
# ConservedDF2<-readRDS("OUTPUT/ConservedList_Lipids.rds")

## Formating longer instead of wider
ConservedDF.longer <- ConservedDF2
ConservedDF.longer$coord <- rownames(ConservedDF.longer)
ConservedDF.longer <-
  pivot_longer(
    ConservedDF.longer,
    cols = 1:length(ConservedDF.longer) - 1,
    names_to = "CELL_TYPE",
    values_to = "COLOR"
  )
ConservedDF.longer$CELL_TYPE<-factor(ConservedDF.longer$CELL_TYPE, 
                                     levels = levels(metadata$CELL_TYPE))
#levels(ConservedDF.longer$CELL_TYPE)
ConservedDF.longer <- transform(ConservedDF.longer, ID = as.numeric(factor(CELL_TYPE)))
ConservedDF.longer <- tidyr::separate(ConservedDF.longer, col = "coord", into = c("chr", "region"))


# Formating coordintaes for the segment plot

mydata = data.frame(
  variable = ConservedDF.longer$CELL_TYPE,
  pos = ConservedDF.longer$ID,
  start = (as.numeric(ConservedDF.longer$region) * 200) - 199,
  end = as.numeric(ConservedDF.longer$region) * 200,
  SymbState = ConservedDF.longer$COLOR
)

#-----------------------------------------------------------------------------------------
# Seting colors for the chromatin states in the plot
#-----------------------------------------------------------------------------------------
chrlabels$SymbState <- as.character(chrlabels$SymbState)
chrlabels$NameState <- as.character(chrlabels$NameState)
chrlabels[8, ] <- c("0", "Not conserved")
chrlabels[9, ] <- c("8", "Gene")

chrcolor <- brewer.pal(n = 7, name = 'Dark2')
chrcolor[6] <- c("grey")
chrcolor[8] <- "#3f4963"
chrcolor[9] <- "black"

level_order <- c("gene", "Not conserved","TransElo","HetChrom",
                 "PolyComb","AProm","AeloE","AconE","PoisE")
names(chrcolor) <- rev(level_order)

## Seting the snp positions as vertical lines
vertical.lines<-SNPsToEnhM0123.summary[SNPsToEnhM0123.summary$chr==chromosome,"CHR_POS"]

##Adding labels to the formated data for the plot
mydata$SymbState <- as.character(mydata$SymbState)
mydata <- left_join(mydata, chrlabels)

InRegionGenes<-lipid.GENES.DF[lipid.GENES.DF$seqnames==chromosome,]
InRegionGenes<-  
  data.frame(
    variable = InRegionGenes$gene_name,
    pos = max(mydata$pos) +1,
    start = as.numeric(InRegionGenes$start),
    end = as.numeric(InRegionGenes$end),
    SymbState ="8",
    NameState = "gene"
  )

mydata.genes<- rbind(mydata,InRegionGenes)
#-------------------------------------------------------------------------------------------
##Graph for all the celltypes
#-------------------------------------------------------------------------------------------
Chromosomename<-sprintf("Position in %s", chromosome)

ggplot(mydata.genes,
       aes(x = start,
           xend = end,
           y = variable,
           yend = variable,
           color = NameState), 
       size = 3) +
  geom_segment(size = 5) +
  scale_color_manual(values = chrcolor) +
  theme_bw() +
  xlab(Chromosomename) +
  ggtitle("Conserved chromatin states by celltype") +
  ylab("cell type")+
  geom_vline(xintercept = vertical.lines, linetype="dotted", 
             color = "white", size=0.5)
# Size of the plot 1180 * 500


#-----------------------------------------------------------------------------------------
#Graph only macrophages
#-----------------------------------------------------------------------------------------
mydata2<-mydata.genes[grepl("macrophage", mydata.genes$variable),]
mydata2$pos<-rep(c(3,2,4), length(mydata2$variable)/3)
mydata2.genes<- rbind(mydata2,InRegionGenes)
mydata2.genes[mydata2.genes$variable=="Genes", "pos"]<-1

ggplot(mydata2.genes,
       aes(
         x = start,
         xend = end,
         y = variable,
         yend = variable,
         color = NameState
       ),
       size = 12) +
  geom_segment(size = 7) +
  scale_color_manual(values = chrcolor) +
  theme_bw() +
  xlab(Chromosomename) +
  ggtitle("Conserved chromatin states in region") +
  ylab("Celltype")+
  geom_vline(xintercept = vertical.lines, linetype="dotted", 
             color = "white", size=0.5)

# Size of the plot 1180 * 200
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

rm(InRegionGenes.sample)

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
  geom_segment(size = 5) +
  scale_color_manual(values = chrcolor )+
  theme_bw()+ 
  xlab(Chromosomename) +
  ggtitle("Conserved chromatin states by sample")+ 
  ylab("Cell type")+
  geom_vline(xintercept = vertical.lines, linetype="dotted", 
                        color = "white", size=0.5)+
  theme(axis.text.y = element_text(colour = coloraxis))



rm(mydata3.genes)
#Size 1180 * 450
# Size of the plot 1180 * 200
# Size of the plot 1180 * 550

