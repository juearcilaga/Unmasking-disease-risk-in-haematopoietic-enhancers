setwd("/mnt/nocode/juliana/PID/ThemeII")

chrlabels <-
  read.table("INPUTS/labelsChr.txt", header = TRUE, sep = "\t")

#Import chromatin states merged
All_data <- readRDS("INPUTS/states_220422_BACKGROUND_merged.rds")
dim(All_data)#15441337      108
head(All_data)

#NR1H3
#Chromosome 11: 47,248,300-47,269,033 forward strand.
#GRCh38:CM000673.2

start <- 47248300
end <- 47269033

Chrom <- rownames(All_data)[grepl("chr11", rownames(All_data))]
Region <- Chrom[ceiling(start / 200):ceiling(end / 200)]

RegionStates <- All_data[rownames(All_data) %in% Region, ]

# #Import metadata
metadata <- read.table("OUTPUT/METADATA_paper.txt", sep = "\t")
# length(unique(metadata$SAMPLE_NAME))#107
# #One sample no metadata

# mycelltypes <- unique(metadata$CELL_TYPE)#23
# 
# #Cell types Biola PC-HiC
# Biola_PCHiC_cells <- mycelltypes[c(5, 6, 7, 8, 10, 11, 12, 13, 16, 19, 23)]
# 
# mycelltypes <- metadata %>% filter(., CELL_TYPE %in% Biola_PCHiC_cells)

library(dplyr)
library(tidyverse)

###Consensus chromatin states 

# 75% of the samples should  be consistent
Sum <- function(x) {
  sum(x) >= round(length(x) * 0.75)
}

FindConserved <- function(chromstate) {
  # chromstate<-chrlabels$SymbState[2]
  
  celltype.chromstate <- celltype == chromstate
  
  ifelse(as.numeric(dim(data.frame(celltype))[2]) > 1,
         Suma <- apply(celltype.chromstate, 1, Sum),
         Suma <- celltype.chromstate)
  
  celltype2 <- data.frame(celltype, Suma)
  celltype2$conserved <- 0
  celltype2[celltype2$Suma == TRUE, "conserved"] <-  chromstate
  
  print(table(celltype2$conserved))
  conserved <-
    data.frame(Conserved = celltype2$conserved,
               row.names = rownames(celltype.chromstate))
  colnames(conserved) <- chrlabels$NameState[chromstate]
  return(conserved)
}

library(tidyr)
library(parallel)
library(MASS)

Conserved = list()

for (cell in unique(metadata$CELL_TYPE)) {
  #cell<-unique(metadata$CELL_TYPE)[27]#"unswitched memory B cell"
  print(cell)
  celltype <-
    RegionStates[, colnames(RegionStates) %in%
                   metadata[metadata$CELL_TYPE == cell, "SAMPLE_NAME"]]
  
  Conservedlist<-mclapply(chrlabels$SymbState, FindConserved, mc.cores = 35 )
  Conserved[[as.character(cell)]] <- do.call(cbind, Conservedlist)
}
# str(Conserved)

#Save this object Conserved
ConservedDF<-mclapply(Conserved, rowSums, mc.cores = 35 )
ConservedDF2<-data.frame(do.call(cbind, ConservedDF))
saveRDS(ConservedDF2, "OUTPUT/ConservedList_LiverX.rds")
ConservedDF2<-readRDS("OUTPUT/ConservedList_LiverX.rds")


# setwd("/mnt/nocode/juliana/PAR1")
# # saveRDS(ConservedDF2, "ConservedDF_Female_Par1.rds")

library(tidyverse)
ensayo <- ConservedDF2
ensayo$coord <- rownames(ensayo)
ensayo <-
  pivot_longer(
    ensayo,
    cols = 1:length(ensayo) - 1,
    names_to = "CELL_TYPE",
    values_to = "COLOR"
  )
ensayo <- transform(ensayo, ID = as.numeric(factor(CELL_TYPE)))
ensayo <- separate(ensayo, col = "coord", into = c("chr", "region"))

mydata = data.frame(
  variable = ensayo$CELL_TYPE,
  pos = ensayo$ID,
  start = (as.numeric(ensayo$region) * 200) - 199,
  end = as.numeric(ensayo$region) * 200,
  SymbState = ensayo$COLOR
)


# celltype_Abb <- data.frame(
#   variable = mydata$variable[1:11],
#   acronym = c(
#     "nCD4",
#     "nCD8",
#     "Neu",
#     "Mon",
#     "Ery",
#     "MK",
#     "M0",
#     "M2",
#     "M1",
#     "nB",
#     "EndP"
#   )
# )


# mydata <- left_join(mydata, celltype_Abb)
chrlabels$SymbState <- as.character(chrlabels$SymbState)
chrlabels$NameState <- as.character(chrlabels$NameState)
chrlabels[8, ] <- c("0", "Not conserved")
mydata$SymbState <- as.character(mydata$SymbState)
mydata <- left_join(mydata, chrlabels)

library(RColorBrewer)
library(scales)

chrcolor <- brewer.pal(n = 7, name = 'Dark2')
chrcolor[c(8)] <- c("grey")
chrcolor[6] <- "white"
level_order <-
  c(
    "Not conserved",
    "TransElo",
    "HetChrom",
    "PolyComb",
    "AProm",
    "AeloE",
    "AconE",
    "PoisE"
  )

names(chrcolor) <- rev(level_order)


ggplot(mydata,
       aes(
         x = start,
         xend = end,
         y = variable,
         yend = variable,
         color = NameState
       ),
       size = 3) +
  geom_segment(size = 5) +
  scale_color_manual(values = chrcolor) +
  theme_bw() +
  xlab("Position in chr") +
  ggtitle("Conserved chromatin states in region") +
  ylab("Celltype")

# saveRDS(mydata, "Region_chrstates_bycelltype.rds")
# write.table(mydata, "Region_chrstates_bycelltype.txt", sep = "\t")


mydata2<-mydata[grepl("macrophage", mydata$variable),]

mydata2$pos<-rep(c(2,1,3), length(mydata2$variable)/3)


ggplot(mydata2,
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
  xlab("Position in chr") +
  ggtitle("Conserved chromatin states in region") +
  ylab("Celltype")
#-----------------------------------------------------------

# by sample
#-----------------------------------------------------------

ensayo<-RegionStates
ensayo$coord<-rownames(ensayo)
ensayo<-pivot_longer(ensayo, cols = 1:length(ensayo)-1, names_to = "SAMPLES", values_to = "COLOR")
ensayo<-transform(ensayo, ID = as.numeric(factor(SAMPLES)))
ensayo<-separate(ensayo,col = "coord", into = c("chr", "region"))

mydata = data.frame(
  variable = ensayo$SAMPLES,
  pos = ensayo$ID,
  start = (as.numeric(ensayo$region)*200)-199,
  end = as.numeric(ensayo$region)*200,
  SymbState=ensayo$COLOR
)


# chrlabels$SymbState<-as.character(chrlabels$SymbState)
# chrlabels$NameState<-as.character(chrlabels$NameState)
# chrlabels[8,]<-c("0", "Not conserved")
mydata$SymbState<-as.character(mydata$SymbState)
mydata<-left_join(mydata, chrlabels)

# library(RColorBrewer)
# library(scales)

# chrcolor <- brewer.pal(n = 7, name = 'Dark2')
# chrcolor[c(8,9,10)]<- c("#D95F02", "darkorchid", "grey")
# chrcolor[6]<-"white"
# level_order <-
#   c("Other",
#     "protein_coding",
#     "Not conserved",
#     "TransElo",
#     "HetChrom",
#     "PolyComb",
#     "AProm",
#     "AeloE",
#     "AconE",
#     "PoisE")

# names(chrcolor) <- rev(level_order)

ggplot(
  mydata,
  aes(
    x = start,
    xend = end,
    y = variable,
    yend =  variable,
    color = NameState), size = 3) + 
  geom_segment(size = 5) +
  scale_color_manual(values = chrcolor )+
  theme_bw()+ 
  xlab("Position in chrX") +
  ggtitle("Conserved chromatin states in Par1 region")+ 
  ylab("Celltype")


M012<-unique(metadata[grepl("macrophage", metadata$CELL_TYPE), "SAMPLE_NAME"])

ggplot(
  mydata[mydata$variable %in% M012, ],
  aes(
    x = start,
    xend = end,
    y = variable,
    yend =  variable,
    color = NameState), size = 3) + 
  geom_segment(size = 5) +
  scale_color_manual(values = chrcolor )+
  theme_bw()+ 
  xlab("Position in chr 11") +
  ggtitle("Conserved chromatin states in LiverX receptor region")+ 
  ylab("Celltype")+
  geom_vline(xintercept = 47253513, linetype="dotted", 
             color = "black", size=0.5)#+
  # geom_vline(xintercept = 68159740, linetype="dotted", 
  #            color = "black", size=0.5)


47253513
