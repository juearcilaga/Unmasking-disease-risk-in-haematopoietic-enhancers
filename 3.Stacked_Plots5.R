#--------------------------------------------------------------------------------
#     INPUTS
#--------------------------------------------------------------------------------
setwd("/mnt/nocode/juliana/PID/ThemeII")

#Import chromatin states merged
All_data<-readRDS("INPUTS/states_220422_BACKGROUND_merged.rds")

#Import metadata
metadata<-read.table("OUTPUT/METADATA_paper.txt", sep="\t")

#Import celltypes ontology
celltypes<-read.delim("INPUTS/celltype_data.txt", sep="\t")

#Import chromatin states labels and identifiers
chrlabels<-read.table("INPUTS/labelsChr.txt", header = TRUE,sep="\t")
#--------------------------------------------------------------------------------
library(dplyr)
library(tidyr)

mylineage<-All_data

my_list <- list()
for (i in 1:nrow(chrlabels)) {
  Chromstate <- as.character(chrlabels$NameState[i])
  SymbState <- chrlabels$SymbState[i]
  my_list[[Chromstate]] <- colSums(mylineage == SymbState)
}

StackObj <- do.call(cbind, my_list)
StackObj <- data.frame(SAMPLE_NAME = rownames(StackObj), StackObj)

StackObj.long <-
  pivot_longer(
    StackObj,
    cols = 2:ncol(StackObj),
    names_to = "NameState",
    values_to = "NBins",
  )

StackObj.long<-
  left_join(StackObj.long, metadata[, c("SAMPLE_NAME", "CELL_TYPE", "TISSUE_TYPE")])

StackObj.long <- left_join(StackObj.long, chrlabels)
StackObj.long$SymbState <- as.integer(StackObj.long$SymbState)
StackObj.long<- left_join(StackObj.long,
                      celltypes[, c("MAIN_CELL_TYPE", "CELL_TYPE_Abb", "CELL_TYPE")])

StackObj.long$MAIN_CELL_TYPE <-
  factor(StackObj.long$MAIN_CELL_TYPE, 
  level = c("mesenchymal","endothelial","T",
          "B","plasma","NK","neutrophilic",
          "neutrophil","eosinophil","monocyte", 
          "osteoclast","dendritic", "macrophage",
          "megakaryocyte", "erythroblast"))


StackObj.long$TISSUE_TYPE <-
  factor(StackObj.long$TISSUE_TYPE,
         level = rev(unique(StackObj.long$TISSUE_TYPE)))

StackObj.long$TISSUE_TYPE <-
  factor(StackObj.long$TISSUE_TYPE,
         level = rev(unique(StackObj.long$TISSUE_TYPE)))

StackObj.long$SIZE_bp<-as.numeric(StackObj.long$NBins*200)

StackObj.long<-StackObj.long[StackObj.long$SAMPLE_NAME!="S01Y2UH1",]

#write.table(StackObj.long, "OUTPUT/Enhancerbund.txt", sep="\t")
write.table(StackObj.long, "OUTPUT/Enhancerbund_all.txt", sep="\t")
#-------------------------------------------  -------------------------------------
library(RColorBrewer)
library(scales)
chrcolor <- brewer.pal(n = 7, name = 'Dark2')
level_order <-
  c("TransElo",
    "HetChrom",
    "PolyComb",
    "AProm",
    "AeloE",
    "AconE",
    "PoisE")

names(chrcolor) <- rev(level_order)

data <- StackObj.long
# data <- StackObj.long[StackObj.long$NameState%in% c("AconE","AeloE","PoisE"),]

# Stacked
ggplot(data, aes(fill=factor(NameState, level = rev(level_order)),
                 y=NBins, x=SAMPLE_NAME)) + 
  geom_bar(position="stack", stat="identity") + theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_grid(
    ~ MAIN_CELL_TYPE,
    scales = "free_x", #Let the x axis vary across facets
    space = "free_x", #Let the width of facets vary and force all bars to have the same width.
    switch = "x")  +
  theme(
    strip.placement = "outside", #Place facet labels outside x axis labels.
    strip.background = element_blank(), #Make facet label background white.
    axis.title = element_blank(),
    text = element_text(size = 10)) +
  scale_fill_manual(values = chrcolor, name = "ChrState")


#--------------------------------------------------------------------------------
#                               Gained enhancers 
#--------------------------------------------------------------------------------
monodescend <- mylineage

GainLostStacked <- function(Enhtype) {
  #Lost
  ##All regions where both monocytes have Enhtype enhancers
  MonocytePE <-
    monodescend[apply(monodescend[, Monocytos], 1, function(r)
      all(r %in% c(Enhtype))),]
  dim(MonocytePE)#[1] 193749     12
  
  ##How many AE conservados in the descendants? Each cell type each sample
  PE_lost <- colSums(MonocytePE != Enhtype)
  
  summary(PE_lost)
  # Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
  # 0   75441     107,808   91920  122091  154879
  
  #Gained
  MonocytenoPE <-
    monodescend[apply(monodescend[, Monocytos], 1, function(r)
      all(r %in% c(Enhtype)) == FALSE), ]
  dim(MonocytenoPE)#[1] 15247588       12
  
  PE_gained <- colSums(MonocytenoPE == Enhtype)
  
  summary(PEneweachsample[3:length(PEneweachsample)])
  # Min. 1st Qu.     Median    Mean 3rd Qu.    Max.
  # 116,355  242206  326,596  327122  417170  525,677
  
  Gained_Poise <-
    data.frame(SAMPLE_NAME = names(PE_gained), NBins = PE_gained)
  Lost_Poise <-
    data.frame(SAMPLE_NAME = names(PE_lost), NBins = -PE_lost)
  GainedLost <- rbind(Gained_Poise, Lost_Poise)
  
  GainedLost <-
    left_join(GainedLost, metadata[, c("SAMPLE_NAME", "CELL_TYPE")])
  
  GainedLost <- left_join(GainedLost,
                          celltypes[, c("MAIN_CELL_TYPE", "CELL_TYPE_Abb", "CELL_TYPE")])
  GainedLost$MAIN_CELL_TYPE <-
    factor(GainedLost$MAIN_CELL_TYPE, level = rev(unique(celltypes$MAIN_CELL_TYPE)))
  
  GainedLost$SIZE_bp <- as.numeric(GainedLost$NBins * 200)
  
  return(GainedLost)
}

GainedLost_Plot <- function(GainedLost) {
  GainedLost_Plot_Poised <-
    ggplot(GainedLost[GainedLost$MAIN_CELL_TYPE != "monocyte",],
           aes(
             y = NBins / 10 ^ 3,
             x = SAMPLE_NAME,
             fill = NBins > 0
           )) +
    geom_bar(stat = "identity") + theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    facet_grid(
      ~ MAIN_CELL_TYPE,
      scales = "free_x",
      #  Let the x axis vary across facets.
      space = "free_x",
      # Let the width of facets vary and force all bars to have the same width.
      switch = "x"
    )  +
    theme(
      strip.placement = "outside",
      # Place facet labels outside x axis labels.
      strip.background = element_blank(),
      # Make facet label background white.
      axis.title = element_blank(),
      text = element_text(size = 10)
    ) +
    geom_abline(
      slope = 0,
      intercept = 0,
      col = "grey",
      lty = 2
    ) +
    scale_y_continuous(limits = c(-200, 600), labels = unit_format(unit = "Kb"))
  return(GainedLost_Plot_Poised)
}


GainedLostPoised<-GainLostStacked(5)
write.table(GainedLostPoised, "OUTPUT/GainedLostPoised.txt",sep="\t")
GainedLostPoisedPlot<-GainedLost_Plot(GainedLostPoised)
GainedLostPoisedPlot<-GainedLostPoisedPlot+ 
  scale_fill_manual(values = setNames(c('#1B9E77', '#115c45'), c(T, F))) 
 
GainedLostAconE<-GainLostStacked(6)
write.table(GainedLostAconE, "OUTPUT/GainedLostAconE.txt",sep="\t")
GainedLostAconEPlot<-GainedLost_Plot(GainedLostAconE)
GainedLostAconEPlot<-GainedLostAconEPlot+ 
  scale_fill_manual(values = setNames(c('#D95F02', '#a33c08'), c(T, F))) 
#-------------------------------------------  -------------------------------------
# We found that poised enhancers (median 407,610 bin) are most 
# abundant than active enhancers (median 207,077 bin) in all samples.
# In comparison to monocytes, descendant celltypes win more new enhancers 
# (median new AE:163,866 bins, median new PE:326,596,866 bin) than they loss 
# (median lost Active Enhancers:30,718 bins, median lost Poised Enhancers:107,808 bins)
