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
mycelltypes <-
  c("monocyte", "osteoclast", "immature.*dendritic.*")

Monocyte_lineage <-
  metadata[grepl(paste(mycelltypes, collapse = "|"), metadata$CELL_TYPE) | 
           metadata$CELL_TYPE=="macrophage", ]

# Monocyte_lineage <-metadata


Monocyte_lineage_venous <-
  Monocyte_lineage[Monocyte_lineage$TISSUE_TYPE == "venous blood", ]

Monocyte_lineage_venous$CELL_TYPE <-
  factor(Monocyte_lineage_venous$CELL_TYPE, levels = unique(celltypes$CELL_TYPE))

mylineage<-All_data[,colnames(All_data)%in%Monocyte_lineage_venous$SAMPLE_NAME]
#--------------------------------------------------------------------------------
library(dplyr)
library(tidyr)

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
  left_join(StackObj.long, metadata[, c("SAMPLE_NAME", "CELL_TYPE")])


StackObj.long <- left_join(StackObj.long, chrlabels)
StackObj.long$SymbState <- as.integer(StackObj.long$SymbState)
StackObj.long<- left_join(StackObj.long,
                      celltypes[, c("MAIN_CELL_TYPE", "CELL_TYPE_Abb", "CELL_TYPE")])
StackObj.long$MAIN_CELL_TYPE <-
  factor(StackObj.long$MAIN_CELL_TYPE, level = rev(unique(celltypes$MAIN_CELL_TYPE)))

StackObj.long$SIZE_bp<-as.numeric(StackObj.long$NBins*200)

#write.table(StackObj.long, "OUTPUT/Enhancerbund.txt", sep="\t")
write.table(StackObj.long, "OUTPUT/Enhancerbund_all.txt", sep="\t")
#--------------------------------------------------------------------------------
  Plot <- function(StackedEnh) {
    
    # StackedEnh<-Stacked_PoisE
    # StackedEnh<-Stacked_AeloE
    library(RColorBrewer)
    library(scales)
    chrcolor <- brewer.pal(n = 7, name = 'Dark2')
   
    chrcolor[8]<- "#D95F02"
    
    level_order <-
      c("AE",
      "TransElo",
        "HetChrom",
        "PolyComb",
        "AProm",
        "AeloE",
        "AconE",
        "PoisE")
  
    names(chrcolor) <- rev(level_order)
    
    Stacked_Plot <-
      ggplot(StackedEnh,
             aes(
               fill = factor(NameState, level = rev(level_order)),
               y = SIZE_bp/10^6,
               x = SAMPLE_NAME))+
      geom_bar(stat = "identity") + theme_classic() +
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
        scale_fill_manual(values = chrcolor, name = "ChrState") + 
      scale_y_continuous(
      limits = c(0, 130),
      labels = unit_format(unit = "Mb")
    )
    
    return(Stacked_Plot)
  }

Stacked_AconE<-
  StackObj.long[StackObj.long$NameState%in% c("AconE"),]

Stacked_PoisE<-
  StackObj.long[StackObj.long$NameState%in% c("PoisE"),]

Stacked_AeloE<-
  StackObj.long[StackObj.long$NameState%in% c("AeloE"),]

Stacked_AE<-
  StackObj.long[StackObj.long$NameState%in% c("AconE","AeloE"),]

Stacked_AE2<-aggregate(Stacked_AE$NBins, by=list(SAMPLE_NAME=Stacked_AE$SAMPLE_NAME), FUN=sum)
Stacked_AE<-left_join(Stacked_AE,Stacked_AE2)%>% dplyr::select(c(1,4,6,7,9))
colnames(Stacked_AE)[5]<-"NBins"
Stacked_AE$NameState<-"AE"
Stacked_AE$SIZE_bp<-Stacked_AE$NBins*200
Stacked_AE<-unique(Stacked_AE)

Plot(Stacked_AconE)
Plot(Stacked_PoisE)
Plot(Stacked_AeloE)
Plot(Stacked_AE)
#-------------------------------------------  -------------------------------------
library(RColorBrewer)
library(scales)
chrcolor <- brewer.pal(n = 7, name = 'Dark2')

chrcolor[8]<- "#D95F02"

level_order <-
  c("AE",
    "TransElo",
    "HetChrom",
    "PolyComb",
    "AProm",
    "AeloE",
    "AconE",
    "PoisE")

names(chrcolor) <- rev(level_order)

# data <- StackObj.long

data <- StackObj.long[StackObj.long$NameState%in%c("AconE", "AeloE","PoisE"),]
data <-Stacked_AE
# Stacked
ggplot(data, aes(fill=factor(NameState, level = rev(level_order)),
                 y=NBins, x=SAMPLE_NAME)) + 
  geom_bar(position="stack", stat="identity") + theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_grid(
    ~ CELL_TYPE_Abb,
    scales = "free_x", #Let the x axis vary across facets
    space = "free_x", #Let the width of facets vary and force all bars to have the same width.
    switch = "x")  +
  theme(
    strip.placement = "outside", #Place facet labels outside x axis labels.
    strip.background = element_blank(), #Make facet label background white.
    axis.title = element_blank(),
    text = element_text(size = 10)) +
  scale_fill_manual(values = chrcolor, name = "ChrState")



data <-Stacked_AE

data2<-data[order(data$NBins,decreasing = TRUE),]
data2$SAMPLE_NAME<-factor(data2$SAMPLE_NAME, levels = data2$SAMPLE_NAME)

# Stacked
ggplot(data2, aes(fill=CELL_TYPE_Abb,
                 y=NBins, x=SAMPLE_NAME)) + 
  geom_bar(position="stack", stat="identity") + theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

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



