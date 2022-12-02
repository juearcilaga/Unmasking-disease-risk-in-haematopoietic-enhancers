setwd("/mnt/nocode/juliana/PID/ThemeII")

#Import chromatin states merged
All_data<-readRDS("INPUTS/states_220422_BACKGROUND_merged.rds")
dim(All_data)#15441337      108

#Import metadata
metadata<-read.table("INPUTS/METADATA_paper2.txt", sep="\t", header = TRUE)
library(dplyr)

mycelltypes <- metadata %>%
  filter(., TISSUE_TYPE == "venous blood") %>%
  count(., "CELL_TYPE") %>%
  filter(., freq > 1) %>% data.frame() %>% dplyr::select("CELL_TYPE")

mycelltypes <-metadata %>% filter(., CELL_TYPE  %in% mycelltypes$CELL_TYPE) 

Datavenous<-All_data[, colnames(All_data) %in% mycelltypes$SAMPLE_NAME]
dim(Datavenous)#[1] 15441337       73

unique(Datavenous[,1])#1-7

# #Import chromatin states labels and identifiers
chrlabels<-read.table("INPUTS/labelsChr.txt", header = TRUE,sep="\t")

# 75% of the samples should  be consistent
Sum <- function(x) {
  sum(x) >= round(length(x) * 0.75)
}

#This only works for one chromstate  whe need it for all chromstates

FindConserved <- function(chromstate) {
  # chromstate <- chrlabels$SymbState[2]
  # chromstate <-2
  celltype.chromstate <- celltype == chromstate
  
  ifelse(as.numeric(dim(data.frame(celltype))[2]) > 1,
         Suma <- apply(celltype.chromstate, 1, Sum),
         Suma <- TRUE)
  
  celltype <- data.frame(celltype, Suma)
  celltype$conserved <- 0
  celltype[celltype$Suma == TRUE, "conserved"] <-  chromstate
  
  print(table(celltype$conserved))
  conserved <-
    data.frame(Conserved = celltype$conserved,
               row.names = rownames(celltype.chromstate))
  colnames(conserved) <- chrlabels$NameState[chromstate]
  return(conserved)
}

# library(tidyr)
library(parallel)
library(MASS)

Conserved = list()
for (cell in unique(mycelltypes$CELL_TYPE)) {
  # cell<-mycelltypes$CELL_TYPE[13]
  print(cell)
  celltype <-
    Datavenous[, colnames(Datavenous) %in%
                 metadata[metadata$CELL_TYPE == cell, "SAMPLE_NAME"]]
  
  Conservedlist<-mclapply(chrlabels$SymbState, FindConserved, mc.cores = 35 )
  
  Conserved[[as.character(cell)]] <- do.call(cbind, Conservedlist)
  }

str(Conserved)

# [1] "mature eosinophil"
# [2] "mature conventional dendritic cell"
# [3] "immature conventional dendritic cell"
# [4] "osteoclast"
# [5] "CD14-positive, CD16-negative classical monocyte"
# [6] "macrophage"
# [7] "alternatively activated macrophage"
# [8] "inflammatory macrophage"
# [9] "mature neutrophil"
# [10] "CD38-negative naive B cell"
# [11] "class switched memory B cell"
# [12] "naive B cell"
# [13] "CD4-positive, alpha-beta T cell"
# [14] "effector memory CD8-positive, alpha-beta T cell"
# [15] "mesenchymal stem cell of the bone marrow"
# [16] "adult endothelial progenitor cell"


#Save this object Conserved
saveRDS(Conserved, "OUTPUT/ConservedList.rds")
# ConservedDF <- do.call(rbind, Conserved)
# ConservedDF <- do.call(cbind,rowSums(ConservedDF))

Conserved<-readRDS("ConservedList.rds")
ConservedDF<-mclapply(Conserved, rowSums, mc.cores = 35 )
ConservedDF2<-data.frame(do.call(cbind, ConservedDF))
saveRDS(ConservedDF2, "OUTPUT/ConservedDF.rds")

ConservedDF<-ConservedDF2

Monocites<-ConservedDF[[names(ConservedDF)[grepl("monocyte", names(ConservedDF))]]]
names(Monocites)<-rownames(ConservedDF)
table(Monocites)
Osteoclast<-ConservedDF[[names(ConservedDF)[grepl("osteoclast", names(ConservedDF))]]]
names(Osteoclast)<-rownames(ConservedDF)

Differences <-
  function(One,
           Two,
           name_ONETWO_all_shared,
           name_ONETWO_all_different) {
    # One<-Monocites[Monocites!=0]
    # Two<-Osteoclast[Osteoclast!=0]
    # name_ONETWO_all_shared<-"Mono.Osteo.shared.rds"
    # name_ONETWO_all_different<-"Mono.Osteo.different.rds"
    
    One <- One[One != 0]
    Two <- Two[Two != 0]
    
    a <- One[names(One) %in% names(Two)]
    print(length(a))
    
    b <- Two[names(Two) %in% names(One)]
    print(length(b))
    
    d <- cbind(One_states = a, Two_states = b)
    d <- data.frame(d)
    
    ONETWO_all_shared <- d
    saveRDS(ONETWO_all_shared, name_ONETWO_all_shared)
    
    e <- subset(d, d$One_states != d$Two_states)
    print(dim(e))
    saveRDS(e, name_ONETWO_all_different)
  }


Differences(Monocites,
            Osteoclast,
            "OUTPUT/Mono.Osteo.shared.rds",
            "OUTPUT/Mono.Osteo.different.rds")

Immature.DC<-ConservedDF[[names(ConservedDF)[grepl("immature", names(ConservedDF))]]]
names(Immature.DC)<-rownames(ConservedDF)

Differences(Monocites,
            Immature.DC,
            "OUTPUT/Mono.Immature.DC.shared.rds",
            "OUTPUT/Mono.Immature.DC.different.rds")

Macrophage.M0<-ConservedDF[,names(ConservedDF)=="macrophage"]
names(Macrophage.M0)<-rownames(ConservedDF)

Differences(Monocites,
            Macrophage.M0,
            "OUTPUT/Mono.Macrophage.M0.shared.rds",
            "OUTPUT/Mono.Macrophage.M0.different.rds")

# #######################################################################################################
#Making an alluvial plot

#Import data

Mono_ImDendri<-readRDS("OUTPUT/Mono.Immature.DC.different.rds")
colnames(Mono_ImDendri)<-c("Monocyte","ImDendritic")

Cuenta<-dim(Mono_ImDendri[apply(Mono_ImDendri, 1, function(r)
  any(r==0) == FALSE), ])#1057594       2

Mono_M0<-readRDS("OUTPUT/Mono.Macrophage.M0.different.rds")
colnames(Mono_M0)<-c("Monocyte","M0")

dim(Mono_M0[apply(Mono_M0, 1, function(r)
  any(r==0) == FALSE), ])#1057594       2
Mono_Osteoclast<-readRDS("OUTPUT/Mono.Osteo.different.rds")
colnames(Mono_Osteoclast)<-c("Monocyte","Osteoclast")

dim(Mono_Osteoclast[apply(Mono_Osteoclast, 1, function(r)
  any(r==0) == FALSE), ])#231679      2

#Chromatin states annotations
label_states <-chrlabels

library(RColorBrewer)

chrcolor<-brewer.pal(n = 7, name = 'Dark2')

level_order <-
  c("TransElo",
    "HetChrom",
    "PolyComb",
    "AProm",
    "AeloE",
    "AconE",
    "PoisE")

names(chrcolor)<-rev(level_order)

#Import libraries
library(networkD3)
library(dplyr)


Sankey_plot <- function(A_B, label_states, transition_name) {
  #Counts dataframe
  # A_B<-Mono_M0
  # A_B<-Mono_ImDendri
  # A_B<-Mono_Osteoclast
  A_B <- A_B[apply(A_B, 1, function(r)
    any(r == 0) == FALSE),]
  transition_name <- "MonoM0"
  
  Changes_A_B <-
    plyr::count(A_B, vars = colnames(A_B))
  colnames(Changes_A_B)[1] <- "SymbState"
  Changes_A_B <- left_join(Changes_A_B, label_states)
  colnames(Changes_A_B)[c(1, 2, 4)] <-
    c("state1", "SymbState", "source")
  Changes_A_B <- left_join(Changes_A_B, label_states)
  colnames(Changes_A_B)[c(2, 3, 5)] <-
    c("state2", "value", "target")
  
  #file_name<-paste("ChangesFreq/",transition_name, sep="")
  Changes_A_B <- Changes_A_B[, c("source", "target", "value")]
  
  Changes_A_B$source <- paste("A_", Changes_A_B$source, sep = "")
  Changes_A_B$target <- paste("B_", Changes_A_B$target, sep = "")
  
  Changes_A_B$source <-
    factor(Changes_A_B$source, levels = paste("A_", rev(level_order), sep = ""))
  Changes_A_B$target <-
    factor(Changes_A_B$target, levels = paste("B_", rev(level_order), sep = ""))
  
  #DF<-Changes_A_B
  #saveRDS(DF, paste(file_name, "ChangesFreq.rds", sep="_"))
  
  #str(Changes_A_B)
  
  micolorA <-
    data.frame(color = chrcolor,
               SymbState = levels(Changes_A_B$source))
  micolorB <-
    data.frame(color = chrcolor,
               SymbState = levels(Changes_A_B$target))
  micolor <- rbind(micolorA, micolorB)
  
  my_color <-
    'd3.scaleOrdinal(["#D3D3D3","#d95f02","#1B9E77","#A6761D","#E6AB02","#66A61E",
  "#E7298A", "#7570B3","#D95F02","#1B9E77", "#A6761D",
  "#E6AB02","#66A61E","#E7298A", "#7570B3", "#D95F02"])'
  
  # A connection data frame is a list of flows with intensity for each flow
  links <-
    data.frame(Changes_A_B[, 1:2], value = as.numeric(Changes_A_B[, 3]))
  # From these flows we need to create a node data frame: it lists every entities involved in the flow
  
  nodes <-
    data.frame(name = c(
      paste("A_", level_order, sep = ""),
      paste("B_", level_order, sep = "")
    ))
  
  level_nodes <-
    c(rev(levels(Changes_A_B$source)), rev(levels(Changes_A_B$target)))
  nodes$name <-
    factor(nodes$name, levels = level_nodes)
  
  labels <-
    tidyr::separate(nodes,
                    col = name,
                    into = c("celltype", "NameState"))
  # labels <- dplyr::left_join(labels, label_states, by = "NameState")
  nodes$label <- factor(labels$NameState, levels = level_order)##
  
  nodes$node.colr <-
    factor(as.character(nodes$label), levels = level_order)
  
  level_ID <- data.frame(Names = nodes$name, IDS = 0:13)
  
  # With networkD3, connection must be provided using id, not using real name like in the links dataframe.. So we need to reformat it.
  IDsource <- left_join(data.frame(Names = links$source), level_ID)
  IDtarget <- left_join(data.frame(Names = links$target), level_ID)
  
  links$IDsource <- IDsource$IDS
  links$IDtarget <- IDtarget$IDS
  
  linkis <- links
  linkis$link.colr <- "othercol"
  linkis[linkis$source == "A_PoisE",]$link.colr <- "A_PoisEcol"
  linkis[linkis$target == "B_AconE",]$link.colr <- "B_AconEcol"
  linkis$link.colr <-
    factor(linkis$link.colr,
           levels = c("othercol", "A_PoisEcol", "B_AconEcol"))
  # Make the Network
  p <- sankeyNetwork(
    Links = linkis,
    Nodes = nodes,
    Source = "IDsource",
    Target = "IDtarget",
    Value = "value",
    NodeID = "label",
    sinksRight = TRUE,
    fontSize = 20,
    nodeWidth = 50,
    nodePadding = 20,
    units = 'TWh',
    colourScale = my_color,
    LinkGroup = "link.colr",
    NodeGroup = "node.colr"
  )
  return(p)
  #return(DF)
}



Mono_M0_sankey<-Sankey_plot(Mono_M0, label_states, "MonoM0")

#################################################################################################
ConservedDF <- readRDS("/mnt/nocode/juliana/PID/ThemeII/OUTPUT/ConservedDF.rds")

#Import metadata
metadata<-read.table("OUTPUT/METADATA_paper.txt", sep="\t")

#Import celltypes ontology
celltypes<-read.delim("INPUTS/celltype_data.txt", sep="\t")

#Import chromatin states labels and identifiers
chrlabels<-read.table("INPUTS/labelsChr.txt", header = TRUE,sep="\t")

mylineage<-ConservedDF 
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

StackObj.long$SAMPLE_NAME<-gsub("\\."," ", StackObj.long$SAMPLE_NAME)
StackObj.long$SAMPLE_NAME<-gsub("  "," ", StackObj.long$SAMPLE_NAME)
StackObj.long <- left_join(StackObj.long, chrlabels)
StackObj.long$SymbState <- as.integer(StackObj.long$SymbState)

celltypes$CELL_TYPE2<-gsub("\\-"," ", celltypes$CELL_TYPE)
celltypes$CELL_TYPE2<-gsub("\\,"," ", celltypes$CELL_TYPE2)
celltypes$CELL_TYPE2<-gsub("  "," ", celltypes$CELL_TYPE2)
colnames(StackObj.long)[1]<-"CELL_TYPE2"

StackObj.long<- left_join(StackObj.long,
                          celltypes[, c("MAIN_CELL_TYPE", "CELL_TYPE_Abb", "CELL_TYPE2")])
StackObj.long$MAIN_CELL_TYPE <-
  factor(StackObj.long$MAIN_CELL_TYPE, level = rev(unique(celltypes$MAIN_CELL_TYPE)))

StackObj.long$SIZE_bp<-as.numeric(StackObj.long$NBins*200)

#write.table(StackObj.long, "OUTPUT/Enhancerbund.txt", sep="\t")
write.table(StackObj.long, "OUTPUT/Enhconservbund_all.txt", sep="\t")

#----------------------------------------------------------------------------
data <-StackObj.long[StackObj.long$NameState=="AconE",]
data2<-data[order(data$NBins,decreasing = TRUE),]
data2$CELL_TYPE_Abb<-factor(data2$CELL_TYPE_Abb, levels = data2$CELL_TYPE_Abb)
data3<-data2
data3[data3$NameState=="AeloE","NameState"]<-"AconE"
group(NBins, by=c(CELL_TYPE2,NameState))

data3<-data3 %>% group_by(CELL_TYPE2,NameState)  %>%
  summarise(total_sales = sum(NBins))

ggplot(data3, aes(y=total_sales, x=CELL_TYPE2)) + 
  geom_bar(position="stack", stat="identity") + theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

# Stacked
ggplot(data2, aes(y=NBins, x=CELL_TYPE_Abb)) + 
  geom_bar(position="stack", stat="identity") + theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 