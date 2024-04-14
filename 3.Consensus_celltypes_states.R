# Author: Juliana Arcila Galvis
# Description: This script performs analysis on chromatin states and their conservation across different cell types. 
# The analysis includes identifying conserved enhancers, finding differences between pairs of cell types, 
# and visualizing these differences using Sankey plots and bar plots.

#--------------------------------------------------------------------------------
# INPUTS
#--------------------------------------------------------------------------------
# 1. Chromatin states merged data
# 2. Metadata
# 3. Celltypes ontology
# 4. Chromatin states labels and identifiers

setwd("/mnt/nocode/juliana/PID/ThemeII")

# Import chromatin states merged
All_data <- readRDS("INPUTS/states_220422_BACKGROUND_merged.rds")

# Import metadata
metadata <- read.table("OUTPUT/METADATA_paper.txt", sep="\t")

# Import celltypes ontology
celltypes <- read.delim("INPUTS/celltype_data2.txt", sep="\t")

# Import chromatin states labels and identifiers
chrlabels <- read.table("INPUTS/labelsChr.txt", header = TRUE, sep="\t")

#--------------------------------------------------------------------------------
# OUTPUTS
#--------------------------------------------------------------------------------
# 1. Conserved enhancers list: ConservedList.rds
# 2. Conserved enhancers dataframe:ConservedDF.rds
# 3. Differences between pairs of cell types saved as RDS files
# 4. Alluvial plots and bar plots

#--------------------------------------------------------------------------------
# Data processing
#--------------------------------------------------------------------------------

library(dplyr)

# Filter out unique cell types and samples
mycelltypes <- metadata %>% 
  dplyr::count(., CELL_TYPE) %>% 
  filter(., n > 1) %>% 
  data.frame() %>% 
  dplyr::select(CELL_TYPE)

mysamples <- metadata %>% 
  filter(., CELL_TYPE %in% mycelltypes$CELL_TYPE) %>% 
  dplyr::select(CELL_TYPE, SAMPLE_NAME)

Datavenous <- All_data[, colnames(All_data) %in% mysamples$SAMPLE_NAME]

# Define a function to find conserved enhancers
Sum <- function(x) {
  sum(x) >= round(length(x) * 0.75)
}

FindConserved <- function(chromstate) {
  celltype.chromstate <- celltypes == chromstate
  
  ifelse(as.numeric(dim(data.frame(celltypes))[2]) > 1,
         Suma <- apply(celltype.chromstate, 1, Sum),
         Suma <- TRUE)
  
  celltype <- data.frame(celltypes, Suma)
  celltype$conserved <- 0
  celltype[celltype$Suma == TRUE, "conserved"] <-  chromstate
  
  conserved <- data.frame(Conserved = celltype$conserved,
                          row.names = rownames(celltype.chromstate))
  colnames(conserved) <- chrlabels$NameState[chromstate]
  return(conserved)
}

# Identify conserved enhancers by cell type
library(parallel)

Conserved = list()

for (cell in unique(mycelltypes$CELL_TYPE)) {
  celltype <- Datavenous[, colnames(Datavenous) %in%
                           metadata[metadata$CELL_TYPE == cell, "SAMPLE_NAME"]]
  Conservedlist <- mclapply(chrlabels$SymbState, FindConserved, mc.cores = 35)
  Conserved[[as.character(cell)]] <- do.call(cbind, Conservedlist)
}

# Save conserved enhancers list
saveRDS(Conserved, "OUTPUT/ConservedList.rds")

#--------------------------------------------------------------------------------
# CONSERVED ENHANCERS DATAFRAME
#--------------------------------------------------------------------------------

Conserved <- readRDS("ConservedList.rds")
ConservedDF <- mclapply(Conserved, rowSums, mc.cores = 35 )
ConservedDF2 <- data.frame(do.call(cbind, ConservedDF))
saveRDS(ConservedDF2, "OUTPUT/ConservedDF.rds")

#--------------------------------------------------------------------------------
# DYNAMICS OF ENHANCERS BETWEEN PAIRS 
#--------------------------------------------------------------------------------

# Define a function to find differences between two enhancer lists
Differences <- function(One, Two, name_ONETWO_all_shared, name_ONETWO_all_different) {
  One <- One[One != 0]
  Two <- Two[Two != 0]
  
  a <- One[names(One) %in% names(Two)]
  b <- Two[names(Two) %in% names(One)]
  
  d <- cbind(One_states = a, Two_states = b)
  d <- data.frame(d)
  
  ONETWO_all_shared <- d
  saveRDS(ONETWO_all_shared, name_ONETWO_all_shared)
  
  e <- subset(d, d$One_states != d$Two_states)
  saveRDS(e, name_ONETWO_all_different)
}

# Find differences between pairs of cell types
Monocites <- ConservedDF[[names(ConservedDF)[grepl("monocyte", names(ConservedDF))]]]
Osteoclast <- ConservedDF[[names(ConservedDF)[grepl("osteoclast", names(ConservedDF))]]]
Differences(Monocites, Osteoclast, "OUTPUT/Mono.Osteo.shared.rds", "OUTPUT/Mono.Osteo.different.rds")

#--------------------------------------------------------------------------------
# Alluvial plot
#--------------------------------------------------------------------------------

# Import data
Mono_Osteoclast <- readRDS("OUTPUT/Mono.Osteo.different.rds")
colnames(Mono_Osteoclast) <- c("Monocyte","Osteoclast")

# Chromatin states annotations
label_states <- chrlabels

# Import libraries
library(networkD3)
library(dplyr)

# Define Sankey plot function
Sankey_plot <- function(A_B, label_states, transition_name) {
  A_B <- A_B[apply(A_B, 1, function(r) any(r == 0) == FALSE),]
  
  Changes_A_B <- plyr::count(A_B, vars = colnames(A_B))
  colnames(Changes_A_B)[1] <- "SymbState"
  
  Changes_A_B <- left_join(Changes_A_B, label_states)
  
  colnames(Changes_A_B)[c(1, 2, 4)] <- c("state1", "SymbState", "source")
  Changes_A_B <- left_join(Changes_A_B, label_states)
  colnames(Changes_A_B)[c(2, 3, 5)] <- c("state2", "value", "target")
  
  Changes_A_B <- Changes_A_B[, c("source", "target", "value")]
  Changes_A_B$source <- paste("A_", Changes_A_B$source, sep = "")
  Changes_A_B$target <- paste("B_", Changes_A_B$target, sep = "")
  
  Changes_A_B$source <- factor(Changes_A_B$source, levels = paste("A_", rev(level_order), sep = ""))
  Changes_A_B$target <- factor(Changes_A_B$target, levels = paste("B_", rev(level_order), sep = ""))
  
  micolorA <- data.frame(color = chrcolor, SymbState = levels(Changes_A_B$source))
  micolorB <- data.frame(color = chrcolor, SymbState = levels(Changes_A_B$target))
  
  micolor <- rbind(micolorA, micolorB)
  
  my_color <- 'd3.scaleOrdinal(["#D3D3D3","#1B9E77","#d95f02","#A6761D","#E6AB02","#66A61E","#E7298A", "#7570B3","#D95F02","#1B9E77", "#A6761D","#E6AB02","#66A61E","#E7298A", "#7570B3", "#D95F02"])'
  
  links <- data.frame(Changes_A_B[, 1:2], value = as.numeric(Changes_A_B[, 3]))
  nodes <- data.frame(name = c(paste("A_", level_order, sep = ""), paste("B_", level_order, sep = "")))
  
  level_nodes <- c(rev(levels(Changes_A_B$source)), rev(levels(Changes_A_B$target)))
  nodes$name <- factor(nodes$name, levels = level_nodes)
  
  labels <- tidyr::separate(nodes, col = name, into = c("celltype", "NameState"))
  nodes$label <- factor(labels$NameState, levels = level_order)
  nodes$node.colr <- factor(as.character(nodes$label), levels = level_order)
  
  level_ID <- data.frame(Names = nodes$name, IDS = 0:13)
  
  IDsource <- left_join(data.frame(Names = links$source), level_ID)
  IDtarget <- left_join(data.frame(Names = links$target), level_ID)
  
  links$IDsource <- IDsource$IDS
  links$IDtarget <- IDtarget$IDS
  
  linkis <- links
  linkis$link.colr <- "othercol"
  
  linkis[linkis$target == "B_PoisE",]$link.colr <- "B_PoisEcol"
  linkis[linkis$target == "B_AconE",]$link.colr <- "B_AconEcol"
  
  linkis$link.colr <- factor(linkis$link.colr, levels = c("othercol", "B_PoisEcol","B_AconEcol"))
  
  sankeyNetwork(
    Links = linkis,
    Nodes = nodes,
    Source = "IDsource",
    Target = "IDtarget",
    Value = "value",
    NodeID = "label",
    sinksRight = TRUE,
    fontSize = 20,
    nodeWidth = 50,
    nodePadding = 30,
    units = 'TWh',
    colourScale = my_color,
    LinkGroup = "link.colr",
    NodeGroup = "node.colr"
  )
}

# Generate Sankey plot for Mono_M0
Mono_M0_sankey <- Sankey_plot(Mono_M0, label_states, "MonoM0")
