#-----------------------------------------------------------------
# DESCRIPTION:
# This script identifies  and analyzes enhancer activity patterns across different
# cell types. It also visualizes them through heatmaps
# and bar plots.
#
# AUTHOR:
# Juliana Arcila Galvis
#
# INPUTS:
# - ConservedDF.rds: Dataframe containing conserved enhancer regions.
# - BActiveE.rds: Dataframe containing active enhancer IDs.
#
# OUTPUTS:
# - Table_PatternsAE.txt: Table containing enhancer activity patterns.
# - gr_AE_list.rds: GenomicRanges object of enhancer regions.
# - AE_abundance_celltype.txt: Table containing enhancer abundance per cell type.
# - cellclust.txt: File containing clustered cell types based on enhancer activity.
# - patternclust.txt: File containing clustered enhancer activity patterns.
#
#------------------------------------------------------------------
#               IMPORTING CONSENSUS ENHANCERS 
#------------------------------------------------------------------
setwd("/mnt/nocode/juliana/PID/ThemeII")

# Importing conserved enhancers
conserved <- readRDS("OUTPUT/ConservedDF.rds")

# Importing IDs of active enhancers in at least 2 samples
AllAE <- readRDS("OUTPUT/BActiveE.rds")

# Selecting only active enhancers
All_dataAE <- conserved[rownames(conserved) %in% AllAE,]

# Modifying data to label enhancers as active (1) or inactive (0)
data_modif <- as.matrix(All_dataAE)
data_modif[data_modif %in% c(1,2,3,7,5)] <- 0  # Inactive enhancers
data_modif[data_modif %in% c(4,6)] <- 1  # Active enhancers
colnames(data_modif) <- gsub("\\.", "", colnames(data_modif))

#------------------------------------------------------------------
# Counting enhancer activity patterns
#------------------------------------------------------------------
library(plyr)
Patterns <- plyr::count(df = data.frame(data_modif), vars = colnames(data_modif))
Patterns <- Patterns[order(Patterns$freq, decreasing = TRUE),]

# Visualizing distribution of enhancer activity patterns
length(Patterns$freq)  # Total number of unique enhancer patterns
summary(Patterns[2:length(Patterns$freq), "freq"])

# Creating a boxplot to visualize enhancer activity patterns
library(ggplot2)
ggplot(Patterns[2:length(Patterns$freq),], aes(x = "Enhancer Activity Patterns", y = freq)) +
  geom_boxplot(outlier.size = 0.5) +
  theme_classic() + 
  scale_y_continuous(limits = c(1, 80000)) +
  labs(y = "Number of Active Enhancer bins") +
  theme(axis.title.x = element_blank())

#------------------------------------------------------------------
# Identifying enhancer activity patterns
#------------------------------------------------------------------
pattern <- function(x){
  x <- which(x == 1)
  return(x)
}

patterns <- list()
for (i in 1:length(Patterns$freq)){
  patterns <- apply(Patterns[i,], 1, pattern)
  patterns <- paste(colnames(Patterns[i,])[patterns], collapse = "|")
  patterns[[i]] <- patterns
}

Patterns$pattern <- do.call(rbind, patterns)
write.table(Patterns, "OUTPUT/Table_PatternsAE.txt", sep = "\t")

#------------------------------------------------------------------
# Creating GenomicRanges object and saving
#------------------------------------------------------------------
Mylist <- data.frame(data_modif)
Milista <- as.matrix(Mylist)
colnames(Milista) <- c("endprog", "M2", "bfneut", "mono", "megkaryo", "CD38negnaiveB",
                       "CD4T", "CD8T", "cswitmB", "NK", "effmemCD8T", "endumbprolif",
                       "eryth", "germcenterB", "imDC", "M1", "M0", "mDC", "meos", 
                       "mneut", "mesenBM", "naiveB", "neutmetmyelo", "neutmyelo", "osteo",
                       "plasma", "segneutBM")

pattern2 <- function(y) {
  patterns <- pattern(y)
  patterns <- colnames(Milista)[patterns]
  patterns <- paste(patterns, collapse = "|")
  return(patterns)
}

patterns <- apply(Milista, 1, pattern2)
Mylist2 <- cbind(Milista, patterns)
Mylist <- data.frame(Mylist2)
head(Mylist)

# Creating GRanges object
CtypeChrStates <- Mylist
CtypeChrStates$regionID <- rownames(CtypeChrStates)
CtypeChrStates$coord <- rownames(CtypeChrStates)

library(GenomicRanges)
CtypeChrStates <- tidyr::separate(CtypeChrStates, coord, c("chr", "region"))
CtypeChrStates$end <- as.numeric(CtypeChrStates$region) * 200
CtypeChrStates$start <- CtypeChrStates$end - 200

gr_CtypeChrStates <- makeGRangesFromDataFrame(
  CtypeChrStates,
  keep.extra.columns = TRUE,
  ignore.strand = TRUE,
  seqinfo = NULL,
  seqnames.field = c("seqnames", "seqname", "chromosome", "chrom", "chr",
                     "chromosome_name", "seqid"),
  start.field = "start",
  end.field = "end",
  strand.field = NULL,
  starts.in.df.are.0based = TRUE
)

saveRDS(gr_CtypeChrStates, "OUTPUT/gr_AE_list.rds")

#------------------------------------------------------------------
# Visualization of enhancer activity patterns
#------------------------------------------------------------------
Patterns <- read.table("OUTPUT/Table_PatternsAE.txt", sep = "\t")

# Visualizing the heatmap of enhancer activity patterns
library(superheat)
library(RColorBrewer)

data <- Patterns[100:2, 1:27]
rownames(data) <- Patterns$Pattern[100:2]

# Identifying specific and non-specific enhancer patterns
nonspec <- grepl("\\|", rownames(data))
spec <- grepl("\\|", rownames(data)) == FALSE

logic <- data == 1 & nonspec
data[logic] <- 2

# Creating the heatmap
superheat(data, scale = F, pretty.order.rows = FALSE, pretty.order.cols = FALSE,
          grid.vline = TRUE, smooth.heat = FALSE, smoothing.method = "lm",
          force.grid.hline = TRUE, heat.pal = c("light grey", "#8c1b2c", "#293266"),
          heat.pal.values = c(0, 2, 1), left.label.col = "white",
          left.label.text.size = 3, bottom.label.text.angle = 90,
          bottom.label.text.size = 3, yr.bar.col = "light grey",
          yr.obs.col = ifelse(rownames(data) %in% rownames(data)[(nonspec)],
                             "#293266", "#8c1b2c"),
          yr = Patterns$freq[100:2], yr.axis.name = "Num. AE bins",
          yr.plot.type = "bar", yr.plot.size = 1, legend.width = 0.5,
          legend.text.size = 10, bottom.label.col = "white")

# Biplot clustering rows and columns of the patterns
library(pheatmap)
data_spec_non <- rbind(data_spec, data_nonspec)

paleta <- ifelse(grepl("\\|", rownames(data_spec_non)),
                 palleta <- c("white", "#542788"), 
                 palleta <- c("white", "red"))

anotcols <- data.frame(Freq = data2[rownames(data_spec_non), 28],
                       row.names = data2[rownames(data_spec_non), "Pattern"])

anotcols$Pattern_Frequency <- "Between min & 1st Qu."
anotcols[anotcols$Freq >= summary(anotcols$Freq)["1st Qu."], "Pattern_Frequency"] <- "Between 1st Qu. & Median"
anotcols[anotcols$Freq >= summary(anotcols$Freq)["Median"], "Pattern_Frequency"] <- "Between Median & 3rd Qu."
anotcols[anotcols$Freq >= summary(anotcols$Freq)["3rd Qu."], "Pattern_Frequency"] <- "Between 3rd Qu. & Max."

anotcols$Pattern_Frequency <- factor(anotcols$Pattern_Frequency,
                                     levels = c("Between 3rd Qu. & Max.",
                                                "Between Median & 3rd Qu.",
                                                "Between 1st Qu. & Median",
                                                "Between min & 1st Qu."))

levelsColors <- c("#581845", "#900C3F", "#C70039", "#FF5733")

ann_colors <- list(Pattern_Frequency = levelsColors)

pheatmap(t(data_spec_non), color = c("white", "#293266", "#8c1b2c"),
         cellwidth = 6, cellheight = 8, border_color = "#3b2942",
         scale = "none", cluster_rows = TRUE, cluster_cols = TRUE,
         treeheight_col = 30, treeheight_row = 25, show_colnames = FALSE,
         fontsize_row = 10, annotation_col = anotcols, annotation_colors = ann_colors)