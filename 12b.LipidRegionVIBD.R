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
GENES.DF<-read.table("OUTPUT/IBD.GENES.DF.txt", sep="\t")
GENES.summary<-GENES.DF %>% group_by(seqnames)%>% 
  summarise(genes=paste(unique(gene_name), collapse = ", "))
ggtexttable(GENES.summary, rows = NULL)
# cromosome 2
# CXCR1: https://www.proteinatlas.org/ENSG00000163464-CXCR1/immune+cell
#file:///H:/Downloads/ijms-23-05076.pdf, https://pubmed.ncbi.nlm.nih.gov/34893315/, https://www.sciencedirect.com/science/article/pii/S1074761320300893
# https://www.proteinatlas.org/ENSG00000163464-CXCR1/single+cell+type
# ARPC2 https://www.proteinatlas.org/ENSG00000163466-ARPC2/immune+cell
# in NET https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8488397/pdf/fimmu-12-675315.pdf
#Enhancer coordinates (enhancers snps regions)
SNPsToEnh.summary.bychr<-read.table("OUTPUT/SNPsToEnhsigIBD.summary.bychr.txt", sep="\t")
colnames(SNPsToEnh.summary.bychr)<-c("chr", "start20k", "end20k")

#Snps coordinates
SNPsToEnh.summary<-read.table("OUTPUT/SNPsToEnhsigIBDsummary.txt", sep="\t")

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
#---------------------------TABLES---------------------------------------------
GENES.DF$seqnames <- factor(GENES.DF$seqnames,
                            levels = c("chr2", "chr6", "chr8", "chr11", "chr19", "chr22", "chrX"))

GENES.DF <- with(GENES.DF, GENES.DF[order(seqnames), ])
dev.off()
ggtexttable(GENES.DF[, c(1, 6:11)], rows = NULL)


SNPsToEnh.summary$chr <- factor(
  SNPsToEnh.summary$chr,
  levels = c("chr2", "chr6", "chr8", "chr11", "chr19", "chr22", "chrX")
)
SNPsToEnh.summary <-
  with(SNPsToEnh.summary, SNPsToEnh.summary[order(chr), ])

SNPsToEnh.summary$Enhancer_Size_kb <- SNPsToEnh.summary$NEnh * 200 / 1000
SNPsToEnh.summary <-
  SNPsToEnh.summary %>% arrange(DISEASE.TRAIT,-size)

ggtexttable(SNPsToEnh.summary[, c(1, 3, 4, 9, 2, 8)], rows = NULL)





#-----------------------------------------------------------------------------------------
# Chr 2 Neutrophils
# Chr 7 Eryth y ends
# chr9  chr7  chr6  chr5  chr4  chr21 chr2 
#-----------------------------------------------------------------------------------------
SetofEnh <- 2# chr7
chromosome <- as.character(unique(SNPsToEnh.summary.bychr$chr)[SetofEnh])
start <- SNPsToEnh.summary.bychr$start20k[SetofEnh]
end <- SNPsToEnh.summary.bychr$end20k[SetofEnh]

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
  ifelse(
    as.numeric(dim(data.frame(celltype))[2]) > 1,
    Suma <- apply(celltype.chromstate, 1, Sum),
    Suma <- celltype.chromstate
  )
  
  celltype2 <- data.frame(celltype, Suma)
  celltype2$conserved <- 0
  celltype2[celltype2$Suma == TRUE, "conserved"] <-  chromstate
  print(table(celltype2$conserved))
  conserved <-
    data.frame(
      Conserved = as.numeric(celltype2$conserved),
      row.names = rownames(celltype.chromstate)
    )
  colnames(conserved) <-
    chrlabels$NameState[as.numeric(chromstate)]
  return(conserved)
}



Conserved = list()
for (cell in unique(metadata$CELL_TYPE)) {
  print(cell)
  celltype <-
    RegionStates[, colnames(RegionStates) %in%
                   metadata[metadata$CELL_TYPE == cell, "SAMPLE_NAME"]]
  
  Conservedlist <-
    mclapply(chrlabels$SymbState, FindConserved, mc.cores = 35)
  Conserved[[as.character(cell)]] <- do.call(cbind, Conservedlist)
}
ConservedDF <- mclapply(Conserved, rowSums, mc.cores = 35)
ConservedDF2 <- data.frame(do.call(cbind, ConservedDF))

# filename <- sprintf("OUTPUT/ConservedList_Lipids_%s.rds", chromosome)
# saveRDS(ConservedDF2, filename)

#-----------------------------------------------------------------------------------------
# Formating data for the plot
#-----------------------------------------------------------------------------------------
# filename <- sprintf("OUTPUT/ConservedList_Lipids_%s.rds", chromosome)
# ConservedDF2 <- readRDS(filename)

##Formating Y labels
celltypes_Abb <- read.delim("OUTPUT/celltype_Abb.txt")
celltypes <- read.delim("INPUTS/celltype_data2.txt")
celltypes_Abb$variable <- FixLabels(celltypes_Abb$variable)

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

ConservedDF.longer$CELL_TYPE <-
  factor(ConservedDF.longer$CELL_TYPE,
         levels = levels(metadata$CELL_TYPE))

#levels(ConservedDF.longer$CELL_TYPE)
ConservedDF.longer <-
  transform(ConservedDF.longer, ID = as.numeric(factor(CELL_TYPE)))

###Time expensive step
ConservedDF.longer <-
  tidyr::separate(ConservedDF.longer,
                  col = "coord",
                  into = c("chr", "region"))

#-----------------------------------------------------------------------------------------
# Formating coordintaes for the segment plot
#-----------------------------------------------------------------------------------------
mydata = data.frame(
  variable = ConservedDF.longer$CELL_TYPE,
  pos = ConservedDF.longer$ID,
  start = (as.numeric(ConservedDF.longer$region) * 200) - 199,
  end = as.numeric(ConservedDF.longer$region) * 200,
  SymbState = ConservedDF.longer$COLOR
)

mydata <- left_join(mydata, celltypes_Abb)
mydata$acronym <-
  factor(mydata$acronym, levels = rev(unique(celltypes$PATTERN_NAME)))

#checking unique(mydata$acronym)
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

level_order <- c(
  "gene",
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

## Seting the snp positions as vertical lines

# SNPsToEnh.summary$TRAITAbb <-
#   c(
#     "LDL",
#     "HDL_TG_metT",
#     "TG",
#     "LipMet",
#     "HDLCTG",
#     "metSyndbivT",
#     "TG_BP",
#     "HDL",
#     "metSyndbivT",
#     "metSynd",
#     "metSyndbivT",
#     "TG_HDL",
#     "LipMet",
#     "WaistHip",
#     "MetLevels",
#     "MetLevels",
#     "HDL_HDLCTG_TG",
#     "resp_statins_LDLchange_MetLevels",
#     "HDL"
#   )
#

# SNPsToEnh.summary<-SNPsToEnh.summary[SNPsToEnh.summary$chr==chromosome,] %>%
#   unite(labely, SNPS,TRAITAbb, sep=".")


vertical.lines <-
  SNPsToEnh.summary[SNPsToEnh.summary$chr == chromosome, "CHR_POS"]

# names(vertical.lines)<-
#   SNPsToEnh.summary[SNPsToEnh.summary$chr==chromosome,"labely"]

names(vertical.lines) <-
  SNPsToEnh.summary[SNPsToEnh.summary$chr == chromosome, "SNPS"]

##Adding labels to the formated data for the plot
mydata$SymbState <- as.character(mydata$SymbState)
mydata <- left_join(mydata, chrlabels)

InRegionGenes <- GENES.DF[GENES.DF$seqnames == chromosome,]%>% data.frame()
InRegionGenes$pos<-1:length(InRegionGenes$seqnames)
InRegionGenes <-
  data.frame(
    variable = InRegionGenes$gene_name,
    pos = max(InRegionGenes$pos) + InRegionGenes$pos,
    start = as.numeric(InRegionGenes$start),
    end = as.numeric(InRegionGenes$end),
    SymbState = "8",
    NameState = "gene",
    acronym = InRegionGenes$gene_name
  )

mydata.genes <- rbind(mydata, InRegionGenes)
#-------------------------------------------------------------------------------------------
##Graph for all the celltypes
#-------------------------------------------------------------------------------------------
Chromosomename <- sprintf("Position in %s", chromosome)

# Size of the plot 1200 * 600
dev.off()

# parameter controling segsize
# parameter controling graphize

dim(mydata.genes)#8941552
colnames(mydata.genes)
unique(mydata.genes$pos)#32
dim(mydata)


bycelltype <- 
  ggplot(mydata.genes,
         aes(x = start,
             xend = end,
             y = acronym,
             yend = acronym,
             color = NameState), 
         size = 3) +
  geom_segment(size = 3.5) +
  scale_color_manual(values = chrcolor) +
  theme_bw() +
  xlab(Chromosomename) +
  # ggtitle("Conserved chromatin states by celltype") +
  ylab("cell type")+
  theme(text = element_text(size = 12),
        axis.text = element_text(color = "black"))+
  geom_vline(xintercept = vertical.lines, linetype="dashed", 
             color = "dark red", size=0.2)+
  annotate(geom = "text",
           label = names(vertical.lines),
           x = vertical.lines,
           y = 33,
           size = 1*.pt,
           angle = 40, 
           hjust = 0,
           colour="dark red"
  ) + # Adding text
  theme(plot.margin = unit(c(4, 1, 1, 1), "lines")) +
  coord_cartesian(clip = "off")



#OPTIONAL

#-----------------------------------------------------------
# by sample This must change acording to the celltype in
# which the enchancers are active
#-----------------------------------------------------------
ConservedDF.longer2 <- RegionStates
ConservedDF.longer2$coord <- rownames(ConservedDF.longer2)
ConservedDF.longer2 <-
  pivot_longer(
    ConservedDF.longer2,
    cols = 1:length(ConservedDF.longer2) - 1,
    names_to = "SAMPLE_NAME",
    values_to = "COLOR"
  )
ConservedDF.longer2 <- left_join(ConservedDF.longer2, metadata)
ConservedDF.longer2 <-
  with(ConservedDF.longer2, ConservedDF.longer2[order(CELL_TYPE),])

ConservedDF.longer2 <-
  transform(ConservedDF.longer2, ID = as.numeric(factor(SAMPLE_NAME)))
ConservedDF.longer2 <-
  separate(ConservedDF.longer2,
           col = "coord",
           into = c("chr", "region"))


mydata3 = data.frame(
  variable = ConservedDF.longer2$SAMPLE_NAME,
  # pos = ConservedDF.longer2$ID,
  start = (as.numeric(ConservedDF.longer2$region) * 200) - 199,
  end = as.numeric(ConservedDF.longer2$region) * 200,
  SymbState = ConservedDF.longer2$COLOR,
  celltype = ConservedDF.longer2$CELL_TYPE
)

mydata3$SymbState <- as.character(mydata3$SymbState)
mydata3 <- left_join(mydata3, chrlabels)

M012 <-
  unique(metadata[grepl("neut", metadata$CELL_TYPE), "SAMPLE_NAME"])

mydata3.genes <- mydata3[mydata3$variable %in% M012, ]


InRegionGenes.sample <- data.frame(InRegionGenes[, c(1, 3:5)],
                                   celltype = "NA",
                                   NameState = InRegionGenes[, 6])

mydata3.genes <- rbind(mydata3.genes, InRegionGenes.sample)

rm(InRegionGenes.sample)

positions <-
  mydata3.genes %>% dplyr::select(variable, celltype) %>% unique()
positions$pos <- 0
positions[positions$celltype != "NA", "pos"] <-
  2:(length(positions[positions$celltype != "NA", "pos"]) + 1)
positions[positions$celltype == "NA", "pos"] <- 1

mydata3.genes <- left_join(mydata3.genes, positions)
# mydata3.genes<-with(mydata3.genes, mydata3.genes[order(celltype),])
mydata3.genes$variable <- factor(mydata3.genes$variable,
                                 levels = unique(mydata3.genes$variable))

mydata3.genes %>% dplyr::select(variable, celltype, pos) %>% unique()

coloraxis <- c(
  rep("#6b2520", 13),
  rep("#3a5431", 3),
  rep("#1f5387", 3),
  rep("yellow", 3),
  rep("pink", 3),
  rep("black", length(positions[positions$celltype == "NA", "pos"]))
)
names(coloraxis) <- positions$variable

by_sample <- ggplot(
  mydata3.genes,
  aes(
    x = start,
    xend = end,
    y = variable,
    yend =  variable,
    color = NameState
  ),
  size = 3
) +
  geom_segment(size = 3.5) +
  scale_color_manual(values = chrcolor) +
  theme_bw() +
  xlab(Chromosomename) +
  # ggtitle("Conserved chromatin states by sample")+
  ylab("Cell type") +
  theme(text = element_text(size = 12),
        axis.text.y = element_text(colour = coloraxis)) +
  geom_vline(
    xintercept = vertical.lines,
    linetype = "dashed",
    color = "dark red",
    size = 0.2
  ) +
  annotate(
    geom = "text",
    label = names(vertical.lines),
    x = vertical.lines,
    y = 23,
    size = 1 * .pt,
    angle = 40,
    hjust = 0,
    colour = "dark red"
  ) + # Adding text
  theme(plot.margin = unit(c(3, 8.7, 6, 2.1), "lines")) +
  coord_cartesian(clip = "off") +
  theme(legend.position = "none")


#Size 845*518

unique(mydata3.genes[, c("variable", "celltype")])

library(ggpubr)

ggarrange(
  bycelltype,
  by_sample,
  # bp + rremove("x.text"),
  labels = c("A", "B"),
  ncol = 1,
  nrow = 2
)

Tabla <-
  SNPsToEnh.summary[SNPsToEnh.summary$chr == chromosome,] %>% select(SNPS, DISEASE.TRAIT)

Table.plot <-
  ggtexttable(Tabla, rows = NULL)

Table.plot

#639 * 417


