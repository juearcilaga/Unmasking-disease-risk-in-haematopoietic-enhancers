#-------------------------------------------------------------------------------
# Script Name: Blueprint_Data_Analysis.R
# Author: Juliana Arcila
# Description: 
# This script performs data analysis on ChiPseq and RNAseq samples from BLUEPRINT.
# It imports the data, filters the metadata, plots sample sizes, and subsets 
# Monocyte-derived samples. The output includes metadata files and sample size 
# plots for further analysis.
#
# Inputs:
# - ChiPseq and RNAseq data from BLUEPRINT
# - Metadata files: 20160816.data, METADATA_paper2.txt, celltype_data.txt
#
# Outputs:
# - Metadata files for ChiPseq and RNAseq: Metadata_ChiPseq.txt, Metadata_RNAseq.txt
# - Sample size plot: Number of Biological Replicates by CELL_TYPE

#-------------------------------------------------------------------------------
# Step 1. Importing data:
# Data from the analysis of the ChiPseq samples from BLUEPRINT
# Analysis done by Juan David and Enrique Carrillo in Spain
#-------------------------------------------------------------------------------

# Directory with original data
setwd("/mnt/nocode/data_private/blueprint_segmentations_build38/")
setwd("DRico/chromStates_allBP_hg38_Sep2016/chromHMMresults/")
setwd("Learnmodel_12/healthy_MODEL/")

# Importing bed files related to chromatin segments
bed_files <- dir("SEGMENTS/",  full.names = TRUE)
bed_files <- bed_files[grepl("_segments.bed", bed_files)]
names <- gsub('SEGMENTS//', '', bed_files)
names <- gsub('_12_BPhg38healthyAllsep2016', '', names)
names(bed_files) <- gsub('_segments.bed', '', names)

# Importing state by line files
statebyline_files <- dir("STATEBYLINE/",  full.names = TRUE)
statebyline_files <- statebyline_files[grepl("_statebyline.txt", statebyline_files)]
names <- gsub('STATEBYLINE//', '', statebyline_files)
names <- gsub('_12_BPhg38healthyAllsep2016', '', names)
names(statebyline_files) <- gsub('_statebyline.txt', '', names)

#-------------------------------------------------------------------------------
# Step 2. Importing metadata:
# Metadata of samples from BLUEPRINT
#-------------------------------------------------------------------------------
all_bluprint_metadata <- read.delim(
  "~/../../mnt/nocode/juliana/PID/Chromatin_States_Blueprint/20160816.data",
  sep = "\t",
  header = TRUE
)

# Filtering metadata based on bed file names
df <- all_bluprint_metadata[all_bluprint_metadata$SAMPLE_NAME %in% names(bed_files), ]

# Selecting relevant columns for further analysis
SubDf <- dplyr::select(
  df,
  "EXPERIMENT_ID",
  "STUDY_ID",
  "BIOMATERIAL_TYPE",
  "SAMPLE_NAME",
  "CELL_TYPE",
  "TISSUE_TYPE",
  "DONOR_SEX",
  "DONOR_ID",
  "DONOR_HEALTH_STATUS",
  "EXPERIMENT_TYPE"
)

# Filtering unique subset
df_sub <- unique(SubDf)

# Filtering OTHER_DATASETS based on selected criteria
OTHER_DATASETS <- all_bluprint_metadata %>% 
  dplyr::filter(
    EXPERIMENT_TYPE %in% c("total-RNA-Seq", "mRNA-Seq", "Chromatin Accessibility") &
    CELL_TYPE %in% SubDf$CELL_TYPE &
    DONOR_HEALTH_STATUS %in% unique(SubDf$DONOR_HEALTH_STATUS)
  )

# Filtering RNAseq samples
RNAseq <- OTHER_DATASETS %>% dplyr::filter(FILE_TYPE == "RNA_GENE_QUANT_STAR_CRG") %>%
  dplyr::select(
    "EXPERIMENT_ID",
    "STUDY_ID",
    "BIOMATERIAL_TYPE",
    "SAMPLE_NAME",
    "CELL_TYPE",
    "TISSUE_TYPE",
    "DONOR_SEX",
    "DONOR_ID",
    "EXPERIMENT_TYPE"
  ) %>% unique

# Filtering ChiPseq samples
ChiPseq <- SubDf %>% dplyr::select(-c("DONOR_HEALTH_STATUS")) %>% unique

# Saving metadata to output files
setwd("/mnt/nocode/juliana/PID/ThemeII")
write.table(ChiPseq, "OUTPUT/TABLES/Metadata_ChiPseq.txt", sep = "\t")
write.table(RNAseq, "OUTPUT/TABLES/Metadata_RNAseq.txt", sep = "\t")

#-------------------------------------------------------------------------------
# Step 3. Sample size plot:
# Plotting sample sizes by cell type
#-------------------------------------------------------------------------------
mymetadata <- read.delim("INPUTS/METADATA_paper2.txt", sep = "\t")
celltypes <- read.delim("INPUTS/celltype_data.txt")

# Counting sample sizes
Sample_number <- table(mymetadata$CELL_TYPE)
Sample_number <- data.frame(Sample_number[Sample_number > 0])
names(Sample_number) <- c("CELL_TYPE", "SAMPLE_SIZE")

# Joining with additional cell type data
Sample_number <- left_join(Sample_number, celltypes)
Sample_number <- left_join(
  Sample_number,
  unique(mymetadata[,c("CELL_TYPE", "CELL_TYPE_Abb")])
)

# Plotting sample sizes
p <- ggdotchart(
  Sample_number,
  y = "SAMPLE_SIZE",
  x = "CELL_TYPE_Abb",
  color = "Linage",
  palette = c("#00AFBB", "#E7B800", "#FC4E07"),
  add = "segments",
  rotate = TRUE,
  group = "CELL_TYPE_Abb",
  dot.size = 2,
  font.label = list(
    color = "white",
    size = 2,
    vjust = 0.5
  ),
  xlab = "",
  ylab = "Number of Biological Replicates",
  legend.title = "",
  ggtheme = theme_pubr()
)

# Adjusting plot margins
p <- p + theme(plot.margin = margin(t = 0, r = 20, b = 0, l = 40))

