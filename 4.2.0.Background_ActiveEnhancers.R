#-------------------------------------------------------------------------------------
# Author: Juliana Arcila Galvis
# Description: 
#   This script identifies enhancers that are conserved across all cell types.
#   It takes enhancer data as input, filters it to retain only enhancers,
#   and then calculates the count of samples with specific enhancer states.
#   Finally, it generates GRanges objects for any active enhancers present 
#   in at least two samples of any of the 31 cell types.
#
# Inputs:
#   - "states_220422_BACKGROUND_merged.rds": RDS file containing chromatin states data
#
# Outputs:
#   - "BPoised.rds": RDS file containing identified poised enhancers
#   - "BActiveE.rds": RDS file containing identified active enhancers
#   - "gr_AEBackg.rds": RDS file containing GRanges object for identified enhancers
#
#-------------------------------------------------------------------------------------

# Load required libraries
library(dplyr)
library(tidyverse)

# Set working directory and import data
setwd("/mnt/nocode/juliana/PID/ThemeII")
All_data <- readRDS("INPUTS/states_220422_BACKGROUND_merged.rds")

# Filter data to retain only enhancers
All_dataEnhancers <- All_data %>% filter_all(any_vars(. %in% c(4, 5, 6)))

#### Define functions to count samples with specific enhancer states
AeloE <- function(x) {x != 4}
PoisE <- function(x) {x != 5}
AcanE <- function(x) {x != 6}

# Replace enhancer states based on defined functions and calculate row sums
dataAeloE <- All_dataEnhancers %>% replace(AeloE(.), 0)
dataAeloE$sum <- rowSums(dataAeloE)

dataAcanE <- All_dataEnhancers %>% replace(AcanE(.), 0)
dataAcanE$sum <- rowSums(dataAcanE)

dataPoisE <- All_dataEnhancers %>% replace(PoisE(.), 0)
dataPoisE$sum <- rowSums(dataPoisE)

# Count and store enhancers with specific states
length(dataPoisE$sum[(dataPoisE$sum / 5) > 2])  # 2473087 - Background poised
length(dataAcanE$sum[(dataAcanE$sum / 6) > 2])   # 1145017 - Background active conventional
length(dataAeloE$sum[(dataAeloE$sum / 4) > 2])   # 608584  - Background AeloE

# Extract enhancers with specific states
BAcanE <- rownames(dataAcanE[(dataAcanE$sum / 6) > 2, ])
BAeloE <- rownames(dataAeloE[(dataAeloE$sum / 4) > 2, ])
BPoised <- rownames(dataPoisE[(dataPoisE$sum / 5) > 2, ])
BActiveE <- union(BAcanE, BAeloE)  # 1526184/4403214 

# Save identified enhancers
saveRDS(BPoised, "OUTPUT/BPoised.rds")
saveRDS(BActiveE, "OUTPUT/BActiveE.rds")

#-----------------------------------------------------------------
# Generate GRanges objects for any active enhancers 
# (Present in at least two samples of any of the 31 cell types)
#-----------------------------------------------------------------

# Read saved enhancer data
BkGdE <- readRDS("OUTPUT/BActiveE.rds")  # 1526184
# BkGdE <- readRDS("OUTPUT/BPoised.rds")  # 2473987

# Convert enhancer data to dataframe and separate coordinates
BkGdE <- data.frame(coord = BkGdE, regionID = BkGdE)
BkGdE <- tidyr::separate(BkGdE, coord, c("chr", "region"))

# Define start and end coordinates
BkGdE$end <- as.numeric(BkGdE$region) * 200
BkGdE$start <- BkGdE$end - 200

# Create GRanges object
gr_BkGdE <- makeGRangesFromDataFrame(
  BkGdE,
  keep.extra.columns = TRUE,
  ignore.strand = TRUE,
  seqinfo = NULL,
  seqnames.field = c(
    "seqnames", "seqname", "chromosome", "chrom", "chr", "chromosome_name", "seqid"
  ),
  start.field = "start",
  end.field = "end",
  strand.field = NULL,
  starts.in.df.are.0based = TRUE
)

# Save GRanges object
saveRDS(gr_BkGdE, "OUTPUT/gr_AEBackg.rds")
# saveRDS(gr_BkGdE, "OUTPUT/gr_PoisEBackg.rds")
