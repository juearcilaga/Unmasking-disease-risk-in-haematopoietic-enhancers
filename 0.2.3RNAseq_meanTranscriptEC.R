######################################
# Script Description:
# This script imports processed isoform counts data from BLUEPRINT project and performs 
# data transformation and visualization. It adds metadata columns, calculates 
# summary statistics, and generates heatmaps for specific transcript expressions.
#
# Inputs:
# - dfcountsEC_transcripts.rds: RDS file containing the processed isoform counts data 
#                               from BLUEPRINT.
# - 20160816.data: Metadata file containing information about various samples.
# - celltype_Abb.txt: File mapping cell type abbreviations to full names.
# - celltype_data2.txt: Additional cell type data file.
#
# Outputs:
# - ECAvg.rds: RDS file containing summary statistics of transcript expression.
# - ECAvg.txt: Text file containing summary statistics.
# - Heatmap visualizations.
#
# Author:
# Juliana Arcila Galvis
######################################

# Set working directory to the location of the theme data
setwd("/mnt/nocode/juliana/PID/ThemeII")

# Import processed isoform counts data
dfcountsEC <- readRDS("/mnt/nocode/juliana/PID/Chromatin_States_Blueprint/dfcountsEC_transcripts.rds")

# Add transcript_ID column using row names
dfcountsEC$transcript_ID <- rownames(dfcountsEC)

# Extract transcript_ID list
dfcountsEC_list <- dfcountsEC$transcript_ID

# Load tidyverse library
library(tidyverse)

# Format the data to long format and add metadata
dfcountsECRNAEC_long <- dfcountsEC %>%
  pivot_longer(cols = 1:(length(dfcountsEC) - 1), names_to = "SAMPLE_NAME", values_to = "EC") %>%
  left_join(read.delim("~/../../mnt/nocode/juliana/PID/Chromatin_States_Blueprint/20160816.data", sep = "\t", header = TRUE) %>%
              dplyr::select(BIOMATERIAL_TYPE, SAMPLE_NAME, CELL_TYPE, TISSUE_TYPE, DONOR_SEX, DONOR_ID) %>%
              unique())

# Calculate summary statistics
ECAvg <- dfcountsECRNAEC_long %>% 
  group_by(transcript_ID, CELL_TYPE) %>% 
  summarise(VALS = paste0(EC, collapse = ", "),
            MEAN = mean(EC),
            SD = sd(EC),
            Transcript_ID = paste(unique(transcript_ID), collapse=", "))

# Save summary statistics
saveRDS(ECAvg, "../../ThemeII/OUTPUT/transcriptECAvg.rds")
write.table(ECAvg, "../../ThemeII/OUTPUT/transcriptECAvg.txt", sep="\t")

