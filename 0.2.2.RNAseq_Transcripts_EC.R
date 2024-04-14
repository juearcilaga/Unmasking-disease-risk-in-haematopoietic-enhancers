######################################
# Script Description:
# This script imports gene counts data from BLUEPRINT project transcript quantification 
# files. It reads in expected counts from these files using the edgeR package, 
# processes the data, and saves the resulting gene counts as an RDS file.
#
# Inputs:
# - transcript_quant directory: Directory containing transcript quantification files 
#                              downloaded from BLUEPRINT.
#
# Outputs:
# - dfcountsEC_transcripts.rds: RDS file containing the processed gene counts data.
#
# Author:
# Juliana Arcila Galvis
######################################

# Set working directory to the location of the transcript quantification files
setwd("/mnt/nocode/juliana/PID/Chromatin_States_Blueprint/")

# Step 1. Importing gene counts data from BLUEPRINT:

# List all transcript quantification files in the directory
files <- dir("transcript_quant/", full.names = TRUE)

# Extract sample names from file names
names(files) <- gsub('\\.transcript_quantification.rsem_grape2_crg.GRCh38.[0-9]*.results', '', dir("transcript_quant/"))

# Load the edgeR library for handling RNA-seq data
library(edgeR)

# Read in expected counts from transcript quantification files
countsEC <- readDGE(files, columns=c(1,5))  # tags samples (6 TPM, 5 Expected counts)
# countsTPM <- readDGE(files, columns=c(1,6))  # TPM

# Assign sample names to the columns of the counts matrix
colnames(countsEC) <- names(files)
# colnames(countsTPM) <- names(files)

# Convert counts to a data frame
dfcountsEC <- data.frame(countsEC$counts)
# dfcountsTPM <- data.frame(countsTPM$counts)

# Save processed gene counts as an RDS file
saveRDS(dfcountsEC, "dfcountsEC_transcripts.rds")
# saveRDS(dfcountsTPM, "dfcountsTPM.rds")  # Uncomment if you want to save TPM data
