######################################
# Script Description:
# This script imports and organizes metadata and data files related to ChIP-seq 
# samples from the BLUEPRINT project. It imports various metadata files, filters 
# and merges them, and prepares a list of transcript quantification files for download.
#
# Inputs:
# - 20160816.data: Metadata file containing information about various samples.
# - METADATA_paper2.txt: Additional metadata file for paper-related information.
#
# Outputs:
# - transcript_quant_BLUEPTINT.txt: List of URLs for transcript quantification files.
# - transcript_quant directory: Directory containing downloaded transcript quantification files.
#
# Author:
# Juliana Arcila Galvis
######################################

# Set working directory to the location of the original data
setwd("/mnt/nocode/juliana/PID/Chromatin_States_Blueprint/")

# Uncomment and modify as needed if you want to process these files

# Step 1. Importing all metadata:
# Import metadata for all samples from BLUEPRINT
all_bluprint_metadata <- read.delim("20160816.data", 
                                    sep = "\t", 
                                    header = TRUE)

# Step 2. Importing paper metadata:
# Import metadata for samples from the paper
all_paper_metadata <- read.delim("../ThemeII/INPUTS/METADATA_paper2.txt", 
                                 sep = "\t", 
                                 header = TRUE)

# Filter metadata to include only cell types with chromatin states
df <- all_bluprint_metadata[
  all_bluprint_metadata$CELL_TYPE %in% all_paper_metadata$CELL_TYPE,
]

# Extract transcript quantification files and prepare URLs for download
files_transcript <- df$FILE[grep("transcript_quantification", df$FILE)]
files_transcript <- paste("https://ftp.ebi.ac.uk/pub/databases/", files_transcript, sep = "")

# Write URLs to file
write.table(files_transcript, "transcript_quant_BLUEPTINT.txt", 
            col.names = FALSE, 
            row.names = FALSE, 
            quote = FALSE)

# Download transcript quantification files
system("mkdir transcript_quant")
system("cd transcript_quant")
system("wget -i ../transcript_quant_BLUEPTINT.txt")
