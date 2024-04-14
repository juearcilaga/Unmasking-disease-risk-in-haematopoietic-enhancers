######################################
# Script Description:
# This script processes RNAseq count data and associated metadata from BLUEPRINT
# haematopoietic cells. It imports an RNAseq count matrix, adds cell type 
# metadata, joins it with additional metadata, performs summary statistics by 
# cell type, and saves the results as text files and an RDS file.
#
# Inputs:
# - dfcountsEC.rds:matrix with expected counts from BLUEPRINT RNAseq stored as an RDS file.
# - 20160816.data: Metadata file containing information about cell types, 
#                  tissues, donors, etc., stored as a tab-delimited file.
# - gene_Anotations38_lcincluded.rds: Gene annotations stored as an RDS file.
#
# Outputs:
# - RNAEC_long.txt: Long-format RNAseq EC counts per sample with metadata.
# - Counts_sex_RNAseq.txt: Summary counts by cell type and donor sex.
# - ECAvg.rds: Summary statistics (mean, sd, etc.) by gene and cell type.
# - ECAvg.txt: Summary statistics in text format.
#
# Author:
# Juliana Arcila Galvis
######################################

# Set working directory
setwd("/mnt/nocode/juliana/PID/ThemeII")

# Import RNAseq count matrix
dfcountsEC <- readRDS("/mnt/nocode/juliana/PID/Chromatin_States_Unpublished/results/2022/RNAseq_matrix/dfcountsEC.rds")

# Add GENE_ID column using row names
dfcountsEC$GENE_ID <- rownames(dfcountsEC)

# Extract Gene_ID column
dfcountsEC_list <- dfcountsEC$GENE_ID

# Import additional metadata
all_bluprint_metadata <- read.delim("~/../../mnt/nocode/juliana/PID/Chromatin_States_Blueprint/20160816.data", 
                                    sep = "\t", 
                                    header = TRUE)

# Select and filter relevant metadata columns
metadata <- all_bluprint_metadata %>% 
  select("BIOMATERIAL_TYPE", "SAMPLE_NAME", "CELL_TYPE", "TISSUE_TYPE", "DONOR_SEX", "DONOR_ID") %>% 
  unique

# Format RNAseq count matrix to long form
dfcountsECRNAEC_long <- pivot_longer(data = dfcountsEC,
                                    cols = 1:(length(dfcountsEC) - 1),
                                    names_to = "SAMPLE_NAME",
                                    values_to = "EC")

# Add metadata to the long-form matrix
dfcountsECRNAEC_long <- left_join(dfcountsECRNAEC_long, metadata)

# Write long-form RNAseq count data with metadata to file
write.table(dfcountsECRNAEC_long, "OUTPUT/RNAEC_long.txt", sep="\t")

# Perform summary counts by cell type and donor sex
Counts_sex <- dfcountsECRNAEC_long %>%
  select("BIOMATERIAL_TYPE", "SAMPLE_NAME", "CELL_TYPE", "TISSUE_TYPE", "DONOR_SEX", "DONOR_ID") %>%
  unique %>%
  count(CELL_TYPE, DONOR_SEX)

Counts_sex <- pivot_wider(Counts_sex, names_from = "DONOR_SEX", values_from = "n")

# Write summary counts to file
write.table(Counts_sex, "OUTPUT/Counts_sex_RNAseq.txt", sep="\t")

# Import gene annotations
Genes <- readRDS("extData/gene_Anotations38_lcincluded.rds") %>% data.frame()

# Format gene annotations
Genes$Gene_names <- ifelse(Genes$gene_name == "", 
                           as.character(Genes$ID), 
                           as.character(Genes$gene_name))

colnames(Genes)[4] <- "GENE_ID"
Genes$GENE_ID <- gsub("\\..*", "", Genes$GENE_ID)
dfcountsECRNAEC_long$GENE_ID <- gsub("\\..*", "", dfcountsECRNAEC_long$GENE_ID)

# Join gene annotations with long-form RNAseq count data
Genes_annot <- Genes %>% 
  select(GENE_ID, Gene_names, gene_type) %>%
  unique()

dfcountsECRNAEC_long <- left_join(dfcountsECRNAEC_long, Genes_annot)

# Calculate summary statistics by gene and cell type
ECAvg <- dfcountsECRNAEC_long %>% 
  group_by(Gene_names, CELL_TYPE) %>% 
  summarise(VALS = paste0(EC, collapse = ", "), 
            MEAN = mean(EC), 
            SD = sd(EC), 
            GENE_ID = paste(unique(GENE_ID), collapse=", "), 
            gene_type = unique(gene_type))

# Save summary statistics
saveRDS(ECAvg, "OUTPUT/ECAvg.rds")
write.table(ECAvg, "OUTPUT/ECAvg.txt", sep="\t")
