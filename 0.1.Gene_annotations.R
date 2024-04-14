
######################################
# Script Description:
# This script organizes gene annotations obtained from ENCODE release 38.
# It filters the annotations to retain only those related to protein-coding genes,
# selects relevant columns, and creates a GRanges object to store the genomic ranges
# and gene information. The resulting data is saved as RDS files for further analysis.
#
# Inputs:
# - gencode.v38.basic.annotation.gff3: Gene annotation table downloaded from ENCODE release 38.
#
# Outputs:
# - gene_Anotations38_lcincluded.rds: Filtered and selected gene annotations saved as an RDS file.
# - gr_gene_Anotations38_lcincluded.rds: GRanges object containing gene information saved as an RDS file.
#
# Author:
# Juliana Arcila Galvis
######################################

# Step 0: Organising gene annotations

# Set working directory to the location of the gene annotation file
setwd("/mnt/nocode/juliana/PID/ThemeII/extData/")

# Load required libraries
library(rtracklayer)

# Import gene annotation table downloaded from ENCODE release 38 in GFF3 format
GENE_CODE_38 <- readGFF("gencode.v38.basic.annotation.gff3", 
                        version=3, 
                        columns=NULL, 
                        tags=NULL, 
                        filter=NULL, 
                        nrows=-1,
                        raw_data=FALSE)

# Filter out annotations for protein-coding genes
# Uncomment below lines if you want to filter
# GENES <- GENE_CODE_38[GENE_CODE_38$gene_type == "protein_coding" & 
#                       GENE_CODE_38$type == "gene",]

# Assign filtered annotations to a new variable
# annotation_genes <- GENES
annotation_genes <- GENE_CODE_38

# Display dimensions of the annotation_genes data frame
dim(annotation_genes) # [1] 20082    35

# Filter out annotations with gene_name present
genes <- annotation_genes[annotation_genes$gene_name != "",]

# Display dimensions of the genes data frame
dim(genes) # [1] 20082    35

# Select relevant columns from genes data frame
genes <- genes[,c(1,4:5,9,11,12)]

# Save filtered and selected genes data frame as an RDS file
saveRDS(genes, "gene_Anotations38_lcincluded.rds")

###############################################################

# Creating a granges object for gene information

# Read previously saved genes data frame from RDS file
genes <- readRDS("gene_Anotations38_lcincluded.rds")

# Create a GRanges object to store genomic ranges and gene information
gr_gene <- GRanges(
  seqnames = genes$seqid,
  ranges = IRanges(genes$start, end = genes$end, names = genes$ID),
  gene_name = genes$gene_name,
  gene_ID = genes$ID,
  gene_start = genes$start,
  gene_end = genes$end,
  gene_size = genes$end - genes$start,
  gene_type = genes$gene_type
)

# Save the GRanges object as an RDS file
saveRDS(gr_gene, "gr_gene_Anotations38_lcincluded.rds")

# Uncomment below line to display the gr_gene object
# gr_gene
