#-----------------------------------------------------------------
# DESCRIPTION:
# This script maps SNPs to enhancers of interest based on overlapping genomic regions.
# It imports filtered data from the GWAS catalog and active enhancer datasets. 
# It then filters and processes the data to map SNPs to enhancers.
#
# AUTHOR:
# Juliana Arcila Galvis
#
# INPUTS:
# - gr_gwas_cat_reduc.rds: GenomicRanges object containing filtered 
#   data from GWAS catalog.
# - gr_AE_list.rds: GenomicRanges object containing active enhancers 
#   across different cell types.
#
# OUTPUTS:
# - SNPsToBkGdE.rds: Filtered dataset mapping SNPs to enhancers.
#
#-----------------------------------------------------------------
#               SETTING DIRECTORY
#------------------------------------------------------------------
setwd("/mnt/nocode/juliana/PID/ThemeII")

#------------------------------------------------------------------
#              Mapping SNPs to Background and Enhancers of interest
#------------------------------------------------------------------

#------------------------------------------------------------------
# IMPORTING DATA FROM GWAS CATALOG DB
#------------------------------------------------------------------
gr_gwas_cat_reduc <- readRDS("OUTPUT/gr_gwas_cat_reduc.rds")

#------------------------------------------------------------------
#               IMPORTING ENHANCERS 
#   (Present in at least two samples of any of the 31 cell types)
#------------------------------------------------------------------
gr_CtypeChrStates <- readRDS("OUTPUT/gr_AE_list.rds")
length(gr_CtypeChrStates$region)

gr_CtypeChrStates_NoY <- gr_CtypeChrStates[seqnames(gr_CtypeChrStates) != "chrY"]

#------------------------------------------------------------------
#              Mapping SNPs to All enhancers
#------------------------------------------------------------------
SNPsToBkGdAE <- na.omit(data.frame(plyranges::join_overlap_intersect(
  gr_CtypeChrStates_NoY, gr_gwas_cat_reduc)))

seqnames(gr_gwas_cat_reduc)
seqnames(gr_CtypeChrStates_NoY)

dim(SNPsToBkGdAE)

#------------------------------------------------------------------
# Filtering unwanted SNPs and enhancers
#------------------------------------------------------------------
unwanted <- c("missense", "STOP-GAIN", "frameshift")

SNPsToBkGdAE <- SNPsToBkGdAE[grepl(paste(unwanted, collapse = "|"),
                                    SNPsToBkGdAE$CONTEXT) == FALSE, ]

summary(SNPsToBkGdAE$width)

# Filtering by width (assuming width represents length in this context)
SNPsToBkGdAE <- SNPsToBkGdAE[SNPsToBkGdAE$width >= 200, ]

dim(SNPsToBkGdAE)

#------------------------------------------------------------------
# Saving the processed dataset
#------------------------------------------------------------------
saveRDS(SNPsToBkGdAE, "OUTPUT/SNPsToBkGdE.rds")
