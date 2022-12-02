#-------------------------------------------------------------------------------------
##BACKGROUND
#Step 1. Importing data:
#Data from the analysisis of the ChiPseqsamples from bluprint
#Analysis done by Juan Daviad and Enrique in Spain
#-------------------------------------------------------------------------------------
#directory with original data
setwd("/mnt/nocode/data_private/blueprint_segmentations_build38/")
setwd("DRico/chromStates_allBP_hg38_Sep2016/chromHMMresults/")
setwd("Learnmodel_12/healthy_MODEL/")

#bed files
bed_files <- dir("SEGMENTS/",  full.names = TRUE)
bed_files <- bed_files[grepl("_segments.bed", bed_files)]
names <- gsub('SEGMENTS//', '', bed_files)
names <- gsub('_12_BPhg38healthyAllsep2016', '', names)
names(bed_files) <- gsub('_segments.bed', '', names)

#state by line files
statebyline_files <- dir("STATEBYLINE/",  full.names = TRUE)
statebyline_files <-
  statebyline_files[grepl("_statebyline.txt", statebyline_files)]
names <- gsub('STATEBYLINE//', '', statebyline_files)
names <- gsub('_12_BPhg38healthyAllsep2016', '', names)
names(statebyline_files) <- gsub('_statebyline.txt', '', names)

#-------------------------------------------------------------------------------------
# EXAMPLE CELLTYPES
#-------------------------------------------------------------------------------------
#Step 4. Bulding Matrix data
#directory with original data
files<-statebyline_files

chromosomes <- c(1:22, "X","Y")

#For each chromosome
Matrix <- function(chr) {
  chr_name <- sprintf("states_chr%s.txt", chr)
  chr_name2 <- sprintf("chr%s$", chr)
  chr_name3 <- sprintf("_chr%s$", chr)
  chromosome <- files[grepl(chr_name2, names(files))]
  
  states_chr <-
    do.call(cbind, lapply(chromosome, function(fn)
      read.delim(
        fn,
        header = FALSE,
        sep = "\t",
        skip = 2
      )))
  
  colnames(states_chr) <- gsub(chr_name3, '', names(chromosome))
  
  rownames(states_chr) <-
    paste(gsub("\\$", '', chr_name2), rownames(states_chr), sep = "_")
  return(states_chr)
}

#Joining states of all the chromosomes
states_chr_sex<- do.call(rbind, lapply(chromosomes, Matrix))
ncol(states_chr_sex)

saveRDS(
  states_chr_sex,
  "/mnt/nocode/juliana/PID/Chromatin_States_Unpublished/states_220422_BACKGROUND_raw.rds"
)

#-------------------------------------------------------------------------------------
# Each column is a sample and each row is a genomic region, rownames give
# the chromosome name and the 200pb bin
#-------------------------------------------------------------------------------------
#Step 5. Merging similar states
setwd("/mnt/nocode/juliana/PID/Chromatin_States_Unpublished/")
#States_chr_all <-states_chr_sex
States_chr_all <- readRDS("/mnt/nocode/juliana/PID/Chromatin_States_Unpublished/states_220422_BACKGROUND_raw.rds")

dim(States_chr_all)#[1] [1] 15441337      108

States_chr_all[(States_chr_all == 1) | (States_chr_all == 2)] <- 1
States_chr_all[(States_chr_all == 3) | (States_chr_all == 4)] <- 2
States_chr_all[(States_chr_all == 5) | (States_chr_all == 6)] <- 3
States_chr_all[States_chr_all == 7] <- 4
States_chr_all[States_chr_all == 8] <- 5
States_chr_all[(States_chr_all == 9) | (States_chr_all == 10)] <- 6
States_chr_all[(States_chr_all == 11) | (States_chr_all == 12)] <- 7

saveRDS(States_chr_all, "states_220422_BACKGROUND_merged.rds")
#-------------------------------------------------------------------------------------

saveRDS(
  states_chr_sex,
  "/mnt/nocode/juliana/PID/Chromatin_States_Unpublished/states_220422_BACKGROUND_raw.rds"
)
saveRDS(States_chr_all, "/mnt/nocode/juliana/PID/Chromatin_States_Unpublished/states_220422_BACKGROUND_merged.rds")


