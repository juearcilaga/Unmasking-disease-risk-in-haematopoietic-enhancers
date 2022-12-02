#-------------------------------------------------------------------------------------
#Which enhacers are conserved in all celltypes

library(dplyr)
library(tidyverse)

#Import my data samples in columns and chrstates in bin rows
setwd("/mnt/nocode/juliana/PID/ThemeII")
All_data<-readRDS("INPUTS/states_220422_BACKGROUND_merged.rds")

#Only enhancers
All_dataEnhancers<- All_data %>% filter_all(any_vars(. %in% c(4,5,6)))
#dim(All_dataEnhancers)#4403214     108
#400982800pb--401Mb

####Count samples with particular state 
AeloE<-function(x){x!=4}
PoisE<-function(x){x!=5}
AcanE<-function(x){x!=6}

dataAeloE<-All_dataEnhancers %>% replace(AeloE(.), 0)
dataAeloE$sum<-rowSums(dataAeloE)

dataAcanE<-All_dataEnhancers %>% replace(AcanE(.), 0)
dataAcanE$sum<-rowSums(dataAcanE)

dataPoisE<-All_dataEnhancers %>% replace(PoisE(.), 0)
dataPoisE$sum<-rowSums(dataPoisE)

###################################################################
length(dataPoisE$sum[(dataPoisE$sum/5)>2])#2473087 #Background poised
length(dataAcanE$sum[(dataAcanE$sum/6)>2])#1145017 #Background active conventional 
length(dataAeloE$sum[(dataAeloE$sum/4)>2])#608584 #Background AeloE

BAcanE<-rownames(dataAcanE[(dataAcanE$sum/6)>2,])
BAeloE<-rownames(dataAeloE[(dataAeloE$sum/4)>2,])
BPoised<-rownames(dataPoisE[(dataPoisE$sum/5)>2,])
BActiveE<-union(BAcanE,BAeloE)#1526184/4403214 

saveRDS(BPoised, "OUTPUT/BPoised.rds")
saveRDS(BActiveE, "OUTPUT/BActiveE.rds")

#-----------------------------------------------------------------
#               GR OBJECTS ANY ACTIVE ENHANCERS 
#   (Present in at least two samples of any of the 31 cell types)
#-----------------------------------------------------------------
BkGdE<-readRDS("OUTPUT//BActiveE.rds")#1526184
# BkGdE<-readRDS("OUTPUT/BPoised.rds")#2473987
BkGdE<- data.frame(coord=BkGdE, regionID=BkGdE)
BkGdE <-tidyr::separate(BkGdE, coord, c("chr", "region"))
BkGdE$end <- as.numeric(BkGdE$region) * 200
BkGdE$start <- BkGdE$end - 200

gr_BkGdE<- makeGRangesFromDataFrame(
  BkGdE,
  keep.extra.columns = TRUE,
  ignore.strand = TRUE,
  seqinfo = NULL,
  seqnames.field = c(
    "seqnames",
    "seqname",
    "chromosome",
    "chrom",
    "chr",
    "chromosome_name",
    "seqid"
  ),
  start.field = "start",
  end.field = "end",
  strand.field = NULL,
  starts.in.df.are.0based = TRUE
)
saveRDS(gr_BkGdE, "OUTPUT/gr_AEBackg.rds")
# saveRDS(gr_BkGdE, "OUTPUT/gr_PoisEBackg.rds")
