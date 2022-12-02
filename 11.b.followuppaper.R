
#------------------------------------------------------------------------------------
#                           IMPORTING PVALS
#------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------
#               Building Matrix Disease*Pattern filling with NA
#------------------------------------------------------------------------------------

GWAS_filename<-read.delim("OUTPUT/GWAS_results_0001.txt")
colnames(GWAS_filename)
#[1] "disease"           "EnhancerGroup"     "padj"              "AEinPATTERN"       "AEinPATTERNToSNPs"   


summary(GWAS_filename$padj)
#0.000e+00 0.000e+00 3.560e-09 6.891e-06 1.913e-06 9.874e-05 

library(tidyverse)

summary(GWAS_filename$AEinPATTERNToSNPs)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 3.00    8.00   11.00   13.66   16.00  102.00

GWAS.sign<-GWAS_filename

# GWAS.sign<-GWAS.sign %>% filter(AEinPATTERNToSNPs>=mean(AEinPATTERNToSNPs))
# summary(GWAS.sign$AEinPATTERNToSNPs)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 14.00   16.00   20.00   23.23   25.25  102.00 

GWAS.sign<-GWAS.sign %>% separate(disease, into=c("DISEASE", "PUBMED.ID"), sep="_")


M012.cells <- c("M0","M1","M2")


lipid.traits <- c(
  "LDL cholesterol",
  "Metabolic syndrome (bivariate traits)",
  "Metabolite levels",
  "Waist-hip ratio",
  "Lipid metabolism phenotypes",
  "HDL cholesterol",
  "HDL Cholesterol - Triglycerides (HDLC-TG)",
  "Triglycerides",
  "Triglycerides-Blood Pressure (TG-BP)",
  "Metabolic traits",
  "Metabolic syndrome"
)

M012.lipid.enrich <-
  GWAS.sign[grepl(paste0(lipid.traits, collapse = "|"), GWAS.sign$DISEASE)  &
                   grepl(paste0(M012.cells, collapse = "|"),
                         GWAS.sign$EnhancerGroup), ]



M012.lipid.enrich$TRAIT.PUBMED.ENhList<-M012.lipid.enrich %>% 
  dplyr::select(DISEASE, PUBMED.ID, EnhancerGroup) %>% 
  unite(., TRAIT.PUBMED.ENhList, sep="_")
  
SNPsToEnhInterest <- readRDS("OUTPUT/SNPsToBkGdE.rds")

colnames(SNPsToEnhInterest)[1:5] <-c("chr", "Start", "End", "Width", "Strand")

SNPsToEnhM0123 <- 
  SNPsToEnhInterest[,c("chr", "Start", "End", "Width","regionID", 
                    "Patrons","PUBMEDID", "CHR_ID", "CHR_POS", 
                    "SNPS", "DISEASE.TRAIT", "CONTEXT")]
SNPsToEnhM0123<-
  SNPsToEnhM0123[SNPsToEnhM0123$Patrons%in%M012.lipid.enrich$EnhancerGroup,]

SNPsToEnhM0123<-
  SNPsToEnhM0123[SNPsToEnhM0123$DISEASE.TRAIT%in%M012.lipid.enrich$DISEASE,]

SNPsToEnhM0123<-
  SNPsToEnhM0123[SNPsToEnhM0123$PUBMEDID%in%M012.lipid.enrich$PUBMED.ID,]

SNPsToEnhM0123$TRAIT.PUBMED.ENhList <- SNPsToEnhM0123 %>%
  dplyr::select(DISEASE.TRAIT, PUBMEDID, Patrons) %>%
  unite(., TRAIT.PUBMED.ENhList, sep = "_")

SNPsToEnhM0123<-
  SNPsToEnhM0123[
    SNPsToEnhM0123$TRAIT.PUBMED.ENhList$TRAIT.PUBMED.ENhList %in% 
    M012.lipid.enrich$TRAIT.PUBMED.ENhList$TRAIT.PUBMED.ENhList,]


SNPsToEnhM0123.summary<-
  SNPsToEnhM0123 %>%
  group_by(SNPS) %>%
  summarise(
    DISEASE.TRAIT = paste(unique(DISEASE.TRAIT), collapse = ", "),
    chr = paste(unique(chr), collapse = ", "),
    NEnh = length(unique(regionID)),
    CHR_POS= unique(CHR_POS),
    start= min(Start),
    end= max(End),
    EnhGroup = paste(unique(Patrons), collapse = ", ")
  )

SNPsToEnhM0123.summary <-
  SNPsToEnhM0123.summary[order(SNPsToEnhM0123.summary$chr, decreasing = TRUE),]

SNPsToEnhM0123.summary$size<-SNPsToEnhM0123.summary$end-SNPsToEnhM0123.summary$start


table(SNPsToEnhM0123.summary$chr) %>% data.frame %>% arrange(desc(Freq))
#     CHROM NSNPS
# 1  chr8   16
# 2  chr19    6
# 3  chr15    4
# 4  chr11,2,6,9    2
# 8  chr10,12,20,22,5,X,    1



#Only those with more than 10 SNPS

SNPsToEnhM0123.summary<-SNPsToEnhM0123.summary[SNPsToEnhM0123.summary$NEnh>10,]
# rs5031002 16 snps and is in the X chromosome, should be relevant by sex, check for Ines
write.table(SNPsToEnhM0123.summary,"OUTPUT/SNPsToEnhM0123summary.txt", sep="\t")
#In chromosome X is only the LDL cholesterol 16 ehancers. That should be sex specific.
#In chromosome 11 is only the HDL cholesterol 13 enhancers 
# Chromosome 6 waist hip radio
# Metabolite levels  chr 2 and 22
#-----------------------------------------------------------------------------------------------------

SNPsToEnhM0123.summary.bychr<-
  SNPsToEnhM0123.summary %>% group_by(chr) %>%
  summarise(
    startReg = min(start),
    endReg = max(end),
    sizeReg = endReg- startReg
  ) 


SNPsToEnhM0123.summary.bychr$middle <-
  SNPsToEnhM0123.summary.bychr$startReg + 
  (SNPsToEnhM0123.summary.bychr$sizeReg /2)


SNPsToEnhM0123.summary.bychr$start10K<-SNPsToEnhM0123.summary.bychr$middle - 5000
SNPsToEnhM0123.summary.bychr$end10K<-SNPsToEnhM0123.summary.bychr$middle + 5000

SNPsToEnhM0123.summary.bychr<-
  SNPsToEnhM0123.summary.bychr %>% 
  mutate(start20k=ifelse(start10K>startReg, startReg,start10K),
         end20k=ifelse(end10K<endReg, endReg, end10K) )%>% 
  dplyr::select(chr,start20k,end20k)

# SNPsToEnhM0123.summary$start50K<-SNPsToEnhM0123.summary$middle-50000
# SNPsToEnhM0123.summary$end50K<-SNPsToEnhM0123.summary$middle+50000
write.table(SNPsToEnhM0123.summary.bychr,"../OUTPUT/SNPsToEnhM0123.summary.bychr.txt", sep="\t")



   
##Anotate the region for plotting   
library(plyranges)
gene_Anotations<-readRDS("gr_gene_Anotations38.rds")
gene_Anotations[gene_Anotations$gene_name=="LPL",]

SNPS.GENES<-list()

for(chr in 1:length(SNPsToEnhM0123.summary.bychr$chr)) {
  CHR <- SNPsToEnhM0123.summary.bychr[chr, ]
  gr.CHR <- GRanges(seqnames = CHR$chr ,
                    ranges = IRanges(start = CHR$start20k, end = CHR$end20k ))
  CHR.GENES <-
    na.omit(data.frame(
      plyranges::join_overlap_intersect(gene_Anotations, gr.CHR)))
  
  SNPS.GENES[[CHR$chr]] <- CHR.GENES
}

SNPS.GENES.DF<-do.call(rbind, SNPS.GENES)

# APOE, chr 19 : Mutations in this gene result in familial dysbetalipoproteinemia, 
# or type III hyperlipoproteinemia (HLP III), in which increased plasma 
# cholesterol and triglycerides are the consequence of impaired clearance 
# of chylomicron and VLDL remnants. 

# NR1H3, chr 11

# APOC1, chr 19 : This gene is expressed primarily in the liver, and it is activated when 
# monocytes differentiate into macrophages. The encoded protein plays a 
# central role in high density lipoprotein (HDL) and very low density lipoprotein 
# (VLDL) metabolism. 

# PNPLA3, chr 22: The protein encoded by this gene is a triacylglycerol lipase that mediates triacylglycerol 
# hydrolysis in adipocytes. The encoded protein, which appears to be membrane bound, may be 
# involved in the balance of energy usage/storage in adipocytes

# LPL, chr 8, # LPL lipoprotein lipasa
# The Lipoprotein lipase (LPL) gene is a significant contributor to dyslipidemia. 
# It has shown associations with several conditions including atherosclerosis, obesity, 
# and metabolic syndrome (MetS). We assessed the interactive association between MetS and 
# rs3779788 of the LPL gene based on aerobic exercise.


# AR: The androgen receptor gene, Upon binding the hormone ligand, 
# the receptor dissociates from accessory proteins, translocates into 
# the nucleus, dimerizes, and then stimulates transcription of androgen 
# responsive genes 

write.table(SNPS.GENES.DF,"../OUTPUT/Lipid.GENES.DF.txt", sep="\t")

#Plot
# Consesnsus M0, M1, M2 chromstates
# vertical lines with the SNPs
# Color of vertical lines depends on the diseases associated
# Panel with coordinates of the LPL gene

