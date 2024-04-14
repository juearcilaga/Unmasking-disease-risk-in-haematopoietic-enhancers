#------------------------------------------------------------------
# DESCRIPTION:
# This script performs enrichment analysis to identify relationships 
# between disease-associated SNPs and enhancer activity patterns. It imports 
# SNPs, enhancers, and disease trait data, processes them, and then 
# performs Fisher's exact tests to determine enrichment. Finally, 
# the significant results are filtered and saved.
#
# AUTHOR:
# Juliana Arcila Galvis
#
# INPUTS:
# - gr_AE_list.rds: GenomicRanges object containing enhancers.
# - SNPsToBkGdE.rds: Dataset mapping SNPs to enhancers.
#
# OUTPUTS:
# - EnrichmentEnhparalel_2_100.rds: Enrichment results.
# - EnrichmentEnhparalel_2_100_filt.rds: Filtered enrichment results.
# - GWAS_results_2_100.txt: Initial GWAS results.
# - GWAS_results_2_100_padj.txt: Filtered GWAS results with adjusted p-values.
# - GWAS_results_0001_V2.txt: Significant GWAS results.
#
#------------------------------------------------------------------
#               Load Libraries
#------------------------------------------------------------------
setwd("/mnt/nocode/juliana/PID/ThemeII")

library(dplyr)
library(tidyr)
library(parallel)
library(MASS)

# OBJECTS REQUIRED FOR THE ENRICHMENT ANALYSIS
#-----------------------------------------------------------------
#               IMPORTING ALL ACTIVE ENHANCERS
#------------------------------------------------------------------
gr_TotalAE <- readRDS("OUTPUT/gr_AE_list.rds")

#------------------------------------------------------------------
#              Enhancers mapped to SNPs
#------------------------------------------------------------------
SNPsToBkGdAE <- readRDS("OUTPUT/SNPsToBkGdE.rds")
#------------------------------------------------------------------
#                 OUTPUT FILE'S NAMES
#------------------------------------------------------------------
filenamepvalrow <- "OUTPUT/EnrichmentEnhparalel_2_100.rds"
filenamepvalrowfilt <- "OUTPUT/EnrichmentEnhparalel_2_100_filt.rds"
GWAS_filename <- "OUTPUT/GWAS_results_2_100.txt"

#------------------------------------------------------------------
#       List of most abundant patterns in the SNPS/AE overlaps
#------------------------------------------------------------------
PatternsListSNPS <- data.frame(table(SNPsToBkGdAE$Patrons))

# Further processing and filtering of PatternsListSNPS...
dim(PatternsListSNPS)# 107562      2
summary(PatternsListSNPS$Freq)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.00     0.00     0.00     1.39     0.00 33536.00 

dim(PatternsListSNPS[PatternsListSNPS$Freq>mean(PatternsListSNPS$Freq),])
# 6867    2
write.table(PatternsListSNPS, "OUTPUT/PattersListSNPS.txt", sep="\t")

MyPatternsListSNPS <- PatternsListSNPS[2:100, "Var1"]
summary(PatternsListSNPS[2:100,"Freq"])
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 103.0   143.0   285.0   664.3   858.0  6322.0 

write.table(PatternsListSNPS[2:100,], "OUTPUT/MyPattersListSNPS.txt", sep="\t")

Low_patterns <-c("CD38negnaiveB", "CD4T","meos","mneut")

MyPatternsListSNPS[MyPatternsListSNPS %in% Low_patterns]
###00 

PatternsListSNPS[PatternsListSNPS$Var1 %in% Low_patterns, ]
###04
 #          Var1 Freq
 #         mneut  101
 # CD38negnaiveB  100
 #          meos  100
 #          CD4T   45

MyPatternsListSNPS<-c(as.character(MyPatternsListSNPS), Low_patterns)
# There are no disease associated SNPs mapping to this enhancers

#------------------------------------------------------------------
#              Mapped SNPs to Enhancers of interest
#------------------------------------------------------------------
SNPsTopPatterns <- SNPsToBkGdAE[SNPsToBkGdAE$Patrons %in% MyPatternsListSNPS, ]
dim(SNPsTopPatterns)#[1] 66115    44

diseases <- unique(SNPsTopPatterns[,"Study.Trait"])
print(length(diseases))#1006

#Traits excluded from the analysis
unwanted <-
  unique(
    c(
      "Height",
      "Economic",
      "Age at smoking initiation in chronic obstructive pulmonary disease",
      "Self-rated health",
      "Self-rated health",
      "Educational attainment",
      "Self-reported",
      "Alcohol consumption",
      "Illicit"
    )
  )

diseases <-diseases [grepl(paste(unwanted,
                           collapse = "|"),
                           diseases) == FALSE]#973 Total traits included

length(unique(SNPsToBkGdAE$Study.Trait))#[1] 1029

SNPsToBkGdAE_red <-
  SNPsToBkGdAE[SNPsToBkGdAE$Study.Trait %in% diseases, ]

length(unique(SNPsToBkGdAE_red$Study.Trait))#973
#------------------------------------------------------------------
#                ENRICHMENT
#------------------------------------------------------------------
Enrichment <- function(disease) {
	#Emodule is group of enhancers with same activity pattern)
    #How many enhancers in the Emodule?
    EGroup <- gr_TotalAE[gr_TotalAE$Patrons == Emodule,]
    NEGroup <- length(unique(EGroup$regionID))
    
    
    #How many Active enhancers in general------------------------
    IndBkGdE <-
      gr_TotalAE[(gr_TotalAE$regionID %in%
                    EGroup$regionID) == FALSE,]
    
    NIndBkGdE <- length(unique(IndBkGdE$regionID))
    
    
    #How many enhancers in the Emodule with at least 1 SNP of disease X?
    EGroupToSNPs <-
      SNPsTopPatterns[SNPsTopPatterns$Study.Trait == disease &
                        SNPsTopPatterns$Patrons == Emodule, ]
    NEGroupToSNPs <- length(unique(EGroupToSNPs$regionID))
    #1009
    
    #How many Active enhancers in general with at least 1 SNP?
    IndBkGdEToSNPs <-
      SNPsToBkGdAE_red[SNPsToBkGdAE$Study.Trait == disease,]
    NIndBkGdEToSNPs <- length(unique(IndBkGdEToSNPs$regionID))
    #564
    
    #How many Active enhancers with 0 SNPs?
    NEGroupOutSNPs <- NEGroup - NEGroupToSNPs#782525
    IndNBkEOutSNPs <- NIndBkGdE - NIndBkGdEToSNPs#742086
    
    #FISHER TEST
    Comparison <- matrix(
      c(
        NEGroupToSNPs,
        NEGroupOutSNPs,
        NIndBkGdEToSNPs,
        IndNBkEOutSNPs
      ),
      nrow = 2,
      dimnames = list(
        Lists = c("MyEnh", "IndEnh"),
        Truth = c("WithSNPs", "NoSNPs")
      )
    )
    #Example:
    # Truth
    # Lists    WithSNPs NoSNPs
    # MyEnh      1009    564
    # IndEnh   782525 742086
    
    test_greater <- fisher.test(Comparison, alternative = "greater")
    #p-value < 2.2e-16
    
    ValuesDF <- data.frame(
      DISEASE = as.character(disease),
      PATTERN = as.character(Emodule),
      AEinPATTERN = NEGroup,
      TotalAEInd = NIndBkGdE,
      AEinPATTERNToSNPs = NEGroupToSNPs,
      TotalAEIndToSNPs = NIndBkGdEToSNPs,
      AEinPATTERNoutSNPs = NEGroupOutSNPs,
      TotalAEIndoutSNPs = IndNBkEOutSNPs,
      p.value = test_greater$p.value
    )
    
    return(ValuesDF)
    # return(test_greater$p.value)
    
}

#-----------------------------------------------------------------
#               FOR EACH LIST OF ENHANCERS
#                     DO ENRICHMENT
#------------------------------------------------------------------
AllPatherns <- list()
for (Emodule in MyPatternsListSNPS) {
  All_Pval <- mclapply(diseases, Enrichment, mc.cores = 35)
  All_Pvallist <- do.call(rbind, All_Pval)
  AllPatherns[[Emodule]] <- All_Pvallist
}

AllPathernsDF <- do.call(rbind, AllPatherns)

saveRDS(AllPathernsDF, filenamepvalrow)

#------------------------------------------------------------------
#                  CORRECT FOR MULTIPLE TESTING
#					bonferroni Padj<0.0001
#------------------------------------------------------------------

All_PvalDF<- readRDS("/mnt/nocode/juliana/PID/ThemeII/OUTPUT/EnrichmentEnhparalel_2_100.rds")

dim(All_PvalDF)#100219      9

filenamepvalrowfilt<-"OUTPUT/EnrichmentEnhparalel_2_100_filt.rds"
GWAS_filename<-"OUTPUT/GWAS_results_2_100.txt"

wide_PvalDF<-pivot_wider(All_PvalDF[, c("DISEASE","PATTERN","p.value")],
                         names_from = "PATTERN",values_from ="p.value" )

PadJ <- function(x) {
  padj <- p.adjust(x, method = "bonferroni")
  return(padj)
}

filtered <-wide_PvalDF
PADJS <- lapply(filtered[2:length(filtered)], PadJ)
PADJS <- do.call(cbind, PADJS)
colnames(PADJS) <- paste("padj", colnames(PADJS), sep = "_")

filtered <- cbind(filtered, PADJS)

filtered.long <-
  pivot_longer(
    filtered,
    cols = 2:length(filtered),
    names_to = "PATTERN",
    values_to = "pValue"
  )


Padjs<-filtered.long[grepl("padj", filtered.long$PATTERN), ]
colnames(Padjs)[3]<-"Padj"
Padjs$PATTERN<-gsub("padj_","", Padjs$PATTERN)
Pvals<-filtered.long[grepl("padj", filtered.long$PATTERN)==FALSE, ]

AllResults<-left_join(All_PvalDF,Padjs)

write.table(AllResults[AllResults$Padj<0.0001,], "OUTPUT/GWAS_results_2_100_padj.txt", sep= "\t")

Significant_Results<-AllResults[AllResults$Padj<0.0001,]%>% arrange(DISEASE,PATTERN, desc( Padj))
Significant_Results<-separate(Significant_Results,col="DISEASE",into=c("TRAIT", "PUBMED.ID"), sep="_")

# Additional analysis, filtering, and saving of final results...

#Num of traits
length(unique(Significant_Results$TRAIT))#338
#Num of Pubmed ID
length(unique(Significant_Results$PUBMED.ID))#478
#Num of Patterns
length(unique(Significant_Results$PATTERN))#95

write.table(Significant_Results, "OUTPUT/GWAS_results_0001_V2.txt", sep= "\t")

