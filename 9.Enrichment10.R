#------------------------------------------------------------------
#               Load Libraries
#------------------------------------------------------------------
#Change lines 34 166  and 112 
setwd("/mnt/nocode/juliana/PID/ThemeII")

library(dplyr)
library(tidyr)
library(parallel)
library(MASS)

# OBJECTS REQUIRED FOR THE ENRICHMENT ANALYSIS
#-----------------------------------------------------------------
#               IMPORTING ANY ENHANCERS
#   (Present in at least two samples of any of the 31 cell types)
#------------------------------------------------------------------
#Genomic Ranges object with the presence/absence information [1,0]
# Of the AC in the consensus celltype chromatin states and 
# pattern name

gr_TotalAE <- readRDS("OUTPUT/gr_AE_list.rds")#49309
length(gr_TotalAE$region)#1526184
#------------------------------------------------------------------
#              Mapped SNPs to Enhancers
#------------------------------------------------------------------
SNPsToBkGdAE <- readRDS("OUTPUT/SNPsToBkGdE.rds")
dim(SNPsToBkGdAE)#150048     34
#------------------------------------------------------------------
#                 OUTPUT FILE'S NAMES
#------------------------------------------------------------------
filenamepvalrow<-"OUTPUT/EnrichmentEnhparalel_2_100.rds"
filenamepvalrowfilt<-"OUTPUT/EnrichmentEnhparalel_2_100_filt.rds"
GWAS_filename<-"OUTPUT/GWAS_results_2_100.txt"
#------------------------------------------------------------------
#       List of most abundant patterns in the total AE list
#------------------------------------------------------------------
# PatternsList <- data.frame(table(gr_TotalAE$Patrons))
# PatternsList <-
#   PatternsList[order(PatternsList$Freq, decreasing = TRUE), ]
# PatternsList$index <- 1:nrow(PatternsList)
# dim(PatternsList)
# #107562      3
# summary(PatternsList$Freq)
# # Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# # 1.0      1.0      1.0     14.2      2.0 397225.0 
# dim(PatternsList[PatternsList$Freq>mean(PatternsList$Freq),])
# # 3609    3
# MyPatternsList <- PatternsList[2:100, "Var1"]
# write.table(PatternsList, "OUTPUT/ALLPattersList.txt", sep="\t")
#------------------------------------------------------------------
#       List of most abundant patterns in the SNPS/AE overlaps
#------------------------------------------------------------------
PatternsListSNPS <- data.frame(table(SNPsToBkGdAE$Patrons))
PatternsListSNPS <-
  PatternsListSNPS[order(PatternsListSNPS$Freq, decreasing = TRUE), ]

# PatternsListSNPS$index <- 1:nrow(PatternsListSNPS)
dim(PatternsListSNPS)
# 107562      2
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

#Less than 2 samples
# falta una B unswitched
# cent mem CD8T
# cent mem CD4T
#resting endothelial

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
print(length(diseases))
#1006

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
                           diseases) == FALSE]
#973

length(unique(SNPsToBkGdAE$Study.Trait))#
#[1] 1029

SNPsToBkGdAE_red <-
  SNPsToBkGdAE[SNPsToBkGdAE$Study.Trait %in% diseases, ]

length(unique(SNPsToBkGdAE_red$Study.Trait))
#973
#------------------------------------------------------------------
#                ENRICHMENT
#------------------------------------------------------------------
Enrichment <-
  function(disease) {
    # disease = diseases[1] #"sarampion"
    # Emodule = MyPatternsListSNPS[4]
    
    #How many enhancers in the Emodule?
    # gr_TotalAE$Patrons<-as.character(gr_TotalAE$Patrons)
    # Emodule<-as.character(Emodule)
    
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
#                     ENRICHMENT
#------------------------------------------------------------------
AllPatherns <- list()
for (Emodule in MyPatternsListSNPS) {
  print(paste("Processing pathern", Emodule))
  All_Pval <- mclapply(diseases, Enrichment, mc.cores = 35 )
  All_Pvallist <- do.call(rbind, All_Pval)
  AllPatherns[[Emodule]] <- All_Pvallist
}


AllPathernsDF <- do.call(rbind, AllPatherns)
# AllPathernsDF <- do.call(cbind, AllPatherns)
# All_PvalDF <- data.frame(disease = diseases, AllPathernsDF)
# 
# colnames(All_PvalDF) <- c("disease", as.character(MyPatternsListSNPS))

saveRDS(AllPathernsDF, filenamepvalrow)
#filenamepvalrow<-"OUTPUT/EnrichmentEnhparalel_2_100.rds"
# ---------------------------------------------------------
#       Disease *enhancerModules fill=pvalue enrichmen
# https://davemcg.github.io/post/lets-plot-scrna-dotplots/
# -----------------------------------------------------------

All_PvalDF<- readRDS("/mnt/nocode/juliana/PID/ThemeII/OUTPUT/EnrichmentEnhparalel_2_100.rds")
dim(All_PvalDF)#100219      9

filenamepvalrowfilt<-"OUTPUT/EnrichmentEnhparalel_2_100_filt.rds"
GWAS_filename<-"OUTPUT/GWAS_results_2_100.txt"
#------------------------------------------------------

wide_PvalDF<-pivot_wider(All_PvalDF[, c("DISEASE","PATTERN","p.value")],
                         names_from = "PATTERN",values_from ="p.value" )

PadJ <- function(x) {
  padj <- p.adjust(x, method = "bonferroni")
  return(padj)
}

##Did it by list of enhancers or by all together??
filtered <-wide_PvalDF
PADJS <- lapply(filtered[2:length(filtered)], PadJ)
PADJS <- do.call(cbind, PADJS)
colnames(PADJS) <- paste("padj", colnames(PADJS), sep = "_")

filtered <- cbind(filtered, PADJS)

filtered.long <-
  pivot_longer(
    filtered,
    cols = 2:length(filtered),
    names_to = "Pathern",
    values_to = "pValue"
  )

# ###################################################################################
filtered.long<-
  filtered.long[grepl("padj", filtered.long$Pathern) & filtered.long$pValue < 0.0001, ]

signiflist <- list()
for (pathern in unique(filtered.long$Pathern)) {
  PadjPathern <- filtered.long[filtered.long$Pathern == pathern,]
  PadjPathern <- PadjPathern[order(PadjPathern$pValue),]
  PadjPathern <- PadjPathern[PadjPathern$pValue < 0.05,]
  print(pathern)
  print(dim(PadjPathern))
  signiflist[[pathern]] <- PadjPathern
}


Stable_GWAS <- do.call(rbind, signiflist)
Stable_GWAS$Pathern <- gsub("padj_", "", Stable_GWAS$Pathern)
colnames(Stable_GWAS) <- c("disease", "EnhancerGroup", "padj")
# A tibble: 1,982 × 3
# A tibble: 1,117 × 3

Stable_GWAS$ID<-paste(Stable_GWAS$disease, Stable_GWAS$EnhancerGroup, sep="|")
METADATA<-unique(All_PvalDF[, c(1,2,3,5)])
METADATA$ID<-paste(METADATA$DISEASE, METADATA$PATTERN, sep="|")

Stable_GWA<-left_join(Stable_GWAS, METADATA)
Stable_GWA<-Stable_GWA[,c("disease", "EnhancerGroup",
                          "padj", "AEinPATTERN", "AEinPATTERNToSNPs")]
summary(Stable_GWA$AEinPATTERNToSNPs)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 2.00    6.00    9.00   11.03   13.00  102.00 
GWAS_filename<-"OUTPUT/GWAS_results_0001.txt"
write.table(Stable_GWA, GWAS_filename, sep = "\t")


