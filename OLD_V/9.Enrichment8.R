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

gr_TotalAE$Patrons[grepl("imDC",unique(gr_TotalAE$Patrons))]
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
PatternsList <- data.frame(table(gr_TotalAE$Patrons))
PatternsList <-
  PatternsList[order(PatternsList$Freq, decreasing = TRUE), ]
PatternsList$index <- 1:nrow(PatternsList)
dim(PatternsList)
#107562      3
summary(PatternsList$Freq)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 1.0      1.0      1.0     14.2      2.0 397225.0 
dim(PatternsList[PatternsList$Freq>mean(PatternsList$Freq),])
# 3609    3
MyPatternsList <- PatternsList[2:100, "Var1"]
write.table(PatternsList, "OUTPUT/ALLPattersList.txt", sep="\t")
#------------------------------------------------------------------
#       List of most abundant patterns in the SNPS/AE overlaps
#------------------------------------------------------------------
PatternsListSNPS <- data.frame(table(SNPsToBkGdAE$Patrons))
PatternsListSNPS <-
  PatternsListSNPS[order(PatternsListSNPS$Freq, decreasing = TRUE), ]

# PatternsListSNPS$index <- 1:nrow(PatternsListSNPS)
dim(PatternsListSNPS)
# 12901     3
summary(PatternsListSNPS$Freq)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.00     0.00     0.00    11.63     1.00 61771.00 
dim(PatternsListSNPS[PatternsListSNPS$Freq>mean(PatternsListSNPS$Freq),])
# 642   3

MyPatternsListSNPS <- PatternsListSNPS[2:100, "Var1"]
write.table(PatternsListSNPS, "OUTPUT/PattersListSNPS.txt", sep="\t")

#FROM SCRIT 7
Low_patterns <- c(
  "bfneut",
  "megkaryo",
  "CD8T",
  "NK",
  "endumbprolif",
  "eryth",
  "germcenterB",
  "neutmetmyelo",
  "neutmyelo",
  "plasma",
  "segneutBM"
)

#Less than 2 samples
# falta una B unswitched
# cent mem CD8T
# cent mem CD4T

PatternsListSNPS[PatternsListSNPS$Var1 %in% Low_patterns, ]
###00 
# There are no disease associated SNPs mapping to this enhancers
#------------------------------------------------------------------
#              Mapped SNPs to Enhancers of interest
#------------------------------------------------------------------
SNPsTopPatterns <- SNPsToBkGdAE[SNPsToBkGdAE$Patrons %in% MyPatternsListSNPS, ]
dim(SNPsTopPatterns)#[1] 63992    34

diseases <- unique(SNPsTopPatterns[,"Study.Trait"])
print(length(diseases))
#979

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
#946

length(unique(SNPsToBkGdAE$Study.Trait))#
#[1] 1029

SNPsToBkGdAE_red <-
  SNPsToBkGdAE[SNPsToBkGdAE$Study.Trait %in% diseases, ]

length(unique(SNPsToBkGdAE_red$Study.Trait))
#946
#------------------------------------------------------------------
#                ENRICHMENT
#------------------------------------------------------------------
# Enrichment <-
#   function(disease) {
    disease = "Inflammatory bowel disease" #"sarampion"
    Emodule = MyPatternsListSNPS[4]
    
    #How many enhancers in the Emodule?
    gr_TotalAE$Patrons<-as.character(gr_TotalAE$Patrons)
    Emodule<-as.character(Emodule)
    
    EGroup <- gr_TotalAE[gr_TotalAE$Patrons == Emodule, ]
    NEGroup <- length(unique(EGroup$regionID))
    
    data.frame(DISEASE=as.character(disease),
               PATTERN=as.character(Emodule),
               NAEPATTERN <-NEGroup)###OJO, A LISTA
    
    #How many Active enhancers in general------------------------
    IndBkGdE <-
      gr_TotalAE[(gr_TotalAE$regionID %in%
                           EGroup$regionID) == FALSE, ]
    
    NIndBkGdE <- length(unique(IndBkGdE$regionID))
    
    #How many enhancers in the Emodule with at least 1 SNP of disease X?
    EGroupToSNPs <-
      SNPsTopPatterns[SNPsTopPatterns$Study.Trait == disease &
                     SNPsTopPatterns$Patrons == Emodule,]
    NEGroupToSNPs <- length(unique(EGroupToSNPs$regionID))
    #1009
    
    #How many Active enhancers in general with at least 1 SNP?
    IndBkGdEToSNPs <-
      SNPsToBkGdAE_red[SNPsToBkGdAE$Study.Trait == disease, ]
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
    return(test_greater$p.value)
    
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

AllPathernsDF <- do.call(cbind, AllPatherns)
All_PvalDF <- data.frame(disease = diseases, AllPathernsDF)

colnames(All_PvalDF) <- c("disease", as.character(MyPatternsListSNPS))

saveRDS(All_PvalDF, filenamepvalrow)

# ---------------------------------------------------------
#       Disease *enhancerModules fill=pvalue enrichmen
# https://davemcg.github.io/post/lets-plot-scrna-dotplots/
# -----------------------------------------------------------

filtered <-
  All_PvalDF[rowSums(All_PvalDF[, 2:length(MyPatternsListSNPS)] < 0.05) > 0, ]###Question here

saveRDS(filtered, filenamepvalrowfilt)

# filtered <- readRDS("RDS/EnrichmentEnhparalel.rds")
filtered2 <- filtered

PadJ <- function(x) {
  padj <- p.adjust(x, method = "bonferroni")
  return(padj)
}

PADJS <- lapply(filtered2[2:length(filtered2)], PadJ)
PADJS <- do.call(cbind, PADJS)
colnames(PADJS) <- paste("padj", colnames(PADJS), sep = "_")

filtered3 <- cbind(filtered2, PADJS)

df.long <-
  pivot_longer(
    filtered3,
    cols = 2:length(filtered3),
    names_to = "Pathern",
    values_to = "pValue"
  )

###################################################################################
Padjusteds <-
  df.long[grepl("padj", df.long$Pathern) & df.long$pValue < 0.05, ]

signiflist <- list()
for (pathern in unique(Padjusteds$Pathern)) {
  PadjPathern <- Padjusteds[Padjusteds$Pathern == pathern,]
  PadjPathern <- PadjPathern[order(PadjPathern$pValue),]
  PadjPathern <- PadjPathern[PadjPathern$pValue < 0.05,]
  print(pathern)
  print(dim(PadjPathern))
  signiflist[[pathern]] <- PadjPathern
}


Stable_GWAS <- do.call(rbind, signiflist)
Stable_GWAS$Pathern <- gsub("padj_", "", Stable_GWAS$Pathern)
colnames(Stable_GWAS) <- c("disease", "EnhancerGroup", "padj")

write.table(Stable_GWAS,
            GWAS_filename,
            sep = "\t")

# 
# 
# listDisease <- list()
# for (modules in names(signiflist)) {
#   print(modules)
#   top10 <- signiflist[[modules]]
#   top10 <- top10[order(top10$disease, decreasing = TRUE), ]
#   lengTH <- length(top10$disease)
#   if (lengTH > 10) {
#     top10 <- top10$disease[1:10]
#   }
#   else{
#     top10 <- top10$disease
#   }
#   listDisease[[modules]] <- data.frame(top10)
# }
# 
# df_diseases <- do.call(rbind, listDisease)
# df_diseases <- data.frame(disease = na.omit(unique(df_diseases$top10)))
# 
# 
# # unique(df.long$Pathern)
# 
# Mygraph <- df.long[df.long$disease %in% df_diseases$disease &
#                      df.long$Pathern %in% colnames(PADJS), ]
# 
# Mygraph$disease <- factor(Mygraph$disease,
#                           levels = rev(unique(df_diseases$disease)))
# 
# 
# coloring <- function(x)
#   if_else(x < 0.05, "<0.05", ">0,05")
# Pvalue <- coloring(Mygraph$pValue)
# Mygraph$color <- Pvalue
# 
# Mygraph <- Mygraph[with(Mygraph, order(Pathern, color)),]
# 
# Mygraph$color <- factor(Mygraph$color,
#                         levels = unique(Mygraph$color))
# 
# Mygraph$Pathern <- gsub("padj_", "", Mygraph$Pathern)
# 
# Mygraph$Pathern <- factor(Mygraph$Pathern,
#                           levels = gsub("padj_", "", colnames(PADJS)))
# 
# 
# Mygraph %>%
#   ggplot(aes(
#     x = Pathern,
#     y = disease,
#     fill = color,
#     size = I(5)
#   )) +
#   geom_point(shape = 21, colour = "black") +
#   #scale_color_viridis_c(name = 'P value') +
#   scale_fill_manual(values = c("#5c1e3d", "#fcf6f6")) +
#   cowplot::theme_cowplot() +
#   theme(axis.line  = element_blank()) +
#   theme(axis.text.x = element_text(
#     angle = 90,
#     vjust = 0.5,
#     hjust = 1
#   )) +
#   ylab('') +
#   theme(axis.ticks = element_blank())

#--------------------------------------------------------------------------------
#Empty list

# For each disease in diseases
# Table SNPS/Enhancers overlap column pattern
# colnames()<-c("Pattern","Frequency") 
# Add to list

#Combine list into a dataframe
#Export results and Contrast with diseases with significant results
#Maybe left join?
#--------------------------------------------------------------------------------
