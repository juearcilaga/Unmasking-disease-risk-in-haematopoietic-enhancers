
#------------------------------------------------------------------------------------
#                           IMPORTING PVALS
#------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------
#               Building Matrix Disease*Pattern filling with NA
#------------------------------------------------------------------------------------
EnrichmentEnhparalel_2_100 <- readRDS("/mnt/nocode/juliana/PID/ThemeII/OUTPUT/EnrichmentEnhparalel_2_100.rds")
colnames(EnrichmentEnhparalel_2_100)
# [1] "DISEASE"            "PATTERN"            "AEinPATTERN"        "TotalAEInd"        
# [5] "AEinPATTERNToSNPs"  "TotalAEIndToSNPs"   "AEinPATTERNoutSNPs" "TotalAEIndoutSNPs" 
# [9] "p.value"       

GWAS_results<-EnrichmentEnhparalel_2_100
summary(GWAS_results$p.value)

library(tidyverse)
GWAS.sign<-GWAS_results %>% filter(p.value<0.05)

summary(GWAS.sign$AEinPATTERNToSNPs)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.000   3.000   4.000   6.691   8.000 123.000 

GWAS.sign<-GWAS.sign %>% filter(AEinPATTERNToSNPs>=mean(AEinPATTERNToSNPs))
summary(GWAS.sign$AEinPATTERNToSNPs)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 7.00    8.00   11.00   13.37   15.00  123.00

GWAS.sign<-GWAS.sign %>% separate(DISEASE, into=c("DISEASE", "PUBMED.ID"), sep="_")

GWAS.collapsed<-GWAS.sign %>% group_by(DISEASE,PATTERN) %>% 
  summarise(sig.PUBMED.IDs=paste0(unique(PUBMED.ID), collapse = ", "),
            N.PUBMED.IDs= length(unique(PUBMED.ID)))
#------------------------------------------------------------------------------------
Patterns <- read.table("OUTPUT/Table_PatternsAE.txt", sep="\t")
Patterns.names <- colnames(Patterns)
Patterns <- Patterns[,c(29,1:27)]
Patterns$matrix <- apply(Patterns[,2:28], 1, paste, collapse = ",")
Patterns <- Patterns[, c("Patron", "matrix")]
colnames(Patterns)[1] <- "PATTERN"

GWAS.collapsed <- left_join(GWAS.collapsed,Patterns)
GWAS.collapsed <- GWAS.collapsed %>% 
  separate(matrix, into=Patterns.names[1:27], sep=",")

GWAS.collapsed <-data.frame(GWAS.collapsed[,1:4], 
           apply(GWAS.collapsed[,5:length(GWAS.collapsed)], 2, as.numeric))

GWAS.collapsed.disease <- 
  GWAS.collapsed %>% group_by(DISEASE) %>%
  dplyr::select(-c("sig.PUBMED.IDs","PATTERN","N.PUBMED.IDs"))%>%
  summarise(across(everything(), ~ sum(., is.na(.), 0)))
  


rownames(MAT)<-GWAS.collapsed$DISEASE
#------------------------------------------------------------------------------------
GWAS_results<-pivot_wider(Stable_GWA[,c(1,2,3)], names_from = "EnhancerGroup",values_from = "padj" )
MAT<-GWAS_results[,2:length(GWAS_results)]
rownames(MAT)<-GWAS_results$disease
MAT$sum<-rowSums(MAT)
MAT<-MAT[order(MAT$sum, decreasing = FALSE),]
#Top are the most complex diseases
#Top more signifcant
  
# ordercols<-read.delim("OUTPUT/patternclust.txt")
# ordercols$x<-factor(ordercols$x, levels= ordercols$x)
# ordercols$x[1:5]
# levels(ordercols$x)[1:5]
# 
# MATordered<-MAT[,as.character(ordercols$x)]
# colnames(MAT)[1:5]
# colnames(MATordered)

# MATbinary<-ifelse(MATordered[,2:length(MATordered)]<0.001, 1, 0)#Esta ajustado o no??
# MATbinary<-MATbinary[rowSums(MATbinary)!=0,]
# dim(MATbinary)# 798  98
# sort(rowSums(MATbinary))

MATbinary<-ifelse(MAT[,2:length(MAT)]<0.001, 1, 0)#Esta ajustado
MATbinary<-MATbinary[rowSums(MATbinary)!=0,]
dim(MATbinary)# 798  98
sort(rowSums(MATbinary))



hist(table(MAT$sum))
summary(MAT$sum)

few_patterns<-MAT[MAT$sum>90 & MAT$sum<98,]#Enriquecida en 1 a 10 patrones
dim(few_patterns)#480 100

few_patterns_binary<-ifelse(few_patterns[,2:99]<0.001, 1, 0)
few_patterns_binary<-few_patterns_binary[rowSums(few_patterns_binary)!=0,]
dim(few_patterns_binary)#[1] 379  98

pheatmap(few_patterns_binary[1:50,], border_color = "black", 
         color = c("white", "#ad2487"), 
         cluster_cols = FALSE, 
         show_colnames = FALSE, 
         fontsize_row = 5, treeheight_row = 20,scale = "none")

mostsignif<-MAT[1:50,]
mostsignif_binary<-ifelse(mostsignif[,2:length(mostsignif)-1]<0.001, 1, 0)
mostsignif_binary<-mostsignif_binary[rowSums(mostsignif_binary)!=0,]

sort(rowSums(mostsignif_binary))
dim(few_patterns_binary)#[1] 379  98

# [1] mDC|imDC|osteo|cswitmB|naiveB           
# mDC|imDC|cswitmB|naiveB                
# [3] mDC|imDC|cswitmB                       
# mDC|cswitmB                            
# [5] mDC|cswitmB|naiveB 

colnames(MATordered)





#Number of patterns with no significant enrichment per disease
#Frequency is number of diseases 
Diseases1<-rowSums(GWAS_results==1)
summary(Diseases1)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#19.0    83.0    90.0    87.3    95.0    98.0 

#Number of diseases with no significant enrichment per pattern
#Frequency is number of patterns

Pattern1<-colSums(GWAS_results==1)
summary(Pattern1)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.0   832.0   868.0   818.9   889.5   913.0 
head(GWAS_results==1)
#-PLOTS--------------------yes
-----------------------------------------------------
hist(Diseases1, breaks = 100 )
hist(Pattern1, breaks = 100)

summary(Diseases1)
summary(Pattern1)
#------------------------------------------------------------------------------------
#                            ADJUST PVALS
#------------------------------------------------------------------------------------
PadJ <- function(x) {
  padj <- p.adjust(x, method = "bonferroni")
  return(padj)
}

library(tidyr)
GWAS_long<-tidyr::pivot_longer(GWAS_results, cols=2:length(GWAS_results), names_to= "pathern", values_to="pval")

PADJS <- PadJ(GWAS_long$pval)
GWAS_long$Padj<-PADJS

newdf<-GWAS_long 
newdf$pathern<-gsub("\\|","_",newdf$pathern)
newdf$abrev<-abbreviate(newdf$disease, minlength = 5, use.classes = TRUE,
                        dot = FALSE, strict = FALSE,
                        method = "both.sides", named = TRUE)

write.table(newdf, "OUTPUT/EnrichmentEnhparalel_2_100_long.txt", sep="\t")
newdf[newdf$pval<0.05,]
min(newdf$pval)
#-PLOTS-------------------------------------------------------------------------

hist(GWAS_long$Padj, breaks = 100)
hist(GWAS_long$pval, breaks = 100)

hist(-log(GWAS_long$Padj), breaks = 100)
hist(-log(GWAS_long$pval), breaks = 100)

hist(-log(GWAS_long$Padj), breaks = 100, ylim = c(0,10))
hist(-log(GWAS_long$pval), breaks = 100,ylim = c(0,10))

summary(-log(GWAS_long$pval))


plot(-log(GWAS_long$pval),-log(GWAS_long$Padj))
plot(GWAS_long$pval, -log(GWAS_long$pval))
##########################################






