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
gr_gwas_cat_reduc<-readRDS("OUTPUT/gr_gwas_cat_reduc.rds")

# gwas_catalog<-gr_gwas_cat_reduc
# SNP.STUDY.CONTEXT<-gwas_catalog%>% 
#   dplyr::select(c("SNPS", "CONTEXT")) %>% 
#   unique()
# 
# SNP.STUDY <-
#   data.frame(
#     STUDY = paste(gwas_catalog$PUBMEDID, gwas_catalog$DISEASE.TRAIT, sep = "|"),
#     SNPS = gwas_catalog$SNPS,
#     CONTEXT=gwas_catalog$CONTEXT
#   )
# 
# N.SNP.STUDY<-table(SNP.STUDY$STUDY) %>% data.frame()
# N.SNP.STUDY<-N.SNP.STUDY[order(-N.SNP.STUDY$Freq),]
# colnames(N.SNP.STUDY)<-c("Study", "NumSNPS")
# min(N.SNP.STUDY$NumSNPS)

# #-----------------------------------------------------------------
# #               IMPORTING ANY ENHANCERS 
# #   (Present in at least two samples of any of the 31 cell types)
# #------------------------------------------------------------------
gr_CtypeChrStates <- readRDS("OUTPUT/gr_AE_list.rds")#49309
length(gr_CtypeChrStates$region)#1526184

gr_CtypeChrStates_NoY <-gr_CtypeChrStates[seqnames(gr_CtypeChrStates)!="chrY"]
# seqlevelsStyle(gr_CtypeChrStates_NoY)#"UCSC"

#------------------------------------------------------------------
#              Mapping SNPs to All enhancers
#------------------------------------------------------------------
SNPsToBkGdAE <-na.omit(data.frame(plyranges::join_overlap_intersect(
  gr_CtypeChrStates_NoY,gr_gwas_cat_reduc)))

seqnames(gr_gwas_cat_reduc)###23
seqnames(gr_CtypeChrStates_NoY)###x y Y remove Y 

dim(SNPsToBkGdAE)#165714     44

#------------------------------------------------------------------
# NonORF<-c("intergenic","","ncRNA")
# ORFNoncode<-c("cds-synon", "intron","UTR-5", "UTR-3", "splice-3",
#               "nearGene-3",  "nearGene-5",  "splice-5")

unwanted<-c("missense", "STOP-GAIN","frameshift")

SNPsToBkGdAE<-SNPsToBkGdAE[grepl(paste(unwanted,collapse="|"),
                                 SNPsToBkGdAE$CONTEXT)==FALSE,]

summary(SNPsToBkGdAE$width)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.0   200.0   200.0   196.2   200.0   200.0 

SNPsToBkGdAE<-SNPsToBkGdAE[SNPsToBkGdAE$width>=200,]
dim(SNPsToBkGdAE)#150048     44

saveRDS(SNPsToBkGdAE, "OUTPUT/SNPsToBkGdE.rds")
#------------------------------------------------------------------
#              Counting SNPs of interest
#------------------------------------------------------------------
# SNP.ENh <-
#   data.frame(
#     STUDY = paste(SNPsToBkGdAE$PUBMEDID, SNPsToBkGdAE$DISEASE.TRAIT, sep = "|"),
#     SNPS = SNPsToBkGdAE$SNPS,
#     CONTEXT=SNPsToBkGdAE$CONTEXT,
#     PATRON=SNPsToBkGdAE$Patrons,
#     ENHANCER=SNPsToBkGdAE$regionID)
# 
# SNP.ENh.context<-unique(SNP.ENh[,c("SNPS","CONTEXT")])
# SNP.ENh.context$CONTEXT<-gsub("\\;.*", "",SNP.ENh.context$CONTEXT)
# N.SNPE.STUDY<-table(SNP.ENh.context$CONTEXT) %>% data.frame()
# N.SNPE.STUDY<-N.SNPE.STUDY[order(-N.SNPE.STUDY$Freq),]
# colnames(N.SNPE.STUDY)<-c("Context", "NumSNPS")
# N.SNPE.STUDY$perc<-(N.SNPE.STUDY$NumSNPS/sum(N.SNPE.STUDY$NumSNPS))*100
# 
# colnames(N.SNPE.STUDY)<-c("Context", "NumSNPSEnh", "perc")
# 
# N.SNP.CONTEXT<-gwas_catalog%>% 
#   dplyr::select(c("SNPS", "CONTEXT")) %>% 
#   unique()
# 
# #Fix the context column, sometimes there is more than one category and they are separated by semicolons
# #I will consider only the first annotation
# 
# SNP.STUDY.CONTEXT$CONTEXT<-gsub("\\;.*", "",SNP.STUDY.CONTEXT$CONTEXT)
# N.SNP.CONTEXT<-table(SNP.STUDY.CONTEXT$CONTEXT) %>% data.frame()
# N.SNP.CONTEXT<-N.SNP.CONTEXT[order(-N.SNP.CONTEXT$Freq),]
# colnames(N.SNP.CONTEXT)<-c("Context", "NumSNPS")
# 
# 
# 
# PLOT<-left_join(N.SNPE.STUDY,N.SNP.CONTEXT, by="Context")
# PLOT<-PLOT[,c("Context", "NumSNPS", "NumSNPSEnh")]
# PLOT$perc<-PLOT$NumSNPSEnh/PLOT$NumSNPS*100
# 
# 
# ggplot(PLOT, aes(x = Context , y = perc)) +
#   geom_bar(color = "black", stat = "identity") +
#   theme_classic() +
#   theme(axis.text.x = element_text(
#     size = 12,
#     angle = 90,
#     vjust = 0.5,
#     hjust = 1
#   )) +
#   theme(axis.text.y = element_text(size = 10)) +
#   geom_text(
#     aes(label = NumSNPS, sep="\n") ,vjust= -0.2,
#     size = 3.5,
#     fontface = "bold",
#     family = "Fira Sans"
#   ) +
#   scale_fill_manual(values = c("#E69F00", "darkorchid", "#56B4E9", "#999999"))
# 
# 
# 
# N.Enh.STUDY<-table(unique(SNP.ENRICH[,c("STUDY","ENHANCER")])$STUDY) %>% data.frame()
# N.Enh.STUDY<-N.Enh.STUDY[order(-N.Enh.STUDY$Freq),]
# colnames(N.Enh.STUDY)<-c("Study", "NumEnhE")
# 
# Perc.SNPSE<-left_join(N.SNPE.STUDY,N.SNP.STUDY)
# Perc.SNPSE$Perc<-Perc.SNPSE$NumSNPSE/Perc.SNPSE$NumSNPS
# Perc.SNPSE<-Perc.SNPSE[order(-Perc.SNPSE$Perc),]
# Perc.SNPSE$PUBMED<-gsub("\\|.*","",Perc.SNPSE$Study)
# 
# Enriched<-read.delim("OUTPUT/EnrichmentEnhparalel_2_100_long.txt", sep="\t")
# Enriched<-Enriched[Enriched$Padj<0.01,]
# Enriched<-Enriched[order(Enriched$Padj),]
# 
# NPatrones<-table(Enriched$pathern) %>% data.frame()
# NPatrones<-NPatrones[order(-NPatrones$Freq),]      
# colnames(NPatrones)<-c("Pattern", "Enriched in N studies")
# write.table(NPatrones, "OUTPUT/NumEnrichedStudysbyPatterns.txt",sep="\t")
# 
# Enriched[Enriched$pathern=="mDC_imDC", ]
# 
# List.Enriched.Study<-Perc.SNPSE[Perc.SNPSE$PUBMED %in% unique(gsub(".*_","",Enriched$disease) ),]
# dim(List.Enriched.Study)#594 studies
# List.Enriched.Study<-List.Enriched.Study[,c("Study", "NumSNPS", "NumSNPSE", "Perc")]
# write.table(List.Enriched.Study, "OUTPUT/NSNPs.Enriched.Study.txt", sep="\t")

#=----------------------------------------------------------------------

