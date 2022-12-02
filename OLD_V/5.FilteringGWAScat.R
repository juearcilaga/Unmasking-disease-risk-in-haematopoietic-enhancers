#------------------------------------
# IMPORTING DATA FROM GWAS CATALOG DB
#------------------------------------
setwd("/mnt/nocode/juliana/PID/ThemeII")
#We will use the GWAS catalog database
library(gwascat)
data(ebicat38)
gwas_catalog <- data.frame(ebicat38, stringsAsFactors=TRUE)
# summary(as.factor(gwas_catalog$SNPS))
#How many STUDIES in the catalog?

length(unique(gwas_catalog$PUBMEDID))#2001

#How many TRAITS/DISEASEs in the catalog?
length(unique(gwas_catalog$DISEASE.TRAIT))#1320

#How many SNPS in the catalog? 
length(unique(gwas_catalog$SNPS)) #17728

##################################################################################################################
#Filter disease by N SNPs
SNP_TRAIT<-gwas_catalog%>% 
  dplyr::select(c("SNPS", "DISEASE.TRAIT")) %>% 
  unique()

#How many SNPS per TRAITS/DISEASEs in the catalog?
SNP_TRAIT<-table(SNP_TRAIT$DISEASE.TRAIT) %>% data.frame()
SNP_TRAIT<-SNP_TRAIT[order(-SNP_TRAIT$Freq),]
colnames(SNP_TRAIT)<-c("Trait", "NumSNPS")

summary(SNP_TRAIT$NumSNPS)
# Min.  1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.00    3.00    6.00   15.09   15.00  836.00

# ggplot(SNP_TRAIT, aes(NumSNPS)) + theme_classic()+ geom_histogram( binwidth= 10) 

gwascatbymedianSNP <- gwas_catalog[
  gwas_catalog$DISEASE.TRAIT %in% 
    SNP_TRAIT[SNP_TRAIT$NumSNPS>5, "Trait"], ]


#How many SNPS in the catalog?
length(unique(gwascatbymedianSNP $SNPS))#16439

#How many STUDIES in the catalog?
length(unique(gwascatbymedianSNP$PUBMEDID))#1480

#How many TRAITS/DISEASEs in the catalog?
length(unique(gwascatbymedianSNP$DISEASE.TRAIT))#719
###################################################################
#Construct a data frame of studies
studies = gwas_catalog[, c("PUBMEDID", "INITIAL.SAMPLE.DESCRIPTION")] %>% unique()
colnames(studies) = c("pubmed_id", "description")

#Extract sample size from the gwas catalog
sample_sizes = lapply(as.list(studies$description), function(x) {
  stringr::str_extract_all(x, "(\\d+,\\d+)|(\\d+)") %>%
    unlist() %>%
    stringr::str_replace(",", "") %>%
    as.numeric() %>%
    sum()
})

studies_df = dplyr::mutate(studies, sample_size = unlist(sample_sizes)) 
summary(studies_df$sample_size)
# Min. 1st Qu.  Median  Mean  3rd Qu. Max. 
# 55    1133    2950    8874    7802  270930 

pubmed_id_2K<-studies_df[studies_df$sample_size>2000,]
length(unique(pubmed_id_2K$pubmed_id))#1206

gwascat_2K<-gwas_catalog[gwas_catalog$PUBMEDID%in%pubmed_id_2K$pubmed_id,]
length(unique(gwascat_2K$SNPS))#11689
length(unique(gwascat_2K$DISEASE.TRAIT))#810

# #colnames(gwas_cat_10K)
# ggplot(studies_df, aes(sample_size)) +
#   theme_classic()+
#   geom_histogram(binwidth= 100) +
#   scale_x_continuous(breaks = seq(10000,300000, by = 50000),limits = (x=c(10000,300000)))
# 
# ggplot(studies_df, aes(sample_size)) +
#   theme_classic()+
#   geom_histogram(binwidth= 1000) +
#   ylim(c(0,45))+
#   scale_x_continuous(limits = (x=c(0,3E5)), 
#                      labels = unit_format(unit = "e^5", scale = 1 / 1e+5, digits = 0))
##############################################################################
gwascatreduc <-
  gwas_catalog[
    gwas_catalog$DISEASE.TRAIT %in% gwascatbymedianSNP$DISEASE.TRAIT &
      gwas_catalog$PUBMEDID%in%pubmed_id_2K$pubmed_id, ]

length(unique(gwascatreduc$SNPS))# 11060
length(unique(gwascatreduc$DISEASE.TRAIT))#518
length(unique(gwascatreduc$PUBMEDID))#985

###############################################################################
#My snps of interest are outside genes
#gwas_catalog$CONTEXT they are semicolon separated
###############################################################################
contexto <- 
  lapply(
    as.list(gwas_catalog$CONTEXT),
    function(x) { 
      stringr::str_split(x, ";") %>% unlist()
      }) %>% 
  unlist() %>% 
  unique()

# [1] "intron"     "intergenic" "cds-synon"  "missense"   "ncRNA"      ""           "nearGene-5"
# [8] "UTR-5"      "UTR-3"      "splice-3"   "nearGene-3" "STOP-GAIN"  "splice-5"   "frameshift"

#"missense", "STOP-GAIN","frameshift" 

NonORF<-c("intergenic","","ncRNA")
ORFNoncode<-c("cds-synon", "intron","UTR-5", "UTR-3", "splice-3", "nearGene-3",  "nearGene-5",  "splice-5")

gwascatreduc2<-gwascatreduc[gwascatreduc$CONTEXT%in%NonORF,]
length(unique(gwascatreduc2$SNPS))# 11060#4700
length(unique(gwascatreduc2$DISEASE.TRAIT))#518#502
length(unique(gwascatreduc2$PUBMEDID))#985#809

gwascatreduc2<-gwascatreduc[gwascatreduc$CONTEXT%in%ORFNoncode,]
length(unique(gwascatreduc2$SNPS))# 11060#5617
length(unique(gwascatreduc2$DISEASE.TRAIT))#518#510
length(unique(gwascatreduc2$PUBMEDID))#985#894


#Extract snps form catalog
gwascatred= gwascatreduc[,c(
  "seqnames",
  "start",
  "PUBMEDID",
  "CHR_ID",
  "CHR_POS",
  "SNPS",
  "DISEASE.TRAIT",
  "CONTEXT")]

colnames(gwascatred)<-c(
  "Seqnames",
  "Start",
  "PUBMEDID",
  "CHR_ID",
  "CHR_POS",
  "SNPS",
  "DISEASE.TRAIT",
  "CONTEXT")

gwascatred$chromosome <-
  paste("chr", gwascatred$CHR_ID, sep = "")

gwascatred$left <- gwascatred$CHR_POS-5000

gwascatred$right <-gwascatred$CHR_POS+5000


gwascatred<-
  gwascatred%>% unite("Study.Trait",
                   c("DISEASE.TRAIT","PUBMEDID"),
                   remove = FALSE)


annot_disease<-read.delim("../SecondYearAssesment/RESULTADOS/disease_mapping_to_attributes_DISGENET.txt")
annot_disease_target <- annot_disease[grepl("bowel|infectious|inflamm|immune|disease|syndrome",
                                            annot_disease$name) |
                                        annot_disease$type %in% c("disease", "group"), ]


gwascatred$is.disease<-"FALSE"
gwascatred[gwascatred$DISEASE.TRAIT%in%annot_disease_target$name,]$is.disease<-"TRUE"

table(gwascatred$is.disease)
# #FALSE  TRUE 
# 13838  1087 

#unique(gwascatred$Seqnames)
gwascatred$chromosome[gwascatred$chromosome=="chr23"]<-"chrX"

gr_gwas_cat_reduc <- makeGRangesFromDataFrame(
  gwascatred,
  keep.extra.columns = TRUE,
  ignore.strand = TRUE,
  seqinfo = NULL,
  seqnames.field ="chromosome",
  start.field = "left",
  end.field = "right",
  strand.field = NULL,
  starts.in.df.are.0based = FALSE
)

seqlevels(gr_gwas_cat_reduc)

saveRDS(gr_gwas_cat_reduc,
        "OUTPUT/gr_gwas_cat_reduc.rds")