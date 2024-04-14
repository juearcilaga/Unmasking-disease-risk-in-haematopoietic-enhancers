#------------------------------------
# IMPORTING DATA FROM GWAS CATALOG DB
#------------------------------------
setwd("/mnt/nocode/juliana/PID/ThemeII")

#We will use the GWAS catalog database
library(gwascat)
data(ebicat38)
gwas_catalog <- data.frame(ebicat38, stringsAsFactors=TRUE)

citation("gwascat")
#----OPTIONAL----------------------------------------------
#How many STUDIES in the catalog?
length(unique(gwas_catalog$PUBMEDID))#2001

#How many TRAITS/DISEASEs in the catalog?
length(unique(gwas_catalog$DISEASE.TRAIT))#1320

#How many SNPS in the catalog? 
length(unique(gwas_catalog$SNPS)) #17728

#"STRONGEST.SNP.RISK.ALLELE"

#How many studies per trait? 
gwas_catalog %>% group_by(DISEASE.TRAIT)%>% summarise(num.Studies=length(unique(PUBMEDID)))
#---------------------------------------------------------
#Filter disease by N SNPs

SNP_TRAIT<-gwas_catalog%>% 
  dplyr::select(c("SNPS", "DISEASE.TRAIT")) %>% 
  unique()

#How many SNPS per TRAITS/DISEASEs in the catalog?
SNP_TRAIT<-table(SNP_TRAIT$DISEASE.TRAIT) %>% data.frame()
SNP_TRAIT<-SNP_TRAIT[order(-SNP_TRAIT$Freq),]
colnames(SNP_TRAIT)<-c("Trait", "NumSNPS")

#----OPTIONAL----------------------------------------------
summary(SNP_TRAIT$NumSNPS)
# Min.  1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.00    3.00    6.00   15.09   15.00  836.00

# ggplot(SNP_TRAIT, aes(NumSNPS)) + theme_classic()+ geom_histogram( binwidth= 10) 
#---------------------------------------------------------
#Filter traits with more than the meadian N SNPs in the catalogue
gwascatbymedianSNP <- gwas_catalog[
  gwas_catalog$DISEASE.TRAIT %in% 
    SNP_TRAIT[SNP_TRAIT$NumSNPS>5, "Trait"], ]

#--------------------------------------------------------------

#How many SNPs per study?
SNP_STUDY<-table(gwascatbymedianSNP$PUBMEDID) %>% data.frame()
SNP_STUDY<-SNP_STUDY%>% arrange(-Freq)
colnames(SNP_STUDY)<-c("PUBMED.ID", "NumSNPS")
summary(SNP_STUDY$NumSNPS)
summary(SNP_STUDY$NumSNPS)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.00    2.00    6.00   14.22   15.00  949.00 

#But One study can involve more than 1 trait
#How many SNPs per Study-Trait?

SNP_STUDY.TRAIT<-with(gwascatbymedianSNP,table(PUBMEDID,DISEASE.TRAIT)) %>% data.frame() %>% arrange(PUBMEDID,-Freq)

SNP_STUDY.TRAIT<-SNP_STUDY.TRAIT[SNP_STUDY.TRAIT$Freq>0,]
summary(SNP_STUDY.TRAIT$Freq)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.00    3.00    7.00   11.96   13.00  949.00 
SNP_STUDY.TRAIT<-SNP_STUDY.TRAIT[SNP_STUDY.TRAIT$Freq>6,]
colnames(SNP_STUDY.TRAIT)[3]<-"NumSNPS"

gwascatbymedianSNP$PUBMEDID<-factor(gwascatbymedianSNP$PUBMEDID)
left_join(SNP_STUDY.TRAIT,gwascatbymedianSNP)
#----OPTIONAL----------------------------------------------

#How many STUDIES in the catalog?
length(unique(gwascatbymedianSNP$PUBMEDID))#1480

#How many SNPS in the catalog?
length(unique(gwascatbymedianSNP $SNPS))#16439

#How many TRAITS/DISEASEs in the catalog?
length(unique(gwascatbymedianSNP$DISEASE.TRAIT))#719


#----------------------------------------------------------------------------
#Construct a data frame of studies description
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

#----OPTIONAL----------------------------------------------
summary(studies_df$sample_size)
# Min. 1st Qu.  Median  Mean  3rd Qu. Max. 
# 55    1133    2950    8874    7802  270930 
#----------------------------------------------------------------------------
pubmed_id_2K<-studies_df[studies_df$sample_size>2000,]
length(unique(pubmed_id_2K$pubmed_id))#1206

gwascat_2K<-gwas_catalog[gwas_catalog$PUBMEDID%in%pubmed_id_2K$pubmed_id,]

#----OPTIONAL----------------------------------------------
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

#----------------------------------------------------------------------------
gwascatreduc <-
  gwas_catalog[
    gwas_catalog$DISEASE.TRAIT %in% gwascatbymedianSNP$DISEASE.TRAIT &
      gwas_catalog$PUBMEDID%in%pubmed_id_2K$pubmed_id, ]

#----OPTIONAL----------------------------------------------
length(unique(gwascatreduc$SNPS))# 11060
length(unique(gwascatreduc$DISEASE.TRAIT))#518
length(unique(gwascatreduc$PUBMEDID))#985

FILTER.STATS <- 
  data.frame(
           ORIGINAL       = c(length(unique(gwas_catalog$PUBMEDID)),
                              length(unique(gwas_catalog$DISEASE.TRAIT)),
                              length(unique(gwas_catalog$SNPS))), 
           
           by.N.SNPs      = c(length(unique(gwascatbymedianSNP$PUBMEDID)),
                              length(unique(gwascatbymedianSNP$DISEASE.TRAIT)),
                              length(unique(gwascatbymedianSNP$SNPS))), 
           
           by.SAMPLE.SIZE = c(length(unique(gwascat_2K$PUBMEDID)),
                              length(unique(gwascat_2K$DISEASE.TRAIT)),
                              length(unique(gwascat_2K$SNPS))),
           
           by.BOTH.FILTERS = c(length(unique(gwascatreduc$PUBMEDID)),
                              length(unique(gwascatreduc$DISEASE.TRAIT)), 
                              length(unique(gwascatreduc$SNPS)))
             )

rownames(FILTER.STATS)<-c("PUBMEDID", "DISEASE.TRAIT", "SNPS")

FILTER.STATS$type<-rownames(FILTER.STATS)
t(FILTER.STATS)
relig_income<-FILTER.STATS

relig_income <-relig_income %>%
  pivot_longer(!type, names_to = "filter", values_to = "count")

relig_income$perc<-0

relig_income[relig_income$type == "PUBMEDID", "perc"] <-
  relig_income[relig_income$type == "PUBMEDID", "count"] /
  as.numeric(relig_income[relig_income$type == "PUBMEDID" &
                            relig_income$filter == "ORIGINAL",
                          "count"])

relig_income[relig_income$type == "DISEASE.TRAIT", "perc"] <-
  relig_income[relig_income$type == "DISEASE.TRAIT", "count"] /
  as.numeric(relig_income[relig_income$type == "DISEASE.TRAIT" &
                            relig_income$filter == "ORIGINAL",
                          "count"])

relig_income[relig_income$type == "SNPS", "perc"] <-
  relig_income[relig_income$type == "SNPS", "count"] /
  as.numeric(relig_income[relig_income$type == "SNPS" &
                            relig_income$filter == "ORIGINAL",
                          "count"])
library(forcats)
library(RColorBrewer)

ggplot(relig_income, aes(x=type,y=perc*100, fill=fct_rev(filter))) +
          geom_bar(width=0.8, color = "black", stat = "identity", position = position_dodge(width=0.9))  + 
          theme_classic()+ scale_fill_brewer(palette="Dark2")

 
#----------------------------------------------------------------------------------
##Analysing the context of the disease associated snps in original dataset
SNP.STUDY.CONTEXT<-gwas_catalog%>% 
  dplyr::select(c("SNPS", "CONTEXT")) %>% 
  unique()

#Fix the context column, sometimes there is more than one category and they are separated by semicolons
#I will consider only the first annotation

SNP.STUDY.CONTEXT$CONTEXT<-gsub("\\;.*", "",SNP.STUDY.CONTEXT$CONTEXT)
N.SNP.CONTEXT<-table(SNP.STUDY.CONTEXT$CONTEXT) %>% data.frame()
N.SNP.CONTEXT<-N.SNP.CONTEXT[order(-N.SNP.CONTEXT$Freq),]
colnames(N.SNP.CONTEXT)<-c("Context", "NumSNPS")

N.SNP.CONTEXT$Perc.NumSNPS=round((N.SNP.CONTEXT$NumSNPS/(sum(N.SNP.CONTEXT$NumSNPS)))*100, digits = 1)
N.SNP.CONTEXT$Context.Perc=paste(N.SNP.CONTEXT$Context, N.SNP.CONTEXT$Perc.NumSNPS, sep ="\n")

# Package
library(treemap)



# Plot
treemap(N.SNP.CONTEXT,
        # data
        index="Context.Perc",
        vSize="NumSNPS",
        type="index",
        
        # Main
        title="",
        palette="Dark2",
        
        # Borders:
        border.col=c("black"),             
        border.lwds=1,                         
        
        # Labels
        fontsize.labels=10,
        fontcolor.labels="white",
        fontface.labels=1,            
        bg.labels=c("transparent"),              
        align.labels=c("left", "top"),                                  
        overlap.labels=0.5,
        inflate.labels=F                      # If true, labels are bigger when rectangle is bigger.
)
treemap
#----------------------------------------------------------------------------------
##Analysing the context of the disease associated snps in filtered dataset
SNP.STUDY.CONTEXT<-gwascatreduc%>% 
  dplyr::select(c("SNPS", "CONTEXT")) %>% 
  unique()

#Fix the context column, sometimes there is more than one category and they are separated by semicolons
#I will consider only the first annotation

SNP.STUDY.CONTEXT$CONTEXT<-gsub("\\;.*", "",SNP.STUDY.CONTEXT$CONTEXT)
N.SNP.CONTEXT<-table(SNP.STUDY.CONTEXT$CONTEXT) %>% data.frame()
N.SNP.CONTEXT<-N.SNP.CONTEXT[order(-N.SNP.CONTEXT$Freq),]
colnames(N.SNP.CONTEXT)<-c("Context", "NumSNPS")


# Package
library(treemap)

# Plot
treemap(N.SNP.CONTEXT,
        # data
        index="Context",
        vSize="NumSNPS",
        type="index",
        
        # Main
        title="",
        palette="Dark2",
        
        # Borders:
        border.col=c("black"),             
        border.lwds=1,                         
        
        # Labels
        fontsize.labels=10,
        fontcolor.labels="white",
        fontface.labels=1,            
        bg.labels=c("transparent"),              
        align.labels=c("left", "top"),                                  
        overlap.labels=0.5,
        inflate.labels=T                       # If true, labels are bigger when rectangle is bigger.
)





# N.SNP.CONTEXT$perc<-(N.SNP.CONTEXT$NumSNPS/sum(N.SNP.CONTEXT$NumSNPS))*100

# ggplot(N.SNP.CONTEXT, aes(x = Context , y = perc, fill = color)) +
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

#-------------------------------------------------------------------------------------------------------------

#My snps of interest are outside genes
#gwas_catalog$CONTEXT they are semicolon separated
###############################################################################
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
  gwascatred %>% unite("Study.Trait",
                   c("DISEASE.TRAIT","PUBMEDID"),
                   remove = FALSE)

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

gr_gwas_cat_reduc<-readRDS("OUTPUT/old_OUTPUT/gr_gwas_cat_reduc.rds")

gr_gwas_cat_reduc<-gr_gwas_cat_reduc %>%
  data.frame() %>%
  group_by(DISEASE.TRAIT) %>%
  summarise(num.Studies = length(unique(PUBMEDID)),
            PUBMEDIDs = paste(unique(PUBMEDID), collapse = ", "),
            num.SNPS=paste(length(unique(SNPS)),collapse=", ")) 

write.table(gr_gwas_cat_reduc, "OUTPUT/gwas_studies_trait.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
            
