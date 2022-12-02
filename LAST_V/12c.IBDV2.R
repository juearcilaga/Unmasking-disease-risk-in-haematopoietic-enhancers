
#------------------------------------------------------------------------------------
#                           IMPORTING PVALS
#------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------
#               Building Matrix Disease*Pattern filling with NA
#------------------------------------------------------------------------------------

GWAS_filename<-read.delim("OUTPUT/GWAS_results_0001.txt")
colnames(GWAS_filename)
#[1] "disease"           "EnhancerGroup"     "padj"              "AEinPATTERN"       "AEinPATTERNToSNPs"   

library(tidyverse)

GWAS.sign<-GWAS_filename

GWAS.sign<-GWAS.sign %>% filter(AEinPATTERNToSNPs>=mean(AEinPATTERNToSNPs))
summary(GWAS.sign$AEinPATTERNToSNPs)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 14.00   16.00   20.00   23.23   25.25  102.00 

GWAS.sign<-GWAS.sign %>% separate(disease, into=c("DISEASE", "PUBMED.ID"), sep="_")
GWAS.sign<-GWAS.sign %>% tidyr::unite("TEST.SIG", c(DISEASE, PUBMED.ID, EnhancerGroup), sep =", ", remove=FALSE)

GWAS.collapsed<-GWAS.sign %>% group_by(DISEASE,EnhancerGroup) %>% 
  summarise(sig.PUBMED.IDs=paste0(unique(PUBMED.ID), collapse = ", "),
            Nsig.PUBMED.IDs= length(unique(PUBMED.ID)))

DISEASE.traits <- c(
"Ulcerative colitis",
"Crohn's disease",
"Pyoderma gangrenosum in inflammatory bowel disease"
)


# DISEASE.traits <- c(
#   "Ulcerative colitis")

DISEASE.enrich <- GWAS.sign[GWAS.sign$DISEASE %in% DISEASE.traits,]
DISEASE.enrich %>% arrange(DISEASE, EnhancerGroup)
DISEASE.enrich$EnhancerGroup<-
  factor(
    DISEASE.enrich$EnhancerGroup,
    levels =
      c(
        "effmemCD8T",
        "germcenterB|naiveB",
        "naiveB",
        "mDC",
        "imDC|mDC",
        "mneut|segneutBM",
        "eryth",
        "endprog|endumbprolif",
        "endprog|endumbprolif|mesenBM"
      )
  )

DISEASE.enrich<-with(DISEASE.enrich, DISEASE.enrich[order(DISEASE,EnhancerGroup),])
library(ggpubr)
ggtexttable(DISEASE.enrich[, c(2,3,4,6,7,5)], rows = NULL)
#------------------------------------------------------------------------------
#                   BLUEPRINT Active enhancers containing SNPs
#------------------------------------------------------------------------------
EnhSNPsig <- readRDS("OUTPUT/SNPsToBkGdE.rds")
colnames(EnhSNPsig)[1:5] <-
  c("chr", "Start", "End", "Width", "Strand")

library(GenomicRanges)
gr_EnhSNP <- makeGRangesFromDataFrame(
  data.frame(EnhSNPsig),
  keep.extra.columns = TRUE,
  ignore.strand = TRUE,
  seqinfo = NULL,
  seqnames.field = "chr",
  start.field = "Start",
  end.field = "End",
  strand.field = "Strand",
  starts.in.df.are.0based = FALSE
)
#------------------------------------------------------------------------------
#   FILTERING Enh OVERLAPPING WITH TRAIT SNPS IN CELLTYPES 
#   FILTERING Enh LIST ENRICHED IN TRAIT's SNPs 
#------------------------------------------------------------------------------
EnhSNPsig<-gr_EnhSNP %>% data.frame()%>% tidyr::unite("TEST.SIG", c(DISEASE.TRAIT, PUBMEDID, Patrons), sep =", ", remove=FALSE)
EnhSNPsig<-EnhSNPsig[EnhSNPsig$TEST.SIG%in% DISEASE.enrich$TEST.SIG,]

EnhSNP.enriched.summary<-EnhSNPsig %>% 
  data.frame()%>%
  group_by(TEST.SIG) %>% 
  summarise(
    CHR_ID=paste(unique(CHR_ID), collapse=", "),
    CONTEXT= paste(unique(CONTEXT), collapse=", "),
    Num.SNPS=length(unique(SNPS)),
    Num.Enh=length(unique(regionID))) %>% 
  arrange(-Num.Enh)


ggtexttable(EnhSNP.enriched.summary, rows = NULL)

colnames(EnhSNPsig)[1:5] <-c("chr", "Start", "End", "Width", "Strand")

SNPsToEnhsig <- 
  EnhSNPsig[,c("chr", "Start", "End", "Width","regionID", 
                       "Patrons","PUBMEDID", "CHR_ID", "CHR_POS", 
                       "SNPS", "DISEASE.TRAIT", "CONTEXT")]
SNPsToEnhsig.summary<-
  SNPsToEnhsig %>%
  group_by(SNPS) %>%
  summarise(
    DISEASE.TRAIT = paste(unique(DISEASE.TRAIT), collapse = ", "),
    chr = paste(unique(chr), collapse = ", "),
    NEnh = length(unique(regionID)),
    CHR_POS= unique(CHR_POS),
    startSNPreg= min(Start),
    endSNPreg= max(End),
    EnhGroup = paste(unique(Patrons), collapse = ", ")
  )

SNPsToEnhsig.summary <-
  SNPsToEnhsig.summary[order(SNPsToEnhsig.summary$chr, decreasing = TRUE),]

SNPsToEnhsig.summary$size<-SNPsToEnhsig.summary$endSNPreg-SNPsToEnhsig.summary$startSNPreg


#Only those with more than 10 enh

SNPsToEnhsig.summary<-SNPsToEnhsig.summary[SNPsToEnhsig.summary$NEnh>10,]
# rs11676348 14 Enhs in the  chromosome 2  Ilcerative colitis
# Identification of a common variant with potential pleiotropic 
# effect on risk of inflammatory bowel disease and colorectal cancer

write.table(SNPsToEnhsig.summary,"OUTPUT/SNPsToEnhsigIBDsummary.txt", sep="\t")

#-----------------------------------------------------------------------------------------------------
# 
# SNPsToEnhsig.summary.bychr<-
#   SNPsToEnhsig.summary %>% group_by(chr) %>%
#   summarise(
#     startReg = min(start),
#     endReg = max(end),
#     sizeReg = endReg- startReg
#   ) 

# if siReg > 100, split region

SNPsToEnhsig.summary.bychrReg<-SNPsToEnhsig.summary
SNPsToEnhsig.summary.bychrReg$Regstart<-SNPsToEnhsig.summary.bychrReg$CHR_POS-10000
SNPsToEnhsig.summary.bychrReg$Regend<-SNPsToEnhsig.summary.bychrReg$CHR_POS+10000


gr_SNPsToEnhsig.summary.bychrReg<-makeGRangesFromDataFrame(
  data.frame(SNPsToEnhsig.summary.bychrReg),
  keep.extra.columns = TRUE,
  ignore.strand = TRUE,
  seqinfo = NULL,
  seqnames.field = "chr",
  start.field = "Regstart",
  end.field = "Regend",
  strand.field = "Strand",
  starts.in.df.are.0based = FALSE
)


gr_SNPsToEnhsig.summary.bychrReg<-plyranges::join_overlap_self(gr_SNPsToEnhsig.summary.bychrReg,  minoverlap=10000) %>% data.frame

SNPsToEnhsig.summary.bychr<-
  transform(gr_SNPsToEnhsig.summary.bychrReg, ID = as.numeric(factor(start)))%>% 
  group_by(ID)%>%
  summarise(
    chr=unique(seqnames),     
    start20k=unique(start),       
    end20k=unique(end),
    SNPS=paste(unique(SNPS), collapse=","))


write.table(SNPsToEnhsig.summary.bychr,"OUTPUT/SNPsToEnhsigIBD.summary.bychr.txt", sep="\t")
#------------------------------------------------------------------------------------------
##Anotate the region for plotting   
library(plyranges)
gene_Anotations<-readRDS("extData/gr_gene_Anotations38_lcincluded.rds")
gene_Anotations<-gene_Anotations[grepl("exon:|CDS:|stop_codon:|UTR",gene_Anotations$gene_ID)==FALSE,]
gene_Anotations<-gene_Anotations %>% select(-gene_ID)%>%  unique()


SNPS.GENES<-list()
for(chr in 1:length(SNPsToEnhsig.summary.bychr$chr)) {
  CHR <- SNPsToEnhsig.summary.bychr[chr, ]
  gr.CHR <- GRanges(seqnames = CHR$chr ,
                    ranges = IRanges(start = CHR$start20k, end = CHR$end20k ))
  CHR.GENES <-
    na.omit(data.frame(
      plyranges::join_overlap_intersect(gene_Anotations, gr.CHR)))
  SNPS.GENES[[CHR$chr]] <- CHR.GENES
}

SNPS.GENES.DF<-do.call(rbind, SNPS.GENES)
SNPS.GENES.DF<-SNPS.GENES.DF%>%group_by(gene_name) %>% 
  filter(gene_size==max(gene_size))

write.table(SNPS.GENES.DF,"OUTPUT/IBD.GENES.DF.txt", sep="\t")

SNPsToEnhsigGene.summary<-
  SNPsToEnhsig.summary.bychrReg %>% group_by(chr)%>% 
  summarise(SNPS=paste(unique(SNPS), collapse=", "),
            DISEASE.TRAIT=paste(unique(DISEASE.TRAIT), collapse=", "),
            EnhGroup=paste(unique(EnhGroup), collapse=", ")) %>% 
  data.frame() %>% arrange(DISEASE.TRAIT)

ggtexttable(SNPsToEnhsigGene.summary, rows = NULL)
#------------------------------------------------------------------------------
#                           NETWORK INTERACTIONS
#   File with the Chicago results of Javierre et al 2016
#------------------------------------------------------------------------------
#/home/maninder/Documents/PhD_Project/Immune_Cells/Hi-C_Data/EDIT4.txt
PCHIC.Intractions <-
  read.table(
    "/mnt/nocode/data_public/pchic_interaction_networks/merged_samples_12Apr2015_full.txt",
    header = TRUE,
    skip = 5
  )
#------------------------------------------------------------------------------
#                           NETWORK FRAGMENTS 
#------------------------------------------------------------------------------
# ggranges object with the coordinates of the Biola Javierre 2016 network 
# fragments and their degree and other statistics in M0,M1,M2

gr_metadata_net_38 <-
  readRDS("/mnt/nocode/juliana/PID/integration/RDS/gr_metadata_net_38.rds")

#------------------------------------------------------------------------------
#                           NETWORK FRAGMENTS GENE ANNOTATION
#------------------------------------------------------------------------------
##Genes in the fragment's of the network
##Corregir esto porque no coinciden las coordenadas
Genes_mapedtonet <-
  readRDS("/mnt/nocode/juliana/PID/integration/RDS/Genes_mapedtonet.rds")

#------------------------------------------------------------------------------
#   MAPPING "ENHANCERS FROM ENHANCERS LISTS ENRICHED ON TRAIT SNPS"
#   TO "NETWORK FRAGMENTS"
#------------------------------------------------------------------------------
# GRanges object with 142 ranges and 40 metadata columns:
# unique(EnhSNP.enriched$Patrons) %in% Enriched.Enh.List

gr_EnhSNPsig<-gr_EnhSNP[(elementMetadata(gr_EnhSNP)[, "regionID"] %in% EnhSNPsig$regionID)]

Enh.frags <-
  plyranges::join_overlap_intersect(gr_EnhSNPsig, gr_metadata_net_38)


Enh.frags.reduc <- unique(Enh.frags[, c("Patrons", "regionID","PUBMEDID", 
                                        "CHR_ID", "CHR_POS", "SNPS", "DISEASE.TRAIT", 
                                        "CONTEXT", "fragment_ID", "fragment_size", 
                                        "fragment_genes")])

## GRanges object with 51 ranges and 11 metadata columns

Enh.frags.reduc.summary <-
  Enh.frags.reduc %>% 
  data.frame()%>%
  filter(width==200)%>%
  group_by(SNPS) %>% 
  summarise(
    PUBMED.ID= paste(unique(PUBMEDID), sep=", "),
    CONTEXT= paste(unique(CONTEXT), sep=", "),
    CHR_ID= paste(unique(CHR_ID), sep=", "),
    Num.Enh=length(unique(regionID)),
    Enh.Activity = paste(unique(Patrons), collapse=", "),
    fragment_ID= paste(unique(fragment_ID), collapse=", "),
    fragment_size= paste(unique(fragment_size), collapse=", "),
    Enh.ID= paste(unique(regionID), collapse=", "),
  )


Enh.frags.reduc %>% 
  data.frame()%>%
  filter(width==200)%>%
  group_by(DISEASE.TRAIT) %>% 
  summarise(
    CHR_ID=paste(unique(CHR_ID), collapse=", "),
    CONTEXT= paste(unique(CONTEXT), collapse=", "),
    PUBMEDID= paste(unique(PUBMEDID), collapse=", "),
    Enh.Activity= paste(unique(Patrons), collapse=", "),
    Num.SNPS=length(unique(SNPS)),
    Num.Enh=length(unique(regionID)),
    fragment_ID= paste(unique(fragment_ID), collapse=", "),
    fragment_size= paste(unique(fragment_size), collapse=", "),) %>% 
  arrange(-Num.Enh)

# Only two frags one for each celltype net eryth
# Chromosomes 2: 1 enhancers 1 SNP window rs11676348
# Chromosomes 7: 3  enhancer  1 snp  rs4730273
# All in introns
# Only one pubmend  ID:  21297633 ,19122664


# Fragment genes: AC098820.2;SMARCAL1
# no info

#------------------------------------------------------------------------------
# WHICH GENES ARE LOCATED INSIDE ENHANCER FRAGMENTS?
#------------------------------------------------------------------------------
Genes.Enh.Frag<-Genes_mapedtonet$fragment_ID %in% Enh.frags.reduc$fragment_ID
unique(Genes_mapedtonet[Genes.Enh.Frag, ])
data.frame(unique(Genes_mapedtonet[Genes.Enh.Frag, "gene_name"]))

# chr16:NOD2
# chr16:SBNO2
# chr2:DNMT3A
# chr2:CXCR2: CXCR2 knockout mice are protected against DSS-colitis-induced acute kidney injury and inflammation
# chr4:BANK1
# chr6:LTA :Association of TNF-alpha/LTA polymorphisms with Crohn's disease in Koreans
# chr6:TAGAP-AS1:TAGAP instructs Th17 differentiation by bridging Dectin activation to EPHB2 signaling in innate antifungal response
# chr7:GNA12:IBD Candidate Genes and Intestinal Barrier Regulation, IBD risk loci are enriched in multigenic regulatory modules encompassing putative causative genes
# chr7:IRF5:A two-marker haplotype in the IRF5 gene is associated with inflammatory bowel disease in a North American cohort

##All are associated

Genes.Enh.Frag.summary<-
  unique(Genes_mapedtonet[Genes.Enh.Frag, ]) %>% 
  data.frame%>% 
  dplyr::select(gene_name,  fragment_ID, gene_ID, gene_size, gene_type)


# Chr2 CXCR2: CD182; IL8R2; IL8RA; IL8RB; CMKAR2; WHIMS2; CDw128b
# motif chemokine receptor 2 
# Chr7 NA

----------------------------------------------------------------------
# WHICH GENES ARE INTERACTING WITH ENHANCER FRAGMENTS?
#------------------------------------------------------------------------------
##Fragments with enhancers associated to Trait
unique(Enh.frags.reduc$fragment_ID)
# Two fragments 
# 396503 #bait
# 692374
Enh.frags.Interactions <-
  PCHIC.Intractions[PCHIC.Intractions$baitID %in% Enh.frags.reduc$fragment_ID,]

# Three interactions (oe) in Chr 11
# 116345
# 116349  >5 in M1 & M2
# 116430    >5 in M0, M1, M2, Neutrophils, Megakaryocytes, Endothelial_precursors, Erythroblasts

###Remember thta chicago score can be less than 5 in some ofthe celltypes
unique(Genes_mapedtonet[Genes_mapedtonet$fragment_ID%in%Enh.frags.Interactions$oeID,])

unique(Genes_mapedtonet[Genes_mapedtonet$fragment_ID%in%Enh.frags.Interactions$oeID, "gene_name"])


Enh.frags.Interactions.summary<-
  unique(
    Genes_mapedtonet[Genes_mapedtonet$fragment_ID %in% Enh.frags.Interactions$oeID,]
  )%>% 
  data.frame%>% 
  dplyr::select(fragment_ID, gene_name,  gene_ID, gene_size, gene_type)

unique(Enh.frags.Interactions.summary$gene_name)#

# MAP4K4  : MAP4K4 negatively regulates CD8 T cell–mediated antitumor and antiviral immunity
# LINC01127 
# IL1RL2:Inhibiting Interleukin 36 Receptor Signaling Reduces Fibrosis in Mice With Chronic Intestinal Inflammation
# NBEAL1   
# ERBB4:ErbB4 signaling stimulates pro-inflammatory macrophage apoptosis and limits colonic inflammation 
# ABCA12:Transcriptional Signatures That Define Ulcerative Colitis in Remission 
# RPL37A    <NA>      
# DIRC3    
# TNS1      
# ARPC2    something 
# LINC02832 INPP5D    LFNG      BRAT1     IQCE     
# TTYH3     AMZ1      
# GNA12  


#-----------------------------------------------------------------------------------
# PLOT FOR PAPER
#------------------------------------------------------------------------------------
Patterns <- read.table("OUTPUT/Table_PatternsAE.txt", sep="\t")
Patterns.names <- colnames(Patterns)
Patterns <- Patterns[,c(29,1:27)]
Patterns$matrix <- apply(Patterns[,2:28], 1, paste, collapse = ",")
Patterns <- Patterns[, c("Patron", "matrix")]
colnames(Patterns)[1] <- "EnhancerGroup"

DISEASE.enrich <- left_join(DISEASE.enrich,Patterns)

annot<-DISEASE.enrich[, c("DISEASE", "EnhancerGroup", "matrix")]
annot$values<-1

annot <-pivot_wider(annot, names_from = DISEASE, values_from = values)


annotlefth<-annot[, 1:2]
annotlefth <-annotlefth %>% 
  separate(matrix, into=Patterns.names[1:27], sep=",")

annotlefth <-data.frame(annotlefth[,1], 
           apply(annotlefth[,2:length(annotlefth)], 2, as.numeric))

rownames(annotlefth)<-annotlefth$EnhancerGroup

numericAnnotLeft<-annotlefth[,2:length(annotlefth)]

numericAnnotLeft<-numericAnnotLeft[,colSums(numericAnnotLeft)>0]

numericAnnotLeft <- numericAnnotLeft[c(
  "effmemCD8T",
  "germcenterB|naiveB",
  "naiveB",
  "mDC",
  "imDC|mDC",
  "mneut|segneutBM",
  "eryth",
  "endprog|endumbprolif",
  "endprog|endumbprolif|mesenBM"
), c(
  "effmemCD8T",
  "germcenterB",
  "naiveB",
  "mDC",
  "imDC",
  "segneutBM",
  "mneut",
  "eryth",
  "mesenBM",
  "endumbprolif",
  "endprog"
)]


library(pheatmap)
pheatmap(
  numericAnnotLeft,
  color = c("#642194", "yellow"),
  border_color = "black",
  scale = "none",
  cellwidth = 10 ,
  cellheight = 10,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  treeheight_col = 50,
  treeheight_row = 25,
  show_colnames = TRUE,
  fontsize_col = 10,
  fontsize_row = 10,
)


annotrigth<-annot[, c(1, 3:5)]
colnames(annotrigth)<- c("Pattern", "crohn disease", 
                         "pyoderma gangrenosum", "ulcerative colitis")
annotrigth <- data.frame(annotrigth[, 1],
                         apply(annotrigth[, 2:length(annotrigth)], 2, as.numeric))
rownames(annotrigth)<-annotrigth$Pattern
annotrigth[is.na(annotrigth)]<-0



annotrigth <- annotrigth[c(
  "effmemCD8T",
  "germcenterB|naiveB",
  "naiveB",
  "mDC",
  "imDC|mDC",
  "mneut|segneutBM",
  "eryth",
  "endprog|endumbprolif",
  "endprog|endumbprolif|mesenBM"
), c(2, 4, 3)]

pheatmap(
  annotrigth,
  color = c("#642194", "yellow"),
  cellwidth = 10 ,
  cellheight = 10,
  border_color = "black",
  scale = "none",
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  treeheight_col = 50,
  treeheight_row = 25,
  show_colnames = TRUE,
  fontsize_col = 10,
  fontsize_row = 10,
)
