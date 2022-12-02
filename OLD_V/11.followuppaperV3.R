
setwd("/mnt/nocode/juliana/PID/ThemeII")
#------------------------------------------------------------------------------
#                   SETTING PARAMETERS
#------------------------------------------------------------------------------
# TRAIT OF INTEREST
trait<-"Metabolic syndrome"

# CELL TYPES OF INTEREST
celltypes<-c("M0", "M1", "M2")
#----------------------------------INPUTS--------------------------------------
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

##Subseting only fragments IDs of non-coding fragments (without genes)
NonGenic_frags <-
  Genes_mapedtonet[is.na(Genes_mapedtonet$gene_type), "fragment_ID"]

#------------------------------------------------------------------------------
#                   BLUEPRINT Active enhancers containing SNPs
#------------------------------------------------------------------------------
SNPsToEnhInterest <- readRDS("OUTPUT/SNPsToBkGdE.rds")
colnames(SNPsToEnhInterest)[1:5] <-
  c("chr", "Start", "End", "Width", "Strand")

gr_EnhSNP <- makeGRangesFromDataFrame(
  data.frame(SNPsToEnhInterest),
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
#             BLUEPRINT Active enhancers Enrichment results
#------------------------------------------------------------------------------
GWAS.ENRICH.RESULTS <- read.delim("OUTPUT/GWAS_results_0001.txt")
GWAS.ENRICH.RESULTS <-
  GWAS.ENRICH.RESULTS %>% 
  separate(disease, into = c("TRAIT", "STUDY"), sep = "_")


#----------------------------------ANALISIS--------------------------------------

#------------------------------------------------------------------------------
#   FILTERING Enh OVERLAPPING WITH TRAIT SNPS IN CELLTYPES 
#------------------------------------------------------------------------------
gr_EnhSNP[ gr_EnhSNP$DISEASE.TRAIT== trait &
             (gr_EnhSNP$M2==1 | gr_EnhSNP$M1==1 |gr_EnhSNP$M0==1),]

EnhSNP.celltypes<-
  gr_EnhSNP[,celltypes] %>% 
  data.frame() %>% 
  dplyr::select(celltypes) %>% 
  apply(2, as.numeric) %>% 
  rowSums()>0

EnhSNP.trait.celltypes<- 
  gr_EnhSNP[gr_EnhSNP$DISEASE.TRAIT== trait &
               EnhSNP.celltypes,]


length(unique(EnhSNP.trait.celltypes$SNPS))#8

#------------------------------------------------------------------------------
#   FILTERING Enh LIST ENRICHED IN TRAIT's SNPs 
#------------------------------------------------------------------------------
GWAS.ENRICH.TRAIT <-GWAS.ENRICH.RESULTS[GWAS.ENRICH.RESULTS$TRAIT==trait,]
Enriched.Enh.List<-unique(GWAS.ENRICH.TRAIT$EnhancerGroup)

EnhSNP.enriched<-EnhSNP.trait.celltypes[
  EnhSNP.trait.celltypes$Patrons %in% Enriched.Enh.List,]
# GRanges object with 42 ranges and 40 metadata columns:
# unique(EnhSNP.enriched$Patrons) %in% Enriched.Enh.List
data.frame(chr = unique(EnhSNP.enriched$CHR_ID) )

EnhSNP.enriched.summary<-
  EnhSNP.enriched %>% 
  data.frame()%>%
  group_by(SNPS) %>% 
  summarise(
    PUBMED.ID= paste(unique(PUBMEDID), sep=", "),
    CONTEXT= paste(unique(CONTEXT), sep=", "),
    CHR_ID= paste(unique(CHR_ID), sep=", "),
    Num.Enh=length(unique(regionID)),
    Enh.Activity = paste(unique(Patrons), collapse=", "),
    Enh.ID= paste(unique(regionID), collapse=", "))

table(EnhSNP.enriched$SNPS) %>% data.frame
#------------------------------------------------------------------------------
#   MAPPING "ENHANCERS FROM ENHANCERS LISTS ENRICHED ON TRAIT SNPS"
#   TO "NETWORK FRAGMENTS"
#------------------------------------------------------------------------------
Enh.frags <-
  plyranges::join_overlap_intersect(EnhSNP.enriched, gr_metadata_net_38)

## 10 ranges overlaps
Enh.frags.reduc <- unique(Enh.frags[, c("Patrons", "regionID","PUBMEDID", 
                            "CHR_ID", "CHR_POS", "SNPS", "DISEASE.TRAIT", 
                            "CONTEXT", "fragment_ID", "fragment_size", 
                            "fragment_genes")])

Enh.frags.reduc.summary <-
Enh.frags.reduc %>% 
  data.frame()%>%
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

# Only mapped the ones in "M1|M0" 
# Chromosomes 11: 9 enhancers 1 SNP window rs10838681
# Chromosomes 19: 1 enhancer 1 snp rs157582
# All in introns
# Only one pubmend  ID:  22399527
# Fragment genes: DDB2;RN7SL772P;RP11-17G12.2
# DDB2: required for UV-B tolerance and genomic integrity
# RN7SL772P: pseudogene 
# RP11-17G12.2 no info

#------------------------------------------------------------------------------
# WHICH GENES ARE LOCATED INSIDE ENHANCER FRAGMENTS?
#------------------------------------------------------------------------------
Genes.Enh.Frag<-Genes_mapedtonet$fragment_ID %in% Enh.frags.reduc$fragment_ID
unique(Genes_mapedtonet[Genes.Enh.Frag, ])
unique(Genes_mapedtonet[Genes.Enh.Frag, "gene_name"])

Genes.Enh.Frag.summary<-
  unique(Genes_mapedtonet[Genes.Enh.Frag, ]) %>% 
  data.frame%>% 
  dplyr::select(fragment_ID, gene_name,  gene_ID, gene_size, gene_type)


# NR1H3: El receptor hepático X alfa o en inglés 
# es un receptor nuclear de transcripción codificado en humanos por el gen RN1H3. 
# Este receptor está regulado por formas específicas oxidadas del colesterol, 
# los oxisteroles, y por productos que median la biosíntesis del colesterol.

# TOMM40: codes for a protein that is embedded into outer membranes of mitochondria 
# and is required for the movement of proteins into mitochondria. More precisely, 
# TOMM40 is the channel-forming subunit of a translocase of the mitochondrial outer 
# membrane (TOM) that is essential for protein transport into mitochondria.

#------------------------------------------------------------------------------
# WHICH GENES ARE INTERACTING WITH ENHANCER FRAGMENTS?
#------------------------------------------------------------------------------
##Fragments with enhancers associated to Trait
unique(Enh.frags.reduc$fragment_ID)
# Two fragments 
# 116432 (bait) NR1H3 
# 329318 (oe) TOMM40

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


# LRP4: La proteína 4 relacionada con el receptor de lipoproteínas de baja densidad, 
# también conocida como dominios 7 similares al factor de crecimiento epidérmico múltiple,
# es una proteína que en humanos está codificada por el gen LRP4.

# ACP2:
# La fosfatasa ácida lisosomal es una enzima que en humanos está codificada por el gen ACP2. 
# La fosfatasa ácida lisosomal se compone de dos subunidades, alfa y beta, y es química 
# y genéticamente distinta de la fosfatasa ácida de glóbulos rojos.
