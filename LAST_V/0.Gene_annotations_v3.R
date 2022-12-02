#Step 3: Organising gene annotations 
setwd("/mnt/nocode/juliana/PID/ThemeII/extData/")
library(rtracklayer)
#Import gene annotation table downloaded from ENCODE release 38
GENE_CODE_38<-readGFF("gencode.v38.basic.annotation.gff3", version=3,  columns=NULL, tags=NULL, filter=NULL, nrows=-1,raw_data=FALSE)

# GENES<-GENE_CODE_38[GENE_CODE_38$gene_type =="protein_coding" & 
#                       GENE_CODE_38$type =="gene",]

# annotation_genes<-GENES
annotation_genes<-GENE_CODE_38
dim(annotation_genes)#[1] 20082    35
#Using only annotations with gene_name
genes<-annotation_genes[annotation_genes$gene_name!="",]
dim(genes)#[1] 20082    35
genes<-genes[,c(1,4:5,9,11,12)]
saveRDS(genes, "gene_Anotations38_lcincluded.rds")
###############################################################
#Creating a granges object for gene information

# genes<-readRDS("gene_Anotations.rds")

genes<-readRDS("gene_Anotations38_lcincluded.rds")

gr_gene<- GRanges(
  seqnames =genes$seqid,
  ranges = IRanges(genes$start, end = genes$end,
                   names = genes$ID),
  gene_name=genes$gene_name,
  gene_ID=genes$ID,
  gene_start= genes$start,
  gene_end=genes$end,
  gene_size=genes$end-genes$start,
  gene_type=genes$gene_type)
#gene_chr=genes$Chr

saveRDS(gr_gene, "gr_gene_Anotations38_lcincluded.rds")

#gr_gene
