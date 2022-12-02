#-------------------------------------------------------------------------------
#Step 1. Importing data:
#Data from the analysis of the ChiPseqsamples from BLUEPRINT
#Analysis done by Juan David and Enrique in Spain
#-------------------------------------------------------------------------------

#Directory with original data
setwd("/mnt/nocode/data_private/blueprint_segmentations_build38/")
setwd("DRico/chromStates_allBP_hg38_Sep2016/chromHMMresults/")
setwd("Learnmodel_12/healthy_MODEL/")

#bed files
bed_files <- dir("SEGMENTS/",  full.names = TRUE)
bed_files <- bed_files[grepl("_segments.bed", bed_files)]
names <- gsub('SEGMENTS//', '', bed_files)
names <- gsub('_12_BPhg38healthyAllsep2016', '', names)
names(bed_files) <- gsub('_segments.bed', '', names)

#state by line files
statebyline_files <- dir("STATEBYLINE/",  full.names = TRUE)
statebyline_files <-
  statebyline_files[grepl("_statebyline.txt", statebyline_files)]
names <- gsub('STATEBYLINE//', '', statebyline_files)
names <- gsub('_12_BPhg38healthyAllsep2016', '', names)
names(statebyline_files) <- gsub('_statebyline.txt', '', names)

#-------------------------------------------------------------------------------
#Step 2. Importing metadata:
#The file with all metadata of samples from BLUEPRINT
#-------------------------------------------------------------------------------
all_bluprint_metadata <-
  read.delim(
    "~/../../mnt/nocode/juliana/PID/Chromatin_States_Blueprint/20160816.data",
    sep = "\t",
    header = TRUE
  )


df <-
  all_bluprint_metadata[all_bluprint_metadata$SAMPLE_NAME %in% names(bed_files), ]

#Subset only columns of interest in this step
SubDf <-
  dplyr::select(
    df,
    "BIOMATERIAL_TYPE",
    "SAMPLE_NAME",
    "CELL_TYPE",
    "TISSUE_TYPE",
    "DONOR_SEX",
    "DONOR_ID"
  )

df_sub <- unique(SubDf)

setwd("/mnt/nocode/juliana/PID/ThemeII")

write.table(df_sub, "OUTPUT/METADATA_paper.txt", sep = "\t")
#df_sub is the metadata file with the sample names of the files with unpublished chromatin states


#-------------------------------------------------------------------------------
# PLOT with sample sizes 
#-------------------------------------------------------------------------------
mymetadata <- read.delim("INPUTS/METADATA_paper2.txt", sep = "\t")
celltypes<-read.delim("INPUTS/celltype_data.txt")

Sample_number <- table(mymetadata$CELL_TYPE)
Sample_number <- data.frame(Sample_number[Sample_number > 0])
names(Sample_number) <- c("CELL_TYPE", "SAMPLE_SIZE")
Sample_number <- left_join(Sample_number,celltypes)
Sample_number <- left_join(
  Sample_number,
  unique(mymetadata[,c("CELL_TYPE", "CELL_TYPE_Abb")]) 
  )

Sample_number<-Sample_number %>% tidyr::unite("z",
                              c("Linage", "Progen_1", "Progen_2",
                                "Progen_3", "MAIN_CELL_TYPE",
                                "Detail_1", "Detail_2"),
                              sep = "-",
                              remove = FALSE)

Sample_number <- Sample_number[order(Sample_number$z, decreasing = TRUE), ]

Sample_number$Linage <-
  factor(Sample_number$Linage , levels = rev(unique(celltypes$Linage)))

Sample_number$z <-
  factor(Sample_number$z , levels = rev(levelsSample_number$z))

Sample_number$CELL_TYPE <-
  factor(Sample_number$CELL_TYPE, levels = unique(celltypes$CELL_TYPE))

Sample_number$MAIN_CELL_TYPE <-
  factor(Sample_number$MAIN_CELL_TYPE, 
         levels = rev(unique(celltypes$MAIN_CELL_TYPE)))

Sample_number$Progen_1 <-
  factor(Sample_number$Progen_1, 
         levels = unique(celltypes$Progen_1))


Sample_number$CELL_TYPE_Abb <-
  factor(Sample_number$CELL_TYPE_Abb, levels = rev(unique(celltypes$CELL_TYPE_Abb)))

library(ggpubr)

p <-
  ggdotchart(
    Sample_number,
    y = "SAMPLE_SIZE",
    x = "CELL_TYPE_Abb",
    color = "Linage",
    #"MAIN_CELL_TYPE",                                 # Color by groups
    palette = c("#00AFBB", "#E7B800", "#FC4E07"),
    # Custom color palette
    add = "segments",
    # Add segments from y = 0 to dots
    rotate = TRUE,
    # Rotate vertically
    group = "CELL_TYPE_Abb" ,
    # Order by groups
    dot.size = 2,
    # Large dot size                     # Add mpg values as dot labels
    font.label = list(
      color = "white",
      size = 2,
      vjust = 0.5
    ),
    # Adjust label parameters
    xlab = "",
    ylab = "Number of Biological Replicates",
    legend.title = "",
    ggtheme = theme_pubr()                        # ggplot2 theme
  )


p<-p+ theme(plot.margin = margin(t = 0,  # Top margin
                              r = 20,  # Right margin
                              b = 0,  # Bottom margin
                              l = 40)) # Left margin

# library(cowplot)
# 
# ggdraw(p) +
#   draw_line(#"DC"
#     x = c(0.19, 0.19),
#     y = c(0.73, 0.77),
#     color = "black", size = 1
#   ) +
#   annotate("text", x = 0.17, y = 0.75, label ="II.2" , size = 4, angle =90) +
#   draw_line(#"MF"
#     x = c(0.19, 0.19),
#     y = c(0.785, 0.85),
#     color = "black", size = 1
#   ) +
#   annotate("text", x = 0.17, y = 0.82, label ="II.3" , size = 4, angle =90)+
#   draw_line(#"Neutrophils"
#     x = c(0.19, 0.19),
#     y = c(0.545, 0.67),
#     color = "black", size = 1
#   ) +
#   annotate("text", x = 0.17, y = 0.61, label = "II.1", size = 4, angle =90)+
#   draw_line(#"T cell"
#     x = c(0.14, 0.14),
#     y = c(0.375, 0.495),
#     color = "black", size = 1
#   ) +
#   annotate("text", x = 0.12, y = 0.435, label = "I.3", size = 4, angle =90)+
#   draw_line(#B cell
#     x = c(0.14, 0.14),
#     y = c(0.215, 0.36),
#     color = "black", size = 1
#   ) +
#   annotate("text", x = 0.12, y = 0.28, label = "I.2", size = 4, angle =90)+
#   draw_line(#Granulocyte
#     x = c(0.14, 0.14),
#     y = c(0.51, 0.67),
#     color = "black", size = 1
#   ) +
#   annotate("text", x = 0.12, y = 0.6, label = "I.4", size = 4, angle =90)+
#   draw_line(#Mono derived
#     x = c(0.14, 0.14),
#     y = c(0.68, 0.85),
#     color = "black", size = 1
#   ) +
#   annotate("text", x = 0.12, y = 0.76, label = "I.5", size = 4, angle =90)+
#   draw_line(#MEP
#     x = c(0.14, 0.14),
#     y = c(0.865, 0.90),
#     color = "black", size = 1
#   ) +
#   annotate("text", x = 0.12, y = 0.885, label = "I.6", size = 4, angle =90)+
#   draw_line(#"End"
#     x = c(0.14, 0.14),
#     y = c(0.11, 0.17),
#     color = "black", size = 1
#   ) +
#   annotate("text", x = 0.12, y = 0.138, label = "I.1", size = 4, angle =90)
# 
#-------------------------------------------------------------------------------
#Step 3. This block of code has to be modified
# according to the samples of interest
#-------------------------------------------------------------------------------
df_sub <- read.delim("INPUTS/METADATA_paper.txt", sep = "\t")

mycelltypes <-
  c("monocyte", "macrophage", "osteoclast", "dendritic")

Monocyte_derived <-
  df_sub[grepl(paste(mycelltypes, collapse = "|"), df_sub$CELL_TYPE), ]

Monocyte_derived_venous <-
  Monocyte_derived[Monocyte_derived$TISSUE_TYPE == "venous blood", ]

Sample_number <- table(Monocyte_derived_venous$CELL_TYPE)
Sample_number <- data.frame(Sample_number[Sample_number > 0])
names(Sample_number) <- c("CELL_TYPE", "SAMPLE_SIZE")


write.table(Monocyte_derived_venous,
            "OUTPUT/METADATA_Monolinage_venous.txt",
            sep = "\t")
write.table(Sample_number, "OUTPUT/SAMSIZE_Monolinage_venous_male.txt", sep = "\t")

# #Subset of files to use
# files  <-
#   statebyline_files[grepl(
#     paste(Monocyte_derived_venous$SAMPLE_NAME, collapse = "|"),
#     names(statebyline_files)
#   )]
