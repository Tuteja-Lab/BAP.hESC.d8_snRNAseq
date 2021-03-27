setwd("~/TutejaLab/snn-results_20210308_testing-installation")
# load the modules
library(Seurat)
library(SeuratWrappers)
library(knitr)
library(kableExtra)
library(ggplot2)
library(cowplot)
library(patchwork)
library(metap)
library(multtest)
library(gridExtra)
library(dplyr)
library(stringr)
library(TissueEnrich)
library(gprofiler2)
library(tidyverse)
library(enhancedDimPlot)
library(calibrate)
library(ggrepel)
library(dittoSeq)
library(ComplexHeatmap)

# create seurat object
experiment_name = "BAP"
dataset_loc <- "Z:/Documents/snRNAseq-placenta-project/expression-data"
ids <- c("5pcO2_r1", "5pcO2_r2", "20pcO2_r1", "20pcO2_r2")
d10x.data <- sapply(ids, function(i){
  d10x <- Read10X(file.path(dataset_loc,i,"filtered_feature_bc_matrix"))
  colnames(d10x) <- paste(sapply(strsplit(colnames(d10x),split="-"), '[[' , 1L ), i, sep="-")
  d10x
})
experiment.data <- do.call("cbind", d10x.data)
bapd8.combined <- CreateSeuratObject(
  experiment.data,
  project = "BAPd8",
  min.cells = 10,
  min.genes = 200,
  names.field = 2,
  names.delim = "\\-")
# add the replicate information in the metadata table
head(bapd8.combined@meta.data)
bapd8.combined$log10GenesPerUMI <- log10(bapd8.combined$nFeature_RNA) / log10(bapd8.combined$nCount_RNA)
bapd8.combined$mitoRatio <- PercentageFeatureSet(object = bapd8.combined, pattern = "^MT-")
bapd8.combined$mitoRatio <- bapd8.combined@meta.data$mitoRatio / 100
metadata <- bapd8.combined@meta.data
metadata$cells <- rownames(metadata)
metadata <- metadata %>%
    dplyr::rename(seq_folder = orig.ident,
                  nUMI = nCount_RNA,
                  nGene = nFeature_RNA, 
				  sample = Condition)
mt <- ggplot(dat = metadata, aes(x=nUMI, y=nGene, color=mitoRatio)) + 
  geom_point(alpha = 0.5) + 
  scale_colour_gradient(low = "gray90", high = "black") + labs(colour="MT ratio") +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    panel.border = element_rect(colour = "black")) + 
  xlab("RNA counts") + ylab("Gene counts") +
  stat_smooth(method=lm) +
  facet_wrap(~seq_folder, labeller = labeller(seq_folder = 
                                                c("20pcO2_r1" = "20% Oxygen (rep1)",
                                                  "20pcO2_r2" = "20% Oxygen (rep2)",
                                                  "5pcO2_r1" = "5% Oxygen (rep1)",
                                                  "5pcO2_r2" = "5% Oxygen (rep2)"))) +
  scale_y_continuous(label=comma) +
  scale_x_continuous(label=comma) 