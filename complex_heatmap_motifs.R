# COMPLEX HEATMAP OF MOTIFS AND RNA-seq

#install.packages("installr")
#library(installr)
#updateR()
#install.packages("dplyr")
library(dplyr)
#install.packages('Seurat')
library(Seurat)
#install.packages("patchwork")
library(patchwork)

# Loading the important repositories # 
require("devtools")
library(ggplot2)
library(dplyr)
library(tidyr)
library(mvoutlier)
#BiocManager::install("limma")
library(limma)
# error no limma?

library(knitr)
library(SingleCellExperiment)

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 

library(scater)
library(RColorBrewer)
library(plot3D)
library(stringr)
library(SAVER)

#library(DoubletFinder) #remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')
#BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)
#BiocManager::install("org.Hs.eg.db")
library(org.Hs.eg.db)
#BiocManager::install("clusterProfiler")
#library(clusterProfiler) this one does not work correctly
#BiocManager::install("ggpubr")
library(ggpubr)
#install.packages("circlize")
library(circlize)
library(cowplot)
#install.packages("clustree")
library(clustree)
#install.packages("grid")
library(grid)

#library('biomaRt') #BiocManager::install("biomaRt") (update none)
#library('bitr') #remotes::install_github("GuangchuangYu/bitr")

#install.packages('ape')
library(ape)
#install.packages("ggplot2")
library("ggplot2")

# install.packages("scran")
# library(scran)
# scran not available for R version 3.5.1

## load in the datasets:
lakorna <- read.table(file = '/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/20210712/Zscoretable.tsv', sep = '\t', header = TRUE, row.names=1)
lakorna <- na.omit(lakorna)
# importing files and selecting unique row names and genes
lakoqnormmotifs <- read.table(file = '/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/data/lako2021/motif_analysis/quantile_sigregions/joined.txt', sep = ' ', header = TRUE, row.names=1)
lakoqnormmotifs <- lakoqnormmotifs[,1:12]
selector <- unique(lakoqnormmotifs[c("genes")])
names <- rownames(selector)
x <- subset(lakoqnormmotifs, rownames(lakoqnormmotifs) %in% names)

# because the motif that contains TP63 is selected with TP73 and TP53 (that are not expressed), this value is manually changed
x$genes[x$genes == "TP73"] <- "TP63"
x$genes == "TP63"

lakovstmotifs <- read.table(file = '/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/data/lako2021/motif_analysis/vst_sigregions/joined.txt', sep = ' ', header = TRUE, row.names=1)
lakovstmotifs <- lakovstmotifs[,1:12]
selector2 <- unique(lakovstmotifs[c("genes")])
names2 <- rownames(selector2)
a <- subset(lakovstmotifs, rownames(lakovstmotifs) %in% names2)

# because the motif that contains TP63 is selected with TP73 and TP53 (that are not expressed), this value is manually changed
a$genes[a$genes == "TP73"] <- "TP63"
a$genes == "TP63"

# subselect overlapping genes in the two datasets
z <- subset(lakorna, rownames(lakorna) %in% (x$genes))
z <- subset(z, rownames(z) %in% (a$genes))

x <- subset(x, x$genes %in% rownames(z))
a <- subset(a, a$genes %in% rownames(z))

rownames(x) <- x$genes
rownames(a) <- a$genes
a <- a[,-1]
x <- x[,-1]
z$ESC <- NULL

matx <- as.matrix(x)
matz <- as.matrix(z)
maty <- as.matrix(a)

# making the complex heatmap
f1 = colorRamp2(c(-5, 0, 5), c("blue", "#EEEEEE", "yellow"), space = "RGB")
f2 = colorRamp2(c(-1, 0, 1), c("blue", "#EEEEEE", "red"), space = "RGB")
f3 = colorRamp2(c(-5, 0, 5), c("blue", "#EEEEEE", "green"), space = "RGB")

ht3 = Heatmap(maty, col = f3, cluster_columns = T, name = "motifs z-score qnorm")

roword <- row_order(ht3)
maty <- maty[roword,]
colord <- column_order(ht3)
maty <- maty[, colord]

matx <- matx[rownames(maty),colnames(maty)]
matz <- matz[rownames(maty),colnames(maty)]


ht1 = Heatmap(matz, col= f2, cluster_rows = F, cluster_columns = F, name = "Z-score RNA count")
ht2 = Heatmap(matx, col = f1, cluster_columns = F, cluster_rows = F, name = "motifs z-score vst")
ht3 = Heatmap(matx, col = f3, cluster_columns = F, cluster_rows=F, name = "motifs z-score qnorm")

## Storing results ##
workdir <- "/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/complexheatmap/"

## setting up results directory
dateoftoday <- gsub("-", "", as.character(Sys.Date()))
resultsdir <- paste0(workdir, dateoftoday)
system(paste("mkdir -p ", resultsdir))
colnames(maty)

mattest <- data.frame(maty)
sel <- mattest[order(mattest$CjS,decreasing = T),]
rownames(sel)
mattest <- data.frame(maty)


# plotting the heatmaps
pdf(paste(resultsdir,'/complexheatmap2.pdf',sep="/") ,width=15,height=8,paper='special')
for (i in colnames(mattest)) {
  print(i)
  sel <- mattest[order(mattest[[i]],decreasing = T),]
  vec <- rownames(sel)[1:15]
  print(vec)
  vec2 <- as.matrix(sel)
  vec3 <- which(rownames(matz) %in% rownames(vec2)[1:15], arr.ind=TRUE)
  htlist <- ht3 + ht2 + ht1 + rowAnnotation(link = anno_mark(at =  vec3,labels = vec))
  draw(htlist, column_title = i)
}
dev.off()

# plotting venndiagram of motifs
library(VennDiagram)
venn.diagram(
  x = list(rownames(selector),rownames(selector2)),
  category.names = c("Set vst","Set qnorm "),
  filename = '#14_venn_diagramm.png',
  output=TRUE
)
