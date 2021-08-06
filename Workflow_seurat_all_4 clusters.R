### SEURAT WORKFLOW ###
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

# load in all datasets
lako.data1 <- Read10X(data.dir = "/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/scRNA358/outs/filtered_feature_bc_matrix/") ## change dataset here
# Initialize the Seurat object with the raw (non-normalized data).
lako1 <- CreateSeuratObject(counts = lako.data1, project = "lako1", min.cells = 3, min.features = 200)
lako1

lako.data2 <- Read10X(data.dir = "/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/data/lako2021/SRR12386359/outs/filtered_feature_bc_matrix/") ## change dataset here
# Initialize the Seurat object with the raw (non-normalized data).
lako2 <- CreateSeuratObject(counts = lako.data2, project = "lako2", min.cells = 3, min.features = 200)
lako2

lako.data3 <- Read10X(data.dir = "/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/data/lako2021/SRR12386360/outs/filtered_feature_bc_matrix/") ## change dataset here
# Initialize the Seurat object with the raw (non-normalized data).
lako3 <- CreateSeuratObject(counts = lako.data3, project = "lako3", min.cells = 3, min.features = 200)
lako3

lako.data4 <- Read10X(data.dir = "/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/data/lako2021/SRR12386361/outs/filtered_feature_bc_matrix/") ## change dataset here
# Initialize the Seurat object with the raw (non-normalized data).
lako4 <- CreateSeuratObject(counts = lako.data4, project = "lako4", min.cells = 3, min.features = 200)
lako4

lako.data5 <- Read10X(data.dir = "/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/data/lako2021/SRR12386362/outs/filtered_feature_bc_matrix/") ## change dataset here
# Initialize the Seurat object with the raw (non-normalized data).
lako5 <- CreateSeuratObject(counts = lako.data5, project = "lako5", min.cells = 3, min.features = 200)
lako5

lako.data6 <- Read10X(data.dir = "/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/data/lako2021/SRR12386363/outs/filtered_feature_bc_matrix/") ## change dataset here
# Initialize the Seurat object with the raw (non-normalized data).
lako6 <- CreateSeuratObject(counts = lako.data6, project = "lako6", min.cells = 3, min.features = 200)
lako6

lako.data7 <- Read10X(data.dir = "/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/data/lako2021/SRR12386367/outs/filtered_feature_bc_matrix/") ## change dataset here
# Initialize the Seurat object with the raw (non-normalized data).
lako7 <- CreateSeuratObject(counts = lako.data7, project = "lako7", min.cells = 3, min.features = 200)
lako7

lako.data8 <- Read10X(data.dir = "/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/data/lako2021/SRR12386368/outs/filtered_feature_bc_matrix/") ## change dataset here
# Initialize the Seurat object with the raw (non-normalized data).
lako8 <- CreateSeuratObject(counts = lako.data8, project = "lako8", min.cells = 3, min.features = 200)
lako8

# merge all files together
# library(SeuratData)
# InstallData("pbmc3k")
# pbmc3k <- LoadData("pbmc3k", type = "pbmc3k.final")
# pbmc3k

lako <- merge(lako1, y = c(lako2, lako3, lako4, lako5, lako6, lako7, lako8), add.cell.ids = c("Lako1", "Lako2", "Lako3", "Lako4", "Lako_cencor", "Lako_limbring", "lako7", "lako8"), project = "LakoAll")
lako

#######################
# converting the data for making a seur_obj_all
#######################

seur_obj_all <- lako
######################################################################
## Filtering of the dataset ##

## Storing results ##
workdir <- "/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq/"

## setting up results directory
dateoftoday <- gsub("-", "", as.character(Sys.Date()))
resultsdir <- paste0(workdir, dateoftoday)
system(paste("mkdir -p ", resultsdir))

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
seur_obj_all[["percent.mt"]] <- PercentageFeatureSet(seur_obj_all, pattern = "^MT-")


### QC tresholds based on Jos smits his tresholds###
seur_obj <- subset(seur_obj_all, subset = nCount_RNA > 2000 & nFeature_RNA > 1000 & percent.mt < 30)

seur_obj <- NormalizeData(seur_obj, normalization.method = "LogNormalize", scale.factor = 10000)

seur_obj <- NormalizeData(seur_obj)

seur_obj <- FindVariableFeatures(seur_obj, selection.method = "vst", nfeatures = 2000)

####################################################
# percentage MTs before and after
#####




pdf(paste(resultsdir,'3.QC_depth_ERCC_MT.pdf',sep="/") ,width=10,height=6,paper='special')
VlnPlot(object = seur_obj_all, features = ("nCount_RNA"), assay = 'RNA') + scale_y_continuous(limits = c(0,20000)) + ggtitle('total UMI counts all cells') 
VlnPlot(object = seur_obj, features = c("nCount_RNA"), assay = 'RNA') + scale_y_continuous(limits = c(0,20000)) + ggtitle('total UMI counts filtered cells') 
VlnPlot(object = seur_obj_all, features = c("nFeature_RNA"), assay = 'RNA') + scale_y_continuous(limits = c(0,5000)) + ggtitle('genes measured all cells') 
VlnPlot(object = seur_obj, features = c("nFeature_RNA"), assay = 'RNA')  + scale_y_continuous(limits = c(0,5000)) + ggtitle('genes measured filtered cells') 
VlnPlot(object = seur_obj_all, features = c("percent.mt"))+ scale_y_continuous(limits = c(0,40)) + ggtitle('% mitochondrial reads all cells') 
VlnPlot(object = seur_obj, features = c("percent.mt"))+ scale_y_continuous(limits = c(0,40)) + ggtitle('% mitochondrial reads filtered cells') 

#VlnPlot(object = seur_obj_all, features = c("pct_counts_ERCC"))+ scale_y_continuous(limits = c(0,1.4)) + ggtitle('% ERCC reads all cells') 
#VlnPlot(object = seur_obj, features = c("pct_counts_ERCC"))+ scale_y_continuous(limits = c(0,1.4)) + ggtitle('% ERCC reads filtered cells')

FeatureScatter(seur_obj_all, feature1 = ("nCount_RNA"), feature2 = ("nFeature_RNA")) + ggtitle('% counts vs genes measured all cells') 
FeatureScatter(seur_obj, feature1 = ("nCount_RNA"), feature2 = ("nFeature_RNA")) + ggtitle('% counts vs genes measured filtered cells') 
dev.off()

#table(seur_obj_all$cell_type)




#################################################################################

Idents(seur_obj) <- "cell_type"
Idents(seur_obj_all) <- "cell_type"

seur_obj <- NormalizeData(
  object = seur_obj, assay = "RNA",
  normalization.method = "LogNormalize", 
  scale.factor = 10000
)

seur_obj <- FindVariableFeatures(seur_obj, selection.method = "vst")
seur_obj <- ScaleData(seur_obj, features = rownames(seur_obj), assay = 'RNA')

seur_obj <- CellCycleScoring(seur_obj, s.features = cc.genes.updated.2019$s.genes, g2m.features = cc.genes.updated.2019$g2m.genes, set.ident = TRUE)
seur_obj <- RunPCA(seur_obj, features = c(cc.genes.updated.2019$s.genes, cc.genes.updated.2019$g2m.genes))

pdf(paste(resultsdir,'5.cell_cycle_markers.pdf',sep="/") ,width=12,height=6,paper='special')
Idents(seur_obj) <- "Phase"
PCAPlot(seur_obj)

#only for multiple cells
#Idents(seur_obj) <- "cell_type"
#PCAPlot(seur_obj)

RidgePlot(seur_obj, features = c("PCNA", "TOP2A", "MCM6", "MKI67"), assay= 'RNA', ncol = 2, slot = "data")
RidgePlot(seur_obj, features = c("S.Score", "G2M.Score"))
          #, ncol = 2, group.by = "cell_type")
dev.off()

##### PCA plots for all genes

library(grid)
library(gridExtra)

seur_obj <- RunPCA(seur_obj, verbose = FALSE, npcs = 50)
mat <- Seurat::GetAssayData(seur_obj, assay = "RNA", slot = "scale.data")
pca <- seur_obj[["pca"]]

# Get the total variance:
total_variance <- sum(matrixStats::rowVars(mat))

eigValues = (pca@stdev)^2  ## EigenValues
varExplained = eigValues / total_variance
Stdev(object = seur_obj[["pca"]])[1]

Idents(seur_obj) <- "cell_type"
pdf(paste0(resultsdir,'/6.Principle_components.pdf') ,width=12,height=6,paper='special')
print(ElbowPlot(seur_obj))
for (pcs in c(1:11)){
  pc1_viz <- pcs*2-1
  pc2_viz <- pcs*2
  y_label = paste0(paste0(paste0("PC",pc2_viz), ' stdev:  '),round(Stdev(object = seur_obj[["pca"]])[pc2_viz],3))
  x_label = paste0(paste0(paste0("PC",pc1_viz), ' stdev:  '),round(Stdev(object = seur_obj[["pca"]])[pc1_viz],3))
  PC_dimred <- DimPlot(seur_obj, reduction = "pca", dims = c(pc1_viz,pc2_viz))+labs(y= y_label, x = x_label)
  PC1_genes <- DimHeatmap(seur_obj, dims = c(pc1_viz), fast = FALSE)
  PC2_genes <- DimHeatmap(seur_obj, dims = c(pc2_viz), fast = FALSE)
  final_plot <- grid.arrange(as_grob(PC_dimred),as_grob(PC2_genes),as_grob(PC1_genes),
                             ncol=2,
                             as.table=TRUE)
  #heights=c(3,1))
  print(final_plot)}
dev.off()

# see which genes contribute the most
#Loadings(seur_obj[["pca"]])

##### generating the UMAP #####
# chosen 30 PCAs was the most variable untill that number
seur_obj <- RunUMAP(seur_obj, dims = 1:10)
pdf(paste(resultsdir,'7a.norm_umap.pdf',sep="/") ,width=8,height=8,paper='special') 
#DimPlot(seur_obj, label = TRUE)
DimPlot(seur_obj, label=FALSE)+ ggtitle("original identity") 
FeaturePlot(seur_obj, features = 'nCount_RNA', cols =c('white','purple'))
FeaturePlot(seur_obj, features = 'nFeature_RNA', cols =c('white','purple'))
FeaturePlot(seur_obj, features = 'percent.mt', cols =c('white','purple'))
DimPlot(seur_obj, group.by = 'Phase')
FeaturePlot(seur_obj, features = 'S.Score', cols =c('white','dodgerblue3'))
FeaturePlot(seur_obj, features = 'G2M.Score', cols =c('white','green4'))
dev.off()

# run this the second time once the resolution is set
pdf(paste(resultsdir,'clusdendro.pdf',sep="/") ,width=12,height=8,paper='special') 
seur_obj <- FindNeighbors(seur_obj, dims = 1:10)
res <- 0.01
cluster_variable_name <- paste0("RNA_snn_res.", res)
seur_obj <- FindClusters(seur_obj, verbose = FALSE, resolution = res, graph.name = "RNA_snn")
seur_obj <- BuildClusterTree(seur_obj)
PlotClusterTree(seur_obj, direction = "downwards")
DimPlot(seur_obj, label=TRUE, label.size = 6, group.by = 'seurat_clusters') + ggtitle("Louvain Clustering") + ggtitle(paste0("cluster resolution ", res))
dev.off()

#lets first perform clustering, after which we can start to quantify gene expression in all the clusters.
pdf(paste(resultsdir, "8.clustering_resolution.pdf", sep = '/'), width = 15, height = 8)
seur_obj <- FindNeighbors(seur_obj, dims = 1:10) 
for (i in 1:20){
  res <- i/20
  cluster_variable_name <- paste0("RNA_snn_res.", res)
  seur_obj <- FindClusters(seur_obj, verbose = FALSE, resolution = res, graph.name = "RNA_snn")
  seur_obj <- BuildClusterTree(seur_obj)
  #p1 <- DimPlot(seur_obj, label=TRUE) + ggtitle(paste0("cluster resolution ", res))
  p2 <- DimPlot(seur_obj, label=TRUE, group.by = 'seurat_clusters') + ggtitle("Louvain Clustering") + ggtitle(paste0("cluster resolution ", res))
  p3 <- ggplot(seur_obj@meta.data, aes(eval(parse(text= cluster_variable_name))))+geom_bar(stat="count")
  
  if (i != 1){
  print(clustree(
    seur_obj,
    prefix = "RNA_snn_res.",
    exprs = c("data", "counts", "scale.data"),
    assay = NULL,
    node_colour = "sc3_stability"
  ))}
  #print(p1)
  print(p2)
  print(p3)
  #PlotClusterTree(object = seur_obj)
  
  
}
dev.off()


############################################################################
#Use this only if it needs to be merged; cell cycle related and such!!!!!!!
############################################################################

#lets cheat and merge the iPSC and LSC cell cycle related clusters:
current.cluster.ids <- c(0, 1, 2, 3)
new.cluster.ids <- c('cell1','cell2','cell3','cell4')
seur_obj$costum_clustering <- plyr::mapvalues(x = seur_obj$seurat_clusters, from = current.cluster.ids, to = new.cluster.ids)
cluster_order <- c('cell1','cell2','cell3','cell4') # change if nessesary
# Re-level object@ident
seur_obj$costum_clustering <- factor(x = seur_obj$costum_clustering, levels = cluster_order)

# add custom labels to the phylogenetic tree
data.tree <- Tool(object = seur_obj, slot = "BuildClusterTree")
data.tree$tip.label <- levels(seur_obj$costum_clustering)
# plot the tree

pdf(paste(resultsdir, "8b.clustering_costum_tree.pdf", sep = '/'), width = 12, height = 8)
#p1 <- DimPlot(seur_obj, label=TRUE, group.by = 'cell_type')
p1 <- DimPlot(seur_obj, label=TRUE,label.size = 8) 
p2 <- DimPlot(seur_obj, label=TRUE,label.size = 8, group.by = 'costum_clustering') + ggtitle("Louvain Clustering") 
p3 <- DimPlot(seur_obj, label=F,label.size = 5, group.by = 'orig.ident') + ggtitle("Louvain Clustering") 
ape::plot.phylo(x = data.tree, direction = "downwards",)
# 
p1
p2
p3
#print(multiplot(p1, p2, p3, cols = 3))
dev.off()

Idents(seur_obj) <- "costum_clustering"

marker_gene_file <- "/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/LSC_Marker_Genes3.csv"
marker_genes_df <- read.table(marker_gene_file, header = TRUE, sep = ',',comment.char = "#", stringsAsFactors = F)
library(BiocParallel)
SnowParam(workers = 2)
for (paper in unique(marker_genes_df$abreviation)){
  marker_dir <- paste0(paste0(resultsdir, '/marker_genes/'),paper)
  print(marker_dir)
  system(paste("mkdir -p ", marker_dir))
  sub_marker_df <- marker_genes_df[marker_genes_df$abreviation == paper,]
  
  for (cell_type in unique(sub_marker_df$name_plot)){
    #print(cell_type)
    marker_genes <- sub_marker_df[sub_marker_df$name_plot == cell_type,]$marker_genes
    marker_genes <- unlist(strsplit(marker_genes, ','))
    
    n <- length(marker_genes)
    k <- 10 ## your LEN
    marker_gene_sets <- split(marker_genes, rep(1:ceiling(n/k), each=k)[1:n])
    
    for (count_type in c('raw_counts', 'seurat_norm_counts','SAVER_imputed_counts')){        #loop  over the different type of normalized data
      #for (count_type in c('seurat_norm_counts')){        #loop  over the different type of normalized data
      
      print(count_type)
      
      pdf(paste0(paste(paste(marker_dir, cell_type, sep = '/'), count_type, sep = '_'),'.pdf'), width = 20, height = 15)
      
      for (genes in  marker_gene_sets){#loop over subsets of max 12 genes to vizualize
        print(genes)
        plot_list1 <- DimPlot(seur_obj, label=TRUE, combine = F, pt.size = 0.1)
        #plot_list1 <- append(plot_list1, DimPlot(seur_obj, label=TRUE, group.by = 'custom_clustering', combine = F, pt.size = 0.1))
        count_plots <- plot_list1
        count_plots_clusters <- plot_list1
        
        if (count_type == 'raw_counts'){
          count_plots <- append(count_plots, try(VlnPlot(seur_obj, combine=F, features = genes, ncol = 4, group.by = 'cell_type', assay = 'RNA', slot = 'counts')))
          count_plots_clusters <- append(count_plots_clusters, try(VlnPlot(seur_obj, combine=F, features = genes, ncol = 4, group.by = 'costum_clustering', assay = 'RNA', slot = 'counts')))}
        if(count_type == 'seurat_norm_counts'){
          count_plots <- append(count_plots, try(VlnPlot(seur_obj, combine=F, features = genes, ncol = 4, group.by = 'cell_type', assay = 'RNA')))
          count_plots_clusters <- append(count_plots_clusters, try(VlnPlot(seur_obj, combine=F, features = genes, ncol = 4, group.by = 'costum_clustering', assay = 'RNA')))}
        
        if(count_type == 'SAVER_imputed_counts'){
          count_plots <- append(count_plots, try(VlnPlot(seur_obj, combine=F, features = genes, ncol = 4, group.by = 'cell_type', assay = 'RNA_SAVER')))
          count_plots_clusters <- append(count_plots_clusters, try(VlnPlot(seur_obj, combine=F, features = genes, ncol = 4, group.by = 'costum_clustering', assay = 'RNA_SAVER')))}
        
        
        plot_list1 <- append(plot_list1,FeaturePlot(seur_obj, features = genes, pt.size = 0.2, ncol = 4, combine=F))
        print(cowplot::plot_grid(plotlist = plot_list1))
        #print(FeaturePlot(seur_obj, features = genes, pt.size = 0.2, ncol = 4))
        
        p1 <- cowplot::plot_grid(plotlist = count_plots)
        title <- ggdraw() + draw_label(paste0(count_type,"of marker genes per timepoint"), fontface = 'bold')
        print(cowplot::plot_grid(title, p1, ncol = 1, rel_heights = c(0.1, 1)))
        p1 <- cowplot::plot_grid(plotlist = count_plots_clusters)
        title <- ggdraw() + draw_label(paste0(count_type,"of marker genes per cluster"), fontface = 'bold')
        print(cowplot::plot_grid(title, p1, ncol = 1, rel_heights = c(0.1, 1)))
      }
      
      dev.off()}
    
    
  }
}

#Lets find marker genes for each cluster
cluster.markers <- FindAllMarkers(seur_obj, only.pos = TRUE)

# making a nice heatmap of expression in each cluster not significant yet
n_genes <- 40
heatmap.markers <- cluster.markers %>% group_by(cluster) %>% top_n(n_genes, avg_log2FC)
DoHeatmap(seur_obj, features = heatmap.markers$gene) + NoLegend()

#filter on only sigi
cluster.markers$gene_name_shorter <- str_sub(cluster.markers$gene,end = -1)
cluster.markers_fc <- cluster.markers[cluster.markers$p_val_adj < 0.01,]

cluster.markers_fc <- cluster.markers_fc[(cluster.markers_fc$avg_log2FC > 0.58) ,]
# avg_logFC...
table(cluster.markers_fc$cluster) # number of significantly expressed genes with a high log2FC (>0.58)

# GO terms
plot_list <- list()
cluster_GOs <- c() 
#expressed_genes <- row.names(counts(sce_qc)[rowSums(counts(sce_qc))>20,])
expressed_genes <- rownames(cluster.markers_fc)

unique(expressed_genes)

#saveRDS(seur_obj, file = "/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq/lakoRNAunannotated2.rds") # 10 PC
saveRDS(seur_obj, file = "/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq/lakoRNAunannotated_4clus.rds") 
#seur_obj <- readRDS("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq/lakoRNAunannotated2.rds") # 10 PC
seur_obj <- readRDS("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq/lakoRNAunannotated_4clus.rds")

###########################################################################
# Appending the GO terms of the unannotated clusters
plot_list <- list()
cluster_GOs <- c() 
expressed_genes <- rownames(cluster.markers_fc)

for (cluster in unique(cluster.markers_fc$cluster)){
  print(cluster)
  df_subset <- cluster.markers_fc[cluster.markers_fc$cluster == cluster,]
  cluster_genes <- unique(str_sub(df_subset$gene,end = -1))
  PC1_ego <- enrichGO(gene = cluster_genes,
                      #universe = expressed_genes,
                      OrgDb         = 'org.Hs.eg.db',
                      keyType       = "SYMBOL",
                      ont           = "BP",
                      pAdjustMethod = "BH",
                      #pvalueCutoff  = 0.01,
                      qvalueCutoff  = 0.05)
  if(is.null(PC1_ego)){
    go_plot1 <- rectGrob(gp=gpar(col=NA))
  }else{
    go_plot1 <- barplot(PC1_ego, showCategory=10)
  }
  plot_list <- c(plot_list, list(go_plot1))
  cluster_GOs <- append(cluster_GOs, (paste('cluster of cell type ', cluster, sep = '')))
}
dev.off()

# plotting the GO-terms
grob2 <- ggarrange(plotlist= plot_list, nrow = 18, ncol = 1, labels = cluster_GOs) #, align = 'v')
pdf(paste(resultsdir,'/all_DEGS_heatmap.pdf',sep="") ,width=10,height=40,paper='special')
print(ggarrange(grob2, ncol =1 , nrow = 18, widths= c(6, 4)))
print(grob2)
dev.off()

# Heatmap cluster marking genes:
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

heatmap.markers <- cluster.markers_fc
cluster.markers_fc$log2FC <- cluster.markers_fc$avg_log2FC 
top10 <- cluster.markers_fc %>%  group_by(cluster) %>%  top_n(n = 5)
top5_genes <- DoHeatmap(seur_obj, features = top10$gene) + NoLegend()
top10 <- cluster.markers_fc %>% group_by(cluster) %>% top_n(n = 10, wt = 'avg_log2FC')
marker_top5 <- DoHeatmap(seur_obj, features = top10$gene) + NoLegend()

n_genes <- 100
heatmap.markers <- cluster.markers_fc %>% group_by(cluster) #%>% top_n(n_genes, avg_logFC)

cell_type_plot <- DimPlot(seur_obj, label=TRUE,pt.size = 2) + ggtitle("Timepoint") 
cluster_plot <- DimPlot(seur_obj, label=TRUE, group.by = 'costum_clustering', pt.size = 2) + ggtitle("Louvain Clustering") 
grob_clustering = ggarrange(plotlist = list(cell_type_plot, cluster_plot) , ncol =1) #labels = cluster_GOs) #align = 'v')
nclust<- length(unique(as.numeric(cluster.markers_fc$cluster)))
RGB_colours_ggplot <- as.list(as.character(gg_color_hue(nclust)))
names <- as.numeric(unique(heatmap.markers$cluster))
col_fun = colorRamp2(names, RGB_colours_ggplot)#make sure the rows correspond to ggplot colour mapping

mat <- as.matrix(seur_obj@assays$RNA@scale.data[heatmap.markers$gene,])
mat <- rbind(mat,seur_obj@meta.data$costum_clustering)
mat <- mat[,order(mat[nrow(mat),])]

column_ha <- HeatmapAnnotation(cluster = mat[nrow(mat),], col =list(cluster = col_fun))
breaks <- mat[nrow(mat),] 
mat <- mat[-nrow(mat),] 
row_ha = rowAnnotation(adj_p_vallue = heatmap.markers$p_val_adj)
clust_heatmap <- grid.grabExpr(draw(Heatmap(mat, column_split = breaks,row_split = as.numeric(heatmap.markers$cluster), cluster_columns = F, cluster_rows = F,show_row_names = F,show_column_names = F,row_names_gp = gpar(fontsize = 6), top_annotation = column_ha, left_annotation = row_ha, row_names_rot = -40)))

###########################################################################
# Plotting the GO terms of the un-annotated clusters and the heatmap
#plot_list <- list()
#cluster_GOs <- c() 
#expressed_genes <- rownames(cluster.markers_fc)

pdf(paste0(resultsdir,'/9.GO_terms_clusters.pdf') ,width=120,height=4,paper='special')
cowplot::plot_grid(plotlist =grob2, ncol = 1, nrow = 1)#, rel_heights = c(1, 4,2), rel_widths = c(3,2,2))
dev.off()

pdf(paste0(resultsdir,'/9.complex_heatmap.pdf') ,width=40,height=24,paper='special') 
cowplot::plot_grid(marker_top5,grob_clustering, ncol = 3, nrow = 1, rel_heights = c(1, 4,2), rel_widths = c(3,2,2))
dev.off()

# Write the markers as a table
write.table(cluster.markers_fc, file = paste0(resultsdir,'/9b.cluster_markers.csv'), sep = ',')

###########################################################################
# Determining the cell types based upon markers
# Annotation notes: see github file regarding annotation

cluster_order <- c('LiCo','StCSC',
                   'FCVes','MECIC')

# Generate new cluster labeling:
current.cluster.ids <- c('cell1','cell2',
                         'cell3','cell4')
new.cluster.ids <- c('LiCo','StCSC',
                     'FCVes','MECIC')

seur_obj$costum_clustering <- plyr::mapvalues(x = seur_obj$costum_clustering, from = current.cluster.ids, to = new.cluster.ids)

# Re-level object@ident
seur_obj$costum_clustering <- factor(x = seur_obj$costum_clustering, levels = cluster_order)
seur_obj$costum_clustering

# Add custom labels to the phylogenetic tree
data.tree <- Tool(object = seur_obj, slot = "BuildClusterTree")
data.tree$tip.label <- levels(seur_obj$costum_clustering)

# Plot the tree and the final clustering
pdf(paste(resultsdir, "13.clustering_costum_final.pdf", sep = '/'), width = 10, height = 8)
#p1 <- DimPlot(seur_obj, label=TRUE, group.by = 'cell_type')
p1 <- DimPlot(seur_obj, label=TRUE,label.size = 8) 
p2 <- DimPlot(seur_obj, label=TRUE,label.size = 8, group.by = 'costum_clustering') + ggtitle("Louvain Clustering") 
p3 <- DimPlot(seur_obj, label=F,label.size = 5, group.by = 'orig.ident') + ggtitle("Louvain Clustering") 
ape::plot.phylo(x = data.tree, direction = "downwards",)
# 
p1
p2
p3
#print(multiplot(p1, p2, p3, cols = 3))
dev.off()

###########################################################################
# Saving the annotated file
saveRDS(seur_obj, file = "/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq/lakoRNAannotated_4cells.rds")
seur_obj <- readRDS("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq/lakoRNAannotated_4cells.rds")

###########################################################################
# Generating metadata file for the cell counts
RNA_metadata <- seur_obj@meta.data

# pre-allocate space
a <- c()
b <- character()

for(i in unique(RNA_metadata$costum_clustering)){
  b <- append(b,i)
  a <- append(a, nrow(RNA_metadata[RNA_metadata$costum_clustering == i,]))
}
popdf <- data.frame(b,a)
names(popdf) <- c("population","number_of_cells")
write.table(popdf, file = paste0(resultsdir,'/','cellcounts_4cells.tsv'), sep = '\t',quote = F,row.names = F)

#################################### generate a pseudo bulk expression plot
# Make a pseudobulk table
seur_obj@meta.data$costum_clustering

sce_qc <- SingleCellExperiment(assays = list(counts = GetAssayData(object = seur_obj, slot = "data")))

pseudobulk_df <- as.data.frame(row.names(counts(sce_qc)))

sce_qc <- SingleCellExperiment(assays = list(counts = GetAssayData(object = seur_obj, slot = "data")))

pseudobulk_df <- as.data.frame(row.names(counts(sce_qc)))

sce_qc$sample <- seur_obj@meta.data$costum_clustering
sce_qc$sample

#sce_qc$rep <- seur_obj@meta.data$orig.ident
#sce_qc$rep

#sce_qc$reps <- "rep1"
#sce_qc$reps[sce_qc$rep == "lako1" | sce_qc$rep == "lako3" | sce_qc$rep == "lako5" | sce_qc$rep == "lako7"] <- "rep2"

# putting values of artificial rep2 and rep1 into the pseudobulkdataframecolumns
for (sample in unique(sce_qc$sample)){
  pseudobulk_df[[sample]] <- rowSums(counts(sce_qc)[,sce_qc$sample == sample])
}

row.names(pseudobulk_df) <- pseudobulk_df$`row.names(counts(sce_qc))`
pseudobulk_df$`row.names(counts(sce_qc))` <- NULL
write.table(data.frame("ID"=rownames(pseudobulk_df),pseudobulk_df), file = paste0(resultsdir,'/pseudobulk_4cells_TPM.tsv'), sep = '\t',quote = F, row.names = F)# this is only if you do not want the replicates
###########################################################################
# Splitting the pseudobulk datasets unbiased into two populations for DEG
sce_qc$sample <- seur_obj@meta.data$costum_clustering
sce_qc$sample

sce_qc$rep <- seur_obj@meta.data$orig.ident
sce_qc$rep

sce_qc$reps <- "rep1"
sce_qc$reps[sce_qc$rep == "lako1" | sce_qc$rep == "lako3" | sce_qc$rep == "lako5" | sce_qc$rep == "lako7"] <- "rep2"

# putting values of artificial rep2 and rep1 into the pseudobulkdataframecolumns
for (sample in unique(sce_qc$sample)){
  pseudobulk_df[[sample]] <- rowSums(counts(sce_qc)[,sce_qc$sample == sample & sce_qc$reps == "rep2"])
  pseudobulk_df[ , paste0(sample,"1")] <- rowSums(counts(sce_qc)[,sce_qc$sample == sample & sce_qc$reps == "rep1"])
}

row.names(pseudobulk_df) <- pseudobulk_df$`row.names(counts(sce_qc))`
pseudobulk_df$`row.names(counts(sce_qc))` <- NULL
write.table(data.frame("ID"=rownames(pseudobulk_df),pseudobulk_df), file = paste0(resultsdir,'/pseudobulk_4cells.tsv'), sep = '\t',quote = F, row.names = F)

###########################################################################
# Additional: interesting genes
Markergenes <- c("IRF1", "IRF5", "FOXC1", "NR2F2", "PAX6", "ELF3", "TNF", "CXCL1", "CXCL2", "CXCL3", "PTGS2", "NFKBIA", "OTX1", "ELK3", "HES4")
pdf(paste(resultsdir,'interesting.pdf'), width = 20, height = 15)
FeaturePlot(seur_obj, features = Markergenes)
dev.off()