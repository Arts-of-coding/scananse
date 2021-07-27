# setting up libs for installing SnapATAC package
library(doSNOW)
library(devtools)
#install.packages("plot3D")
library(plot3D)
library(rhdf5)
library(Rhdf5lib)
library(BiocParallel)
library(GenomicRanges)
#devtools::install_github("r3fang/SnapATAC") #skip all updates
library(SnapATAC)

## Storing results ##
workdir <- "/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scATAC"

## setting up results directory
dateoftoday <- gsub("-", "", as.character(Sys.Date()))
resultsdir <- paste0(workdir, dateoftoday)
system(paste("mkdir -p ", resultsdir))

# loading in and joining the snap files
x.sp1 = createSnap(file="/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/data/lako2021/snap_files/ataclako4.snap", sample="atac_lako4",num.cores=1)
summarySnap(x.sp1)
x.sp2 = createSnap(file="/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/data/lako2021/snap_files/ataclako3.snap", sample="atac_lako3",num.cores=1)
summarySnap(x.sp2)
x.sp3 = createSnap(file="/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/data/lako2021/snap_files/ataclako1.snap", sample="atac_lako1",num.cores=1)
summarySnap(x.sp3)
x.sp0 = createSnap(file="/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/data/lako2021/snap_files/ataclako2.snap", sample="atac_lako2",num.cores=1)
summarySnap(x.sp0)


x.sp = snapRbind(x.sp1,x.sp2)
x.sp = snapRbind(x.sp,x.sp3)
x.sp = snapRbind(x.sp,x.sp0)
# x.sp = x.sp0 #test pmat
#barcodes = read.csv("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/data/lako2021/single/singlecell2.csv", header = T) #test pmat
summarySnap(x.sp)
x.sp1@pmat
addPmatToSnap(x.sp1)

# loading in and joining the barcode files
library(plyr)
library(readr)
mydir = "/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/data/lako2021/single/"
myfiles = list.files(path=mydir, pattern="*.csv", full.names=TRUE)
myfiles

barcodes = ldply(myfiles, read_csv)
barcodes = barcodes[2:nrow(barcodes),]
head(barcodes)

barcodes$logUMI = log10(barcodes$passed_filters + 1)
head(barcodes$logUMI)

#########################################
# QC

pdf(paste(resultsdir,'barcodeQC1_all.pdf',sep=""),width=6,height=6,paper='special') 
for (name in unique(x.sp@sample)){
  plotBarcode(
    main.title = name,
    obj=x.sp[x.sp@sample==name,], 
    pdf.file.name=NULL, 
    pdf.width=7, 
    pdf.height=7, 
    col="grey",
    border="grey",
    breaks=50)
  mtext(name, outer=TRUE,  cex=1, line=-1)
}
dev.off()

# figure with barcodes?
pdf(paste(resultsdir,'barcodeQC_UMI1_all.pdf',sep=""),width=6,height=6,paper='special') 
hist(barcodes$logUMI, main= "frequency of unique counts barcodes",xlab = "log10 unique counts")
dev.off()

idd = which(barcodes$logUMI>=3.75 & barcodes$logUMI <= 5)
barcodes <- barcodes[idd,]

#new figure with barcodes?
pdf(paste(resultsdir,'barcodeQC_UMI2_all.pdf',sep=""),width=6,height=6,paper='special') 
hist(barcodes$logUMI, main= "frequency of unique counts barcodes",xlab = "log10 unique counts")
dev.off()

### QC
x.sp = filterCells(
  obj=x.sp, 
  subset.names=c("fragment.num"),
  low.thresholds=c(5000),
  high.thresholds=c(500000)
)

x.sp = filterCells(
  obj=x.sp, 
  subset.names=c('mito.ratio'),
  low.thresholds=c(0),
  high.thresholds=c(0.2)
)

x.sp = filterCells(
  obj=x.sp, 
  subset.names=c('dup.ratio'),
  low.thresholds=c(0),
  high.thresholds=c(0.8)
)
summary(x.sp)

pdf(paste(resultsdir,'barcodeQC2_all.pdf',sep=""),width=6,height=6,paper='special') 
for (name in unique(x.sp@sample)){
  plotBarcode(
    main.title = name,
    obj=x.sp[x.sp@sample==name,], 
    pdf.file.name=NULL, 
    pdf.width=7, 
    pdf.height=7, 
    col="grey",
    border="grey",
    breaks=50)
  mtext(name, outer=TRUE,  cex=1, line=-1)
}
plotBarcode(
  main.title = "all",
  obj=x.sp, 
  pdf.file.name=NULL, 
  pdf.width=7, 
  pdf.height=7, 
  col="grey",
  border="grey",
  breaks=50)
mtext(name, outer=TRUE,  cex=1, line=-1)
dev.off()

barcode.shared = intersect(x.sp@barcode, barcodes$barcode);
x.sp = x.sp[match(barcode.shared, x.sp@barcode),];
barcodes = barcodes[match(barcode.shared, barcodes$barcode),];
x.sp@metaData = barcodes;
summary(x.sp)

# matrix binarisation
x.sp = addBmatToSnap(x.sp, bin.size=5000)
x.sp = makeBinary(x.sp, mat="bmat")

# adding the promoter ratio to the barcodes file
promoter.df = read.table("promoter.bed")
promoter.gr = GRanges(promoter.df[,1], IRanges(promoter.df[,2], promoter.df[,3]))
ov = findOverlaps(x.sp@feature, promoter.gr)
x.sp@feature

# find promoter ratio and exclude overlapping bins 
idy = queryHits(ov)
idy = unique(idy)

log_cov = log10(SnapATAC::rowSums(x.sp, mat="bmat")+1)

x.sp@metaData$log_cov = log_cov
promoter_ratio = Matrix::rowSums(x.sp@bmat[,idy]) / Matrix::rowSums(x.sp@bmat);

pdf(paste(resultsdir,'barcodeQC3_promcount_all.pdf',sep=""),width=6,height=6,paper='special') 
plot(log_cov, promoter_ratio, cex=0.5, col="grey", xlab="log(count)", ylab="FIP Ratio", ylim=c(0,1 ))
dev.off()


idz = which(promoter_ratio > 0.2 & promoter_ratio < 0.8 & log_cov > 3.6)
x.sp = x.sp[idz,]
x.sp

promoter_ratio_filt = Matrix::rowSums(x.sp@bmat[,idy]) / Matrix::rowSums(x.sp@bmat);

pdf(paste(resultsdir,'barcodeQC4_promcount_all.pdf',sep=""),width=6,height=6,paper='special') 
plot(x.sp@metaData$log_cov, promoter_ratio_filt, cex=0.5, col="grey", xlab="log(count)", ylab="FIP Ratio", ylim=c(0,1 ))
dev.off()

#bin filtering
system("wget http://mitra.stanford.edu/kundaje/akundaje/release/blacklists/hg38-human/hg38.blacklist.bed.gz");
library(GenomicRanges);
black_list = read.table("hg38.blacklist.bed.gz")
black_list.gr = GRanges(
  black_list[,1], 
  IRanges(black_list[,2], black_list[,3])
);
idy = queryHits(findOverlaps(x.sp@feature, black_list.gr));
if(length(idy) > 0){x.sp = x.sp[,-idy, mat="bmat"]};
x.sp
#remove unwanted chromosomes
chr.exclude = seqlevels(x.sp@feature)[grep("random|chrM", seqlevels(x.sp@feature))];
idy = grep(paste(chr.exclude, collapse="|"), x.sp@feature);
if(length(idy) > 0){x.sp = x.sp[,-idy, mat="bmat"]};
x.sp

# remove invariant features
bin.cov = log10(Matrix::colSums(x.sp@bmat)+1);
binhist <- hist(
  bin.cov[bin.cov > 0], 
  xlab="log10(bin cov)", 
  main="log10(Bin Cov)", 
  col="lightblue", 
  xlim=c(0, 5)
);

pdf(paste(resultsdir,'bins_all.pdf',sep=""),width=6,height=6,paper='special') 
plot(binhist)
dev.off()

# was this the wierd step??
idx = which(Matrix::rowSums(x.sp@bmat) > 1000);
x.sp = x.sp[idx,];
x.sp

bin.cutoff = quantile(bin.cov[bin.cov > 0], 0.95);
idy = which(bin.cov <= bin.cutoff & bin.cov > 0);
x.sp = x.sp[, idy, mat="bmat"];
x.sp

# remove low mapq cells
x.sp = x.sp[x.sp@metaData$lowmapq<5000,]
x.sp

sum(x.sp@sample == "atac_lako1")
# setting core usage for the server below is intensive
library(BiocParallel)
SnowParam(workers = 2)

#dimensionality reduction
x.sp = runDiffusionMaps(
  obj=x.sp,
  input.mat="bmat", 
  num.eigs=30
)

pdf(paste(resultsdir,'numPCAs_all.pdf',sep=""),width=6,height=6,paper='special') 
plotDimReductPW(
  obj=x.sp, 
  eigs.dims=1:30,
  point.size=0.3,
  point.color="grey",
  point.shape=19,
  point.alpha=0.6,
  down.sample=500,
  pdf.file.name=NULL, 
  pdf.height=7, 
  pdf.width=7
)
dev.off()

plotDimReductElbow(x.sp)
################################ annotating clusters
# selected number of 11 dims:
x.sp = runKNN(
  obj=x.sp,
  eigs.dims=1:10,
  k=10
)

#install.packages("leiden")
# library(reticulate)
# library(igraph)
# path_to_python <- "/vol/mbconda/julian/envs/kb_snapatac/bin/python"
# use_python(path_to_python, required = T)
# py_config()

library("leiden")

x.sp=runCluster(
  obj=x.sp,
  resolution=0.2,
  tmp.folder=tempdir(),
  louvain.lib="leiden",
  seed.use=10
)
x.sp@metaData$cluster = x.sp@cluster

#######################################
## UMAP VIZ
library(umap)

x.sp@cluster

x.sp = runViz(
  obj=x.sp, 
  tmp.folder=tempdir(),
  dims=2,
  eigs.dims=1:10, 
  method="umap",
  seed.use=10
)
pdf(paste(resultsdir,'visualisation_UMAP_all.pdf',sep=""),width=6,height=6,paper='special') 
plotViz(
  obj=x.sp, 
  method='umap', 
  point.size=0.1, 
  point.shape=19, 
  point.alpha=0.8, 
  point.color=x.sp@cluster,
  text.add=TRUE,
  text.size=1.2,
  text.color="black",
  text.halo.add=TRUE,
  text.halo.color="white",
  text.halo.width=0.2,
  down.sample=10000,
  pdf.file.name=NULL,
  pdf.width=7, 
  pdf.height=7
)

plotViz(
  obj= x.sp,
  method="umap", 
  main="Cluster",
  point.color=x.sp@cluster, 
  point.size=0.1, 
  point.shape=19, 
  text.add=TRUE,
  text.size=0.5,
  text.color="black",
  down.sample=10000,
  legend.add=FALSE
)

plotViz(
  obj= x.sp,
  method="umap", 
  main="Sample",
  point.size=0.1, 
  point.shape=19, 
  point.color=x.sp@sample, 
  text.add=FALSE,
  text.size=0.5,
  text.color="black",
  down.sample=10000,
  legend.add=TRUE
)
dev.off()

count(x.sp@sample == "atac_lako1")
count(x.sp@sample == "atac_lako2")
count(x.sp@sample == "atac_lako3")
count(x.sp@sample == "atac_lako4")

# calculate the ensemble signals for each cluster
ensemble.ls = lapply(split(seq(length(x.sp@cluster)), x.sp@cluster), function(x){
  SnapATAC::colMeans(x.sp[x,], mat="bmat");
})

# cluster using 1-cor as distance  
hc = hclust(as.dist(1 - cor(t(do.call(rbind, ensemble.ls)))), method="ward.D2")
plot(hc, hang=-1, xlab="")

pdf(paste(resultsdir,'cluster_tree_all_unannotated.pdf',sep=""),width=6,height=6,paper='special') 
plot(hc, hang=-1, xlab="")
dev.off()

##########################################################################################################3
# gene annotation
gene_df <- read.table("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/gencode.v30.gene.annotation.gtf")
genes.gr = GRanges(gene_df[,1], 
                   IRanges(gene_df[,2], gene_df[,3]), 
                   name=gene_df[,4]
)
genes.gr2  <- genes.gr[!duplicated(genes.gr$name),]

marker.genes = c(
  "CPVL","PAX6","TP63","ACKR1","POSTN",'KERA')

genes.sel.gr2 <- genes.gr2[which(genes.gr2$name %in% marker.genes)]
x.sp4 = createGmatFromMat(
  obj=x.sp, 
  genes= genes.sel.gr2,
  do.par=TRUE#,
  #num.cores=10)
)

# as.integer(SnapATAC::rowSums(x.sp4, mat="gmat"))
# 
# boxPlotFeature(x.sp4,feature = as.integer(SnapATAC::rowSums(x.sp4, mat="gmat")),)
# 
# 
# 
# as.integer(SnapATAC::rowSums(x.sp4@gmat[,"TP63"]))
#x.sp4@gmat[,"TP63"]

# test boxplot function
#boxPlotFeature(x.sp5,feature = x.sp5@gmat[,"PMEL"])

x.sp4 = scaleCountMatrix(
  obj=x.sp4, 
  cov=as.integer(SnapATAC::rowSums(x.sp, mat="bmat")),
  mat="gmat",
  method="logRPM")

x.sp4 = runMagic(
   obj=x.sp4,
   input.mat="gmat",
   step.size=3
 )

pdf(paste(resultsdir,'Marker_genes_sample_all.pdf',sep=""),width=6,height=6,paper='special') 
plotViz(
  obj=x.sp4,
  method="umap", 
  main="lako_atac_all",
  point.color=x.sp4@cluster, 
  point.size=0.5, 
  point.shape=19, 
  point.alpha=0.8, 
  text.add=FALSE,
  text.size=1.5,
  text.color="black",
  text.halo.add=TRUE,
  text.halo.color="white",
  text.halo.width=0.2,
  down.sample=10000,
  legend.add=TRUE
)

for(i in 1:length(marker.genes)){
  plotFeatureSingle(
    obj=x.sp4,
    feature.value=x.sp4@gmat[, marker.genes[i]],
    method="umap", 
    main=marker.genes[i],
    point.size=0.1, 
    point.shape=19, 
    down.sample=10000,
    quantiles=c(0, 1)
  )}
dev.off()

##############################################################################
# now all marker genes
####################################################3
# gene expression upon the UMAP
library(BiocParallel)
SnowParam(workers = 2)

marker_gene_file <- "/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/LSC_Marker_Genes3.csv"
marker_genes_df <- read.table(marker_gene_file, header = TRUE, sep = ',',comment.char = "#", stringsAsFactors = F)

#Add an empty list where markers can be stored
mylist <- c()

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
    
    for (genes in  marker_gene_sets){ 
      #if (genes != 0){#loop over subsets of max 12 genes to vizualize
      print(genes)
      
      for (z in genes){
        print(z)
        mylist <- c(mylist, z)
      }
    }
  }
}

# remove duplicates from the list
mylist <- unique(mylist)

genes.sel.gr2 <- genes.gr2[which(genes.gr2$name %in% mylist)]
x.sp5 = createGmatFromMat(
  obj=x.sp, 
  genes= genes.sel.gr2,
  do.par=TRUE#,
  #num.cores=10)
)

x.sp5 = scaleCountMatrix(
  obj=x.sp5, 
  cov=as.integer(SnapATAC::rowSums(x.sp, mat="bmat")),
  mat="gmat",
  method="logRPM")

# x.sp5 = runMagic(
#    obj=x.sp5,
#    input.mat="gmat",
#    step.size=3
#  )

# all lako 2021 marker genes
x.sp5@gmat[,]

for (paper in unique(marker_genes_df$abreviation)){
  marker_dir <- paste0(paste0(resultsdir, '/marker_genes_nomagic/'),paper)
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
    
    
    pdf(paste0(paste(paste(marker_dir, cell_type, sep = '/'), sep = '_'),'.pdf'), width = 10, height = 6)
    
    for (genes in  marker_gene_sets){ 
      #if (genes != 0){#loop over subsets of max 12 genes to vizualize
      print(genes)
      
      for (z in genes){
        print(z)
        plotFeatureSingle(
           obj=x.sp5,
           feature.value=x.sp5@gmat[,z],
           method="umap", 
           main=z,
           point.size=0.5, 
           point.shape=19, 
           down.sample=10000,
           quantiles=c(0.01, 0.99)
         )
        boxPlotFeature(x.sp5,feature = x.sp5@gmat[,z],main = z,ylab = "logRPM norm accesibility")
      }
    }
    dev.off()}
}

saveRDS(x.sp,file="/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/lakoatacxspall.rds")
x.sp <- readRDS("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/lakoatacxspall.rds")

# plotting QC value of each cluster
library(ggplot2)

# QC of the clusters
pdf(paste(resultsdir,'QC clusters.pdf',sep=""),width=6,height=6,paper='special')
data <- data.frame(name= as.factor(x.sp@metaData$cluster), value= x.sp@metaData$log_cov)
ggplot(data, aes(x=as.factor(data$name), y= data$value)) +
geom_boxplot() +
labs(x= "cluster",y="log10 coverage")

data <- data.frame(name= as.factor(x.sp@metaData$cluster), value= x.sp@metaData$logUMI)
ggplot(data, aes(x=as.factor(data$name), y= data$value)) +
geom_boxplot()+
labs(x= "cluster",y="log10 UMI count")

data <- data.frame(name= as.factor(x.sp@metaData$cluster), value= x.sp@metaData$mitochondrial)
ggplot(data, aes(x=as.factor(data$name), y= data$value)) +
geom_boxplot()+
labs(x= "cluster",y="counts mitochondrial")

data <- data.frame(name= as.factor(x.sp@metaData$cluster), value= x.sp@metaData$passed_filters)
ggplot(data, aes(x=as.factor(data$name), y= data$value)) +
geom_boxplot()+
labs(x= "cluster",y="value passed filters")
dev.off()

data <- data.frame(name= as.factor(x.sp@metaData$cluster), value= x.sp@metaData$duplicate)
ggplot(data, aes(x=as.factor(data$name), y= data$value)) +
  geom_boxplot()+
  labs(x= "cluster",y="duplicate value")

# more quality control
pdf(paste(resultsdir,'QC clusters 2 on UMAP.pdf',sep=""),width=6,height=6,paper='special')
plotFeatureSingle(
  obj=x.sp,
  feature.value=x.sp@metaData$lowmapq,
  method="umap", 
  main="QC low mapq",
  point.size=0.5, 
  point.shape=19, 
  down.sample=10000,
  quantiles=c(0.01, 0.99)
)
# is high in cluster 4; perhaps regions that overlap within cell populations?

plotFeatureSingle(
  obj=x.sp,
  feature.value=x.sp@metaData$duplicate,
  method="umap", 
  main="QC duplicate",
  point.size=0.5, 
  point.shape=19, 
  down.sample=10000,
  quantiles=c(0.01, 0.99)
)

plotFeatureSingle(
  obj=x.sp,
  feature.value=x.sp@metaData$passed_filters,
  method="umap", 
  main="QC passed filters",
  point.size=0.5, 
  point.shape=19, 
  down.sample=10000,
  quantiles=c(0.01, 0.99)
)

plotFeatureSingle(
  obj=x.sp,
  feature.value=x.sp@metaData$chimeric,
  method="umap", 
  main="QC chimeric",
  point.size=0.5, 
  point.shape=19, 
  down.sample=10000,
  quantiles=c(0.01, 0.99)
)
# little bit cluster 4

plotFeatureSingle(
  obj=x.sp,
  feature.value=x.sp@metaData$unmapped,
  method="umap", 
  main="QC unmapped",
  point.size=0.5, 
  point.shape=19, 
  down.sample=10000,
  quantiles=c(0.01, 0.99)
)
# little bit cluster 4

plotFeatureSingle(
  obj=x.sp,
  feature.value=x.sp@metaData$logUMI,
  method="umap", 
  main="QC log unique counts",
  point.size=0.5, 
  point.shape=19, 
  down.sample=10000,
  quantiles=c(0.01, 0.99)
)
# little bit cluster 4

plotFeatureSingle(
  obj=x.sp,
  feature.value=x.sp@metaData$log_cov,
  method="umap", 
  main="QC log coverage",
  point.size=0.5, 
  point.shape=19, 
  down.sample=10000,
  quantiles=c(0.01, 0.99)
)
# little bit cluster 4
dev.off()
#############################################################
#### recluster according to accessibility of marker genes
current.cluster.ids <-  c( "1","2","3","4","5","6")
new.cluster.ids <- c( "cell1","cell2","cell3","cell4","cell5","cell6")
# 
x.sp@cluster <- plyr::mapvalues(
   x = x.sp@cluster, 
   from = current.cluster.ids, 
   to = new.cluster.ids
)

x.sp@cluster

pdf(paste(resultsdir,'visualisation_annotated_UMAP_all.pdf',sep=""),width=6,height=6,paper='special') 
plotViz(
  obj= x.sp,
  method="umap", 
  main="Cluster",
  point.color=x.sp@cluster, 
  point.size=0.1, 
  point.shape=19, 
  text.add=TRUE,
  text.size=1,
  text.color="black",
  down.sample=10000,
  legend.add=FALSE
)
dev.off()
####################################
# script from Rebecca and Jelle

# Create narrowpeak per cluster
#Create a narrowpeak file for each of the clusters using MACS2
#```{r pseudo_bulk_cluster, echo = TRUE, results="hide"}
#Create directory for narrowpeak files of cluster

dir.create(paste0(resultsdir, "clusters/"))
setwd(paste0(resultsdir, "clusters/"))

#Get amount of cells per cluster:
#summary(x.sp@metaData$cluster) for all 8 earlier found
summary(x.sp@cluster)
SnowParam(workers = 2)
# call peaks for all cluster with more than 80 cells
clusters.sel = names(table(x.sp@cluster))#[which(table(x.sp@cluster) > 80)]

system("which snaptools")
system("which macs2")
SnowParam(workers = 2)

peaks.ls = mclapply(seq(clusters.sel), function(i){
  print(clusters.sel[i])
  peaks = runMACS(
    obj=x.sp[which(x.sp@cluster==clusters.sel[i]),], 
    output.prefix=paste0("LAKO.", gsub(" ", "_", clusters.sel)[i]),
    path.to.snaptools="/vol/mbconda/julian/envs/kb_snapatac/bin/snaptools",
    path.to.macs="/vol/mbconda/julian/envs/kb_snapatac/bin/macs2",
    gsize="hs", # mm, hs, etc
    buffer.size=500, 
    num.cores=1,
    macs.options="--nomodel --shift 100 --ext 200 --qval 5e-2 -B --SPMR",
    tmp.folder="/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/data/lako2021/R_peaks_scATAC"
  )
  peaks
})#, mc.cores=5)

# peaks.ls = mclapply(seq(clusters.sel), function(i){
#   print(clusters.sel[i])
#   runMACS(
#     obj=x.sp[which(x.sp@cluster==clusters.sel[i]),],
#     output.prefix=paste0(pseudo_bulk_name, "_cluster", clusters.sel[i]),
#     path.to.snaptools="/vol/mbconda/julian/envs/kb_snapatac/bin/snaptools",
#     path.to.macs = "/vol/mbconda/julian/envs/kb_snapatac/bin/macs2",
#     gsize="hs",
#     buffer.size=500,
#     num.cores=1,
#     macs.options="--nomodel --shift 37 --ext 73 --qval 5e-2 -B --SPMR --call-summits",
# #macs.options="--nomodel --shift 100 --ext 200 --qval 5e-2 -B --SPMR", earlier paramets
#     tmp.folder=tempdir()
#   )
# }, mc.cores=5)
#Set parameter to create bigwigs
Sys.setenv(bdg_location=paste0("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/peaks_scATAC_20210603/"))

# Settings v3: --nomodel --shift 37 --ext 73 --qval 5e-2 -B --SPMR --call-summits
# Settings BAMPE (v2): "--format BAMPE --qval 5e-2 -B --SPMR --call-summits"
# Settings Jos (change macs): "--nomodel --shift 100 --ext 200 --qval 5e-2 -B --SPMR"
# Settings Original: "--nomodel --shift 37 --ext 73 --qval 1e-2 -B --SPMR --call-summits"
#``

###############################################
#continue from below
#Create bigwigs
#Create a bigwig from the processed data based on the bedgraph produced above. Folder containing date should be changed to the correct name.

#################################################################
## bash commands
#```{bash create_bigwigs, echo= TRUE}
#source activate /mbshome/jverploegen/anaconda3/envs/genome_tools

# ##set directory where output bedgraph is located
# cd $bdg_location
# 
#
# For all bdg files: exclude GL000*, KI* and chrM (not in the hg38.chrom.sizes file) and remove for each cluster the b and '.
# grep -Ev 'GL000|KI|chrM' $bdg_file > $ bdg_file2
# nice -30 sed -i ' s/b//g' $bdg_file2 && sed -i " s/[']//g" $bdg_file2
# 
# conda activate bigwig
#
# ##Loop to create bigwig from all bedgraphs
# for bdg_file in *bdg
# do
# ##sort bedgraph and remove random chromosomes etc.
# sort -k1,1 -k2,2n $bdg_file > conversion_file.bdg
# 
# nice -30 bedGraphToBigWig conversion_file.bdg /ceph/rimlsfnwi/data/moldevbio/zhou/jarts/cellrangeratac/hg38.chrom.sizes $(basename -s .bdg $bdg_file).bw
# 
# ##remove conversion file
# rm conversion_file.bdg
# done
# ```
####################################################################

#back to R
#Make bed file with bins for all cells
#
#```{r get_bins, echo=TRUE}
#Extract binarized counts
binarized_counts <- Matrix::colSums(x.sp@bmat)

#Scale scoring for bed file using the maximum n of reads in a bin and create a name for the bins
scaled_binarized_counts <- binarized_counts * (1000/max(binarized_counts))
names_binarized_counts <- paste0(binarized_counts, "of", length(x.sp@barcode))

#Create bed datafame with regions of bins extracted from x.sp@feature and scores based on sum of all cells (colSums of x.sp@bmat)
bed_df <- data.frame(seqnames=seqnames(x.sp@feature),
                     starts=start(x.sp@feature)-1,
                     ends=end(x.sp@feature),
                     names=names_binarized_counts,
                     scores=scaled_binarized_counts,
                     strands=strand(x.sp@feature))

#Rename strands
levels(bed_df$strands)[levels(bed_df$strands)=="*"] <- "."

#Write out table
write.table(bed_df, file=paste0("LAKO.", gsub(" ", "_", clusters.sel)[i], "_binary_matrix.bed"), quote=FALSE, sep= "\t", row.names=FALSE, col.names=FALSE)
#```

#Make bed file with bins for seperate clusters
#```{r get_bins2, echo=TRUE}
#Create variable with cluster names
clusters.sel = names(table(x.sp@cluster))

#Create bed_df for every cluster
for (cluster_name in seq(clusters.sel)) {
  binarized_counts_cluster <- Matrix::colSums(x.sp@bmat[which(x.sp@cluster==clusters.sel[cluster_name]),])
  
  #Scale scoring for bed file using the maximum n of reads in a bin and create a name for the bins
  scaled_binarized_counts <- binarized_counts_cluster * (1000/max(binarized_counts))
  names_binarized_counts <- paste0(binarized_counts_cluster, "of", length(x.sp@barcode))
  
  #Create bed datafame with regions of bins extracted from x.sp@feature and scores based on sum of all cells (colSums of x.sp@bmat)
  bed_df <- data.frame(seqnames=seqnames(x.sp@feature),
                       starts=start(x.sp@feature)-1,
                       ends=end(x.sp@feature),
                       names=names_binarized_counts,
                       scores=scaled_binarized_counts,
                       strands=strand(x.sp@feature))
  #Rename strands
  levels(bed_df$strands)[levels(bed_df$strands)=="*"] <- "."
  
  #Write out table
  write.table(format(bed_df, scientific=FALSE), file=paste0(resultsdir, "_cluster", cluster_name, "_binary_matrix.bed"), quote=FALSE, sep= "\t", row.names=FALSE, col.names=FALSE)
  
}
#```

#Peak matrix generation
#```{r get_peaks, echo=TRUE, results="hide"}
#Set working directory to directory containing narrowpeaks of clusters
# with session

#Create GRanges object with all peak locations from clusters
peaks.names = system("ls | grep narrowPeak", intern = TRUE)
peak.gr.ls = lapply(peaks.names, function(x){
  peak.df = read.table(x)
  GRanges(peak.df[,1], IRanges(peak.df[,2], peak.df[,3]))
})
peak.gr = suppressWarnings(reduce(Reduce(c, peak.gr.ls)))
##2 times reduce? Do I need this.
#Create a cell-by-peak matrix
peaks.df = as.data.frame(peak.gr)[,1:3]

# delete wierd characters
# peaks.df$seqnames<-gsub("b'","",as.character(peaks.df$seqnames))
# peaks.df
# peaks.df$seqnames<-gsub("'","",as.character(peaks.df$seqnames))
# peaks.df

write.table(peaks.df, file=paste0(resultsdir, "_peaks_combined.bed"), append = FALSE, quote=FALSE, sep="\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE, qmethod = c("escape", "double"), fileEncoding = "")

# delete wierd characters
peaks.df$seqnames<-gsub("b'","",as.character(peaks.df$seqnames))
peaks.df
peaks.df$seqnames<-gsub("'","",as.character(peaks.df$seqnames))
peaks.df

write.table(peaks.df, file=paste0(resultsdir, "_peaks_combined_nochar.bed"), append = FALSE, quote=FALSE, sep="\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE, qmethod = c("escape", "double"), fileEncoding = "")

#Save the current snap object into an RDS
x.sp@feature # seems right
saveRDS(x.sp, file=paste0(resultsdir, ".snap.rds"))

#Set variable for bash chunk

#delete below?
snap_files <- c("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/data/lako2021/snap_files/")
Sys.setenv(snap_name=paste0(resultsdir, snap_files))
Sys.setenv(bed_name=paste0(resultsdir, "_peaks_combined.bed"))

# ```{bash addpmat, echo=TRUE}
# source activate /mbshome/jverploegen/anaconda3/envs/snaptools_py2
# snaptools snap-del --snap-file $snap_name --session-name PM
# snaptools snap-add-pmat --snap-file $snap_name --peak-file $bed_name
# ```

#```{r addpmat2, echo=TRUE}
x.sp = readRDS(x.sp, file= "/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/20210518.snap.rds")
x.sp = addPmatToSnap(x.sp,do.par = T,num.cores = 10)
x.sp@pmat

#Create directory and save unbinarized cell by peak matrix for scRNA integration (part3)
scRNA_int_dir <- paste0(resultsdir, "scRNA_integration/")
dir.create(scRNA_int_dir)

#Get peak matrix and add peak locations as column names, transpose
peak_matrix <- x.sp@pmat
peak_matrix@Dimnames[2] =list(x.sp@peak@elementMetadata@listData$name)
peak_matrix <- t(peak_matrix)
#Save for later use with Seurat integration with scRNA-seq
saveRDS(peak_matrix, file=paste0(scRNA_int_dir, "cell_peak_matrix.rds"))

#Extract metadata and save for later use
ATAC_metadata <- x.sp@metaData
write.csv(ATAC_metadata, paste0(scRNA_int_dir, "ATAC_metadata.csv"), row.names=FALSE)

#Binarize pmat for further usage in snap object
x.sp = makeBinary(x.sp, mat="pmat")
#```

#```{r find_DARs, echo=TRUE}
# Define the differentially accessible regions for the clusters
idy.ls = lapply(levels(x.sp@cluster), function(cluster_i){
  DARs = findDAR(
    obj=x.sp,
    input.mat="pmat",
    cluster.pos=cluster_i,
    cluster.neg=NULL,
    cluster.neg.method="knn",
    test.method="exactTest",
    bcv=0.01, #Only one plate, so bcv = 0.01
    seed.use=10
  )
  DARs$FDR = p.adjust(DARs$PValue, method="BH")
  idy = which(DARs$FDR < 0.05 & DARs$logFC > 0)
  #For small clusters, we here select the top 2000 peaks
  if((x=length(idy)) < 2000L){
    PValues = DARs$PValue
    PValues[DARs$logFC < 0] = 1
    idy = order(PValues, decreasing=FALSE)[1:2000]
    rm(PValues); # free memory
  }
  idy
})

covs = Matrix::rowSums(x.sp@pmat)
names(idy.ls) = levels(x.sp@cluster)

pdf(paste(resultsdir,'DAR_clusters.pdf',sep=""),width=6,height=6,paper='special') 
for(cluster_i in levels(x.sp@cluster)){
  idy = idy.ls[[cluster_i]]
  vals = Matrix::rowSums(x.sp@pmat[,idy]) / covs
  vals.zscore = (vals - mean(vals)) / sd(vals)
  plotFeatureSingle(
    obj=x.sp,
    feature.value=vals.zscore,
    method="umap",
    main=paste0("Cluster ", cluster_i, " z-score"),
    point.size=0.5,
    point.shape=19,
    down.sample=5000,
    quantiles=c(0.01, 0.99)
  )
}
dev.off()
#```

#```{r annotation_info}
#Create a list with all the directories and file names for part 2 and 3
folder_info <- list("working_dir" = workdir, 
                    "results_folder" = resultsdir)

#Save the snap object for further analysis (part2, fluff, gimme)
saveRDS(x.sp, file=paste0(resultsdir, ".snap.rds"))

#Save the folder_info list for further analysis
saveRDS(folder_info, paste0(workdir, "_folder_locations.rds"))

######################################################################## GO TO SEURAT ENV for determining clusters from scRNA-seq
library(Seurat)

lako.rna = readRDS("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq/lakoRNAannotated.rds")
variable.genes = VariableFeatures(object = lako.rna) 

################################################################ 
#return to snapatac env
## reload the bmat, this is optional but highly recommanded
library(SnapATAC)
library(GenomicRanges)
x.sp = readRDS("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/20210609.snap.rds")

genes.df <- read.table("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/gencode.v30.gene.annotation.gtf") #already ran this
genes.gr3 = GRanges(genes.df[,1], IRanges(genes.df[,2], genes.df[,3]), name=genes.df[,4])
genes.sel.gr3 = genes.gr3[which(genes.gr3$name %in% variable.genes)]

#x.sp = addBmatToSnap(x.sp)
x.sp@bmat

x.sp3 = createGmatFromMat(
  obj=x.sp, 
  input.mat="bmat",
  genes=genes.sel.gr3,
  do.par=TRUE#,
  #num.cores=10
)

x.sp3

### convert snap to seurat object
lakoall.atac <- snapToSeurat(obj=x.sp3,
                             eigs.dims=1:10, 
                             norm=TRUE,
                             scale=TRUE
)

saveRDS(lakoall.atac,file="/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/lakoatacall.rds")

################################################################
### from here on different environment kb_scrna_R_seurat4
library(Seurat)
lako.rna = readRDS("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq/lakoRNAannotated.rds")
lako.rna@meta.data
#variable.genes = VariableFeatures(object = lako.rna) 
#variable.genes2 = variable.genes

lakoall.atac = readRDS("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/lakoatacall.rds")
lakoall.atac@meta.data$cluster

transfer.anchors <- FindTransferAnchors(
  reference = lako.rna, 
  query = lakoall.atac, 
  features = variable.genes, 
  reference.assay = "RNA", 
  query.assay = "ACTIVITY", 
  reduction = "cca", 
    
)

lako.rna$costum_clustering

celltype.predictions <- TransferData(
  anchorset = transfer.anchors, 
  refdata = lako.rna$costum_clustering,
  weight.reduction = lakoall.atac[["SnapATAC"]],
  dims = 1:10
)

lakoall.atac@meta.data$predict.id = celltype.predictions$predicted.id;
lakoall.atac@meta.data$predict.max.score = apply(celltype.predictions[,-1], 1, max);

saveRDS(lakoall.atac,file="/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/lakoatacall.rds")

refdata <- GetAssayData(
   object = lako.rna, 
   assay = "RNA", 
   slot = "data"
 )
imputation <- TransferData(
   anchorset = transfer.anchors, 
refdata = refdata, 
   weight.reduction = lakoall.atac[["SnapATAC"]], 
   dims = 1:10
 )

saveRDS(imputation,file="/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/imputation.rds")
rm(lako.rna)

################################### go back to snapatac env
library(SnapATAC)
library(Seurat)
library(Matrix)
# reload the gmat upon the original x.sp (not the seurat converted file)
x.sp = readRDS("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/20210609.snap.rds")
#imputation = readRDS("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/imputation.rds") # does not work due to seurat boject missing

genes.df <- read.table("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/gencode.v30.gene.annotation.gtf") #already ran this
genes.gr3 = GRanges(genes.df[,1], IRanges(genes.df[,2], genes.df[,3]), name=genes.df[,4])
genes.sel.gr3 = genes.gr3[which(genes.gr3$name %in% variable.genes)]

#x.sp = addBmatToSnap(x.sp)
x.sp@bmat

x.sp = createGmatFromMat(
  obj=x.sp, 
  input.mat="bmat",
  genes=genes.sel.gr3,
  do.par=TRUE#,
  #num.cores=10
)

#x.sp@gmat = t(imputation@data);

### DELETE
#celltype.predictions$predicted.id <- gsub("Mel","MEC", celltype.predictions$predicted.id)

x.sp@metaData$predict.id = celltype.predictions$predicted.id;
x.sp@metaData$predict.max.score = apply(celltype.predictions[,-1], 1, max);


pdf(paste(resultsdir,'predictscore_all.pdf',sep=""),width=18,height=18,paper='special') 
hist(
  x.sp@metaData$predict.max.score, 
  xlab="prediction score", 
  col="lightblue", 
  xlim=c(0, 1),
  main="Lako 10X"
)
abline(v=0.5, col="red", lwd=2, lty=2)
table(x.sp@metaData$predict.max.score > 0.5)
dev.off()

x.sp2 = x.sp[x.sp@metaData$predict.max.score > 0.2,];
x.sp2@metaData
x.sp2@metaData[,"predict.id"]

# own annotation of the clusters
current.cluster.ids <-  c("cell1","cell2","cell3","cell4","cell5","cell6")
new.cluster.ids <- c( "cell1","cornealcells","cell3","cell4","stromalcells","LPCs")
# 
x.sp@cluster <- plyr::mapvalues(
  x = x.sp@cluster, 
  from = current.cluster.ids, 
  to = new.cluster.ids
)

x.sp@cluster

saveRDS(x.sp,file="/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/lakopredicted.rds")

x.sp <- readRDS("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/lakopredicted.rds")

x.sp3 = x.sp[x.sp@metaData$predict.max.score > 0.2,];
x.sp3@metaData
x.sp3@metaData[,"predict.id"]

# cross referencing barcodes from own annotated clusters and imputated clusters for QC
# x.sp2[x.sp2@metaData$predict.id == "CjB",]
# x.sp2[x.sp2@cluster == "cell1" & x.sp2@metaData$predict.id == "CjB",]
# 
# # selecting LNPCs
# x.sp2[x.sp2@metaData$predict.id == "LNPCs",]
# x.sp2[x.sp2@cluster == "cell3" & x.sp2@metaData$predict.id == "LNPCs",]
# 
# # multi selection
# x.sp2
# x.sp4 <- x.sp2[(x.sp2@cluster == "cell1" & x.sp2@metaData$predict.id == "CjB") | 
#         (x.sp2@cluster == "cell3" & x.sp2@metaData$predict.id == "LNPCs") |
#         (x.sp2@cluster == "cell2" & x.sp2@metaData$predict.id == "CjS") |
#         (x.sp2@cluster == "cell1" | x.sp2@cluster == "cell3" & x.sp2@metaData$predict.id == "CB") |
#         (x.sp2@cluster == "cell4" & x.sp2@metaData$predict.id == "Mel" | x.sp2@metaData$predict.id == "IC") |
#         (x.sp2@cluster == "cell6" & x.sp2@metaData$predict.id == "LPCs") |
#         (x.sp2@cluster == "cell5" & x.sp2@metaData$predict.id == "Ves" | x.sp2@metaData$predict.id == "FCECs" | x.sp2@metaData$predict.id == "LSC" | x.sp2@metaData$predict.id == "CSSCs"),]

library(umap)
pdf(paste(resultsdir,'visualisation_imputated_UMAP_final.pdf',sep=""),width=10,height=8,paper='special') 
plotViz(obj=x.sp2, method="umap",text.size = 1.5, point.color=x.sp2@metaData$predict.id, legend.add = T)
plotViz(obj=x.sp, method="umap",text.size = 1.5, point.color=x.sp@cluster, legend.add = T)
plotViz(obj=x.sp3, method="umap",text.size = 1.5, point.color=x.sp3@metaData$predict.id, legend.add = T)
plotViz(obj=x.sp4, method="umap",text.size = 1.5, point.color=x.sp4@metaData$predict.id, legend.add = T)
dev.off()

x.sp4[x.sp4@metaData$predict.id == "LNPCs",]
#########################################################################################################################
#########################################################################################################################
# generating bigwig tracks from imputated data
#########################################################################################################################
#########################################################################################################################

# script from Rebecca and Jelle

# Create narrowpeak per cluster
#Create a narrowpeak file for each of the clusters using MACS2
#```{r pseudo_bulk_cluster, echo = TRUE, results="hide"}
#Create directory for narrowpeak files of cluster

dir.create(paste0(resultsdir, "clustersimputated/"))
setwd(paste0(resultsdir, "clustersimputated/"))

#Get amount of cells per cluster:
#summary(x.sp@metaData$cluster) for all 8 earlier found
x.sp3@metaData$predict.id
SnowParam(workers = 2)
# call peaks for all cluster with more than 80 cells
clusters.sel = names(table(x.sp3@metaData$predict.id))[which(table(x.sp3@metaData$predict.id) > 80)]

system("which snaptools")
system("which macs2")
SnowParam(workers = 6)

peaks.ls = mclapply(seq(clusters.sel), function(i){
  print(clusters.sel[i])
  peaks = runMACS(
    obj=x.sp3[which(x.sp3@metaData$predict.id==clusters.sel[i]),], 
    output.prefix=paste0("LAKOIMP02.", gsub(" ", "_", clusters.sel)[i]),
    path.to.snaptools="/vol/mbconda/julian/envs/kb_snapatac/bin/snaptools",
    path.to.macs="/vol/mbconda/julian/envs/kb_snapatac/bin/macs2",
    gsize="hs", # mm, hs, etc
    buffer.size=500, 
    num.cores=1,
    macs.options="--nomodel --shift 100 --ext 200 --qval 5e-2 -B --SPMR",
    tmp.folder="/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/data/lako2021/R_peaks_scATAC"
  )
  peaks
})#, mc.cores=5)

# peaks.ls = mclapply(seq(clusters.sel), function(i){
#   print(clusters.sel[i])
#   runMACS(
#     obj=x.sp[which(x.sp@cluster==clusters.sel[i]),],
#     output.prefix=paste0(pseudo_bulk_name, "_cluster", clusters.sel[i]),
#     path.to.snaptools="/vol/mbconda/julian/envs/kb_snapatac/bin/snaptools",
#     path.to.macs = "/vol/mbconda/julian/envs/kb_snapatac/bin/macs2",
#     gsize="hs",
#     buffer.size=500,
#     num.cores=1,
#     macs.options="--nomodel --shift 37 --ext 73 --qval 5e-2 -B --SPMR --call-summits",
# #macs.options="--nomodel --shift 100 --ext 200 --qval 5e-2 -B --SPMR", earlier paramets
#     tmp.folder=tempdir()
#   )
# }, mc.cores=5)
#Set parameter to create bigwigs
Sys.setenv(bdg_location=paste0("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/peaks_scATAC_20210614/"))

# Settings v3: --nomodel --shift 37 --ext 73 --qval 5e-2 -B --SPMR --call-summits
# Settings BAMPE (v2): "--format BAMPE --qval 5e-2 -B --SPMR --call-summits"
# Settings Jos (change macs): "--nomodel --shift 100 --ext 200 --qval 5e-2 -B --SPMR"
# Settings Original: "--nomodel --shift 37 --ext 73 --qval 1e-2 -B --SPMR --call-summits"
#``

###############################################
#continue from below
#Create bigwigs
#Create a bigwig from the processed data based on the bedgraph produced above. Folder containing date should be changed to the correct name.

#################################################################
## bash commands
#```{bash create_bigwigs, echo= TRUE}
#source activate /mbshome/jverploegen/anaconda3/envs/genome_tools

# ##set directory where output bedgraph is located
# cd $bdg_location
# 
#
# For all bdg files: exclude GL000*, KI* and chrM (not in the hg38.chrom.sizes file) and remove for each cluster the b and '.
# grep -Ev 'GL000|KI|chrM' $bdg_file > $ bdg_file2
# nice -30 sed -i ' s/b//g' $bdg_file2 && sed -i " s/[']//g" $bdg_file2
# 
# conda activate bigwig
#
# ##Loop to create bigwig from all bedgraphs
# for bdg_file in *bdg
# do
# ##sort bedgraph and remove random chromosomes etc.
# sort -k1,1 -k2,2n $bdg_file > conversion_file.bdg
# 
# nice -30 bedGraphToBigWig conversion_file.bdg /ceph/rimlsfnwi/data/moldevbio/zhou/jarts/cellrangeratac/hg38.chrom.sizes $(basename -s .bdg $bdg_file).bw
# 
# ##remove conversion file
# rm conversion_file.bdg
# done
# ```
####################################################################

#back to R
#Make bed file with bins for all cells
#
#```{r get_bins, echo=TRUE}
#Extract binarized counts
binarized_counts <- Matrix::colSums(x.sp3@bmat)

#Scale scoring for bed file using the maximum n of reads in a bin and create a name for the bins
scaled_binarized_counts <- binarized_counts * (1000/max(binarized_counts))
names_binarized_counts <- paste0(binarized_counts, "of", length(x.sp@barcode))

#Create bed datafame with regions of bins extracted from x.sp@feature and scores based on sum of all cells (colSums of x.sp@bmat)
bed_df <- data.frame(seqnames=seqnames(x.sp3@feature),
                     starts=start(x.sp3@feature)-1,
                     ends=end(x.sp3@feature),
                     names=names_binarized_counts,
                     scores=scaled_binarized_counts,
                     strands=strand(x.sp3@feature))

#Rename strands
levels(bed_df$strands)[levels(bed_df$strands)=="*"] <- "."

#Write out table
for (i in clusters.sel){
write.table(bed_df, file=paste0("LAKOIMP02.", gsub(" ", "_", clusters.sel[i]), "_binary_matrix.bed"), quote=FALSE, sep= "\t", row.names=FALSE, col.names=FALSE)
#```
}

#Make bed file with bins for seperate clusters
#```{r get_bins2, echo=TRUE}


#Create bed_df for every cluster
for (cluster_name in seq(clusters.sel)) {
  binarized_counts_cluster <- Matrix::colSums(x.sp@bmat[which(x.sp@cluster==clusters.sel[cluster_name]),])
  
  #Scale scoring for bed file using the maximum n of reads in a bin and create a name for the bins
  scaled_binarized_counts <- binarized_counts_cluster * (1000/max(binarized_counts))
  names_binarized_counts <- paste0(binarized_counts_cluster, "of", length(x.sp@barcode))
  
  #Create bed datafame with regions of bins extracted from x.sp@feature and scores based on sum of all cells (colSums of x.sp@bmat)
  bed_df <- data.frame(seqnames=seqnames(x.sp3@feature),
                       starts=start(x.sp3@feature)-1,
                       ends=end(x.sp3@feature),
                       names=names_binarized_counts,
                       scores=scaled_binarized_counts,
                       strands=strand(x.sp3@feature))
  #Rename strands
  levels(bed_df$strands)[levels(bed_df$strands)=="*"] <- "."
  
  #Write out table
  write.table(format(bed_df, scientific=FALSE), file=paste0(resultsdir, "_cluster", cluster_name, "_binary_matrix.bed"), quote=FALSE, sep= "\t", row.names=FALSE, col.names=FALSE)
  
}
#```

#Peak matrix generation
#```{r get_peaks, echo=TRUE, results="hide"}
#Set working directory to directory containing narrowpeaks of clusters
# with session

#Create GRanges object with all peak locations from clusters
peaks.names = system("ls | grep narrowPeak", intern = TRUE)
peak.gr.ls = lapply(peaks.names, function(x){
  peak.df = read.table(x)
  GRanges(peak.df[,1], IRanges(peak.df[,2], peak.df[,3]))
})
peak.gr = suppressWarnings(reduce(Reduce(c, peak.gr.ls)))
##2 times reduce? Do I need this.
#Create a cell-by-peak matrix
peaks.df = as.data.frame(peak.gr)[,1:3]

# delete wierd characters
# peaks.df$seqnames<-gsub("b'","",as.character(peaks.df$seqnames))
# peaks.df
# peaks.df$seqnames<-gsub("'","",as.character(peaks.df$seqnames))
# peaks.df

write.table(peaks.df, file=paste0(resultsdir, "_peaks_combined.bed"), append = FALSE, quote=FALSE, sep="\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE, qmethod = c("escape", "double"), fileEncoding = "")

# delete wierd characters
peaks.df$seqnames<-gsub("b'","",as.character(peaks.df$seqnames))
peaks.df
peaks.df$seqnames<-gsub("'","",as.character(peaks.df$seqnames))
peaks.df

write.table(peaks.df, file=paste0(resultsdir, "_peaks_combined_nochar.bed"), append = FALSE, quote=FALSE, sep="\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE, qmethod = c("escape", "double"), fileEncoding = "")

#Save the current snap object into an RDS
x.sp3@feature # seems right
saveRDS(x.sp3, file=paste0(resultsdir, ".snap.rds"))

#Set variable for bash chunk

#delete below?
snap_files <- c("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/data/lako2021/snap_files/")
Sys.setenv(snap_name=paste0(resultsdir, snap_files))
Sys.setenv(bed_name=paste0(resultsdir, "_peaks_combined.bed"))

# ```{bash addpmat, echo=TRUE}
# source activate /mbshome/jverploegen/anaconda3/envs/snaptools_py2
# snaptools snap-del --snap-file $snap_name --session-name PM
# snaptools snap-add-pmat --snap-file $snap_name --peak-file $bed_name
# ```

#```{r addpmat2, echo=TRUE}
SnowParam(workers = 6)
x.sp3 = addPmatToSnap(x.sp3,do.par = T,num.cores = 6)
x.sp3@pmat

saveRDS(x.sp3,file="/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/lakopredicted_wpmat.rds")

x.sp3 <- readRDS("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/lakopredicted_wpmat.rds")

#Create directory and save unbinarized cell by peak matrix for scRNA integration (part3)
scRNA_int_dir <- paste0(resultsdir, "scRNA_integration_all/")
dir.create(scRNA_int_dir)

tracks_int_dir <- paste0(resultsdir, "tracks_integration_all/")
split_int_dir <- paste0("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/split_integration_all/")
dir.create(tracks_int_dir)

#Get peak matrix and add peak locations as column names, transpose
peak_matrix <- x.sp3@pmat
peak_matrix@Dimnames[2] =list(x.sp3@peak@elementMetadata@listData$name)
peak_matrix <- t(peak_matrix)



###################################################

#Binarize pmat for further usage in snap object
x.sp3 = makeBinary(x.sp3, mat="pmat")

# add new labels for correct DAG calculation
x.sp3@metaData$predict.id

current.cluster.ids <- c('CB','CjB',
                         'CjS','LPCs','IC','LNPCs','Mel',
                         'Ves','CSSCs',
                         'LSC','FCECs')
new.cluster.ids <- c('CB','CSB',
                     'CjS','LPCs','IC','LNPCs','MEC',
                     'Ves','CSSCs',
                     'StC','FCECs')
x.sp3@metaData$predict.id
x.sp3@metaData$predict.id <- plyr::mapvalues(x = x.sp3@metaData$predict.id, from = current.cluster.ids, to = new.cluster.ids)

# for filtering the cells
x.sp4 <- x.sp3[x.sp3@metaData$predict.id != "IC" & x.sp3@metaData$predict.id != "FCECs" & x.sp3@metaData$predict.id != "Ves",]
x.sp4@cluster <- as.factor(x.sp4@metaData$predict.id)
x.sp4@cluster

# for filtering

# generating a cluster tree
ensemble.ls = lapply(split(seq(length(x.sp4@cluster)), x.sp4@cluster), function(x){
  SnapATAC::colMeans(x.sp4[x,], mat="bmat");
})

# cluster using 1-cor as distance  
hc = hclust(as.dist(1 - cor(t(do.call(rbind, ensemble.ls)))), method="ward.D2")
plot(hc, hang=-1, xlab="")

pdf(paste(resultsdir,'/cluster_tree_all_imputated_filt.pdf',sep=""),width=6,height=6,paper='special') 
plot(hc, hang=-1, xlab="")
dev.off()
############################################################
#Save for later use with Seurat integration with scRNA-seq
saveRDS(peak_matrix, file=paste0(scRNA_int_dir, "cell_peak_matrix.rds"))

#Extract metadata and save for later use
ATAC_metadata <- x.sp3@metaData
write.csv(ATAC_metadata, paste0(scRNA_int_dir, "ATAC_metadata.csv"), row.names=FALSE)

############################################################
# generate barcode files for each cell population to be used in "bamslice" for each sample
x.sp3@sample
ATAC_metadata$sample <- x.sp3@sample
ATAC_metadata$sample
for (z in unique(ATAC_metadata$sample)){
  for (i in unique(ATAC_metadata$predict.id)){
    local <- ATAC_metadata[i==ATAC_metadata$predict.id & ATAC_metadata$sample == z,]
    write.csv(local$barcode, paste0(tracks_int_dir, i,"_", z,"_metadata.csv"), row.names=F) 
  }
}  

# split the dataset randomly based upon each sample; less clean
for (z in unique(ATAC_metadata$sample)){
  for (i in unique(ATAC_metadata$predict.id)){
    dflako <- ATAC_metadata[i==ATAC_metadata$predict.id & ATAC_metadata$sample == z,]
    ind <- sample(c(TRUE, FALSE), replace=TRUE, prob=c(0.5, 0.5))
    data1 <- dflako[ind, ]
    data2 <- dflako[!ind, ]
    write.csv(data1$barcode, paste0(resultsdir, i,"_", z,"_metadata1.csv"), row.names=F) 
    write.csv(data2$barcode, paste0(resultsdir, i,"_", z,"_metadata2.csv"), row.names=F) 
  }
}  
for (z in unique(ATAC_metadata$sample)){
  for (i in unique(ATAC_metadata$predict.id)){
    dflako <- ATAC_metadata[i==ATAC_metadata$predict.id & ATAC_metadata$sample == z,]
    #dflako$sampname <- sample(factor(rep(1:3, length.out=nrow(dflako)), 
     #                                labels=paste0("samp", 1:3)))
    #print(dflako$sampname)
    #data1 <- dflako[dflako$sample == "atac_lako2",]
    #data2 <- dflako[dflako$sample == "atac_lako3",]
    #data3 <- dflako[dflako$sample == "atac_lako4",]
    write.csv(dflako$barcode, paste0(scRNA_int_dir, i,"_", z,"_metadata1.csv"), row.names=F) 
    #write.csv(data2$barcode, paste0(scRNA_int_dir, i,"_", z,"_metadata2.csv"), row.names=F) 
    #write.csv(data2$barcode, paste0(scRNA_int_dir, i,"_", z,"_metadata3.csv"), row.names=F) 
  }
}  
# randomly split dataset into three equal amount of barcodes (best way for vst normalization)
for (i in unique(ATAC_metadata$predict.id)){
  dflako <- ATAC_metadata[i==ATAC_metadata$predict.id,]
  dflako$rep <- sample(factor(rep(1:3, length.out=nrow(dflako)), labels=paste0("rep", 1:3)))
  for (z in unique(dflako$rep)){
    write.csv(dflako$barcode, paste0(resultsdir, i,"_", z,"_metadata.csv"), row.names=F) 
  }
}

# split the dataset into blob
ATAC_metadata <- x.sp4@metaData
unique(ATAC_metadata$predict.id)

for (i in unique(ATAC_metadata$predict.id)){
  dflako <- ATAC_metadata[i!=ATAC_metadata$predict.id,]
  write.table(dflako$barcode, paste0(resultsdir,"/", i,"_blob_metadata.csv"), row.names=F,col.names = F,quote = F) 
}

# split the dataset into blob (INCLUDE A NORMALISATOIN AND RANDOMIZATION FACTOR HERE)
ATAC_metadata <- x.sp4@metaData
unique(ATAC_metadata$predict.id)

for (i in unique(ATAC_metadata$predict.id)){
  dflako <- ATAC_metadata[i!=ATAC_metadata$predict.id,]
  write.table(dflako$barcode, paste0(resultsdir,"/", i,"_blob_metadata.csv"), row.names=F,col.names = F,quote = F) 
}

########################################################################################
# Define the differentially accessible regions for the clusters
idy.ls = lapply(levels(x.sp3@cluster), function(cluster_i){
  DARs = findDAR(
    obj=x.sp3,
    input.mat="pmat",
    cluster.pos=cluster_i,
    cluster.neg=NULL,
    cluster.neg.method="knn",
    test.method="exactTest",
    bcv=0.01, #Only one plate, so bcv = 0.01
    seed.use=10
  )
  DARs$FDR = p.adjust(DARs$PValue, method="BH")
  idy = which(DARs$FDR < 0.05 & DARs$logFC > 0)
  #For small clusters, we here select the top 2000 peaks
  if((x=length(idy)) < 2000L){
    PValues = DARs$PValue
    PValues[DARs$logFC < 0] = 1
    idy = order(PValues, decreasing=FALSE)[1:2000]
    rm(PValues); # free memory
  }
  idy
})

covs = Matrix::rowSums(x.sp3@pmat)
names(idy.ls) = levels(x.sp3@cluster)

pdf(paste(resultsdir,'DAR_clusters.pdf',sep=""),width=6,height=6,paper='special') 
for(cluster_i in levels(x.sp3@cluster)){
  idy = idy.ls[[cluster_i]]
  vals = Matrix::rowSums(x.sp3@pmat[,idy]) / covs
  vals.zscore = (vals - mean(vals)) / sd(vals)
  plotFeatureSingle(
    obj=x.sp3,
    feature.value=vals.zscore,
    method="umap",
    main=paste0(cluster_i, " z-score"),
    point.size=0.5,
    point.shape=19,
    down.sample=5000,
    quantiles=c(0.01, 0.99)
  )
}
dev.off()

# make pmat binairy for motif analysis in R if you like (see tutorial)