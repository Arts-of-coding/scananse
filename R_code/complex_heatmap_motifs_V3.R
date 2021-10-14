###########################################################################
# Workflow COMPLEX HEATMAP OF MOTIFS AND RNA-seq

# Loading in important packages
library(chorddiag)
library(circlize)
library(viridis)
library(ggvenn)
#install.packages('venn') #run this the first time when initializing
library(venn)
library(dplyr)
library(tidyr)
library(tidyselect)
library(ComplexHeatmap)

# Setting up the working directory
workdir <- "/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/complexheatmap_motifs/"

# Setting up results directory
dateoftoday <- gsub("-", "", as.character(Sys.Date()))
resultsdir <- paste0(workdir, dateoftoday)
system(paste("mkdir -p ", resultsdir))

# Read in transcription factor annotation files if they are available
expected <- read.csv("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/data/lako2021/ANANSE/Expected_tfs.csv",header = T, comment.char = '#')
exgeneral <- read.csv("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/data/lako2021/ANANSE/Expectedgeneral_tfs.csv",header = T, comment.char = '#')

# Load in the Z-score normalized dataset of gene expression:
lakorna <- read.table(file = '/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/20210712/Zscoretable.tsv', sep = '\t', header = TRUE, row.names=1) #for all individual pops

# For filtering; the three clusters were not used in motif analysis, so excluded from gene expression
lakorna <- lakorna[, !names(lakorna) %in% c("IC", "FCECs", "Ves")]

# Setting the colors of the general transcription factors expected to be found
eyeexpected <- unlist(strsplit(exgeneral[exgeneral$cell_type == "eye",]$expected_tfs,","))
epiexpected <- unlist(strsplit(exgeneral[exgeneral$cell_type == "epidermal",]$expected_tfs,","))

lakorna <- na.omit(lakorna)

###########################################################################
# Load in the q-quantile normalized files of the motifs where one TF is associated with one motif

lakoqnormmotifs <- read.table(file = "/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/data/lako2021/motif_analysis/quantile_sigregions/joined.txt",sep = ' ', header = TRUE, row.names=1)
lakoqnormmotifs <- lakoqnormmotifs[,1:12]
lakoqnormmotifs<- lakoqnormmotifs[, !colnames(lakoqnormmotifs) %in% c("IC", "FCECs", "Ves")]

# Select only one motif for each gene
selector <- unique(lakoqnormmotifs[c("genes")])
names <- rownames(selector)
y <- subset(lakoqnormmotifs, rownames(lakoqnormmotifs) %in% names)

# Subselect overlapping genes in the two datasets
z <- subset(lakorna, rownames(lakorna) %in% (x$genes))
z <- subset(lakorna, rownames(lakorna) %in% (y$genes))
y <- subset(y, y$genes %in% rownames(z))

rownames(y) <- y$genes
y$genes <- NULL

###########################################################################
# Converting the dataframes to matrices for complex heatmap of all TFs

matz <- as.matrix(z)
maty <- as.matrix(y)

# Selecting the top 10 factors of all populations
vec2 <- NULL
for (i in colnames(y)) {
sel <- maty[order(maty[,i],decreasing = T),]
vec <- rownames(sel)[1:10]
print(vec)
vec2 <- c(vec2,vec)
}
vec2 <- unique(vec2)

z <- subset(z, rownames(z) %in% vec2)
y <- subset(y, rownames(y) %in% vec2)

###########################################################################
# Converting the dataframes to matrices for complex heatmap of only the top 10 TFs

matz <- as.matrix(z)
maty <- as.matrix(y)

# Setting color ramps for the two heatmaps
f1 = colorRamp2(c(-3, 0, 3), c("blue", "#EEEEEE", "red"), space = "RGB")
f2 = colorRamp2(c(-5, 0, 5), c("blue", "#EEEEEE", "yellow"), space = "RGB")


###########################################################################
# Making the complex heatmap of gene expression

ht1 = Heatmap(matz, col = f1, cluster_columns = T,cluster_rows = T, name = "Z-score RNA count")

# Order the rows based on heatmap clustering in both matrices
roword <- row_order(ht1)
matz <- matz[roword,]
colord <- column_order(ht1)
matz <- matz[, colord]
maty <- maty[rownames(matz),]

# Setting up the vector that will give the names at the side of the heatmap
vec2 <- rownames(matz)

###########################################################################
# Making the complex heatmap of motifs

ht2 = Heatmap(maty, col = f2, cluster_columns = F,cluster_rows = F, name = "Z-score qnorm motifs ")


###########################################################################
# Plotting both heatmaps and colors if they are found in literature

pdf(paste(resultsdir,'/complexheatmap_motifs.pdf',sep="/") ,width=15,height=8,paper='special')

# Color annotation of the TFs
vec3 <- which(rownames(maty) %in% vec2, arr.ind=TRUE)
fontcolors <- rep('black', length(vec3))

# Linking the celltype specific factors
row_idx <- which(vec2 %in% eyeexpected)
fontcolors[row_idx] <- 'red'

# linking the celltype specific factors
row_idx <- which(vec2 %in% epiexpected)
fontcolors[row_idx] <- 'orange'

htlist <- ht1 + ht2 + rowAnnotation(link = anno_mark(at =  vec3,labels = unique(vec2),labels_gp = gpar(col = fontcolors)))
draw(htlist, column_title = "RNAmotifscomparison")
dev.off()
