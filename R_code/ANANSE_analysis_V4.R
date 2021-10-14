# Loading in important packages
library(chorddiag)
library(circlize)
library(viridis)
library(ggvenn)
#install.packages('venn')
library(venn)
library(dplyr)
library(tidyr)
library(tidyselect)
#devtools::install_github("hms-dbmi/UpSetR")
library(UpSetR)
library(ComplexHeatmap)
library(ggplot2)
library(reshape)
#install.packages("dtwclust")
library(dtwclust)
library(clustree)
#####################################
# Storing results
workdir <- "/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/chordESC"

# Setting up results directory
dateoftoday <- gsub("-", "", as.character(Sys.Date()))
resultsdir <- paste0(workdir, dateoftoday)
system(paste("mkdir -p ", resultsdir))

# Read in the files to analyse with chord
files <- read.csv("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/data/lako2021/ANANSE/chord_files.csv",header = T, comment.char = '#') 

# Read in transcription factor annotation files if they are available
expected <- read.csv("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/data/lako2021/ANANSE/Expected_tfs.csv",header = T, comment.char = '#')
exgeneral <- read.csv("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/data/lako2021/ANANSE/Expectedgeneral_tfs.csv",header = T, comment.char = '#')

# Generate the interaction dataframe
factorsheat1 <- NULL
fullheat <- list()
for (i in unique(files$cell_id)){
  factorsheat1 <- read.table(files$network[files$cell_id==i],header = T)
  factorsheat2 <- factorsheat1[,c("factor","targetScaled")]
  colnames(factorsheat2) <- c("factor",i)
  fullheat[[i]] <- factorsheat2
}

merged.data.frame = Reduce(function(...) merge(..., all=T), fullheat)
merged.data.frame[merged.data.frame == 0] <- NA

merged.data.frame2 <- merged.data.frame
merged.data.frame2[!is.na(merged.data.frame2)] <- 1
merged.data.frame2[is.na(merged.data.frame2)] <- 0
merged.data.frame2$factor <- merged.data.frame$factor
merged.data.frame2 <-merged.data.frame2[rowSums(merged.data.frame2[, -1])>0, ]

for(i in 2:ncol(merged.data.frame2)){ merged.data.frame2[ , i] <- as.integer(merged.data.frame2[ , i]) }
colnoms <- colnames(merged.data.frame2)[colnames(merged.data.frame2) != "factor"]

fill = viridis(length(files$cell_id),alpha = 0.8,option = "H")

pdf(paste(resultsdir,'/intersectionblob.pdf',sep="/") ,width=10,height=5,paper='special')
upset(merged.data.frame2,nintersects = NA,
      sets = colnoms, 
      order.by="freq", matrix.color="black", point.size=1,
      sets.bar.color=fill)
dev.off()

# add this!: set_order = c("a", "b", "c")

# Extracting unique factors
lstnames <- c()
fullheatfacts= list()
for (i in unique(files$cell_id)){
  lstnames <- c(lstnames,paste0(merged.data.frame2$i))
  factorsheat3 <- read.table(files$network[files$cell_id==i],header = T)
  factorsheat4 <- factorsheat3[,"factor"]
  fullheatfacts[[i]] <- factorsheat4
}

# Extracting the shared factors of all populations
sharedall <- Reduce(intersect, fullheatfacts)

###################
# Extracting unique factors for each cell population to list
uniq <- list()
for (v in unique(files$cell_id)){
  fullheatfacts2 <- fullheatfacts
  fullheatfacts2[[v]] <- NULL
  fac_interest <- fullheatfacts[[v]]
  `%notin%` <- Negate(`%in%`)
  uniq[[v]] <- fac_interest[fac_interest %notin% unlist(fullheatfacts2)]
}

###############################################
# Generating the complex heatmaps quantitation for all unique factors + shared factors and influence scores

## load in the z-score normalized dataset:
lakorna <- read.table(file = '/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/20210712/Zscoretable.tsv', sep = '\t', header = TRUE, row.names=1) #for all individual pops
#lakorna <- read.table(file = '/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/Z-score20210825/Zscoretable.tsv', sep = '\t', header = TRUE, row.names=1)

lakorna <- na.omit(lakorna)

# for filtering
lakorna <- lakorna[, !names(lakorna) %in% c("IC", "FCECs", "Ves")]

# Combine all lists into one giant dataframe
merged.data.frame = Reduce(function(...) merge(..., all=T), fullheat)

heatdfshared <- subset(merged.data.frame, merged.data.frame$factor %in% sharedall)
uniqtfs <- unlist(uniq)
heatdfuniq <- subset(merged.data.frame, merged.data.frame$factor %in% uniqtfs)
rownames(heatdfuniq) <- heatdfuniq$factor
heatdfuniq$factor <- NULL

rownames(heatdfshared) <- heatdfshared$factor
heatdfshared$factor <- NULL

# Say influence score is 0 for uniq factors that have NA
heatdfuniq[is.na(heatdfuniq)] <- 0

matuniq <- as.matrix(heatdfuniq)
matshared <- as.matrix(heatdfshared)

f1 = colorRamp2(c(0, 1), c("white","purple"), space = "RGB")
f2 = colorRamp2(c(0, 1), c("white", "red"), space = "RGB")

ht1 = Heatmap(matuniq, col= f1,cluster_rows = T, cluster_columns = T, name = "influence score")
ht2 = Heatmap(matshared, col= f1,cluster_rows = T, cluster_columns = T, name = "influence score")

dev.off()
pdf(paste(resultsdir,'/complexheatmap.pdf',sep="/") ,width=15,height=8,paper='special')
ht1
ht2
dev.off()

# Generating a dataframe from influence scores
heatdfall <- merged.data.frame
rownames(heatdfall) <- heatdfall$factor
heatdfall$factor <- NULL
heatdfall[is.na(heatdfall)] <- 0

# Selecting a cutoff for the influence values
maxvec <- unname(apply(heatdfall[,], 1, max))
heatdfall$maxscore <- maxvec
heatdfall= heatdfall[heatdfall$maxscore > 0.5,]
heatdfall$maxscore <- NULL

# Selecting tfs based on the cutoff in the Z-score table
z <- subset(lakorna, rownames(lakorna) %in% (merged.data.frame$factor))
z <- z[rownames(z) %in% rownames(heatdfall),]

# Generate the matrix 
matall <- as.matrix(heatdfall)
matz <- as.matrix(z)

ht4 = Heatmap(matall, col= viridis(100), cluster_rows = T, cluster_columns = T, name = "Influence score")

# Order the rows and columns based on Z-score gene expression table
roword <- row_order(ht4)
matall <- matall[roword,]
colord <- column_order(ht4)
matall <- matall[, colord]
matz <- matz[rownames(matall),colnames(matall)]

col_mat = viridis(3,alpha = 0.8)
f3 = colorRamp2(c(-3, 0, 3), c("blue", "#EEEEEE","red"), space = "RGB")
ht3 = Heatmap(matz, col= f3, cluster_rows = F, cluster_columns = F, name = "Z-score RNA count")

vec3 <- sharedall[sharedall %in% rownames(heatdfall)]
vec3 <- rownames(matz[rownames(matz) %in% vec3 ,])

vec4 <- unname(uniqtfs[uniqtfs %in% rownames(heatdfall)])
vec4 <- rownames(matz[rownames(matz) %in% vec4 ,])

vec1 <- which(rownames(matz) %in% vec3, arr.ind = T)
vec2 <- which(rownames(matz) %in% vec4, arr.ind = T)

# Setting the colors of the general transcription factors expected to be found
#neurexpected <- unlist(strsplit(exgeneral[exgeneral$cell_type == "neural",]$expected_tfs,","))
eyeexpected <- unlist(strsplit(exgeneral[exgeneral$cell_type == "eye",]$expected_tfs,","))
epiexpected <- unlist(strsplit(exgeneral[exgeneral$cell_type == "epidermal",]$expected_tfs,","))

if (!length(vec1) == 0) {
  # linking the celltype specific factors
  fontcolors <- rep('black', length(vec3))
  #row_idx <- which(vec3 %in% neurexpected)
  #fontcolors[row_idx] <- 'purple'
  
  # linking the celltype specific factors
  row_idx <- which(vec3 %in% eyeexpected)
  fontcolors[row_idx] <- 'red'
  
  # linking the celltype specific factors
  row_idx <- which(vec3 %in% epiexpected)
  fontcolors[row_idx] <- 'orange'
  
  # linking the expected factors as a blue color
  expected2 <- unlist(strsplit(expected[expected$cell_type == i,]$expected_tfs,","))
  
  # linking the celltype specific factors
  row_idx <- which(vec3 %in% expected2)
  fontcolors[row_idx] <- 'blue'
  htlist1 <- ht4 + ht3 + rowAnnotation(link = anno_mark(at =  vec1,labels = vec3,labels_gp = gpar(col = fontcolors)))
}

# linking the expected factors as a blue color for the single cell populations
expected2 <- unlist(strsplit(expected$expected_tfs,","))

# linking the celltype specific factors
fontcolors <- rep('black', length(vec4))
#row_idx <- which(vec4 %in% neurexpected)
#fontcolors[row_idx] <- 'purple'

# linking the celltype specific factors
row_idx <- which(vec4 %in% eyeexpected)
fontcolors[row_idx] <- 'red'

# linking the celltype specific factors
row_idx <- which(vec4 %in% epiexpected)
fontcolors[row_idx] <- 'orange'

# linking the celltype specific factors
row_idx <- which(vec4 %in% expected2)
fontcolors[row_idx] <- 'blue'

htlist2 <- ht4 + ht3 + rowAnnotation(link = anno_mark(at =  vec2,labels = vec4,labels_gp = gpar(col = fontcolors)))

pdf(paste(resultsdir,'/complexheatmap_rna.pdf',sep="/") ,width=10,height=10,paper='special')
if (!length(vec1) == 0) {
  draw(htlist1, column_title = "shared_factors")
}
draw(htlist2, column_title = "unique_factors")
dev.off()

# Top 25 selection and intersection

# Selecting subpopulations
cellids <- c("LNPCs","LPCs","CB","CSB","CjS")

for (i in unique(files$cell_id)){
  print(i)
  #diffnetwork <- read.table(files$diffnetwork[files$cell_id==i])
  factors <- read.table(files$network[files$cell_id==i],header = T)
  factors <- factors[order(factors$sumScaled,decreasing = T),]
  # selecting the top 25 factors
  top <- factors$factor[1:25]
  vec5 <- top
  vec5 <- vec5[vec5 %in% rownames(heatdfall)]
  print(vec5)
  vec25 <- c(vec5,vec5)
}

# Selecting subpopulations of interest
#cellids <- c("CSSCs","StC")

for (i in cellids){
  print(i)
  #diffnetwork <- read.table(files$diffnetwork[files$cell_id==i])
  factors <- read.table(files$network[files$cell_id==i],header = T)
  factors <- factors[order(factors$sumScaled,decreasing = T),]
  # selecting the top 25 factors
  top <- factors$factor[1:25]
  vec5 <- top
  vec5 <- vec5[vec5 %in% rownames(heatdfall)]
  print(vec5)
  vec25 <- c(vec5,vec5)
}

vec25 <- unique(vec25)

# Subselecting the matrix
matall25 <- matall[rownames(matall) %in% unique(vec25),colnames(matall) %in% cellids]
ht7 = Heatmap(matall25, col= viridis(100), cluster_rows = T, cluster_columns = T, name = "Influence score")

# Order the rows and columns based on Z-score gene expression table
matz25 <- matz[rownames(matz) %in% unique(vec25),colnames(matz) %in% cellids]

roword <- row_order(ht7)
matall25 <- matall25[roword,]
colord <- column_order(ht7)
matall25 <- matall25[, colord]
matz25 <- matz25[rownames(matall25),colnames(matall25)]

###
col_mat = viridis(3,alpha = 0.8)
f3 = colorRamp2(c(-3, 0, 3), c("blue", "#EEEEEE","red"), space = "RGB")
ht8 = Heatmap(matz25, col= f3, cluster_rows = F, cluster_columns = F, name = "Z-score RNA count")

vec6 <- rownames(matz25)
vec5 <- which(rownames(matz25) %in% vec6, arr.ind = T)

# linking the celltype specific factors
fontcolors <- rep('black', length(vec6))
#row_idx <- which(vec6 %in% neurexpected)
#fontcolors[row_idx] <- 'purple'

# linking the celltype specific factors
row_idx <- which(vec6 %in% eyeexpected)
fontcolors[row_idx] <- 'red'

# linking the celltype specific factors
row_idx <- which(vec6%in% epiexpected)
fontcolors[row_idx] <- 'orange'
htlist3 <- ht7 + ht8 + rowAnnotation(link = anno_mark(at =  vec5,labels = vec6,labels_gp = gpar(col = fontcolors)))

pdf(paste(resultsdir,'/complexheatmap_rna25_corneal.pdf',sep="/") ,width=10,height=5,paper='special')
draw(htlist3, column_title = "shared_factors_top25")
dev.off()

# Single the cell populations for a clear overview
pdf(paste(resultsdir,'/complexheatmap_rna_singled.pdf',sep="/") ,width=10,height=10,paper='special')
for (i in unique(files$cell_id)) {
  vec5 <- uniq[[i]]
  vec5 <- vec5[vec5 %in% rownames(heatdfall)]
  if(length(vec5) == 0) {
    next
  }
  vec5 <- rownames(matz[rownames(matz) %in% vec5 ,])
  vec6 <- which(rownames(matz) %in% vec5, arr.ind = T)
  
  # linking the expected factors as a blue color
  expected2 <- unlist(strsplit(expected[expected$cell_type == i,]$expected_tfs,","))
  
  # linking the celltype specific factors
  fontcolors <- rep('black', length(vec6))
  #row_idx <- which(vec6 %in% neurexpected)
  #fontcolors[row_idx] <- 'purple'
  
  # linking the celltype specific factors
  row_idx <- which(vec6 %in% eyeexpected)
  fontcolors[row_idx] <- 'red'
  
  # linking the celltype specific factors
  row_idx <- which(vec6 %in% epiexpected)
  fontcolors[row_idx] <- 'orange'
  
  # linking the celltype specific factors
  row_idx <- which(vec6 %in% expected2)
  fontcolors[row_idx] <- 'blue'
  
  htlist2 <- ht4 + ht3 + rowAnnotation(link = anno_mark(at =  vec6,labels = vec5,labels_gp = gpar(col = fontcolors)))
  draw(htlist2, column_title = i)
}
dev.off()

# Generate a heatmap of the top 20 found factors for each cell population
pdf(paste(resultsdir,'/', paste0('complexheatmap_top25.pdf'),sep="/") ,width=10,height=10,paper='special')
for (i in unique(files$cell_id)){
  print(i)
  diffnetwork <- read.table(files$diffnetwork[files$cell_id==i])
  factors <- read.table(files$network[files$cell_id==i],header = T)
  factors <- factors[order(factors$sumScaled,decreasing = T),]
  # selecting the top 25 factors
  top <- factors$factor[1:25]
  vec5 <- top
  vec5 <- vec5[vec5 %in% rownames(heatdfall)]
  if(length(vec5) != 0) {
    vec5 <- rownames(matz[rownames(matz) %in% vec5 ,])
    vec6 <- which(rownames(matz) %in% vec5, arr.ind = T)
    
    # linking the expected factors as a blue color
    expected2 <- unlist(strsplit(expected[expected$cell_type == i,]$expected_tfs,","))
    
    # linking the celltype specific factors
    fontcolors <- rep('black', length(vec5))
    #row_idx <- which(vec5 %in% neurexpected)
    #fontcolors[row_idx] <- 'purple'
    
    # linking the celltype specific factors
    row_idx <- which(vec5 %in% eyeexpected)
    fontcolors[row_idx] <- 'red'
    
    # linking the celltype specific factors
    row_idx <- which(vec5 %in% epiexpected)
    fontcolors[row_idx] <- 'orange'
    
    # linking the celltype specific factors
    row_idx <- which(vec5 %in% expected2)
    fontcolors[row_idx] <- 'blue'
    
    # plotting the list
    htlist2 <- ht4 + ht3 + rowAnnotation(link = anno_mark(at =  vec6,labels = vec5,labels_gp = gpar(col = fontcolors)))
    draw(htlist2, column_title = paste0(i,' top 25'))
  }
}
dev.off()



##############################
# Extracting unique factors shared across populations of interest
list_of_cells = c("LPCs","CSB","CjS")

list_of_vecs <- list()
fullheatfacts2 <- fullheatfacts
for (v in list_of_cells){
  list_of_vecs[[v]]<-fullheatfacts[[v]]
  fullheatfacts2[[v]] <- NULL
}

fac_interest <- Reduce(intersect, list_of_vecs)
`%notin%` <- Negate(`%in%`)

vecname <- paste("int_",paste(list_of_cells, collapse = ''),sep="")
assign(vecname,fac_interest[fac_interest %notin% unlist(fullheatfacts2)])
facs <- fac_interest[fac_interest %notin% unlist(fullheatfacts2)]

# Generate a heatmap of interesting factors
pdf(paste(resultsdir,'/', paste0(vecname,'complexheatmap_interesting.pdf'),sep="/") ,width=10,height=10,paper='special')
vec5 <- facs
vec5 <- vec5[vec5 %in% rownames(heatdfall)]
if(length(vec5) != 0) {
  vec5 <- rownames(matz[rownames(matz) %in% vec5 ,])
  vec6 <- which(rownames(matz) %in% vec5, arr.ind = T)
  
  # linking the expected factors as a blue color
  expected2 <- unlist(strsplit(expected$expected_tfs,","))
  
  # linking the celltype specific factors
  fontcolors <- rep('black', length(vec5))
  #row_idx <- which(vec5 %in% neurexpected)
  #fontcolors[row_idx] <- 'purple'
  
  # linking the celltype specific factors
  row_idx <- which(vec5 %in% eyeexpected)
  fontcolors[row_idx] <- 'red'
  
  # linking the celltype specific factors
  row_idx <- which(vec5 %in% epiexpected)
  fontcolors[row_idx] <- 'orange'
  
  # linking the celltype specific factors
  row_idx <- which(vec5 %in% expected2)
  fontcolors[row_idx] <- 'blue'
  htlist2 <- ht4 + ht3 + rowAnnotation(link = anno_mark(at =  vec6,labels = vec5,labels_gp = gpar(col = fontcolors)))
  draw(htlist2, column_title = vecname)
}
dev.off()

# Subselecting interesting cells
colsub <- c("LPCs","LNPCs","CjS","CSB","CB","CSSCs","StC")
list_of_cells1 = c("LPCs","LNPCs","CjS","CSB","CB")
list_of_cells4 = c("LPCs","CjS","CSB","CB")
list_of_cells5 = c("LNPCs","CjS","CSB","CB")
list_of_cells6 = c("LNPCs","CB")
list_of_cells3 = c("LPCs","CjS")
list_of_cells9 = c("LPCs","CB")
list_of_cells10 = c("LNPCs","CjS")
list_of_cells2 = c("CjS","CSB","CB")
list_of_cells7 = c("LNPCs","LPCs")
list_of_cells8 = c("LNPCs","LPCs","CB")
list_of_cells11 = c("LNPCs","LPCs","CjS")
list_of_cells12 = c("CSSCs","StC")
list_of_cells13 = c("LPCs","CSB","CjS")
list_of_cells14 = c("LNPCs","CB","CSB")

#list_of_cells1 = c("LiCo")
#list_of_cells2 = c("StCSC")
#list_of_cells3 = c("CSSCs","StC")

veclst <- list(list_of_cells1,list_of_cells2,list_of_cells3,list_of_cells4,list_of_cells5,list_of_cells6,
               list_of_cells7,list_of_cells8,list_of_cells9,list_of_cells10,list_of_cells11,list_of_cells12,list_of_cells13,list_of_cells14)#,

facs2 <- NULL

for (i in veclst){
  print(i)
  fullheatfacts2 <- fullheatfacts
  list_of_vecs <- list()
  for (v in i){
    print(v)
    list_of_vecs[[v]]<-fullheatfacts[[v]]
    fullheatfacts2[[v]] <- NULL
  }
  
  fac_interest <- Reduce(intersect, list_of_vecs)
  `%notin%` <- Negate(`%in%`)
  
  vecname <- paste("int_",paste(i, collapse = ''),sep="")
  assign(vecname,fac_interest[fac_interest %notin% unlist(fullheatfacts2)])
  facs <- fac_interest[fac_interest %notin% unlist(fullheatfacts2)]
  facs2 <- c(facs2,facs)
}

# Add all expected factors in the subselection
#expected2 <- unlist(strsplit(expected$expected_tfs,","))
#facs2 <- c(facs2,epiexpected,eyeexpected,neurexpected,expected2)

# Selecting a cutoff for the influence values
heatdfall2 <- heatdfall[rownames(heatdfall) %in% (facs2),]

# subselecting columns based on interesting cells
heatdfall2 <- heatdfall2[, colnames(heatdfall2) %in% colsub]


# Selecting tfs based on the cutoff in the Z-score table
#z <- subset(lakorna, rownames(lakorna) %in% (merged.data.frame$factor))
#zz <- z[rownames(z) %in% rownames(heatdfall2),]
zz <- z[rownames(z) %in% rownames(heatdfall2),]

# subselecting columns based on interesting cells
zz <- zz[, colnames(zz) %in% colsub]

# Generate the matrix 
#matall2 <- as.matrix(heatdfall2)
matall2 <- as.matrix(heatdfall2)
matz2 <- as.matrix(zz)

col.order <- c("LPCs","LNPCs","CB","CSB","CjS","CSSCs","StC")
matall2 <- matall2[,col.order]

ht5 = Heatmap(matall2, col= viridis(100), cluster_rows = T, cluster_columns = F, name = "Influence score")

# Order the rows and columns based on Z-score gene expression table
roword <- row_order(ht5)
matall2 <- matall2[roword,]
colord <- column_order(ht5)
matall2 <- matall2[, colord]
matz2 <- matz2[rownames(matall2),colnames(matall2)]

col_mat = viridis(3,alpha = 0.8)
f3 = colorRamp2(c(-3, 0, 3), c("blue", "#EEEEEE","red"), space = "RGB")
ht6 = Heatmap(matz2, col= f3, cluster_rows = F, cluster_columns = F, name = "Z-score RNA count")

vec3 <- rownames(matz2)
vec1 <- which(rownames(matz2) %in% vec3, arr.ind = T)

#pdf(paste(resultsdir,'/complexheatmap_rna_subselected.pdf',sep="/") ,width=10,height=13.5,paper='special')
pdf(paste(resultsdir,'/complexheatmap_rna_subselected3.pdf',sep="/") ,width=10,height=10,paper='special')
if (!length(vec1) == 0) {
  # linking the celltype specific factors
  fontcolors <- rep('black', length(vec3))
  #row_idx <- which(vec3 %in% neurexpected)
  #fontcolors[row_idx] <- 'purple'
  
  # linking the celltype specific factors
  row_idx <- which(vec3 %in% eyeexpected)
  fontcolors[row_idx] <- 'red'
  
  # linking the celltype specific factors
  row_idx <- which(vec3 %in% epiexpected)
  fontcolors[row_idx] <- 'orange'
  
  # linking the celltype specific factors
  row_idx <- which(vec3 %in% expected2)
  fontcolors[row_idx] <- 'red'
  htlist1 <- ht5 + ht6 + rowAnnotation(link = anno_mark(at =  vec1,labels = vec3,labels_gp = gpar(col = fontcolors)))
}
draw(htlist1, column_title = "subselected")#, row_km = 2, cluster_rows = TRUE)
dev.off()

# Line plots of interesting factors together (WIP)
datalineplot <- NULL
# classifications of time series data

datalineplot2 <- NULL
datalineplot <-as.data.frame(as.table(matall))
matall4 <- matall2
for (res in 4:20){
dtw_cluster = tsclust(heatdfall, type="partitional",k=res,
                      distance="dtw_basic",centroid = "pam",seed=1234,trace=T,
                      args = tsclust_args(dist = list(window.size = 5)))


dfall3 <- NULL
matall3 <- as.data.frame(unname(dtw_cluster@cluster),rownames(matall))

datalineplot[[paste0("dtw", res)]]<- rep(matall3[rownames(matall3)%in%datalineplot$Var1,], length(unique(datalineplot$Var2)))


}


pdf(paste(resultsdir,'/grouped_tfs_pseudotime_clustree.pdf',sep="/") ,width=10,height=20,paper='special')
datalineplot2 <- datalineplot[1:46,4:16]

print(clustree(datalineplot2,prefix = "dtw"))
print(clustree(
  datalineplot2,
  prefix = "dtw",
  #exprs = c("data", "counts", "scale.data"),
  assay = NULL,
  node_colour = "sc3_stability"
))
dev.off()

########### corresponding expression in scRNA-seq of the factors (violin plots)
# function for stacked plots was based on https://divingintogeneticsandgenomics.rbind.io/post/stacked-violin-plot-for-visualizing-single-cell-data-in-seurat/
library(Seurat)
library(patchwork)
library(ggplot2)

unique(datalineplot$Var1)
## remove the x-axis text and tick
## plot.margin to adjust the white space between each plot.
## ... pass any arguments to VlnPlot in Seurat
modify_vlnplot<- function(obj, 
                          feature, 
                          pt.size = 0, 
                          plot.margin = unit(c(-0.05, 0, -0.05, 0), "cm"),
                          ...) {
  p<- VlnPlot(obj, features = feature, pt.size = pt.size, ... )  + 
    xlab("") + ylab(feature) + ggtitle("") + 
    theme(legend.position = "none", 
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(), 
          axis.title.y = element_text(size = rel(1), angle = 0), 
          axis.text.y = element_text(size = rel(1)), 
          plot.margin = plot.margin ) 
  return(p)
}

## extract the max value of the y axis
extract_max<- function(p){
  ymax<- max(ggplot_build(p)$layout$panel_scales_y[[1]]$range$range)
  return(ceiling(ymax))
}


## main function
StackedVlnPlot<- function(obj, features,
                          pt.size = 0, 
                          plot.margin = unit(c(-0.05, 0, -0.05, 0), "cm"),
                          ...) {
  
  plot_list<- purrr::map(features, function(x) modify_vlnplot(obj = obj,feature = x, ...))
  
  # Add back x-axis title to bottom plot. patchwork is going to support this?
  plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
    theme(axis.text.x=element_text(), axis.ticks.x = element_line())
  
  # change the y-axis tick to only max value 
  ymaxs<- purrr::map_dbl(plot_list, extract_max)
  plot_list<- purrr::map2(plot_list, ymaxs, function(x,y) x + 
                            scale_y_continuous(breaks = c(y)) + 
                            expand_limits(y = y))
  
  p<- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
  return(p)
}

seur_obj <- readRDS('/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq/lakoRNAannotated.rds')
seur_obj2 <-subset(x=seur_obj, subset = (costum_clustering == "LPCs"|costum_clustering == "LNPCs"|costum_clustering == "CB"|costum_clustering == "CSB"|costum_clustering == "CjS"))
seur_obj2$costum_clustering
seur_obj2 <- droplevels(x = seur_obj2)
seur_obj2$costum_clustering <- droplevels(x = seur_obj2$costum_clustering)

features<- unique(datalineplot$Var1)

seur_obj2@meta.data$seurat_clusters <- seur_obj2$costum_clustering
seur_obj2@active.ident <- seur_obj2$costum_clustering

new_order <- c("LPCs","LNPCs","CB","CSB","CjS")
seur_obj2@active.ident <- factor(seur_obj2@active.ident, levels = new_order)

# Selected 13 clusters from clustree
pdf(paste(resultsdir,'/grouped_tfs_pseudotime.pdf',sep="/") ,width=6,height=8,paper='special')
for (i in 1:13){
datalineplot5 <- datalineplot[datalineplot$dtw13==i,]
print(ggplot(datalineplot5, aes(x=Var2, y=Freq, group=Var1)) +
  geom_line(aes(color=Var1))+
  geom_point() +
  labs(color = "Transcription factors",y = "Influence score",x= "Ident"))
  print(StackedVlnPlot(obj = seur_obj2, features = unique(datalineplot5$Var1)))
}
dev.off()

#### interesting subselection according to clustering
int_nfkb <- c("NFKB1","RELA","ESRRA","HES4")

pdf(paste(resultsdir,'/grouped_tfs_NFKB.pdf',sep="/") ,width=6,height=8,paper='special')
datalineplot5 <- datalineplot[datalineplot$Var1 %in%int_nfkb,]
print(ggplot(datalineplot5, aes(x=Var2, y=Freq, group=Var1)) +
          geom_line(aes(color=Var1))+
          geom_point() +
          labs(color = "Transcription factors",y = "Influence score",x= "Ident"))
print(StackedVlnPlot(obj = seur_obj2, features = unique(datalineplot5$Var1)))
dev.off()


#int_cornea <- c("ZBTB33","ELK4","PAX6","USF1","DPB","KLF5","FOSL2","EHF","ETV2","GRHL1","ELF3","GRHL2") first approach; all intersections
int_cornea <- c("KLF5","EHF","ELF3","GRHL1","HSF4","HES4","KLF6","FOXA1","MAFK","MAFF") #first approach; all intersections

pdf(paste(resultsdir,'/grouped_tfs_cornea.pdf',sep="/") ,width=6,height=5,paper='special') 
datalineplot5 <- datalineplot[(datalineplot$Var1 %in%int_cornea),]#& (datalineplot$Var2 %notin% c("CSSCs","StC"))
datalineplot5 <- datalineplot5[datalineplot5$Var2 %in% new_order,]
print(ggplot(datalineplot5, aes(x=factor(Var2, level = new_order), y=Freq,group = Var1)) +
        geom_line(aes(color=Var1))+
        geom_point() +
        labs(color = "Transcription factors",y = "Influence score",x= "Ident"))
dev.off()

datalineplot5 <- datalineplot5[datalineplot5$Var2 %in% new_order,]
#neworder2 <- rep(new_order,5)
#group_by(.data = datalineplot5,datalineplot5$Var2,new_order)
datalineplot5[match(neworder2, datalineplot5$Var2),]
pdf(paste(resultsdir,'/grouped_tfs_cornea_expr.pdf',sep="/") ,width=6,height=12,paper='special') 
print(StackedVlnPlot(obj = seur_obj2, features = unique(datalineplot5$Var1)))
dev.off()

fontcolors <- rep('black', length(datalineplot$Var1))
#row_idx <- which(datalineplot$Var1 %in% neurexpected)
#fontcolors[row_idx] <- 'purple'

# linking the celltype specific factors
row_idx <- which(datalineplot$Var1 %in% eyeexpected)
fontcolors[row_idx] <- 'red'

# linking the celltype specific factors
row_idx <- which(datalineplot$Var1 %in% epiexpected)
fontcolors[row_idx] <- 'orange'

# linking the celltype specific factors
row_idx <- which(datalineplot$Var1 %in% expected2)
fontcolors[row_idx] <- 'blue'

datalineplot$type <- fontcolors

datalineplot1 <- datalineplot[datalineplot$type =="red",]
ggplot(datalineplot1, aes(x=Var2, y=Freq, group=Var1)) +
  geom_line(aes(color=Var1))+
  geom_point() +
  labs(color = "Transcription factors")

datalineplot2 <- datalineplot[datalineplot$type =="blue",]
ggplot(datalineplot2, aes(x=Var2, y=Freq, group=Var1)) +
  geom_line(aes(color=Var1))+
  geom_point() +
  labs(color = "Transcription factors")

datalineplot3 <- datalineplot[datalineplot$type =="black",]
ggplot(datalineplot3, aes(x=Var2, y=Freq, group=Var1)) +
  geom_line(aes(color=Var1))+
  geom_point() +
  labs(color = "Transcription factors")

datalineplot4 <- datalineplot[datalineplot$type =="orange",]
ggplot(datalineplot4, aes(x=Var2, y=Freq, group=Var1)) +
  geom_line(aes(color=Var1))+
  geom_point() +
  labs(color = "Transcription factors")




StackedVlnPlot(obj = seur_obj2, features = features)


#####################################
#rating chord diagrams of the top 10 factors
# note: perhaps based on the shared and unique factors?

########################## WIP
# subselecting factors based on the plot (by hand first)


circos.par(gap.degree=0.5)

myFun <- function(vector) {
  ind <- ave(rep(1, length(vector)), vector, FUN = length)
  print(ind)
  thresh <- max(ind)
  vector[ind > thresh - 1] ## added "+1" to match your terminology
}

top <- c("ZBTB33","ELK4","PAX6","FOXQ1","TFAP2C","OTX1","IRF6","E2F4","ELK1","TGIF1","RFX5","TP63","ASCL2")
length(top)
cells <- c("LPCs","LNPCs","CB","CjS","CSB")

# add a list in the code for the upset R plot
fullheat <- list()

pdf(paste(resultsdir,"/",cells,'chordsselected.pdf',sep=""),width=6,height=6,paper='special') 
for (i in cells){
  print(i)
  diffnetwork <- read.table(files$diffnetwork[files$cell_id==i])
  #factors <- read.table(files$network[files$cell_id==i],header = T)
  #factors <- factors[order(factors$sumScaled,decreasing = T),]
  
  # selecting the top 20 factors
  #top <- factors$factor[1:20]
  
  colnames(diffnetwork) <- c("from","to","value")
  df <- data.frame(diffnetwork)
  
  # selecting
  df3 <- df[df$from %in% top,]
  df3 <- df3[order(df3$value,decreasing = T) & (df3$value >0.8),]
  vec <- df3$to
  vec2 <- unique(myFun(vec))
  df3 <- df3[(df3$to %in% vec2),]
  vec3 <- df3$from
  
  #ind <- ave(rep(1, length(vec)), vec, FUN = length)
  #thresh <- max(ind)
  
  vec4 <- unique(myFun(vec3))
  df3 <- df3[(df3$from %in% vec4),]
  #vec <- df3$to
  #vec2 <- unique(myFun(vec,length(top)))
  # for (z in unique(top)){
  #   print(z)
  #   df4 <- df3
  #   df5 <- df4[df4$from == z,]
  #   df6 <- df5[df5$to == head(df5$to,10),]
  #   df7 = NULL
    
    # Determine if one of the top 10 factors also influences the factor of interest
    # for (q in unique(head(df5$to,10))){
    #   print(q)
    #   df7 = rbind(df7, df4[df4$from == q & df4$to == z,])
    #   df7 = rbind(df7, df4[df4$from == q & df4$to %in% head(df5$to,10),])
    #   
    #   # If a factor is found influencing the factor of interest or another factor in the top 10 then add to original dataframe
    #   if (!is.null(df7)){
    #     df6 = rbind(df6, df7)#[df4$from == q & df4$to %in% head(df5$to,10) | df4$to == z,])
    #   }
    # }
    
    # making the dataframe consist of only unique rows
    df6 <- unique(df3)
    
    # Setting the right color scale
    u1 <- unique(unlist(df6[, c("from", "to")]))
    
    col_mat = viridis(length(u1),alpha = 0.8,option = "H")
    
    chordDiagram(df3, directional = 1, direction.type = c("diffHeight", "arrows"),
                 link.arr.type = "big.arrow", grid.col = col_mat)
    title(paste(i , "joint regulation factors"))
  #}
}

dev.off()

install.packages("alluvial")
library(alluvial)



df3$freq <- sum(df3$from)
alluvial(df3[,1:2], freq = df3$value,col=ifelse( df3$from == "PAX6" & df3$to == df3$to, "red", "gray") )#, freq=df3$from, border=NA)

## alluvial from top 10s
top <- c("ZBTB33","ELK4","PAX6","FOXQ1","TFAP2C","OTX1","IRF6","E2F4","ELK1","TGIF1","RFX5","TP63","ASCL2")
length(top)
cells <- c("LPCs","LNPCs","CB","CjS","CSB")

# add a list in the code for the upset R plot
fullheat <- list()

pdf(paste(resultsdir,"/",cells,'chordsselected.pdf',sep=""),width=6,height=6,paper='special') 
for (i in cells){
  print(i)
  diffnetwork <- read.table(files$diffnetwork[files$cell_id==i])
  #factors <- read.table(files$network[files$cell_id==i],header = T)
  #factors <- factors[order(factors$sumScaled,decreasing = T),]
  
  # selecting the top 20 factors
  #top <- factors$factor[1:20]
  
  colnames(diffnetwork) <- c("from","to","value")
  df <- data.frame(diffnetwork)
  
  # selecting
  df3 <- df[df$from %in% top,]
  df3 <- df3[order(df3$value,decreasing = T) & (df3$value >0.8),]
  vec <- df3$to
  vec2 <- unique(myFun(vec))
  df3 <- df3[(df3$to %in% vec2),]
  vec3 <- df3$from
  
  #ind <- ave(rep(1, length(vec)), vec, FUN = length)
  #thresh <- max(ind)
  
  vec4 <- unique(myFun(vec3))
  df3 <- df3[(df3$from %in% vec4),]
  #vec <- df3$to
  #vec2 <- unique(myFun(vec,length(top)))
  # for (z in unique(top)){
  #   print(z)
  #   df4 <- df3
  #   df5 <- df4[df4$from == z,]
  #   df6 <- df5[df5$to == head(df5$to,10),]
  #   df7 = NULL
  
  # Determine if one of the top 10 factors also influences the factor of interest
  # for (q in unique(head(df5$to,10))){
  #   print(q)
  #   df7 = rbind(df7, df4[df4$from == q & df4$to == z,])
  #   df7 = rbind(df7, df4[df4$from == q & df4$to %in% head(df5$to,10),])
  #   
  #   # If a factor is found influencing the factor of interest or another factor in the top 10 then add to original dataframe
  #   if (!is.null(df7)){
  #     df6 = rbind(df6, df7)#[df4$from == q & df4$to %in% head(df5$to,10) | df4$to == z,])
  #   }
  # }
  
  # making the dataframe consist of only unique rows
  df6 <- unique(df3)
  
  # Setting the right color scale
  u1 <- unique(unlist(df6[, c("from", "to")]))
  
  col_mat = viridis(length(u1),alpha = 0.8,option = "H")
  
  chordDiagram(df3, directional = 1, direction.type = c("diffHeight", "arrows"),
               link.arr.type = "big.arrow", grid.col = col_mat)
  title(paste(i , "joint regulation factors"))
  #}
}

chordDiagram(df3)

### UPSET intersect of factors in different cell types  WIP
factorsheat1 <- NULL
fullheat <- list()
for (i in unique(files$cell_id)){
  factorsheat1 <- read.table(files$network[files$cell_id==i],header = T)
  factorsheat2 <- factorsheat1[,c("factor","targetScaled")]
  colnames(factorsheat2) <- c("factor",i)
  fullheat[[i]] <- factorsheat2
}

merged.data.frame = Reduce(function(...) merge(..., all=T), fullheat)
merged.data.frame[merged.data.frame == 0] <- NA

merged.data.frame2 <- merged.data.frame
merged.data.frame2[!is.na(merged.data.frame2)] <- 1
merged.data.frame2[is.na(merged.data.frame2)] <- 0
merged.data.frame2$factor <- merged.data.frame$factor
merged.data.frame2 <-merged.data.frame2[rowSums(merged.data.frame2[, -1])>0, ]

for(i in 2:ncol(merged.data.frame2)){ merged.data.frame2[ , i] <- as.integer(merged.data.frame2[ , i]) }
colnoms <- colnames(merged.data.frame2)[colnames(merged.data.frame2) != "factor"]

fill = viridis(length(files$cell_id),alpha = 0.8,option = "H")

pdf(paste(resultsdir,'/intersectionblob.pdf',sep="/") ,width=10,height=5,paper='special')
upset(merged.data.frame2,nintersects = NA,
      sets = colnoms, 
      order.by="freq", matrix.color="black", point.size=1,
      sets.bar.color=fill)
dev.off()

########################################################3




df5 <- NULL
cells <- c("LPCs","LNPCs")#,"CB","CjS","CSB")
col_mat = viridis(length(cells),alpha = 0.8,option = "H")
dfcol <- as.data.frame(col_mat,cells)

for (i in cells){
  print(i)
  diffnetwork <- read.table(files$diffnetwork[files$cell_id==i])
  factors <- read.table(files$network[files$cell_id==i],header = T)
  factors <- factors[order(factors$sumScaled,decreasing = T),]
  # selecting the top factors or iteresting factors
  #top <- factors$factor[1:20]
  top <- c("PAX6","TP63","TGIF1","SP3","ELK4","ZBTB33")
  colnames(diffnetwork) <- c("from","to","value")
  df <- data.frame(diffnetwork)
  df <- df[df$value>0.80,]
  df3 <- df[df$from %in% top,]
  df3 <- df3[order(df3$value,decreasing = T),]
  
  # Selecting the top 20 target genes
  top2 <- df3$to[1:20]
  df4 <- df3[df3$to %in% top2,]
  df4$cell <- i
  df4$col <- dfcol[rownames(dfcol) == i,]
  df5 <- rbind(df5, df4)

}

col_mat = viridis(length(cells),alpha = 0.8,option = "H")
ord <- list(NULL, with(df5, order(df5$cell,df5$to), NULL))

pdf(paste(resultsdir,"/",'top20top20sankeyblob.pdf',sep=""),width=8,height=10,paper='special') 
alluvial(df5[,1:2], freq = df5$value,col=df5$col,border = NA,blocks = TRUE,ordering=ord)#,col=ifelse( df3$from == "PAX6" & df3$to == df3$to, "red", "gray") )#, freq=df3$from, border=NA)
dev.off()




pdf(paste(resultsdir,"/",'chordssingleESC.pdf',sep=""),width=6,height=6,paper='special') 
for (i in unique(files$cell_id)){
  print(i)
  diffnetwork <- read.table(files$diffnetwork[files$cell_id==i])
  factors <- read.table(files$network[files$cell_id==i],header = T)
  factors <- factors[order(factors$sumScaled,decreasing = T),]
  
  # selecting the top 20 factors
  top <- factors$factor[1:20]
  
  colnames(diffnetwork) <- c("from","to","value")
  df <- data.frame(diffnetwork)
  
  # selecting
  df3 <- df[df$from %in% top,]
  df3 <- df3[order(df3$value,decreasing = T),]

  for (z in unique(top)){
    print(z)
    df4 <- df3
    df5 <- df4[df4$from == z,]
    df6 <- df5[df5$to == head(df5$to,10),]
    df7 = NULL
    
    # Determine if one of the top 10 factors also influences the factor of interest
    for (q in unique(head(df5$to,10))){
      print(q)
      df7 = rbind(df7, df4[df4$from == q & df4$to == z,])
      df7 = rbind(df7, df4[df4$from == q & df4$to %in% head(df5$to,10),])
      
      # If a factor is found influencing the factor of interest or another factor in the top 10 then add to original dataframe
      if (!is.null(df7)){
        df6 = rbind(df6, df7)#[df4$from == q & df4$to %in% head(df5$to,10) | df4$to == z,])
      }
    }
    
    # making the dataframe consist of only unique rows
    df6 <- unique(df6)
    
    # Setting the right color scale
    u1 <- unique(unlist(df6[, c("from", "to")]))
    col_mat = viridis(length(u1),alpha = 0.8,option = "H")
  
    chordDiagram(df6, directional = 1, direction.type = c("diffHeight", "arrows"),
               link.arr.type = "big.arrow", grid.col = col_mat)
    title(paste(i , "top 10 regulated genes by", z))
  }
}

dev.off()
