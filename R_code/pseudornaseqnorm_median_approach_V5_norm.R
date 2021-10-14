###########################################################################
# Deseq2 significance calculation DE genes (of the median approach)
# Based upon http://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html

# Setting up the working directory
workdir <- "/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/blob"

# Setting up results directory
dateoftoday <- gsub("-", "", as.character(Sys.Date()))
resultsdir <- paste0(workdir, dateoftoday)
system(paste("mkdir -p ", resultsdir))

# Loading in important packages
library("DESeq2")

###########################################################################
# Loading in the count matrix and coldata files
# Excluding cell populations that were not used in ANANSE

lakocountfile <- read.table(file = '/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq/20210710/joinedcounts2.tsv', sep = '\t', header = TRUE, row.names = 1)
lakocountfile <- lakocountfile[,colnames(lakocountfile) != "ESC1" & colnames(lakocountfile) != "ESC2" & colnames(lakocountfile) != "IC1" & colnames(lakocountfile) != "IC2" & colnames(lakocountfile) != "FCECs1" & colnames(lakocountfile) != "FCECs2" & colnames(lakocountfile) != "Ves1" & colnames(lakocountfile) != "Ves2"]
coldata <- read.table(file = '/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq/20210710/col2.tsv', sep = '\t', header = TRUE, row.names = 1)
coldata <- coldata[rownames(coldata) != "ESC1" & rownames(coldata) != "ESC2" & rownames(coldata) != "IC1" & rownames(coldata) != "IC2" & rownames(coldata) != "FCECs1" & rownames(coldata) != "FCECs2" & rownames(coldata) != "Ves1" & rownames(coldata) != "Ves2",]
coldata <- coldata[,c("condition","type")]

# Setting the correct columns for the coldata
coldata$condition <- factor(coldata$condition)
coldata$type <- factor(coldata$type)

###########################################################################
# Median cell count normalization

# Loading in the cell count file for median normalization of large populations
cellcounts <- read.table(file = '/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq/20210730/cellcounts.tsv', sep = '\t', header = TRUE, row.names = 1)

# Normalizing the counts for the median amount of cells
 coldata$count <- cellcounts[match(coldata$condition,rownames(cellcounts)),]
 coldata$ratio <- median(coldata$count)/coldata$count
 coldata$ratio[coldata$ratio>1] <-1
 for (i in colnames(lakocountfile)){
  print(i)
  lakocountfile[,i] <- lakocountfile[,i] * coldata[rownames(coldata) == i,]$ratio
 }

lakoinf <- as.matrix(lakocountfile,row.names="ID")

# Checking the data
#head(lakoinf,2)
#coldata
#coldata$condition

###########################################################################
# Calculating DE genes
for (i in unique(coldata$condition)){
lakoinf2 <- lakoinf[, rownames(coldata)]
all(rownames(coldata) == colnames(lakoinf))

coldata2 <- coldata
coldata2$condition <- as.character(coldata2$condition)   
coldata2$condition[coldata2$condition != i] <- "blob"

print(coldata2$condition)

dds <- DESeqDataSetFromMatrix(countData = round(lakoinf2),
                              colData = coldata2,
                              design = ~ condition)
dds

dds2 <- DESeq(dds)
dds2

lakoinfcalc <- results(dds2)
lakoinfcalc3 <- as.data.frame(lakoinfcalc@listData$log2FoldChange,row.names = rownames(lakoinf))
lakoinfcalc3

print(i)
lakoinfcalc <- results(dds2,contrast = c("condition", i,"blob"))
lakoinfcalc3 <- as.data.frame(lakoinfcalc@listData$log2FoldChange,row.names = rownames(lakoinf))
lakoinfcalc3$padj <- lakoinfcalc@listData$padj
colnames(lakoinfcalc3) <- c("log2FoldChange","padj")
lakoinfcalc3 <- lakoinfcalc3[!is.na(lakoinfcalc3$log2FoldChange) & !is.na(lakoinfcalc3$padj) ,]

# Writing out the values of DE expression and the adjoined Z-score

write.table(data.frame("resid"=rownames(lakoinfcalc3),lakoinfcalc3), file = paste0(resultsdir, "/", i, 'blobpseudobulkpadj.tsv'),sep = "\t", quote = F, row.names = F, col.names = T)
}
