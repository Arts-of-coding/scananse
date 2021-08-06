# Deseq2 significance calculation peaks
##based upon http://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
workdir <- "/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/blob"

## setting up results directory
dateoftoday <- gsub("-", "", as.character(Sys.Date()))
resultsdir <- paste0(workdir, dateoftoday)
system(paste("mkdir -p ", resultsdir))
# loading in the count matrix and coldata files

lakocountfile <- read.table(file = '/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq/20210710/joinedcounts2.tsv', sep = '\t', header = TRUE, row.names = 1)

lakocountfile <- lakocountfile[,colnames(lakocountfile) != "ESC1" & colnames(lakocountfile) != "ESC2" & colnames(lakocountfile) != "IC1" & colnames(lakocountfile) != "IC2" & colnames(lakocountfile) != "FCECs1" & colnames(lakocountfile) != "FCECs2" & colnames(lakocountfile) != "Ves1" & colnames(lakocountfile) != "Ves2"]
coldata <- read.table(file = '/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq/20210710/col2.tsv', sep = '\t', header = TRUE, row.names = 1)
coldata <- coldata[rownames(coldata) != "ESC1" & rownames(coldata) != "ESC2" & rownames(coldata) != "IC1" & rownames(coldata) != "IC2" & rownames(coldata) != "FCECs1" & rownames(coldata) != "FCECs2" & rownames(coldata) != "Ves1" & rownames(coldata) != "Ves2",]

coldata <- coldata[,c("condition","type")]

cellcounts <- read.table(file = '/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq/20210730/cellcounts.tsv', sep = '\t', header = TRUE, row.names = 1)

# setting the correct columns for the coldata
coldata$condition <- factor(coldata$condition)
coldata$type <- factor(coldata$type)

# Normalizing the counts for the minimal amount of cells
coldata$count <- cellcounts[match(coldata$condition,rownames(cellcounts)),]
coldata$ratio <- min(coldata$count)/coldata$count

library("DESeq2")

for (i in colnames(lakocountfile)){
  print(i)
  lakocountfile[,i] <- lakocountfile[,i] * coldata[rownames(coldata) == i,]$ratio
}

# Generating the matrix
lakoinf <- as.matrix(lakocountfile[,colnames(lakocountfile) != "ESC1" & colnames(lakocountfile) != "ESC2" & colnames(lakocountfile) != "IC1" & colnames(lakocountfile) != "IC2" & colnames(lakocountfile) != "FCECs1" & colnames(lakocountfile) != "FCECs2" & colnames(lakocountfile) != "Ves1" & colnames(lakocountfile) != "Ves2"],row.names="ID")
lakoinf <- as.matrix(lakocountfile,row.names="ID")

# checking the data
head(lakoinf,2)

coldata
coldata$condition

for (i in coldata$condition){
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
  lakoinfcalc3 <- as.data.frame(lakoinfcalc@listData$log2FoldChange,row.names = rownames(lakoinf2))
  lakoinfcalc3

  print(i)
  lakoinfcalc <- results(dds2,contrast = c("condition", i,"blob"))
  lakoinfcalc3 <- as.data.frame(lakoinfcalc@listData$log2FoldChange,row.names = rownames(lakoinf2))
  lakoinfcalc3$padj <- lakoinfcalc@listData$padj
  colnames(lakoinfcalc3) <- c("log2FoldChange","padj")
  lakoinfcalc3 <- lakoinfcalc3[!is.na(lakoinfcalc3$log2FoldChange) & !is.na(lakoinfcalc3$padj) ,]
  write.table(data.frame("resid"=rownames(lakoinfcalc3),lakoinfcalc3), file = paste0(resultsdir, "/", i, 'blobpseudobulkpadj.tsv'),sep = "\t", quote = F, row.names = F)
}