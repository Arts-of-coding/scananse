# Deseq2 significance calculation peaks
##based upon http://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html

# loading in the count matrix and coldata files

workdir <- "/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/"

## setting up results directory
dateoftoday <- gsub("-", "", as.character(Sys.Date()))
resultsdir <- paste0(workdir, dateoftoday)
system(paste("mkdir -p ", resultsdir))


# For the deseq2 matrix vs ESC
lakocountfile <- read.table(file = '/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq/20210710/joinedcounts2.tsv', sep = '\t', header = TRUE, row.names = 1)
lakoinf <- as.matrix(lakocountfile,row.names="ID")
coldata <- read.table(file = '/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq/20210710/col2.tsv', sep = '\t', header = TRUE, row.names = 1)

coldata <- coldata[,c("condition","type")]

# setting the correct columns for the coldata
coldata$condition <- factor(coldata$condition)
coldata$type <- factor(coldata$type)

# checking the data
head(lakoinf,2)

coldata

# setting the rowdata in coldata similar to the coldata in the count matrix
rownames(coldata) <- sub("fb", "", rownames(coldata))
all(rownames(coldata) %in% colnames(lakoinf))

all(rownames(coldata) == colnames(lakoinf))

lakoinf <- lakoinf[, rownames(coldata)]
all(rownames(coldata) == colnames(lakoinf))

library("DESeq2")
dds <- DESeqDataSetFromMatrix(countData = round(lakoinf),
                              colData = coldata,
                              design = ~ condition)
dds

dds2 <- DESeq(dds)
dds2

lakoinfcalc <- results(dds2)
lakoinfcalc3 <- as.data.frame(lakoinfcalc@listData$log2FoldChange,row.names = rownames(lakoinf))

for (i in dds2$condition){
  if (i != "ESC"){
    print(i)
    lakoinfcalc <- results(dds2,contrast = c("condition", i,"ESC"))
    lakoinfcalc3 <- as.data.frame(lakoinfcalc@listData$log2FoldChange,row.names = rownames(lakoinf))
    lakoinfcalc3$padj <- lakoinfcalc@listData$padj
    colnames(lakoinfcalc3) <- c("log2FoldChange","padj")
    lakoinfcalc3 <- lakoinfcalc3[!is.na(lakoinfcalc3$log2FoldChange) & !is.na(lakoinfcalc3$padj) ,]
    write.table(data.frame("resid"=rownames(lakoinfcalc3),lakoinfcalc3), file = paste0(resultsdir, "/", i, 'ESCpseudobulkpadj.tsv'),sep = "\t", quote = F, row.names = F)
    }
}
