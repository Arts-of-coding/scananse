# Deseq2 matrix vst normalisation
##based upon http://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html

# loading in the count matrix and coldata files

# test matrix file
lakocountfile <- read.table(file = '/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/data/lako2021/merge_peaks/2021-27-06/tmp/newtable.tsv', sep = '\t', header = TRUE, row.names = 1)

#lakocountfile <- read.table(file = '/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/20210701scRNA_integration/col.tsv', sep = '\t', header = TRUE, row.names = 1)
lakovst <- as.matrix(lakocountfile,row.names="loc")

#test coldata
coldata <- read.table(file = '/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/data/lako2021/merge_peaks/2021-27-06/tmp/col.tsv', sep = '\t', header = TRUE, row.names = 1)

#coldata <- read.table(file = '/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/20210701scRNA_integration/col.tsv', sep = '\t', header = TRUE, row.names = 1)

coldata <- coldata[,c("condition","type")]

# setting the correct columns for the coldata
coldata$condition <- factor(coldata$condition)
coldata$type <- factor(coldata$type)

# checking the data
head(lakovst,2)

coldata

# setting the rowdata in coldata similar to the coldata in the count matrix
rownames(coldata) <- sub("fb", "", rownames(coldata))
all(rownames(coldata) %in% colnames(lakovst))

all(rownames(coldata) == colnames(lakovst))

lakovst <- lakovst[, rownames(coldata)]
all(rownames(coldata) == colnames(lakovst))

library("DESeq2")
dds <- DESeqDataSetFromMatrix(countData = lakovst,
                              colData = coldata,
                              design = ~ condition)
dds

lakovstcalc <- vst(dds,blind = T)
lakovstcalc2 <- as.data.frame(assay(lakovstcalc), 3,row.names = rownames(lakovst))
lakovstcalc2

# calculating the mean of the samples to have one column for each population again
colnames(lakovstcalc2) <- coldata$condition

lakovstcalc3 <- sapply(split.default(lakovstcalc2, names(lakovstcalc2)), rowMeans)
lakovstcalc3 <- as.data.frame(lakovstcalc3,row.names = rownames(lakovst))
lakovstcalc3

write.csv(lakovstcalc3, 
          file="/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/data/lako2021/merge_peaks/2021-27-06/tmp/vstnosplit.csv")
