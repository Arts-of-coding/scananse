# Deseq2 significance calculation peaks
##based upon http://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html

# loading in the count matrix and coldata files

# test matrix file
#lakocountfile <- read.table(file = '/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq/20210710/joinedcounts.tsv', sep = '\t', header = TRUE, row.names = 1)
lakocountfile <- read.table(file = '/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/data/lako2021/merge_peaks/2021-08-25_4cellreps/tmp/joinedcovtable.tsv', sep = '\t', header = TRUE, row.names = 1)
lakopeaks <- as.matrix(lakocountfile,row.names="loc")

#coldata <- read.table(file = '/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/20210701scRNA_integration/col.tsv', sep = '\t', header = TRUE, row.names = 1)
coldata <- read.table(file = '/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/data/lako2021/merge_peaks/2021-08-25_4cellreps/tmp/col2.tsv', sep = '\t', header = TRUE, row.names = 1)
coldata <- coldata[,c("condition","type")]

# setting the correct columns for the coldata
coldata$condition <- factor(coldata$condition)
coldata$type <- factor(coldata$type)

# checking the data
head(lakopeaks,2)

coldata

# setting the rowdata in coldata similar to the coldata in the count matrix
rownames(coldata) <- sub("fb", "", rownames(coldata))
all(rownames(coldata) %in% colnames(lakopeaks))

all(rownames(coldata) == colnames(lakopeaks))

lakopeaks <- lakopeaks[, rownames(coldata)]
all(rownames(coldata) == colnames(lakopeaks))

library("DESeq2")
dds <- DESeqDataSetFromMatrix(countData = lakopeaks,
                              colData = coldata,
                              design = ~ condition)
dds

dds2 <- DESeq(dds)
dds2
lakopeakscalc <- results(dds2)

lakopeakscalc3 <- as.data.frame(lakopeakscalc@listData$padj,row.names = rownames(lakopeaks))
lakopeakscalc3
names(lakopeakscalc3) <- "padj"

lakopeakscalc3 <- na.omit(lakopeakscalc3)
rownames(lakopeakscalc3)
lakopeakscalc3 <- row.names(lakopeakscalc3)[apply(lakopeakscalc3, 1, function(u) any(u < 0.05))]
# calculating the mean of the samples to have one column for each population again

lakopeakscalc3 <- as.data.frame(lakopeakscalc3,row.names = rownames(lakopeaks))
lakopeakscalc3

write.table(lakopeakscalc3, file = "/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/data/lako2021/merge_peaks/2021-08-25_4cellreps/tmp/sigregions.txt", append = FALSE, quote = F, sep = " ",
            eol = "\n", na = "NA", dec = ".", row.names = F,
            col.names = F, qmethod = c("escape", "double"),
            fileEncoding = "")