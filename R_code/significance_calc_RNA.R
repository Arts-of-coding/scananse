# Deseq2 significance calculation peaks
##based upon http://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html

# loading in the count matrix and coldata files

# For the deseq2 matrix vs ESC
#lakocountfile <- read.table(file = '/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq/20210710/joinedcounts2.tsv', sep = '\t', header = TRUE, row.names = 1)
lakocountfile <- read.table(file = '/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq/20210802/pseudobulk_4cells.tsv', sep = '\t', header = TRUE, row.names = 1)
workdir <- "/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/Z-score"

## setting up results directory
dateoftoday <- gsub("-", "", as.character(Sys.Date()))
resultsdir <- paste0(workdir, dateoftoday)
system(paste("mkdir -p ", resultsdir))

collist<- c("LiCo","StCSC")

# for complex heatmap
#lakocountfile <- read.table(file = '/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq/20210710/pseudobulk.tsv', sep = '\t', header = TRUE, row.names = 1)
lakovst <- as.matrix(lakocountfile,row.names="ID")

#test coldata
#coldata <- read.table(file = '/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq/20210710/col2.tsv', sep = '\t', header = TRUE, row.names = 1)
coldata <- read.table(file = '/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/scRNA-seq/20210802/col.tsv', sep = '\t', header = TRUE, row.names = 1)

# for complex heatmap
#coldata <- coldata[1:22,]

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
dds <- DESeqDataSetFromMatrix(countData = round(lakovst),
                              colData = coldata,
                              design = ~ condition)
dds

###################################################################
# for complex heatmap
dds2 <- DESeq(dds)

# transpose and scale the matrix per row Z-score per gene per sample is scaled upon 
vsd <- assay(vst(dds2,blind = T))
Z <- t(scale(t(vsd)))
Z
Z_score <- as.data.frame(Z,row.names = rownames(lakovst))

# joining the columns on the means sequential
n <- 2
Z-joined <- t(rowMeans(t(Z-score), as.integer(gl(ncol(Z-score), n, ncol(Z-score))))) / n

df %>% mutate(mean_all = rowMeans(.),
              mean_sel = rowMeans(select(., select_vars)))

#  generate list of factors
vec <- rownames(coldata)
x <- split(vec, ceiling(seq_along(vec)/2))

scoretable <- as.data.frame(do.call(cbind, lapply(x, function(i) rowMeans(Z_score[, i]))), row.names = rownames(lakovst))
colnames(scoretable) <- unique(coldata$condition)
scoretable <- scoretable[,collist]
write.table(scoretable, file = paste0(resultsdir, '/Zscoretable.tsv'),sep = "\t", quote = F,row.names = T,col.names = T)
##################################################################
dds2 <- DESeq(dds)
dds2

lakovstcalc <- results(dds2)
lakovstcalc3 <- as.data.frame(lakovstcalc@listData$log2FoldChange,row.names = rownames(lakovst))

for (i in dds2$condition){
  if (i != "ESC"){
    print(i)
    lakovstcalc <- results(dds2,contrast = c("condition", i,"ESC"))
    lakovstcalc3 <- as.data.frame(lakovstcalc@listData$log2FoldChange,row.names = rownames(lakovst))
    lakovstcalc3$padj <- lakovstcalc@listData$padj
    colnames(lakovstcalc3) <- c("log2FoldChange","padj")
    write.table(lakovstcalc3, file = paste0(resultsdir, i, 'ESCpseudobulkpadj.tsv'),sep = "\t", quote = F)
    }
}
