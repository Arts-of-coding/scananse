#devtools::install_github("mattflor/chorddiag")
library(chorddiag)
library(circlize)
library(viridis)

## Storing results ##
workdir <- "/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/R/chord"

## setting up results directory
dateoftoday <- gsub("-", "", as.character(Sys.Date()))
resultsdir <- paste0(workdir, dateoftoday)
system(paste("mkdir -p ", resultsdir))

files <- read.csv("/ceph/rimlsfnwi/data/moldevbio/zhou/jarts/data/lako2021/ANANSE/chord_files.csv",header = T) 

pdf(paste(resultsdir,"/",'chordsESC.pdf',sep=""),width=6,height=6,paper='special') 
for (i in unique(files$cell_id)){
  print(i)
  diffnetwork <- read.table(files$diffnetwork[files$cell_id==i])
  factors <- read.table(files$network[files$cell_id==i],header = T)
  factors <- factors[order(factors$sumScaled,decreasing = T),]
  # selecting the top 10 factors
  top <- factors$factor[1:10]
  
  colnames(diffnetwork) <- c("from","to","value")
  df <- data.frame(diffnetwork)
  
  # cross referencing the top 15 factors against each other
  df3 <- df[df$from %in% top,]
  df3 <- df3[df3$to %in% top,]
  df3 <- df3[order(df3$value,decreasing = T),]
  
  # cutoff value of influence score >0.95
  df4 <- df3[df3$value>0.95,]
  df4 <- df4[order(df4$from),]
  u1 <- unique(unlist(df4[, c("from", "to")]))
  
  col_mat = viridis(length(u1),alpha = 0.8,option = "H")
  
  
  chordDiagram(df4, directional = 1, direction.type = c("diffHeight", "arrows"),
               link.arr.type = "big.arrow", grid.col = col_mat)
  title(paste(i ,"top 15 tfs >0.95 influence network vs ESCs"))
}

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
