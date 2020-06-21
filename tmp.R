library(doParallel)
library(tidyverse)

numCores=5
numChunks=100

if(numCores==1){
  # Split string for each row, this will return a wide DT so
  # transpose back into long DT.
  formatDat <- t(vcfDT[, strsplit(DATA, ':')])

  # Give missing values NA
  formatDat[which(formatDat=='.')] <- NA

  # Label the columns
  colnames(formatDat) <- unlist(strsplit(vcfDT$FORMAT[1], ':'))

  # Combine with the main data table, drop the $FORMAT and $DATA columns
  vcfDT <- cbind(vcfDT[, !c('FORMAT','DATA'), with=FALSE], as.data.table(formatDat))

} else if(numCores>1){
  # Register the number of cores for parallel processing
  clust <- makeCluster(numCores)
  registerDoParallel(cl=clust)

  # The start position of the data chunks
  chunkSt <- seq(1, nrow(vcfDT), by=numChunks)

  # The max number of rows
  maxRows <- nrow(vcfDT)

  # Make row IDs
  vcfDT[, ROW.ID:=1:nrow(vcfDT)]

  # Iterate over the chunks
  formatDat <- foreach(i=chunkSt, .combine='rbind') %dopar% {
    require(data.table); require(tidyverse)
    i2j <- i:(i + numChunks - 1)
    i2j <- i2j[i2j<=maxRows]
    xx <- vcfDT[i2j,]
    yy <- t(xx[, strsplit(DATA, ':')])
    yy[which(yy=='.')] <- NA
    colnames(yy) <- unlist(strsplit(vcfDT$FORMAT[1], ':'))
    yy <- as.data.table(yy)
    yy$ROW.ID <- i2j
    yy
  }

  # End cluster
  stopCluster(clust)

  # Combine data table and the parsed format data
  vcfDT <- left_join(
    vcfDT[, !c('FORMAT','DATA'), with=FALSE]
    , formatDat, by='ROW.ID') %>%
    vcfDT[, !'ROW.ID']
}



detectCores()
