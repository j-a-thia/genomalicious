#' Plot missing data, by samples
#'
dat <- genomalicious_4pops
missIdx <- sample(1:nrow(dat), size=0.15*nrow(dat), replace=FALSE)
dat$GT[missIdx] <- NA

sampCol <- 'SAMPLE'
respCol <- 'GT'
popCol <- NA
locusCol <- 'LOCUS'

type <- 'hist'

stats <- FALSE

############## plot_missingBYsamples <- function(){}
colnames(dat)[
  match(c(locusCol, respCol, sampCol), colnames(dat))
  ] <- c('LOCUS', 'RESP', 'SAMPLE')

if(is.na(popCol)==FALSE){
  colnames(dat)[which(colnames(dat)==popCol)] <- 'POP'
}

if(type=='hist'){
  if(is.na(popCol)){
    stats <- dat[, sum(is.na(RESP)), by=SAMPLE]
  } else{
    stats <- dat[, sum(is.na(RESP)), by=c('SAMPLE', 'POP')]
  }
  hist(stats$V1)
}
