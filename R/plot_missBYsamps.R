#' Plot missing data, by samples
#'
data("genomalicious_4pops")
dat <- genomalicious_4pops
missIdx <- c(sample(1:nrow(dat), size=0.1*nrow(dat), replace=FALSE)
             , 100:500, 800:1200, 8000:8888, 200000:200500)
dat$GT[missIdx] <- NA

type <- 'hist'
look <- 'ggplot'
sampCol <- 'SAMPLE'
respCol <- 'GT'
locusCol <- 'LOCUS'

popCol <- 'POP'
plotColour <- 'white'




stats <- FALSE

############## plot_missingBYsamples <- function(){}
colnames(dat)[
  match(c(locusCol, respCol, sampCol), colnames(dat))
  ] <- c('LOCUS', 'RESP', 'SAMPLE')

if(is.na(popCol)==FALSE){
  colnames(dat)[which(colnames(dat)==popCol)] <- 'POP'
}

# Set the plot theme by look
if(look=='ggplot'){
  plotTheme <- theme_gray() + theme(legend.position='top'
                                    , text=element_text(size=12))
} else if(look=='classic'){
  plotTheme <- theme_bw() + theme(panel.grid.major=element_blank()
                                  , panel.grid.minor=element_blank()
                                  , text=element_text(colour='black', size=12)
                                  , legend.position='top')
}

if(type=='hist'){
  if(is.na(popCol)){
    stats <- dat[, sum(is.na(RESP)), by=SAMPLE]
    gg <- (ggplot(stats, aes(x=V1))
         + plotTheme
         + geom_histogram(fill=plotColour[1], colour='black')
         + labs(x='Missing data', y='Number of samples'))
  } else{
    stats <- dat[, sum(is.na(RESP)), by=c('SAMPLE', 'POP')]
    gg <- (ggplot(stats, aes(x=V1))
         + plotTheme + theme(strip.text.x=element_text(face='bold'))
         + geom_histogram(fill=plotColour[1], colour='black')
         + labs(x='Missing data', y='Number of samples')
         + facet_wrap(~ POP))
  }

  # Heat map
  # https://learnr.wordpress.com/2010/01/26/ggplot2-quick-heatmap-plotting/
  if(type=='heatmap'){
    if(length(plotColour)==1){
      plotColour <- c('white', 'royalblue')
    }

    dat[, MISS:=as.integer(!is.na(RESP))]
    gg <- (ggplot(dat, aes(x=SAMPLE, y=LOCUS))
      + theme(
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank()
        )
      + geom_tile(aes(fill=as.factor(MISS)), colour=NA)
      + scale_fill_manual(values=c('0'=plotColour[1], '1'=plotColour[2]))
      + facet_wrap(~ POP)
      )
    if(is.na(popCol)==FALSE){
      gg <- gg + facet_wrap(~ POP)
    }
  }

}
