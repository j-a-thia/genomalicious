#' Plot missing data, by samples
#'
#' Use to visualise missing data with respect to samples and
#' their associated populations.
#'
#' @param dat Data table: Contains genetic information and must have
#' the following columns,
#' \enumerate{
#'   \item (1) The sampled individuals (see param \code{sampCOl}).
#'   \item (2) The locus ID (see param \code{locusCol}).
#'   \item (3) The response column, e.g. a genotype of allele frequency,
#'   (see param \code{respCol}).
#' }
#'
#' @param type Character: The type of plot to make. A histogram depicting
#' the frequency distribution of samples missing information at each
#' locus (\code{'hist'}, the default)? Or a heatmap illustrating the
#' missing loci for each sample (\code{'heatmap'})?
#'
#' @param look Character: The look of the plot. Default = \code{'ggplot'}, the
#' typical gray background with gridlines produced by \code{ggplot2}. Alternatively,
#' when set to \code{'classic'}, produces a base R style plot.
#'
#' @param sampCol Character: The column name with the sampled
#' individual information. Default = \code{'SAMPLE'}.
#'
#' @param locusCol Character: The column name with the locus information.
#' Default = \code{'LOCUS'}.
#'
#' @param respCol Character: The column name with the response information.
#' Default = \code{'GT'}. Missing data should be represeted by \code{NA}.
#'
#' @param popCol Character: The column name with the population information.
#' Optional parameter. Default = \code{NA}.
#'
#' @param plotColours Character: Vector of colours to use in plotting.
#' Size depends on values specified for parameters \code{type} and
#' \code{popCol}, see Details.
#'
#' @param plotNCol Integer: The number of columns to arrange indiviudal
#' population plots into. Only takes effect when \code{popCol} is specified.
#' Default = 2.
#'
#' @details 1) Return a plot, 2) display plot 3) Break down the different
#' plot types and population inclusions
#'
#' @examples
#' data(genomalicious_4pops)
#' dat <- genomalicious_4pops
#'
#' # Add missing values.
#' missIdx <- c(sample(1:nrow(dat), size=0.05*nrow(dat), replace=FALSE)
#'              , 100:500, 800:1200, 8000:8888, 200000:200500)
#' dat$GT[missIdx] <- NA
#'
#' @export
plot_missBYsamps <- function(dat, type='hist', look='ggplot'
                             , sampCol='SAMPLE', locusCol='LOCUS'
                             , respCol='GT', popCol=NA
                             , plotColours='white', plotNCol=2){
  # --------------------------------------------+
  # Libraries, assertions, and setup
  # --------------------------------------------+
  for(lib in c('ggplot2', 'data.table','gridExtra')){ require(lib, character.only = TRUE)}

  # Rename columns
  colnames(dat)[
    match(c(locusCol, respCol, sampCol), colnames(dat))
    ] <- c('LOCUS', 'RESP', 'SAMPLE')

  # Rename the population column, if it was specified
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

  # --------------------------------------------+
  # Code
  # --------------------------------------------+

  # Histogram
  if(type=='hist'){
    # If
    if(is.na(popCol)){
      # If no population column is specified (NA)
      stats <- dat[, sum(is.na(RESP)), by=SAMPLE]
      gg <- (ggplot(stats, aes(x=V1))
             + plotTheme
             + geom_histogram(fill=plotColours[1], colour='black')
             + labs(x='Missing data', y='Number of samples'))
    } else{
      # If population column is specified, facet plot by population
      stats <- dat[, sum(is.na(RESP)), by=c('SAMPLE', 'POP')]
      gg <- (ggplot(stats, aes(x=V1))
             + plotTheme + theme(strip.text.x=element_text(face='bold'))
             + geom_histogram(fill=plotColours[1], colour='black')
             + labs(x='Missing data', y='Number of samples')
             + facet_wrap(~ POP, ncol=plotNCol))
    }
  }
  # Heat map
  if(type=='heatmap'){
    # Make two colours, if only ine specified
    if(length(plotColours)<=1){
      plotColours <- c('white', 'royalblue')
    }

    # Assign a new column to record missing data
    dat[, MISS:=as.integer(!is.na(RESP))]

    if(is.na(popCol)){
      # If no population column is specified (NA)
      gg <- (ggplot(dat, aes(x=SAMPLE, y=LOCUS))
             + theme(
               axis.text.x = element_blank()
               ,axis.text.y = element_blank()
               ,axis.ticks = element_blank()
               ,panel.border=element_rect(fill=NA, colour='black')
               ,legend.position='none'
             )
             + geom_tile(aes(fill=as.factor(MISS)), colour=NA)
             + scale_fill_manual(values=c('0'=plotColours[1], '1'=plotColours[2]))
             + labs(x='Samples', y='Locus')
      )
    } else {
      # If population column is specified, make indiviudal plot
      # for each population
      ggLs <- lapply(unique(dat$POP), function(pop){
        g <- (ggplot(dat[POP==pop], aes(x=SAMPLE, y=LOCUS))
              + theme(
                axis.text.x = element_blank()
                ,axis.text.y = element_blank()
                ,axis.ticks = element_blank()
                ,panel.border=element_rect(fill=NA, colour='black')
                ,plot.title=element_text(hjust=0.5, size=12, face='bold')
                ,legend.position='none'
              )
              + geom_tile(aes(fill=as.factor(MISS)), colour=NA)
              + scale_fill_manual(values=c('0'=plotColours[1], '1'=plotColours[2]))
              + labs(x='Samples', y='Locus', title=pop)
        )
        return(g)
      })
      gg <- do.call('grid.arrange', c(ggLs, ncol=plotNCol))
    }
  }

  # FInish up
  plot(gg)
  return(gg)
}
