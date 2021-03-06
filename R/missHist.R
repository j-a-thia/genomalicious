#' Plot missing genotypes, by samples
#'
#' Use to visualise missing data with respect to samples and
#' their associated populations.
#'
#' @param dat Data table: Contains genetic information and must have
#' the following columns,
#' \enumerate{
#'   \item The sampled individuals (see param \code{sampCol}).
#'   \item The locus ID (see param \code{locusCol}).
#'   \item The genotype column, e.g. a genotype of allele frequency,
#'   (see param \code{genoCol}).
#' }
#'
#' @param look Character: The look of the plot. Default = \code{'ggplot'}, the
#' typical gray background with gridlines produced by \code{ggplot2}. Alternatively,
#' when set to \code{'classic'}, produces a base R style plot.
#'
#' @param sampCol Character: The column name with the sampled
#' individual ID. Default = \code{'SAMPLE'}.
#'
#' @param locusCol Character: The column name with the locus ID.
#' Default = \code{'LOCUS'}.
#'
#' @param genoCol Character: The column name with the genotype info.
#' Default = \code{'GT'}. Missing data should be represeted by \code{NA}.
#'
#' @param popCol Character: The column name with the population ID.
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
#' @details When \code{popCol} is unspecified, then all samples are used to create the plots.
#' If it is specified, then that column name is used to make one plot for
#' each population. These are arranged in rows and columns, and the
#' user can specify the number of columns with the argument \code{plotNCol}.
#'
#' The number of colours that need to be specified in \code{plotColours}
#' depends on which kind of plot is being created (\code{type}) and
#' whether the samples are to be grouped by population.
#' \itemize{
#'    \item \code{type=='hist'}: This is a histogram with the number
#'       of samples with missing data at a locus (x-axis) vs the
#'       frequency of such observations (y-axis). Only a single colour
#'       needs to be specified, \code{plotColours[1]}, which is the
#'       fill colour of the histogram bars.
#'
#'    \item \code{type=='heatmap'}: This is a heatmap of missing data
#'       with samples (x-axis) vs loci (y-axis), with missing and present
#'       values highlighted in two different colours. The first colour,
#'       \code{plotColours[1]}, should be the colour for missing data.
#'       The second colour, \code{plotColours[2]}, therefore is the colour
#'       for present data. If there are <2 colours, will default to
#'       'white' and 'royalblue'.
#'
#' @return Displays plot and returns the plot object.
#'
#' @examples
#' ####   MISSING GENOTYPE DATA   ####
#' data(data_4pops)
#' datGt <- data_4pops
#'
#' # Add missing values
#' datGt <- do.call('rbind'
#'                 , lapply(split(datGt, datGt$SAMPLE), function(x){
#'                   if(x$POP[1]=='Pop1'){ pr <- 0.1
#'                   } else if(x$POP[1]=='Pop2'){ pr <- 0.2
#'                   } else{ pr <- 0.05}
#'                   numMiss <- rnbinom(1, size=8, prob=pr)
#'                   idxMiss <- sample(1:nrow(x), size=numMiss, replace=FALSE)
#'                   x$GT[idxMiss] <- NA
#'                  return(x)
#'                 }))
#'
#' head(datGt, 10)
#'
#' ####   PLOT MISSING BY SAMPLES   ####
#' # Histograms, ggplot and classic looks
#' missHist(datGt, plotBy='samples', look='ggplot')
#' missHist(datGt, plotBy='samples',, look='classic')
#'
#' # Histograms, by population, specifying colour
#' missHist(datGt, plotBy='samples',, look='ggplot'
#'                  , popCol='POP' , plotColours='plum2')
#'
#' ####   PLOT MISSING BY LOCI   ####
#' missHist(datGt, plotBy='loci',, look='classic'
#'                  , popCol='POP' , plotColours='plum2')
#'
#' ####   CATCH PLOT OUTPUT FOR LATER USE   ####
#' gg4pops <- missHist(datGt, plotBy='samples',, popCol='POP')
#' plot(gg4pops)
#'
#' @export
missHist <- function(dat, plotBy, look='ggplot'
                             , sampCol='SAMPLE', locusCol='LOCUS'
                             , genoCol='GT', popCol=NA
                             , plotColours='white', plotNCol=2){

  # --------------------------------------------+
  # Libraries, assertions, and setup
  # --------------------------------------------+
  for(lib in c('ggplot2', 'data.table','gridExtra')){ require(lib, character.only = TRUE)}

  # Check that plotBy is specified properly
  if(sum(c('samples', 'loci') %in% plotBy)!=1){
    stop("Argument plotBy must be of length==1, and be either 'samples' or 'loci'.")
  }

  # A column specifier to focus data
  if(plotBy=='samples'){ focus <- 'SAMPLE'
  } else{ focus <- 'LOCUS' }

  # Rename columns
  colnames(dat)[
    match(c(locusCol, genoCol, sampCol), colnames(dat))
    ] <- c('LOCUS', 'GT', 'SAMPLE')

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

  if(is.na(popCol)){
    # If no population column is specified (NA)
    stats <- dat[, sum(is.na(GT)), by=focus]
    gg <- (ggplot(stats, aes(x=V1))
           + plotTheme
           + geom_histogram(fill=plotColours[1], colour='black')
           + labs(x='Missing genotypes'
                  , y=paste0('Number of ', plotBy)
                  ))
  } else{
    # If population column is specified, facet plot by population
    stats <- dat[, sum(is.na(GT)), by=c(focus, 'POP')]
    gg <- (ggplot(stats, aes(x=V1))
           + plotTheme + theme(strip.text.x=element_text(face='bold'))
           + geom_histogram(fill=plotColours[1], colour='black')
           + labs(x='Missing genotypes'
                  , y=paste0('Number of ', plotBy))
           + facet_wrap(~ POP, ncol=plotNCol))
  }

  # Finish up
  plot(gg)
  return(gg)
}
