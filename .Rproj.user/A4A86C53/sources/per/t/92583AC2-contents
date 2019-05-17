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
#' @param type Character: The type of plot to make. A histogram depicting
#' the frequency distribution of samples missing information at each
#' locus (\code{'hist'}, the default)? Or a heatmap illustrating the
#' missing loci for each sample (\code{'heatmap'})? See Details.
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
#' @details Whenever \code{popCol} is unspecififed (Default = \code{NA}),
#' then any plot is created from all samples. If \code{popCol} is specified,
#' i.e. to include a column containing population IDs, then individual
#' plots are created for each population. These are arranged in rows by
#' columns, with the number of columns specified in \code{plotNCol}.
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
#' data(genomalicious_4pops)
#' datGt <- genomalicious_4pops
#'
#' # Add missing values
#' missIdx <- c(sample(1:nrow(datGt), size=0.05*nrow(datGt), replace=FALSE)
#'              , 100:500, 800:1200, 8000:8888, 200000:200500)
#'
#' datGt$GT[missIdx] <- NA
#'
#' head(datGt, 10)
#'
#' # Histograms, ggplot and classic looks
#' plot_missBYsamps(datGt, type='hist', look='ggplot')
#' plot_missBYsamps(datGt, type='hist', look='classic')
#'
#' # Histograms, by population, specifying colour
#' plot_missBYsamps(datGt, type='hist', look='ggplot'
#'                  , popCol='POP' , plotColours='plum2')
#'
#' # Heatmaps, default and specified colours
#' plot_missBYsamps(datGt, type='heatmap')
#' plot_missBYsamps(datGt, type='heatmap', plotColours=c('black', 'plum2'))
#'
#' # Heatmaps, by population
#' plot_missBYsamps(datGt, type='heatmap', popCol='POP')
#'
#' ####   CATCH PLOT OUTPUT FOR LATER USE   ####
#' gg4pops <- plot_missBYsamps(datGt, popCol='POP')
#' plot(gg4pops)
#'
#' @export
plot_missBYsamps <- function(dat, type='hist', look='ggplot'
                             , sampCol='SAMPLE', locusCol='LOCUS'
                             , genoCol='GT', popCol=NA
                             , plotColours='white', plotNCol=2){

  # --------------------------------------------+
  # Libraries, assertions, and setup
  # --------------------------------------------+
  for(lib in c('ggplot2', 'data.table','gridExtra')){ require(lib, character.only = TRUE)}

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

  # Histogram
  if(type=='hist'){
    # If
    if(is.na(popCol)){
      # If no population column is specified (NA)
      stats <- dat[, sum(is.na(GT)), by=SAMPLE]
      gg <- (ggplot(stats, aes(x=V1))
             + plotTheme
             + geom_histogram(fill=plotColours[1], colour='black')
             + labs(x='Missing genotypes', y='Number of samples'))
    } else{
      # If population column is specified, facet plot by population
      stats <- dat[, sum(is.na(GT)), by=c('SAMPLE', 'POP')]
      gg <- (ggplot(stats, aes(x=V1))
             + plotTheme + theme(strip.text.x=element_text(face='bold'))
             + geom_histogram(fill=plotColours[1], colour='black')
             + labs(x='Missing genotypes', y='Number of samples')
             + facet_wrap(~ POP, ncol=plotNCol))
    }
  }
  # Heat map
  if(type=='heatmap'){
    # Make two colours, if only ine specified
    if(length(plotColours)<2){
      plotColours <- c('white', 'royalblue')
    }

    # Assign a new column to record missing data
    dat[, MISS:=as.integer(!is.na(GT))]

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
