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
#' @param plotBy Character: One of 'samples' or 'loci', the focus of missing data.
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
#' @param plotColours Character: The fill colour for histogram bars.
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
#' @return Returns a ggplot object.
#'
#' @examples
#' library(genomalicious)
#'
#' ####   MISSING GENOTYPE DATA   ####
#' data(data_Genos)
#' datGt <- data_Genos
#'
#' # Add missing values
#' datGt <- do.call(
#'  'rbind',
#'  # Split data table by sample, and iterate through samples, X
#'  split(datGt, by='POP') %>%
#'    lapply(., function(Dpop){
#'      pop <- Dpop$POP[1]
#'
#'      if(pop=='Pop1'){
#'        pr <- 0.1
#'      } else if(pop=='Pop2'){
#'        pr <- 0.2
#'      } else if(pop %in% c('Pop3','Pop4')){
#'        pr <- 0.05
#'      }
#'
#'      # Numbers and unique loci and samples
#'      num.loc <- Dpop$LOCUS %>% unique %>% length
#'      uniq.loc <- Dpop$LOCUS %>% unique
#'      num.samp <- Dpop$SAMPLE %>% unique %>% length
#'      uniq.samp <- Dpop$SAMPLE %>% unique
#'
#'      # Vector of missingness
#'      num.miss <- rbinom(n=num.samp, size=num.loc, prob=pr)
#'
#'      # Iterate through samples and add unique loci
#'      for(i in 1:num.samp){
#'        locs <- sample(uniq.loc, size=num.miss[i], replace=FALSE)
#'        Dpop[SAMPLE==uniq.samp[i] & LOCUS%in%locs, GT:=NA]
#'      }
#'
#'      # Return
#'      return(Dpop)
#'    }
#'    )
#' )
#'
#' head(datGt, 10)
#'
#' ####   PLOT MISSING BY SAMPLES   ####
#' # Histograms, ggplot and classic looks
#' miss_plot_hist(datGt, plotBy='samples', look='ggplot')
#' miss_plot_hist(datGt, plotBy='samples',, look='classic')
#'
#' # Histograms, by population, specifying colour
#' miss_plot_hist(datGt, plotBy='samples',, look='ggplot'
#'                  , popCol='POP' , plotColours='deeppink2')
#'
#' ####   PLOT MISSING BY LOCI   ####
#' miss_plot_hist(datGt, plotBy='loci',, look='classic'
#'                  , popCol='POP' , plotColours='deeppink2')
#'
#' ####   CATCH PLOT OUTPUT FOR LATER USE   ####
#' gg4pops <- miss_plot_hist(datGt, plotBy='samples', popCol='POP')
#' plot(gg4pops)
#'
#' @export
miss_plot_hist <- function(
  dat, plotBy, look='ggplot', sampCol='SAMPLE', locusCol='LOCUS', genoCol='GT',
  popCol=NA, plotColours='white', plotNCol=2
  ){

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
    stats <- dat[, 100 * sum(is.na(GT))/length(GT), by=focus]
    gg <- (ggplot(stats, aes(x=V1))
           + plotTheme
           + geom_histogram(fill=plotColours[1], colour='black')
           + labs(x='Missing genotypes (%)'
                  , y=paste0('Number of ', plotBy)
                  ))
  } else{
    # If population column is specified, facet plot by population
    stats <- dat[, 100 * sum(is.na(GT))/length(GT), by=c(focus, 'POP')]
    gg <- (ggplot(stats, aes(x=V1))
           + plotTheme + theme(strip.text.x=element_text(face='bold'))
           + geom_histogram(fill=plotColours[1], colour='black')
           + labs(x='Missing genotypes (%)'
                  , y=paste0('Number of ', plotBy))
           + facet_wrap(~ POP, ncol=plotNCol))
  }

  # Finish up
  return(gg)
}
