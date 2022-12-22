#' Plot a heatmap of missing information
#'
#' Use to visualise the presence or absence of genetic information
#' for all sample and locus combinations. Can also be subset
#' by population.
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
#' @param plotColours Character: Vector of colours to use in plotting with a
#' length of 2. The first colour is the missing colour, and the second colour
#' is the non-missing colour.
#'
#' @param plotNCol Integer: The number of columns to arrange indiviudal
#' population plots into. Only takes effect when \code{popCol} is specified.
#' Default = 2.
#'
#' @param showPlot Logical: Should the plot be shown automatically? Default is
#' \code{TRUE}, otherwise \code{FALSE}.
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
#' # Heatmaps, default and specified colours
#' missHeatmap(datGt)
#' missHeatmap(datGt, plotColours=c('black', 'deeppink2'))
#'
#' # Heatmaps, by population
#' missHeatmap(datGt, popCol='POP', plotNCol=4)
#'
#' ####   CATCH PLOT OUTPUT FOR LATER USE   ####
#' gg4pops <- missHeatmap(datGt, popCol='POP')
#' plot(gg4pops)
#'
#' @export
missHeatmap <- function(dat
                             , sampCol='SAMPLE', locusCol='LOCUS'
                             , genoCol='GT', popCol=NA
                             , plotColours='white', plotNCol=2, showPlot=TRUE){

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

  # --------------------------------------------+
  # Code
  # --------------------------------------------+

  # Make two colours, if only one specified
  if(length(plotColours)<2){
    plotColours <- c('white', 'royalblue')
  } else if(length(plotColours)>2){
    plotColours <- plotColours[1:2]
  }

  # Assign a new column to record missing data
  dat[, MISS:=as.integer(!is.na(GT))]

  if(is.na(popCol)){
    # If no population column is specified (NA)
    gg <- (ggplot(dat, aes(x=SAMPLE, y=LOCUS))
           + geom_tile(aes(fill=as.factor(MISS)), colour=NA)
           + scale_fill_manual(values=c('0'=plotColours[1], '1'=plotColours[2]))
           + labs(x='Samples', y='Locus')
           + theme(
             axis.text.x = element_blank()
             ,axis.text.y = element_blank()
             ,axis.ticks = element_blank()
             ,panel.border=element_rect(fill=NA, colour='black')
             ,panel.background=element_blank()
             ,panel.grid=element_blank()
             ,legend.position='none'
           )
    )
  } else {
    # If population column is specified, make indiviudal plot
    # for each population
    gg <- (ggplot(dat, aes(x=SAMPLE, y=LOCUS))
           + geom_tile(aes(fill=as.factor(MISS)), colour=NA)
           + scale_fill_manual(values=c('0'=plotColours[1], '1'=plotColours[2]))
           + labs(x='Samples', y='Locus')
           + facet_wrap(~POP, scales='free_x', ncol=plotNCol)
           + theme(
             axis.text.x = element_blank()
             ,axis.text.y = element_blank()
             ,axis.ticks = element_blank()
             ,legend.position='none'
             ,strip.text.x=element_text(face='bold')
           )
    )
  }

  # Finish up
  if(showPlot==TRUE){plot(gg)}
  return(gg)
}
