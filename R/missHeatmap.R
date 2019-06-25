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
#' @param plotColours Character: Vector of colours to use in plotting.
#' Size depends on values specified for parameters \code{type} and
#' \code{popCol}, see Details.
#'
#' @param plotNCol Integer: The number of columns to arrange indiviudal
#' population plots into. Only takes effect when \code{popCol} is specified.
#' Default = 2.
#'
#' @examples
#' ####   MISSING GENOTYPE DATA   ####
#' data(genomalicious_4pops)
#' datGt <- genomalicious_4pops
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
#' # Heatmaps, default and specified colours
#' missHeatmap(datGt)
#' missHeatmap(datGt, plotColours=c('black', 'plum2'))
#'
#' # Heatmaps, by population
#' missHeatmap(datGt, popCol='POP')
#'
#' ####   CATCH PLOT OUTPUT FOR LATER USE   ####
#' gg4pops <- missHeatmap(datGt, popCol='POP')
#' plot(gg4pops)
#'
#' @export
missHeatmap <- function(dat
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

  # --------------------------------------------+
  # Code
  # --------------------------------------------+

  # Make two colours, if only one specified
  if(length(plotColours)<2){
    plotColours <- c('white', 'royalblue')
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
    ggLs <- lapply(unique(dat$POP), function(pop){
      g <- (ggplot(dat[POP==pop], aes(x=SAMPLE, y=LOCUS))
            + geom_tile(aes(fill=as.factor(MISS)), colour=NA)
            + scale_fill_manual(values=c('0'=plotColours[1], '1'=plotColours[2]))
            + labs(x='Samples', y='Locus', title=pop)
            + theme(
              axis.text.x = element_blank()
              ,axis.text.y = element_blank()
              ,axis.ticks = element_blank()
              ,panel.border=element_rect(fill=NA, colour='black')
              ,plot.title=element_text(hjust=0.5, size=12, face='bold')
              ,legend.position='none'
            )
      )
      return(g)
    })
    gg <- do.call('grid.arrange', c(ggLs, ncol=plotNCol))
  }

  # Finish up
  plot(gg)
  return(gg)
}
