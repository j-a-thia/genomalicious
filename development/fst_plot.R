#' Plot pairwise FST estimates
#'
#' Takes the output \code{genomalicious::fst_calc} and produces a pairiwse plot of FST values.
#'
#' @param fstList List: Output from \code{genomalicious::fst_calc(..., pairwise=TRUE)}.
#'
#' @param plotColours Character: A vector of colours to use for the gradient of
#' FST values.
#'
#' @param legendPos Character: Where should the legend be positioned? Default is
#' \code{'top'}, but could also be one of, \code{'right'}, \code{'bottom'},
#' \code{'left'}, or \code{'none'}.
#'
#' @param scaleLims Numeric: The limits on the FST gradient. Default is NULL,
#' which will place a limit between 0 and the max FST.
#'
#' @return Returns a ggplot object.
#'
#' @examples
#' library(genomalicious)
#'
#' data(data_Genos)
#'
#' # Pairwise FST
#' genoPairs <- fst_calc(
#'    data_Genos,
#'    type='genos',
#'    global=FALSE,
#'    pairwise=TRUE
#'    )
#'
#' # Default pairwise FST plot
#' plot.fst.def <- fst_plot(genoPairs)
#' plot.fst.def
#'
#' # Custom pairwise FST plot, colours and legend positions.
#' plot.fst.custom.1 <- fst_plot(
#'    genoPairs,
#'    legendPos='right',
#'    plotColours=c('white', 'grey20')
#' )
#' plot.fst.custom.1
#'
#' # Custom pairwise FST plot, as above, but changing the scale.
#' plot.fst.custom.2 <- fst_plot(
#'    genoPairs,
#'    legendPos='right',
#'    plotColours=c('white', 'grey20'),
#'    scaleLims=c(0,1)
#' )
#' plot.fst.custom.2
#'
#' @export
fst_plot <- function(
    fstList, plotColours=NULL, legendPos='top', scaleLims=NULL){
  require(data.table); require(tidyverse)

  # Check arguments
  if(!(legendPos%in%c('top','right','left','bottom','none'))){
    stop('Argument `legendPos` must be a character value, one of
       "top", "right", "left", "bottom", or "none". See ?dapc_plot.')
  }

  if(is.null(plotColours)){
    plotColours <- c('white', '#0030C1', '#D6012C')
  }

  if(is.null(scaleLims)){
    scaleLims <- c(0,max(fstList$genome$FST))
  }

  # Sort the data
  fst.pairs.tab <- fstList$genome %>%
    pairwiseMat2DT(., flip=TRUE, X1='POP1', X2='POP2', Y='FST') %>%
    pairwiseMat2DT(., flip=FALSE, X1='POP1', X2='POP2', Y='FST')

  # Plot
  gg <- ggplot(fst.pairs.tab, aes(x=POP1, y=POP2, fill=FST)) +
    theme(
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.ticks.length.x = unit(1.5, 'mm'),
      axis.ticks.length.y = unit(1.5, 'mm'),
      legend.position=legendPos,
      axis.title=element_blank()
    ) +
    geom_tile(colour='grey20') +
    scale_x_discrete(expand=c(0,0)) +
    scale_y_discrete(expand=c(0,0)) +
    scale_fill_gradientn(
      colours=plotColours, limit=scaleLims,
      guide = guide_colorbar(frame.colour = "grey20", ticks.colour = "grey20")
    ) +
    labs(x='', y='', fill=bquote(italic('F')[ST]))

  # Output
  return(gg)
}


