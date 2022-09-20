#' Plot DAPC results
#'
#' Plot the results of a DAPC. Includes options to produce scatterplot,
#' probability plots, and correct assignment heat maps. Expects output produced
#' from the function \code{genomalicious::dapc_fit}.
#'
#' @param dapcList List: A list object generated from \code{genomalicious::dapc_fit}.
#' Depending on what plot is requested, the contents of \code{dapcList} will change.
#' See details.
#'
#' @param type Character: The type of plot to produce. If \code{'scatter'}, a
#' scatter plot is produced. If \code{'probs'}, a plot of posterior probability
#' for populations identity is produced. If \code{'assign'}, a heatmap of
#' assignment rates among population pairs is produced.
#'
#' @param scatterLook Character: The plot theme, only applicable when
#' \code{type=='scatter'}. Default = \code{'ggplot'}, the typical gray
#' background with gridlines produced by \code{ggplot2}. Alternatively,
#' when set to \code{'classic'}, produces a base R style plot.
#'
#' @param axisIndex Integer: The LD axes to plot, only applicable when
#' \code{type=='scatter'}. There must be exactly 2 values, the two LD axes
#' to plot as a scatterplot.
#'
#' @param popBarScale Numeric: A scaling value for the population guide bar,
#' only applicable when \code{type=='probs'}. A guide bar is used to delineate
#' the original designated populations. The thickness of this bar is controlled
#' with \code{popBarScale}. Default is 1, and the value must be >= 1.
#'
#' @param sampleShow Logical: Should the sample names be displayed in the
#' probability plot? Only applicable when \code{type=='probs'}.
#' Default is \code{TRUE}.
#'
#' @param plotColours Character: A vector of colours to use. If \code{type=='scatter'}
#' or \code{type=='probs'}, then this needs to be a named vector, where the values
#' are colours and their indexed names are populations. All populations need
#' to be represented with a colour. If \code{type=='assign'}, a series of colours
#' needs to be specified to parameterise the colour gradient in the assignment
#' rate heatmap.
#'
#' @examples
#' library(genomalicious)
#'
#' data("data_4pops")
#'
#' # DAPC fit on all samples
#' DAPC.fit <- dapc_fit(data_4pops, pcPreds=3, method='fit')
#'
#' # DAPC using training and testing partitions
#' DAPC.tt <- dapc_fit(data_4pops, pcPreds=3, method='train_test')
#'
#' # Scatterplot, LD1 and LD2, with default colours, and ggplot look
#' dapc_plot(DAPC.fit, type='scatter', axisIndex=c(1,2))
#'
#' # Scatterplot, LD2 and LD3, with manual colours, and classic look
#' dapc_plot(DAPC.fit, type='scatter', axisIndex=c(2,3), scatterLook='classic',
#' plotColours=c(Pop1='#08c7e0', Pop2='#4169e1', Pop3='#e46adf', Pop4='#ce0073')
#' )
#'
#' # Probability plot, default colours
#' dapc_plot(DAPC.fit, type='probs')
#'
#' # Probability plot, manual colours, population bar rescaled, sample names
#' turned off, and legend turned off.
#' dapc_plot(DAPC.fit, type='probs', popBarScale=5,
#' plotColours=c(Pop1='#08c7e0', Pop2='#4169e1', Pop3='#e46adf', Pop4='#ce0073'),
#' sampleShow=FALSE, legendPos='none'
#' )
#'
#' # Assignment heatmap, default colours
#' dapc_plot(DAPC.tt, type='assign')
#'
#' # Assignment heatmap, manual colours, and legend repositioned to top
#' dapc_plot(DAPC.tt, type='assign', plotColours=c('white', 'grey50', 'grey20'),
#' legendPos='top')
#'
#' @export
dapc_plot <- function(
  dapcList, type, scatterLook='ggplot', axisIndex=c(1,2), popBarScale=1,
  sampleShow=TRUE, plotColours=NULL, legendPos='right', showPlot=TRUE
){
  # --------------------------------------------+
  # Libraries and assertions
  # --------------------------------------------+
  for(lib in c('data.table', 'ggplot2', 'tidyverse')){ require(lib, character.only = TRUE)}

  # Check that the input is a list
  if('list'!=class(dapcList)){
    stop(
    'Argument `dapcList` must be a list, and specifically generated using
    the function, dapc_fit. See ?dapc_plot.')
  }

  # Parameter checks
  if(!scatterLook %in% c('ggplot','classic')){
    stop('Argument `scatterLook` must be one of "ggplot" or "classic". See ?dapc_plot.')
  }

  if(class(popBarScale)!='numeric'){
    stop('Argument `popBarScale` must be a numeric value. See ?dapc_plot.')
  }

  if(popBarScale<1){
    stop('Argument `popBarScale` must be >=1. See ?dapc_plot.')
  }

  if(class(sampleShow)!='logical'){
    stop('Argument `sampleShow` must be a logical value. See ?dapc_plot.')
  }

  if(!(legendPos%in%c('top','right','left','bottom','none'))){
    stop('Argument `legendPos` must be a character value, one of
       "top", "right", "left", "bottom", or "none". See ?dapc_plot.')
  }

  if(scatterLook=='ggplot'){
    plotTheme <- theme_gray() + theme(legend.position='top', axis.ticks.length = unit(0.2, 'cm'))
  } else if(scatterLook=='classic'){
    plotTheme <- theme_bw() + theme(
      panel.grid.major=element_blank()
      , panel.grid.minor=element_blank()
      , text=element_text(colour='black')
      , legend.position='top'
      , axis.ticks.length=unit(0.2, 'cm'))
  }

  # Specific checks for each type
  if(type=='scatter'){
    # Check table
    if(!'da.tab'%in%names(dapcList)){
      stop(
        'For `type=="scatter"`, the argument `dapcList` needs the index
      `dapcList$da.tab`. See ?dapc_plot.'
      )
    }
    # Colours
    if(!is.null(plotColours)){
      pops.uniq <- dapcList$da.tab$POP %>% unique
      pops.k <- length(pops.uniq)
      if(sum(names(plotColours) %in% pops.uniq)!=pops.k){
        stop(
          'For `type=="scatter"`, argument `plotColours` must be a named
        character vector with all names matching the populations in
        `dapcList$da.tab$POP`. See ?dapc_plot.'
        )
      }
    }
  }

  if(type=='probs'){
    # Check table
    if(!'da.prob'%in%names(dapcList)){
      stop(
        'For `type=="probs"`, the argument `dapcList` needs the index
      `dapcList$da.prob`. See ?dapc_plot.'
      )
    }
    # Colours
    if(!is.null(plotColours)){
      pops.uniq <- dapcList$da.prob$POP %>% unique
      pops.k <- length(pops.uniq)
      if(sum(names(plotColours) %in% pops.uniq)!=pops.k){
        stop(
          'For `type=="probs"`, argument `plotColours` must be a named
        character vector with all names matching the populations in
        `dapcList$da.prob$POP`. See ?dapc_plot.'
        )
      }
    }
  }

  if(type=='assign'){
    if(!'pairs.long' %in% names(dapcList)){
      stop('For `type==assign`, there must be index `dapcList$pairs.long`. See dapc_plot?')
    }
    if(!is.null(plotColours) & length(plotColours)<2){
      stop('Argument `plotColours` requires 2 or more colours.')
    }
    if(is.null(plotColours)){
      plotColours <- c('#0030C1', '#EA49DF', '#D6012C')
    }
  }

  # --------------------------------------------+
  # Code
  # --------------------------------------------+
  ### Plot scatter
  if(type=='scatter'){
    plot.tab <- dapcList$da.tab

    # Get axes
    axX <- paste0('LD',axisIndex[1])
    axY <- paste0('LD',axisIndex[2])

    # Percent explained variance
    eigvals <- dapcList$da.fit$svd^2
    varX <- round(eigvals[axisIndex[1]]/sum(eigvals) * 100, 2)
    varY <- round(eigvals[axisIndex[2]]/sum(eigvals) * 100, 2)

    # Create skeleton of plot
    gg <- ggplot(plot.tab, aes_string(x=axX, y=axY, colour='POP')) +
      plotTheme +
      geom_point() +
      stat_ellipse(type='norm') +
      labs(
        x=paste0('LD', axisIndex[1], ' (', varX, '%)')
        , y=paste0('LD', axisIndex[2], ' (', varY, '%)')
        , colour=NULL
      )

    # Add points and population colours if specified
    if(is.null(plotColours)==FALSE){
      gg <- gg + scale_colour_manual(values=plotColours) + labs(colour=NULL)
    }
  }

  ### Plot probabilities
  if(type=='probs'){
    plot.tab <- dapcList$da.prob

    samp.order <- plot.tab %>%
      .[, c('POP','SAMPLE')] %>%
      unique %>%
      setorder(., POP, SAMPLE) %>%
      .[, ORDER:=(1:.N)]

    plot.tab <- left_join(plot.tab, samp.order) %>% as.data.table

    plot.pops <- left_join(
      plot.tab[, .(MIN=(min(as.integer(ORDER)))-0.5), by=POP],
      plot.tab[, .(MAX=(max(as.integer(ORDER)))+0.5), by=POP],
    )

    # Create skeleton of plot
    gg <- ggplot() +
      theme(
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_text(angle=30, hjust=1),
        axis.ticks.length.x = unit(1.5, 'mm'),
        axis.ticks.length.y = unit(1.5, 'mm'),
        legend.position=legendPos
      ) +
      geom_col(
        data=plot.tab,
        mapping=aes(x=ORDER, y=PROB, fill=POP.PRED)
      ) +
      geom_rect(
        data=plot.pops,
        mapping=aes(xmin=MIN, xmax=MAX, fill=POP),
        ymax=-0.01,
        ymin=-0.025*popBarScale,
        inherit.aes=FALSE
      ) +
      scale_x_continuous(
        expand=c(0,0),
        breaks=1:max(samp.order$ORDER),
        labels=samp.order$SAMPLE
      ) +
      scale_y_continuous(
        expand=c(0,0),
        limits=c(-0.025*popBarScale, 1.1)
      ) +
      labs(x='Samples', y='Probability', fill=NULL)

    # Add points and population colours if specified
    if(is.null(plotColours)==FALSE){
      gg <- gg + scale_fill_manual(values=plotColours) + labs(colour=NULL)
    }

    # Remove sample names if specified
    if(sampleShow==FALSE){
      gg <- gg + theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
      )
    }
  }

  ### Plot assignment rates
  if(type=='assign'){
    plot.tab <- dapcList$pairs.long

    gg <- ggplot(plot.tab, aes(x=POP, y=POP.PRED, fill=ASSIGN)) +
      theme(
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.ticks.length.x = unit(1.5, 'mm'),
        axis.ticks.length.y = unit(1.5, 'mm'),
        legend.position=legendPos
      ) +
      geom_tile(colour='grey20') +
      scale_x_discrete(expand=c(0,0)) +
      scale_y_discrete(expand=c(0,0)) +
      scale_fill_gradientn(colours=plotColours) +
      labs(x='Observed', y='Predicted', fill='Assignment rate')
  }

  ### Output
  plot(gg)
  return(gg)
}


