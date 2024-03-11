#' Plot DAPC results
#'
#' Plot the results of a DAPC. Includes options to produce scatterplot,
#' probability plots, and correct assignment heat maps. Expects output produced
#' from the function \code{genomalicious::dapc_fit}.
#'
#' @param dapcList List: A list object generated from \code{genomalicious::dapc_fit}.
#' See details.
#'
#' @param type Character: The type of plot to produce. If \code{'scatter'}, a
#' scatter plot is produced. If \code{'scree'}, a screeplot of the explained
#' among-population variation is produced. If \code{'cumvar'}, a plot of the
#' cumulative explained among-population variation is produced.
#' If \code{'probs'}, a plot of posterior probability
#' for populations identity is produced. If \code{'assign'}, a heatmap of
#' assignment rates among population pairs is produced.
#'
#' @param plotLook Character: The plot theme, only applicable when
#' \code{type=='scatter'}. Default = \code{'ggplot'}, the typical gray
#' background with gridlines produced by \code{ggplot2}. Alternatively,
#' when set to \code{'classic'}, produces a base R style plot.
#'
#' @param axisIndex Integer: The LD axes to plot. If \code{type=='scatter'},
#' there must be exactly 2 values, the two LD axes to plot as a scatterplot.
#' If \code{type=='scree'} or \code{type=='cumvar'}, these are the axes for
#' which you want to see the explained or cumulative among-population variance,
#' respectively. Default is NULL.
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
#' @param sampleOrder Logical: Should the samples be ordered by their IDs
#' ('by_id') or by their estimated  probability ('by_probs') in the probability
#' plot? Only applicable when \code{type=='probs'}. Default is \code{'by_id'}.
#'
#' @param plotColours Character: A vector of colours to use. If \code{type=='scatter'}
#' or \code{type=='probs'}, then this needs to be a named vector, where the values
#' are colours and their indexed names are populations. All populations need
#' to be represented with a colour. If \code{type=='assign'}, a series of colours
#' needs to be specified to parameterise the colour gradient in the assignment
#' rate heatmap. If \code{type=='scree'} or \code{type=='cumvar'}, a single value.
#'
#' @param legendPos Character: Where should the legend be positioned? Default is
#' \code{'top'}, but could also be one of, \code{'right'}, \code{'bottom'},
#' \code{'left'}, or \code{'none'}.
#'
#' @details If you want to produce a DAPC scatterplot (\code{type=='scatter'}) or
#' a probability plot (\code{type=='probs'}), then this function receives the
#' output of \code{dapc_fit(..., type='fit')}. If instead, you have performed
#' as assignment analysis with \code{dapc_fit(..., type='loo_cv')} or
#' \code{dapc_fit(..., type='traint_test')}, then you want to parameterise with
#' \code{type=='assign'}.
#'
#' In the probability plot is requested (\code{type=='probs'}), you can choose
#' to order samples by their estimated probabilities for their designated
#' population by setting \code{sampleOrder='by_probs'}. The default is
#' \code{sampleOrder='by_id'}, in which case, samples are ordered
#' alpha-numerically by populations and their sample ID, i.e., the command:
#' \code{data.table::setorder(.., POP, SAMPLE)}.
#'
#' @return Returns a ggplot object.
#'
#' @examples
#' library(genomalicious)
#'
#' data("data_Genos")
#'
#' # DAPC fit on all samples
#' DAPC.fit <- dapc_fit(data_Genos, pcPreds=3, method='fit')
#'
#' # DAPC using training and testing partitions
#' DAPC.tt <- dapc_fit(data_Genos, pcPreds=3, method='train_test')
#'
#' # Scatterplot, LD1 and LD2, with default colours, and ggplot look
#' dapc_plot(DAPC.fit, type='scatter', axisIndex=c(1,2))
#'
#' # Scatterplot, LD2 and LD3, with manual colours, and classic look
#' dapc_plot(DAPC.fit, type='scatter', axisIndex=c(2,3), plotLook='classic',
#' plotColours=c(Pop1='#08c7e0', Pop2='#4169e1', Pop3='#e46adf', Pop4='#ce0073')
#' )
#'
#' # Screeplot
#' dapc_plot(DAPC.fit, type='scree', )
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
    dapcList, type, plotLook='ggplot', axisIndex=NULL, popBarScale=1,
    sampleShow=TRUE, sampleOrder='by_id', plotColours=NULL, legendPos='top'
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
  if(!plotLook %in% c('ggplot','classic')){
    stop('Argument `plotLook` must be one of "ggplot" or "classic". See ?dapc_plot.')
  }

  if(!(legendPos%in%c('top','right','left','bottom','none'))){
    stop('Argument `legendPos` must be a character value, one of
       "top", "right", "left", "bottom", or "none". See ?dapc_plot.')
  }

  if(plotLook=='ggplot'){
    plotTheme <- theme_gray() + theme(legend.position=legendPos, axis.ticks.length = unit(0.2, 'cm'))
  } else if(plotLook=='classic'){
    plotTheme <- theme_bw() + theme(
      panel.grid.major=element_blank()
      , panel.grid.minor=element_blank()
      , text=element_text(colour='black')
      , legend.position=legendPos
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

  # Scatterplot axis
  if(type=='scatter'){
    if(is.null(axisIndex)==TRUE){
      axisIndex <- 1:2
    }
  }

  # Scree colours and axes
  if(type=='scree'){
    if(is.null(plotColours)==TRUE){
      plotColours <- 'grey20'
    }
    if(length(plotColours)!=1){
      stop('For `type`=="scree", argument `plotColours` must be a single value. See ?dapc_plot.')
    }
  }

  if(type=='scree'){
    if(is.null(axisIndex)==TRUE){
      axisIndex <- 1:length(dapcList$da.fit$svd)
    }
  }

  # Cumulative variance colours and axes
  if(type=='cumvar'){
    if(is.null(plotColours)==TRUE){
      plotColours <- 'grey20'
    }
    if(length(plotColours)!=1){
      stop('For `type`=="cumvar", argument `plotColours` must be a single value. See ?dapc_plot.')
    }
  }

  if(type=='cumvar'){
    if(is.null(axisIndex)==TRUE){
      axisIndex <- 1:length(dapcList$da.fit$svd)
    }
  }

  # Probability plots
  if(type=='probs'){
    # Population scale bar class
    if(class(popBarScale)!='numeric'){
      stop('Argument `popBarScale` must be a numeric value. See ?dapc_plot.')
    }
    # Population scale bar values
    if(popBarScale<1){
      stop('Argument `popBarScale` must be >=1. See ?dapc_plot.')
    }
    # Sample label class
    if(class(sampleShow)!='logical'){
      stop('Argument `sampleShow` must be a logical value. See ?dapc_plot.')
    }
    # Sample order class
    if(!sampleOrder %in% c('by_id','by_probs')){
      stop('Argument `sampleOrder` must be either "by_id" or "by_probs". See ?dapc_plot.')
    }
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

  # Assignment plots
  if(type=='assign'){
    if(!'pairs.long' %in% names(dapcList)){
      stop('For `type==assign`, there must be index `dapcList$pairs.long`. See dapc_plot?')
    }
    if(!is.null(plotColours) & length(plotColours)<2){
      stop('Argument `plotColours` requires 2 or more colours.')
    }
    if(is.null(plotColours)){
      plotColours <- c('white', '#0030C1', '#D6012C')
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
      theme(legend.position=legendPos) +
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

  ### Plot scree or cumulative variance
  if(type %in% c('scree', 'cumvar')){
    # Vector of number PCs for X axis
    S <- dapcList$da.fit$svd^2
    X <- 1:length(S)

    # If explained variance, divide eigenvalues by sum,
    # also create Y axis label
    if(type=='cumvar'){
      Y <- unlist(lapply(1:length(S), function(i){
        sum(S[1:i])/sum(S) * 100
      }))
      axY <- 'Cumulative variance (%)'
    } else if(type=='scree'){
      Y <- S/sum(S) * 100
      axY <- 'Explained variance (%)'
    }

    # The plot
    gg <- (data.frame(X=X[axisIndex], Y=Y[axisIndex]) %>%
             ggplot(., aes(x=X, y=Y))
           + plotTheme
           + geom_col(fill=plotColours)
           + scale_x_continuous(breaks = ~round(unique(pretty(.))))
           + labs(x='LD axes', y=axY)
    )
  }

  ### Plot probabilities
  if(type=='probs'){
    plot.tab <- dapcList$da.prob

    # Should the samples be ordered by their posterior probabilities?
    if(sampleOrder=='by_probs'){
      samp.order <- plot.tab %>%
        # Split by pop
        split(., by='POP') %>%
        # Iterate through pop, and order by probability
        lapply(., function(D){
          # The focal population
          pop <- D$POP[1]
          # Get the predictions for the focal population
          D[POP.PRED==pop] %>%
            # Organise in descending order
            setorder(., -PROB) %>%
            # Output the designated population and sample ID
            .[, c('POP','SAMPLE')]
        }) %>%
        # Combine
        do.call('rbind', .) %>%
        # Add in an order ID
        .[, ORDER:=1:.N]
    } else if(sampleOrder=='by_id'){
      samp.order <- plot.tab[, c('POP','SAMPLE')] %>%
        unique() %>%
        setorder(., POP, SAMPLE) %>%
        .[, ORDER:=1:.N]
    }

    # Combine data with sample order
    plot.tab <- left_join(plot.tab, samp.order) %>%
      as.data.table %>%
      setorder(., ORDER)

    # Plot the populations
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
        ymax=-0.02,
        ymin=-0.04*popBarScale,
        inherit.aes=FALSE
      ) +
      scale_x_continuous(
        expand=c(0,0),
        breaks=1:max(samp.order$ORDER),
        labels=samp.order$SAMPLE
      ) +
      scale_y_continuous(
        expand=c(0,0),
        limits=c(-0.04*popBarScale, 1.01)
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
      scale_fill_gradientn(
        colours=plotColours,
        guide = guide_colorbar(frame.colour = "grey20", ticks.colour = "grey20")
      ) +
      labs(x='Observed', y='Predicted', fill='Assignment rate')
  }

  ### Output
  return(gg)
}
