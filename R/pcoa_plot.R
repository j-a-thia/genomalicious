#' Plot PCoA results
#'
#' Plots results of a PCoA, e.g., scatterplot, screeplot, and cumulative
#' explained variance plots. Takes \code{pcoa} object as the main input.
#'
#' @param pcoaObj Prcomp object: A PCA of genotype data fitted using the
#' \code{prcomp} function. Either manually fitted, or using \code{genomalicious::pca_genos}.
#'
#' @param type Character: What type of plot to make: a scatterplot (\code{'scatter'}),
#' a screeplot of explained variances (\code{'scree'}), or the cumulative explained
#' variance (\code{'cumvar'}).
#'
#' @param axisIndex Integer: The PC axes to plot. If \code{type=='scatter'},
#' then must be exactly 2 values, the two PC axes to plot as a scatterplot.
#' If either \code{type=='scree'} or \code{type=='cumvar'}, then can be of
#' length from 1 to p, where p is the number of PC axes, and values again
#' represent the desired PC axes to plot.
#'
#' @param pops Character: A vector of population IDs, should match the
#' rows in \code{pcoaObj$x}, but is an optional argument. Default = \code{NULL}.
#' The function will search for \code{pcoaObj$pops} to assign to this argument
#' if not specified. Only valid when \code{type=='scatter'}.
#'
#' @param plotColours Character: A vector of colours to use for plotting,
#' but is an optional argument. Default = \code{NULL}.
#' When \code{type=='scatter'}, this must be a named list with one colour
#' per population. When \code{type=='scree'} or \code{type=='cumvar'}, only
#' a single colour is required, which is the colour of bars in the screeplot
#' or cumulative variance plot, respectively, and will default to 'grey20'
#' if unspecified.
#'
#' @param look Character: The look of the plot. Default = \code{'ggplot'}, the
#' typical gray background with gridlines produced by \code{ggplot2}. Alternatively,
#' when set to \code{'classic'}, produces a base R style plot.
#'
#' @param legendPos Character: Where should the legend be positioned? Default is
#' \code{'top'}, but could also be one of, \code{'right'}, \code{'bottom'},
#' \code{'left'}, or \code{'none'}.
#'
#' @return Returns a ggplot object.
#'
#' @examples
#' library(genomalicious)
#' data(data_PoolFreqs)
#' data(data_PoolInfo)
#'
#' # Note columns in data_PoolFreqs
#' colnames(data_PoolFreqs)
#'
#' # We need to add in the number of diploid individuals, $INDS
#' newFreqData <- left_join(data_PoolFreqs, data_PoolInfo)
#' head(newFreqData)
#'
#' # Fit the PCoA
#' PCOA <-pcoa_freqs(newFreqData)
#'
#' # Plot scatter with default settings
#' pcoa_plot(PCOA, type='scatter')
#'
#' # Plot scatter with custom colours and a classic look
#' pcoa_plot(
#'    PCOA,
#'    type='scatter',
#'    plotColours=c(Pop1='gray30', Pop2='royalblue', Pop3='palevioletred3', Pop4='plum2'),
#'    look='classic'
#' )
#'
#' # Explained variance
#' pcoa_plot(PCOA, type='scree')
#'
#' # Cumulative variance with custom colour
#' pcoa_plot(PCOA, type='cumvar', plotColours='royalblue')
#'
#' @export
pcoa_plot <- function(
    pcoaObj, type='scatter', axisIndex=NULL,
    plotColours=NULL, look='ggplot', legendPos='top'){

  # --------------------------------------------+
  # Libraries and assertions
  # --------------------------------------------+
  for(lib in c('data.table', 'ggplot2')){ require(lib, character.only = TRUE)}

  # Check the pcoaObj is the correct data class
  if(!'pcoa' %in% class(pcoaObj)){
    stop("Argument `pcoaObj` must be a prcomp class object.")
  }

  # Check that type is specified correctly
  if(!type %in% c('scatter', 'scree', 'cumvar')){
    stop("Argument `type` must be either: 'scatter', 'scree', or 'cumvar'.")
  }

  # Check that axisIndex is only length == 2
  if(type=='scatter' & length(axisIndex)>2){
    stop("Argument `axisIndex` should only contain two integer values for type=='scatter'.")
  }

  if(type%in%c('scree','cumvar') & sum(!axisIndex %in% 1:nrow(pcoaObj$values))){
    stop("Argument `axisIndex` should only contain values for the number of eigenvalues
         stored in pcoaObj$values for type=='scree' or type=='cumvar'.")
  }

  # Check that look is ggplot or classic.
  if(!look%in%c('ggplot', 'classic')){
    stop("Argument `look` is not one of: 'ggplot' or 'classic'.")
  }

  # Check that specified populations in plotColours are all in pops.
  if(type=='scatter' & is.null(plotColours)==FALSE){
    names.colours <- names(plotColours)
    names.pops <- pcoaObj$vectors %>% rownames()

    if(sum(names.pops %in% names.colours)!=length(names.pops)){
      stop("Argument `plotColours` misspecified: names of colours must be in
           rownames(pcoaObj$vectors).")
    }
  }

  # Specify axes if unassigned
  if(type=='scatter' & is.null(axisIndex)){
    axisIndex <- c(1,2)
  }

  if(type%in%c('scree','cumvar') & is.null(axisIndex)){
    axisIndex <- 1:nrow(pcoaObj$values)
  }

  # Assign colour if unspecified for scree and cumulative variance plots.
  if(type%in%c('scree','cumvar') & is.null(plotColours)){
    plotColours <- 'grey20'
  }

  # Set the plot theme by look
  if(look=='ggplot'){
    plotTheme <- theme_gray() + theme(legend.position=legendPos, axis.ticks.length = unit(0.2, 'cm'))
  } else if(look=='classic'){
    plotTheme <- theme_bw() + theme(
      panel.grid.major=element_blank()
      , panel.grid.minor=element_blank()
      , text=element_text(colour='black')
      , legend.position=legendPos
      , axis.ticks.length=unit(0.2, 'cm'))
  }

  # Make pcoaObj a data table of PC scores
  if(class(pcoaObj)=='pcoaObj'){
    plot.tab <- pcoaObj$vectors %>%
    as.data.frame() %>%
    rownames_to_column(., 'POP') %>%
    as.data.table()
  }

  # --------------------------------------------+
  # Code
  # --------------------------------------------+

  if(type=='scatter'){
    # Get axes
    axX <- paste0('Axis.', axisIndex[1])
    axY <- paste0('Axis.', axisIndex[2])

    # Percent explained variance
    eigvals <- pcoaObj$values$Eigenvalues
    varX <- round(eigvals[axisIndex[1]]/sum(eigvals) * 100, 2)
    varY <- round(eigvals[axisIndex[2]]/sum(eigvals) * 100, 2)

    # Create skeleton of plot
    gg <- ggplot(plot.tab, aes_string(x=axX, y=axY)) +
      plotTheme +
      labs(
        x=paste0('Axis ', axisIndex[1], ' (', varX, '%)')
        , y=paste0('Axis ', axisIndex[2], ' (', varY, '%)')
      )

    # Add points and population colours if specified
    if(is.null(plotColours)==TRUE){
      gg <- gg + geom_point(aes(colour=POP)) + labs(colour=NULL)
    } else if(is.null(plotColours)==FALSE){
      gg <- gg + geom_point(aes(colour=POP)) + scale_colour_manual(values=plotColours) + labs(colour=NULL)
    }
  }

  if(type %in% c('scree', 'cumvar')){
    # Vector of number PCs for X axis
    S <- pcoaObj$values$Eigenvalues
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
           + labs(x='PCo axes', y=axY)
    )
  }

  # Plot and return
  return(gg)
}
