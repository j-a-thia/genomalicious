#' Plot PCA results
#'
#' Plots different results of a PCA, e.g.: scatter, scree, and cumulative
#' explained vairance plots. See also \code{pca_DTinds()} to perform PCA on SNP genotypes.
#'
#' @param dat Prcomp/Data table/Data frame/Matrix or Numeric: The class of this
#' object depends on what plot you want to make, as specified by \code{type}.
#' If making a...
#' \itemize{
#' \item \strong{Scatter plot}: \code{dat} contains the PCA scores for
#' each individual. If not using a \code{prcomp} object, the tabular data must
#' be wide-format: rows are indivdiuals, columns are PC axes.
#' \item \strong{Scree plot} of eigenvalues, or plotting the \strong{cumulative
#' explained variance}: \code{dat} is simply a numeric vector containing the
#' eigenvalues of the PCA.
#' }
#'
#' @param type Character: What type of plot to make: scatter (\code{'scatter'})
#' the scree plot of eigenvalues (\code{'scree'}), or the cumulative explained
#' variance (\code{'expvar'}).
#'
#' @param axisIndex Integer: Vector of lenght = 2. Corresponds to the column
#' index to plots, i.e. the PC axes. Default = \code{c(1,2)}, the first and second
#' PC axis column tabulated in \code{dat}. Only valid when \code{type=='scatter'}.
#'
#' @param pops Character: A vector of population ID, should match the
#' rows in \code{dat}, but is an optional argument. Default = \code{NULL}.
#' If \code{dat} is a \code{prcomp} object, function will search for \code{dat$pops}
#' to assign to this argument. Only valid when \code{type=='scatter'}.
#'
#' @param popCols Character: A vector of colours to use for each unique population
#' in \code{pops}, but is an optional argument. Default = \code{NULL}.
#' The name of each colour must correspond to a population in \code{pops}.
#' Only valid when \code{type=='scatter'}.
#'
#' @param look Character: The look of the plot. Default = \code{'ggplot'}, the
#' typical gray background with gridlines produced by \code{ggplot2}. Alternatively,
#' when set to \code{'classic'}, produces a base R style plot.
#'
#' @return Prints the plot to screen, and also returns an object of
#' \code{gg}/\code{ggplot} class.
#'
#' @examples
#' # Data
#' data(genomalicious4pops)
#'
#' # Conduct the PCA with Patterson et al.'s (2006) normalisation, and
#' # population specified
#' pca <- pca_DTinds(dat=genomalicious4pops, scaling='patterson', popCol='POP')
#'
#' # Plot PCA scatter
#' pca_plot(pca)
#'
#' # Get more specific on scatter
#' pca_plot(pca
#'             , axisIndex=c(2,3)
#'             , pops=pca$pops
#'             , popCols=c(Pop1='gray30', Pop2='royalblue', Pop3='palevioletred3', Pop4='plum2')
#'             , look='classic')
#'
#' # Plot scree of eigenvalues
#' pca_plot(pca$sdev^2, type='scree')
#'
#' # Plot cumulative explained variance
#' pca_plot(pca$sdev^2, type='expvar', look='classic')
#'
#' @export

pca_plot <- function(dat, type='scatter', axisIndex=c(1,2)
                     , pops=NULL, popCols=NULL, look='ggplot'){

  # --------------------------------------------+
  # Libraries and assertions
  # --------------------------------------------+
  for(lib in c('data.table', 'ggplot2')){ require(lib, character.only = TRUE)}

  # Check the dat is the correct data class
  if(type=='scatter'){
    if(class(dat)!='prcomp' & !class(dat)[1]%in%c('data.table','data.frame','matrix')){
      stop('Argument dat must be one of the following object classes for making
            scatter plots: prcomp, data.table, data.frame, or matrix.')
    }
  } else if(type %in% c('scree', 'expvar')){
    if(class(dat)!='numeric'){
      stop('Argument dat must a numeric to make a scree or cumulative
           explained variance plot.')
    }
  }

  # Check that type is specified correctly
  if(!type %in% c('scatter', 'scree', 'expvar')){
    stop("Argument type must be either: 'scatter', 'scree', or 'expvar'.")
  }

  # Check that axisIndex is only length == 2
  if(length(axisIndex)>2){
    stop('Argument axisIndex should only contain two integer values.')
  }

  # Check that look is ggplot or classic.
  if(!look%in%c('ggplot', 'classic')){
    stop("Argument look is not one of: 'ggplot' or 'classic'.")
  }

  # If dat is a prcomp obect, and if there is a $pops index in dat,
  # assign the pops variable the $pops.
  if(class(dat)=='prcomp'){
    if(is.null(dat$pops)==FALSE){ pops <- dat$pops}
  }

  # Check that specified populations in popCols are all in pops.
  if(is.null(pops)==FALSE & is.null(popCols)==FALSE &
     !sum(names(popCols)%in%unique(pops))==length(unique(pops))){
    stop("Argument popCols misspecified: names of colours must be in argument pops.")
  }

  # Set the plot theme by look
  if(look=='ggplot'){
    plotTheme <- theme_gray() + theme(legend.position='top')
  } else if(look=='classic'){
    plotTheme <- theme_bw() + theme(panel.grid.major=element_blank()
                                     , panel.grid.minor=element_blank()
                                     , text=element_text(colour='black')
                                     , legend.position='top')
  }

  # Make dat a data table of PC scores
  if(class(dat)=='prcomp'){ plot.tab <- as.data.table(dat$x)
  } else{ plot.tab <- as.data.table(dat) }

  # If pops has been assigned, add this as a column to new dat
  if(is.null(pops)==FALSE){
    plot.tab$POPS <- pops
  }

  # --------------------------------------------+
  # Code
  # --------------------------------------------+

  if(type=='scatter'){
    # Get axes
    axX <- colnames(plot.tab)[axisIndex[1]]
    axY <- colnames(plot.tab)[axisIndex[2]]

    # Create skeleton of plot
    gg <- ggplot(plot.tab, aes_string(x=axX, y=axY)) + plotTheme

    # Add points and population colours if specified
    if(is.null(pops)==TRUE){ gg <- gg + geom_point()
    } else if(is.null(pops)==FALSE & is.null(popCols)==TRUE){
      gg <- gg + geom_point(aes(colour=POPS)) + labs(colour=NULL)
    } else if(is.null(pops)==FALSE & is.null(popCols)==FALSE){
      gg <- gg + geom_point(aes(colour=POPS)) + scale_colour_manual(values=popCols) + labs(colour=NULL)
    }
  }


  if(type %in% c('scree', 'expvar')){
    # Vector of number PCs for X axis
    X <- 1:length(dat)

    # If explained variance, divide eigenvalues by sum,
    # also create Y axis label
    if(type=='expvar'){
      Y <- 1- dat/sum(dat); axY <- 'Cumulative explained variance'
    } else{ Y <- dat; axY <- 'Eigenvalue'}

    # The plot
    gg <- (ggplot(data.frame(X, Y), aes(x=X, y=Y))
           + plotTheme
           + geom_line()
           + geom_point()
           + labs(x='PC axis', y=axY)
           )
  }

  # Plot and return
  plot(gg)
  return(gg)
}
