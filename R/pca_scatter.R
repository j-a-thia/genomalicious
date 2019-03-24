#' Create a scatter plot of a PCA on genotypes
#'
#' @export


pca_scatter <- function(dat, axisIndex=c(1,2), pops=NULL, popCols=NULL
                     , look='ggplot'){

  # --------------------------------------------+
  # Libraries and assertions
  # --------------------------------------------+
  for(lib in c('data.table', 'ggplot')){ require(lib, character.only = TRUE)}

  if(class(dat)!='prcomp' & !class(dat)[1]%in%c('data.table','data.frame','matrix')){
    stop('Argument dat must be one of the following object classes:
         prcomp, data.table, data.frame, or matrix')
  }

  if(!type%in%c('scatter', 'eigen')){
    stop("Argument scatter is not one of: 'scatter' or 'scree'.")
  }

  if(!look%in%c('ggplot', 'classic')){
    stop("Argument look is not one of: 'ggplot' or 'classic'.")
  }

  if(class(dat)=='prcomp'){
    if(is.null(dat$pops)==FALSE){ pops <- dat$pops}
  }

  if(is.null(pops)==FALSE &
     !sum(names(popCols)%in%unique(pops))==length(unique(pops))){
    stop("Argument popCols misspecified: names of colours must be in argument pops.")
  }

  if(look=='ggplot'){
    plot.theme <- theme_gray() + theme(legend.position='top')
  } else if(look=='classic'){
    plot.theme <- theme_bw() + theme(panel.grid.major=element_blank()
                                     , panel.grid.minor=element_blank()
                                     , text=element_text(colour='black')
                                     , legend.position='top')
  }

  if(class(dat)=='prcomp'){ plot.tab <- as.data.table(dat$x)
  } else{ plot.tab <- as.data.table(dat) }

  if(is.null(pops)==FALSE){
    plot.tab$POPS <- pops
  }

  # --------------------------------------------+
  # Code
  # --------------------------------------------+

  # Get axes
  axX <- colnames(plot.tab)[axisIndex[1]]
  axY <- colnames(plot.tab)[axisIndex[2]]

  # Create skeleton of plot
  gg <- ggplot(plot.tab, aes_string(x=axX, y=axY)) + plot.theme

  # Add points and population colours if specified
  if(is.null(pops)==TRUE){ gg <- gg + geom_point()
  } else if(is.null(pops)==FALSE & is.null(popCols)==TRUE){
    gg <- gg + geom_point(aes(colour=POPS)) + labs(colour=NULL)
  } else if(is.null(pops)==FALSE & is.null(popCols)==FALSE){
    gg <- gg + geom_point(aes(colour=POPS)) + scale_colour_manual(values=popCols) + labs(colour=NULL)
  }

  return(gg)
}
