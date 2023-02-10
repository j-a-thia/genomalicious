#' Convert a \code{prcomp} object to a long-format \code{data.table}
#'
#' @param pcaObj Prcomp: An object produced by \code{prcomp}. If produced from
#' \code{genomalicious::pca_genos}, it may contain a \code{$pops} indexed
#' item, which contains the population designations for each sample in
#' \code{pcaObj$x}.
#'
#' @param pops Character: An optional vector of population IDs. Default is \code{NULL}.
#' Be sure to check that they match samples in rows of \code{pcaObj$x}.
#'
#' @param subAxes Integer: An optional vector of values to subset the PC axes.
#' Default is NULL.
#'
#' @returns Returns a long-format data.table with columns \code{$POP}, the
#' population IDs (if data available), \code{$SAMPLE}, the sample IDs,
#' \code{$PC[x]}, the individual PC axes, with \code{[x]} denoting the axis number.
#'
#' @examples
#' library(genomalicious)
#' data(data_Genos)
#'
#' # With populations
#' PCA.pops <- pca_genos(data_Genos, popCol='POP')
#'
#' pca2DT(PCA.pops) %>% print
#'
#' # Without populations and all axes (default)
#' PCA.nopops <- pca_genos(data_Genos)
#'
#' pca2DT(PCA.nopops) %>% print
#'
#' # Manual populations, with just the first three axes
#' pops <- rownames(PCA.nopops$x) %>%
#'    sub('Ind', '', .) %>%
#'    sub('_.*', '', .)
#'
#' pops
#'
#' pca2DT(PCA.nopops, pops=pops, subAxes=1:3) %>% print
#'
#' @export

pca2DT <- function(pcaObj, pops=NULL, subAxes=NULL){
  require(data.table); require(tidyverse)

  if(!c('prcomp')%in%class(pcaObj)){
    stop('Argument `pcaObj` must be of "prcomp" class.')
  }

  if(is.null(subAxes)){
    subAxes <- 1:ncol(pcaObj$x)
  }

  # Convert the PCA matrix into a data.frame
  Xtab <- pcaObj$x[, subAxes] %>% as.data.frame() %>%
    rownames_to_column(., 'SAMPLE')

  # Check if $pops is in pcaObj, add in, and then convert to long-format
  if("pops" %in% names(pcaObj)){
    Xtab$POP <- pcaObj$pops
    Xtab %>%
      as.data.table %>%
      return()
  } else if(!is.null(pops)){
    Xtab %>%
      as.data.table %>%
      .[, POP:=pops] %>%
      return()
  } else{
    Xtab %>%
      as.data.table %>%
      return()
  }
}
