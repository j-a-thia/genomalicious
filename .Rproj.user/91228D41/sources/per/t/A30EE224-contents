#' Convert a \code{prcomp} object to a long-format \code{data.table}
#'
#' @param pcaObj Prcomp: An object produced by \code{prcomp}. If produced from
#' \code{genomalicious::pca_genos}, it may contain a \code{$pops} indexed
#' item, which contains the population designations for each sample in
#' \code{pcaObj$x}.
#'
#' @returns Returns a long-format data.table with columns \code{$POP}, the
#' population IDs, \code{$SAMPLE}, the sample IDs, \code{$AXIS}, the PC axis ID,
#' and \code{$SCORE}, the score on the PC axis.
#'
#' @examples
#' data(data_4pops)
#'
#' PCA <- pca_genos(data_4pops, popCol='POP')
#'
#' pca2DT(PCA)
#'
#' @export

pca2DT <- function(pcaObj){
  require(data.table); require(tidyverse)

  if(!class('prcomp')%in%class(pcaObj)){
    stop('Argument `pcaObj` must be of "prcomp" class.')
  }

  # Convert the PCA matrix into a data.frame
  Xtab <- pcaObj$x %>% as.data.frame() %>%
    rownames_to_column(., 'SAMPLE')

  # Check if $pops is in pcaObj, add in, and then convert to long-format
  if("pops" %in% names(pcaObj)){
    Xtab$POP <- pcaObj$pops
    Xtab %>%
      as.data.table %>%
      melt(., id.vars=c('POP','SAMPLE'), variable.name='AXIS', value.name='SCORE') %>%
      return()
  } else{
    Xtab %>%
      as.data.table %>%
      melt(., id.vars=c('POP','SAMPLE'), variable.name='AXIS', value.name='SCORE') %>%
      return()
  }
}
