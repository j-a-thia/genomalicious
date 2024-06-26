#' Convert a \code{lda} object to a long-format \code{data.table}
#'
#' @param daObj Lda: An object produced by \code{lda}.
#'
#' @param sampVec Character: The sample IDs, which need to match the sample
#' order in \code{predict(daObj)$x}.
#'
#' @param obsPops Character/Factor: The observed population IDs, which need to
#' match sample order in \code{predict(daObj)$x}.
#'
#' @return Returns a wide-format data.table with columns \code{$POP}, the
#' population IDs, \code{$SAMPLE}, the sample IDs, \code{$LD[x]}, the individual
#' LD axes, with \code{[x]} denoting the axis number.
#'
#' @examples
#' library(genomalicious)
#'
#' data(data_Genos)
#'
#' PCA <- pca_genos(data_Genos, popCol='POP')
#' DA <- lda(PCA$pops ~ PCA$x[,1:3])
#'
#' da2DT(DA, sampVec=rownames(PCA$x), obsPops=PCA$pops)
#'
#' @export

da2DT <- function(daObj, sampVec, obsPops){
  require(data.table); require(tidyverse)

  if(!c('lda')%in%class(daObj)){
    stop('Argument `daObj` must be of "lda" class.')
  }

  # Get DA predictions
  DA.pred <- predict(daObj)

  # Add in populations then convert to long-format
  data.table(POP=obsPops, SAMPLE=sampVec, DA.pred$x) %>%
    return()
}
