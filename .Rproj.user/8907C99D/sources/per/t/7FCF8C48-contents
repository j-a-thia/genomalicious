#' Convert a pairwise matrix into a pairwise data table
#'
#' This function takes pairwise matrix (all possible combinations in rows and
#' columns) and converts it into a long-format data table. This function can
#' also be flipped, i.e., going from a pairwise data table to a pairwise matrix.
#'
#' @param dat Matrix or Data table: The default parameterisation (\code{flip==FALSE})
#' is expecting a pairwise matrix. This matrix needs to be symmetrical, with a
#' zero on the diagonal and the off-diagonal containing the measured variable
#' values between pairs. There should be row and column names. If instead you
#' want to go in the flipped parameterisation (\code{flip==TRUE}), then a
#' data table is expected, where pair combinations are represented into two
#' columns and the measured variable is in a third column.
#'
#' @param flip Logical: Should the function be flipped? Default is \code{FALSE},
#' a pairwise matrix should be converted into a pairwise data table. If set to
#' \code{TRUE}, then a pairwise data table is converted into a matrix.
#'
#' @param X1 Character: The column name for the first member of the pair.
#' If \code{flip==FALSE}, then \code{X1} will take the row names from the matrix
#' and put them into the first column of the output data table.
#' If \code{flip==TRUE}, then \code{X1} will be take the column in the data table and
#' put it as row names in the output matrix.
#'
#' @param X2 Character: The column name for the first member of the pair.
#' If \code{flip==FALSE}, then \code{X2} will take the column names from the matrix
#' and put them into the second column of the output data table.
#' If \code{flip==TRUE}, then \code{X2} will be take the column in the data table and
#' put it as column names in the output matrix.
#'
#' @param Y Character: The column name or the measured variable between pairs.
#' If \code{flip==FALSE}, then \code{Y} will take the values in the off-diagonal
#' in the matrix and put them into the third column of the output data table.
#' If \code{flip==TRUE}, then \code{Y} will be take the column in the data table and
#' put it as values in the off-diagonal in the output matrix.
#'
#' @examples
#' data(data_Genos)
#'
#' # Calculate pairwise FST
#' pairFst <- fst_calc(data_Genos, type='genos', global=FALSE, pairwise=TRUE)
#'
#' # Convert into a pairwise matrix
#' pairMat <- pairwiseMat2DT(
#'    dat=pairFst$genome, flip=TRUE, X1='POP1', X2='POP2', Y='FST'
#' )
#' pairMat
#'
#' # Convert the pairwise matrix into a data table
#' pairDT <- pairwiseMat2DT(
#'    dat=pairMat, flip=FALSE, X1='POP1', X2='POP2', Y='FST'
#' )
#' pairDT
#'
#' @export

pairwiseMat2DT <- function(dat, flip=FALSE, X1, X2, Y){

  # --------------------------------------------+
  # Libraries and assertions
  # --------------------------------------------+

  require(data.table); require(tidyverse)

  if(flip==FALSE){
    if(!'matrix' %in% class(dat)){
      stop('Argument `dat` is not a matrix, but `type` is "FALSE", i.e.,
           turn a pairwise matrix into a pairwise data table. See ?pairwiseMat2DT.')
    }
  }

  if(flip==TRUE){
    dat <- as.data.table(dat)
    if(sum(c(X1, X2, Y) %in% colnames(dat))!=3){
      stop('Argument `flip` is "TRUE", i.e., turn a pairwise data table into a
           matrix, but not all the the required columns appear to be in `dat`.
           Check parameterisation of `X1`, `X2`, and `Y`. See ?pairwiseMat2DT.')
    }
  }

  # --------------------------------------------+
  # Code
  # --------------------------------------------+

  # Convert a pairwise matrix into a data table
  if(flip==FALSE){
    samp.uniq <- rownames(dat)

    output <- data.table(CJ(X1=samp.uniq, X2=samp.uniq), Y=0)

    for(i in 1:nrow(output)){
      x1 <- output$X1[i]
      x2 <- output$X2[i]
      output$Y[i] <- dat[x1, x2]
    }

    setnames(output, old=c('X1','X2','Y'), new=c(X1,X2,Y))
  }

  # Flipped: convert a pairiwse data table into matrix.
  if(flip==TRUE){
    samp.uniq <- c(dat[[X1]], dat[[X2]]) %>%  unique() %>%  sort()

    n <- length(samp.uniq)

    output <- matrix(0, ncol=n, nrow=n, dimnames=list(samp.uniq, samp.uniq))

    for(i in 1:nrow(dat)){
      x1.i <- dat[[X1]][i]
      x2.i <- dat[[X2]][i]
      measure.i <- dat[[Y]][[i]]
      output[x1.i, x2.i] <- measure.i
      output[x2.i, x1.i] <- measure.i
    }
  }

  # Output
  return(output)
}
