#' Convert a pairwise matrix into a pairwise data table
#'
#' This function takes pairwise matrix (all possible combinations in rows and
#' columns) and converts it into a long-format data table. This function can
#' also be flipped, i.e., going from a pairwise data table to a pairwise matrix.
#'
#' @param dat Matrix or Data table: The default parameterisation (\code{flip==FALSE})
#' is expecting a pairwise matrix. This matrix needs to be symmetrical, with a
#' values on the diagonal representing comparisons within a subject, and values on
#' the off-diagonal containing the measures between subject pairs.
#' There should be row and column names. If instead you want to go in the
#' flipped parameterisation (\code{flip==TRUE}), then a data table is expected,
#' where pair combinations are represented into two columns and the measured
#' variable is in a third column. In this second parameterisation, the data table
#' does not need to have within-subject comparisons (diagonals; you can add these in).
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
#' @param diagAdd Logical: Should the diagonal be added? Default is
#' \code{TRUE}. Only takes effect if \code{flip==TRUE}.
#'
#' @param diagVal Numeric: The value to add along the diagonal when
#' \code{diagAdd==TRUE}. Default is 0.
#'
#' @examples
#' library(genomalicious)
#' data(data_Genos)
#'
#' # Calculate pairwise FST
#' pairFst <- fstat_calc(data_Genos, type='genos', fstat='FST', global=FALSE, pairwise=TRUE)
#'
#' # Convert into a pairwise matrix
#' pairMat <- pairwiseMat2DT(
#'    dat=pairFst$genome, flip=TRUE, X1='POP1', X2='POP2', Y='FST',
#'    diagAdd=TRUE, diagVal=0
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

pairwiseMat2DT <- function(dat, flip=FALSE, X1, X2, Y, diagAdd=TRUE, diagVal=0){

  # --------------------------------------------+
  # Libraries and assertions
  # --------------------------------------------+

  require(data.table); require(tidyverse)

  if(flip==FALSE){
    if(!'matrix' %in% class(dat)){
      stop('Argument `dat` is not a matrix, but `type` is "FALSE", i.e.,
           turn a pairwise matrix into a pairwise data table. See ?pairwiseMat2DT.')
    }
    if(sum(is.na(dat))>0){
      stop('Argument `dat` has NAs. If there is missing data, please encode it in a different way.')
    }
  }

  if(flip==TRUE){
    dat <- as.data.table(dat)

    # Check for appropriate columns
    if(sum(c(X1, X2, Y) %in% colnames(dat))!=3){
      stop('Argument `flip` is "TRUE", i.e., turn a pairwise data table into a
           matrix, but not all the the required columns appear to be in `dat`.
           Check parameterisation of `X1`, `X2`, and `Y`. See ?pairwiseMat2DT.')
    }

    # Check to make sure diagonal parameterisation is correct. If the data.table
    # already contains information for the diagonal, then choosing `diagAdd==TRUE`
    # will cause an error.
    within <- sum(dat[[X1]]==dat[[X2]])
    if(within>0 & diagAdd==TRUE){
      stop_msg <- paste0(
        'Argument `diagAdd` is "TRUE", but there are ', within,
        ' within-subject comparisons: ' ,
        '`diagAdd` should only be "TRUE" if there are no within-subject ',
        'measurements. See ?pairwiseMat2DT'
      )
      stop(stop_msg)
    }
  }

  # --------------------------------------------+
  # Code
  # --------------------------------------------+

  # Convert a pairwise matrix into a data table
  if(flip==FALSE){
    dat[upper.tri(dat, diag=FALSE)] <- NA
    output <- dat %>%
      as.data.frame() %>%
      rownames_to_column(., 'X1') %>%
      as.data.table %>%
      melt(., id.vars='X1', variable.name='X2', value.name='Y') %>%
      .[!is.na(Y)] %>%
      setnames(., old=c('X1','X2','Y'), new=c(X1,X2,Y))
  }

  # Flipped: convert a pairiwse data table into matrix.
  if(flip==TRUE){
    dat <- dat %>%
      copy %>%
      setnames(., old=c(X1,X2,Y), new=c('X1','X2','Y'))

    samp.uniq <- dat[,c('X1','X2')] %>%
      unlist %>%
      unique %>%
      sort

    if(diagAdd==TRUE){
      dat2cast <- rbind(
        # Data in X1-X2 orientation
        dat[, c('X1','X2','Y')],
        # Data in X2-X1 orientation
        data.table(X1=dat$X2, X2=dat$X1, Y=dat$Y),
        # The diagonal
        data.table(X1=samp.uniq, X2=samp.uniq, Y=diagVal)
      )
    } else if(diagAdd==FALSE){
      dat2cast <- dat[, c('X1','X2','Y')]
    }

    output <- dat2cast %>%
      data.table::dcast(., X1 ~ X2, value.var = 'Y') %>%
      as.data.frame %>%
      column_to_rownames(., 'X1') %>%
      as.matrix()
  }

  # Output
  return(output)
}
