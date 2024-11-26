#' Remove gaps from an alignment matrix
#'
#' Provide a sequence alignment in matrix format and remove gaps ('-').
#'
#' @param align Matrix: A sequence alignment converted into a matrix, with
#' positions in columns, samples in rows, and cells being DNA nucleotides or
#' amino acids. Rowsnames should be names of samples.
#'
#' @returns Returns a matrix with same number of rows as `align`, but with any
#' columns with `'-'` being removed. Columns are relabelled as `'pos.X'`, with
#' 'X', being the relative position in the original alignment.
#'
#' @export
#'
align_mat_rm_gaps <- function(align){
  n <- ncol(align)
  m <- nrow(align)

  align.rm.gaps <- matrix(
    NA, nrow=m, ncol=0,
    dimnames=list(rownames(align),NULL)
  )

  for(i in 1:n){
    x <- align[,i]
    if(sum(grepl('-',x))==0){
      new.col <- matrix(
        x, nrow=m, ncol=1,
        dimnames=list(names(x), paste0('pos.',i))
      )

      align.rm.gaps <- cbind(align.rm.gaps, new.col)
    }
  }

  return(align.rm.gaps)
}
