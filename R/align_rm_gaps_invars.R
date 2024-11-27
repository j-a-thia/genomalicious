#' Remove gaps and invariant sites from an alignment matrix
#'
#' Provide a sequence alignment in matrix format and remove gaps and/or
#' invariant sites.
#'
#' @param alignMat Matrix: A sequence alignment converted into a matrix, with
#' positions in columns, samples in rows, and cells being DNA nucleotides or
#' amino acids. Row names should be names of samples.
#'
#' @param removeGaps Logical: Should gaps be removed? `TRUE` or `FALSE`.
#'
#' @param removeInvars Logical: Should invariant sites be removed? `TRUE` or `FALSE`.
#'
#' @details Gaps will be removed before invariant sites. Removal of gaps
#' occurs for any columns in `alignMat` that possess a `'_'`. Removal of invariant
#' sites occurs for any column in `alignMat` where `length(unique(x))` does
#' not evaluate to >2. Note, this will also include gaps if you do not remove them.
#'
#' @returns Returns a matrix with same number of rows as `alignMat`, but with any
#' columns removed that are gaps and/or invariant sites.
#' Columns are relabeled as `'pos.X'`, with 'X', being the relative
#' position in the original alignment.
#'
#' @export
#'
align_rm_gaps_invars <- function(alignMat, removeGaps, removeInvars){
  # Parameters
  n <- ncol(alignMat)
  m <- nrow(alignMat)
  colnames(alignMat) <- paste0('pos.',1:n)

  # Output alignment
  alignOut <- copy(alignMat)

  # Should gaps be removed?
  if(removeGaps==TRUE){
    gap.test <- apply(
      alignOut, 2,
      function(x) sum(x=='-')
    )

    alignOut <- alignOut[, names(gap.test)[gap.test==0]]
  }

  # Should invariant sites be removed?
  if(removeInvars==TRUE){
    var.test <- apply(
      alignOut, 2,
      function(x) length(unique(x))
    )

    alignOut <- alignOut[, names(var.test)[var.test>1]]
  }

  # Output
  return(alignOut)
}
