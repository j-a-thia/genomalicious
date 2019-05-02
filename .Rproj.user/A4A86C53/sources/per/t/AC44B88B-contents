#' Calculate allele counts
#'
#' @param dat Character/Integer: A vector of genotypes. If the class
#' is \code{'character'}, then assumes alleles are separated
#' with a '/'. Doesn't have to be biallelic (see param \code{biallelic}).
#' If the class is \code{'integer'}, then assumes counts of
#' Ref alleles, in which case, it does assumes biallelic.
#'
#' @param biallelic Logical: Is the data biallelic? If alleles separated
#' by '/' (character vector), expect \code{dat} to take form: '0/0', '0/1', or '1/1'.
#' If alleles coded as counts of Ref allele, expect \code{dat} to take form:
#' 0, 1, or 2. Default =\code{TRUE}.
#'
#' @return A vector of alleles counts.
#'
#' @examples
#' allele_counts(c('1/1', '0/1', '0/1', '0/0', '0/0'))
#' allele_counts(c('1/1', '2/3', '1/3', '0/0', '2/2'), biallelic=FALSE)
#' allele_counts(c(2, 1, 1, 0, 2))
#'
#' @export
allele_counts <- function(dat, biallelic=TRUE){

  # --------------------------------------------+
  # Libraries and assertins
  # --------------------------------------------+
  # Get the alleles
  als <- unlist(strsplit(dat, split='/'))
  # Get unique allelles
  uniqAls <- unique(als)

  if(class(dat)=='numeric'){ dat <- as.integer(dat)}

  if(length(uniqAls)==2 & sum(c('1', '0') %in% uniqAls)==1 & biallelic==FALSE){
    warning('Argument biallelic==FALSE, but argument dat looks biallelic.')
  }

  if(length(uniqAls)>2){ biallelic <- FALSE }

  if(class(dat)=='integer' & biallelic==FALSE){
    stop('Argument dat is an integer, and argument biallelic==FALSE.
    Cannot count alleles in dat if there are >2 alleles.')
  }

  # --------------------------------------------+
  # Code
  # --------------------------------------------+
  # Sample size
  n <- length(dat)*2

  # If given a character vector:
  if(class(dat)=='character'){
    # If data is biallelic
    if(biallelic==TRUE){
      ref <- sum(genoscore_converter(dat))
      alt <- n - ref
      return(c(ref=ref, alt=alt))
    # Else, if >2 alleles
    } else{
      tab <- table(als)
      vec <- as.vector(tab)
      names(vec) <- names(tab)
      return(vec)
    }
  # If given counts of bialleles
  } else if(class(dat)=='integer'){
    ref <- sum(dat)
    alt <- n - ref
    return(c(ref=ref, alt=alt))
  }

}
