#' Calculate allele counts
#'
#' @param dat Character/Integer: A vector of genotypes. If the class
#' is \code{'character'}, then assumes alleles are separated
#' with a '/'. Doesn't have to be biallelic (see param \code{biallelic}).
#' If the class is \code{'integer'}, then assumes counts of
#' Alt alleles, in which case, it assumes data is biallelic.
#'
#' @param biallelic Logical: Is the data biallelic? Affects processing and
#' output (see also "Value" section).
#'
#' @return A vector of alleles counts. If \code{biallelic==TRUE}, returns
#' vector with names \code{'ref'} and \code{'alt'}. If \code{biallelic==FALSE},
#' returns a vector with names as alleles.
#'
#' @examples
#' library(genomalicious)
#'
#' # Genotypes as separated alleles, biallelic
#' allele_counts(c('1/1', '0/1', '0/1', '0/0', '0/0'))
#'
#' # Genotypes as separated alleles, not biallelic
#' allele_counts(c('1/1', '2/3', '1/3', '0/0', '2/2'), biallelic=FALSE)
#'
#' # Genotypes as counts of the Alt allele
#' allele_counts(c(2, 1, 1, 0, 0))
#'
#' @export
allele_counts <- function(dat, biallelic=TRUE){

  # --------------------------------------------+
  # Libraries and assertins
  # --------------------------------------------+
  if(class(dat)=='character'){
    # Get the alleles
    als <- unlist(strsplit(dat, split='/'))
    # Get unique allelles
    uniqAls <- unique(als)
    # Check to see if number of unique alleles matches biallelic spec
    if(length(uniqAls)==2 & sum(c('1', '0') %in% uniqAls)==1 & biallelic==FALSE){
      warning('Argument biallelic==FALSE, but argument dat looks biallelic.')
    }
  }

  if(class(dat)=='numeric' | class(dat)=='integer'){
    # Make sure it's an integer
    dat <- as.integer(dat)
    # Set biallelic to TRUE
    biallelic <- TRUE
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
      alt <- sum(genoscore_converter(dat))
      ref <- n - alt
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
    alt <- sum(dat)
    ref <- n - alt
    return(c(ref=ref, alt=alt))
  }

}
