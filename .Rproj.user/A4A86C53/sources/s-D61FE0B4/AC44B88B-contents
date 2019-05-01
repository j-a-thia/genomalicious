#' Calculate allele counts
#'
#' @param dat Character/Integer: A vector of genotypes. If the class
#' is \code{'character'}, then assumes alleles are separated
#' with a '/'. If theclass is \code{'integer'}, then assumes counts of
#' Ref alleles.
#'
#' @return A vector of alleles counts for Ref ('0') and Alt ('1') alleles.
#'
#' @examples
#' allele_counts(c('1/1', '0/1', '0/1', '0/0', '0/0'))
#' allele_counts(c(2, 1, 1, 0, 2))
#'
#' @export
allele_counts <- function(dat){
  if(class(dat)=='numeric'){ dat <- as.integer(dat)}
  n <- length(dat)*2
  if(class(dat)=='character'){
    ref <- sum(genoscore_converter(dat))
    alt <- n - ref
  } else if(class(dat)=='integer'){
    ref <- sum(dat)
    alt <- n - ref
  }
  return(c(ref=ref, alt=alt))
}
