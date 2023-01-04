#' Calculate variance components for F-statistics
#'
#' Takes a list of values of allele frequencies, sample sizes, number of
#' populations, and heterozygosity, and returns variance components, as
#' per Weir & Cockerham (1984). Assumes that all values are for a single
#' biallelic SNP locus.
#'
#' Note, this function is not exported.
#'
#' @param pi Numeric: A vector of allele frequencies for each population.
#'
#' @param ni Numeric: A vector of sample sizes (number of diploid individuals)
#' for each population.
#'
#' @param r Integer: A single value, the number of populations.
#'
#' @param hetStand Logical: Should the estimates be standardised for observed
#' heterozygosity?
#'
#' @param hi Numeric: A vector of observed heterozygosities for each population.
#'
#' @returns Returns a data table. If \code{hetStand==FALSE}, then the list has three
#' columns: \code{$MSP}, mean sqaures for populations; \code{$MSG}, mean squares
#' for gametes; \code{$Nc}, a sample size constant. If \code{hetStand==TRUE},
#' then the three columns are \code{$A}, \code{$B}, \code{$C}, which
#' correspond to the 'a', 'b' and 'c' variance components described in
#' Weir & Cockheram (1984).
#'
#' @references
#' Weir & Cockerham (1984) Evolution. DOI: 10.1111/j.1558-5646.1984.tb05657.x
#' Weir et al. (2002) Annals of Human Genetics. DOI: 10.1146/annurev.genet.36
#'
#' @examples
#' fstat_varcomps(pi=c(0.5,0.25), ni=c(20,20), r=2, hetStand=FALSE)
#'
#' fstat_varcomps(pi=c(0.5,0.25), ni=c(20,20), r=2, hi=c(0.05,0.375), hetStand=TRUE)

fstat_varcomps <- function(pi, ni, r, hi=NULL, hetStand=FALSE){
  require(data.table); require(tidyverse)

  if(hetStand==FALSE){
    # Mean weighted allele frequency
    p.mean <- sum(ni * pi)/sum(ni)

    # Mean squares variance components
    msp <- (1/(r-1)) * sum(ni * (pi - p.mean)^2)
    msg <- (1/sum(ni-1)) * sum(ni * pi * (1-pi))

    # Sample size correction factor
    nc <- (1/(r-1)) * ( sum(ni) - (sum(ni^2)/sum(ni)) )

    # Output as data.table
    return(data.table(MSP=msp, MSG=msg, Nc=nc))
  }

  if(hetStand==TRUE){
    # The mean sample size
    n.mean <- sum(ni/r)

    # The sample size scaling parameter
    nc <- (r*n.mean - sum((ni^2)/(r*n.mean))) / (r-1)

    # The average sample allele frequency
    p.mean <- sum((ni*pi)/(r*n.mean))

    # The variance in allele frequencies
    s2 <- sum( (ni*(pi-p.mean)^2)/((r-1)*n.mean) )

    # The average heterozygosity
    h.mean <- sum( (ni*hi)/(r*n.mean) )

    # The a, b, and c components
    a <- (nc/n.mean) * (s2 - (1/(n.mean-1))*((p.mean*(1-p.mean)) - (s2*(r-1)/r) - (0.25*h.mean)))

    b <- (n.mean/(n.mean-1)) * ((p.mean*(1-p.mean)) - (s2*(r-1)/r) - (h.mean*((2*n.mean-1)/(4*n.mean))))

    c <- 0.5 * h.mean

    # Return as list
    return(data.table(A=a, B=b, C=c))
  }
}
