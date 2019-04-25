#' Calculate the variance components of Weir & Cockerham's FST
#'
#' @param freqMat Matrix: Ref allele counts. Rows = populations,
#' columns = loci; make sure both are named. Row names used to label output FST matrix.
#'
#' @param sampMat Matrix: Number of sampled individuals. Rows = populations,
#' columns = loci.
fstWC_varcomps <- function(freqMat, sampleMat){
  lociNames <- colnames(freqMat)
  numPops <- nrow(freqMat)

  lociVar <- lapply(lociNames, function(locus){
    # Allele frequency and sample size for each ith population
    pi <- freqMat[,locus]
    ni <- sampMat[,locus]

    # Mean weighted allele frequency
    p.mean <- sum( pi * (ni/sum(ni)) )

    # Sum squared deviations, population
    ssp <- sum(ni * (pi - p.mean)^2)

    # Variance, gametes
    varg <- sum(ni * pi * (1-pi))

    # Mean squares variance components
    msp <- (1/(length(ni)-1)) * ssp
    msg <- (1/sum(ni-1)) * varg

    # Sample size correction factor
    nc <- (1/(length(ni)-1)) * ( sum(ni) - (sum(ni^2)/sum(ni)) )

    # Return the locus specific parameters
    return(data.table(LOCUS=locus
                      , NUMER=msp - msg
                      , DENOM=msp + (nc-1) * msg))
  })

  lociVar <- do.call('rbind', lociVar)

  return(lociVar)
}

