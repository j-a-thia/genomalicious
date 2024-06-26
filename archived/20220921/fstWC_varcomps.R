#' Calculate the variance components of Weir & Cockerham's FST
#'
#' Calculates the numerator and denominator variance components for
#' Weir and Cockerham's FST (theta) for each locus from allele
#' frequency or individual genotype data.
#'
#' @param dat Matrix: Allele frequencies for populations, or biallelic genotypes
#' of individuals scored as integer counts of the Alt allele (0, 1, 2).
#' Populations or individuals are in rows, loci are in columns.
#' Column names are loci IDs and must exactly match those in argument \code{samp_size}.
#'
#' @param input_type Character: One of two possible values: 'genos', calcualte
#' variance components from genotype matrix, or 'freqs', calculate variance
#' components from an allele frequency matrix.
#'
#' @param samp_size Matrix: The sample size for each locus in each population.
#' Default is \code{NULL} and is required if \code{input_type=='freqs'}.
#' Allows for different sample sizes at each locus, for example, if there
#' is missing data. Populations in rows, loci in columns. Rows must be
#' in the same order as rows in \code{dat}. Column names are loci IDs
#' and must all occur in \code{dat}.
#'
#' @param pop_id  Charater: A vector of population IDs. Default is NULL and
#' is only required if \code{input_type=='genos'}. Must be the same length
#' as \code{nrow(dat)} and the order of values must match the order of rows in \code{dat}.
#'
#' @param num.cores Integer: Number of cores to use for estimation. Default is 1.
#'
#' @details
#' Estimation of the variance components from allele frequencies versus genotypes will
#' give slightly different values because of the different ways these are calculated.
#' Variance components for allele frequencies are calculated as described in
#' Weir and Hill (2002), where the numerator = MSP - MSG, and the
#' denominator = MSP + (nc - 1)MSG. Variance components for individual genotypes
#' are caculated as described in Weir and Cockerham (1984), where the numberator = a,
#' and the denominator = a + b + c. The methods are fundamentally similar, except in
#' estimation of variance components from individual genotypes, inclusion of observed
#' heterozygosity occurs in the calculation. Heterozygosity is not available when
#' only frequencies are assayed (e.g. in pool-seq), so its exclusion from calculations
#' produces these slight divergences in the variance components from those
#' observed if indivdiual genotypes were used. \cr\cr
#' The \code{samp_size} matrix can be used to allow for different sample sizes at
#' each locus, within a population, when \code{input_type=='freqs'}.
#' This facilitates calculation of variance components when there is missing data
#' at certain loci. In practise, for pool-seq experiments, it is hard to disentangle
#' the number of effectively contributing indivdiuals versus the number of pooled
#' diploids used to estimate the frequency at a specific locus in a specific population.
#'
#' @return A data table with columns \code{$LOCUS} (locus ID),
#' \code{$NUMER} (the numerator for FST) and \code{$DENOM} (the denominator for FST).
#'
#' @references
#' Weir, Cockerham (1984) Estimating F-statistics for the analysis of population structure. Evol. \cr
#' Weir, Hill (2002) Estimating F-statistics. Annu. Rev. Genet
#'
#' @examples
#' # FROM ALLELE FREQUENCIES
#' data(data_FreqsMat)
#' freqDat <- data_FreqsMat
#'
#' # Create a sample size matrix (order is important!)
#' freqSamps <- matrix(30, nrow=nrow(freqDat), ncol=ncol(freqDat), dimnames=list(NULL, colnames(freqDat)))
#'
#' # Variance components
#' freqVarcomps <- fstWC_varcomps(dat=freqDat, input_type='freqs', samp_size=freqSamps)
#' freqVarcomps
#'
#' # FROM INDIVIDUAL GENOTYPES
#' data(data_4pops)
#'
#' # Convert to matrix and get the first 8 loci
#' genoDat <- DT2Mat_genos(
#'     dat=data_4pops
#'     , sampCol='SAMPLE'
#'     , locusCol='LOCUS'
#'     , genoCol='GT')[, 1:8]
#'
#' # Must convert genotypes to integers
#' genoDat <- apply(genoDat, 2, genoscore_converter)
#'
#' # Create a population ID vector (order is important!)
#' # See that individuals in rows follow naming convention: Ind[pop].[sample].
#' # We can strip the characters after the '.' to get the population specific value.
#' genoPops <- gsub("\\..*", "", rownames(genoDat))
#'
#' # Variance components
#' genoVarcomps <- fstWC_varcomps(
#'     dat=genoDat
#'     , input_type='genos'
#'     , pop_id=genoPops)
#' genoVarcomps
#'
#' @export
fstWC_varcomps <- function(dat, input_type, samp_size=NULL, pop_id=NULL, num.cores=1){
  # BEGIN ..........

  # --------------------------------------------+
  # Assertions and environment
  # --------------------------------------------+
  # Environment checks
  require(data.table)

  # Check input type
  if(!input_type %in% c('genos', 'freqs')){
    stop("Argument input type is incorrectly specified, see ?fstWC_varcomps")
  }

  # If there are no locus names, stop.
  if(length(colnames(dat))!=ncol(dat)){
    stop("All colnames of argument dat (the loci) need to named, see ?fstWC_varcomps")
  }

  # If using genotypes, but no population IDs specified...
  if(input_type=='genos' & is.null(pop_id)==TRUE){
    stop("Argument input_type=='genos', but pop_id is not specified, see ?fstWC_varcomps")
  }

  # If using frequencies, but no sample size matrix
  if(input_type=='freqs' & is.null(samp_size)==TRUE){
    stop("Argument input_type=='freqs', but samp_size is not specified, see ?fstWC_varcomps")
  }

  # If using frequencies:
  # (1) Check columns are equal for dat and samp_size,
  # (2) Make sure column names are the same for dat and samp_size,
  # (3) For good measure, make columns be in the same order.
  if(input_type=='freqs'){
    if(ncol(dat)!=ncol(samp_size)){
      stop("The number columns (loci) must be equal in dat and samp_size, see ?fstWC_varcomps")
    }
    if(sum(colnames(samp_size) %in% colnames(dat))!=ncol(dat)){
      stop("All column names in dat must be in samp_size, see ?fstWC_varcomps")
    }

    samp_size <- samp_size[,colnames(dat)]
  }

  if(class(dat[,1])=='character'){
    warning('Argument dat contains genotypes coded as characters, attempting to conversion to integer allele counts, see ?fstWC_varcomps')
    dat <- apply(dat, 2, genoscore_converter)
  }

  # Get the locus names
  lociNames <- colnames(dat)

  # --------------------------------------------+
  # For allele frequencies
  # --------------------------------------------+
  FUN_freqs_varcomp <- function(pi, ni, r){
    # Mean weighted allele frequency
    p.mean <- sum(ni * pi)/sum(ni)

    # Mean squares variance components
    msp <- (1/(r-1)) * sum(ni * (pi - p.mean)^2)
    msg <- (1/sum(ni-1)) * sum(ni * pi * (1-pi))

    # Sample size correction factor
    nc <- (1/(r-1)) * ( sum(ni) - (sum(ni^2)/sum(ni)) )

    # Output as data.table
    data.table(NUMER=msp-msg, DENOM=msp+(nc-1)*msg)
  }

  if(input_type=='freqs'){
    # The number of populations
    r <- nrow(dat)

    # Iterate over each locus
    # ... Single core
    if(num.cores==1){
      lociVar <- lapply(lociNames, function(locus){
        data.table(
          LOCUS=locus,
          FUN_freqs_varcomp(pi=dat[, locus], ni=samp_size[, locus], r=r)
        )
      })
    }

    # .... Multiple cores
    if(num.cores>1){
      lociVar <- foreach(locus=lociNames) %dopar% {
        require(data.table)

        data.table(
          LOCUS=locus,
          FUN_freqs_varcomp(pi=dat[, locus], ni=samp_size[, locus], r=r)
        )
      }
    }
  }

  # --------------------------------------------+
  # For individual genotypes
  # --------------------------------------------+
  FUN_genos_varcomps <- function(ni, pi, hi, r){
    # Reorder, just to be sure
    ni <- ni[names(pi)]
    hi <- hi[names(pi)]

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

    # Return as numerator and denominator
    return(data.table(NUMER=a, DENOM=a+b+c))
  }

  if(input_type=='genos'){
    popNames <- sort(unique(pop_id))

    # The number of populations
    r <- length(unique(pop_id))

    # Create allele frequency, sample size, and heterozygosity matrices.
    # .... Get the values for each population
    allpopVals <- lapply(popNames, function(pop){
      popGenos <- dat[which(pop_id==pop),]

      # Get values for each locus for the population
      locpopVals <- apply(popGenos, 2, function(xx){
        yy <- na.omit(xx)
        n <- length(yy)
        p <- sum(yy)/(length(yy)*2)
        h <- sum(yy==1)/length(yy)
        return(list(n=n, p=p, h=h))
      })

      # Combine information from all loci into single row
      n <- do.call('cbind', lapply(locpopVals, function(xx){ xx$n }))
      p <- do.call('cbind', lapply(locpopVals, function(xx){ xx$p }))
      h <- do.call('cbind', lapply(locpopVals, function(xx){ xx$h }))

      rownames(n) <- pop
      rownames(p) <- pop
      rownames(h) <- pop

      return(list(n=n, p=p, h=h))
    })

    # .... Now combine populations
    sampMat <- do.call('rbind', lapply(allpopVals, function(xx){ xx$n }))
    freqMat <- do.call('rbind', lapply(allpopVals, function(xx){ xx$p }))
    hetMat <- do.call('rbind', lapply(allpopVals, function(xx){ xx$h }))

    # .... Locus names
    lociNames <- colnames(freqMat)

    # Iterate over each locus
    # ... Single core
    if(num.cores==1){
      lociVar <- lapply(lociNames, function(locus){
        data.table(
          LOCUS=locus,
          FUN_genos_varcomps(
            ni = sampMat[, locus],
            pi = freqMat[, locus],
            hi = hetMat[, locus],
            r = r
          )
        )
      })
    }

    # ... Multiple cores
    if(num.cores>1){
      lociVar <- foreach(locus=lociNames) %dopar% {
        require(data.table)

        data.table(
          LOCUS=locus,
          FUN_genos_varcomps(
            ni = sampMat[, locus],
            pi = freqMat[, locus],
            hi = hetMat[, locus],
            r = r
          )
        )
      }
    }

  }

  # --------------------------------------------+
  # Return the variance components
  # --------------------------------------------+
  # Merge list items together
  lociVar <- do.call('rbind', lociVar)

  # Return data table
  return(lociVar)

  # .......... END
}

