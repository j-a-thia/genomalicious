#' Calculate the variance components of Weir & Cockerham's FST

#' Calculates the numerator and denominator variance components for
#' Weir and Cockerham's FST (theta) for each locus from allele
#' frequency or individual genotype data.

#' @param dat Matrix: Allele frequencies for populations, or biallelic genotypes of
#' of individuals scores as integer counts of the Alt allele (0, 1, 2).
#' Populations or individuals are in rows, loci are in columns.
#' Column names are loci IDs and must exactly match those in argument \code{samp_size}.

#' @param input_type Character: One of two possible values: 'genos', calcualte
#' variance components from genotype matrix, or 'freqs', calculate variance
#' components from an allele frequency matrix.

#' @param samp_size Matrix: The sample size for each locus in each population.
#' Allows for different sample sizes at each locus, for example, if there
#' is missing data. Populations in rows, loci in columns. Rows must be
#' in the same order as rows in \code{dat}. Column names are loci IDs
#' and must all occur in \code{dat}.

#' @param pop_id  Charater: A vector of population IDs. Default is NULL and
#' is only required if \code{input_type=='genos'}. Must be the same length
#' as \code{nrow(dat)} and the order of values must match the order of rows in \code{dat}.
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
#' observed if indivdiual genotypes were used. \cr
#' The \code{samp_size} matrix can be used to allow for different sample sizes at
#' each locus, within a population, when \code{input_type=='freqs'}.
#' This facilitates calculation of variance components when there is missing data
#' at certain loci. In practise, for pool-seq experiments, it is hard to disentangle
#' the number of effectively contributing indivdiuals versus the number of pooled
#' diploids used to estimate the frequency at a specific locus in a specific population.
#'
#' @return A data table with columns $LOCUS (locus ID), $NUMER (the numerator for FST)
#' and $DENOM (the denominator for FST).
#'
#' @references
#' Weir, Cockerham (1984) Estimating F-statistics for the analysis of population structure. Evol. \cr
#' Weir, Hill (2002) Estimating F-statistics. Annu. Rev. Genet
#'
#' @examples
#' # FROM ALLELE FREQUENCIES
#' data(genomalicious_Freqs)
#' freqDat <- genomalicious_Freqs
#'
#' # Create a sample size matrix (order is important!)
#' freqSamps <- matrix(30, nrow=nrow(freqDat), ncol=ncol(freqDat), dimnames=list(NULL, colnames(freqDat)))
#'
#' # Variance components
#' freqVarcomps <- fstWC_varcomps(dat=freqDat, input_type='freqs', samp_size=freqSamps)
#' freqVarcomps
#'
#' # FROM INDIVIDUAL GENOTYPES
#' data(genomalicious_4pops)
#'
#' # Convert to matrix and get the first 8 loci
#' genoDat <- DT2Mat_genos(
#'     dat=genomalicious_4pops
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
fstWC_varcomps <- function(dat, input_type, samp_size=NULL, pop_id=NULL){
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

  # Get the locus names
  lociNames <- colnames(dat)

  # --------------------------------------------+
  # For allele frequencies
  # --------------------------------------------+
  if(input_type=='freqs'){
    # The number of populations
    r <- nrow(dat)

    # Iterate over each locus
    lociVar <- lapply(lociNames, function(locus){
      # Allele frequency and sample size for each ith population
      pi <- dat[, locus]
      ni <- samp_size[, locus]

      # Mean weighted allele frequency
      p.mean <- sum(ni * pi)/sum(ni)

      # Mean squares variance components
      msp <- (1/(r-1)) * sum(ni * (pi - p.mean)^2)
      msg <- (1/sum(ni-1)) * sum(ni * pi * (1-pi))

      # Sample size correction factor
      nc <- (1/(r-1)) * ( sum(ni) - (sum(ni^2)/sum(ni)) )

      # Return the locus specific parameters
      return(data.table(LOCUS=locus, NUMER=msp-msg, DENOM=msp+(nc-1)*msg))
    })
  }

  # --------------------------------------------+
  # For individual genotypes
  # --------------------------------------------+
  if(input_type=='genos'){
    popNames <- sort(unique(pop_id))

    # The number of populations
    r <- length(unique(pop_id))

    # Iterate over each locus
    lociVar <- lapply(lociNames, function(locus){
      # Get the allele frequencies, observed heterozygosity,
      # and sample size for each ith population
      gtDT <- na.omit(data.table(POP=pop_id, GT=dat[, locus]))
      piDT <- gtDT[, .( FREQ=sum(GT)/(length(GT)*2) ), by=POP]
      niDT <- gtDT[, .( SAMPS=length(GT) ), by=POP]
      hiDT <- gtDT[, .( HO=sum(GT==1)/length(GT) ), by=POP]
      pi <- piDT$FREQ
      names(pi) <- piDT$POP
      ni <- niDT$SAMPS
      names(ni) <- niDT$POP
      hi <- hiDT$HO
      names(hi) <- hiDT$POP

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
      return(data.table(LOCUS=locus, NUMER=a, DENOM=a+b+c))
    })
  }

  # Merge list items together
  lociVar <- do.call('rbind', lociVar)

  # Return data table
  return(lociVar)

  # .......... END
}

